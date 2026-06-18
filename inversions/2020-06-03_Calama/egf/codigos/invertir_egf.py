"""
invertir_egf.py
===============
FASE 2: inversion cinematica de la elipse usando la EGF validada
(2020-06-16 M4.9) en la banda 0.1-0.3 Hz.

Cadena:
  modelo 7 params -> FaultGeometry (kdellipspy, plano de input.ctl)
    -> momentos m_k + tiempos de ruptura t_r(k)  (to_axitra_hist)
    -> posiciones lat/lon/z                       (to_axitra_sources)
    -> EGFForward.predict  (shift-and-sum con directividad)
    -> MisfitCalculator.l2_misfit  (ventanas P/S de 20 s, igual que AXITRA)

Observado : mainshock desde registros_20200603 (acel -> vel causal 0.1-0.3).
EGF       : mseed FDSN, respuesta removida a VEL, mismo filtro causal,
            alineada por diferencia de tiempo de viaje S (taup); los lags
            residuales (directividad del mainshock) son la SENAL a invertir.
M0_egf    : de Mw catalogo, 10^(1.5*Mw + 9.05).
Estaciones: las CC_S >= 0.5 de estaciones_validadas.csv (7).

Salidas en egf/inversion[<sufijo>]/: inversion_result.joblib, ajuste_ondas.png,
appraisal_corner_egf.png, INFORME_fase2.md
Uso: python codigos/invertir_egf.py [--smoke] [--seed N] [--excluir STA]
                                    [--tag TAG] [--out SUFIJO] [--dmax M]
                                    [--banda F1,F2] [--csv NOMBRE]
                                    [--lagmax S] [--sinP STA1,STA2]

--banda  : banda del filtro causal (default 0.1,0.3).
--csv    : CSV de estaciones validadas (default segun tag/banda).
--lagmax : lag estatico maximo (s) buscado por ventana y estacion dentro del
           misfit (P y S independientes; 0 = sin lags, misfit clasico).
           Remueve errores de alineacion taup/localizacion; la directividad
           debe ajustarse entonces por forma/duracion, no por lag absoluto.
--sinP   : estaciones cuyas ventanas P (R,Z) se excluyen del misfit
           (CC_P baja: solo inyectan ruido).
"""
import os
import sys
import csv
import glob
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import joblib
import obspy
from obspy import UTCDateTime, read_inventory
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)                       # .../egf
EV2020 = os.path.dirname(BASE)                     # .../2020-06-03
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(EV2020, "codigos"))

from io_sismica import read_txt_trace                          # noqa: E402
from egf_forward import EGFForward                             # noqa: E402
from validar_egf import tt                                     # noqa: E402

from kdellipspy import ConfigParser                            # noqa: E402
from kdellipspy.core.forward_model import AxitraForwardModel   # noqa: E402
from kdellipspy.inversion.base import MisfitCalculator, NAModel, NAResult  # noqa: E402
from neighpy import NASearcher                                 # noqa: E402

# ------------------------------------------------------------------ config
def _arg(flag, default=None):
    if flag in sys.argv:
        return sys.argv[sys.argv.index(flag) + 1]
    return default


TAG = _arg("--tag", "20200616T175441")
EGF_EVENTS = {
    "20200616T175441": {"time": UTCDateTime("2020-06-16T17:54:41.000"),
                        "lat": -23.307, "lon": -68.578, "depth_km": 108.0,
                        "mw": 4.89},
    "20201025T213242": {"time": UTCDateTime("2020-10-25T21:32:42.030"),
                        "lat": -23.241, "lon": -68.537, "depth_km": 118.0,
                        "mw": 4.6},
}
EGF_EV = EGF_EVENTS[TAG]
M0_EGF = 10 ** (1.5 * EGF_EV["mw"] + 9.05)         # ~2.4e16 Nm (M4.9)
SEED = int(_arg("--seed", "11"))
EXCLUIR = (_arg("--excluir", "") or "").split(",") if _arg("--excluir") else []

MAIN = {"time": UTCDateTime("2020-06-03T07:35:34.000"),
        "lat": -23.247, "lon": -68.530, "depth_km": 123.4}

F1, F2 = (float(v) for v in _arg("--banda", "0.1,0.3").split(","))
CSV_VALIDADAS = _arg("--csv")   # None -> default segun tag
LAG_MAX = float(_arg("--lagmax", "0.0"))
SIN_P = (_arg("--sinP", "") or "").split(",") if _arg("--sinP") else []
DT = 0.25                       # FS=4 Hz (entero para MisfitCalculator)
NPTS = 512                      # 0..128 s
WIN_MISFIT = 20.0               # s, igual que las corridas AXITRA
C_PHASE = 4.6                   # km/s, vel. de fase efectiva (S en la fuente)

DMAX_MAX = float(_arg("--dmax", "4.0"))
BOUNDS = [(2.0, 13.0), (2.0, 13.0), (0.0, 1.0), (0.0, 1.0),
          (0.0, 1.0), (0.5, DMAX_MAX), (1.5, 4.5)]
PARAM_NAMES = ["a1 (km)", "a2 (km)", "theta (x pi)", "np (frac)",
               "tp (x 2pi)", "dmax (m)", "vr (km/s)"]

TXT_ROOT = os.path.join(EV2020, "registros_20200603")
OUTDIR = os.path.join(BASE, "inversion" + (_arg("--out") or ""))


# ------------------------------------------------------------------ datos
def _preparar(tr, t0, integrar, shift_s=0.0):
    """detrend(+integra)+filtro causal y rejilla comun NPTS x DT desde t0."""
    tr = tr.copy()
    tr.detrend("linear"); tr.detrend("demean")
    tr.taper(max_percentage=0.05, type="cosine")
    if integrar:
        tr.integrate(method="cumtrapz")
        tr.detrend("linear")
    tr.filter("bandpass", freqmin=F1, freqmax=F2, corners=4, zerophase=False)
    t_ini = t0 - shift_s                     # shift>0 retrasa la senal
    tr.trim(t_ini - 5.0, t_ini + NPTS * DT + 5.0, pad=True, fill_value=0.0)
    tr.interpolate(sampling_rate=1.0 / DT, starttime=t_ini, npts=NPTS)
    return tr.data[:NPTS].astype(np.float64)


def cargar_estaciones():
    """Lee estaciones_validadas[_tag].csv y arma observado + EGF + azi_times."""
    csv_name = CSV_VALIDADAS or ("estaciones_validadas.csv"
                                 if TAG == "20200616T175441"
                                 else f"estaciones_validadas_{TAG}.csv")
    with open(os.path.join(BASE, csv_name)) as f:
        rows = [r for r in csv.DictReader(f)
                if int(r["valida"]) == 1 and r["sta"] not in EXCLUIR]

    inv = read_inventory(os.path.join(BASE, TAG, "RAW", "inventory.xml"))
    obs, egf, azi_rows, keep = [], [], [], []
    for r in rows:
        sta, la, lo = r["sta"], float(r["lat"]), float(r["lon"])
        # mainshock (acel txt -> vel)
        comps_m = {}
        for c in "NEZ":
            m = glob.glob(os.path.join(TXT_ROOT, sta, f"*-{sta}-*{c}.txt"))
            if not m:
                break
            comps_m[c] = _preparar(read_txt_trace(m[0]), MAIN["time"], True)
        if len(comps_m) != 3:
            print(f"  [!] {sta}: mainshock incompleto"); continue
        # EGF (mseed -> VEL fisica) alineada por dt de viaje S
        fn = glob.glob(os.path.join(BASE, TAG, "RAW", f"*.{sta}.mseed"))
        if not fn:
            print(f"  [!] {sta}: sin EGF"); continue
        st = obspy.read(fn[0])
        for pref in ("HH", "HL", "HN", "BH"):
            sel = st.select(channel=f"{pref}?")
            if len(sel) >= 3:
                st = sel; break
        try:
            st.remove_response(inventory=inv, output="VEL",
                               pre_filt=(0.02, 0.05, 8.0, 9.0))
        except Exception as ex:  # noqa: BLE001
            print(f"  [!] {sta}: respuesta fallo ({str(ex)[:50]})"); continue
        shift = float(r["dt_align_s"])       # ts_main - ts_egf (taup)
        comps_e = {}
        for tr in st:
            c = tr.stats.channel[-1].replace("1", "N").replace("2", "E")
            comps_e[c] = _preparar(tr, EGF_EV["time"], False, shift_s=-shift)
        if not all(c in comps_e for c in "NEZ"):
            print(f"  [!] {sta}: EGF incompleta"); continue

        obs.append([comps_m["N"], comps_m["E"], comps_m["Z"]])
        egf.append([comps_e["N"], comps_e["E"], comps_e["Z"]])
        d_m, az, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], la, lo)
        tp, ts = tt(kilometers2degrees(d_m / 1000.0), MAIN["depth_km"])
        azi_rows.append([math.radians(az), tp - 1.0, ts - 1.0])
        keep.append({"sta": sta, "lat": la, "lon": lo,
                     "dist_km": d_m / 1000.0, "ccS": float(r["ccS_T"])})
        print(f"  {sta}: OK (dist {d_m/1000:.0f} km, CC_S {r['ccS_T']}, "
              f"align {shift:+.2f} s)")
    return (np.array(obs), np.array(egf), np.array(azi_rows), keep)


# ------------------------------------------------------------------ misfit
class MisfitLagEstatico:
    """L2 normalizado en ventanas P (R+Z) y S (T), mismas convenciones de
    indices que kdellipspy MisfitCalculator, pero con lag estatico por
    estacion y ventana (P y S independientes) buscado en +-lag_max dentro de
    cada evaluacion: el sintetico se desplaza al lag que minimiza el residuo.
    Con lag_max=0 reproduce el misfit clasico. use_p[j]=False excluye las
    ventanas P de la estacion j (CC_P baja)."""

    def __init__(self, observed, time_array, azi_times_array, time_window_s,
                 lag_max_s=0.0, use_p=None):
        self.obs = np.asarray(observed, float)
        nsta, _, npts = self.obs.shape
        dt = float(time_array[1] - time_array[0])
        samp = max(1, int(np.rint(1.0 / dt)))
        self.dt = dt
        win = max(1, int(np.rint(time_window_s * samp)))
        self.nlag = int(np.rint(lag_max_s / dt))
        self.use_p = (np.ones(nsta, bool) if use_p is None
                      else np.asarray(use_p, bool))
        azi = np.asarray(azi_times_array, float)
        self.az = azi[:, 0]
        self.win_idx, self.obs_win, self.den = [], [], 0.0
        for j in range(nsta):
            start_p = int(np.rint(azi[j, 1])) * samp
            start_s = int(np.rint(azi[j, 2])) * samp
            kp0, kp1 = max(0, start_p - 1), min(npts - 1, start_p + win - 1)
            ks0, ks1 = max(0, start_s - 1), min(npts - 1, start_s + win - 1)
            ca, sa = math.cos(self.az[j]), math.sin(self.az[j])
            r_o = self.obs[j, 0] * ca + self.obs[j, 1] * sa
            t_o = self.obs[j, 1] * ca - self.obs[j, 0] * sa
            self.win_idx.append((kp0, kp1, ks0, ks1))
            self.obs_win.append((r_o[kp0:kp1 + 1], self.obs[j, 2, kp0:kp1 + 1],
                                 t_o[ks0:ks1 + 1]))
            if self.use_p[j]:
                self.den += float(np.sum(r_o[kp0:kp1 + 1] ** 2))
                self.den += float(np.sum(self.obs[j, 2, kp0:kp1 + 1] ** 2))
            self.den += float(np.sum(t_o[ks0:ks1 + 1] ** 2))

    @staticmethod
    def _ventana(trace, k0, k1, lag):
        """trace[k0-lag : k1-lag+1] con relleno de ceros fuera de rango."""
        npts = len(trace)
        s0, s1 = k0 - lag, k1 - lag
        out = np.zeros(k1 - k0 + 1)
        a, b = max(0, s0), min(npts - 1, s1)
        if b >= a:
            out[a - s0:b - s0 + 1] = trace[a:b + 1]
        return out

    def l2_misfit(self, synth, return_lags=False):
        nsta = self.obs.shape[0]
        num = 0.0
        lags = np.zeros((nsta, 2))
        rango = range(-self.nlag, self.nlag + 1)
        for j in range(nsta):
            ca, sa = math.cos(self.az[j]), math.sin(self.az[j])
            r_s = synth[j, 0] * ca + synth[j, 1] * sa
            t_s = synth[j, 1] * ca - synth[j, 0] * sa
            kp0, kp1, ks0, ks1 = self.win_idx[j]
            r_o, z_o, t_o = self.obs_win[j]
            if self.use_p[j]:
                mejor, ml = None, 0
                for L in rango:
                    nj = float(np.sum((r_o - self._ventana(r_s, kp0, kp1, L)) ** 2)
                               + np.sum((z_o - self._ventana(synth[j, 2], kp0,
                                                             kp1, L)) ** 2))
                    if mejor is None or nj < mejor:
                        mejor, ml = nj, L
                num += mejor
                lags[j, 0] = ml * self.dt
            mejor, ml = None, 0
            for L in rango:
                nj = float(np.sum((t_o - self._ventana(t_s, ks0, ks1, L)) ** 2))
                if mejor is None or nj < mejor:
                    mejor, ml = nj, L
            num += mejor
            lags[j, 1] = ml * self.dt
        mis = num / max(self.den, 1e-30)
        return (mis, lags) if return_lags else mis


# ------------------------------------------------------------------ main
def main():
    smoke = "--smoke" in sys.argv
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"EGF Fase 2 — {TAG} M{EGF_EV['mw']} (M0_egf={M0_EGF:.2e} Nm), "
          f"banda {F1}-{F2} Hz, dt={DT}")

    obs, egf, azi, stas = cargar_estaciones()
    nsta = len(stas)
    print(f"\n{nsta} estaciones | observado {obs.shape} | egf {egf.shape}")
    time_arr = np.arange(NPTS) * DT

    use_p = np.array([s["sta"] not in SIN_P for s in stas])
    if LAG_MAX > 0 or SIN_P:
        mc = MisfitLagEstatico(obs, time_arr, azi, WIN_MISFIT,
                               lag_max_s=LAG_MAX, use_p=use_p)
        print(f"misfit: lags estaticos +-{LAG_MAX} s; sin ventana P: "
              f"{SIN_P or 'ninguna'}")
    else:
        mc = MisfitCalculator(observed_waveforms=obs, time_array=time_arr,
                              azi_times_array=azi, time_window_s=WIN_MISFIT)

    cfg = ConfigParser(os.path.join(EV2020, "input.ctl"))
    fm = AxitraForwardModel.from_config(cfg)
    fwd = EGFForward(egf, DT, M0_EGF,
                     [s["lat"] for s in stas], [s["lon"] for s in stas],
                     (MAIN["lat"], MAIN["lon"], MAIN["depth_km"]),
                     c_phase=C_PHASE)

    n_eval = [0]

    def objective(x):
        n_eval[0] += 1
        try:
            geom = fm.build_geometry_with_ellipse_slip(np.asarray(x, float))
            hist = geom.to_axitra_hist()
            srcs = geom.to_axitra_sources(latlon=True)
            if len(hist) == 0:
                return 1.0e3
            synth = fwd.predict(srcs[:, 1], srcs[:, 2], srcs[:, 3] / 1000.0,
                                hist[:, 1], x[6], t_rup=hist[:, 7])
            mis = mc.l2_misfit(synth)
        except Exception:  # noqa: BLE001
            return 1.0e3
        if n_eval[0] % 200 == 0:
            print(f"  eval {n_eval[0]:5d}  misfit {mis:.4f}", flush=True)
        return float(mis)

    if smoke:
        ni, ns, n, nr = 50, 20, 3, 4
    else:
        ni, ns, n, nr = 500, 50, 50, 10        # 3000 evals
    np.random.seed(SEED)
    print(f"NA: ni={ni} ns={ns} n={n} nr={nr} -> {ni + ns*n} evals "
          f"(seed={SEED}, excluidas={EXCLUIR or 'ninguna'})")
    searcher = NASearcher(objective, ns=ns, nr=nr, ni=ni, n=n,
                          bounds=tuple(BOUNDS))
    searcher.run()

    samples = np.asarray(searcher.samples)
    objs = np.asarray(searcher.objectives)
    ib = int(np.argmin(objs))
    best, best_mis = samples[ib], float(objs[ib])
    print(f"\nMejor misfit: {best_mis:.4f}")
    for nm, v in zip(PARAM_NAMES, best):
        print(f"  {nm:14s} = {v:7.3f}")

    # mejor sintetico + momento
    geom = fm.build_geometry_with_ellipse_slip(best)
    hist = geom.to_axitra_hist(); srcs = geom.to_axitra_sources(latlon=True)
    synth = fwd.predict(srcs[:, 1], srcs[:, 2], srcs[:, 3] / 1000.0,
                        hist[:, 1], best[6], t_rup=hist[:, 7])
    m0 = float(hist[:, 1].sum())
    mw = (2.0 / 3.0) * (math.log10(m0) - 9.05)
    print(f"M0={m0:.3e} Nm  Mw={mw:.2f}")
    best_lags = None
    if isinstance(mc, MisfitLagEstatico):
        _, best_lags = mc.l2_misfit(synth, return_lags=True)
        print("lags estaticos del mejor modelo (s):")
        for s, (lp, ls) in zip(stas, best_lags):
            pp = f"{lp:+.2f}" if s["sta"] not in SIN_P else "  -- "
            print(f"  {s['sta']:5s}  P {pp}  S {ls:+.2f}")

    models = [NAModel(model=np.asarray(s, float), misfit=float(o), iteration=0)
              for s, o in zip(samples, objs)]
    result = NAResult(models, param_names=PARAM_NAMES,
                      extra_metadata={"algorithm": "NA-EGF", "egf_tag": TAG,
                                      "m0_egf_nm": M0_EGF, "band": [F1, F2],
                                      "stations": [s["sta"] for s in stas],
                                      "lag_max_s": LAG_MAX, "sin_p": SIN_P,
                                      "best_lags": best_lags},
                      best_synthetics=synth, observed=obs, time=time_arr,
                      azi_times_array=azi)
    joblib.dump(result, os.path.join(OUTDIR, "inversion_result.joblib"))

    # ---------------- figura de ajuste (R/T/Z con ventanas) ----------------
    fig, axes = plt.subplots(nsta, 3, figsize=(13, 1.9 * nsta), sharex=True)
    axes = np.atleast_2d(axes)
    for i, s in enumerate(stas):
        az = azi[i, 0]
        r_o = obs[i, 0] * np.cos(az) + obs[i, 1] * np.sin(az)
        t_o = obs[i, 1] * np.cos(az) - obs[i, 0] * np.sin(az)
        r_s = synth[i, 0] * np.cos(az) + synth[i, 1] * np.sin(az)
        t_s = synth[i, 1] * np.cos(az) - synth[i, 0] * np.sin(az)
        for j, (o, sy, lbl, t_on) in enumerate(
                [(r_o, r_s, "R", azi[i, 1]), (t_o, t_s, "T", azi[i, 2]),
                 (obs[i, 2], synth[i, 2], "Z", azi[i, 1])]):
            ax = axes[i, j]
            ax.plot(time_arr, o, "k", lw=0.9)
            ax.plot(time_arr, sy, "r", lw=0.9, alpha=0.85)
            ax.axvspan(t_on, t_on + WIN_MISFIT, color="y", alpha=0.18, lw=0)
            if best_lags is not None:
                lag = best_lags[i, 1] if lbl == "T" else best_lags[i, 0]
                txt = (f"lag {lag:+.2f}s" if (lbl == "T" or use_p[i])
                       else "P excl.")
                ax.text(0.99, 0.93, txt, transform=ax.transAxes, ha="right",
                        va="top", fontsize=6, color="0.4")
            ax.set_ylabel(f"{s['sta']}\n{lbl}", fontsize=7)
            ax.tick_params(labelsize=6)
    for j in range(3):
        axes[-1, j].set_xlabel("t (s)")
    fig.suptitle(f"EGF inversion {TAG} — misfit {best_mis:.3f} "
                 f"(negro=obs, rojo=synth; banda {F1}-{F2} Hz)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(OUTDIR, "ajuste_ondas.png"), dpi=140)
    plt.close(fig)

    # ---------------- appraisal ----------------
    if not smoke:
        try:
            result.run_appraisal(n_resample=20000, verbose=False, save=False,
                                 bounds=np.array(BOUNDS, float))
            result.plot_appraisal(
                save_path=os.path.join(OUTDIR, "appraisal_corner_egf.png"),
                show=False)
            print("appraisal -> inversion/appraisal_corner_egf.png")
        except Exception as ex:  # noqa: BLE001
            print(f"appraisal fallo: {ex}")

    # ---------------- informe ----------------
    with open(os.path.join(OUTDIR, "INFORME_fase2.md"), "w") as f:
        f.write(f"# Inversion EGF — Calama 2020 (banda {F1}-{F2} Hz)\n\n")
        f.write(f"EGF: {TAG} M{EGF_EV['mw']} (M0_egf={M0_EGF:.2e} Nm, "
                f"de Mw catalogo).\n")
        f.write(f"Estaciones ({nsta}): "
                + ", ".join(s["sta"] for s in stas) + "\n\n")
        f.write(f"**Mejor misfit: {best_mis:.4f}**  (M0={m0:.2e} Nm, "
                f"Mw={mw:.2f})\n\n")
        f.write("| param | valor |\n|---|---|\n")
        for nm, v in zip(PARAM_NAMES, best):
            f.write(f"| {nm} | {v:.3f} |\n")
        if LAG_MAX > 0 and best_lags is not None:
            f.write(f"Lags estaticos (max +-{LAG_MAX} s) del mejor modelo; "
                    f"ventanas P excluidas: {', '.join(SIN_P) or 'ninguna'}\n\n")
            f.write("| sta | lag P (s) | lag S (s) |\n|---|---|---|\n")
            for s, (lp, ls) in zip(stas, best_lags):
                pp = f"{lp:+.2f}" if s["sta"] not in SIN_P else "--"
                f.write(f"| {s['sta']} | {pp} | {ls:+.2f} |\n")
            f.write("\n")
        f.write("\nNotas: EGF alineada solo por dt de viaje S (taup); los lags"
                " de directividad quedan como senal. M0_egf incierto ->"
                " trade-off directo con dmax (geometria no afectada). P (R,Z)"
                " con CC moderada; la senal dominante es S (T).\n")
        f.write("\n![ajuste](ajuste_ondas.png)\n")
        if not smoke:
            f.write("\n![appraisal](appraisal_corner_egf.png)\n")
    print(f"informe -> {os.path.join(OUTDIR, 'INFORME_fase2.md')}")


if __name__ == "__main__":
    main()
