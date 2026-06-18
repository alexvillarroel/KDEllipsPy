"""
validar_egf_full.py
===================
Validacion EGF extendida a TODAS las estaciones comunes entre la descarga
FDSN (CX PB* + C1.AF01) y el mainshock de registros_20200603 (CSN, acel txt).

Mismo criterio que validar_egf.py (CC mainshock<->EGF, P en Z / S en T,
banda 0.1-0.3 Hz causal), pero el mainshock se toma de registros_20200603.

Salidas:
  egf/estaciones_validadas.csv  (sta, lat, lon, snr, cc, dt_align para Fase 2)
  egf/validacion_full_<tag>.png
Uso:
  python codigos/validar_egf_full.py [tag]    (default 20200616T175441, el M4.9)
"""
import os
import sys
import csv
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)                                  # .../egf
EV2020 = os.path.dirname(BASE)                                # .../2020-06-03
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(EV2020, "codigos"))

from io_sismica import read_txt_trace                          # noqa: E402
from validar_egf import (MAIN, F1, F2, FS, WIN_P, WIN_S, MAX_LAG,   # noqa: E402
                         tt, preparar, rotar, cc_max, cargar_egf)

TXT_ROOT = os.path.join(EV2020, "registros_20200603")
# Estaciones comunes registros_20200603 <-> descarga EGF.
STATIONS = ["PB01", "PB02", "PB03", "PB05", "PB06", "PB09", "PB19", "AF01"]
CC_MIN = 0.5     # umbral para pasar a la inversion (Fase 2)


def cargar_main_txt(station):
    """Mainshock desde registros (acel 200 Hz) -> velocidad filtrada comun."""
    comps, coords = {}, None
    for c in "NEZ":
        m = glob.glob(os.path.join(TXT_ROOT, station, f"*{station}*N{c}.txt")) or \
            glob.glob(os.path.join(TXT_ROOT, station, f"*-{station}-*{c}.txt"))
        if not m:
            return None, None
        tr = read_txt_trace(m[0])
        coords = (tr.stats.coordinates["latitude"], tr.stats.coordinates["longitude"])
        comps[c] = preparar(tr, MAIN["time"], integrar=True)
    return comps, coords


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "20200616T175441"
    with open(os.path.join(BASE, "candidatos_ranking.csv")) as f:
        rows = list(csv.DictReader(f))
    row = next(r for r in rows
               if UTCDateTime(r["time"]).strftime("%Y%m%dT%H%M%S") == tag)
    ev = {"time": UTCDateTime(row["time"]), "lat": float(row["lat"]),
          "lon": float(row["lon"]), "depth_km": float(row["depth_km"])}
    print(f"Validacion FULL del EGF {tag} (M{row['mag']}) en {len(STATIONS)} "
          f"estaciones, banda {F1}-{F2} Hz")

    res, figs = [], []
    for sta in STATIONS:
        mn, coords = cargar_main_txt(sta)
        egf = cargar_egf(tag, ev, sta)
        if mn is None or egf is None:
            print(f"  {sta}: faltan datos (main={mn is not None}, egf={egf is not None})")
            continue
        la, lo = coords
        d_m, az, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], la, lo)
        deg = kilometers2degrees(d_m / 1000.0)
        tp_m, ts_m = tt(deg, MAIN["depth_km"])
        d2, az2, _ = gps2dist_azimuth(ev["lat"], ev["lon"], la, lo)
        tp_e, ts_e = tt(kilometers2degrees(d2 / 1000.0), ev["depth_km"])

        def win(data, t_on, dur):
            i0 = int((t_on + 60.0) * FS)
            return data[max(0, i0):i0 + int(dur * FS)]

        rn_m, tt_m = rotar(mn["N"].data, mn["E"].data, az)
        rn_e, tt_e = rotar(egf["N"].data, egf["E"].data, az2)
        noise = egf["Z"].data[:int(50 * FS)]
        nrm = max(np.sqrt(np.mean(noise ** 2)), 1e-20)
        snr_p = np.sqrt(np.mean(win(egf["Z"].data, tp_e, WIN_P) ** 2)) / nrm
        noise_t = tt_e[:int(50 * FS)]
        snr_s = np.sqrt(np.mean(win(tt_e, ts_e, WIN_S) ** 2)) / \
            max(np.sqrt(np.mean(noise_t ** 2)), 1e-20)

        ccP, lagP = cc_max(win(mn["Z"].data, tp_m, WIN_P),
                           win(egf["Z"].data, tp_e, WIN_P), FS, MAX_LAG)
        ccS, lagS = cc_max(win(tt_m, ts_m - 5.0, WIN_S),
                           win(tt_e, ts_e - 5.0, WIN_S), FS, MAX_LAG)

        res.append({"sta": sta, "lat": la, "lon": lo,
                    "dist_km": round(d_m / 1000, 1),
                    "snrP": round(float(snr_p), 1), "snrS": round(float(snr_s), 1),
                    "ccP_Z": round(ccP, 2), "lagP_s": round(lagP, 1),
                    "ccS_T": round(ccS, 2), "lagS_s": round(lagS, 1),
                    "dt_align_s": round(ts_m - ts_e, 2),
                    "valida": int(ccS >= CC_MIN)})
        figs.append((sta, mn["Z"].data, egf["Z"].data, tp_m, tp_e,
                     tt_m, tt_e, ts_m, ts_e))
        print(f"  {sta}: SNR_P={snr_p:.1f} SNR_S={snr_s:.1f}  CC_P={ccP:.2f}"
              f" ({lagP:+.1f}s)  CC_S={ccS:.2f} ({lagS:+.1f}s)"
              f"  {'OK' if ccS >= CC_MIN else 'descartada'}")

    out = os.path.join(BASE, "estaciones_validadas.csv"
                       if tag == "20200616T175441"
                       else f"estaciones_validadas_{tag}.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader(); w.writerows(res)
    n_ok = sum(r["valida"] for r in res)
    print(f"\n{n_ok}/{len(res)} estaciones validas (CC_S>={CC_MIN}) -> {out}")

    # figura
    t = np.arange(len(figs[0][1])) / FS - 60.0
    fig, axes = plt.subplots(len(figs), 2, figsize=(13, 2.2 * len(figs)),
                             sharex=False)
    axes = np.atleast_2d(axes)
    for i, (sta, zm, ze, tpm, tpe, tm, te_, tsm, tse) in enumerate(figs):
        for j, (m_, e_, on_m, on_e, lbl) in enumerate(
                [(zm, ze, tpm, tpe, "Z (P)"), (tm, te_, tsm, tse, "T (S)")]):
            ax = axes[i, j]
            ax.plot(t, m_ / max(np.max(np.abs(m_)), 1e-20), "k", lw=0.8)
            ax.plot(t + (on_m - on_e), e_ / max(np.max(np.abs(e_)), 1e-20),
                    "r", lw=0.8, alpha=0.8)
            ax.axvline(on_m, color="b", ls=":", lw=0.7)
            ax.set_ylabel(f"{sta}\n{lbl}", fontsize=8)
            ax.set_xlim(on_m - 20, on_m + 80)
            ax.tick_params(labelsize=7)
    axes[-1, 0].set_xlabel("t desde origen (s)"); axes[-1, 1].set_xlabel("t (s)")
    fig.suptitle(f"Validacion FULL EGF {tag} (negro=mainshock, rojo=EGF alineada)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fp = os.path.join(BASE, f"validacion_full_{tag}.png")
    fig.savefig(fp, dpi=140)
    print(f"figura -> {fp}")


if __name__ == "__main__":
    main()
