"""
comparar_cc_bandas.py
=====================
Correlacion cruzada mainshock <-> EGF (M4.9, 20200616T175441) en dos bandas:
la actual de inversion (0.1-0.3 Hz) y la candidata (0.15-0.4 Hz), motivada
por el SNR espectral (espectros_snr.pdf: el microsismo cae justo en 0.1-0.3
y el SNR de la EGF mejora ~3x en 0.15-0.4).

Mismo criterio que validar_egf_full.py (P en Z, S en T, mainshock txt CSN
integrado a velocidad, EGF respuesta removida a VEL, filtro causal), pero
parametrizado en banda. NO sobreescribe estaciones_validadas.csv.

Salidas:
  egf/cc_bandas.csv                       (tabla larga: banda x estacion)
  egf/validacion_cc_0.15-0.40.png         (formas de onda en la banda nueva)
Uso:
  python codigos/comparar_cc_bandas.py
"""
import os
import sys
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)                                   # .../egf
EV2020 = os.path.dirname(BASE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(EV2020, "codigos"))

import validar_egf as ve                                        # noqa: E402
from validar_egf import (MAIN, FS, WIN_P, WIN_S, MAX_LAG,       # noqa: E402
                         tt, rotar, cc_max)
import validar_egf_full as vf                                   # noqa: E402

TAG = "20200616T175441"
BANDS = [(0.10, 0.30), (0.15, 0.40)]
CC_MIN = 0.5


def evaluar_banda(ev, f1, f2):
    """CC P/S por estacion con la banda (f1, f2). Devuelve filas + figura."""
    ve.F1, ve.F2 = f1, f2          # preparar() lee estos globales del modulo
    filas, figs = [], []
    for sta in vf.STATIONS:
        mn, coords = vf.cargar_main_txt(sta)
        egf = ve.cargar_egf(TAG, ev, sta)
        if mn is None or egf is None:
            continue
        la, lo = coords
        d_m, az, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], la, lo)
        tp_m, ts_m = tt(kilometers2degrees(d_m / 1000.0), MAIN["depth_km"])
        d2, az2, _ = gps2dist_azimuth(ev["lat"], ev["lon"], la, lo)
        tp_e, ts_e = tt(kilometers2degrees(d2 / 1000.0), ev["depth_km"])

        def win(data, t_on, dur):
            i0 = int((t_on + 60.0) * FS)
            return data[max(0, i0):i0 + int(dur * FS)]

        rn_m, tt_m = rotar(mn["N"].data, mn["E"].data, az)
        rn_e, tt_e = rotar(egf["N"].data, egf["E"].data, az2)
        ccP, lagP = cc_max(win(mn["Z"].data, tp_m, WIN_P),
                           win(egf["Z"].data, tp_e, WIN_P), FS, MAX_LAG)
        ccS, lagS = cc_max(win(tt_m, ts_m - 5.0, WIN_S),
                           win(tt_e, ts_e - 5.0, WIN_S), FS, MAX_LAG)
        filas.append({"f1": f1, "f2": f2, "sta": sta,
                      "ccP_Z": round(ccP, 2), "lagP_s": round(lagP, 1),
                      "ccS_T": round(ccS, 2), "lagS_s": round(lagS, 1),
                      "valida": int(ccS >= CC_MIN)})
        figs.append((sta, mn["Z"].data, egf["Z"].data, tp_m, tp_e,
                     tt_m, tt_e, ts_m, ts_e, ccP, ccS))
        print(f"  {sta}: CC_P={ccP:.2f} ({lagP:+.1f}s)  CC_S={ccS:.2f} "
              f"({lagS:+.1f}s)  {'OK' if ccS >= CC_MIN else 'descartada'}")
    return filas, figs


def figura(figs, f1, f2):
    t = np.arange(len(figs[0][1])) / FS - 60.0
    fig, axes = plt.subplots(len(figs), 2, figsize=(13, 2.2 * len(figs)))
    axes = np.atleast_2d(axes)
    for i, (sta, zm, ze, tpm, tpe, tm, te_, tsm, tse, ccP, ccS) in \
            enumerate(figs):
        for j, (m_, e_, on_m, on_e, lbl, cc) in enumerate(
                [(zm, ze, tpm, tpe, "Z (P)", ccP),
                 (tm, te_, tsm, tse, "T (S)", ccS)]):
            ax = axes[i, j]
            ax.plot(t, m_ / max(np.max(np.abs(m_)), 1e-20), "k", lw=0.8)
            ax.plot(t + (on_m - on_e), e_ / max(np.max(np.abs(e_)), 1e-20),
                    "r", lw=0.8, alpha=0.8)
            ax.axvline(on_m, color="b", ls=":", lw=0.7)
            ax.text(0.98, 0.93, f"CC={cc:.2f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8,
                    color="forestgreen" if cc >= CC_MIN else "crimson",
                    fontweight="bold")
            ax.set_ylabel(f"{sta}\n{lbl}", fontsize=8)
            ax.set_xlim(on_m - 20, on_m + 80)
            ax.tick_params(labelsize=7)
    axes[-1, 0].set_xlabel("t desde origen (s)")
    axes[-1, 1].set_xlabel("t (s)")
    fig.suptitle(f"Validacion EGF {TAG} en {f1}-{f2} Hz "
                 "(negro=mainshock, rojo=EGF alineada)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fp = os.path.join(BASE, f"validacion_cc_{f1:.2f}-{f2:.2f}.png")
    fig.savefig(fp, dpi=140)
    plt.close(fig)
    print(f"  figura -> {fp}")


def main():
    with open(os.path.join(BASE, "candidatos_ranking.csv")) as f:
        rows = list(csv.DictReader(f))
    row = next(r for r in rows
               if UTCDateTime(r["time"]).strftime("%Y%m%dT%H%M%S") == TAG)
    ev = {"time": UTCDateTime(row["time"]), "lat": float(row["lat"]),
          "lon": float(row["lon"]), "depth_km": float(row["depth_km"])}

    todas = []
    for f1, f2 in BANDS:
        print(f"\n=== banda {f1}-{f2} Hz ===")
        filas, figs = evaluar_banda(ev, f1, f2)
        todas += filas
        if (f1, f2) != (0.10, 0.30):           # la 0.1-0.3 ya tiene figura
            figura(figs, f1, f2)
        ccs = [r["ccS_T"] for r in filas]
        print(f"  CC_S mediana = {np.median(ccs):.2f}  "
              f"({sum(r['valida'] for r in filas)}/{len(filas)} validas)")

    out = os.path.join(BASE, "cc_bandas.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(todas[0].keys()))
        w.writeheader(); w.writerows(todas)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
