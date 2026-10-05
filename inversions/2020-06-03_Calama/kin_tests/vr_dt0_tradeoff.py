"""Mapa del misfit en el plano (Vr, dt0), con el resto del modelo fijo.

El NA encuentra misfits casi iguales con Vr sub-Rayleigh (~4.3 km/s) y
supershear (~6.6 km/s), y cada familia viene con su propio dt0 (1.37 y 1.78 s).
Eso huele a intercambio: una ruptura más rápida llega antes y se compensa
atrasando la hora de origen. Este script lo mide directamente, barriendo solo
esos dos parámetros con los otros seis congelados en el mejor modelo de un caso.

Si sale un valle alargado y plano en diagonal, Vr no está resuelto y no se
puede afirmar nada sobre supershear desde el cinemático.

Uso:  KIN_CASE=np1_isola_potin_b0.04-0.15_n20_wide python vr_dt0_tradeoff.py
"""

import os
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kdellipspy as kde

from run_kin_multiseed import CASE, DT0_RANGE, KinematicDt0NA

VR = np.linspace(1.0, 7.5, 66)     # km/s, cubre sub-Rayleigh y supershear
DT0 = np.linspace(-1.0, 5.0, 61)   # s
BETA = 4.770                       # km/s en el hipocentro (Potin, capa de 112.5 km)


def best_model(case_dir):
    """Mejor modelo (8 parámetros) entre las semillas de multiseed."""
    d = joblib.load(case_dir / "multiseed" / "multiseed.joblib")
    return np.asarray(d["models"])[int(np.argmin(d["misfits"]))], d["param_names"]


def main():
    cfg = kde.ConfigParser(filepath=str(CASE / "input.ctl"))
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE / "input.ctl"), data_dir=str(CASE / "DATA"))
    inv = KinematicDt0NA(config=cfg, observed_waveforms=observed, time_array=time_array,
                         azi_times_array=kde.build_azi_times_array(config=cfg, model_name="iasp91"))
    inv.use_green_cache = True
    inv.param_ranges = np.vstack([inv.param_ranges, DT0_RANGE])
    inv.param_names = list(inv.param_names) + ["dt0 (s)"]

    m0, names = best_model(CASE)
    print(f"modelo base ({CASE.name}): " + "  ".join(f"{n.split()[0]}={v:.3f}" for n, v in zip(names, m0)),
          flush=True)

    mis = np.full((DT0.size, VR.size), np.nan)
    for j, vr in enumerate(VR):
        for i, dt0 in enumerate(DT0):
            m = m0.copy()
            m[6], m[7] = vr, dt0
            mis[i, j] = inv._evaluate_model(m)[0]
        print(f"  Vr {vr:5.2f} km/s ({vr / BETA:4.2f} beta)  misfit min {np.nanmin(mis[:, j]):.4f}", flush=True)
    inv.clear_green_cache()

    out = CASE / "multiseed"
    np.savez(out / "vr_dt0_tradeoff.npz", vr=VR, dt0=DT0, misfit=mis, model=m0, beta=BETA)

    best = np.unravel_index(np.nanargmin(mis), mis.shape)
    mmin = mis[best]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    lv = mmin * np.array([1.0, 1.02, 1.05, 1.10, 1.20, 1.40, 1.80, 2.5])
    cf = ax.contourf(VR, DT0, mis, levels=np.linspace(mmin, min(np.nanmax(mis), mmin * 2.5), 40), cmap="viridis_r")
    cs = ax.contour(VR, DT0, mis, levels=lv[1:], colors="w", linewidths=0.7)
    ax.clabel(cs, fmt="%.3f", fontsize=7)
    fig.colorbar(cf, label="L2 misfit")
    ax.plot(VR[best[1]], DT0[best[0]], "r*", ms=14, label=f"global min {mmin:.4f}")
    ax.axvline(BETA, color="k", ls="--", lw=1)
    ax.axvline(0.92 * BETA, color="k", ls=":", lw=1)
    ax.axvline(np.sqrt(2) * BETA, color="k", ls="--", lw=1)
    ax.text(0.92 * BETA, DT0[-1], " 0.92$\\beta$ (Rayleigh)", rotation=90, va="top", fontsize=7)
    ax.text(BETA, DT0[-1], r" $\beta$", rotation=90, va="top", fontsize=7)
    ax.text(np.sqrt(2) * BETA, DT0[-1], r" $\sqrt{2}\beta$", rotation=90, va="top", fontsize=7)
    ax.axvspan(BETA, np.sqrt(2) * BETA, color="0.5", alpha=0.25)
    ax.set_xlabel(r"Rupture velocity $V_r$ (km s$^{-1}$)")
    ax.set_ylabel(r"Origin-time correction $\delta t_0$ (s)")
    ax.set_title(f"Calama 2020 · {CASE.name}\n"
                 r"misfit with the other six parameters fixed at the best model")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    png = out / "vr_dt0_tradeoff.png"
    fig.savefig(png, dpi=140)
    print(f"\nminimo global: Vr {VR[best[1]]:.2f} km/s ({VR[best[1]] / BETA:.2f} beta), "
          f"dt0 {DT0[best[0]]:+.2f} s, misfit {mmin:.4f}", flush=True)
    # ¿cuánto se degrada el ajuste al imponer sub-Rayleigh?
    sub = VR <= 0.92 * BETA
    msub = np.nanmin(mis[:, sub])
    print(f"mejor sub-Rayleigh (Vr <= 0.92 beta = {0.92 * BETA:.2f}): misfit {msub:.4f} "
          f"({100 * (msub - mmin) / mmin:+.1f}% respecto del mínimo global)", flush=True)
    print(f"figura: {png}", flush=True)


if __name__ == "__main__":
    main()
