#!/usr/bin/env python3
"""
test_source_type.py
===================
Compara source type 4 (triangulo, actual) vs 5 (ramp, como el legacy de Calama)
con el MISMO mejor modelo y rise time t0=3 s, SIN re-invertir. Calcula el misfit
con cada uno y grafica la componente T (la que domina) obs vs sint-4 vs sint-5
para unas estaciones, para ver cuanto cambia a la banda 0.06-0.15 Hz.

Uso:
    python codigos/test_source_type.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import joblib

import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "output" / "test_source_type.png"
T0 = 3.0


def main():
    r = joblib.load(BASE / "output" / "inversion_result.joblib")
    best7 = np.asarray(r.best_model.model, float)
    cfg = r.config

    obs, t = kde.load_and_filter_observed_data(
        input_ctl_path=str(BASE / "input.ctl"), data_dir=str(BASE / "DATA"),
        prefer_raw=False, freq1=cfg.ellipse.freq1, freq2=cfg.ellipse.freq2,
    )
    azt = np.asarray(kde.build_azi_times_array(config=cfg, model_name="iasp91"), float)
    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs,
                                 time_array=t, azi_times_array=azt)
    model.use_green_cache = True
    cfg.ellipse.t0 = T0

    res = {}
    for st in (4, 5):
        cfg.ellipse.source_type = st
        mf, syn = model._evaluate_model(best7)
        res[st] = (mf, np.asarray(syn, float))
        print(f"source_type={st} (t0={T0}s): misfit = {mf:.6f}")

    mf4, syn4 = res[4]
    mf5, syn5 = res[5]
    dpc = 100 * (mf5 - mf4) / mf4
    print(f"\nΔmisfit (5 vs 4) = {mf5-mf4:+.6f}  ({dpc:+.2f}%)")
    rms_rel = np.sqrt(np.mean((syn5 - syn4)**2)) / (np.sqrt(np.mean(syn4**2)) or 1)
    print(f"diferencia RMS sinteticas (5 vs 4) = {100*rms_rel:.2f}% de la amplitud")

    # --- figura: componente T (rotada) obs vs sint4 vs sint5 por estacion ---
    azi = azt[:, 0]
    names = [s.name for s in cfg.stations.stations]
    flags = [(s.use_n, s.use_e, s.use_z) for s in cfg.stations.stations]
    used_T = [j for j in range(len(names)) if flags[j][1]]   # use_e -> T
    n = len(used_T)
    fig, ax = plt.subplots(n, 1, figsize=(10, 1.6 * n), sharex=True)
    if n == 1:
        ax = [ax]
    for row, j in enumerate(used_T):
        az = float(azi[j])
        def rotT(a): return a[j, 1] * np.cos(az) - a[j, 0] * np.sin(az)
        ax[row].plot(t, rotT(obs), "k", lw=0.8, label="obs")
        ax[row].plot(t, rotT(syn4), "tab:blue", lw=0.9, label="type 4 (triang)")
        ax[row].plot(t, rotT(syn5), "tab:red", lw=0.9, ls="--", label="type 5 (ramp)")
        ax[row].set_yticks([])
        ax[row].text(0.01, 0.85, f"{names[j]} (T)", transform=ax[row].transAxes, fontsize=8)
        if row == 0:
            ax[row].legend(fontsize=7, ncol=3, loc="upper right")
    ax[-1].set_xlabel("Tiempo [s]")
    fig.suptitle(f"Source type 4 vs 5 (t0={T0}s) — comp. T | "
                 f"misfit 4={mf4:.4f}  5={mf5:.4f} ({dpc:+.1f}%)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT, dpi=150)
    print(f"\nfigura -> {OUT}")


if __name__ == "__main__":
    main()
