#!/usr/bin/env python3
"""
plot_fit_amplitude.py
=====================
Igual que plot_fit_rot_win pero SIN normalizar entre paneles: todas las
estaciones comparten el mismo eje Y, asi se compara visualmente la amplitud
real (estaciones cercanas grandes, lejanas chicas). Anota el pico |obs| por
panel. Rota a R/T/Z para ser consistente con el misfit.

Lee output/inversion_result.joblib. No toca all_plots.pdf ni el dashboard.

Uso:
    python codigos/plot_fit_amplitude.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import joblib
import numpy as np
import matplotlib.pyplot as plt

import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
JOBLIB = BASE / "output" / "inversion_result.joblib"
OUT = BASE / "output" / "waveform_fit_amplitude.png"


def build_fig(obs, syn, time, names, azi):
    nsta = obs.shape[0]
    R_o, R_s = np.zeros_like(obs[:, 0]), np.zeros_like(obs[:, 0])
    T_o, T_s = np.zeros_like(obs[:, 0]), np.zeros_like(obs[:, 0])
    for j in range(nsta):
        c, s = np.cos(azi[j]), np.sin(azi[j])
        R_o[j] = obs[j, 0] * c + obs[j, 1] * s
        R_s[j] = syn[j, 0] * c + syn[j, 1] * s
        T_o[j] = obs[j, 1] * c - obs[j, 0] * s
        T_s[j] = syn[j, 1] * c - syn[j, 0] * s
    rot_o = [R_o, T_o, obs[:, 2]]
    rot_s = [R_s, T_s, syn[:, 2]]
    lbl = ["R (P)", "T (S)", "Z (P)"]

    # sharey=True -> misma escala de amplitud en TODOS los paneles (sin normalizar)
    fig, axs = plt.subplots(
        nsta, 3, figsize=(11, 1.6 * nsta), squeeze=False, sharex=True, sharey=True
    )
    for j in range(nsta):
        for c in range(3):
            ax = axs[j, c]
            ax.plot(
                time,
                rot_o[c][j],
                color="black",
                lw=1.4,
                label="Obs" if j == 0 and c == 0 else None,
            )
            ax.plot(
                time,
                rot_s[c][j],
                color="red",
                lw=1.1,
                ls="--",
                label="Syn" if j == 0 and c == 0 else None,
            )
            pk = float(np.max(np.abs(rot_o[c][j])))
            ax.text(
                0.97,
                0.9,
                f"{pk:.2e}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7,
                color="dimgray",
            )
            if j == 0:
                ax.set_title(lbl[c])
            if c == 0:
                ax.set_ylabel(names[j], fontweight="bold")
            if j == nsta - 1:
                ax.set_xlabel("Tiempo (s)")
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelleft=(c == 0), left=(c == 0), right=False, top=False)
    fig.legend(loc="upper right", bbox_to_anchor=(1.02, 1.0))
    fig.text(
        0.005,
        0.5,
        "Desplazamiento (m) — escala comun",
        rotation=90,
        va="center",
        fontsize=9,
    )
    return fig, axs


def main():
    r = joblib.load(JOBLIB)
    obs = np.asarray(r.observed, dtype=float)
    syn = np.asarray(r.best_synthetics, dtype=float)
    time = np.asarray(r.time, dtype=float)
    cfg = r.config
    names = [s.name for s in cfg.stations.stations]
    azi = kde.build_azi_times_array(config=cfg, model_name="iasp91")[:, 0]

    fig, _ = build_fig(obs, syn, time, names, azi)
    fig.suptitle(
        f"Ajuste R/T/Z — amplitudes reales (eje Y comun)   "
        f"misfit={r.best_model.misfit:.4f}",
        fontsize=13,
        y=1.005,
    )
    plt.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"figura -> {OUT}")


def _check():
    # dos estaciones, una 10x mas grande: los paneles deben compartir ylim
    # (no normalizado) -> la chica se ve chica.
    t = np.linspace(0, 50, 200)
    obs = np.stack(
        [
            np.array([np.sin(t), np.sin(t), np.sin(t)]),
            10 * np.array([np.sin(t), np.sin(t), np.sin(t)]),
        ]
    )
    fig, axs = build_fig(obs, obs * 0.9, t, ["S1", "S2"], np.array([0.0, 0.0]))
    yl = {ax.get_ylim() for ax in axs.ravel()}
    assert len(yl) == 1, "eje Y no compartido -> estaria normalizando"
    print("check OK: eje Y compartido, amplitudes sin normalizar")


if __name__ == "__main__":
    import sys

    _check() if "--check" in sys.argv else main()
