#!/usr/bin/env python3
"""
plot_fit_rot_win.py
===================
Plot del ajuste de formas de onda CONSISTENTE con el misfit: rota las trazas a
R/T/Z por estacion (usando el azimut) y SOMBREA la ventana usada en el misfit:
    columna R (radial)     -> ventana P  [tP, tP+TWIN]
    columna T (transversal)-> ventana S  [tS, tS+TWIN]
    columna Z (vertical)   -> ventana P  [tP, tP+TWIN]

(El plot estandar plot_waveform_fit dibuja N/E/Z crudas y sin ventanas, lo que
NO refleja como se calcula el misfit; este script lo corrige.)

Lee output/inversion_result.joblib (observed, best_synthetics, time, config).

Uso:
    python codigos/plot_fit_rot_win.py
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
OUT = BASE / "output" / "waveform_fit_rotwin.png"
TWIN = 20.0   # ventana del misfit (s)


def main():
    r = joblib.load(JOBLIB)
    obs = np.asarray(r.observed, dtype=float)        # (nsta,3,npts) N,E,Z
    syn = np.asarray(r.best_synthetics, dtype=float)
    time = np.asarray(r.time, dtype=float)
    cfg = r.config
    names = [s.name for s in cfg.stations.stations]
    flags = np.array([[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations], bool)
    azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")  # [azi,tP,tS]
    azi, tP, tS = azt[:, 0], azt[:, 1], azt[:, 2]
    nsta = obs.shape[0]

    # Rotar a R/T/Z.
    R_o = np.zeros_like(obs[:, 0]); R_s = np.zeros_like(obs[:, 0])
    T_o = np.zeros_like(obs[:, 0]); T_s = np.zeros_like(obs[:, 0])
    for j in range(nsta):
        c, s = np.cos(azi[j]), np.sin(azi[j])
        R_o[j] = obs[j, 0]*c + obs[j, 1]*s
        R_s[j] = syn[j, 0]*c + syn[j, 1]*s
        T_o[j] = obs[j, 1]*c - obs[j, 0]*s
        T_s[j] = syn[j, 1]*c - syn[j, 0]*s
    rot_o = [R_o, T_o, obs[:, 2]]   # R, T, Z
    rot_s = [R_s, T_s, syn[:, 2]]
    comp_lbl = ["R (P)", "T (S)", "Z (P)"]
    # Ventana por componente: R->P(tP), T->S(tS), Z->P(tP)
    win_start = [tP, tS, tP]

    fig, axs = plt.subplots(nsta, 3, figsize=(11, 1.6*nsta), squeeze=False, sharex=True)
    fig.suptitle(f"Ajuste rotado R/T/Z + ventana {TWIN:.0f}s "
                 f"(P->R,Z ; S->T)   misfit={r.best_model.misfit:.4f}",
                 fontsize=13, y=1.005)

    for j in range(nsta):
        for c in range(3):
            ax = axs[j, c]
            ax.plot(time, rot_o[c][j], color="black", lw=1.4,
                    label="Obs" if j == 0 and c == 0 else None)
            ax.plot(time, rot_s[c][j], color="red", lw=1.1, ls="--",
                    label="Syn" if j == 0 and c == 0 else None)
            # Ventana sombreada [t0, t0+TWIN].
            t0 = float(win_start[c][j])
            ax.axvspan(t0, t0 + TWIN, color="orange", alpha=0.18, lw=0)
            ax.axvline(t0, color="orange", lw=0.8, alpha=0.7)
            # Marca si la componente esta desactivada por flags.
            comp_flag = flags[j, {0: 0, 1: 1, 2: 2}[c]]
            if not comp_flag:
                ax.text(0.5, 0.5, "off", transform=ax.transAxes, ha="center",
                        color="gray", fontsize=8)
            if j == 0:
                ax.set_title(comp_lbl[c])
            if c == 0:
                ax.set_ylabel(names[j], fontweight="bold")
            if j == nsta - 1:
                ax.set_xlabel("Tiempo (s)")
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelleft=False, left=False, right=False, top=False)

    fig.legend(loc="upper right", bbox_to_anchor=(1.02, 1.0))
    plt.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"figura -> {OUT}")
    print(f"  ventana P (naranja) en R y Z desde tP; ventana S en T desde tS; TWIN={TWIN:.0f}s")


if __name__ == "__main__":
    main()
