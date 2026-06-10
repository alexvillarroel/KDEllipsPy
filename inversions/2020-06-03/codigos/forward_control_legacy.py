#!/usr/bin/env python3
"""
forward_control_legacy.py
=========================
Forward de CONTROL: corre KDEllipsPy con el MODELO del legacy y sus MISMAS
convenciones de forward (STF tipo 5 = rampa causal, aw=2, desplazamiento, banda
0.02-0.1 causal) y lo SUPERPONE con el sintetico del legacy
(Event/kine_files/best_disp_*) para ver si calzan.

Si calzan -> las diferencias previas eran STF (5 vs 4) + aw (2 vs 0.5), no las
Green en si.

Uso:
    python codigos/forward_control_legacy.py
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
KF = BASE / "legacy" / "Kinematic_inversion" / "Event" / "kine_files"
OUT = BASE / "output" / "forward_control_legacy.png"

LEGACY_MODEL = np.array([8.06, 4.81, 0.31, 1.00, 0.81, 3.64, 1.26])  # a1,a2,theta,np,tp,dmax,vr
AW = 0.5
SRC_TYPE = 4  # prueba con la STF default de KDEllipsPy


def load_flat(prefix):
    """Carga best/real _disp_{x,y,z} -> (nsta,3,npts)."""
    comps = []
    for c in "xyz":
        v = np.loadtxt(KF / f"{prefix}_{c}")
        comps.append(v)
    npts = 512
    nsta = comps[0].size // npts
    arr = np.stack([c.reshape(nsta, npts) for c in comps], axis=1)  # (nsta,3,npts)
    return arr


def main():
    synth_leg = load_flat("best_disp")     # sintetico legacy (desplazamiento)
    obs_leg = load_flat("real_disp")       # observado legacy
    nsta, _, npts = obs_leg.shape
    print(f"legacy synth {synth_leg.shape}, obs {obs_leg.shape}")

    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    cfg.ellipse.source_type = SRC_TYPE     # rampa causal
    cfg.observed_data.units = 1            # desplazamiento
    dt = float(cfg.observed_data.delta)
    time = np.arange(npts) * dt
    azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")

    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs_leg,
                                 time_array=time, azi_times_array=azt)
    model.use_green_cache = False
    model.axitra_aw = AW
    print(f"Forward KDEllipsPy: modelo legacy, ics={SRC_TYPE}, aw={AW}, disp, 0.02-0.1 causal...")
    misfit, synth_kde = model._evaluate_model(LEGACY_MODEL)
    print(f"misfit (vs obs legacy) = {misfit:.4f}")

    # Correlacion sintetico-KDE vs sintetico-legacy por traza.
    names = [s.name for s in cfg.stations.stations]
    corrs = []
    for j in range(nsta):
        for c in range(3):
            a, b = synth_kde[j, c], synth_leg[j, c]
            if a.std() > 0 and b.std() > 0:
                corrs.append(np.corrcoef(a, b)[0, 1])
    print(f"correlacion KDE-synth vs legacy-synth: media {np.mean(corrs):.3f}, "
          f"min {np.min(corrs):.3f}")

    # Figura: obs (negro) / legacy-synth (azul) / KDE-synth (rojo)
    comp_lbl = ["N", "E", "Z"]
    fig, axs = plt.subplots(nsta, 3, figsize=(12, 1.5 * nsta), squeeze=False, sharex=True)
    fig.suptitle(f"Forward de control — modelo legacy · ics=5 · aw=2 · disp\n"
                 f"(negro=obs, azul=synth LEGACY, rojo=synth KDEllipsPy)  "
                 f"corr media {np.mean(corrs):.2f}", fontsize=12, y=1.005)
    for j in range(nsta):
        for c in range(3):
            ax = axs[j, c]
            ax.plot(time, obs_leg[j, c], color="black", lw=1.0, alpha=0.6)
            ax.plot(time, synth_leg[j, c], color="tab:blue", lw=1.6)
            ax.plot(time, synth_kde[j, c], color="red", lw=1.1, ls="--")
            if j == 0:
                ax.set_title(comp_lbl[c])
            if c == 0:
                ax.set_ylabel(names[j], fontweight="bold")
            if j == nsta - 1:
                ax.set_xlabel("Tiempo (s)")
            ax.grid(alpha=0.3); ax.tick_params(labelleft=False)
    plt.tight_layout()
    fig.savefig(str(OUT), dpi=180, bbox_inches="tight")
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
