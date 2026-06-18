#!/usr/bin/env python3
"""
plot_trazas_usadas.py
=====================
Inspecciona las trazas que REALMENTE entran al misfit de la inversion, para
decidir que estaciones conservar. Carga `DATA/real_disp_{x,y,z}` (mismo loader
que run_inversion.py -> exactamente lo que se invierte), rota a R/T/Z por
estacion con el azimut, y dibuja una grilla:

    filas   = estaciones (ordenadas por azimut)
    columnas= R / T / Z

Sombrea la VENTANA usada en el misfit (igual que misfit_breakdown.py):
    R, Z -> ventana P  [tP, tP+TWIN]
    T    -> ventana S  [tS, tS+TWIN]

Y calcula un SNR por traza (RMS en la ventana de senal / RMS en ruido pre-P),
imprimiendo una tabla para identificar trazas de baja calidad.

NO corre ningun forward ni toca la inversion en curso (solo lee real_disp).

Uso:
    python codigos/plot_trazas_usadas.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
DATA_DIR = BASE / "DATA"
OUT = BASE / "output" / "trazas_usadas.png"

# Margen de ruido pre-P para el SNR (segundos antes de tP).
NOISE_PRE_S = 8.0


def rms(x):
    return float(np.sqrt(np.mean(x**2))) if x.size else 0.0


def main():
    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    f1, f2 = cfg.ellipse.freq1, cfg.ellipse.freq2
    twin = float(cfg.inversion_process.misfit_time_window)

    observed, time = kde.load_and_filter_observed_data(
        input_ctl_path=str(INPUT_CTL), data_dir=str(DATA_DIR),
        prefer_raw=False, freq1=f1, freq2=f2,
    )
    azt = np.asarray(kde.build_azi_times_array(config=cfg, model_name="iasp91"), float)
    azi, tP, tS = azt[:, 0], azt[:, 1], azt[:, 2]

    stas = cfg.stations.stations
    names = [s.name for s in stas]
    flags = np.array([[s.use_n, s.use_e, s.use_z] for s in stas], bool)
    nsta = len(stas)

    dt = float(time[1] - time[0])
    samp = max(1, int(np.rint(1.0 / dt)))
    win = max(1, int(np.rint(twin * samp)))
    npts = observed.shape[-1]

    # Orden por azimut (refleja la cobertura).
    order = list(np.argsort(azi))

    # --- Figura: nsta filas x 3 columnas (R, T, Z) --------------------------- #
    comps = ["R (vent. P)", "T (vent. S)", "Z (vent. P)"]
    fig, axes = plt.subplots(nsta, 3, figsize=(13, 2.0 * nsta), sharex=True)
    if nsta == 1:
        axes = axes[None, :]

    snr_table = []
    for row, j in enumerate(order):
        az = float(azi[j])
        N, E, Z = observed[j, 0], observed[j, 1], observed[j, 2]
        R = N * np.cos(az) + E * np.sin(az)
        T = E * np.cos(az) - N * np.sin(az)
        traces = [R, T, Z]

        # Indices de ventana (igual que el misfit).
        sp = int(np.rint(tP[j])) * samp
        ss = int(np.rint(tS[j])) * samp
        kp0, kp1 = max(0, sp - 1), min(npts - 1, sp + win - 1)       # P -> R, Z
        ks0, ks1 = max(0, ss - 1), min(npts - 1, ss + win - 1)       # S -> T
        win_idx = [(kp0, kp1), (ks0, ks1), (kp0, kp1)]
        # Ruido pre-P (comun a las 3 componentes).
        n0, n1 = max(0, sp - int(NOISE_PRE_S * samp)), max(1, sp - 1)

        snr_row = [names[j], f"{np.degrees(az):.0f}"]
        for c in range(3):
            ax = axes[row, c]
            tr = traces[c]
            k0, k1 = win_idx[c]
            used = bool(flags[j, c])
            color = "k" if used else "0.6"
            ax.plot(time, tr, color=color, lw=0.6)
            # sombrear ventana de senal
            ax.axvspan(time[k0], time[k1], color="tab:orange", alpha=0.20)
            # marca de llegada
            tarr = tP[j] if c in (0, 2) else tS[j]
            ax.axvline(tarr, color="tab:red", ls="--", lw=0.7)
            # SNR
            sig = tr[k0:k1 + 1]; noi = tr[n0:n1 + 1]
            snr = rms(sig) / rms(noi) if rms(noi) > 0 else np.inf
            snr_row.append(f"{snr:5.1f}" + ("" if used else "·off"))
            ax.set_yticks([])
            txt = f"{names[j]}  az={np.degrees(az):.0f}°  SNR={snr:.1f}"
            ax.text(0.01, 0.92, txt, transform=ax.transAxes, fontsize=7,
                    va="top", color=color,
                    bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))
            if not used:
                ax.text(0.5, 0.5, "NO USADA", transform=ax.transAxes,
                        ha="center", va="center", fontsize=11, color="0.6", alpha=0.5)
            if row == 0:
                ax.set_title(comps[c], fontsize=10)
        snr_table.append(snr_row)

    for c in range(3):
        axes[-1, c].set_xlabel("Tiempo desde el origen [s]")
    fig.suptitle(f"Trazas usadas (real_disp R/T/Z, {f1:.3f}-{f2:.3f} Hz, "
                 f"vent. {twin:.0f}s) — naranja=ventana misfit, rojo=llegada",
                 fontsize=11, y=1.002)
    fig.tight_layout(rect=[0, 0, 1, 0.995])
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"figura -> {OUT}")

    # --- Tabla de SNR por estacion/componente ------------------------------- #
    print(f"\nSNR (RMS senal en ventana / RMS ruido pre-P {NOISE_PRE_S:.0f}s):")
    print(f"  {'Est':<6}{'az':>5}{'SNR_R':>9}{'SNR_T':>9}{'SNR_Z':>9}")
    print("  " + "-" * 38)
    for r in snr_table:
        print(f"  {r[0]:<6}{r[1]:>5}{r[2]:>9}{r[3]:>9}{r[4]:>9}")
    print("\n  (SNR bajo ~<2-3 => traza ruidosa, candidata a descartar)")


if __name__ == "__main__":
    main()
