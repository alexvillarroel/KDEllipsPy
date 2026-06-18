#!/usr/bin/env python3
"""
forward_units_check.py
======================
Verifica que el flag de unidades realmente cambia el forward: corre el MISMO
modelo (mejor solucion) con unit=2 (velocidad) y unit=1 (desplazamiento), sin
filtrar, y comprueba que  desplazamiento ~ integral temporal de la velocidad
(es decir, velocidad ~ d/dt desplazamiento). Si se cumple, el units funciona.
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "output" / "units_check.png"

# Mejor solucion (corrida de velocidad previa)
BEST = np.array([9.8454, 4.5649, 0.4062, 0.6339, 0.9116, 3.4557, 1.5239])


def to_nsta3npts(sx, sy, sz):
    a = np.array([sx, sy, sz])
    a = np.transpose(a, (1, 2, 0))
    a = np.transpose(a, (0, 2, 1))
    return a  # (nsta,3,npts)


def main():
    cfg = kde.ConfigParser(filepath=str(BASE / "input.ctl"))
    dt = float(cfg.observed_data.delta)
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.build_geometry_with_ellipse_slip(BEST)
    ap = fm.build_axitra(geom, latlon=False, freesurface=True, aw=0.5)
    ap = fm.green(ap, quiet=True)

    _, vx, vy, vz = fm.conv(ap, geom, source_type=4, t0=3.0, unit=2, quiet=True)  # velocidad
    _, dx, dy, dz = fm.conv(ap, geom, source_type=4, t0=3.0, unit=1, quiet=True)  # desplazamiento
    vel = to_nsta3npts(vx, vy, vz)
    dis = to_nsta3npts(dx, dy, dz)
    npts = vel.shape[2]
    t = np.arange(npts) * dt

    print(f"shapes vel {vel.shape} disp {dis.shape}")
    print(f"max |vel|  = {np.abs(vel).max():.3e}")
    print(f"max |disp| = {np.abs(dis).max():.3e}")

    # disp deberia ~ integral de vel (quitando media para la constante)
    intvel = np.zeros_like(vel)
    intvel[:, :, 1:] = cumulative_trapezoid(vel, dx=dt, axis=2)
    # correlacion disp vs integral(vel)
    corr = []
    for j in range(vel.shape[0]):
        for c in range(3):
            a = dis[j, c] - dis[j, c].mean()
            b = intvel[j, c] - intvel[j, c].mean()
            if a.std() > 0 and b.std() > 0:
                corr.append(np.corrcoef(a, b)[0, 1])
    print(f"corr( disp , ∫vel ) media = {np.mean(corr):.4f}  (≈1 => units OK)")

    # Figura: una estacion-componente representativa (la de mayor energia en disp)
    e = np.array([[np.sum(dis[j, c]**2) for c in range(3)] for j in range(dis.shape[0])])
    j, c = np.unravel_index(np.argmax(e), e.shape)
    names = [s.name for s in cfg.stations.stations]
    comp = ["N", "E", "Z"][c]
    fig, axs = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    axs[0].plot(t, vel[j, c], "b"); axs[0].set_title(f"{names[j]} {comp} — velocidad (unit=2)")
    axs[1].plot(t, dis[j, c], "r"); axs[1].set_title("desplazamiento (unit=1)")
    axs[2].plot(t, dis[j, c], "r", label="disp (unit=1)")
    ig = intvel[j, c] - intvel[j, c].mean() + dis[j, c].mean()
    axs[2].plot(t, ig, "k--", label="∫ vel dt")
    axs[2].legend(); axs[2].set_title("disp vs integral de la velocidad")
    axs[2].set_xlabel("Tiempo (s)")
    for ax in axs:
        ax.grid(alpha=0.3)
    plt.tight_layout(); fig.savefig(str(OUT), dpi=160)
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
