#!/usr/bin/env python3
"""
mapa_azimutal.py
================
Diagrama polar (radar / rosa) de la cobertura azimutal de las estaciones usadas
en la inversion, respecto al epicentro. Para cada estacion dibuja un rayo desde
el centro con:
    angulo = azimut evento->estacion (grados, N arriba, sentido horario)
    radio  = distancia epicentral (km)
Sombrea el MAYOR hueco azimutal (azimuthal gap), metrica clave de la calidad de
la cobertura.

Lee event.ctl (epicentro) y input.ctl (estaciones de la inversion).

Uso:
    python codigos/mapa_azimutal.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth

import kdellipspy as kde
from event import load_event

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
CTL_FILE = BASE / "event.ctl"
OUT = BASE / "output" / "mapa_azimutal.png"


def plot_azimuthal_coverage(cfg, ev_lat, ev_lon, save_path=None, show=False):
    """Diagrama polar de cobertura azimutal de las estaciones de ``cfg`` respecto
    al epicentro (ev_lat, ev_lon). Devuelve (fig, gap_max_deg)."""
    az_deg, dist_km, names = [], [], []
    for st in cfg.stations.stations:
        d_m, az, _baz = gps2dist_azimuth(ev_lat, ev_lon, st.latitude, st.longitude)
        az_deg.append(az)
        dist_km.append(d_m / 1000.0)
        names.append(st.name)
    az_deg = np.array(az_deg)
    dist_km = np.array(dist_km)

    order = np.argsort(az_deg)
    a_sorted = az_deg[order]
    gaps = np.diff(np.concatenate([a_sorted, [a_sorted[0] + 360.0]]))
    gmax_i = int(np.argmax(gaps))
    gap_max = float(gaps[gmax_i])
    gap_lo = a_sorted[gmax_i]
    gap_hi = gap_lo + gap_max

    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    rmax = dist_km.max() * 1.15
    ax.set_ylim(0, rmax)

    th = np.deg2rad(np.linspace(gap_lo, gap_hi, 50))
    ax.fill_between(th, 0, rmax, color="red", alpha=0.12, label=f"Gap máx = {gap_max:.0f}°")

    for az, d, nm in zip(az_deg, dist_km, names):
        t = np.deg2rad(az)
        ax.plot([t, t], [0, d], color="0.5", lw=1.0, zorder=2)
        ax.plot(t, d, marker="^", ms=12, color="dodgerblue", markeredgecolor="k", zorder=3)
        ax.text(t, d + rmax*0.05, f"{nm}\n{az:.0f}° · {d:.0f}km",
                ha="center", va="center", fontsize=8, zorder=4)

    ax.plot(0, 0, marker="*", ms=14, color="yellow", markeredgecolor="k", zorder=5)
    ax.set_rlabel_position(135)
    ax.set_title(f"Cobertura azimutal — {len(names)} estaciones\n"
                 f"gap máximo {gap_max:.0f}° ({gap_lo:.0f}°–{gap_hi:.0f}°)",
                 fontsize=12, pad=20)
    ax.legend(loc="lower right", bbox_to_anchor=(1.15, -0.05), fontsize=9)

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return fig, gap_max, list(zip(names, az_deg, dist_km))


def main():
    ev = load_event(str(CTL_FILE))
    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    _fig, gap_max, info = plot_azimuthal_coverage(cfg, ev["lat"], ev["lon"], save_path=OUT)
    print(f"mapa azimutal -> {OUT}")
    print("  azimuts: " + ", ".join(f"{n}={a:.0f}°" for n, a, _d in info))
    print(f"  gap azimutal máximo = {gap_max:.0f}°")


if __name__ == "__main__":
    main()
