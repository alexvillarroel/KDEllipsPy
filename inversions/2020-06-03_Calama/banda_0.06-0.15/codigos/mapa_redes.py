#!/usr/bin/env python3
"""
mapa_redes.py
=============
Compara las DOS redes de estaciones usadas para el terremoto de Calama 2020:
  * RECIENTE (kdellipspy, 9 est): A22F A08F A09F PB19 A05F PB06 A06F A18F A23F
  * LEGACY   (Fortran, 8 est):    A05C AC01 AC02 PB05 PB06 PB08 PB19 T15A
Comunes: PB19, PB06.

Genera:
  output/mapa_redes.png            -> mapa PyGMT (relieve+costa) con ambas redes,
                                       epicentro y mecanismo focal.
  output/cobertura_azimutal.png    -> dos diagramas polares (azimut vs distancia)
                                       con el GAP azimutal máximo sombreado.

Uso:
    python codigos/mapa_redes.py
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
import pygmt

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "output"

EV = dict(lat=-23.247, lon=-68.530, depth=123.4, strike=333, dip=60, rake=-91, mw=6.5)

RECIENTE = {  # 9 estaciones
    "A22F": (-22.909, -68.201), "A08F": (-23.186, -68.006), "A09F": (-23.826, -67.442),
    "PB19": (-23.905, -69.291), "A05F": (-22.891, -69.321), "PB06": (-22.706, -69.572),
    "A06F": (-22.344, -69.658), "A18F": (-22.457, -68.929), "A23F": (-22.337, -68.651)}
LEGACY = {   # 8 estaciones
    "A05C": (-27.3612, -70.3390), "AC01": (-26.1479, -70.5987), "AC02": (-26.8355, -69.1291),
    "PB05": (-22.8528, -70.2023), "PB06": (-22.7058, -69.5717), "PB08": (-20.1411, -69.1535),
    "PB19": (-23.9048, -69.2900), "T15A": (-20.2393, -70.0540)}
COMUNES = set(RECIENTE) & set(LEGACY)


def az_dist(net):
    out = {}
    for n, (la, lo) in net.items():
        d, az, _ = gps2dist_azimuth(EV["lat"], EV["lon"], la, lo)
        out[n] = (d / 1000.0, az)
    return out


def max_gap(azs):
    a = sorted(azs); n = len(a)
    gaps = [((a[(i + 1) % n] - a[i]) % 360) for i in range(n)]
    i = int(np.argmax(gaps))
    return gaps[i], a[i], a[i] + gaps[i]   # tamaño, inicio, fin del hueco


# ----------------------------------------------------------------------------- #
# 1) MAPA PyGMT con ambas redes
# ----------------------------------------------------------------------------- #
def mapa():
    REGION = [-71.2, -67.0, -28.2, -19.6]
    PROJ = "M14c"
    fig = pygmt.Figure()
    try:
        grid = pygmt.datasets.load_earth_relief(resolution="30s", region=REGION)
    except Exception:
        grid = pygmt.datasets.load_earth_relief(resolution="01m", region=REGION)
    pygmt.makecpt(cmap="gray", series=[float(grid.min()), float(grid.max())])
    fig.grdimage(grid=grid, region=REGION, projection=PROJ, cmap=True, shading=True)
    fig.coast(region=REGION, projection=PROJ, shorelines="0.5p,black",
              borders=["1/0.6p,black"], frame=["af", "+tRedes de estaciones - Calama 2020"])

    # lineas evento->estacion (tenues) para visualizar azimutes
    for net, pen in [(RECIENTE, "0.4p,red@60"), (LEGACY, "0.4p,blue@60")]:
        for n, (la, lo) in net.items():
            fig.plot(x=[EV["lon"], lo], y=[EV["lat"], la], pen=pen, region=REGION, projection=PROJ)

    # estaciones
    for net, style, fill, lab in [(LEGACY, "i0.42c", "dodgerblue", "Legacy (8)"),
                                  (RECIENTE, "t0.40c", "red", "Reciente (9)")]:
        xs = [lo for la, lo in net.values()]; ys = [la for la, lo in net.values()]
        fig.plot(x=xs, y=ys, style=style, fill=fill, pen="0.6p,black",
                 region=REGION, projection=PROJ, label=lab)
        for n, (la, lo) in net.items():
            fig.text(x=lo, y=la, text=n, font="6p,Helvetica-Bold,black",
                     offset="0/0.32c", fill="white@40", region=REGION, projection=PROJ)
    # comunes: anillo amarillo
    for n in COMUNES:
        la, lo = RECIENTE[n]
        fig.plot(x=[lo], y=[la], style="c0.55c", pen="1.6p,gold", region=REGION, projection=PROJ)

    # mecanismo focal + epicentro
    fig.meca(spec=dict(strike=EV["strike"], dip=EV["dip"], rake=EV["rake"], magnitude=EV["mw"]),
             scale="1.0c", longitude=EV["lon"], latitude=EV["lat"], depth=EV["depth"],
             compressionfill="black", region=REGION, projection=PROJ)
    fig.plot(x=[EV["lon"]], y=[EV["lat"]], style="a0.5c", fill="yellow", pen="1p,black",
             region=REGION, projection=PROJ)
    fig.legend(region=REGION, projection=PROJ, position="jTR+o0.2c", box="+gwhite+p0.5p")
    fig.savefig(str(OUT / "mapa_redes.png"), dpi=300)
    print("mapa -> output/mapa_redes.png")


# ----------------------------------------------------------------------------- #
# 2) Cobertura azimutal (polar) de ambas redes
# ----------------------------------------------------------------------------- #
def cobertura():
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), subplot_kw={"projection": "polar"})
    for ax, (net, color, titulo) in zip(axes, [(RECIENTE, "red", "Reciente (9 est)"),
                                               (LEGACY, "dodgerblue", "Legacy (8 est)")]):
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        ad = az_dist(net)
        azs = [a for _, a in ad.values()]
        gap, g0, g1 = max_gap(azs)
        rmax = max(d for d, _ in ad.values()) * 1.15
        # sombrear el hueco azimutal maximo
        th = np.radians(np.linspace(g0, g1, 60))
        ax.fill_between(th, 0, rmax, color="red", alpha=0.12)
        ax.plot(np.radians([g0, g0]), [0, rmax], "r--", lw=0.8)
        ax.plot(np.radians([g1, g1]), [0, rmax], "r--", lw=0.8)
        # estaciones
        for n, (d, az) in ad.items():
            mk = "*" if n in COMUNES else ("v" if net is RECIENTE else "^")
            ax.plot(np.radians(az), d, mk, ms=13 if n in COMUNES else 9,
                    color="gold" if n in COMUNES else color,
                    mec="black", mew=0.6, zorder=5)
            ax.text(np.radians(az), d, " " + n, fontsize=7, va="bottom")
        ax.set_rlim(0, rmax); ax.set_rlabel_position(112)
        ax.set_title(f"{titulo}\nGAP máx = {gap:.0f}°", fontsize=11, pad=18)
        ax.grid(alpha=0.4)
    fig.suptitle("Cobertura azimutal — Calama 2020 (estrella dorada = estación común)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "cobertura_azimutal.png", dpi=150, bbox_inches="tight")
    print("cobertura -> output/cobertura_azimutal.png")
    # tabla resumen
    for net, nm in [(RECIENTE, "Reciente"), (LEGACY, "Legacy")]:
        ad = az_dist(net); gap, _, _ = max_gap([a for _, a in ad.values()])
        dd = [d for d, _ in ad.values()]
        print(f"  {nm:9s}: {len(net)} est | gap {gap:5.1f}° | dist {min(dd):.0f}-{max(dd):.0f} km")


if __name__ == "__main__":
    cobertura()
    mapa()
