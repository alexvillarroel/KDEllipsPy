#!/usr/bin/env python3
"""
map.py
======
Mapa de TODAS las estaciones disponibles para el terremoto de Calama 2020:
  * Todas las de ``registros_20200603/`` (lee lat/lon del header de cada .txt).
  * + las de la red LEGACY que no estan en registros (AC01, AC02, PB08).

Categoriza:
  - gris  : disponibles, no usadas
  - rojo  : usadas en la inversion reciente (kdellipspy, 9)
  - azul  : red legacy (8)
  - anillo dorado : estaciones comunes a ambas redes

Relieve + costa, epicentro + mecanismo focal. La region se autoajusta a todas
las estaciones graficadas.

Uso:
    python codigos/map.py
"""
import glob
import os
import re

import numpy as np
import pygmt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REGISTROS = os.path.join(BASE, "registros_20200603")
OUT = os.path.join(BASE, "output", "mapa_todas_estaciones.png")

# Evento Calama 2020
EV = dict(lat=-23.247, lon=-68.530, depth=123.4, strike=333, dip=60, rake=-91, mw=6.5)

USED = {"A22F", "A08F", "A09F", "PB19", "A05F", "PB06", "A06F", "A18F", "A23F"}   # reciente (9)
LEGACY = {"A05C", "AC01", "AC02", "PB05", "PB06", "PB08", "PB19", "T15A"}         # legacy (8)
COMUNES = USED & LEGACY
# coords de las legacy que NO estan en registros_20200603
LEGACY_EXTRA = {"AC01": (-26.1479, -70.5987), "AC02": (-26.8355, -69.1291),
                "PB08": (-20.1411, -69.1535)}


def coords_registros():
    """Lee lat/lon del header (# Latitud: .. Longitud: ..) de un .txt por estacion."""
    out = {}
    for d in sorted(glob.glob(os.path.join(REGISTROS, "*/"))):
        sta = os.path.basename(d.rstrip("/"))
        if sta.lower() == "zips":
            continue
        txts = glob.glob(os.path.join(d, "*.txt"))
        if not txts:
            continue
        with open(txts[0]) as fh:
            for ln in fh:
                if "Latitud" in ln:
                    m = re.findall(r"-?\d+\.\d+", ln)
                    if len(m) >= 2:
                        out[sta] = (float(m[0]), float(m[1]))
                    break
                if not ln.startswith("#"):
                    break
    return out


def main():
    coords = coords_registros()
    coords.update(LEGACY_EXTRA)               # agrega legacy faltantes
    print(f"Estaciones totales en el mapa: {len(coords)}")

    # region autoajustada a todas las estaciones + epicentro
    lons = [lo for la, lo in coords.values()] + [EV["lon"]]
    lats = [la for la, lo in coords.values()] + [EV["lat"]]
    m = 0.4
    region = [min(lons) - m, max(lons) + m, min(lats) - m, max(lats) + m]
    proj = "M15c"
    print(f"Region: {[round(r,2) for r in region]}")

    fig = pygmt.Figure()
    try:
        grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    except Exception:
        grid = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
    pygmt.makecpt(cmap="gray", series=[float(grid.min()), float(grid.max())])
    fig.grdimage(grid=grid, region=region, projection=proj, cmap=True, shading=True)
    fig.coast(region=region, projection=proj, shorelines="0.5p,black",
              borders=["1/0.6p,black"], frame=["af", "+tEstaciones - Calama 2020 (rojo=usadas, azul=legacy, verde=disponibles)"])

    # 1) disponibles no usadas (triangulos verdes, sin etiqueta)
    otras = [s for s in coords if s not in USED and s not in LEGACY]
    if otras:
        fig.plot(x=[coords[s][1] for s in otras], y=[coords[s][0] for s in otras],
                 style="t0.30c", fill="green2", pen="0.4p,black",
                 region=region, projection=proj, label="Disponibles")
    # 2) legacy (azul, con etiqueta)
    leg = [s for s in coords if s in LEGACY]
    fig.plot(x=[coords[s][1] for s in leg], y=[coords[s][0] for s in leg],
             style="i0.40c", fill="dodgerblue", pen="0.6p,black",
             region=region, projection=proj, label="Legacy (8)")
    # 3) usadas reciente (rojo, con etiqueta)
    usa = [s for s in coords if s in USED]
    fig.plot(x=[coords[s][1] for s in usa], y=[coords[s][0] for s in usa],
             style="t0.40c", fill="red", pen="0.6p,black",
             region=region, projection=proj, label="Usadas reciente (9)")
    # 4) comunes (anillo dorado)
    for s in COMUNES:
        fig.plot(x=[coords[s][1]], y=[coords[s][0]], style="c0.55c", pen="1.6p,gold",
                 region=region, projection=proj)
    # etiquetas solo de usadas + legacy (para no saturar)
    for s in (USED | LEGACY):
        if s in coords:
            fig.text(x=coords[s][1], y=coords[s][0], text=s, font="6p,Helvetica-Bold,black",
                     offset="0/0.30c", fill="white@40", region=region, projection=proj)

    # mecanismo focal + epicentro
    fig.meca(spec=dict(strike=EV["strike"], dip=EV["dip"], rake=EV["rake"], magnitude=EV["mw"]),
             scale="1.0c", longitude=EV["lon"], latitude=EV["lat"], depth=EV["depth"],
             compressionfill="black", region=region, projection=proj)
    fig.plot(x=[EV["lon"]], y=[EV["lat"]], style="a0.5c", fill="yellow", pen="1p,black",
             region=region, projection=proj)
    fig.legend(region=region, projection=proj, position="jTR+o0.2c", box="+gwhite+p0.5p")

    fig.savefig(OUT, dpi=300)
    print(f"mapa -> {OUT}")
    print(f"  usadas={len(usa)} legacy={len(leg)} otras={len(otras)} comunes={sorted(COMUNES)}")


if __name__ == "__main__":
    main()
