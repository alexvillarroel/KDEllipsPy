"""
mapa_mecanismos.py
==================
Mapa PyGMT (estilo de map_elipse_pygmt.py de Calama 2026) pero SIN elipses:
dibuja los EPICENTROS de los dos terremotos (Calama 2020-06-03 y 2026-05-25)
con su respectiva PELOTA DE PLAYA (mecanismo focal, convencion Aki-Richards
strike/dip/rake leida de cada event.ctl) y una etiqueta con la FECHA.

Fondo: relieve + costa + fronteras, misma estetica que los mapas de cada evento.

Uso:
    python codigos/mapa_mecanismos.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
from pygmt.params import Box

# ---------------------------------------------------------------------------- #
BASE = Path(__file__).resolve().parent.parent      # comparacion_2020_vs_2026/
INVERSIONS = BASE.parent                            # inversions/
OUT = BASE / "output" / "mapa_mecanismos.png"

EVENTOS = {
    "Calama 2020": dict(dir=INVERSIONS / "2020-06-03", color="dodgerblue"),
    "Calama 2026": dict(dir=INVERSIONS / "2026-05-25_Calama", color="orangered"),
}

# Mw de cada inversion (de comparar_modelos.py / inversion_result.joblib).
MW = {"Calama 2020": 7.0, "Calama 2026": 6.9}

# Region que contiene AMBOS epicentros con margen (lon W/E, lat S/N).
# 2020 esta en lat -23.247 y 2026 en -22.380: el lat debe abarcarlos a los dos.
REGION = [-69.20, -68.20, -23.50, -22.05]

PROJ = "M15c"

# Desplazamiento del beachball respecto al epicentro real (al ESTE), para que
# la estrella del epicentro quede visible; una linea los conecta (offset=True).
DLON, DLAT = 0.16, 0.0


def load_event_ctl(event_dir):
    """Lee event.ctl -> dict lat, lon, depth, strike, dip, rake, date."""
    line = None
    with open(event_dir / "event.ctl") as fh:
        for raw in fh:
            raw = raw.strip()
            if raw and not raw.startswith("#"):
                line = raw
                break
    p = line.split()
    return dict(lat=float(p[0]), lon=float(p[1]), depth=float(p[2]),
                strike=float(p[3]), dip=float(p[4]), rake=float(p[5]),
                date=p[6][:10])   # 'YYYY-MM-DD' del campo UTC


def main():
    evs = {name: {**meta, **load_event_ctl(meta["dir"])}
           for name, meta in EVENTOS.items()}

    fig = pygmt.Figure()
    # --- relieve de fondo + costa ---
    try:
        grid = pygmt.datasets.load_earth_relief(resolution="03s", region=REGION)
    except Exception:
        grid = pygmt.datasets.load_earth_relief(resolution="15s", region=REGION)
    pygmt.makecpt(cmap="gray", series=[float(grid.min()), float(grid.max())])
    fig.grdimage(grid=grid, region=REGION, projection=PROJ, cmap=True, shading=True)
    fig.coast(region=REGION, projection=PROJ, shorelines="0.5p,black",
              borders=["1/0.6p,black"],
              frame=["af", "+tMecanismos focales: Calama 2020 vs 2026"],
              map_scale="jBL+o1c/1c+c-7+w20k+f+lkm+ar",
                  box=Box(fill="white@30", pen="0.5p,gray30,solid", radius="3p"),
               )

    # --- estrella en cada epicentro + pelota de playa (desplazada) + etiqueta ---
    for name, e in evs.items():
        plon, plat = e["lon"] + DLON, e["lat"] + DLAT   # posicion del beachball
        # estrella en el EPICENTRO real
        fig.plot(x=e["lon"], y=e["lat"], style="a0.55c", fill="yellow",
                 pen="1p,black", region=REGION, projection=PROJ)
        # pelota de playa desplazada, unida al epicentro por una linea (offset=True)
        spec = dict(strike=e["strike"], dip=e["dip"], rake=e["rake"],
                    magnitude=MW[name])
        fig.meca(spec=spec, scale="1.2c",
                 longitude=e["lon"], latitude=e["lat"], depth=e["depth"],
                 plot_longitude=plon, plot_latitude=plat, offset=True,
                 compressionfill=e["color"], pen="0.8p,black",
                 region=REGION, projection=PROJ)
        # etiqueta: fecha + Mw + profundidad, bajo la pelota.
        label = f"{e['date']}  Mw {MW[name]:.2f}  z={e['depth']:.0f} km"
        fig.text(x=plon, y=plat, text=label,
                 font="9p,Helvetica-Bold,black", fill="white@30",
                 offset="0/-0.95c", region=REGION, projection=PROJ)

    fig.savefig(str(OUT), dpi=300)
    fig.show()
    print(f"mapa -> {OUT}")
    for name, e in evs.items():
        print(f"  {name}: ({e['lat']:.3f},{e['lon']:.3f}) z={e['depth']:.0f} km  "
              f"sdr={e['strike']:.0f}/{e['dip']:.0f}/{e['rake']:.0f}  Mw {MW[name]:.2f}")


if __name__ == "__main__":
    main()
