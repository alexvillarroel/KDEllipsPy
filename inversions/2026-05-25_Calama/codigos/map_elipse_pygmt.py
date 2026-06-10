#!/usr/bin/env python3
"""
map_elipse_pygmt.py
===================
Mapa PyGMT (estilo del map.py / Calama_Seismicity) con la MISMA region que el
mapa de sismicidad, pero dibujando la ELIPSE de slip de la inversion en vez de
la sismicidad — para comparar la orientacion de la ruptura con la sismicidad.

Dibuja: relieve + costa, el footprint de la elipse (subfallas coloreadas por
slip), el contorno del parche, el EJE MAYOR (orientacion de la ruptura, vía PCA
ponderada por slip), el epicentro y las estaciones dentro del recuadro.

Region inferida del Calama_Seismicity.png (GMT Mercator):
    W=-69.26  E=-68.62  S=-22.87  N=-21.84

Uso:
    python codigos/map_elipse_pygmt.py
"""

from pathlib import Path
import numpy as np
import pygmt

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel
from event import load_event

BASE = Path(__file__).resolve().parent.parent
JOBLIB = BASE / "output" / "inversion_result.joblib"
INPUT_CTL = BASE / "input.ctl"
CTL_FILE = BASE / "event.ctl"
OUT = BASE / "output" / "mapa_elipse_pygmt.png"

# Region del mapa de sismicidad (dada por Alex).
REGION = [-69.25, -68.25, -22.85, -21.85]
PROJ = "M15c"


def main():
    import joblib
    ev = load_event(str(CTL_FILE))
    r = joblib.load(JOBLIB)
    best7 = np.asarray(r.best_model.model, float)
    cfg = r.config if r.config is not None else kde.ConfigParser(filepath=str(INPUT_CTL))

    # Geometria con slip (toda la malla; slip=0 fuera de la elipse).
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), best7, keep_all_sources=True)
    slip = np.array([sp.displacement for sp in geom.source_points], float)
    smax = float(slip.max())
    inside = slip > 1e-6                                  # subfallas DENTRO de la elipse
    src = geom.to_axitra_sources(latlon=True)            # [idx, lat, lon, z]
    lat_in, lon_in, slip_in = src[inside, 1], src[inside, 2], slip[inside]

    # Parametros de la elipse proyectada a superficie por PCA de las subfallas
    # internas (proyeccion afin de una elipse = elipse). Para relleno uniforme,
    # semieje = 2*sqrt(autovalor). (E=East, N=North en metros locales.)
    E = np.array([sp.y_m for sp in geom.source_points])[inside]
    N = np.array([sp.x_m for sp in geom.source_points])[inside]
    cE, cN = E.mean(), N.mean()
    cov = np.cov(np.vstack([E - cE, N - cN]))
    evals, evecs = np.linalg.eigh(cov)
    semi_maj = 2.0 * float(np.sqrt(evals.max()))         # m
    semi_min = 2.0 * float(np.sqrt(evals.min()))         # m
    vmaj = evecs[:, int(np.argmax(evals))]               # (E, N)
    azim = float(np.degrees(np.arctan2(vmaj[0], vmaj[1])) % 180.0)   # azimut CW desde N
    maj_km, min_km = 2 * semi_maj / 1000.0, 2 * semi_min / 1000.0    # ejes COMPLETOS (km)
    c_lat, c_lon = geom._local_to_latlon(cN, cE)         # centro (x_m=North, y_m=East)
    vmin = evecs[:, int(np.argmin(evals))]               # (E, N)
    # Endpoints del eje mayor para la linea de orientacion.
    p_lat, p_lon = [], []
    for s in (+1, -1):
        la, lo = geom._local_to_latlon(cN + s * vmaj[1] * semi_maj,
                                       cE + s * vmaj[0] * semi_maj)
        p_lat.append(la); p_lon.append(lo)
    # Polígono del borde de la elipse (para enmascarar el grid de slip).
    th = np.linspace(0, 2 * np.pi, 180)
    ell_lat, ell_lon = [], []
    for t in th:
        e = cE + semi_maj * np.cos(t) * vmaj[0] + semi_min * np.sin(t) * vmin[0]
        n = cN + semi_maj * np.cos(t) * vmaj[1] + semi_min * np.sin(t) * vmin[1]
        la, lo = geom._local_to_latlon(n, e)
        ell_lat.append(la); ell_lon.append(lo)
    print(f"max slip={smax:.2f} m | elipse: centro ({c_lat:.3f},{c_lon:.3f}) "
          f"ejes {maj_km:.1f}x{min_km:.1f} km  azimut mayor ~{azim:.0f}°")

    # ---------------- Figura ----------------
    fig = pygmt.Figure()
    try:
        grid = pygmt.datasets.load_earth_relief(resolution="03s", region=REGION)
    except Exception:
        grid = pygmt.datasets.load_earth_relief(resolution="15s", region=REGION)
    pygmt.makecpt(cmap="gray", series=[float(grid.min()), float(grid.max())])
    fig.grdimage(grid=grid, region=REGION, projection=PROJ, cmap=True, shading=True)
    fig.coast(region=REGION, projection=PROJ, shorelines="0.5p,black",
              borders=["1/0.6p,black"],
              frame=["af", "+tElipse de slip - Calama 2026-05-25"])

    # --- Slip dentro de la elipse: mallar los puntos e insertarlo como grid ----
    import math
    pad = 0.06
    sp_grid = 0.004
    # bbox alineado a múltiplos exactos del spacing (evita warnings de surface).
    bbox = [math.floor((min(lon_in) - pad) / sp_grid) * sp_grid,
            math.ceil((max(lon_in) + pad) / sp_grid) * sp_grid,
            math.floor((min(lat_in) - pad) / sp_grid) * sp_grid,
            math.ceil((max(lat_in) + pad) / sp_grid) * sp_grid]
    slipgrid = pygmt.surface(x=lon_in, y=lat_in, z=slip_in, region=bbox,
                             spacing=sp_grid, tension=0.35)
    # Máscara: NaN fuera del polígono de la elipse (test punto-en-polígono).
    from matplotlib.path import Path as _Path
    glon = slipgrid["x"].values if "x" in slipgrid.coords else slipgrid["lon"].values
    glat = slipgrid["y"].values if "y" in slipgrid.coords else slipgrid["lat"].values
    GX, GY = np.meshgrid(glon, glat)
    poly = _Path(np.column_stack([ell_lon, ell_lat]))
    inside_grid = poly.contains_points(np.column_stack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    slip_masked = slipgrid.where(inside_grid)
    pygmt.makecpt(cmap="hot", reverse=True, series=[0.0, smax])
    fig.grdimage(grid=slip_masked, region=REGION, projection=PROJ, cmap=True,
                 nan_transparent=True)
    fig.colorbar(frame='x+l"Slip (m)"', position="JMR+o0.6c/0c+w7c")
    # Contorno del borde de la elipse + eje mayor.
    fig.plot(x=ell_lon, y=ell_lat, pen="1.5p,black", region=REGION, projection=PROJ)
    fig.plot(x=p_lon, y=p_lat, pen="1.5p,green2,-", region=REGION, projection=PROJ)

    # Epicentro y estaciones dentro de la region.
    fig.plot(x=[ev["lon"]], y=[ev["lat"]], style="a0.6c", fill="yellow",
             pen="1p,black", region=REGION, projection=PROJ)
    W, Ee, S, Nn = REGION
    for st in cfg.stations.stations:
        if W <= st.longitude <= Ee and S <= st.latitude <= Nn:
            fig.plot(x=[st.longitude], y=[st.latitude], style="t0.45c",
                     fill="dodgerblue", pen="0.8p,black", region=REGION, projection=PROJ)
            fig.text(x=st.longitude, y=st.latitude, text=st.name, font="8p,Helvetica-Bold,black",
                     offset="0/0.35c", region=REGION, projection=PROJ)

    fig.savefig(str(OUT), dpi=300)
    print(f"mapa -> {OUT}")


if __name__ == "__main__":
    main()
