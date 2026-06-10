#!/usr/bin/env python3
"""
map_elipse.py
=============
Mapa con zoom (cartopy) de la elipse de slip del mejor modelo, PROYECTADA a la
superficie. Dibuja:
  * El footprint de la elipse: subfallas coloreadas por slip (tricontourf),
    proyectadas a lat/lon (la geometria ya incorpora strike + dip).
  * El contorno del borde de la elipse a 0.15*max_slip (misma convencion que
    plot_slip_distribution).
  * Hipocentro / epicentro (estrella) y las estaciones usadas.

Usa el plano agrandado (64 km) para que la elipse salga completa, con los
parametros que convergieron (output/inversion_result.joblib). No re-invierte.

Uso:
    python codigos/map_elipse.py
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import joblib
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel
from event import load_event

# ---------------------------------------------------------------------------- #
BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
OUTPUT = BASE / "output"
JOBLIB = OUTPUT / "inversion_result.joblib"
CTL_FILE = BASE / "event.ctl"

# Plano grande (igual que test_plano_grande.py) para no recortar la elipse.
LX = LY = 64000.0
HX = HY = 32000.0
NX = NY = 60


def main():
    ev = load_event(str(CTL_FILE))
    result = joblib.load(JOBLIB)
    best7 = np.asarray(result.best_model.model, dtype=float)

    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    cfg.fault_plane.lx, cfg.fault_plane.ly = LX, LY
    cfg.fault_plane.hx, cfg.fault_plane.hy = HX, HY
    cfg.fault_plane.nx, cfg.fault_plane.ny = NX, NY

    # Geometria con TODA la malla (slip=0 fuera de la elipse) -> footprint completo.
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.apply_ellipse_model_to_geometry(
        fm.build_geometry(), best7, keep_all_sources=True
    )

    # lat/lon por source point (ya proyectado con strike+dip) + slip.
    src = geom.to_axitra_sources(latlon=True)   # [idx, lat, lon, z]
    lat = src[:, 1]
    lon = src[:, 2]
    slip = np.array([sp.displacement for sp in geom.source_points], dtype=float)
    max_slip = float(slip.max())
    print(f"max_slip = {max_slip:.3f} m  |  {len(slip)} subfallas  |  "
          f"slip>0 en {int((slip > 0).sum())}")

    # --- Figura cartopy ----------------------------------------------------- #
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9, 9))
    ax = plt.axes(projection=proj)

    # Extent: footprint de la elipse (slip>0) + margen.
    m = slip > 0.05 * max_slip
    pad = 1.5
    west, east = lon[m].min() - pad, lon[m].max() + pad
    south, north = lat[m].min() - pad, lat[m].max() + pad
    ax.set_extent([west, east, south, north], crs=proj)

    # Fondo (features de NaturalEarth; si no hay red, se omite sin romper).
    for feat, kw in [
        (cfeature.LAND, dict(facecolor="0.92")),
        (cfeature.BORDERS, dict(edgecolor="0.4", linewidth=0.8)),
        (cfeature.COASTLINE, dict(edgecolor="0.3", linewidth=0.8)),
    ]:
        try:
            ax.add_feature(feat, **kw)
        except Exception as exc:  # noqa: BLE001
            print(f"  [aviso] no se pudo dibujar feature ({exc}); se omite.")

    # Footprint de slip (tricontourf sobre la malla proyectada).
    levels = np.linspace(0.15 * max_slip, max_slip, 12)
    tcf = ax.tricontourf(lon, lat, slip, levels=levels, cmap="hot_r",
                         alpha=0.85, transform=proj, extend="max")
    cb = fig.colorbar(tcf, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Slip (m)")

    # Borde de la elipse (0.15*max_slip).
    ax.tricontour(lon, lat, slip, levels=[0.15 * max_slip], colors="lime",
                  linewidths=2.0, transform=proj)

    # Hipocentro / epicentro y estaciones.
    ax.plot(ev["lon"], ev["lat"], marker="*", markersize=11, color="yellow",
            markeredgecolor="k", markeredgewidth=0.8, transform=proj,
            zorder=6, label="Epicentro/Hipocentro")
    # Solo las estaciones que caen dentro del zoom (el resto estan a 50-200 km).
    n_in = 0
    for st in cfg.stations.stations:
        if west <= st.longitude <= east and south <= st.latitude <= north:
            ax.plot(st.longitude, st.latitude, marker="^", markersize=10,
                    color="dodgerblue", markeredgecolor="k", transform=proj, zorder=6)
            ax.text(st.longitude, st.latitude + 0.012, st.name, fontsize=9,
                    ha="center", transform=proj, zorder=7)
            n_in += 1
    print(f"  estaciones dentro del zoom: {n_in}/{len(cfg.stations.stations)}")

    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5,
                      linestyle="--")
    gl.top_labels = gl.right_labels = False

    ax.set_title(f"Elipse de slip (proyectada a superficie)\n"
                 f"Calama 2026-05-25  ·  prof. {ev['depth']:.0f} km  ·  "
                 f"strike/dip/rake {ev['strike']:.0f}/{ev['dip']:.0f}/{ev['rake']:.0f}",
                 fontsize=12)
    ax.legend(loc="upper right", fontsize=9)

    out = OUTPUT / "mapa_elipse.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"mapa -> {out}")


if __name__ == "__main__":
    main()
