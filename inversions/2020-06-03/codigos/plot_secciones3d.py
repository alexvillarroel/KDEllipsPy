#!/usr/bin/env python3
"""
plot_secciones3d.py
===================
Tres vistas del parche de slip de la inversion (subfallas dentro de la elipse),
coloreadas por slip:
    1) Profundidad vs Longitud  (corte E-O)
    2) Profundidad vs Latitud   (corte N-S)
    3) Vista 3D (lon, lat, profundidad)

Cada subfalla tiene lon/lat (proyeccion del plano de falla a la superficie) y su
profundidad real (z), asi se ve como buza el plano. Marca el hipocentro.

Lee output/inversion_result.joblib.

Uso:
    python codigos/plot_secciones3d.py
"""

from pathlib import Path
import numpy as np
import joblib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel

BASE = Path(__file__).resolve().parent.parent
JOBLIB = BASE / "output" / "inversion_result.joblib"
OUT = BASE / "output" / "secciones_profundidad.png"


def main():
    r = joblib.load(JOBLIB)
    best7 = np.asarray(r.best_model.model, float)
    cfg = r.config
    sp = cfg.source_position

    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), best7, keep_all_sources=True)
    src = geom.to_axitra_sources(latlon=True)               # [idx, lat, lon, z_m]
    slip = np.array([spt.displacement for spt in geom.source_points], float)
    m = slip > 1e-6                                          # dentro de la elipse
    lat, lon, dep = src[m, 1], src[m, 2], src[m, 3] / 1000.0   # depth en km
    s = slip[m]
    smax = float(s.max())
    hyp_lon, hyp_lat, hyp_dep = sp.longitude, sp.latitude, float(sp.depth)
    order = np.argsort(s)   # dibujar primero los de menor slip (los altos arriba)

    cmap = "hot_r"
    fig = plt.figure(figsize=(16, 5.2))

    # 1) Profundidad vs Longitud
    ax1 = fig.add_subplot(1, 3, 1)
    sc = ax1.scatter(lon[order], dep[order], c=s[order], cmap=cmap, vmin=0, vmax=smax,
                     s=22, edgecolor="0.3", linewidth=0.2)
    ax1.plot(hyp_lon, hyp_dep, marker="*", ms=18, color="yellow", mec="k", zorder=6, label="Hipocentro")
    ax1.set_xlabel("Longitud (°)"); ax1.set_ylabel("Profundidad (km)")
    ax1.invert_yaxis(); ax1.set_title("Profundidad – Longitud (corte E–O)")
    ax1.grid(alpha=0.3); ax1.legend(fontsize=8)

    # 2) Profundidad vs Latitud
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.scatter(lat[order], dep[order], c=s[order], cmap=cmap, vmin=0, vmax=smax,
                s=22, edgecolor="0.3", linewidth=0.2)
    ax2.plot(hyp_lat, hyp_dep, marker="*", ms=18, color="yellow", mec="k", zorder=6)
    ax2.set_xlabel("Latitud (°)"); ax2.set_ylabel("Profundidad (km)")
    ax2.invert_yaxis(); ax2.set_title("Profundidad – Latitud (corte N–S)")
    ax2.grid(alpha=0.3)

    # 3) Vista 3D
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    ax3.scatter(lon[order], lat[order], dep[order], c=s[order], cmap=cmap, vmin=0, vmax=smax,
                s=18, edgecolor="0.3", linewidth=0.2)
    ax3.scatter([hyp_lon], [hyp_lat], [hyp_dep], marker="*", s=220, color="yellow", edgecolor="k")
    ax3.set_xlabel("Lon (°)"); ax3.set_ylabel("Lat (°)"); ax3.set_zlabel("Prof. (km)")
    ax3.invert_zaxis(); ax3.set_title("Vista 3D")
    ax3.view_init(elev=22, azim=-60)

    cbar = fig.colorbar(sc, ax=[ax1, ax2, ax3], shrink=0.7, pad=0.02)
    cbar.set_label("Slip (m)")
    fig.suptitle(f"{getattr(sp,'event_name','Evento')}  ·  prof. {hyp_dep:.0f} km  ·  "
                 f"misfit {r.best_model.misfit:.3f}", fontsize=13, y=1.02)

    fig.savefig(str(OUT), dpi=200, bbox_inches="tight")
    print(f"profundidad: {dep.min():.1f}–{dep.max():.1f} km | max slip {smax:.2f} m")
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
