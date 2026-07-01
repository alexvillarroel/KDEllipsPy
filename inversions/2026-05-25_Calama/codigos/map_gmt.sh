#!/usr/bin/env python3
"""
map_elipse_3d.py
================
Vista 3D de la interfaz de subduccion (el plano de falla de la inversion) con la
ELIPSE de slip dibujada SOBRE el plano inclinado, en su orientacion real 3D.

A diferencia de map_elipse_pygmt.py (que proyecta todo a superficie), aqui se
conserva la profundidad de cada subfalla, de modo que se ve como buza la
interfaz y como queda la elipse de ruptura apoyada sobre ella.

Dibuja:
  - la interfaz de subduccion como superficie tenue (todas las subfallas),
  - el slip dentro de la elipse (scatter coloreado por slip),
  - el contorno 3D de la elipse y su EJE MAYOR (orientacion de la ruptura),
  - la proyeccion del contorno en superficie (Up=0) como referencia al mapa 2D,
  - el epicentro.

Uso:
    python codigos/map_elipse_3d.py
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registra la proyeccion 3d)

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel
from event import load_event

BASE = Path(__file__).resolve().parent.parent
JOBLIB = BASE / "output" / "inversion_result.joblib"
INPUT_CTL = BASE / "input.ctl"
CTL_FILE = BASE / "event.ctl"
OUT = BASE / "output" / "mapa_elipse_3d.png"

# Exageracion vertical (1.0 = escala real). Subelo (p.ej. 2-3) si el buzamiento
# se ve muy plano por la diferencia de escala horizontal vs. vertical.
VE = 1.0


def _source_depth_m(geom, src):
    """Profundidad (m, positiva hacia abajo) por subfalla, robusto a unidades."""
    sp0 = geom.source_points[0]
    for attr in ("z_m", "depth_m"):
        if hasattr(sp0, attr):
            return np.array([getattr(sp, attr) for sp in geom.source_points], float)
    z = np.asarray(src[:, 3], float)
    # to_axitra_sources suele entregar z en km; si parece km, pasar a m.
    if np.nanmax(np.abs(z)) < 1000.0:
        z = z * 1000.0
    return z


def main():
    import joblib
    ev = load_event(str(CTL_FILE))
    r = joblib.load(JOBLIB)
    best7 = np.asarray(r.best_model.model, float)
    cfg = r.config if r.config is not None else kde.ConfigParser(filepath=str(INPUT_CTL))

    # --- Geometria con slip (toda la malla; slip=0 fuera de la elipse) ---
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), best7, keep_all_sources=True)
    slip = np.array([sp.displacement for sp in geom.source_points], float)
    smax = float(slip.max())
    inside = slip > 1e-6
    src = geom.to_axitra_sources(latlon=True)            # [idx, lat, lon, z]

    # --- Coordenadas locales 3D (km). E=Este, N=Norte, Up = -profundidad ---
    E_m = np.array([sp.y_m for sp in geom.source_points], float)   # Este (m)
    N_m = np.array([sp.x_m for sp in geom.source_points], float)   # Norte (m)
    z_m = _source_depth_m(geom, src)                                # prof (m, +abajo)
    Ek, Nk, Up = E_m / 1000.0, N_m / 1000.0, -z_m / 1000.0

    Ei, Ni, Ui, si = Ek[inside], Nk[inside], Up[inside], slip[inside]

    # --- PCA 3D de las subfallas internas: la elipse vive en el plano de falla ---
    # (la subfalla esta sobre un plano => el 3er autovalor ~ 0 = espesor nulo;
    #  los dos mayores definen los ejes mayor/menor dentro del plano inclinado.)
    P = np.vstack([Ei, Ni, Ui]).T
    c = P.mean(axis=0)
    cov = np.cov((P - c).T)
    evals, evecs = np.linalg.eigh(cov)               # ascendente
    order = np.argsort(evals)[::-1]                  # descendente
    evals, evecs = evals[order], evecs[:, order]
    v1, v2 = evecs[:, 0], evecs[:, 1]                # ejes en el plano (mayor, menor)
    semi_maj = 2.0 * float(np.sqrt(max(evals[0], 0.0)))   # km
    semi_min = 2.0 * float(np.sqrt(max(evals[1], 0.0)))   # km

    # Contorno 3D de la elipse y eje mayor (sobre el plano inclinado)
    th = np.linspace(0, 2 * np.pi, 200)
    ell = c[:, None] + semi_maj * np.cos(th) * v1[:, None] + semi_min * np.sin(th) * v2[:, None]
    axis_pts = np.column_stack([c - semi_maj * v1, c + semi_maj * v1])

    # Azimut del eje mayor (CW desde N), proyectando v1 a la horizontal (E,N)
    azim = float(np.degrees(np.arctan2(v1[0], v1[1])) % 180.0)
    print(f"max slip={smax:.2f} m | elipse 3D: centro E,N,Up=({c[0]:.1f},{c[1]:.1f},{c[2]:.1f}) km | "
          f"ejes {2*semi_maj:.1f} x {2*semi_min:.1f} km | azimut mayor ~{azim:.0f}° | "
          f"prof. centro ~{-c[2]:.1f} km")

    # --- Epicentro: lon/lat -> E,N por ajuste afin con las subfallas (robusto) ---
    epi = None
    try:
        lat_s, lon_s = src[:, 1], src[:, 2]
        A = np.column_stack([lon_s, lat_s, np.ones_like(lon_s)])
        cE, *_ = np.linalg.lstsq(A, E_m, rcond=None)
        cN, *_ = np.linalg.lstsq(A, N_m, rcond=None)
        eE = float(np.array([ev["lon"], ev["lat"], 1.0]) @ cE) / 1000.0
        eN = float(np.array([ev["lon"], ev["lat"], 1.0]) @ cN) / 1000.0
        edep = float(ev.get("depth", -c[2]))         # km (fallback: prof. del centro)
        epi = (eE, eN, -edep)
    except Exception as e:
        print("epicentro: no se pudo ubicar en 3D:", e)

    # ---------------- Figura 3D ----------------
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Interfaz de subduccion: superficie tenue con todas las subfallas
    try:
        ax.plot_trisurf(Ek, Nk, Up, color="0.6", alpha=0.20,
                        linewidth=0, antialiased=True, shade=False)
    except Exception:
        ax.scatter(Ek, Nk, Up, s=2, c="0.7", alpha=0.4)

    # Slip dentro de la elipse (alto=brillante)
    p = ax.scatter(Ei, Ni, Ui, c=si, cmap="hot", vmin=0, vmax=smax,
                   s=18, depthshade=False, edgecolors="none")
    cb = fig.colorbar(p, ax=ax, shrink=0.6, pad=0.08)
    cb.set_label("Slip (m)")

    # Contorno y eje mayor de la elipse en 3D
    ax.plot(ell[0], ell[1], ell[2], color="black", lw=2.0, label="Elipse de slip")
    ax.plot(axis_pts[0], axis_pts[1], axis_pts[2], color="lime", lw=2.0, ls="--",
            label="Eje mayor (ruptura)")

    # Proyeccion del contorno a superficie (Up=0), como puente con el mapa 2D
    ax.plot(ell[0], ell[1], np.zeros_like(ell[0]), color="0.4", lw=1.0, ls=":")

    # Epicentro
    if epi is not None:
        ax.scatter([epi[0]], [epi[1]], [epi[2]], marker="*", s=260,
                   c="yellow", edgecolors="black", linewidths=1.0, label="Epicentro")

    ax.set_xlabel("Este (km)")
    ax.set_ylabel("Norte (km)")
    ax.set_zlabel("Profundidad (km)  [Up = -z]")
    ax.set_title("Elipse de slip sobre la interfaz de subduccion - Calama 2026-05-25")
    ax.legend(loc="upper left")

    # Aspecto a escala real (VE para exagerar la vertical si hace falta)
    xr, yr = np.ptp(Ek), np.ptp(Nk)
    zr = max(np.ptp(Up), 1.0) * VE
    ax.set_box_aspect((xr, yr, zr))
    ax.view_init(elev=22, azim=-60)

    fig.tight_layout()
    fig.savefig(str(OUT), dpi=300)
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
