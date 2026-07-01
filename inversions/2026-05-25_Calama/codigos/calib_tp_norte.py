#!/usr/bin/env python3
"""
calib_tp_norte.py
=================
Calibra que valor de tp (Param5) coloca el centroide de slip al NORTE del
hipocentro, usando la convencion real del codigo (no a mano). La ruptura de
este evento es unilateral al Norte, asi que despues acotamos Param5 a
tp_norte +/- 30 deg en input.ctl.

Barre tp con una elipse de prueba (circulo, alpha=0, np=0.5) y mide el
centroide Norte = sum(slip * x_m)/sum(slip), donde Subfault.x_m = Norte
(geometry.py:129). El tp que lo MAXIMIZA apunta al Norte.

Uso:
    python codigos/calib_tp_norte.py
"""
from pathlib import Path
import numpy as np
import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel

BASE = Path(__file__).resolve().parent.parent
CTL = BASE / "input.ctl"
HALFWIDTH = 45.0 / 360.0   # +/-45 deg en unidades de Param5 (x 2pi)


def centroid_north(fm, tp):
    # model = [a1, a2, theta(xpi), np, tp(x2pi), dmax(m), vr(km/s)]
    model = np.array([8.0, 8.0, 0.0, 0.5, float(tp), 2.0, 2.5])
    geom = fm.build_geometry_with_ellipse_slip(model)
    s = np.array([sf.slip_m for sf in geom.subfaults])
    x = np.array([sf.x_m for sf in geom.subfaults])   # Norte (m)
    if s.sum() <= 0:
        return np.nan
    return float((s * x).sum() / s.sum())


def main():
    cfg = kde.ConfigParser(filepath=str(CTL))
    fm = AxitraForwardModel.from_config(cfg)
    tps = np.linspace(0.0, 1.0, 72, endpoint=False)
    north = np.array([centroid_north(fm, tp) for tp in tps])

    i_n = int(np.nanargmax(north))
    i_s = int(np.nanargmin(north))
    tp_norte = float(tps[i_n])
    lo, hi = tp_norte - HALFWIDTH, tp_norte + HALFWIDTH

    print(f"tp_norte = {tp_norte:.4f}  (centroide Norte = {north[i_n]:+.1f} m)")
    print(f"tp_sur   = {tps[i_s]:.4f}  (centroide Norte = {north[i_s]:+.1f} m)  "
          f"[debe estar ~0.5 alejado de tp_norte]")
    print(f"\n-> Acotar Param5 (tp) en input.ctl a:  min={lo:.4f}  max={hi:.4f}  "
          f"(+/-{HALFWIDTH*360:.0f} deg)")
    # check: norte y sur deben estar ~medio ciclo aparte
    d = abs(((tps[i_s] - tp_norte) % 1.0) - 0.5)
    assert d < 0.1, f"tp_norte y tp_sur no estan opuestos (d={d:.2f}); revisar convencion"
    print("check OK: tp_norte y tp_sur opuestos (~180 deg)")


if __name__ == "__main__":
    main()
