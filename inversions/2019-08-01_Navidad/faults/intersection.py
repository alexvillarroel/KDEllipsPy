"""Intersección del plano de la falla de Pichilemu (CHAF, tramo principal) con el
plano de ruptura de Navidad (NP1 por el hipocentro del caso), en coordenadas de
la malla (strike desde el borde inicial, dip desde el borde superior) y en
lat/lon. Escribe pichilemu_intersection.json.

Uso: DYN_CASE=isc_tomo/b0.05-0.30_noM14L python intersection.py
"""

import json
import os
from pathlib import Path

import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
CASE = os.environ.get("DYN_CASE", "isc_tomo/b0.05-0.30_noM14L")
K = 111.19  # km/grado


def frame(strike, dip):
    s, d = np.radians(strike), np.radians(dip)
    along = np.array([np.sin(s), np.cos(s), 0.0])                       # E, N, Up
    down = np.array([np.cos(s) * np.cos(d), -np.sin(s) * np.cos(d), -np.sin(d)])  # buzamiento a la derecha
    return along, down, np.cross(along, down)


def main():
    cfg = kde.ConfigParser(str(NAV / CASE / "input.ctl"))
    sp, fp = cfg.source_position, cfg.fault_plane
    lat0, lon0 = sp.latitude, sp.longitude
    to_en = lambda lo, la: np.array([(lo - lon0) * K * np.cos(np.radians(lat0)), (la - lat0) * K, 0.0])
    to_ll = lambda e, n: (lon0 + e / (K * np.cos(np.radians(lat0))), lat0 + n / K)

    gj = json.loads((HERE / "chaf_navidad.geojson").read_text())
    pich = next(f for f in gj["features"] if f["properties"].get("F_name") == "Pichilemu" and f["properties"].get("FT_name") == "Principal")
    pr = pich["properties"]
    q = np.mean([to_en(lo, la) for lo, la in pich["geometry"]["coordinates"]], axis=0)  # punto de la traza (z=0)
    _, _, n2 = frame(float(pr["strike"]), float(pr["dip"]))

    a1, d1, _ = frame(sp.strike, sp.dip)
    p0 = np.array([0.0, 0.0, -sp.depth])  # hipocentro (km)
    # Puntos de Navidad: p = p0 + L*a1 + W*d1 (L, W en km desde el hipocentro).
    # En el plano de Pichilemu si (p - q)·n2 = 0  ->  L*(a1·n2) + W*(d1·n2) = (q - p0)·n2
    cL, cW, rhs = a1 @ n2, d1 @ n2, (q - p0) @ n2
    hx, hy = fp.hx / 1e3, fp.hy / 1e3
    pts = []
    for L in np.linspace(-hx, fp.lx / 1e3 - hx, 201):
        W = (rhs - L * cL) / cW
        if -hy <= W <= fp.ly / 1e3 - hy:
            p = p0 + L * a1 + W * d1
            lo, la = to_ll(p[0], p[1])
            pts.append({"strike_km": L + hx, "dip_km": W + hy, "lon": lo, "lat": la, "depth_km": -p[2]})
    out = {"case": CASE, "pichilemu": {k: pr.get(k) for k in ("strike", "dip", "dipdir", "sense", "activity", "max_z_km")},
           "line_in_mesh": pts}
    (HERE / "pichilemu_intersection.json").write_text(json.dumps(out, indent=1))
    if pts:
        print(f"intersección dentro de la malla: {len(pts)} pts; strike {pts[0]['strike_km']:.1f}→{pts[-1]['strike_km']:.1f} km, "
              f"dip {pts[0]['dip_km']:.1f}→{pts[-1]['dip_km']:.1f} km desde el borde superior; prof. {pts[0]['depth_km']:.1f}→{pts[-1]['depth_km']:.1f} km")
        print(f"  extremos lat/lon: ({pts[0]['lat']:.3f}, {pts[0]['lon']:.3f}) → ({pts[-1]['lat']:.3f}, {pts[-1]['lon']:.3f})")
    else:
        L = np.linspace(-60, 60, 7)
        print("la intersección no cruza la malla; W(L) =", [round((rhs - l * cL) / cW, 1) for l in L], "para L =", L)


if __name__ == "__main__":
    main()
