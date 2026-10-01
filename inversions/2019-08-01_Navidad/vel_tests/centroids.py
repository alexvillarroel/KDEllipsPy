"""Centroide de momento de los mejores modelos cinemáticos (multiseed) vs agencias.

Centroide = promedio de lat/lon/prof. de las subfallas pesado por momento
(mu*A*slip). Agencias: boletín ISC del evento 616616159 (centroides NEIC Mww,
GFZ, GCMT). Escribe centroids.txt.
"""

from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
AGENCIES = {"NEIC Mww": (-34.2433, -72.3015, 23.5), "GFZ": (-34.2080, -72.1250, 22.0), "GCMT": (-34.3200, -72.2400, 25.0)}
CASES = ["loc_csn", "loc_mid", "loc_isc", "loc_neic"]


def km(a, b):
    dy = (a[0] - b[0]) * 111.19
    dx = (a[1] - b[1]) * 111.19 * np.cos(np.radians(0.5 * (a[0] + b[0])))
    return np.hypot(dx, dy), a[2] - b[2]


def main():
    lines = []
    for case in CASES:
        runs = sorted((NAV / case / "output_multiseed").glob("seed*.joblib"))
        if not runs:
            continue
        best = min((joblib.load(r) for r in runs), key=lambda r: r.best_model.misfit)
        cfg = kde.ConfigParser(str(NAV / case / "input.ctl"))
        fm = kde.AxitraForwardModel.from_config(cfg)
        g = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(best.best_model.model[:7]), keep_all_sources=True)
        w = np.array([sf.mu_pa * sf.area_m2 * sf.slip_m for sf in g.subfaults])
        c = (np.average([sf.lat for sf in g.subfaults], weights=w),
             np.average([sf.lon for sf in g.subfaults], weights=w),
             np.average([sf.z_m for sf in g.subfaults], weights=w) / 1e3)
        hyp = (cfg.source_position.latitude, cfg.source_position.longitude, cfg.source_position.depth)
        d_h, _ = km(c, hyp)
        line = (f"{case:9s} misfit {best.best_model.misfit:.4f} | centroide {c[0]:.3f} {c[1]:.3f} {c[2]:5.1f} km"
                f" (a {d_h:4.1f} km del hipocentro) | " +
                " | ".join(f"{n}: {km(c, a)[0]:4.1f} km, dz {km(c, a)[1]:+5.1f}" for n, a in AGENCIES.items()))
        lines.append(line)
    (HERE / "centroids.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
