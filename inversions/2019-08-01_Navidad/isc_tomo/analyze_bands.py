"""Resumen de la matriz banda x estaciones (cinemático multisemilla).

Por caso: misfit (media±std entre semillas), dt0 y Vr (media±std), centroide
de momento del mejor modelo vs GCMT/NEIC, y misfit por estación del mejor
modelo (sum (obs-sint)^2 / sum obs^2 sobre las 3 componentes, traza completa).
Escribe analysis.txt.
"""

from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
GCMT, NEIC = (-34.3200, -72.2400, 25.0), (-34.2433, -72.3015, 23.5)


def dist(a, b):
    dy = (a[0] - b[0]) * 111.19
    dx = (a[1] - b[1]) * 111.19 * np.cos(np.radians(a[0]))
    return np.hypot(dx, dy), a[2] - b[2]


def main():
    out = []
    for case in sorted(p for p in HERE.glob("b*_*") if (p / "output_multiseed").exists()):
        runs = [joblib.load(r) for r in sorted((case / "output_multiseed").glob("seed*.joblib"))]
        if not runs:
            continue
        ms = np.array([r.best_model.misfit for r in runs])
        mods = np.array([r.best_model.model for r in runs])
        names = list(runs[0].param_names)
        best = runs[int(np.argmin(ms))]
        cfg = kde.ConfigParser(str(case / "input.ctl"))
        fm = kde.AxitraForwardModel.from_config(cfg)
        g = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(best.best_model.model[:7]), keep_all_sources=True)
        w = np.array([sf.mu_pa * sf.area_m2 * sf.slip_m for sf in g.subfaults])
        c = tuple(np.average(v, weights=w) for v in (
            [sf.lat for sf in g.subfaults], [sf.lon for sf in g.subfaults], [sf.z_m / 1e3 for sf in g.subfaults]))
        m0 = w.sum()
        obs, syn = np.asarray(best.observed), np.asarray(best.best_synthetics)
        per_sta = ((obs - syn) ** 2).sum(axis=(1, 2)) / (obs ** 2).sum(axis=(1, 2))
        sta = [s.name for s in cfg.stations.stations]
        i_dt, i_vr = names.index("dt0 (s)"), next(i for i, n in enumerate(names) if n.startswith("vr"))
        dg, dzg = dist(c, GCMT)
        dn, _ = dist(c, NEIC)
        out += [
            f"== {case.name}",
            f"  misfit {ms.mean():.3f} ± {ms.std():.3f} (mejor {ms.min():.3f}) | dt0 {mods[:, i_dt].mean():+.2f} ± {mods[:, i_dt].std():.2f} s"
            f" | Vr {mods[:, i_vr].mean():.2f} ± {mods[:, i_vr].std():.2f} km/s | Mw {(np.log10(m0) - 9.1) / 1.5:.2f}",
            f"  centroide {c[0]:.3f} {c[1]:.3f} {c[2]:.1f} km | a GCMT {dg:.1f} km (dz {dzg:+.1f}) | a NEIC {dn:.1f} km",
            "  misfit por estación: " + "  ".join(f"{n} {v:.2f}" for n, v in zip(sta, per_sta)),
        ]
    (HERE / "analysis.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
