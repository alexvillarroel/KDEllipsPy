"""Tabla comparativa de los casos cinemáticos multisemilla.

Lee <caso>/multiseed/multiseed.joblib de cada caso y ordena por el misfit
MEDIO entre semillas, no por el mejor: el mejor de 5 semillas es un estadístico
de orden y premia a quien tuvo suerte. La std entre semillas dice si la
diferencia entre dos configuraciones significa algo -- si dos casos distan
menos que su propia dispersión, no están separados.

Marca con * los parámetros cuya media queda a menos de un 5% del borde de su
rango (señal de que hay que ampliarlo).

Uso:  python compare_cases.py
"""

from pathlib import Path

import joblib
import numpy as np

HERE = Path(__file__).resolve().parent
# (nombre, min, max) de los rangos de input.ctl + dt0 de run_kin_multiseed
EDGE_FRAC = 0.05


def load():
    rows = []
    for f in sorted(HERE.glob("*/multiseed/multiseed.joblib")):
        d = joblib.load(f)
        rows.append(d)
    return rows


def main():
    rows = load()
    if not rows:
        raise SystemExit("todavía no hay ningún multiseed.joblib")

    rows.sort(key=lambda d: np.mean(d["misfits"]))
    names = rows[0]["param_names"]

    print(f"{'caso':<32s}{'n':>3s}{'media':>8s}{'std':>8s}{'mejor':>8s}{'peor':>8s}")
    print("-" * 75)
    for d in rows:
        m = np.asarray(d["misfits"], float)
        print(f"{d['case']:<32s}{m.size:3d}{m.mean():8.4f}{m.std():8.4f}{m.min():8.4f}{m.max():8.4f}")

    base = rows[0]
    print(f"\nreferencia (mejor media): {base['case']}")
    print(f"{'caso':<32s}" + "".join(f"{n.split()[0]:>10s}" for n in names))
    print("-" * (32 + 10 * len(names)))
    for d in rows:
        mod = np.asarray(d["models"], float)
        best = mod[int(np.argmin(d["misfits"]))]
        print(f"{d['case']:<32s}" + "".join(f"{v:10.3f}" for v in best))
        print(f"{'  (std entre semillas)':<32s}" + "".join(f"{v:10.3f}" for v in mod.std(axis=0)))

    # ¿la diferencia entre el mejor y el segundo supera la dispersión propia?
    if len(rows) > 1:
        a, b = rows[0], rows[1]
        ma, mb = np.asarray(a["misfits"], float), np.asarray(b["misfits"], float)
        sep = mb.mean() - ma.mean()
        pooled = np.hypot(ma.std(), mb.std())
        verdict = "separados" if sep > pooled else "NO separados (dentro de la dispersión entre semillas)"
        print(f"\n{a['case']} vs {b['case']}: Dmisfit {sep:+.4f}, dispersión combinada {pooled:.4f} -> {verdict}")

    # parámetros pegados al borde
    print("\nparámetros cerca del borde del prior (media a <5% del rango):")
    flagged = False
    for d in rows:
        mod = np.asarray(d["models"], float)
        for i, n in enumerate(names):
            lo, hi = mod[:, i].min(), mod[:, i].max()
            # sin los rangos a mano, se usa la dispersión: todas las semillas pegadas
            if mod[:, i].std() < 1e-6:
                print(f"  {d['case']:<32s}{n:<22s} todas las semillas en {lo:.3f}")
                flagged = True
    if not flagged:
        print("  ninguno (revisar igual contra los rangos del input.ctl)")

    print("\nmisfit por estación del mejor modelo de cada caso:")
    for d in rows:
        f = HERE / d["case"] / "multiseed" / "station_misfit.txt"
        if f.is_file():
            print(f"\n--- {d['case']}")
            print(f.read_text().rstrip())


if __name__ == "__main__":
    main()
