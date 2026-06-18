"""
resumen_robustez_egf.py
=======================
Compara los mejores modelos de la bateria de robustez EGF (semillas, jackknife,
EGF cruzada) + el run principal, y escribe el informe final de la noche.

Salida: egf/INFORME_NOCHE.md
Uso:    python codigos/resumen_robustez_egf.py
"""
import os
import glob

import numpy as np
import joblib

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]
BOUNDS = [(2.0, 13.0), (2.0, 13.0), (0.0, 1.0), (0.0, 1.0),
          (0.0, 1.0), (0.5, 4.0), (1.5, 4.5)]

RUNS = [
    ("principal (seed11, 7 est)", "inversion"),
    ("seed22", "inversion_seed22"),
    ("seed33", "inversion_seed33"),
    ("sin AF01", "inversion_sinAF01"),
    ("sin PB19", "inversion_sinPB19"),
    ("EGF M4.6 (3 est)", "inversion_egfM46"),
]


def main():
    filas = []
    for label, d in RUNS:
        fn = os.path.join(BASE, d, "inversion_result.joblib")
        if not os.path.exists(fn):
            print(f"  [!] falta {d}")
            continue
        r = joblib.load(fn)
        b = np.asarray(r.best_model.model, float)
        filas.append((label, float(r.best_model.misfit), b))

    print(f"{'run':28s} {'misfit':>7s} " + " ".join(f"{n:>6s}" for n in NAMES))
    for lb, mis, b in filas:
        print(f"{lb:28s} {mis:7.3f} " + " ".join(f"{v:6.2f}" for v in b))

    # dispersion entre runs del M4.9 (excluye la EGF cruzada para el spread)
    core = np.array([b for lb, _, b in filas if "M4.6" not in lb])
    print("\nDispersion entre runs M4.9 (max-min, % del rango del parametro):")
    veredictos = []
    for i, nm in enumerate(NAMES):
        sp = core[:, i].max() - core[:, i].min()
        b0, b1 = BOUNDS[i]
        frac = 100 * sp / (b1 - b0)
        v = "ROBUSTO" if frac < 15 else ("MEDIO" if frac < 35 else "NO")
        veredictos.append((nm, sp, frac, v))
        print(f"  {nm:6s}: spread {sp:5.2f}  ({frac:3.0f}%)  {v}")

    out = os.path.join(BASE, "INFORME_NOCHE.md")
    with open(out, "w") as f:
        f.write("# Robustez inversion EGF — Calama 2020 (0.1-0.3 Hz)\n\n")
        f.write("| run | misfit | " + " | ".join(NAMES) + " |\n")
        f.write("|---|---|" + "---|" * len(NAMES) + "\n")
        for lb, mis, b in filas:
            f.write(f"| {lb} | {mis:.3f} | "
                    + " | ".join(f"{v:.2f}" for v in b) + " |\n")
        f.write("\n## Dispersion entre runs M4.9 (semillas + jackknife)\n\n")
        f.write("| param | spread | % rango | veredicto |\n|---|---|---|---|\n")
        for nm, sp, frac, v in veredictos:
            f.write(f"| {nm} | {sp:.2f} | {frac:.0f}% | {v} |\n")
        f.write("\nEGF cruzada (M4.6, 3 estaciones) como consistencia "
                "independiente — no entra al spread.\n")
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
