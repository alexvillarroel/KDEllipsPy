"""
generar_validadas_banda.py
==========================
Genera estaciones_validadas_b015-040.csv para la inversion EGF en 0.15-0.4 Hz:
toma lat/lon/dist/snr/dt_align_s de estaciones_validadas.csv (banda 0.1-0.3;
dt_align_s es taup, independiente de banda) y actualiza ccP/lagP/ccS/lagS/
valida desde cc_bandas.csv (filas 0.15-0.4). PB02 queda reincorporada.

Uso: python codigos/generar_validadas_banda.py
"""
import os
import csv

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    with open(os.path.join(BASE, "estaciones_validadas.csv")) as f:
        base_rows = {r["sta"]: r for r in csv.DictReader(f)}
    with open(os.path.join(BASE, "cc_bandas.csv")) as f:
        nuevas = {r["sta"]: r for r in csv.DictReader(f)
                  if float(r["f1"]) == 0.15}

    out_rows = []
    for sta, r in base_rows.items():
        if sta not in nuevas:
            continue
        n = nuevas[sta]
        r = dict(r)
        for k in ("ccP_Z", "lagP_s", "ccS_T", "lagS_s", "valida"):
            r[k] = n[k]
        out_rows.append(r)

    out = os.path.join(BASE, "estaciones_validadas_b015-040.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)
    n_ok = sum(int(r["valida"]) for r in out_rows)
    print(f"{n_ok}/{len(out_rows)} validas -> {out}")
    for r in out_rows:
        print(f"  {r['sta']}: ccS={r['ccS_T']} valida={r['valida']}")


if __name__ == "__main__":
    main()
