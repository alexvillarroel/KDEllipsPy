"""Caso dinámico derivado de un caso cinemático, con la malla de subfallas remuestreada.

El input.ctl de Calama trae Nx = Ny = 40 sobre 40 x 40 km, es decir subfallas de
1 km. Para el dinámico esa malla es el doble gasto por nada: es a la vez la
grilla gruesa del prestress y la malla de fuentes de axitra, y el costo de la
base de respuestas impulsivas va como 2 x Nx x Ny convoluciones (1600 subfallas
-> 3200 llamadas, ~20 min, contra ~5 min con 20 x 20). A 0.15 Hz la longitud de
onda más corta es ~29 km (Vs 4.4 km/s), así que subfallas de 2 km sobran -- y
2 km es justo la resolución validada en Navidad, lo que deja los rangos de a, b
de na_dynamic.py (en puntos de la grilla gruesa) directamente comparables.

Lx/dh y Ly/dh (80 x 80 celdas de 500 m) siguen siendo múltiplos de Nx, Ny.

Uso:  python make_dyn_case.py kin_tests/np1_isola_potin_b0.04-0.15 [N]
      -> ../dyn_cases/<nombre>/{input.ctl, DATA, SAC}
"""

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CALAMA = HERE.parent
NSUB = 20  # subfallas por lado (40 km / 20 = 2 km)


def set_field(text, label, value):
    pat = re.compile(rf"^(\s*{re.escape(label)}\s*(?:\([^)]*\))?\s*:\s*)\S+", re.M)
    text, n = pat.subn(lambda m: m.group(1) + value, text, count=1)
    assert n == 1, f"no se pudo fijar {label!r}"
    return text


def make(kin_case, nsub=NSUB):
    src = CALAMA / kin_case
    assert (src / "input.ctl").is_file(), f"no existe {src}/input.ctl"
    text = (src / "input.ctl").read_text()
    for label in ("Number of subfaults along strike (Nx)", "Number of subfaults along dip (Ny)"):
        text = set_field(text, label, str(nsub))

    out = CALAMA / "dyn_cases" / Path(kin_case).name
    out.mkdir(parents=True, exist_ok=True)
    # kde-prep busca integracion.py en <caso>/codigos o <caso>/../codigos
    cod = out.parent / "codigos"
    if not cod.exists():
        cod.symlink_to(CALAMA / "np1" / "codigos")
    (out / "input.ctl").write_text(text)
    link = out / "SAC"
    if link.is_symlink():
        link.unlink()
    link.symlink_to(CALAMA / "np1" / "SAC")
    # DATA depende solo de la banda y de las estaciones, no de Nx/Ny, pero se
    # regenera con kde-prep para que el caso sea autocontenido y reproducible.
    subprocess.run([sys.executable, "-m", "kdellipspy.prep_data", str(out)],
                   check=True, stdout=subprocess.DEVNULL)

    import kdellipspy as kde
    cfg = kde.ConfigParser(str(out / "input.ctl"))
    fp = cfg.fault_plane
    assert int(round(fp.lx / 500.0)) % fp.nx == 0 and int(round(fp.ly / 500.0)) % fp.ny == 0, \
        "Lx/dh y Ly/dh deben ser múltiplos de Nx, Ny"
    print(f"{out.relative_to(CALAMA)}: {fp.nx}x{fp.ny} subfallas de "
          f"{fp.lx / fp.nx / 1e3:g} x {fp.ly / fp.ny / 1e3:g} km  "
          f"(strike {cfg.source_position.strike:g} / dip {cfg.source_position.dip:g} "
          f"/ rake {cfg.source_position.rake:g})")
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    make(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else NSUB)
