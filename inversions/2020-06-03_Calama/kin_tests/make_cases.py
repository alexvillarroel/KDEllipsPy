"""Casos de prueba del cinematico: plano nodal x hipocentro x modelo 1D x banda.

Escribe kin_tests/<caso>/input.ctl (plantilla ../np1/input.ctl con las secciones
2, 4 y 9 sustituidas), enlaza SAC/ -> ../../np1/SAC y corre ``kde-prep`` en cada
caso para generar su DATA/. kde-prep integra desde la aceleracion CRUDA
(SAC/ACC) con la cadena de ISOLA y la banda del propio input.ctl, que es la
unica preparacion valida: real_disp.py integra dos veces en el tiempo y filtra
solo al final, lo que deja energia antes de la llegada P. Correr despues:

    KIN_CASE=<caso> python run_kin_multiseed.py

Hipocentros
  isola : -23.2536 -68.4991 115.557  (tensor de momento ISOLA, el de np1/np2)
  cat   : -23.2470 -68.5300 123.400  (catalogo, el que cita resumen_trabajo.tex)

Modelos 1D
  potin   : 21 capas hasta 150 km, promedio de la tomografia regional de Potin
            et al. (2024) -- el que documenta el suplementario del paper.
  crust16 : 16 capas hasta 70 km (de ../banda_0.02-0.1/input.ctl). NO cubre la
            profundidad de la fuente (~115 km): todo el plano de falla cae en
            su semiespacio de Vp 8.48 km/s. Se prueba solo para medir cuanto
            cambia el resultado, no como candidato.
"""

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CALAMA = HERE.parent
TEMPLATE = CALAMA / "np1" / "input.ctl"
CRUST16_SRC = CALAMA / "banda_0.02-0.1" / "input.ctl"

PLANES = {"np1": (335.0, 63.0, -85.0), "np2": (143.0, 27.0, -101.0)}
HYPOS = {"isola": (-23.2536, -68.4991, 115.557), "cat": (-23.2470, -68.5300, 123.400)}

# (plano, hipocentro, modelo 1D, banda) -- una sola cosa cambia por caso
CASES = [
    ("np1", "isola", "potin", (0.04, 0.15)),   # base (= configuracion de np1)
    ("np2", "isola", "potin", (0.04, 0.15)),   # base (= configuracion de np2)
    ("np1", "cat", "potin", (0.04, 0.15)),     # hipocentro
    ("np2", "cat", "potin", (0.04, 0.15)),     # hipocentro
    ("np1", "isola", "crust16", (0.04, 0.15)),  # modelo 1D
    ("np1", "isola", "potin", (0.02, 0.10)),   # banda baja
    ("np1", "isola", "potin", (0.05, 0.25)),   # banda ancha
]

LAYERS_RE = r"( Number of layers\s*:\s*)\d+\n(?:\s*[-0-9.eE]+(?:\s+[-0-9.eE]+){5}\s*\n)+"


def layer_block(text):
    """Las lineas de capas (con su conteo) de un input.ctl."""
    m = re.search(LAYERS_RE, text)
    assert m, "no se encontro la seccion 9"
    return m.group(0)[len(m.group(1)):]


def set_field(text, label, value):
    """Reemplaza el valor de una linea ' <label> ... : <valor>' manteniendo el formato."""
    pat = re.compile(rf"^(\s*{re.escape(label)}\s*(?:\([^)]*\))?\s*:\s*)\S+", re.M)
    text, n = pat.subn(lambda m: m.group(1) + value, text, count=1)
    assert n == 1, f"no se pudo fijar {label!r}"
    return text


def case_name(plane, hypo, vel, band):
    return f"{plane}_{hypo}_{vel}_b{band[0]:g}-{band[1]:g}"


def write_case(plane, hypo, vel, band):
    text = TEMPLATE.read_text()
    strike, dip, rake = PLANES[plane]
    lat, lon, dep = HYPOS[hypo]
    for label, val in [("Latitude", f"{lat:.4f}"), ("Longitude", f"{lon:.4f}"),
                       ("Depth", f"{dep:.3f}"), ("Strike", f"{strike:.1f}"),
                       ("Dip", f"{dip:.1f}"), ("Rake", f"{rake:.1f}"),
                       ("Frequency 1 (Freq1)", f"{band[0]:g}"),
                       ("Frequency 2 (Freq2)", f"{band[1]:g}")]:
        text = set_field(text, label, val)
    text = set_field(text, "Event Name", f"Calama2020 {plane} {hypo} {vel}")
    if vel == "crust16":
        text = re.sub(LAYERS_RE, lambda m: m.group(1) + layer_block(CRUST16_SRC.read_text()), text, count=1)

    out = HERE / case_name(plane, hypo, vel, band)
    out.mkdir(exist_ok=True)
    (out / "input.ctl").write_text(text)
    link = out / "SAC"
    if link.is_symlink():
        link.unlink()
    link.symlink_to(CALAMA / "np1" / "SAC")
    subprocess.run([sys.executable, "-m", "kdellipspy.prep_data", str(out)],
                   check=True, stdout=subprocess.DEVNULL)
    nlayers = text[re.search(LAYERS_RE, text).start():].split(":")[1].split()[0]
    print(f"{out.name:34s} strike {strike:5.1f} dip {dip:4.1f} rake {rake:6.1f} "
          f"z {dep:7.3f} km  {nlayers} capas  {band[0]}-{band[1]} Hz")
    return out.name


if __name__ == "__main__":
    wanted = sys.argv[1:]
    names = [write_case(*c) for c in CASES if not wanted or case_name(*c) in wanted]
    print(f"\n{len(names)} casos en {HERE}")
