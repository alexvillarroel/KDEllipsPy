#!/usr/bin/env python3
"""
kde-prep · Genera DATA/real_{vel,disp}_{x,y,z} desde input.ctl + SAC/

Lee la configuración (npts, delta, banda, estaciones, units) de ``input.ctl``
y procesa los SAC de ``SAC/VEL`` (units=vel) o ``SAC/DISP`` (units=disp):
filtra en banda, recorta a la ventana del ctl y exporta SOLO las estaciones
listadas en input.ctl, en su orden.

    kde-prep                      # ./input.ctl + ./SAC → ./DATA
    kde-prep ruta/al/proyecto
    python -m kdellipspy.prep_data
"""
import argparse
from pathlib import Path

from .cli import banner, section, info, ok, warn, err


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="kde-prep", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("project_dir", type=Path, nargs="?", default=Path.cwd())
    p.add_argument("--input-ctl", type=Path, default=None)
    p.add_argument("--sac-dir", type=Path, default=None,
                   help="Carpeta SAC (default <proyecto>/SAC)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Carpeta de salida (default <proyecto>/DATA)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    import kdellipspy as kde

    banner("KDEllipsPy · Generación de datos observados", "input.ctl + SAC → DATA/real_*")

    ctl = args.input_ctl or args.project_dir / "input.ctl"
    sac = args.sac_dir or args.project_dir / "SAC"
    out = args.output or args.project_dir / "DATA"

    section("Validando rutas")
    if not ctl.is_file():
        err(f"No existe input.ctl: {ctl}"); return 1
    if not sac.is_dir():
        err(f"No existe la carpeta SAC: {sac}"); return 1
    ok(f"input.ctl : {ctl}")

    cfg = kde.ConfigParser(filepath=str(ctl))
    od, el = cfg.observed_data, cfg.ellipse
    units = int(od.units)                       # 1: disp, 2: vel
    prefix = "real_disp" if units == 1 else "real_vel"
    # Fuente SIEMPRE aceleracion cruda (SAC/ACC). Integramos aqui, una sola vez:
    # acc->disp = 2 pasos (units=1), acc->vel = 1 paso (units=2). El filtro de
    # banda (del ctl) se aplica tambien una sola vez, con el MISMO orden que veran
    # los sinteticos -> sin doble filtrado asimetrico.
    integrate_n = 3 - units                     # units=1->2, units=2->1
    raw_dir = sac / "ACC"
    if not raw_dir.is_dir():
        err(f"Falta {raw_dir} (kde-prep integra desde SAC/ACC crudo)"); return 1
    ok(f"SAC fuente: {raw_dir} (acc)  →  {out}/{prefix}_{{x,y,z}}  "
       f"[{integrate_n} integracion(es), Butterworth orden {2*integrate_n}]")

    section("Configuración (desde input.ctl)")
    info(f"Banda      : {el.freq1:g}–{el.freq2:g} Hz  |  filtro "
         f"{'acausal' if el.zerophase else 'causal'}")
    info(f"Ventana    : t1={od.t1:g}s  npts={od.npts}  dt={od.delta:g}s  "
         f"({'disp' if units == 1 else 'vel'})")
    info(f"Estaciones : {len(cfg.stations.stations)}  "
         f"[{', '.join(s.name for s in cfg.stations.stations)}]")

    section("Procesando SAC (demean+detrend, taper, integra, filtro, recorte)")
    import numpy as np
    from kdellipspy.core.signal_utils import _load_from_raw

    # Lee los .sac de acc, quita media+tendencia, taper, INTEGRA (integrate_n),
    # filtra en banda (orden 2*integrate_n) una sola vez, recorta a la ventana y
    # apila SOLO las estaciones del input.ctl, en orden. Devuelve (nsta,3,npts) NEZ.
    observed, _ = _load_from_raw(raw_dir, cfg, el.freq1, el.freq2, integrate_n=integrate_n)
    out.mkdir(parents=True, exist_ok=True)
    comps = {"x": 0, "y": 1, "z": 2}   # x=N, y=E, z=Z
    for c, idx in comps.items():
        np.savetxt(out / f"{prefix}_{c}", observed[:, idx, :].flatten(), fmt="%.8e")

    section("Listo")
    for c in comps:
        ok(f"{prefix}_{c}  shape={observed.shape[0]}x{observed.shape[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
