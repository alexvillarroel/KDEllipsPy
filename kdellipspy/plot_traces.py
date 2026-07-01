#!/usr/bin/env python3
"""
kde-traces · Grafica las trazas integradas de DATA/ (salida de kde-prep)

Lee DATA/real_{disp,vel}_{x,y,z} (segun units de input.ctl) y dibuja una grilla
de nsta filas x 3 columnas (N, E, Z), una fila por estacion, en el orden del ctl.
La banda/ventana/estaciones salen de input.ctl (misma fuente que kde-prep).

    kde-traces                      # muestra el plot (ventana interactiva)
    kde-traces -o salida.png        # lo guarda en vez de mostrarlo
    kde-traces ruta/al/proyecto
    python -m kdellipspy.plot_traces
"""
import argparse
from pathlib import Path

from .cli import banner, section, info, ok, err


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="kde-traces", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("project_dir", type=Path, nargs="?", default=Path.cwd())
    p.add_argument("--input-ctl", type=Path, default=None)
    p.add_argument("--data-dir", type=Path, default=None,
                   help="Carpeta DATA (default <proyecto>/DATA)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Guarda el PNG aqui en vez de mostrarlo en pantalla")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    import numpy as np
    import matplotlib
    import kdellipspy as kde          # importa primero: el paquete fuerza Agg al cargar
    # Elegimos backend DESPUES (con force) para que el Agg del paquete no nos pise.
    if args.output is not None:
        matplotlib.use("Agg", force=True)
    else:
        for be in ("QtAgg", "TkAgg"):  # interactivos disponibles en el env
            try:
                matplotlib.use(be, force=True); break
            except Exception:
                continue
    import matplotlib.pyplot as plt

    banner("KDEllipsPy · Trazas observadas", "DATA/real_* → grilla nsta×3 (N,E,Z)")

    ctl = args.input_ctl or args.project_dir / "input.ctl"
    data = args.data_dir or args.project_dir / "DATA"
    out = args.output                 # None -> mostrar; ruta -> guardar

    section("Validando rutas")
    if not ctl.is_file():
        err(f"No existe input.ctl: {ctl}"); return 1
    if not data.is_dir():
        err(f"No existe la carpeta DATA: {data}"); return 1
    ok(f"input.ctl : {ctl}")

    cfg = kde.ConfigParser(filepath=str(ctl))
    od, el = cfg.observed_data, cfg.ellipse
    units = int(od.units)
    prefix = "real_disp" if units == 1 else "real_vel"
    mag = "desplazamiento (m)" if units == 1 else "velocidad (m/s)"
    npts, dt, t1 = int(od.npts), float(od.delta), float(od.t1)
    names = [s.name for s in cfg.stations.stations]
    nsta = len(names)
    t = t1 + np.arange(npts) * dt

    section("Cargando DATA")
    comps = ("x", "y", "z")          # x=N, y=E, z=Z
    arrs = {}
    for c in comps:
        f = data / f"{prefix}_{c}"
        if not f.is_file():
            err(f"Falta {f}"); return 1
        v = np.loadtxt(f)
        if v.size != nsta * npts:
            err(f"{f.name}: {v.size} valores, esperaba {nsta}×{npts}={nsta*npts} "
                f"(¿coincide DATA con input.ctl?)"); return 1
        arrs[c] = v.reshape(nsta, npts)
    ok(f"{prefix}_{{x,y,z}}  {nsta} estaciones × {npts} muestras")

    section("Graficando")
    labels = ["N", "E", "Z"]
    fig, axes = plt.subplots(nsta, 3, figsize=(12, max(2.0, 1.3 * nsta)),
                             sharex=True, squeeze=False)
    for i, sta in enumerate(names):
        # Mismo ylim para las 3 componentes de la estacion (comparables entre si).
        row = np.array([arrs[c][i] for c in comps])
        m = np.abs(row).max() or 1.0
        for j, c in enumerate(comps):
            ax = axes[i][j]
            ax.plot(t, arrs[c][i], lw=0.8, color="C0")
            ax.set_ylim(-1.1 * m, 1.1 * m)
            ax.grid(alpha=0.25)
            if j == 0:
                ax.set_ylabel(sta, rotation=0, ha="right", va="center",
                              fontsize=9, labelpad=18)
            if i == 0:
                ax.set_title(labels[j], fontsize=11, weight="bold")
            if i == nsta - 1:
                ax.set_xlabel("t (s)")
            ax.tick_params(labelsize=7)
    fig.suptitle(f"Trazas observadas — {mag} · banda {el.freq1:g}–{el.freq2:g} Hz "
                 f"· {nsta} estaciones", fontsize=12, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.99))

    section("Listo")
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=140)
        ok(f"guardado en {out}")
    else:
        ok("mostrando ventana (cierra para terminar)")
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
