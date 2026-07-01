#!/usr/bin/env python3
"""
kdellipspy CLI
==============
Corre la inversión cinemática (Neighbourhood Algorithm) de KDEllipsPy.

Lee ``input.ctl`` + ``DATA/`` de la carpeta del proyecto (por defecto el
directorio actual) y guarda la solución y las figuras en ``output/``:

    <proyecto>/
    ├── input.ctl
    ├── DATA/
    └── output/
        ├── inversion_result.joblib   ← la solución
        ├── best_model_live.txt
        ├── misfit_breakdown.txt
        └── figures/                  ← PNGs + all_plots.pdf

Ejemplos:
    kdellipspy                       # usa ./input.ctl + ./DATA → ./output
    kdellipspy ruta/al/proyecto
    kdellipspy --freq1 0.06 --freq2 0.15
    kdellipspy --no-plots
    python -m kdellipspy             # equivalente sin instalar el comando
"""

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

# Backend no interactivo: evita que las figuras intenten abrir ventana en terminal.
import matplotlib
matplotlib.use("Agg")

# --------------------------------------------------------------------------- #
#  Utilidades de impresión estilo CLI
# --------------------------------------------------------------------------- #
_T0 = time.perf_counter()


def _clock():
    return datetime.now().strftime("%H:%M:%S")


def _el():
    return f"{time.perf_counter() - _T0:6.1f}s"


def banner(title, subtitle=""):
    line = "═" * 70
    print(f"\n╔{line}╗")
    print(f"║ {title:<68} ║")
    if subtitle:
        print(f"║ {subtitle:<68} ║")
    print(f"╚{line}╝")


def section(name):
    print(f"\n┌─[{_clock()} | {_el()}] {name}")


def info(msg):
    print(f"│  {msg}")


def ok(msg):
    print(f"│  \033[32m✓\033[0m {msg}")


def warn(msg):
    print(f"│  \033[33m!\033[0m {msg}")


def err(msg):
    print(f"\033[31m✗ ERROR:\033[0m {msg}", file=sys.stderr)


# --------------------------------------------------------------------------- #
#  Argumentos
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="kdellipspy",
        description="Inversión cinemática NA (KDEllipsPy). Corre desde la "
                    "carpeta del proyecto (input.ctl + DATA/ → output/).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("project_dir", type=Path, nargs="?", default=Path.cwd(),
                   help="Carpeta del proyecto (contiene input.ctl y DATA/)")
    p.add_argument("--input-ctl", type=Path, default=None,
                   help="Ruta a input.ctl (por defecto <project_dir>/input.ctl)")
    p.add_argument("--data-dir", type=Path, default=None,
                   help="Carpeta de datos (por defecto <project_dir>/DATA)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Carpeta de salida (por defecto <project_dir>/output)")
    p.add_argument("--freq1", type=float, default=None,
                   help="Sobrescribe Freq1 (Hz); por defecto usa input.ctl")
    p.add_argument("--freq2", type=float, default=None,
                   help="Sobrescribe Freq2 (Hz); por defecto usa input.ctl")
    p.add_argument("--n-jobs", type=int, default=None,
                   help="Workers paralelos para NA; por defecto usa input.ctl")
    p.add_argument("--model", default="iasp91",
                   help="Modelo de Tierra para tiempos P/S (azi_times)")
    p.add_argument("--prefer-raw", action="store_true",
                   help="Cargar desde DATA/RAW en vez de archivos planos")
    p.add_argument("--no-plots", action="store_true",
                   help="No generar figuras (solo guardar el .joblib)")
    args = p.parse_args(argv)
    base = args.project_dir
    # Defaults derivados del project_dir si no se pasaron explícitos.
    args.input_ctl = (args.input_ctl or base / "input.ctl").resolve()
    args.data_dir = (args.data_dir or base / "DATA").resolve()
    args.output = (args.output or base / "output").resolve()
    return args


# --------------------------------------------------------------------------- #
#  Programa principal
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)

    banner("KDEllipsPy · Inversión cinemática (NA)")

    # --- Importar la librería (puede tardar un poco) ------------------------ #
    section("Cargando KDEllipsPy")
    try:
        import numpy as np
        import kdellipspy as kde
    except Exception as exc:  # noqa: BLE001
        err(f"No se pudo importar kdellipspy: {exc}")
        return 1
    ok(f"kdellipspy importado (numpy {np.__version__})")

    # --- Validar rutas ------------------------------------------------------ #
    section("Validando rutas")
    if not args.input_ctl.is_file():
        err(f"No existe input.ctl: {args.input_ctl}")
        return 1
    if not args.data_dir.is_dir():
        err(f"No existe la carpeta de datos: {args.data_dir}")
        return 1
    figdir = args.output / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    ok(f"input.ctl : {args.input_ctl}")
    ok(f"DATA      : {args.data_dir}")
    ok(f"output    : {args.output}  (figuras en figures/)")

    # --- Leer configuración ------------------------------------------------- #
    section("Leyendo configuración")
    cfg = kde.ConfigParser(filepath=str(args.input_ctl))
    if args.n_jobs is not None:
        cfg.inversion_process.n_jobs = args.n_jobs

    sp, el, od, ip = (cfg.source_position, cfg.ellipse,
                      cfg.observed_data, cfg.inversion_process)
    f1 = args.freq1 if args.freq1 is not None else el.freq1
    f2 = args.freq2 if args.freq2 is not None else el.freq2
    fase = "causal" if not el.zerophase else "acausal (fase cero)"
    algo = "Neighbourhood Algorithm" if ip.algorithm_type == 0 else "MCMC"

    info(f"Evento      : {getattr(sp, 'event_name', '?')}  "
         f"({sp.latitude:.3f}, {sp.longitude:.3f}, {sp.depth:.1f} km)")
    info(f"Mecanismo   : strike {sp.strike:g} / dip {sp.dip:g} / rake {sp.rake:g}"
         f"  ({'doble cupla' if cfg.moment_tensor.flag == 0 else 'MT completo'})")
    info(f"Banda       : {f1:g}–{f2:g} Hz   |  filtro {fase}")
    info(f"Ventana     : t1={od.t1:g}s  t2={od.t2:g}s  npts={od.npts}  dt={od.delta:g}s"
         f"  (units={'disp' if od.units == 1 else 'vel'})")
    info(f"Estaciones  : {len(cfg.stations.stations)}  "
         f"[{', '.join(s.name for s in cfg.stations.stations)}]")
    info(f"Mod. veloc. : {len(cfg.velocity_model.layers)} capas")
    info(f"Algoritmo   : {algo}  |  n_jobs={ip.n_jobs}")
    ok("Configuración leída")

    # --- Cargar datos observados ------------------------------------------- #
    section("Cargando datos observados")
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(args.input_ctl),
        data_dir=str(args.data_dir),
        prefer_raw=args.prefer_raw,
        freq1=f1,
        freq2=f2,
    )
    ok(f"observed shape = {observed.shape}  (estaciones, comp, npts)")
    ok(f"eje de tiempo  = {time_array[0]:g} … {time_array[-1]:g} s "
       f"({len(time_array)} pts)")

    # --- Tiempos P/S (azimutes) -------------------------------------------- #
    section("Calculando tiempos P/S por estación (azi_times)")
    azi_times_array = kde.build_azi_times_array(config=cfg, model_name=args.model)
    ok(f"azi_times_array shape = {azi_times_array.shape}  (modelo {args.model})")

    # --- Construir modelo --------------------------------------------------- #
    section("Construyendo modelo de inversión")
    model = kde.NAInversionModel(
        config=cfg,
        observed_waveforms=observed,
        time_array=time_array,
        azi_times_array=azi_times_array,
    )
    # Cachear las funciones de Green de la malla completa: se calculan una vez y
    # se reutilizan en cada eval. ~4-8x más rápido.
    model.use_green_cache = True
    # Checkpoint en vivo: mejor modelo a output/best_model_live.txt al mejorar.
    model.checkpoint_path = args.output / "best_model_live.txt"
    n_fwd = ip.ss1 + ip.num_iterations * ip.ss_other
    ok("NAInversionModel listo  (cache de Green ACTIVADO)")
    info(f"checkpoint en vivo → {model.checkpoint_path}")
    info(f"NA: {ip.num_iterations} iter · SS1={ip.ss1} · "
         f"SS={ip.ss_other} · resample={ip.cells_resample} "
         f"→ ~{n_fwd} modelos directos")

    # --- Correr la inversión ------------------------------------------------ #
    section("Corriendo búsqueda NA  (esto puede tardar)")
    t_inv = time.perf_counter()
    result = model.run_na_search()
    dt_inv = time.perf_counter() - t_inv
    ok(f"Inversión terminada en {dt_inv:.1f} s")

    # --- Mejor solución ----------------------------------------------------- #
    section("Mejor solución")
    best = result.best_model
    print("│  ┌──────────────────────────────────────────────┐")
    for name, val in zip(result.param_names, best.model):
        print(f"│  │ {name:<32s} {val:12.4f} │")
    print("│  └──────────────────────────────────────────────┘")
    ok(f"misfit mínimo = {best.misfit:.6g}")

    # --- Guardar resultados ------------------------------------------------- #
    section("Guardando resultados")
    joblib_path = args.output / "inversion_result.joblib"
    result.save(str(joblib_path))
    ok(f"resultado     → {joblib_path}")

    if not args.no_plots:
        section("Generando figuras")
        figs = [
            ("plot",             "na_results.png"),
            ("plot_convergence", "parameter_convergence.png"),
            ("plot_fit",         "waveform_fit.png"),
            ("plot_elipse",      "ellipse_fit.png"),
            ("plot_azimuthal",   "mapa_azimutal.png"),
            ("plot_misfit_breakdown", "misfit_breakdown.png"),
        ]
        for method, fname in figs:
            try:
                getattr(result, method)(show=False, save_path=str(figdir / fname))
                ok(f"figura        → {figdir / fname}")
            except Exception as exc:  # noqa: BLE001
                warn(f"no se pudo generar {fname}: {exc}")

        # Dashboard de una página con todos los paneles.
        try:
            pdf = result.plot_yolo(save_path=str(figdir / "all_plots.pdf"))
            if pdf:
                ok(f"dashboard     → {pdf}")
        except Exception as exc:  # noqa: BLE001
            warn(f"no se pudo generar all_plots.pdf: {exc}")
    else:
        info("(--no-plots) figuras omitidas")

    # --- Resumen ------------------------------------------------------------ #
    banner("Inversión completada",
           f"tiempo total {_el().strip()}  ·  salida en {args.output}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        err("Interrumpido por el usuario (Ctrl-C).")
        sys.exit(130)
    except Exception as exc:  # noqa: BLE001
        err(str(exc))
        import traceback
        traceback.print_exc()
        sys.exit(1)
