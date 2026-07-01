#!/usr/bin/env python3
"""
run_inversion.py
================
CLI para correr la inversion cinematica (Neighbourhood Algorithm) del evento
de Calama 2026-05-25 con KDEllipsPy.

Lee ``input.ctl`` + ``DATA/real_disp_{x,y,z}`` (generados por real_disp.py),
corre la busqueda NA y guarda resultados y figuras en ``output/``.

Ejemplos:
    python codigos/run_inversion.py
    python codigos/run_inversion.py --n-jobs 4
    python codigos/run_inversion.py --input-ctl input.ctl --data-dir DATA -o output
    python codigos/run_inversion.py --no-plots
"""

import argparse
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# Backend no interactivo: evita que las figuras intenten abrir ventana en terminal.
import matplotlib
matplotlib.use("Agg")

# --------------------------------------------------------------------------- #
#  Utilidades de impresion estilo CLI
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
def parse_args():
    BASE = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(
        description="Inversion cinematica NA (KDEllipsPy) - Calama 2026.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-ctl", type=Path, default=BASE / "input.ctl",
                   help="Ruta al archivo input.ctl")
    p.add_argument("--data-dir", type=Path, default=BASE / "DATA",
                   help="Carpeta con real_disp_{x,y,z}")
    p.add_argument("-o", "--output", type=Path, default=BASE / "output",
                   help="Carpeta de salida")
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
    return p.parse_args()


# --------------------------------------------------------------------------- #
#  Programa principal
# --------------------------------------------------------------------------- #
def main():
    args = parse_args()
    # Rutas absolutas: load_and_filter_observed_data resuelve data_dir relativo
    # al input.ctl, asi que conviene fijarlas para evitar sorpresas.
    args.input_ctl = args.input_ctl.resolve()
    args.data_dir = args.data_dir.resolve()
    args.output = args.output.resolve()

    banner("KDEllipsPy · Inversion cinematica (NA)",
           "Evento Calama  ·  2026-05-25")

    # --- Importar la libreria (puede tardar un poco) ------------------------ #
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
    args.output.mkdir(parents=True, exist_ok=True)
    ok(f"input.ctl : {args.input_ctl}")
    ok(f"DATA      : {args.data_dir}")
    ok(f"output    : {args.output}")

    # --- Leer configuracion ------------------------------------------------- #
    section("Leyendo configuracion")
    cfg = kde.ConfigParser(filepath=str(args.input_ctl))
    if args.n_jobs is not None:
        cfg.inversion_process.n_jobs = args.n_jobs

    sp, el, od, ip = (cfg.source_position, cfg.ellipse,
                      cfg.observed_data, cfg.inversion_process)
    f1 = args.freq1 if args.freq1 is not None else el.freq1
    f2 = args.freq2 if args.freq2 is not None else el.freq2
    fase = "causal" if not el.zerophase else "acausal (fase cero)"
    algo = "Neighbourhood Algorithm" if ip.algorithm_type == 0 else "MCMC"

    info(f"Evento      : {sp.latitude:.3f}, {sp.longitude:.3f}, {sp.depth:.1f} km")
    info(f"Mecanismo   : strike {sp.strike:g} / dip {sp.dip:g} / rake {sp.rake:g}"
         f"  ({'doble cupla' if cfg.moment_tensor.flag == 0 else 'MT completo'})")
    info(f"Banda       : {f1:g}–{f2:g} Hz   |  filtro {fase}")
    info(f"Ventana     : t1={od.t1:g}s  t2={od.t2:g}s  npts={od.npts}  dt={od.delta:g}s"
         f"  (units={'disp' if od.units == 1 else 'vel'})")
    info(f"Estaciones  : {len(cfg.stations.stations)}  "
         f"[{', '.join(s.name for s in cfg.stations.stations)}]")
    info(f"Mod. veloc. : {len(cfg.velocity_model.layers)} capas")
    info(f"Algoritmo   : {algo}  |  n_jobs={ip.n_jobs}")
    ok("Configuracion leida")

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
    section("Calculando tiempos P/S por estacion (azi_times)")
    azi_times_array = kde.build_azi_times_array(config=cfg, model_name=args.model)
    ok(f"azi_times_array shape = {azi_times_array.shape}  (modelo {args.model})")

    # --- Construir modelo --------------------------------------------------- #
    section("Construyendo modelo de inversion")
    model = kde.NAInversionModel(
        config=cfg,
        observed_waveforms=observed,
        time_array=time_array,
        azi_times_array=azi_times_array,
    )
    # Cachear las funciones de Green de la malla completa: se calculan una vez y
    # se reutilizan en cada eval (solo corre conv por modelo). ~4-8x mas rapido.
    model.use_green_cache = True
    # Checkpoint en vivo: vuelca el mejor modelo a output/best_model_live.txt
    # cada vez que mejora el misfit (para ver los valores durante la corrida).
    model.checkpoint_path = args.output / "best_model_live.txt"
    n_fwd = ip.ss1 + ip.num_iterations * ip.ss_other
    ok(f"NAInversionModel listo  (cache de Green ACTIVADO)")
    info(f"checkpoint en vivo → {model.checkpoint_path}")
    info(f"NA: {ip.num_iterations} iter · SS1={ip.ss1} · "
         f"SS={ip.ss_other} · resample={ip.cells_resample} "
         f"→ ~{n_fwd} modelos directos")

    # --- Correr la inversion ------------------------------------------------ #
    section("Corriendo busqueda NA  (esto puede tardar)")
    t_inv = time.perf_counter()
    result = model.run_na_search()
    dt_inv = time.perf_counter() - t_inv
    ok(f"Inversion terminada en {dt_inv:.1f} s")

    # --- Mejor solucion ----------------------------------------------------- #
    section("Mejor solucion")
    best = result.best_model
    print("│  ┌──────────────────────────────────────────────┐")
    for name, val in zip(result.param_names, best.model):
        print(f"│  │ {name:<32s} {val:12.4f} │")
    print("│  └──────────────────────────────────────────────┘")
    ok(f"misfit minimo = {best.misfit:.6g}")

    # --- Stress drop (crack circular, Eshelby) ------------------------------ #
    # r = promedio de los dos semiejes (a1, a2 ya son semiejes en km).
    # M0 sale del propio modelo (mu = rho*Vs^2 del modelo de velocidades).
    # Delta_sigma = (7/16) * M0 / r^3   [Pa]  (Eshelby 1957, crack circular).
    try:
        m0, mw = model.fm.estimate_total_moment_and_mw(best.model)
        a1, a2 = float(best.model[0]), float(best.model[1])   # semiejes (km)
        r = 0.5 * (a1 + a2) * 1000.0                          # radio equiv. (m)
        dsigma = (7.0 / 16.0) * m0 / r**3                     # Pa
        ok(f"M0 = {m0:.3e} N·m  (Mw {mw:.2f})")
        ok(f"radio equiv. r = {r/1000:.2f} km  (semiejes {a1:.2f}, {a2:.2f} km)")
        ok(f"stress drop Δσ = {dsigma/1e6:.3f} MPa")
    except Exception as exc:  # noqa: BLE001
        warn(f"no se pudo calcular el stress drop: {exc}")

    # --- Guardar resultados ------------------------------------------------- #
    section("Guardando resultados")
    joblib_path = args.output / "inversion_result.joblib"
    result.save(str(joblib_path))
    ok(f"resultado     → {joblib_path}")

    if not args.no_plots:
        figs = [
            ("plot",            "na_results.png"),
            ("plot_convergence", "parameter_convergence.png"),
            ("plot_fit",        "waveform_fit.png"),
            ("plot_elipse",     "ellipse_fit.png"),
            ("plot_ellipse_depth", "ellipse_depth.png"),
            ("plot_ellipse_map",   "ellipse_map.png"),
        ]
        for method, fname in figs:
            try:
                getattr(result, method)(show=False,
                                        save_path=str(args.output / fname))
                ok(f"figura        → {args.output / fname}")
            except Exception as exc:  # noqa: BLE001
                warn(f"no se pudo generar {fname}: {exc}")

        # Mapa de cobertura azimutal (radar polar).
        try:
            from mapa_azimutal import plot_azimuthal_coverage
            _f, gap_max, _info = plot_azimuthal_coverage(
                cfg, sp.latitude, sp.longitude,
                save_path=str(args.output / "mapa_azimutal.png"))
            ok(f"figura        → {args.output / 'mapa_azimutal.png'}  (gap azimutal {gap_max:.0f}°)")
        except Exception as exc:  # noqa: BLE001
            warn(f"no se pudo generar mapa_azimutal.png: {exc}")
    else:
        info("(--no-plots) figuras omitidas")

    # --- Desglose de misfit por estacion/componente ------------------------- #
    section("Desglose de misfit por estacion/componente")
    try:
        from misfit_breakdown import breakdown_from_result
        txt = breakdown_from_result(result, window_s="auto")
        bd_path = args.output / "misfit_breakdown.txt"
        bd_path.write_text(txt + "\n")
        print(txt)
        ok(f"desglose      → {bd_path}")
    except Exception as exc:  # noqa: BLE001
        warn(f"no se pudo generar el desglose de misfit: {exc}")

    # --- PDF con todos los plots (plot_yolo) -------------------------------- #
    if not args.no_plots:
        section("Generando PDF con todos los plots (plot_yolo)")
        try:
            pdf = result.plot_yolo(save_path=str(args.output / "all_plots.pdf"))
            if pdf:
                ok(f"PDF           → {pdf}")
        except Exception as exc:  # noqa: BLE001
            warn(f"no se pudo generar all_plots.pdf: {exc}")

    # --- Resumen ------------------------------------------------------------ #
    banner("Inversion completada",
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
