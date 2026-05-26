import sys
from pathlib import Path
import argparse

# --- Bloque de resolución de rutas para que los imports funcionen siempre ---
# Agregamos la raíz del proyecto al path si no está presente
src_root = Path(__file__).resolve().parent.parent
if str(src_root) not in sys.path:
    sys.path.insert(0, str(src_root))

from kdellipspy import (
    ConfigParser,
    load_and_filter_observed_data,
    build_azi_times_array,
    NAInversionModel,
    NAConfig,
    MCMCInversionModel,
    MCMCConfig,
    GraphicsSuite
)

def step(number: int, total: int, title: str) -> None:
    print("\n" + "=" * 72, flush=True)
    print(f"[{number}/{total}] {title}", flush=True)
    print("=" * 72, flush=True)

def validate_paths(input_dir: Path, data_dir: Path) -> Path:
    """Validate that input_dir contains input.ctl and data_dir exists."""
    input_dir = input_dir.resolve()
    data_dir = data_dir.resolve()
    
    input_ctl = input_dir / "input.ctl"
    
    errors = []
    if not input_ctl.is_file():
        errors.append(f"  ✗ {input_ctl} (no se encontró o no es un archivo)")
    if not data_dir.is_dir():
        errors.append(f"  ✗ {data_dir} (no se encontró o no es un directorio)")
    
    if errors:
        print("ERROR: Estructura de directorios inválida", flush=True)
        print(f"Directorio input: {input_dir}", flush=True)
        print(f"Directorio datos: {data_dir}", flush=True)
        print("Faltan archivos o directorios requeridos:", flush=True)
        for error in errors:
            print(error, flush=True)
        print("\nUso: python cli.py <INPUT_DIR> [DATA_DIR]", flush=True)
        print("       <INPUT_DIR> debe contener: input.ctl", flush=True)
        print("       [DATA_DIR] por defecto es <INPUT_DIR>/DATA", flush=True)
        sys.exit(1)
    
    return input_ctl

def get_axitra_dir(input_dir: Path, src_root: Path) -> Path:
    """
    Resolve axitra directory: try input_dir first, then kdellipspy/, then parent directory.
    """
    search_paths = [
        input_dir / "axitra" / "src",           # Check in input directory first
        src_root / "axitra" / "src",          # Check in kdellipspy/ (src_root)
        src_root.parent / "axitra" / "src",   # Check in KDEllipsPy/ (parent)
        input_dir / "axitra",
        src_root / "axitra",
    ]
    
    for path in search_paths:
        if path.is_dir():
            return path
    
    raise RuntimeError(
        f"axitra directory not found in:\n" +
        "\n".join(f"  - {p}" for p in search_paths) +
        "\nPlease ensure axitra is available in one of these locations."
    )

def main() -> None:
    src_root = Path(__file__).resolve().parent
    
    # Parse arguments first
    parser = argparse.ArgumentParser(
        description="Pipeline de inversión cinemática con soporte para directorios de entrada y datos"
    )
    parser.add_argument(
        "input_dir",
        type=str,
        nargs="?",
        default=".",
        help="Carpeta que contiene input.ctl (por defecto: .)"
    )
    parser.add_argument(
        "data_dir",
        type=str,
        nargs="?",
        help="Carpeta que contiene los datos observados (por defecto: <INPUT_DIR>/DATA)"
    )
    parser.add_argument(
        "--run",
        type=str,
        help="Legacy: Carpeta de ejecución (contiene input.ctl y DATA/)"
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Saltar la generación de gráficos"
    )
    
    args = parser.parse_args()
    
    # Resolvimiento de rutas basado en los nuevos parámetros
    if args.run:
        input_dir = Path(args.run).resolve()
        data_dir = input_dir / "DATA"
    else:
        input_dir = Path(args.input_dir).resolve()
        if args.data_dir:
            data_dir = Path(args.data_dir).resolve()
        else:
            data_dir = input_dir / "DATA"

    print("=" * 72, flush=True)
    print("KINEMATIC INVERSION - PYTHON PIPELINE", flush=True)
    print("=" * 72, flush=True)
    
    step(1, 6, "Validando rutas de entrada y datos")
    input_ctl = validate_paths(input_dir, data_dir)
    print(f"✓ Directorio de configuración: {input_dir}", flush=True)
    print(f"✓ Directorio de datos: {data_dir}", flush=True)

    step(2, 6, "Cargando configuración")
    cfg = ConfigParser(str(input_ctl))
    print(f"Evento: {cfg.source_position.event_name}", flush=True)
    print(f"Estaciones: {len(cfg.stations.stations)}", flush=True)
    print(f"Frecuencia: {cfg.ellipse.freq1} - {cfg.ellipse.freq2} Hz", flush=True)

    step(3, 6, "Preparando datos observados")
    observed_waveforms, time_array = load_and_filter_observed_data(
        freq1=cfg.ellipse.freq1,
        freq2=cfg.ellipse.freq2,
        input_ctl_path=input_ctl,
        data_dir=data_dir,
        prefer_raw=False,
    )
    print(f"Datos observados listos con forma: {observed_waveforms.shape}", flush=True)

    algo = int(cfg.inversion_process.algorithm_type)
    if algo == 0:
        step(4, 6, "Inicializando inversión (NA / neighpy)")
    elif algo == 1:
        step(4, 6, "Inicializando inversión MCMC (PyMC + ArviZ)")
    else:
        print(f"ERROR: Algorithm type {algo} not supported (0=NA, 1=MCMC).", flush=True)
        sys.exit(1)

    axitra_dir = get_axitra_dir(input_dir, src_root)

    if algo == 0:
        inversion = NAInversionModel(
            config=cfg,
            axitra_dir=str(axitra_dir),
            observed_waveforms=observed_waveforms,
            time_array=time_array,
        )
    else:
        inversion = MCMCInversionModel(
            config=cfg,
            axitra_dir=str(axitra_dir),
            observed_waveforms=observed_waveforms,
            time_array=time_array,
        )

    step(5, 6, "Ejecutando inversión")
    if algo == 0:
        na_config = NAConfig(
            n_samples_initial=cfg.inversion_process.ss1,
            n_samples_iteration=cfg.inversion_process.ss_other,
            n_iterations=cfg.inversion_process.num_iterations,
            n_cells_resample=cfg.inversion_process.cells_resample,
            random_seed=None,
        )
        result = inversion.run_na_search(na_config)
    else:
        mc_config = MCMCConfig(
            total_steps=cfg.inversion_process.mcmc_total_steps,
            burn_in=cfg.inversion_process.mcmc_burn_in,
            proposal_scale=cfg.inversion_process.mcmc_proposal_scale,
            thin=cfg.inversion_process.mcmc_thin,
            random_seed=None,
            keep_axitra_files=False,
            chains=cfg.inversion_process.mcmc_chains,
        )
        result = inversion.run_mcmc_search(mc_config)

    best = result.best_model
    if best is None:
        print("ERROR: la inversión no devolvió modelos evaluados.", flush=True)
        sys.exit(1)

    print(f"Mejor misfit: {best.misfit:.6e}", flush=True)
    print(f"Parámetros del mejor modelo:", flush=True)
    for name, value in zip(inversion.param_names, best.model):
        print(f"  {name:20s} = {value:8.3f}", flush=True)

    best_geom = inversion._build_geometry_from_parameters(best.model)
    diag = getattr(best_geom, "slip_scale_diagnostics", None)
    if isinstance(diag, dict):
        print("\n--- Slip / MT scaling (mejor modelo) ---", flush=True)
        for key in sorted(diag.keys()):
            print(f"  {key}: {diag[key]}", flush=True)

    step(6, 6, "Exportando resultados y generando gráficos")
    output_dir = input_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "na_results.json"
    csv_path = output_dir / "na_results.csv"

    result.export_results(json_path)
    result.export_csv(csv_path)

    print(f"Resultados guardados en: {json_path}", flush=True)
    print(f"Resultados guardados en: {csv_path}", flush=True)

    try:
        import matplotlib
        # Check if we're in an interactive environment
        is_interactive = hasattr(sys, 'ps1') or 'ipython' in sys.modules or 'pytest' in sys.modules
        show_plots = is_interactive and not args.no_plot
        
        graphics = GraphicsSuite(base_dir=input_dir, show=show_plots)
        graphics.plot_na_results(result)

        # --- Gráfico de sismogramas observados vs sintéticos del mejor modelo ---
        if inversion.best_synthetics is not None:
            station_names = [st.name for st in cfg.stations.stations]
            graphics.plot_waveform_fit(
                observed=inversion.observed_waveforms,
                synthetic=inversion.best_synthetics,
                time_array=inversion.time_array,
                station_names=station_names,
                misfit=result.best_model.misfit if result.best_model else None,
            )
        else:
            print("⚠ No se generaron sintéticos del mejor modelo (best_synthetics es None).", flush=True)
        
        if show_plots:
            print("✓ Gráficos mostrados", flush=True)
        else:
            print("✓ Gráficos generados (sin mostrar en modo no-interactivo)", flush=True)
    except Exception as e:
        print(f"⚠ Advertencia al generar gráficos: {e}", flush=True)
        print("  Los resultados JSON/CSV se han guardado correctamente.", flush=True)

    print("\n✓ Inversión terminada correctamente.", flush=True)


if __name__ == "__main__":
    main()
