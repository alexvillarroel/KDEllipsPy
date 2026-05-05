from pathlib import Path
import sys
import argparse

from config_parser import ConfigParser
from signal_utils import load_and_filter_observed_data, build_azi_times_array
from inversion_na import NAInversionModel, NAConfig
from inversion_mcmc import MCMCInversionModel, MCMCConfig
from graphics_suite import GraphicsSuite

def step(number: int, total: int, title: str) -> None:
    print("\n" + "=" * 72, flush=True)
    print(f"[{number}/{total}] {title}", flush=True)
    print("=" * 72, flush=True)

def validate_run_dir(run_dir: Path) -> None:
    """Validate that run_dir contains required input.ctl and DATA directory."""
    run_dir = run_dir.resolve()
    
    input_ctl = run_dir / "input.ctl"
    data_dir = run_dir / "DATA"
    
    errors = []
    if not input_ctl.is_file():
        errors.append(f"  ✗ {input_ctl} (not found or not a file)")
    if not data_dir.is_dir():
        errors.append(f"  ✗ {data_dir} (not found or not a directory)")
    
    if errors:
        print("ERROR: Invalid run directory structure", flush=True)
        print(f"Run directory: {run_dir}", flush=True)
        print("Missing required files/directories:", flush=True)
        for error in errors:
            print(error, flush=True)
        print("\nUsage: python main.py --run <RUN_DIR>", flush=True)
        print("       <RUN_DIR> must contain: input.ctl and DATA/", flush=True)
        sys.exit(1)

def get_axitra_dir(run_dir: Path, src_root: Path) -> Path:
    """
    Resolve AXITRA directory: try run_dir first, then kdellipspy/, then parent directory.
    """
    search_paths = [
        run_dir / "AXITRA2024",           # Check in run directory first
        src_root / "AXITRA2024",          # Check in kdellipspy/ (src_root)
        src_root.parent / "AXITRA2024",   # Check in KDEllipsPy/ (parent)
    ]
    
    for path in search_paths:
        if path.is_dir():
            return path
    
    raise RuntimeError(
        f"AXITRA2024 directory not found in:\n" +
        "\n".join(f"  - {p}" for p in search_paths) +
        "\nPlease ensure AXITRA2024 is available in one of these locations."
    )

def main() -> None:
    src_root = Path(__file__).resolve().parent
    
    # Parse arguments first
    parser = argparse.ArgumentParser(
        description="Kinematic inversion pipeline with event-run folder support"
    )
    parser.add_argument(
        "--run",
        type=str,
        default=".",
        help="Event run folder (default: current directory). Must contain input.ctl and DATA/"
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plot generation"
    )
    
    args = parser.parse_args()
    run_dir = Path(args.run).resolve()
    
    print("=" * 72, flush=True)
    print("KINEMATIC INVERSION - PYTHON PIPELINE", flush=True)
    print("=" * 72, flush=True)
    
    step(1, 6, "Validando estructura del directorio de ejecución")
    validate_run_dir(run_dir)
    print(f"✓ Directorio de ejecución: {run_dir}", flush=True)
    
    input_ctl = run_dir / "input.ctl"

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
        data_dir=run_dir / "DATA",
        prefer_raw=False,
    )

    azi_times_array = None
    try:
        azi_times_array = build_azi_times_array(input_ctl_path=input_ctl)
        print(f"azi_times generado en memoria: shape={azi_times_array.shape}", flush=True)
    except Exception as e:
        raise RuntimeError(
            "No se pudo generar azi_times en memoria. Instala obspy para calcular tiempos P/S."
        ) from e
    print(f"Datos listos con forma: {observed_waveforms.shape}", flush=True)

    algo = int(cfg.inversion_process.algorithm_type)
    if algo == 0:
        step(4, 6, "Inicializando inversión (NA / neighpy)")
    elif algo == 1:
        step(4, 6, "Inicializando inversión MCMC (PyMC + ArviZ)")
    else:
        print(f"ERROR: Algorithm type {algo} not supported (0=NA, 1=MCMC).", flush=True)
        sys.exit(1)

    axitra_dir = get_axitra_dir(run_dir, src_root)

    if algo == 0:
        inversion = NAInversionModel(
            str(input_ctl),
            axitra_dir=str(axitra_dir),
            observed_waveforms=observed_waveforms,
            time_array=time_array,
            azi_times_array=azi_times_array,
        )
    else:
        inversion = MCMCInversionModel(
            str(input_ctl),
            axitra_dir=str(axitra_dir),
            observed_waveforms=observed_waveforms,
            time_array=time_array,
            azi_times_array=azi_times_array,
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
    output_dir = run_dir / "output"
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
        
        graphics = GraphicsSuite(base_dir=run_dir, show=show_plots)
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
