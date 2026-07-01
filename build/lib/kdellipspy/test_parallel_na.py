
import numpy as np
import os
import sys
from pathlib import Path

# --- Bloque de resolución de rutas ---
src_root = Path(__file__).resolve().parent.parent
if str(src_root) not in sys.path:
    sys.path.insert(0, str(src_root))

from kdellipspy import ConfigParser, NAInversionModel, NAConfig

def create_mock_config():
    """Crea un diccionario de configuración minimal para la prueba."""
    params = {
        'source_position': {
            'event_name': 'TestEvent',
            'latitude': -33.45,
            'longitude': -70.66,
            'depth': 10.0,
            'strike': 0.0,
            'dip': 90.0,
            'rake': 0.0
        },
        'fault_plane': {
            'Lx': 20000.0,
            'Ly': 20000.0,
            'nx': 2,
            'ny': 2,
            'hx': 10000.0,
            'hy': 10000.0
        },
        'ellipse': {
            'n_ellipses': 1,
            'initial_slip': 0,
            'slip_shape': 1,
            'freq1': 0.1,
            'freq2': 0.5,
            't0': 2.0
        },
        'observed_data': {
            't1': 0.0,
            't2': 20.0,
            'npts': 128,
            'delta': 0.15,
            'units': 1
        },
        'velocity_model': [
            {'thickness': 5000.0, 'vp': 5000.0, 'vs': 3000.0, 'rho': 2500.0, 'qp': 500.0, 'qs': 200.0},
            {'thickness': 15000.0, 'vp': 6500.0, 'vs': 3800.0, 'rho': 2800.0, 'qp': 1000.0, 'qs': 500.0}
        ],
        'moment_tensor': {'flag': 0},
        'stations': [
            {'name': 'STA1', 'latitude': -33.40, 'longitude': -70.66, 'height': 0.0, 'use_n': True, 'use_e': True, 'use_z': True},
        ],
        'inversion_params': [
            # Dummy parameters to satisfy config parser
        ],
        'inversion_process': {
            'Algorithm type': 0,
            'Number of iterations': 1,
            'Sample size for first iteration (SS1)': 4,
            'Sample size for other iterations': 4,
            'Cells to resample': 2
        }
    }
    return ConfigParser.from_dict(params)

def run_test():
    print("Iniciando prueba de NA Paralelo (Conda kdellipspy, Python 3.11)...")
    cfg = create_mock_config()
    
    # Mock inversion parameters manually
    from kdellipspy.core.config_parser import InversionParam, InversionParams
    cfg.inversion_params = InversionParams(parameters=[
        InversionParam(name="a1", min_val=1.0, max_val=10.0, flag=1),
        InversionParam(name="a2", min_val=1.0, max_val=10.0, flag=1),
        InversionParam(name="theta", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="np", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="tp", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="dmax", min_val=0.0001, max_val=0.005, flag=1),
        InversionParam(name="vr", min_val=1.0, max_val=4.0, flag=1),
    ])

    # Generate some mock observed data (zeros)
    nsta = len(cfg.stations.stations)
    npts = cfg.observed_data.npts
    obs_data = np.zeros((nsta, 3, npts))
    time_array = np.linspace(cfg.observed_data.t1, cfg.observed_data.t2, npts)
    
    # Dummy azi_times_array
    azi_times = np.array([[0.0, 5.0, 10.0]])

    axitra_dir = (Path(__file__).resolve().parent.parent / "kdellipspy" / "axitra" / "src").resolve()

    model = NAInversionModel(
        config=cfg,
        axitra_dir=str(axitra_dir),
        observed_waveforms=obs_data,
        time_array=time_array,
        azi_times_array=azi_times
    )

    print("\n--- Ejecutando NA con n_jobs=2 ---")
    na_cfg = NAConfig(
        n_samples_initial=4,
        n_samples_iteration=4,
        n_iterations=1,
        n_cells_resample=2,
        n_jobs=2
    )
    
    result = model.run_na_search(na_cfg)
    
    print(f"\n✓ Prueba completada.")
    print(f"Mejor misfit encontrado: {result.best_model.misfit:.6e}")
    print(f"Total de modelos evaluados: {len(result.all_models)}")

if __name__ == "__main__":
    run_test()
