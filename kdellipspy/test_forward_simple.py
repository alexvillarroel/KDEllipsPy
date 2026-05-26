
import numpy as np
import os
import sys
from pathlib import Path

# --- Bloque de resolución de rutas ---
src_root = Path(__file__).resolve().parent.parent
if str(src_root) not in sys.path:
    sys.path.insert(0, str(src_root))

from kdellipspy import ConfigParser, AxitraForwardModel, GraphicsSuite

def create_mock_config():
    """Crea un diccionario de configuración minimal para la prueba."""
    params = {
        'source_position': {
            'event_name': 'TestEvent',
            'latitude': -33.45,
            'longitude': -70.66,
            'depth': 10.0, # KM
            'strike': 0.0,
            'dip': 90.0,
            'rake': 0.0
        },
        'fault_plane': {
            'Lx': 20000.0, # m
            'Ly': 20000.0, # m
            'nx': 2,
            'ny': 2,
            'hx': 10000.0, # m
            'hy': 10000.0  # m
        },
        'ellipse': {
            'n_ellipses': 1,
            'initial_slip': 0,
            'slip_shape': 1,
            'freq1': 0.1,
            'freq2': 1.0,
            't0': 2.0
        },
        'observed_data': {
            't1': -5.0,
            't2': 30.0,
            'npts': 256,
            'delta': 0.1,
            'units': 1
        },
        'velocity_model': [
            {'thickness': 5.0, 'vp': 5.0, 'vs': 3.0, 'rho': 2.5e-6, 'qp': 500.0, 'qs': 200.0},
            {'thickness': 15.0, 'vp': 6.5, 'vs': 3.8, 'rho': 2.8e-6, 'qp': 1000.0, 'qs': 500.0}
        ],
        'moment_tensor': {'flag': 0},
        'stations': [
            {'name': 'STA1', 'latitude': -33.40, 'longitude': -70.66, 'height': 0.0, 'use_n': True, 'use_e': True, 'use_z': True},
            {'name': 'STA2', 'latitude': -33.45, 'longitude': -70.60, 'height': 0.0, 'use_n': True, 'use_e': True, 'use_z': True}
        ],
        'inversion_params': [] # No necesario para forward simple
    }
    return ConfigParser.from_dict(params)

def run_test():
    print("Iniciando prueba de Forward con unidades consistentes (KM y KG/KM3)...")
    cfg = create_mock_config()
    
    # Ruta al binario de axitra
    # Intentamos encontrarlo relativo al paquete si no se especifica
    src_root = Path(__file__).resolve().parent
    axitra_dir = (src_root / "axitra").resolve()
    
    if not axitra_dir.exists():
        # Fallback al directorio actual
        axitra_dir = Path("./axitra").resolve()
    if not (axitra_dir / "axitra").exists():
        print(f"Error: No se encontró el binario en {axitra_dir}")
        return

    # 1. Construir modelo Forward
    fm = AxitraForwardModel.from_config(cfg, axitra_dir=str(axitra_dir))
    
    # 2. Construir Geometría
    print("Construyendo geometría...")
    geom = fm.build_geometry()
    
    # 3. Aplicar un modelo de elipse manual
    # dmax = 1e-3 (1 km, para ver si sale algo grande)
    model_params = np.array([5.0, 3.0, 0.0, 0.5, 0.0, 0.001, 2.5])
    geom_with_slip = fm.apply_ellipse_model_to_geometry(geom, model_params)
    
    print(f"Número de fuentes generadas: {len(geom_with_slip.source_points)}")

    # 4. Correr axitra
    print("Ejecutando axitra...")
    ap = fm.build_axitra(geom_with_slip, latlon=False) 
    ap = fm.green(ap, quiet=False)
    
    # 5. Convolución
    print("Calculando convolución...")
    _, sx, sy, sz = fm.conv(ap, geom_with_slip, source_type=5, t0=cfg.ellipse.t0)
    
    # Reordenar a (nsta, 3, npts)
    synthetics = np.array([sx, sy, sz])
    synthetics = np.transpose(synthetics, (1, 0, 2))
    
    max_amp = np.max(np.abs(synthetics))
    print(f"Amplitud máxima detectada: {max_amp:.6e}")
    
    # 6. Visualización básica
    print("Generando gráficos de prueba...")
    gs = GraphicsSuite(base_dir=".", cfg=cfg, show=False)
    
    npts = synthetics.shape[2]
    dt = float(cfg.observed_data.delta)
    time_array = np.arange(npts) * dt + float(cfg.observed_data.t1)
    
    station_names = [s.name for s in cfg.stations.stations]
    gs.plot_synthetic_components(time_array, synthetics, station_names)
    
    print("\n✓ Prueba completada. Revisa 'Figures/Synthetic_Seismograms.png'.")
    
    # Limpieza
    ap.clean()

if __name__ == "__main__":
    run_test()
