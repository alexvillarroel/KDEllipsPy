import pytest
import numpy as np
from pathlib import Path
import sys

# Ensure the kdellipspy package is importable from the root directory
src_root = Path(__file__).resolve().parent.parent
if str(src_root) not in sys.path:
    sys.path.insert(0, str(src_root))

from kdellipspy import ConfigParser
from kdellipspy.core.config_parser import InversionParam, InversionParams

@pytest.fixture
def axitra_dir():
    """Returns the absolute path to the compiled axitra binaries."""
    return str((Path(__file__).resolve().parent.parent / "kdellipspy" / "axitra" / "src").resolve())

@pytest.fixture
def mock_config():
    """Returns a fully populated ConfigParser instance suitable for testing without external files."""
    params = {
        'source_position': {
            'Event Name': 'TestEvent',
            'Latitude': -33.45,
            'Longitude': -70.66,
            'Depth': 10.0,
            'Strike': 45.0,
            'Dip': 45.0,
            'Rake': 90.0
        },
        'fault_plane': {
            'Length along strike (Lx)': 20000.0,
            'Length along dip (Ly)': 20000.0,
            'Number of subfaults along strike (Nx)': 20,
            'Number of subfaults along dip (Ny)': 20,
            'Hypocenter position strike (Hx)': 10000.0,
            'Hypocenter position dip (Hy)': 10000.0
        },
        'ellipse': {
            'Number of ellipses': 1,
            'Initial slip': 0,
            'Slip shape': 1,
            'Frequency 1 (Freq1)': 0.1,
            'Frequency 2 (Freq2)': 1.0,
            'Time shift (T0)': 2.0
        },
        'observed_data': {
            'Time window start (t1)': 0.0,
            'Time window end (t2)': 25.5,
            'Number of points (Npts)': 256,
            'Delta / Time step': 0.1,
            'Units': 1
        },
        'velocity_model': [
            {'thickness': 5000.0, 'vp': 5210.0, 'vs': 2990.0, 'rho': 2500.0, 'qp': 500.0, 'qs': 200.0},
            {'thickness': 15000.0, 'vp': 6500.0, 'vs': 3800.0, 'rho': 2800.0, 'qp': 1000.0, 'qs': 500.0}
        ],
        'moment_tensor': {'Moment Tensor Flag': 0},
        'stations': [
            {'name': 'STA1', 'latitude': -33.40, 'longitude': -70.66, 'height': 0.0, 'use_n': True, 'use_e': True, 'use_z': True},
        ],
        'inversion_params': [],
        'inversion_process': {
            'Algorithm type': 0,
            'Number of iterations': 1,
            'Sample size for first iteration (SS1)': 4,
            'Sample size for other iterations': 4,
            'Cells to resample': 2,
            'Parallel workers (n_jobs)': 1
        }
    }
    cfg = ConfigParser.from_dict(params)
    
    # Inject inversion parameters
    cfg.inversion_params = InversionParams(parameters=[
        InversionParam(name="a1", min_val=1.0, max_val=10.0, flag=1),
        InversionParam(name="a2", min_val=1.0, max_val=10.0, flag=1),
        InversionParam(name="theta", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="np", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="tp", min_val=0.0, max_val=1.0, flag=1),
        InversionParam(name="dmax", min_val=0.0001, max_val=0.005, flag=1),
        InversionParam(name="vr", min_val=1.0, max_val=4.0, flag=1),
    ])
    return cfg
