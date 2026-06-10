import numpy as np
from kdellipspy import AxitraForwardModel

def test_forward_model_execution(mock_config, axitra_dir):
    """
    Tests the fundamental physics engine execution (Axitra integration).
    Ensures geometry is generated, compiled subroutines run, and convolution produces valid shapes.
    """
    fm = AxitraForwardModel.from_config(mock_config, axitra_dir=axitra_dir)
    
    # 1. Build Base Geometry
    geom = fm.build_geometry()
    
    # 2. Apply Ellipse Model
    # model = [a1, a2, theta, np, tp, dmax, vr]
    model_params = np.array([12.0, 8.0, 0.2, 0.4, 0.3, 3.5, 2.8])
    geom_with_slip = fm.apply_ellipse_model_to_geometry(geom, model_params)
    
    assert len(geom_with_slip.source_points) > 0, "No source points generated in the fault plane"

    # 3. Build Axitra Input and compute Green's functions
    ap = fm.build_axitra(geom_with_slip, latlon=False) 
    ap = fm.green(ap, quiet=True)
    
    # 4. Convolution
    _, sx, sy, sz = fm.conv(ap, geom_with_slip, source_type=1, unit=1, t0=mock_config.ellipse.t0, quiet=True)
    
    synthetics = np.array([sx, sy, sz])
    
    # Check dimensions: (n_components, n_stations, n_points)
    n_sta = len(mock_config.stations.stations)
    n_comp = 3
    n_pts = mock_config.observed_data.npts
    
    assert synthetics.shape == (n_comp, n_sta, n_pts)
    assert np.any(synthetics != 0.0), "Synthetics are all zeros, Axitra might have failed silently"
    assert not np.any(np.isnan(synthetics)), "Synthetics contain NaNs"
    
    ap.clean()
