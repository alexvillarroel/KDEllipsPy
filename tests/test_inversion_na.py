import numpy as np
from kdellipspy import NAInversionModel, NAConfig

def test_na_inversion_parallel(mock_config, axitra_dir):
    """
    Tests the Neighbourhood Algorithm (NA) search process using multiple CPU cores.
    Verifies that the parallel NASearcher yields the correct number of models and handles misfit logic.
    """
    # Create empty observed data matching the mock_config dimensions
    nsta = len(mock_config.stations.stations)
    npts = mock_config.observed_data.npts
    obs_data = np.zeros((nsta, 3, npts))
    
    # Time array mapping the observed data
    time_array = np.linspace(mock_config.observed_data.t1, mock_config.observed_data.t2, npts)
    
    # Dummy azimuth/travel-time table [azimuth, P-time, S-time]
    azi_times = np.array([[0.0, 5.0, 10.0]])
    
    # Initialize the NA model
    model = NAInversionModel(
        config=mock_config,
        axitra_dir=axitra_dir,
        observed_waveforms=obs_data,
        time_array=time_array,
        azi_times_array=azi_times
    )
    
    # Configure a tiny NA search with 2 jobs (Parallel)
    # 4 initial + (4 per iteration * 1 iteration) = 8 total models evaluated
    na_cfg = NAConfig(
        n_samples_initial=4,
        n_samples_iteration=4,
        n_iterations=1,
        n_cells_resample=2,
        n_jobs=2,
        random_seed=42
    )
    
    result = model.run_na_search(na_config=na_cfg)
    
    # Assertions
    assert len(result.all_models) == 8, f"Expected 8 models evaluated, got {len(result.all_models)}"
    assert result.best_model is not None, "Best model was not recorded"
    
    # Since observed is all zeros, the misfit should reflect the sum of squares of synthetics
    # As long as it is a valid float and not infinite, the objective function executed properly.
    assert np.isfinite(result.best_model.misfit), "Misfit is not a valid finite number"
