import sys
from pathlib import Path
import os
import numpy as np

# Simulate what the notebook does
root = Path(os.getcwd()).parent
sys.path.append(str(root))
import kdellipspy as kde

input_ctl = root / 'inversions' / 'copiapo-06-06-2025' / 'input.ctl'
cfg = kde.ConfigParser(filepath=str(input_ctl))

axitra_dir = root / 'kdellipspy' / 'axitra'
output_dir = root / 'inversions' / 'copiapo-06-06-2025' / 'DATA_PROCESSED'

observed_waveforms, time_array = kde.load_and_filter_observed_data(
    config=cfg,
    data_dir=str(output_dir),
    prefer_raw=False,
)
azi_times_array = kde.build_azi_times_array(config=cfg)

inversion = kde.NAInversionModel(
    config=cfg,
    axitra_dir=str(axitra_dir),
    observed_waveforms=observed_waveforms,
    time_array=time_array,
    azi_times_array=azi_times_array,
    axitra_aw=0.5,
    axitra_ikmax=100000
)

result = inversion.run_na_search()
print("Result best misfit:", result.best_model.misfit if result.best_model else "None")
