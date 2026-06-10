# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: kdellipspy
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Inversión Cinemática: Terremoto de Copiapó (06-06-2025)
#
# Este notebook realiza la inversión utilizando los parametros configurandolos desde el mismo codigo.

# %% Liberias
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime

import kdellipspy as kde

# %%
cfg = kde.ConfigParser(
    filepath=Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl")
)
data_dir = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/DATA_ISOLA")
print([i for i in dir(cfg) if not i.startswith("_")])
# %% [markdown]
# ## 1. Visualización Previa de SAC Files
#

print(cfg.ellipse.source_type)
cfg.plot_stations()
# %% [markdown]
# ## 4. Inversión Neighbourhood Algorithm (NA)
observed_waveforms, time_array = kde.load_and_filter_observed_data(
    input_ctl_path=Path(
        "/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl"
    ),
    data_dir=data_dir,
    prefer_raw=False,
    freq1=0.10,
    freq2=0.2,
)

azi_times_array = kde.build_azi_times_array(config=cfg)

model = kde.NAInversionModel(
    config=cfg,
    observed_waveforms=observed_waveforms,
    time_array=time_array,
    azi_times_array=azi_times_array,
)

# Desactivar la paralelización forzando n_jobs=1
cfg.inversion_process.n_jobs = 1

result = model.run_na_search()

# %%
os.makedirs("./output", exist_ok=True)
result.save("./output/inversion_result.joblib")
result.plot(show=True, save_path="./output/na_results.png")
result.plot_convergence(show=True, save_path="./output/parameter_convergence.png")
result.plot_fit(show=True, save_path="./output/waveform_fit.png")
result.plot_elipse(show=True, save_path="./output/ellipse_fit.png")
