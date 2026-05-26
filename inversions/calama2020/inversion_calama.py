# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
# ---

# %% [markdown]
# # Inversión Cinemática: Terremoto de Copiapó (06-06-2025)
#
# Este notebook realiza la inversión utilizando los parámetros optimizados de AXITRA y el pre-procesador de datos con control de tiempo y padding.

# %%
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime


def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "kdellipspy").exists():
            return p
    raise FileNotFoundError("No root found.")


root = find_project_root(Path(os.getcwd()))
sys.path.append(str(root))
import kdellipspy as kde

print(f"Project root: {root}")

# %%
input_ctl = root / "inversions" / "calama2020" / "input.ctl"
cfg = kde.ConfigParser(filepath=str(input_ctl))
print(
    f"Evento: {cfg.source_position.event_name} | Estaciones originales: {len(cfg.stations.stations)}"
)
print("\n=== Configuración Completa ===")
print(f"Nombre del evento: {cfg.source_position.event_name}")
print(f"Tiempo de origen: {cfg.source_position.origin_time}")
print(f"Profundidad: {cfg.source_position.depth}")
print(f"Latitud: {cfg.source_position.latitude}")
print(f"Longitud: {cfg.source_position.longitude}")
print(f"\nDatos observados:")
print(f"  - Número de puntos (npts): {cfg.observed_data.npts}")
print(f"  - Delta (dt): {cfg.observed_data.delta}")
print(f"  - Duración: {cfg.observed_data.npts * cfg.observed_data.delta} s")
print(f"\nEstaciones ({len(cfg.stations.stations)}):")
for s in cfg.stations.stations:
    print(f"  - {s.name}: Lat {s.latitude}, Lon {s.longitude}, Alt {s.height}")

# %% [markdown]
# # 2. Visualizacion de parametros de inversion y de fuente
#


# %%
def mostrar_objeto(obj, nivel=0, max_nivel=2):
    pad = "  " * nivel

    if isinstance(obj, (str, int, float, bool, type(None), np.number)):
        print(f"{pad}{obj!r}")
        return

    if isinstance(obj, (list, tuple)):
        for i, item in enumerate(obj):
            print(f"{pad}[{i}]")
            mostrar_objeto(item, nivel + 1, max_nivel)
        return

    if hasattr(obj, "__dict__") and max_nivel >= 0:
        for k, v in vars(obj).items():
            if k.startswith("_"):
                continue
            if hasattr(v, "__dict__"):
                print(f"{pad}{k}:")
                mostrar_objeto(v, nivel + 1, max_nivel - 1)
            elif isinstance(v, (list, tuple)):
                print(f"{pad}{k}:")
                mostrar_objeto(v, nivel + 1, max_nivel - 1)
            else:
                print(f"{pad}{k}: {v!r}")
    else:
        print(f"{pad}{repr(obj)}")


mostrar_objeto(cfg)

# %% [markdown]
# ## 4. Inversión Neighbourhood Algorithm (NA)

# %%
axitra_dir = root / "kdellipspy" / "axitra" / "src"
output_dir = root / "inversions" / "calama2020" / "DATA"
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
    axitra_ikmax=100000,
)

result = inversion.run_na_search()

# %% [markdown]
# PLOT RESULTS

# %%
graphics = kde.core.GraphicsSuite(
    base_dir="/home/alex/KDEllipsPy/inversions/calama2020/Figures", show=True
)
graphics.plot_na_results(result)
