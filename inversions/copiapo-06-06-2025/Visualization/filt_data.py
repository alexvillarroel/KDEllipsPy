# %% esto es un python interactivo
import matplotlib.pyplot as plt
import obspy
from pathlib import Path

plt.style.use("ggplot")
# Use path relative to this script's parent directory or CWD
base_dir = Path(__file__).resolve().parent.parent
data_folder = str(base_dir / "RAW" / "*.SAC")
data = obspy.read(data_folder)
# Visualize the data
data[0].plot()
