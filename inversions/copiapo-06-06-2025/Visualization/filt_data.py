# %% esto es un python interactivo
import matplotlib.pyplot as plt
import obspy

plt.style.use("ggplot")
data_folder = "/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/RAW/*.SAC"
data = obspy.read(data_folder)
# Visualize the data
data[0].plot()
