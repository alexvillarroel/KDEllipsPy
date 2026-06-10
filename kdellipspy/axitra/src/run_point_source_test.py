import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --- Configuración ---
sys.path.insert(0, str(Path.cwd().parent))
import kdellipspy as kde

# --- 1. Cargar Configuración y Datos Observados ---
INPUT_CTL_PATH = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl")
DATA_DIR = INPUT_CTL_PATH.parent / "SELECTED_DATA"
OUTPUT_IMAGE = INPUT_CTL_PATH.parent / "point_source_comparison.png"

print(f"Cargando configuración desde: {INPUT_CTL_PATH}")
cfg = kde.ConfigParser(filepath=INPUT_CTL_PATH)

print(f"Cargando datos observados desde: {DATA_DIR}")
try:
    observed_waveforms, time_array = kde.load_and_filter_observed_data(
        config=cfg,
        data_dir=DATA_DIR,
        prefer_raw=False
    )
except FileNotFoundError:
    print(f"Error: No se encontraron los archivos 'real_vel_x/y/z' en {DATA_DIR}.")
    sys.exit(1)

# --- 2. Forzar una única subfalla en el hipocentro ---
cfg.fault_plane.nx = 1
cfg.fault_plane.ny = 1
cfg.fault_plane.hx = cfg.fault_plane.lx / 2
cfg.fault_plane.hy = cfg.fault_plane.ly / 2
print(f"Configuración modificada para usar 1x1 subfallas (fuente puntual).")

# Forzar el modo MT_STRICT
cfg.moment_tensor.scaling_mode = "mt_strict"
print(f"Modo de Tensor de Momento: {cfg.moment_tensor.scaling_mode}")

# Modelo de elipse que solo cubre la subfalla central
point_source_model = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.5])

# --- 3. Generar Geometría y Sintéticos ---
fm = kde.AxitraForwardModel.from_config(cfg)
# build_geometry_with_ellipse_slip crea los 6 puntos fuente en el mismo lugar y les asigna momento
geom = fm.build_geometry_with_ellipse_slip(point_source_model)

print(f"Geometría generada con {geom.nsources} puntos fuente (6 por cada subfalla).")
print(f"Momento Sísmico Total (forzado por mt_strict): {geom.total_moment_nm():.2e} Nm")

print("Calculando funciones de Green...")
ap = fm.build_axitra(geom)
fm.green(ap)

print("Realizando convolución (esto suma la contribución de las 6 fuentes)...")
time_syn, sx, sy, sz = fm.conv(ap, geom, quiet=True)

# --- 4. Procesar y Filtrar Sintéticos ---
n_cut = min(len(time_array), sx.shape[1])
synthetics = np.stack([sz[:, :n_cut], sx[:, :n_cut], sy[:, :n_cut]], axis=1) # Z, N, E

print("Filtrando sintéticos...")
synthetics = kde.bandpass_filter_waveforms(
    synthetics,
    time_array[:n_cut],
    freq1=float(cfg.ellipse.freq1),
    freq2=float(cfg.ellipse.freq2)
)

# --- 5. Graficar la Comparación ---
comp_labels = ['Vertical (Z)', 'Norte (N)', 'Este (E)']
n_stations = len(cfg.stations.stations)
fig, axes = plt.subplots(n_stations, 3, figsize=(14, 2.5 * n_stations), sharex=True, constrained_layout=True)
if n_stations == 1: axes = axes[np.newaxis, :]

fig.suptitle(f"Comparación Forward: Fuente Puntual (mt_strict) vs. Observado\nFiltro: {cfg.ellipse.freq1}-{cfg.ellipse.freq2} Hz", fontsize=14)

for i, st in enumerate(cfg.stations.stations):
    row_max = max(np.abs(observed_waveforms[i]).max(), np.abs(synthetics[i]).max()) or 1.0
    for j in range(3):
        ax = axes[i, j]
        ax.plot(time_array[:n_cut], observed_waveforms[i, j, :n_cut] / row_max, 'k-', lw=1.2, label='Observado')
        ax.plot(time_array[:n_cut], synthetics[i, j] / row_max, 'r--', lw=1.0, label='Sintético')
        ax.grid(True, linestyle=':', alpha=0.6)
        if i == 0: ax.set_title(comp_labels[j])
        if j == 0: ax.set_ylabel(st.name, rotation=0, labelpad=30, va='center', ha='right')
        ax.set_yticks([])

axes[0, 0].legend(loc='upper right', fontsize=8)
fig.supxlabel("Tiempo (s)")

plt.savefig(OUTPUT_IMAGE, dpi=150)
print(f"\\nGráfico guardado en: {OUTPUT_IMAGE}")

