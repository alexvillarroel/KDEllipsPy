import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --- Configuración ---
sys.path.insert(0, str(Path.cwd().parent))
import kdellipspy as kde

# --- 1. Cargar Configuración y Datos Observados ---
INPUT_CTL_PATH = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl")
AXITRA_SRC_DIR = Path("/home/alex/KDEllipsPy/kdellipspy/axitra/src")
DATA_DIR = INPUT_CTL_PATH.parent / "SELECTED_DATA"
OUTPUT_IMAGE = INPUT_CTL_PATH.parent / "fortran_direct_comparison.png"
AXITRA_ID = "axi_test" # Usaremos un ID único para estos archivos

print(f"Directorio de trabajo de Axitra: {AXITRA_SRC_DIR}")
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

# --- 2. Forzar Fuente Puntual y MT_STRICT ---
cfg.fault_plane.nx = 1
cfg.fault_plane.ny = 1
cfg.fault_plane.hx = cfg.fault_plane.lx / 2
cfg.fault_plane.hy = cfg.fault_plane.ly / 2
cfg.moment_tensor.scaling_mode = "mt_strict"
print("Configuración: 1x1 subfalla (fuente puntual) con modo 'mt_strict'.")

# --- 3. Generar Geometría y Archivos de Entrada para Axitra ---
fm = kde.AxitraForwardModel.from_config(cfg)
point_source_model = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.5])
geom = fm.build_geometry_with_ellipse_slip(point_source_model)

print("Escribiendo archivos de entrada para Axitra:")
# .stat (Corregido: usar cfg.stations)
stat_path = AXITRA_SRC_DIR / f"{AXITRA_ID}.stat"
np.savetxt(stat_path, cfg.stations.to_axitra_stations(latlon=False), fmt="%d %.3f %.3f %.3f")
print(f"  ✓ {stat_path.name}")
# .source
source_path = AXITRA_SRC_DIR / f"{AXITRA_ID}.source"
np.savetxt(source_path, geom.to_axitra_sources(latlon=False), fmt="%d %.3f %.3f %.3f")
print(f"  ✓ {source_path.name}")
# .hist
hist_path = AXITRA_SRC_DIR / f"{AXITRA_ID}.hist"
np.savetxt(hist_path, geom.to_axitra_hist(), fmt=("%d", "%.6e") + ("%.3f",) * 6)
print(f"  ✓ {hist_path.name}")
# .data
data_content = fm.get_axitra_data_content(source_filename=source_path.name, stat_filename=stat_path.name)
data_path = AXITRA_SRC_DIR / f"{AXITRA_ID}.data"
with open(data_path, "w") as f:
    f.write(data_content)
print(f"  ✓ {data_path.name}")

# --- 4. Ejecutar Binarios de Fortran ---
import subprocess
import importlib.util

# a) axitra
print("\\nEjecutando el binario 'axitra' de Fortran...")
axitra_executable = AXITRA_SRC_DIR / "axitra"
if not axitra_executable.exists():
    print(f"Error: No se encontró el ejecutable '{axitra_executable}'.")
    sys.exit(1)
process = subprocess.run([str(axitra_executable), AXITRA_ID.split('_')[1]], cwd=AXITRA_SRC_DIR, capture_output=True, text=True)
if process.returncode != 0:
    print("--- Salida de Axitra (stdout) ---\\n", process.stdout)
    print("--- Salida de Axitra (stderr) ---\\n", process.stderr)
    sys.exit("¡La ejecución de 'axitra' falló!")
print("  ✓ 'axitra' completado.")

# b) convm
print("Ejecutando la convolución via wrapper 'convmPy'...")
spec = importlib.util.spec_from_file_location("convmPy", AXITRA_SRC_DIR / "convmPy.cpython-311-x86_64-linux-gnu.so")
convmPy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(convmPy)

id_val = int(AXITRA_ID.split('_')[1])
sx, sy, sz = convmPy.moment_conv(id_val, 5, float(cfg.ellipse.t0), 0.0, int(cfg.observed_data.units), n_sta=len(cfg.stations.stations), n_pts=int(cfg.observed_data.npts))
print("  ✓ Convolución completada.")

# --- 5. Cargar datos y graficar ---
synthetics = np.stack([sz.T, sx.T, sy.T], axis=1)
n_cut = min(len(time_array), synthetics.shape[2])

print("Filtrando sintéticos...")
synthetics_filtered = kde.bandpass_filter_waveforms(
    synthetics[:, :, :n_cut], time_array[:n_cut], float(cfg.ellipse.freq1), float(cfg.ellipse.freq2)
)

comp_labels = ['Vertical (Z)', 'Norte (N)', 'Este (E)']
fig, axes = plt.subplots(len(cfg.stations.stations), 3, figsize=(14, 2.5 * len(cfg.stations.stations)), sharex=True, constrained_layout=True)
if len(cfg.stations.stations) == 1: axes = axes[np.newaxis, :]
fig.suptitle(f"Comparación (Fortran Directo): Punto Fuente vs. Observado\nFiltro: {cfg.ellipse.freq1}-{cfg.ellipse.freq2} Hz", fontsize=14)
for i, st in enumerate(cfg.stations.stations):
    row_max = max(np.abs(observed_waveforms[i]).max(), np.abs(synthetics_filtered[i]).max()) or 1.0
    for j in range(3):
        ax = axes[i, j]
        ax.plot(time_array[:n_cut], observed_waveforms[i, j, :n_cut] / row_max, 'k-', lw=1.2, label='Observado')
        ax.plot(time_array[:n_cut], synthetics_filtered[i, j, :] / row_max, 'r--', lw=1.0, label='Sintético (Fortran)')
        ax.grid(True, linestyle=':', alpha=0.6)
        if i == 0: ax.set_title(comp_labels[j])
        if j == 0: ax.set_ylabel(st.name, rotation=0, labelpad=30, va='center', ha='right')
        ax.set_yticks([])

axes[0, 0].legend(loc='upper right', fontsize=8)
fig.supxlabel("Tiempo (s)")
plt.savefig(OUTPUT_IMAGE, dpi=150)
print(f"\\nGráfico guardado en: {OUTPUT_IMAGE}")

