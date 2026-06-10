import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --- Configuración ---
sys.path.insert(0, str(Path.cwd().parent))
import kdellipspy as kde

# --- 1. Cargar Configuración ---
INPUT_CTL_PATH = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl")
DATA_DIR = INPUT_CTL_PATH.parent / "SELECTED_DATA"
OUTPUT_DIR = INPUT_CTL_PATH.parent

print(f"Cargando configuración desde: {INPUT_CTL_PATH}")
cfg = kde.ConfigParser(filepath=INPUT_CTL_PATH)

# --- 2. Forzar una única subfalla en el hipocentro ---
# Al usar Nx=1, Ny=1, el modelo se convierte en una fuente puntual.
cfg.fault_plane.nx = 1
cfg.fault_plane.ny = 1
# Centramos el hipocentro en esta única subfalla
cfg.fault_plane.hx = cfg.fault_plane.lx / 2
cfg.fault_plane.hy = cfg.fault_plane.ly / 2
print(f"Configuración modificada para usar 1x1 subfallas (fuente puntual).")

# --- 3. Cargar Datos Observados ---
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

# --- 4. Generar Geometría Base y Funciones de Green ---
fm = kde.AxitraForwardModel.from_config(cfg)
# Esta geometría base contendrá los 6 puntos fuente (basis_slot 0 a 5) para nuestra única subfalla.
base_geom = fm.geometry_builder.build()

print(f"Geometría base generada con {base_geom.nsources} puntos fuente.")

print("Calculando funciones de Green una sola vez...")
ap = fm.build_axitra(base_geom)
fm.green(ap)

# --- 5. Calcular pesos y momento total del MT ---
# Usamos los métodos de EllipticalSlipMapper para obtener los pesos y el M0
# como si fuera un slip uniforme sobre una única falla.
slip_mapper = kde.EllipticalSlipMapper(cfg)
total_moment_nm = slip_mapper._mt_target_m0_nm()
mt_weights = slip_mapper._mt_component_weights_signed()
print(f"Momento Sísmico Total (M0): {total_moment_nm:.2e} Nm")
print(f"Pesos de las 6 fuentes base: {[f'{w:.3f}' for w in mt_weights]}")


# --- 6. Calcular Sintéticos para cada Fuente Elemental ---
elemental_synthetics = []
for k in range(6):
    print(f"Calculando sintéticos para la fuente base #{k+1}...")
    
    # Crea una geometría temporal que contiene SOLO UN punto fuente
    temp_geom = base_geom.copy()
    point_k = temp_geom.source_points[k]
    
    # Asignar el momento ponderado a este único punto
    # El slip se calcula a partir del momento para ser consistente con axitra
    mu = point_k.mu_pa
    area = point_k.width * point_k.length if (point_k.width > 0 and point_k.length > 0) else base_geom.subfaults[0].area_m2
    
    point_k.moment = total_moment_nm * mt_weights[k]
    if mu > 0 and area > 0:
        point_k.displacement = abs(point_k.moment) / (mu * area)

    temp_geom.source_points = [point_k]
    
    # Realizar convolución solo para este punto
    _, sx, sy, sz = fm.conv(ap, temp_geom, quiet=True)
    
    n_cut = min(len(time_array), sx.shape[1])
    synthetics_k = np.stack([sz[:, :n_cut], sx[:, :n_cut], sy[:, :n_cut]], axis=1) # Z, N, E
    
    # Filtrar
    synthetics_k_filtered = kde.bandpass_filter_waveforms(
        synthetics_k, time_array[:n_cut], float(cfg.ellipse.freq1), float(cfg.ellipse.freq2)
    )
    elemental_synthetics.append(synthetics_k_filtered)

# Sumar todos los sintéticos elementales para obtener el resultado final
total_synthetics = np.sum(elemental_synthetics, axis=0)

# --- 7. Graficar los Resultados ---

# Gráfico 1: Comparación de las 6 fuentes elementales en una estación
fig1, ax1 = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
comp_labels = ['Vertical (Z)', 'Norte (N)', 'Este (E)']
station_idx_to_plot = 0 # Graficar la primera estación (e.g., AC01)
station_name = cfg.stations.stations[station_idx_to_plot].name

for comp_idx in range(3):
    ax1[comp_idx].set_title(f"Componente {comp_labels[comp_idx]} en Estación {station_name}", fontsize=10)
    
    # Graficar cada fuente elemental
    for k in range(6):
        ax1[comp_idx].plot(time_array, elemental_synthetics[k][station_idx_to_plot, comp_idx, :], lw=0.8, alpha=0.7, label=f'Fuente Base {k+1}')
    
    # Graficar la suma
    ax1[comp_idx].plot(time_array, total_synthetics[station_idx_to_plot, comp_idx, :], 'k-', lw=1.5, label='Suma Total')
    ax1[comp_idx].grid(True, linestyle=':', alpha=0.6)
    ax1[comp_idx].set_ylabel("Amplitud (m/s)")

ax1[0].legend(ncol=4, fontsize=8, loc='upper right')
ax1[-1].set_xlabel("Tiempo (s)")
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig1.suptitle("Descomposición del Tensor de Momento en Fuentes Elementales", fontsize=14)
output_fig1 = OUTPUT_DIR / "elemental_sources_comparison.png"
plt.savefig(output_fig1, dpi=150)
print(f"\\nGráfico de fuentes elementales guardado en: {output_fig1}")


# Gráfico 2: Comparación de la suma total vs. datos observados
fig2, axes = plt.subplots(n_stations, 3, figsize=(14, 2.5 * n_stations), sharex=True, constrained_layout=True)
if n_stations == 1: axes = axes[np.newaxis, :]

fig2.suptitle(f"Comparación: Suma de Sintéticos (mt_strict) vs. Observado\nFiltro: {cfg.ellipse.freq1}-{cfg.ellipse.freq2} Hz", fontsize=14)

for i, st in enumerate(cfg.stations.stations):
    row_max = max(np.abs(observed_waveforms[i]).max(), np.abs(total_synthetics[i]).max()) or 1.0
    for j in range(3):
        ax = axes[i, j]
        ax.plot(time_array, observed_waveforms[i, j, :] / row_max, 'k-', lw=1.2, label='Observado')
        ax.plot(time_array, total_synthetics[i, j, :] / row_max, 'r--', lw=1.0, label='Suma Sintéticos')
        if i == 0: ax.set_title(comp_labels[j])
        if j == 0: ax.set_ylabel(st.name, rotation=0, labelpad=30, va='center', ha='right')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.set_yticks([])

axes[0, 0].legend(loc='upper right', fontsize=8)
fig2.supxlabel("Tiempo (s)")

output_fig2 = OUTPUT_DIR / "total_synthetics_vs_observed.png"
plt.savefig(output_fig2, dpi=150)
print(f"Gráfico de comparación total guardado en: {output_fig2}")


