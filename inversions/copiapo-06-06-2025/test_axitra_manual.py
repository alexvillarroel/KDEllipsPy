# %%
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import obspy

REPO_ROOT = Path("/home/alex/KDEllipsPy")
AXITRA_SRC = REPO_ROOT / "kdellipspy" / "axitra" / "src"

# Add repo root and axitra src so we can import the wrapper directly
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(AXITRA_SRC))

from axitra import Axitra, moment

DATA_DIR = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/DATA_ISOLA")

# 1. Manually set velocity model (from input.ctl)
# [thickness/depth, vp, vs, rho, qp, qs]
model = np.array(
    [
        [0.0, 4700.0, 2750.0, 2562.0, 300.0, 300.0],
        [5000.0, 5730.0, 3300.0, 2692.0, 300.0, 300.0],
        [10000.0, 6320.0, 3620.0, 2817.0, 300.0, 300.0],
        [15000.0, 6680.0, 3830.0, 2946.0, 300.0, 300.0],
        [20000.0, 6900.0, 3950.0, 3026.0, 300.0, 300.0],
        [25000.0, 7020.0, 4020.0, 3069.0, 300.0, 300.0],
        [30000.0, 7100.0, 4060.0, 3098.0, 300.0, 300.0],
        [35000.0, 7250.0, 4140.0, 3152.0, 300.0, 300.0],
        [40000.0, 7480.0, 4260.0, 3235.0, 300.0, 300.0],
        [45000.0, 7790.0, 4430.0, 3346.0, 300.0, 300.0],
        [50000.0, 8070.0, 4570.0, 3447.0, 1000.0, 1000.0],
        [55000.0, 8230.0, 4640.0, 3505.0, 1000.0, 1000.0],
        [60000.0, 8360.0, 4680.0, 3551.0, 1000.0, 1000.0],
        [65000.0, 8480.0, 4730.0, 3595.0, 1000.0, 1000.0],
        [70000.0, 8550.0, 4770.0, 3620.0, 1000.0, 1000.0],
        [75000.0, 8570.0, 4780.0, 3627.0, 1000.0, 1000.0],
        [80000.0, 8550.0, 4760.0, 3620.0, 1000.0, 1000.0],
        [85000.0, 8500.0, 4740.0, 3602.0, 1000.0, 1000.0],
        [90000.0, 8440.0, 4720.0, 3580.0, 1000.0, 1000.0],
        [95000.0, 8380.0, 4690.0, 3559.0, 1000.0, 10000.0],
    ],
    dtype=float,
)

# 2. Manually set stations
# [index, lat, lon, z]
#
# -27.5030 -70.8870 0.0000 A18C
# -27.8300 -70.1080 0.0000 A19C
# -23.6100 -70.2620 0.0000 A24F
# -29.1020 -71.4100 0.0000 A32C
# -26.1480 -70.5990 0.0000 AC01
# -28.8360 -70.2740 0.0000 AC05
# -23.9050 -69.2910 0.0000 PB19
stations = np.array(
    [
        [1, -27.5030, -70.8870, 0.0],
        [2, -27.8300, -70.1080, 0.0],
        [3, -23.6100, -70.2620, 0.0],
        [4, -29.1020, -71.4100, 0.0],
        [5, -26.1480, -70.5990, 0.0],
        [6, -28.8360, -70.2740, 0.0],
        [7, -23.9050, -69.2910, 0.0],
    ],
    dtype=float,
)
n_stations = stations.shape[0]
name_stations = [
    "A18C",
    "A19C",
    "A24F",
    "A32C",
    "AC01",
    "AC05",
    "PB19",
]
# 3. Manually set sources
# Only 1 point source at hypocenter
# Lat = -26.639, Lon = -70.404, Depth = 75.0 km
sources = np.array(
    [
        [1, -26.639, -70.404, 75000.0],
        [2, -26.639, -70.404, 75000.0],
        [3, -26.639, -70.404, 75000.0],
        [4, -26.639, -70.404, 75000.0],
        [5, -26.639, -70.404, 75000.0],
        [6, -26.639, -70.404, 75000.0],
    ],
    dtype=float,
)


# 4. Other axitra parameters
fmax = 2  # Need to be large enough for 0.09 Hz
duration = 204.8  # From config
dt = 0.2
npts = int(duration / dt)

# 5. Initialize Axitra manually
# axpath uses the directory where the binaries reside
axpath = str(AXITRA_SRC.resolve())

ap = Axitra(
    model=model,
    stations=stations,
    sources=sources,
    fmax=fmax,
    duration=duration,
    latlon=True,  # MANUALLY SET LATLON
    axpath=axpath,
    id=123,
)

# 6. Run Green's function
print("Computing Green's functions...")
ap = moment.green(ap)


def mt_decomposition(
    mrr, mtt, mpp, mrt, mrp, mtp, iexp, basis_system="NED", axitra_output="."
):
    def vec(mat):
        return np.array(
            [mat[0, 0], mat[1, 1], mat[2, 2], mat[0, 1], mat[0, 2], mat[1, 2]]
        )

    M = np.array([[mrr, mrt, mrp], [mrt, mtt, mtp], [mrp, mtp, mpp]]) * 10**iexp

    m1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    m2 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    m3 = np.array([[0, 0, 0], [0, 0, -1], [0, -1, 0]])
    m4 = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]])
    m5 = np.array([[0, 0, 0], [0, -1, 0], [0, 0, 1]])
    m6 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    basis = [m1, m2, m3, m4, m5, m6]
    M_vec = vec(M)
    B = np.array([vec(b) for b in basis]).T

    coeffs = np.linalg.solve(B, M_vec)
    print("Coeficients of your base: ", coeffs)
    M_rec = sum(c * b for c, b in zip(coeffs, basis))
    print("\nreconstructured Tensor: \n", M_rec)

    M0_components = [
        abs(c) / np.sqrt(2) * np.linalg.norm(b, "fro") for c, b in zip(coeffs, basis)
    ]

    M0_total = (1 / np.sqrt(2)) * np.linalg.norm(
        sum(c * b for c, b in zip(coeffs, basis)), "fro"
    )
    print(f"M0_total: {M0_total:.3e}")
    Mw = (2 / 3) * (np.log10(M0_total) - 9.1)
    print(f"Mw: {Mw:.2f}")

    return coeffs, M0_components


# 7. Moment tensor 6 orthogonal bases
# From input.ctl:
mrr = -2.486000
mtt = 7.932000
mpp = 0.083000
mrt = 0.038000
mrp = -0.739000
mtp = 0.473000
iexp = 18.0

coeffs, M0_components = mt_decomposition(mrr, mtt, mpp, mrt, mrp, mtp, iexp)

b_str = [0.0, 270.0, 0.0, 90.0, 0.0, 0.0]
b_dip = [90.0, 90.0, 90.0, 45.0, 45.0, 0.0]
b_rak = [0.0, -90.0, 90.0, 90.0, 90.0, 0.0]
b_wd = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
b_len = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

time_array_obs = np.arange(npts) * dt

elemental_synthetics = []
hist = np.zeros((6, 8))
for k in range(6):
    print(f"Computing synthetic for orthogonal basis {k + 1}...")
    hist[k, :] = [
        k + 1,
        M0_components[k],
        b_str[k],
        b_dip[k],
        b_rak[k],
        b_wd[k],
        b_len[k],
        0.0,
    ]
# %% run convolution
t, sx, sy, sz = moment.conv(
    ap, hist, source_type=4, t0=4, unit=1
)  # unit=1 for displacement

fs = 1.0 / dt  # Frecuencia de muestreo (1 / 0.25 = 4 Hz)
lowcut = 0.06  # Límite inferior (Hz)
highcut = 0.09  # Límite superior (Hz)
order = 4  # Orden del filtro Butterworth (4 es estándar en sismología)

# Sorting observed waveforms
print("Plotting comparison...")
# Load observed data from Calama (for comparison)
real_disp_x = np.loadtxt(DATA_DIR / "real_disp_x")
real_disp_y = np.loadtxt(DATA_DIR / "real_disp_y")
real_disp_z = np.loadtxt(DATA_DIR / "real_disp_z")

rx = real_disp_x.reshape((7, npts))[:n_stations]
ry = real_disp_y.reshape((7, npts))[:n_stations]
rz = real_disp_z.reshape((7, npts))[:n_stations]

observed_waveforms = np.stack(
    [rx, ry, rz], axis=1
)  # X, Y, Z (observed already filtered)
syntethic_waveforms = np.stack([sx, sy, sz], axis=1)  # X, Y, Z (to be filtered)
import scipy.signal as signal

# 2. Diseñar el filtro pasabanda
b, a = signal.butter(order, [lowcut, highcut], btype="bandpass", fs=fs)

# 3. Crear arreglos vacíos para guardar las señales filtradas
synthetic_filtered = np.zeros_like(syntethic_waveforms)
observed_filtered = observed_waveforms
n_stations = stations.shape[0]
# 4. Aplicar el filtro traza por traza
for i in range(n_stations):
    for j in range(3):
        # Filtramos la sintética (el eje de tiempo es el último índice ':')
        synthetic_filtered[i, j, :] = signal.filtfilt(
            b, a, syntethic_waveforms[i, j, :]
        )

#  Plotting comparison of filtered waveforms
# %% Manual Misfit Calculation (No rotation, L2 Misfit)

misfit = 0.0
misfit_station = np.zeros(n_stations)
for i in range(n_stations):
    for j in range(3):
        obs = observed_filtered[i, j, :]
        syn = synthetic_filtered[i, j, :]
        # Normalize both signals
        obs_max = np.max(np.abs(obs))
        syn_max = np.max(np.abs(syn))
        obs_norm = obs if obs_max == 0 else obs / obs_max
        syn_norm = syn if syn_max == 0 else syn / syn_max
        # Calculate L2 misfit
        comp_misfit = np.sum((obs_norm - syn_norm) ** 2)
        misfit += comp_misfit
        misfit_station[i] += comp_misfit
print(f"Misfit (L2): {misfit:.6f}")
# %%

# Plotting Components comparison
fig, axes = plt.subplots(n_stations, 3, figsize=(10, 0.5 * n_stations), sharex=True)
comps = ["X", "Y", "Z"]

for i in range(n_stations):
    for j in range(3):
        obs_data = observed_filtered[i, j, :]
        syn_data = synthetic_filtered[i, j, :]
        axes[i, j].plot(time_array_obs, obs_data, "k", lw=1.0, label="Observed")
        axes[i, j].plot(
            time_array_obs, syn_data, "r", lw=1.0, label="Synthetic (filtered)"
        )
        if i == 0:
            axes[i, j].set_title(comps[j])
        if j == 0:
            axes[i, j].set_ylabel(f"{name_stations[i]}\nAmp")
        if i == n_stations - 1:
            axes[i, j].set_xlabel("Time (s)")
        if i == 0 and j == 0:
            axes[i, j].legend()

plt.suptitle(
    f"Component Comparison (No rotation) - L2 Misfit: {misfit:.6f}", fontsize=16
)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
