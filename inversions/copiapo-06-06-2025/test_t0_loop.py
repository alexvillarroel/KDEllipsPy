import sys
from pathlib import Path
import numpy as np
import scipy.signal as signal
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

# Add project root to path
ROOT_DIR = Path("/home/alex/KDEllipsPy")
sys.path.insert(0, str(ROOT_DIR))

from kdellipspy.axitra.src.axitra import Axitra, moment

DATA_DIR = ROOT_DIR / "inversions/copiapo-06-06-2025/SELECTED_DATA"

# 1. Velocity model (from input.ctl)
model = np.array([
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
], dtype=float)

# 2. Stations
stations = np.array([
    [1, -27.1040, -70.8490, 0.0000],
    [2, -27.5030, -70.8870, 0.0000],
    [3, -23.6100, -70.2620, 0.0000],
    [4, -26.1480, -70.5990, 0.0000],
    [5, -25.1626, -69.5904, 0.0000],
    [6, -27.5940, -70.2350, 0.0000],
    [7, -23.9050, -69.2910, 0.0000],
], dtype=float)

# 3. Sources
sources = np.array([
    [1, -26.639, -70.404, 75000.0],
    [2, -26.639, -70.404, 75000.0],
    [3, -26.639, -70.404, 75000.0],
    [4, -26.639, -70.404, 75000.0],
    [5, -26.639, -70.404, 75000.0],
    [6, -26.639, -70.404, 75000.0],
], dtype=float)

fmax = 2.0
duration = 128.0
dt = 0.25
npts = int(duration / dt)

axpath = str((ROOT_DIR / "kdellipspy" / "axitra" / "src").resolve())

# Initialize Axitra and Green functions
ap = Axitra(model=model, stations=stations, sources=sources, fmax=fmax, duration=duration, latlon=True, axpath=axpath, id=123)
ap = moment.green(ap)

# Moment tensor
mrr, mtt, mpp, mrt, mrp, mtp, iexp = -2.486, 7.932, 0.083, 0.038, -0.739, 0.473, 18.0

def mt_decomposition(mrr, mtt, mpp, mrt, mrp, mtp, iexp):
    def vec(mat): return np.array([mat[0, 0], mat[1, 1], mat[2, 2], mat[0, 1], mat[0, 2], mat[1, 2]])
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
    M0_components = [abs(c) / np.sqrt(2) * np.linalg.norm(b, "fro") for c, b in zip(coeffs, basis)]
    return M0_components

M0_components = mt_decomposition(mrr, mtt, mpp, mrt, mrp, mtp, iexp)
b_str = [0.0, 270.0, 0.0, 90.0, 0.0, 0.0]
b_dip = [90.0, 90.0, 90.0, 45.0, 45.0, 0.0]
b_rak = [0.0, -90.0, 90.0, 90.0, 90.0, 0.0]
b_wd = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
b_len = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

hist = np.zeros((6, 8))
for k in range(6):
    hist[k, :] = [k + 1, M0_components[k], b_str[k], b_dip[k], b_rak[k], b_wd[k], b_len[k], 0.0]

# Load observed data
real_vel_x = np.loadtxt(DATA_DIR / "real_vel_x")
real_vel_y = np.loadtxt(DATA_DIR / "real_vel_y")
real_vel_z = np.loadtxt(DATA_DIR / "real_vel_z")

n_stations = stations.shape[0]
rx = real_vel_x.reshape((7, 512))[:n_stations]
ry = real_vel_y.reshape((7, 512))[:n_stations]
rz = real_vel_z.reshape((7, 512))[:n_stations]
observed_waveforms = np.stack([rx, ry, rz], axis=1)

# Filter parameters
fs = 1.0 / dt
lowcut, highcut, order = 0.1, 0.2, 4
b, a = signal.butter(order, [lowcut, highcut], btype="bandpass", fs=fs)

observed_filtered = np.zeros_like(observed_waveforms)
for i in range(n_stations):
    for j in range(3):
        observed_filtered[i, j, :] = signal.filtfilt(b, a, observed_waveforms[i, j, :])

# Windowing parameters
taup_model = TauPyModel(model="iasp91")
evla, evlo, evdep_km = -26.639, -70.404, 75.0
tp_list, ts_list = [], []
for i in range(n_stations):
    stla, stlo = stations[i, 1], stations[i, 2]
    dist_m, _, _ = gps2dist_azimuth(evla, evlo, stla, stlo)
    dist_deg = dist_m / (111.19 * 1000.0)
    arrivals_p = taup_model.get_travel_times(source_depth_in_km=evdep_km, distance_in_degree=dist_deg, phase_list=["P", "p"])
    arrivals_s = taup_model.get_travel_times(source_depth_in_km=evdep_km, distance_in_degree=dist_deg, phase_list=["S", "s"])
    tp = arrivals_p[0].time if arrivals_p else 10.0
    ts = arrivals_s[0].time if arrivals_s else 25.0
    tp_list.append(tp)
    ts_list.append(ts)

window_s = 20.0
sampling_rate = 1.0 / dt

best_misfit = float('inf')
best_t0 = None

# Testing t0 from 0 to 20 with step 0.25
for t0 in np.arange(0, 20.1, 0.25):
    _, sx, sy, sz = moment.conv(ap, hist, source_type=4, t0=t0, unit=2)
    syn_waveforms = np.stack([sx, sy, sz], axis=1)
    syn_filtered = np.zeros_like(syn_waveforms)
    for i in range(n_stations):
        for j in range(3):
            syn_filtered[i, j, :] = signal.filtfilt(b, a, syn_waveforms[i, j, :])
    
    num, den = 0.0, 0.0
    for i in range(n_stations):
        tp, ts = tp_list[i], ts_list[i]
        kp0, kp1 = int((tp - 1.0) * sampling_rate), int((tp - 1.0 + window_s) * sampling_rate)
        ks0, ks1 = int((ts - 1.0) * sampling_rate), int((ts - 1.0 + window_s) * sampling_rate)
        kp0, kp1 = max(0, kp0), min(npts - 1, kp1)
        ks0, ks1 = max(0, ks0), min(npts - 1, ks1)

        # Unrotated windowed misfit
        num += np.sum((observed_filtered[i, 0, kp0:kp1] - syn_filtered[i, 0, kp0:kp1]) ** 2)
        den += np.sum(observed_filtered[i, 0, kp0:kp1] ** 2)
        num += np.sum((observed_filtered[i, 2, kp0:kp1] - syn_filtered[i, 2, kp0:kp1]) ** 2)
        den += np.sum(observed_filtered[i, 2, kp0:kp1] ** 2)
        num += np.sum((observed_filtered[i, 1, ks0:ks1] - syn_filtered[i, 1, ks0:ks1]) ** 2)
        den += np.sum(observed_filtered[i, 1, ks0:ks1] ** 2)

    misfit = num / den if den > 0 else 0.0
    if misfit < best_misfit:
        best_misfit = misfit
        best_t0 = t0
    # print(f"t0={t0:.2f}, misfit={misfit:.6f}")

print(f"Best t0: {best_t0:.2f} with misfit: {best_misfit:.6f}")
