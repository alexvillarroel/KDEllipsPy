import sys
from pathlib import Path
import numpy as np
import scipy.signal as signal
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

ROOT_DIR = Path("/home/alex/KDEllipsPy")
sys.path.insert(0, str(ROOT_DIR))
from kdellipspy.axitra.src.axitra import Axitra, moment

DATA_DIR = ROOT_DIR / "inversions/copiapo-06-06-2025/SELECTED_DATA"

# --- Setup ---
model = np.array([
    [0.0, 4700.0, 2750.0, 2562.0, 300.0, 300.0], [5000.0, 5730.0, 3300.0, 2692.0, 300.0, 300.0],
    [10000.0, 6320.0, 3620.0, 2817.0, 300.0, 300.0], [15000.0, 6680.0, 3830.0, 2946.0, 300.0, 300.0],
    [20000.0, 6900.0, 3950.0, 3026.0, 300.0, 300.0], [25000.0, 7020.0, 4020.0, 3069.0, 300.0, 300.0],
    [30000.0, 7100.0, 4060.0, 3098.0, 300.0, 300.0], [35000.0, 7250.0, 4140.0, 3152.0, 300.0, 300.0],
    [40000.0, 7480.0, 4260.0, 3235.0, 300.0, 300.0], [45000.0, 7790.0, 4430.0, 3346.0, 300.0, 300.0],
    [50000.0, 8070.0, 4570.0, 3447.0, 1000.0, 1000.0], [55000.0, 8230.0, 4640.0, 3505.0, 1000.0, 1000.0],
    [60000.0, 8360.0, 4680.0, 3551.0, 1000.0, 1000.0], [65000.0, 8480.0, 4730.0, 3595.0, 1000.0, 1000.0],
    [70000.0, 8550.0, 4770.0, 3620.0, 1000.0, 1000.0], [75000.0, 8570.0, 4780.0, 3627.0, 1000.0, 1000.0],
    [80000.0, 8550.0, 4760.0, 3620.0, 1000.0, 1000.0], [85000.0, 8500.0, 4740.0, 3602.0, 1000.0, 1000.0],
    [90000.0, 8440.0, 4720.0, 3580.0, 1000.0, 1000.0], [95000.0, 8380.0, 4690.0, 3559.0, 1000.0, 10000.0],
], dtype=float)

stations = np.array([
    [1, -27.1040, -70.8490, 0.0000], [2, -27.5030, -70.8870, 0.0000], [3, -23.6100, -70.2620, 0.0000],
    [4, -26.1480, -70.5990, 0.0000], [5, -25.1626, -69.5904, 0.0000], [6, -27.5940, -70.2350, 0.0000],
    [7, -23.9050, -69.2910, 0.0000],
], dtype=float)

sources = np.array([[1, -26.639, -70.404, 75000.0]]*6, dtype=float)
fmax, duration, dt = 2.0, 128.0, 0.25
axpath = str((ROOT_DIR / "kdellipspy" / "axitra" / "src").resolve())
ap = Axitra(model=model, stations=stations, sources=sources, fmax=fmax, duration=duration, latlon=True, axpath=axpath, id=127)
ap = moment.green(ap)

mrr, mtt, mpp, mrt, mrp, mtp, iexp = -2.486, 7.932, 0.083, 0.038, -0.739, 0.473, 18.0
def get_coeffs(mrr, mtt, mpp, mrt, mrp, mtp, iexp):
    M = np.array([[mrr, mrt, mrp], [mrt, mtt, mtp], [mrp, mtp, mpp]]) * 10**iexp
    m = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
         np.array([[0, 0, 0], [0, 0, -1], [0, -1, 0]]), np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]),
         np.array([[0, 0, 0], [0, -1, 0], [0, 0, 1]]), np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])]
    def vec(mat): return np.array([mat[0,0], mat[1,1], mat[2,2], mat[0,1], mat[0,2], mat[1,2]])
    B = np.array([vec(x) for x in m]).T
    c = np.linalg.solve(B, vec(M))
    return [abs(val)/np.sqrt(2)*np.linalg.norm(m[i], 'fro') for i, val in enumerate(c)]

M0_c = get_coeffs(mrr, mtt, mpp, mrt, mrp, mtp, iexp)
b_s, b_d, b_r, b_w = [0, 270, 0, 90, 0, 0], [90, 90, 90, 45, 45, 0], [0, -90, 90, 90, 90, 0], [0, 0, 0, 0, 0, -1]

rx = np.loadtxt(DATA_DIR / "real_vel_x").reshape((7, 512))[:7]
ry = np.loadtxt(DATA_DIR / "real_vel_y").reshape((7, 512))[:7]
rz = np.loadtxt(DATA_DIR / "real_vel_z").reshape((7, 512))[:7]
obs = np.stack([rx, ry, rz], axis=1)

fs, low, high, order = 1.0/dt, 0.1, 0.2, 4
b_filt, a_filt = signal.butter(order, [low, high], btype="bandpass", fs=fs)
obs_f = signal.filtfilt(b_filt, a_filt, obs, axis=-1)

taup = TauPyModel(model="iasp91")
tp_l, ts_l = [], []
for i in range(7):
    dist_m, _, _ = gps2dist_azimuth(-26.639, -70.404, stations[i,1], stations[i,2])
    tp_l.append(taup.get_travel_times(75.0, dist_m / 111190.0, phase_list=["P", "p"])[0].time)
    ts_l.append(taup.get_travel_times(75.0, dist_m / 111190.0, phase_list=["S", "s"])[0].time)

window_s = 20.0
results = []
fixed_t0 = 1.5

print(f"{'Delay (s)':<10} | {'CC':<8} | {'Opt VR':<8}")
print("-" * 35)

for delay in np.arange(-10, 10.1, 0.5):
    hist = np.array([[i+1, M0_c[i], b_s[i], b_d[i], b_r[i], b_w[i], 0, delay] for i in range(6)])
    _, sx, sy, sz = moment.conv(ap, hist, source_type=4, t0=fixed_t0, unit=2)
    syn_f = signal.filtfilt(b_filt, a_filt, np.stack([sx, sy, sz], axis=1), axis=-1)
    
    all_o, all_s = [], []
    for i in range(7):
        kp0, kp1 = int((tp_l[i]-1)*4), int((tp_l[i]-1+window_s)*4)
        ks0, ks1 = int((ts_l[i]-1)*4), int((ts_l[i]-1+window_s)*4)
        for c, st, en in [(0, kp0, kp1), (2, kp0, kp1), (1, ks0, ks1)]:
            all_o.append(obs_f[i, c, max(0,st):min(512,en)])
            all_s.append(syn_f[i, c, max(0,st):min(512,en)])
    
    fo, fs_arr = np.concatenate(all_o), np.concatenate(all_s)
    cc = np.corrcoef(fo, fs_arr)[0, 1]
    scale = np.sum(fo * fs_arr) / np.sum(fs_arr**2)
    vr = (1 - np.sum((fo - scale*fs_arr)**2) / np.sum(fo**2)) * 100
    results.append((delay, cc, vr))

best = max(results, key=lambda x: x[1])
print(f"\nBest Shift: Delay={best[0]:.2f}s, CC={best[1]:.4f}, VR={best[2]:.2f}%")
