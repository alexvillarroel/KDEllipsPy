import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Ensure kdellipspy is in path if run from here
sys.path.insert(0, "/home/alex/KDEllipsPy")
import kdellipspy as kde

# 1. Load configuration
cfg = kde.ConfigParser(filepath=Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/input.ctl"))

# Hardcode the model (Small ellipse at hypocenter)
# Parameters: [a1, a2, theta, np, tp, dmax, vr]
# a1, a2 = 2.0 km (small enough to be almost a point, but >0 so it covers at least 1 subfault)
test_model = np.array([2.0, 2.0, 0.0, 0.0, 0.0, 1.0, 2.5])

# Ensure MT strict is set
cfg.moment_tensor.scaling_mode = "mt_strict"

# 2. Build geometry and apply ellipse
fm = kde.AxitraForwardModel.from_config(cfg)
geom = fm.build_geometry_with_ellipse_slip(test_model)

print(f"Number of subfaults: {geom.nsubfaults}")
print(f"Number of source points (Axitra): {geom.nsources}")
print(f"Total Moment (Nm): {geom.total_moment_nm():.2e}")

# Check slip values
slips = [sp.displacement for sp in geom.source_points if sp.displacement > 0]
if slips:
    print(f"Max slip applied: {max(slips):.4f} m")
else:
    print("WARNING: No slip applied. Check ellipse parameters.")

# 3. Load observed data from DATA (since DATA_PROCESSED didn't exist)
data_dir = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/DATA")

observed_waveforms, time_array = kde.load_and_filter_observed_data(
    config=cfg,
    data_dir=data_dir,
    prefer_raw=False
)

# 4. Compute Synthetics
ap = fm.build_axitra(geom, quiet=False)
fm.green(ap, quiet=False)
time_syn, sx, sy, sz = fm.conv(ap, geom, quiet=False)

# Format synthetics
n_cut = min(len(time_array), sx.shape[1])
synthetics = np.stack([sz[:, :n_cut], sx[:, :n_cut], sy[:, :n_cut]], axis=1) # Z, N, E

# Filter synthetics to match observed
synthetics = kde.bandpass_filter_waveforms(
    synthetics,
    time_array[:n_cut],
    freq1=float(cfg.ellipse.freq1),
    freq2=float(cfg.ellipse.freq2)
)

# 5. Plotting
comp_labels = ['Z', 'N', 'E']
n_stations = len(cfg.stations.stations)
fig, axes = plt.subplots(n_stations, 3, figsize=(12, 2*n_stations), sharex=True)
if n_stations == 1: axes = axes[np.newaxis, :]

for i, st in enumerate(cfg.stations.stations):
    row_max = max(np.abs(observed_waveforms[i]).max(), np.abs(synthetics[i]).max()) or 1.0
    for j in range(3):
        ax = axes[i, j]
        ax.plot(time_array[:n_cut], observed_waveforms[i, j, :n_cut] / row_max, 'k', label='Obs')
        ax.plot(time_array[:n_cut], synthetics[i, j] / row_max, 'r--', label='Syn')
        if i == 0: ax.set_title(comp_labels[j])
        if j == 0: ax.set_ylabel(st.name)
        ax.set_yticks([])
        
axes[0,0].legend(loc='upper right', fontsize=8)
fig.supxlabel("Time (s)")
plt.tight_layout()
plt.savefig("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/single_forward_test.png", dpi=150)
print("Plot saved to /home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/single_forward_test.png")

