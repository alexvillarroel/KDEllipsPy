import sys
from pathlib import Path
import numpy as np

ROOT_DIR = Path("/home/alex/KDEllipsPy")
sys.path.insert(0, str(ROOT_DIR))
from kdellipspy.axitra.src.axitra import Axitra, moment

# --- Setup simple test ---
model = np.array([[0.0, 4700.0, 2750.0, 2562.0, 300.0, 300.0], 
                  [5000.0, 5730.0, 3300.0, 2692.0, 300.0, 300.0]], dtype=float)
stations = np.array([[1, -27.1, -70.8, 0.0]], dtype=float)
sources = np.array([[1, -26.6, -70.4, 75000.0]]*6, dtype=float)
fmax, duration, dt = 1.0, 100.0, 0.5
axpath = str((ROOT_DIR / "kdellipspy" / "axitra" / "src").resolve())

ap = Axitra(model=model, stations=stations, sources=sources, fmax=fmax, duration=duration, latlon=True, axpath=axpath, id=999)
ap = moment.green(ap)

# Moment tensor orthogonal basis coefficients (arbitrary for this test)
M0_c = [1e17] * 6
b_s, b_d, b_r, b_w = [0, 270, 0, 90, 0, 0], [90, 90, 90, 45, 45, 0], [0, -90, 90, 90, 90, 0], [0, 0, 0, 0, 0, -1]
hist = np.array([[i+1, M0_c[i], b_s[i], b_d[i], b_r[i], b_w[i], 0, 0] for i in range(6)])

# 1. Compute Displacement (unit=1)
t1, sx1, sy1, sz1 = moment.conv(ap, hist, source_type=4, t0=1.5, unit=1)

# 2. Compute Velocity (unit=2)
t2, sx2, sy2, sz2 = moment.conv(ap, hist, source_type=4, t0=1.5, unit=2)

print("\n--- Comparison of Units ---")
print(f"Max Amplitude Displacement (Z): {np.max(np.abs(sz1)):.2e}")
print(f"Max Amplitude Velocity (Z):     {np.max(np.abs(sz2)):.2e}")

identical = np.allclose(sz1, sz2)
print(f"\nAre the waveforms identical? {'YES' if identical else 'NO (as expected)'}")

# check ratio (velocity should be approx displacement * omega)
# in time domain, for a dominant frequency f, vel ~ 2*pi*f * disp
if not identical:
    ratio = np.max(np.abs(sz2)) / np.max(np.abs(sz1))
    print(f"Ratio Peak Velocity / Peak Displacement: {ratio:.4f}")
