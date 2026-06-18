"""Remuestrea los archivos planos real_vel_{x,y,z} de 512 -> 128 muestras por
estacion (dt 0.25 s -> 1.0 s), manteniendo la ventana de 128 s.

La data esta filtrada a 0.02-0.1 Hz, muy por debajo del nuevo Nyquist (0.5 Hz),
asi que la decimacion es sin aliasing (lossless). Respaldo: *_512.bak."""
import os
import shutil
import numpy as np
from scipy.signal import resample

HERE = os.path.dirname(os.path.abspath(__file__))
NSTA = 8
NPTS_OLD, NPTS_NEW = 512, 128

for comp in ("x", "y", "z"):
    f = os.path.join(HERE, f"real_vel_{comp}")
    bak = f + "_512.bak"
    if not os.path.exists(bak):
        shutil.copy2(f, bak)
    arr = np.loadtxt(f, dtype=float).ravel()
    assert arr.size == NSTA * NPTS_OLD, f"{f}: {arr.size} != {NSTA*NPTS_OLD}"
    blocks = arr.reshape(NSTA, NPTS_OLD)
    out = np.empty((NSTA, NPTS_NEW), dtype=float)
    for i in range(NSTA):
        out[i] = resample(blocks[i], NPTS_NEW)
    np.savetxt(f, out.ravel(), fmt="%.8e")
    print(f"  real_vel_{comp}: {NSTA}x{NPTS_OLD} -> {NSTA}x{NPTS_NEW} = {out.size} filas")

print("Listo. Backups *_512.bak.")
