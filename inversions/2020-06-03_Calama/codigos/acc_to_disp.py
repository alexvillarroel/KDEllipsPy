"""
acc_to_disp.py
==============
Doble-integra la aceleracion cruda de SAC/ACC a velocidad y desplazamiento y
escribe los SAC correspondientes:

    SAC/ACC/*.acc.sac  ->  SAC/VEL/*.vel.sac   (velocidad,     m/s)
                       ->  SAC/DISP/*.disp.sac (desplazamiento, m)

Usa ``pssa3`` (codigos/integracion.py, port del MATLAB del grupo): Butterworth
pasa-banda + integracion trapezoidal, fase cero. Banda ancha 0.02-40 Hz (bajo la
banda de inversion 0.04-0.15) para dejar el filtrado fino a la inversion.

Se conservan todos los headers del SAC crudo; solo cambia ``idep``.

Uso:
    python codigos/acc_to_disp.py
"""

import os
import glob

import numpy as np
import obspy

from integracion import pssa3

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ACC_DIR = os.path.join(BASE, "SAC", "ACC")
VEL_DIR = os.path.join(BASE, "SAC", "VEL")
DISP_DIR = os.path.join(BASE, "SAC", "DISP")

F1, F2 = 0.04, 40.0  # esquinas del pasa-banda para integrar
IDEP = {"vel": 7, "disp": 6}  # IVEL / IDISP (ACC=8 en txt_to_sac.py)


def escribir(tr, data, out_dir, kind):
    out = tr.copy()
    out.data = np.asarray(data, dtype=np.float32)
    out.stats.sac.idep = IDEP[kind]
    fname = f"{tr.stats.station}.{tr.stats.channel}.{kind}.sac"
    out.write(os.path.join(out_dir, fname), format="SAC")


def main():
    os.makedirs(VEL_DIR, exist_ok=True)
    os.makedirs(DISP_DIR, exist_ok=True)

    files = sorted(glob.glob(os.path.join(ACC_DIR, "*.acc.sac")))
    for f in files:
        tr = obspy.read(f)[0]
        _ac, v, d = pssa3(
            tr.data.astype(np.float64), tr.stats.sampling_rate, f1=F1, f2=F2
        )
        escribir(tr, v, VEL_DIR, "vel")
        escribir(tr, d, DISP_DIR, "disp")

    print(
        f"Listo. {len(files)} trazas -> SAC/VEL, {len(files)} -> SAC/DISP "
        f"(banda {F1}-{F2} Hz)."
    )


if __name__ == "__main__":
    main()
