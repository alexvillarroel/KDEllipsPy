"""
pb20_disp.py
============
Genera el desplazamiento (SAC/DISP) para PB20, que se anadio a mano y solo
trae velocidad (SAC/VEL/PB20.BL*.vel.sac, de deconvolucion).

Replica el paso vel->disp de integracion.pssa3 (integracion trapezoidal +
Butterworth CAUSAL 0.02 Hz, orden 4) y rellena la geometria del evento en el
header igual que txt_to_sac.py. No toca el resto de estaciones.

Uso:
    python codigos/pb20_disp.py
"""

import os
import glob

import numpy as np
import obspy
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

from event import load_event
from integracion import integraf, _design_butter
from scipy.signal import lfilter, detrend

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTL_FILE = os.path.join(BASE, "event.ctl")
VEL_DIR = os.path.join(BASE, "SAC", "VEL")
DISP_DIR = os.path.join(BASE, "SAC", "DISP")

# Mismos parametros que txt_to_sac.py (banda ancha, causal).
FILT_F1 = 0.02
FILT_F2 = 40.0
FILT_ORDER = 4


def vel_to_disp(vel, Fs):
    """vel -> disp: integraf + Butterworth causal, igual que pssa3 (1 integracion)."""
    f2 = min(FILT_F2, 0.45 * Fs)
    b, a = _design_butter(Fs, FILT_F1, f2, FILT_ORDER)
    nfi = len(b) * 10
    x = np.concatenate([np.zeros(nfi), detrend(vel), np.zeros(nfi)])
    d = integraf(x, Fs)
    d = lfilter(b, a, d)          # causal (1 pasada)
    return d[nfi:nfi + len(vel)]


def set_geometria(sac, ev, stla, stlo):
    dist_m, az, baz = gps2dist_azimuth(ev["lat"], ev["lon"], stla, stlo)
    epic_km = dist_m / 1000.0
    sac.evla, sac.evlo, sac.evdp = ev["lat"], ev["lon"], ev["depth"]
    sac.dist = np.sqrt(epic_km ** 2 + ev["depth"] ** 2)   # hipocentral km
    sac.gcarc = kilometers2degrees(epic_km)
    sac.az, sac.baz = az, baz
    sac.user0 = epic_km
    sac.lcalda = 0
    sac.idep = 6                                            # IDISP


def main():
    ev = load_event(CTL_FILE)
    files = sorted(glob.glob(os.path.join(VEL_DIR, "PB20.*.vel.sac")))
    assert files, "no encontre PB20.*.vel.sac en SAC/VEL/"
    for f in files:
        tr = obspy.read(f)[0]
        Fs = tr.stats.sampling_rate
        tr.data = vel_to_disp(tr.data, Fs).astype(np.float32)
        set_geometria(tr.stats.sac, ev, tr.stats.sac.stla, tr.stats.sac.stlo)
        out = os.path.join(DISP_DIR, f"{tr.stats.station}.{tr.stats.channel}.disp.sac")
        tr.write(out, format="SAC")
        print(f"  {os.path.basename(out)}  npts={tr.stats.npts}  "
              f"max|d|={np.max(np.abs(tr.data)):.3e} m")
    print(f"\nListo: {len(files)} componentes en {DISP_DIR}/")


def _selfcheck():
    # vel constante -> disp debe crecer ~lineal (pendiente = valor), el
    # highpass causal la atenua pero la integral no debe ser ~0 ni NaN.
    Fs = 20.0
    vel = np.ones(4000)
    d = vel_to_disp(vel, Fs)
    assert np.all(np.isfinite(d)) and np.max(np.abs(d)) > 0
    print("selfcheck OK")


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        _selfcheck()
    else:
        main()
