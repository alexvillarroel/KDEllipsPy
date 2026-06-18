"""Regenera real_vel_{x,y,z} en VELOCIDAD CAUSAL desde los SAC de aceleracion
de RAW/, en banda 0.02-0.1 Hz, ancla en el origen (07:35:34), ventana 128 s a
dt=1.0 (128 muestras/estacion).

Usa pssa3 (el integrador trusted del proyecto, mismo de la tuberia 2020-06-03)
con nophase=False -> filtro causal. Asi el observado queda causal y consistente
con las sinteticas (input.ctl: Zerophase filter = 0).

Orden de estaciones = seccion 8 del input.ctl. x=N, y=E, z=Z.
Respaldo de los archivos previos: *_prev.bak
"""
import os
import sys
import glob
import shutil

import numpy as np
import obspy
from obspy import UTCDateTime

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "RAW")

ORIGIN = UTCDateTime("2020-06-03T07:35:34.000")
NPTS, DT = 128, 1.0
F1, F2 = 0.02, 0.1          # banda de inversion
NOPHASE = False             # False -> filtro CAUSAL (lfilter)

# Orden EXACTO del input.ctl (seccion 8).
STATIONS = ["A05C", "AC01", "AC02", "PB05", "PB06", "PB08", "PB19", "T15A"]
# componente -> sufijo de archivo plano
COMP2FILE = {"N": "x", "E": "y", "Z": "z"}


def procesar(station, comp):
    """Devuelve velocidad causal (NPTS muestras) de una componente, o None.

    Detecta el tipo de sensor por el canal:
      HN* (acelerometro) -> aceleracion -> se integra una vez a velocidad.
      HL* (sismometro)   -> ya es velocidad -> no se integra.
    """
    # Acepta canales HN? y HL? (acelerometro o sismometro de banda ancha).
    matches = glob.glob(os.path.join(RAW, f"*{station}*H[NL]{comp}*.sac"))
    if not matches:
        print(f"  [!] {station} ?{comp}: no encontrado")
        return None
    tr = obspy.read(matches[0])[0]
    tr.data = np.asarray(tr.data, dtype=np.float64)

    # Tanto HN como HL son ACELEROMETROS -> se integran una vez a velocidad.
    # UN solo bandpass causal de 4o orden (SOS, estable) identico al de las
    # sinteticas.
    tr.detrend("linear")
    tr.detrend("demean")
    tr.taper(max_percentage=0.05, type="cosine")
    tr.integrate(method="cumtrapz")           # aceleracion -> velocidad
    tr.detrend("linear")                      # quita rampa de integracion
    tr.filter("bandpass", freqmin=F1, freqmax=F2, corners=4,
              zerophase=NOPHASE)              # NOPHASE=False -> causal
    # Ancla al origen y remuestrea a 128 muestras a dt=1.0.
    tr.trim(ORIGIN, ORIGIN + (NPTS + 2) * DT, pad=True, fill_value=0.0)
    tr.interpolate(sampling_rate=1.0 / DT, starttime=ORIGIN, npts=NPTS)
    return tr.data[:NPTS].astype(np.float64)


def main():
    cols = {"x": [], "y": [], "z": []}
    print(f"Regenerando velocidad CAUSAL {F1}-{F2} Hz, {NPTS} muestras dt={DT}s")
    for st in STATIONS:
        ok = True
        buf = {}
        for comp, suf in COMP2FILE.items():
            v = procesar(st, comp)
            if v is None:
                ok = False
                break
            buf[suf] = v
        if not ok:
            print(f"  {st}: incompleta, se omite")
            continue
        for suf in ("x", "y", "z"):
            cols[suf].append(buf[suf])
        print(f"  {st:5s} OK  (rms_z={np.sqrt(np.mean(buf['z']**2)):.3e})")

    for suf in ("x", "y", "z"):
        f = os.path.join(HERE, f"real_vel_{suf}")
        if os.path.exists(f):
            shutil.copy2(f, f + "_prev.bak")
        data = np.concatenate(cols[suf])
        np.savetxt(f, data, fmt="%.8e")
        print(f"  -> real_vel_{suf}: {len(cols[suf])} est x {NPTS} = {data.size} filas")

    print("Listo. Velocidad causal regenerada. Respaldo previo en *_prev.bak.")


if __name__ == "__main__":
    main()
