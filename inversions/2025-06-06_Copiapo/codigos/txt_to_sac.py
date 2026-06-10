"""
txt_to_sac.py
=============
Script PRINCIPAL.

Lee los registros de aceleracion (.txt) del evento de Calama, los procesa con
``integracion.pssa3`` (port del MATLAB) para obtener:

    aceleracion corregida  ->  SAC/ACC/
    velocidad              ->  SAC/VEL/
    desplazamiento         ->  SAC/DISP/

Genera UN archivo SAC por componente y por tipo de senal, es decir 9 SAC por
estacion (3 componentes x {acc, vel, desp}).

Cada SAC incluye en su header:
    kstnm  (estacion), kcmpnm (canal),
    stla / stlo        (lat/lon de la estacion),
    evla / evlo / evdp (lat/lon/prof del hipocentro, desde event.ctl),
    dist               (distancia HIPOCENTRAL en km),
    gcarc / az / baz   (distancia epicentral en grados, azimut y back-azimut),
    idep               (tipo de dato: aceleracion / velocidad / desplazamiento).

Uso:
    python codigos/txt_to_sac.py
"""

import os

import numpy as np
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

from event import load_event
from io_sismica import read_all
from integracion import pssa3


# ==============================================================================
# CONFIGURACION
# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CTL_FILE = os.path.join(BASE, "event.ctl")
TXT_ROOT = os.path.join(BASE, "registros_sismicos")
SAC_ROOT = os.path.join(BASE, "SAC")
OUT = {"ACC": os.path.join(SAC_ROOT, "ACC"),
       "VEL": os.path.join(SAC_ROOT, "VEL"),
       "DISP": os.path.join(SAC_ROOT, "DISP")}

# Banda de filtrado para la integracion (Butterworth).
# f1 = esquina baja (quita la deriva al integrar); f2 = esquina alta (anti-alias).
# Para la inversion de momento tensor luego se re-filtra a la banda fina
# (p.ej. 0.06-0.15 Hz); aqui dejamos productos de banda ancha.
FILT_F1 = 0.02      # Hz
FILT_F2 = 40.0      # Hz  (se recorta a 0.999*Nyquist por traza si hace falta)
FILT_ORDER = 4
ZEROPHASE = False   # False = causal (1 pasada, lfilter); True = acausal (filtfilt)

# Codigo SAC 'idep' para el tipo de variable dependiente.
IDEP = {"ACC": 8, "VEL": 7, "DISP": 6}   # IACC, IVEL, IDISP


def asignar_geometria(tr, ev):
    """Rellena los headers SAC de geometria (estacion + evento + distancias)."""
    stla = tr.stats.coordinates["latitude"]
    stlo = tr.stats.coordinates["longitude"]

    # Distancia epicentral, azimut y back-azimut.
    dist_m, az, baz = gps2dist_azimuth(ev["lat"], ev["lon"], stla, stlo)
    epic_km = dist_m / 1000.0
    # Distancia hipocentral: hipotenusa entre epicentral y profundidad.
    hypo_km = np.sqrt(epic_km ** 2 + ev["depth"] ** 2)

    sac = tr.stats.sac  # ya existe porque seteamos format SAC al escribir
    sac.kstnm = tr.stats.station
    sac.kcmpnm = tr.stats.channel
    sac.stla = stla
    sac.stlo = stlo
    sac.evla = ev["lat"]
    sac.evlo = ev["lon"]
    sac.evdp = ev["depth"]           # km
    sac.dist = hypo_km               # km (HIPOCENTRAL, segun pedido)
    sac.gcarc = kilometers2degrees(epic_km)   # grados (epicentral)
    sac.az = az
    sac.baz = baz
    sac.user0 = epic_km              # distancia epicentral [km], por si se necesita
    sac.lcalda = 0                   # no recalcular dist/az automaticamente
    return tr


def procesar():
    ev = load_event(CTL_FILE)
    print("Evento:", ev["origin_time"],
          f"| lat {ev['lat']} lon {ev['lon']} prof {ev['depth']} km")

    for d in OUT.values():
        os.makedirs(d, exist_ok=True)

    stream = read_all(TXT_ROOT)
    print(f"Leidas {len(stream)} trazas "
          f"de {len({tr.stats.station for tr in stream})} estaciones.\n")

    n_ok = 0
    for tr in stream:
        Fs = tr.stats.sampling_rate
        f2 = min(FILT_F2, 0.45 * Fs)   # margen seguro bajo Nyquist
        ac, v, d = pssa3(tr.data, Fs, f1=FILT_F1, f2=f2,
                         N=FILT_ORDER, nophase=ZEROPHASE, izero=True)

        for kind, signal in (("ACC", ac), ("VEL", v), ("DISP", d)):
            out_tr = tr.copy()
            out_tr.data = np.asarray(signal, dtype=np.float32)
            # Forzamos la creacion del header SAC escribiendo a un buffer:
            out_tr.stats.sac = getattr(out_tr.stats, "sac", None) or _empty_sac()
            asignar_geometria(out_tr, ev)
            out_tr.stats.sac.idep = IDEP[kind]
            fname = f"{tr.stats.station}.{tr.stats.channel}.{kind.lower()}.sac"
            out_tr.write(os.path.join(OUT[kind], fname), format="SAC")
        n_ok += 1

    print(f"\nListo. {n_ok} trazas -> {n_ok * 3} archivos SAC en {SAC_ROOT}/")
    print("  ACC/  velocidad? no: aceleracion corregida (m/s^2)")
    print("  VEL/  velocidad (m/s)")
    print("  DISP/ desplazamiento (m)")


def _empty_sac():
    """Crea un AttribDict SAC vacio para poder asignar headers antes de escribir."""
    from obspy.core import AttribDict
    return AttribDict()


if __name__ == "__main__":
    procesar()
