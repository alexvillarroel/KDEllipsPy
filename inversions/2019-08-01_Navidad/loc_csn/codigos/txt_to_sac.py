"""
txt_to_sac.py
=============
Script PRINCIPAL.

Lee los registros de aceleracion (.txt) del evento de Calama y los deja como
SAC CRUDO (sin integrar ni filtrar) con su metadata:

    aceleracion cruda  ->  SAC/ACC/

La integracion (acc->vel->disp) y el filtrado a la banda los hace despues
``kde-prep`` UNA sola vez, desde este SAC crudo, con la banda del input.ctl y
aplicando el MISMO filtro a los sinteticos (evita el doble filtrado asimetrico).
El crudo se mantiene salvo que cambien los .txt.

Genera UN archivo SAC por componente (3 por estacion).

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


# ==============================================================================
# CONFIGURACION
# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CTL_FILE = os.path.join(BASE, "event.ctl")
TXT_ROOT = os.path.join(BASE, "registros_sismicos")
SAC_ROOT = os.path.join(BASE, "SAC")
OUT = {"ACC": os.path.join(SAC_ROOT, "ACC")}

# Codigo SAC 'idep' para el tipo de variable dependiente.
IDEP = {"ACC": 8}  # IACC (aceleracion)


def asignar_geometria(tr, ev):
    """Rellena los headers SAC de geometria (estacion + evento + distancias)."""
    stla = tr.stats.coordinates["latitude"]
    stlo = tr.stats.coordinates["longitude"]

    # Distancia epicentral, azimut y back-azimut.
    dist_m, az, baz = gps2dist_azimuth(ev["lat"], ev["lon"], stla, stlo)
    epic_km = dist_m / 1000.0
    # Distancia hipocentral: hipotenusa entre epicentral y profundidad.
    hypo_km = np.sqrt(epic_km**2 + ev["depth"] ** 2)

    sac = tr.stats.sac  # ya existe porque seteamos format SAC al escribir
    sac.kstnm = tr.stats.station
    sac.kcmpnm = tr.stats.channel
    sac.stla = stla
    sac.stlo = stlo
    sac.evla = ev["lat"]
    sac.evlo = ev["lon"]
    sac.evdp = ev["depth"]  # km
    sac.dist = hypo_km  # km (HIPOCENTRAL, segun pedido)
    sac.gcarc = kilometers2degrees(epic_km)  # grados (epicentral)
    sac.az = az
    sac.baz = baz
    sac.user0 = epic_km  # distancia epicentral [km], por si se necesita
    sac.lcalda = 0  # no recalcular dist/az automaticamente
    return tr


def procesar():
    ev = load_event(CTL_FILE)
    print(
        "Evento:",
        ev["origin_time"],
        f"| lat {ev['lat']} lon {ev['lon']} prof {ev['depth']} km",
    )

    for d in OUT.values():
        os.makedirs(d, exist_ok=True)

    stream = read_all(TXT_ROOT)
    print(
        f"Leidas {len(stream)} trazas "
        f"de {len({tr.stats.station for tr in stream})} estaciones.\n"
    )

    n_ok = 0
    for tr in stream:
        # Guardia: el pipeline asume crudo en ACELERACION (kde-prep integra acc->
        # disp/vel). Si algun .txt llega en velocidad, falla fuerte aqui.
        units = getattr(tr.stats, "units", "")
        assert "seg/seg" in units, (
            f"{tr.stats.station}.{tr.stats.channel}: esperaba aceleracion, "
            f"vino unidades={units!r}. Extender el pipeline si hay velocidad cruda.")

        out_tr = tr.copy()
        out_tr.data = np.asarray(tr.data, dtype=np.float32)   # crudo, sin tocar
        out_tr.stats.sac = getattr(out_tr.stats, "sac", None) or _empty_sac()
        asignar_geometria(out_tr, ev)
        out_tr.stats.sac.idep = IDEP["ACC"]
        fname = f"{tr.stats.station}.{tr.stats.channel}.acc.sac"
        out_tr.write(os.path.join(OUT["ACC"], fname), format="SAC")
        n_ok += 1

    print(f"\nListo. {n_ok} trazas -> {n_ok} archivos SAC crudos en {SAC_ROOT}/ACC/")
    print("  ACC/  aceleracion CRUDA (m/s^2). Integra/filtra kde-prep desde aqui.")


def _empty_sac():
    """Crea un AttribDict SAC vacio para poder asignar headers antes de escribir."""
    from obspy.core import AttribDict

    return AttribDict()


if __name__ == "__main__":
    procesar()
