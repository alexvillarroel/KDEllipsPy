"""
io_sismica.py
=============
Lectura de los registros de aceleracion en formato ``.txt`` del CSN y
construccion de objetos ObsPy (``Trace`` / ``Stream``) con su metadata.

Cabecera de cada archivo::

    # Tiempo de Origen: 2026-05-25T21:51:52.000000Z   <- starttime de la traza
    # Tasa de muestreo: 200.0 muestras/seg
    # Numero total de muestras: 43800
    # Estacion: A02F Componente: HNN
    # Latitud: -23.101 Longitud: -70.445
    # Unidades: m/seg/seg
    <una muestra por linea>

Nota: "Tiempo de Origen" en estos archivos es el tiempo de la PRIMERA muestra
de la traza (no el tiempo de origen del sismo), por eso lo usamos como
``starttime``.
"""

import glob
import os
import re

import numpy as np
from obspy import Stream, Trace, UTCDateTime


# Coordenadas (lat/lon) que vienen en el header las guardamos en el Trace
# usando un sub-diccionario propio para no perderlas antes de escribir el SAC.
def _parse_header(path):
    """Lee solo las lineas de comentario y devuelve un dict con la metadata."""
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            txt = line.lstrip("#").strip()
            if txt.startswith("Tiempo de Origen:"):
                meta["starttime"] = UTCDateTime(txt.split(":", 1)[1].strip())
            elif txt.startswith("Tasa de muestreo:"):
                # Tomamos el numero que va despues de los dos puntos.
                meta["sampling_rate"] = float(
                    re.search(r":\s*([-+]?\d*\.?\d+)", txt).group(1))
            elif txt.startswith("Numero total de muestras:"):
                meta["npts"] = int(re.findall(r"\d+", txt)[0])
            elif txt.startswith("Estacion:"):
                # "Estacion: A02F Componente: HNN"
                m = re.search(r"Estacion:\s*(\S+)\s+Componente:\s*(\S+)", txt)
                if m:
                    meta["station"] = m.group(1)
                    meta["channel"] = m.group(2)
            elif txt.startswith("Latitud:"):
                nums = re.findall(r"[-+]?\d+\.?\d*", txt)
                meta["latitude"] = float(nums[0])
                meta["longitude"] = float(nums[1])
            elif txt.startswith("Unidades:"):
                meta["units"] = txt.split(":", 1)[1].strip()
    return meta


def read_txt_trace(path):
    """Lee un archivo ``.txt`` y devuelve un ``obspy.Trace``.

    Las coordenadas quedan en ``tr.stats.coordinates`` (lat/lon) para que
    ``txt_to_sac`` pueda volcarlas a los headers SAC.
    """
    meta = _parse_header(path)
    data = np.loadtxt(path, comments="#", dtype=np.float64)

    tr = Trace(data=data)
    tr.stats.station = meta.get("station", "")
    tr.stats.channel = meta.get("channel", "")
    tr.stats.sampling_rate = meta.get("sampling_rate", 1.0)
    tr.stats.starttime = meta.get("starttime", UTCDateTime(0))
    # Metadata extra (no estandar de obspy, pero util):
    tr.stats.coordinates = {
        "latitude": meta.get("latitude"),
        "longitude": meta.get("longitude"),
    }
    tr.stats.units = meta.get("units", "m/seg/seg")
    return tr


def read_all(root):
    """Lee todos los ``.txt`` bajo ``root/<ESTACION>/*.txt`` -> ``Stream``."""
    stream = Stream()
    pattern = os.path.join(root, "*", "*.txt")
    for path in sorted(glob.glob(pattern)):
        try:
            stream += read_txt_trace(path)
        except Exception as exc:  # noqa: BLE001
            print(f"[io_sismica] No se pudo leer {path}: {exc}")
    return stream


if __name__ == "__main__":
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    st = read_all(os.path.join(here, "registros_sismicos"))
    print(st)
    print(f"\nTotal trazas: {len(st)}")
    estaciones = sorted({tr.stats.station for tr in st})
    print(f"Estaciones ({len(estaciones)}): {', '.join(estaciones)}")
