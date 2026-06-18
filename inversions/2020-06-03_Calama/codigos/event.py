"""
event.py
========
Lectura del archivo de control del evento (``event.ctl``).

Formato esperado (una linea de comentario + una linea de datos)::

    # Lat Lon Depth Strike Dip Rake UTC-DateTime
    -22.38  -68.76  114  201  20  -70  2026-05-25-T21:52:18.000

El orden de los angulos es el estandar sismologico **Strike Dip Rake**.
Los campos Strike/Dip/Rake son opcionales: si falta alguno se devuelve ``None``
en su lugar y se avisa por pantalla, en vez de fallar.
"""

from obspy import UTCDateTime


def _parse_datetime(token):
    """Convierte '2026-05-25-T21:52:18.000' -> UTCDateTime."""
    # El formato del CSN usa un guion extra antes de la 'T'; lo normalizamos.
    token = token.strip().replace("-T", "T")
    return UTCDateTime(token)


def load_event(ctl_path):
    """Lee ``event.ctl`` y devuelve un diccionario con los parametros del evento.

    Returns
    -------
    dict con llaves:
        lat, lon, depth      : float  (grados, grados, km)
        strike, rake, dip    : float | None
        origin_time          : obspy.UTCDateTime
    """
    with open(ctl_path) as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        raise ValueError(f"No se encontro linea de datos en {ctl_path}")

    fields = lines[0].split()

    # El ultimo token siempre es la fecha-hora UTC.
    origin_time = _parse_datetime(fields[-1])
    numeric = [float(x) for x in fields[:-1]]

    # Orden: Lat Lon Depth [Strike Dip Rake]  (convencion sismologica estandar)
    lat, lon, depth = numeric[0], numeric[1], numeric[2]
    strike = numeric[3] if len(numeric) > 3 else None
    dip = numeric[4] if len(numeric) > 4 else None
    rake = numeric[5] if len(numeric) > 5 else None

    if None in (strike, dip, rake):
        print("[event.py] AVISO: falta algun angulo en event.ctl "
              "(strike/dip/rake = %s/%s/%s). Complétalo para el mecanismo focal."
              % (strike, dip, rake))

    return {
        "lat": lat,
        "lon": lon,
        "depth": depth,
        "strike": strike,
        "rake": rake,
        "dip": dip,
        "origin_time": origin_time,
    }


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ev = load_event(os.path.join(here, "event.ctl"))
    for k, v in ev.items():
        print(f"{k:12s}: {v}")
