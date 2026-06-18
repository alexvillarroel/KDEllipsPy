"""
descargar_egf.py
================
Descarga formas de onda + metadata (respuesta) para los mejores candidatos
EGF del ranking (candidatos_ranking.csv).

Fuentes:
  - GFZ        : red CX (IPOC; PB01..PB16, estaciones que compartimos con el
                 mainshock: PB05, PB06, PB08, PB19).
  - EARTHSCOPE : red C1 (CSN banda ancha) dentro de 2.5 grados del epicentro.

Por cada candidato crea  egf/<YYYYMMDDTHHMMSS>/RAW/  con:
    <NET>.<STA>.mseed   (todas las componentes)
    inventory.xml       (respuestas de todo lo descargado)

Uso:
    python codigos/descargar_egf.py [n_top]     (default 3)
"""
import os
import sys
import csv

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../egf
RANKING = os.path.join(BASE, "candidatos_ranking.csv")

T_PRE, T_POST = 120.0, 400.0          # s antes/despues del origen
RADIUS_DEG = 2.5                      # para estaciones C1

# Solicitudes: (cliente, red, estaciones, canales)
REQUESTS = [
    ("GFZ", "CX", "PB*", "HH?,HL?,HN?"),
    ("EARTHSCOPE", "C1", "*", "HH?,BH?"),
]


def descargar_evento(row, clients):
    t0 = UTCDateTime(row["time"])
    tag = t0.strftime("%Y%m%dT%H%M%S")
    outdir = os.path.join(BASE, tag, "RAW")
    os.makedirs(outdir, exist_ok=True)
    print(f"\n=== {tag}  M{row['mag']}  d3D={row['dist3d_km']} km ===")

    inv_total = None
    n_tr = 0
    for prov, net, sta, cha in REQUESTS:
        c = clients.get(prov)
        if c is None:
            continue
        try:
            inv = c.get_stations(
                network=net, station=sta, channel=cha, level="response",
                latitude=float(row["lat"]), longitude=float(row["lon"]),
                maxradius=RADIUS_DEG,
                starttime=t0 - T_PRE, endtime=t0 + T_POST,
            )
        except Exception as ex:  # noqa: BLE001
            print(f"  [{prov}/{net}] sin estaciones: {str(ex)[:70]}")
            continue
        inv_total = inv if inv_total is None else inv_total + inv

        for network in inv:
            for station in network:
                try:
                    st = c.get_waveforms(
                        network.code, station.code, "*", cha,
                        t0 - T_PRE, t0 + T_POST,
                    )
                except Exception:
                    continue
                if len(st) == 0:
                    continue
                st = st.merge(method=1, fill_value=0)
                fn = os.path.join(outdir, f"{network.code}.{station.code}.mseed")
                st.write(fn, format="MSEED")
                n_tr += len(st)
                print(f"  {network.code}.{station.code}: {len(st)} trazas")

    if inv_total is not None:
        inv_total.write(os.path.join(outdir, "inventory.xml"),
                        format="STATIONXML")
    print(f"  -> {n_tr} trazas en {outdir}")
    return n_tr


def main():
    n_top = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    with open(RANKING) as f:
        rows = list(csv.DictReader(f))[:n_top]

    clients = {}
    for prov in {p for p, *_ in REQUESTS}:
        try:
            clients[prov] = Client(prov, timeout=120)
        except Exception as ex:  # noqa: BLE001
            print(f"[{prov}] no disponible: {ex}")

    total = 0
    for row in rows:
        try:
            total += descargar_evento(row, clients)
        except Exception as ex:  # noqa: BLE001
            print(f"  ERROR evento {row['time']}: {ex}")

    print(f"\nListo: {total} trazas en total para {len(rows)} candidatos.")


if __name__ == "__main__":
    main()
