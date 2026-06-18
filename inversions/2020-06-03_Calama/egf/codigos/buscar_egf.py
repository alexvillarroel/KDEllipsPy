"""
buscar_egf.py
=============
Busca y ranquea candidatos a EGF (Empirical Green's Function) para el
mainshock de Calama 2020-06-03 (Mw 6.8, z=123 km).

Criterios de un buen EGF:
  - co-localizado con el mainshock (misma trayectoria fuente->estacion),
  - 1.5 <= dM <= 2.5 (suficientemente chico para ser ~impulsivo en la banda
    0.1-0.3 Hz, suficientemente grande para tener SNR),
  - idealmente mismo mecanismo (para M~4.5 no suele haber MT publicado;
    la similitud se verifica a posteriori con la correlacion cruzada
    EGF<->mainshock en validar_egf.py).

Consulta USGS y GFZ via FDSN, puntua y guarda el ranking.

Salida: egf/candidatos_ranking.csv
Uso:    python codigos/buscar_egf.py
"""
import os
import csv
import math

import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../egf

# Mainshock (tiempo corregido)
MAIN = {
    "time": UTCDateTime("2020-06-03T07:35:34.000"),
    "lat": -23.247, "lon": -68.530, "depth_km": 123.4, "mw": 6.8,
}

# Ventana de busqueda: desde 2 h despues del mainshock (fuera de la coda
# inmediata) hasta 3 anios despues.
T_INI = MAIN["time"] + 2 * 3600
T_FIN = MAIN["time"] + 3 * 365 * 86400

MAXRADIUS_DEG = 0.7          # ~78 km horizontales
DEPTH_MIN, DEPTH_MAX = 80.0, 160.0
MAG_MIN, MAG_MAX = 3.8, 5.5  # dM objetivo ~1.3-3.0; ideal 1.5-2.5
DM_IDEAL = 2.0


def dist3d_km(lat, lon, dep):
    d_m, _, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], lat, lon)
    dh = d_m / 1000.0
    dz = dep - MAIN["depth_km"]
    return math.sqrt(dh ** 2 + dz ** 2), dh, dz


def consultar(provider):
    c = Client(provider, timeout=60)
    cat = c.get_events(
        starttime=T_INI, endtime=T_FIN,
        latitude=MAIN["lat"], longitude=MAIN["lon"], maxradius=MAXRADIUS_DEG,
        mindepth=DEPTH_MIN, maxdepth=DEPTH_MAX,
        minmagnitude=MAG_MIN, maxmagnitude=MAG_MAX,
    )
    rows = []
    for e in cat:
        o = e.preferred_origin() or e.origins[0]
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        if m is None:
            continue
        dep = (o.depth or 0.0) / 1000.0
        d3, dh, dz = dist3d_km(o.latitude, o.longitude, dep)
        dm = MAIN["mw"] - m.mag
        # Puntaje: distancia 3D domina; penaliza alejarse del dM ideal.
        score = d3 + 15.0 * abs(dm - DM_IDEAL)
        rows.append({
            "provider": provider,
            "time": str(o.time),
            "lat": round(o.latitude, 4), "lon": round(o.longitude, 4),
            "depth_km": round(dep, 1),
            "mag": m.mag, "dM": round(dm, 2),
            "dist3d_km": round(d3, 1), "dist_h_km": round(dh, 1),
            "dz_km": round(dz, 1),
            "score": round(score, 1),
        })
    return rows


def main():
    rows = []
    for prov in ("USGS", "GFZ"):
        try:
            r = consultar(prov)
            print(f"[{prov}] {len(r)} candidatos")
            rows.extend(r)
        except Exception as ex:  # noqa: BLE001
            print(f"[{prov}] fallo: {ex}")

    # Dedup aproximado (mismo evento en ambos catalogos): por tiempo +-10 s.
    rows.sort(key=lambda r: r["score"])
    unicos = []
    for r in rows:
        t = UTCDateTime(r["time"])
        if any(abs(t - UTCDateTime(u["time"])) < 10.0 for u in unicos):
            continue
        unicos.append(r)

    out = os.path.join(BASE, "candidatos_ranking.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(unicos[0].keys()))
        w.writeheader()
        w.writerows(unicos)

    print(f"\nRanking ({len(unicos)} unicos) -> {out}")
    print(f"{'time':26s} {'M':>4s} {'dM':>5s} {'d3D':>6s} {'dh':>5s} {'dz':>6s} score")
    for r in unicos[:10]:
        print(f"{r['time']:26s} {r['mag']:4.1f} {r['dM']:5.2f} "
              f"{r['dist3d_km']:6.1f} {r['dist_h_km']:5.1f} {r['dz_km']:6.1f} {r['score']}")


if __name__ == "__main__":
    main()
