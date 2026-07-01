"""
comparar_cobertura_waveforms.py
===============================
Compara los eventos Calama 2020-06-03 y 2026-05-25:

  (1) Cobertura azimutal de las estaciones de cada inversión, en un único
      diagrama polar central. Cada estación se colorea por la POLARIDAD
      teórica del primer movimiento P (patrón de radiación de doble cupla,
      Aki & Richards 2002 ec. 4.84), usando el ángulo de salida i (TauP,
      ak135) y el mecanismo strike/dip/rake de cada evento:
          rojo  = compresión  (F_P > 0)
          azul  = dilatación  (F_P < 0)
      Forma del marcador = evento:  ▲ 2026   ▼ 2020.
      Borde del marcador = color del evento (orangered 2026 / dodgerblue 2020).

  (2) Alrededor del polar, en su azimut, las formas de onda en DESPLAZAMIENTO
      (componente Z) de los ~5 pares de estaciones co-localizadas (<6 km),
      superponiendo 2026 (orangered) y 2020 (dodgerblue), normalizadas.

Uso:
    conda run -n kdellipspy python codigos/comparar_cobertura_waveforms.py
"""

import glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

# Ventana de las trazas de desplazamiento (s desde el origen), igual que map.py.
TMIN, TMAX = -10.0, 130.0

# ---------------------------------------------------------------------------- #
BASE = Path(__file__).resolve().parent.parent          # comparacion_2020_vs_2026/
INVERSIONS = BASE.parent
OUT = BASE / "output" / "cobertura_waveforms.png"

EVENTS = {
    "2020": dict(dir=INVERSIONS / "2020-06-03_Calama", color="dodgerblue", marker="v"),
    "2026": dict(dir=INVERSIONS / "2026-05-25_Calama", color="orangered",  marker="^"),
}

# Pares co-localizados (<6 km): (estación_2026, estación_2020). Ver análisis de
# vecino más cercano; el resto de estaciones no tiene contraparte cercana.
PAIRS = [("PB06", "PB06"), ("A06F", "A06F"), ("PB19", "PB19"),
         ("A20F", "A18F"), ("AF01", "A22F")]

MODEL = TauPyModel(model="ak135")


# ---------------------------------------------------------------------------- #
def load_event(d):
    """event.ctl -> dict lat lon depth strike dip rake."""
    with open(d / "event.ctl") as fh:
        for raw in fh:
            raw = raw.strip()
            if raw and not raw.startswith("#"):
                p = raw.split()
                return dict(lat=float(p[0]), lon=float(p[1]), depth=float(p[2]),
                            strike=float(p[3]), dip=float(p[4]), rake=float(p[5]),
                            origin=UTCDateTime(p[6].replace("-T", "T")))
    raise ValueError(f"sin línea de datos en {d/'event.ctl'}")


def load_stations(d):
    """Sección 8 de input.ctl -> lista [(name, lat, lon)] en orden (= orden DATA)."""
    lines = (d / "input.ctl").read_text().splitlines()
    n, out, i = None, [], 0
    while i < len(lines):
        if "Number of stations" in lines[i]:
            n = int(lines[i].split(":")[1])
            i += 1
            while len(out) < n and i < len(lines):
                p = lines[i].split()
                if len(p) >= 4 and not lines[i].lstrip().startswith("#"):
                    out.append((p[3], float(p[0]), float(p[1])))
                i += 1
            break
        i += 1
    return out


def load_disp_trace(d, station, origin):
    """Traza de desplazamiento Z (banda ancha, SIN filtrar) de SAC/DISP, recortada
    a [TMIN, TMAX] respecto al origen. Igual fuente que map.py. -> (t[s], data[mm])."""
    fs = glob.glob(str(d / "SAC" / "DISP" / f"{station}*Z.disp.sac"))
    if not fs:
        return None
    tr = obspy.read(fs[0])[0]
    tr.trim(origin + TMIN, origin + TMAX, pad=True, fill_value=0.0)
    return tr.times() + TMIN, tr.data * 1e3


def fp_radiation(strike, dip, rake, az, takeoff):
    """Patrón de radiación P (Aki & Richards 2002, ec. 4.84). Grados.
    az = azimut fuente→estación; takeoff = ángulo de salida desde vertical hacia abajo."""
    s, d, l = map(np.radians, (strike, dip, rake))
    f = np.radians(az) - s
    i = np.radians(takeoff)
    return (np.cos(l)*np.sin(d)*np.sin(i)**2*np.sin(2*f)
            - np.cos(l)*np.cos(d)*np.sin(2*i)*np.cos(f)
            + np.sin(l)*np.sin(2*d)*(np.cos(i)**2 - np.sin(i)**2*np.sin(f)**2)
            + np.sin(l)*np.cos(2*d)*np.sin(2*i)*np.sin(f))


def takeoff_angle(depth_km, dist_km):
    """Ángulo de salida de P (TauP ak135) desde la vertical hacia abajo, en grados."""
    arr = MODEL.get_travel_times(source_depth_in_km=depth_km,
                                 distance_in_degree=kilometers2degrees(dist_km),
                                 phase_list=["P", "p"])
    return arr[0].takeoff_angle if arr else np.nan


def radiation_field(ax, ev, rmax, alpha):
    """Sombrea el fondo polar con el patrón de radiación P del evento:
    rojo=compresión (F_P>0), azul=dilatación (F_P<0). El takeoff depende del
    radio (distancia) vía TauP; la polaridad varía con el azimut."""
    rr = np.linspace(3.0, rmax, 40)
    th = np.linspace(0.0, 2*np.pi, 181)
    to = np.array([takeoff_angle(ev["depth"], r) for r in rr])     # i(r)
    TH, RR = np.meshgrid(th, rr)
    I = np.repeat(to[:, None], th.size, axis=1)
    F = fp_radiation(ev["strike"], ev["dip"], ev["rake"], np.degrees(TH), I)
    m = np.nanmax(np.abs(F)) or 1.0
    return ax.pcolormesh(TH, RR, F, cmap="RdBu_r", vmin=-m, vmax=m,
                         alpha=alpha, shading="auto", zorder=0)


def station_table(ev, stations):
    """Por estación: az, dist(km), takeoff, polaridad(+1 compresión / -1 dilatación)."""
    rows = {}
    for name, la, lo in stations:
        dm, az, _ = gps2dist_azimuth(ev["lat"], ev["lon"], la, lo)
        dk = dm / 1000.0
        to = takeoff_angle(ev["depth"], dk)
        fp = fp_radiation(ev["strike"], ev["dip"], ev["rake"], az, to)
        rows[name] = dict(az=az, dist=dk, takeoff=to, pol=int(np.sign(fp)))
    return rows


# ---------------------------------------------------------------------------- #
def main():
    evs, sts, tab = {}, {}, {}
    for k, meta in EVENTS.items():
        evs[k] = {**meta, **load_event(meta["dir"])}
        sts[k] = load_stations(meta["dir"])
        tab[k] = station_table(evs[k], sts[k])

    fig = plt.figure(figsize=(16, 16))
    # --- polar central de cobertura ---------------------------------------- #
    ax = fig.add_axes([0.32, 0.32, 0.36, 0.36], projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    rmax = max(r["dist"] for k in tab for r in tab[k].values()) * 1.18
    ax.set_ylim(0, rmax)
    # Fondo = campo de radiación P (2026 más fuerte, 2020 pastel; intersección oscura).
    radiation_field(ax, evs["2026"], rmax, alpha=0.55)
    radiation_field(ax, evs["2020"], rmax, alpha=0.35)
    for k in EVENTS:
        ecol = EVENTS[k]["color"]
        for name, r in tab[k].items():
            t = np.deg2rad(r["az"])
            ax.plot([t, t], [0, r["dist"]], color="0.5", lw=0.7, zorder=1)
            ax.scatter(t, r["dist"], marker=EVENTS[k]["marker"], s=120,
                       facecolor="white", edgecolor=ecol, linewidth=1.8, zorder=3)
            ax.text(t, r["dist"] + rmax*0.045, f"{name}\n{r['dist']:.0f} km",
                    ha="center", va="center", fontsize=7, color=ecol, zorder=4)
    ax.plot(0, 0, marker="*", ms=16, color="yellow", markeredgecolor="k", zorder=5)
    ax.set_rlabel_position(135)
    ax.set_title("Cobertura azimutal · fondo = radiación P (rojo=compresión, azul=dilatación)\n"
                 "▲ 2026 (borde naranja)   ▼ 2020 (borde azul)", fontsize=11, pad=22)

    # --- waveforms de los pares: anillo equiespaciado (ordenado por azimut) + líneas guía
    R, w, h = 0.42, 0.18, 0.12
    order = sorted(PAIRS, key=lambda p: tab["2026"][p[0]]["az"])
    for j, (s26, s20) in enumerate(order):
        ang = np.deg2rad(90 - j/len(order)*360.0)       # equiespaciado, sin solape
        cx, cy = 0.5 + R*np.cos(ang), 0.5 + R*np.sin(ang)
        a = fig.add_axes([cx - w/2, cy - h/2, w, h])
        for k, sname in (("2026", s26), ("2020", s20)):
            tr = load_disp_trace(evs[k]["dir"], sname, evs[k]["origin"])
            if tr is None:
                continue
            t, z = tr                                   # mm, banda ancha SIN filtrar
            a.plot(t, z, color=EVENTS[k]["color"], lw=1.0,
                   label=f"{sname} {k}  (pk {np.abs(z).max():.1f} mm)")
        a.set_xlim(TMIN, TMAX)
        a.set_title(f"{s26} ↔ {s20}", fontsize=8)
        a.tick_params(labelsize=6); a.set_ylabel("Z (mm)", fontsize=6)
        a.legend(fontsize=6, loc="upper right", framealpha=0.6)
        a.axhline(0, color="0.8", lw=0.5)
        # línea guía panel -> estación real en el radar
        t = np.deg2rad(tab["2026"][s26]["az"]); rdist = tab["2026"][s26]["dist"]
        xf, yf = fig.transFigure.inverted().transform(ax.transData.transform((t, rdist)))
        fig.add_artist(plt.Line2D([cx, xf], [cy, yf], color="0.6", lw=0.6,
                                  zorder=0, transform=fig.transFigure))

    fig.suptitle("Calama 2020-06-03  vs  2026-05-25 — cobertura azimutal y desplazamiento Z",
                 fontsize=14, y=0.95)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=140)
    print("figura ->", OUT)

    # self-check: pares realmente co-localizados y polaridades binarias
    for s26, s20 in PAIRS:
        d26 = dict((n, (la, lo)) for n, la, lo in sts["2026"])[s26]
        d20 = dict((n, (la, lo)) for n, la, lo in sts["2020"])[s20]
        dk = gps2dist_azimuth(*d26, *d20)[0] / 1000.0
        assert dk < 6.0, f"par {s26}-{s20} a {dk:.1f} km (>6)"
    assert all(r["pol"] in (-1, 1) for k in tab for r in tab[k].values())
    print("self-check OK")


if __name__ == "__main__":
    main()
