"""
map.py
======
Mapa del evento de Calama (2026-05-25) adaptado desde el map.py original de
Diego de Almagro. Dibuja:

  * Relieve regional del Norte de Chile (tierra + mar).
  * Costa y fronteras.
  * Estaciones usadas (leidas de los SAC de aceleracion) con su forma de onda
    (componente Z) superpuesta y normalizada, igual que en el original.
  * Epicentro y mecanismo focal tomados de ``event.ctl``.
  * (Opcional) contornos de profundidad de la losa Slab2, si esta instalado
    ``rockhound``.

DEPENDENCIAS (no estan en el env actual; instalar antes de correr):
    conda install -c conda-forge pygmt
    pip install rockhound        # opcional, solo para los contornos de slab

Uso:
    python codigos/map.py
"""
# %%

import glob
import os
import sys

import numpy as np
import obspy

# sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pygmt  # noqa: E402
from event import load_event  # noqa: E402
from pygmt.datasets import load_earth_relief  # noqa: E402

# %%

# ==============================================================================
# 0. Parametros del evento y rutas
# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ev = load_event(os.path.join(BASE, "event.ctl"))

lon_event = ev["lon"]
lat_event = ev["lat"]
depth_event = ev["depth"]

MAGNITUDE = 6.9  # ajusta si tienes el valor oficial

# Mecanismo focal en convencion Aki (strike/dip/rake) desde event.ctl.
if None in (ev["strike"], ev["dip"], ev["rake"]):
    raise SystemExit("Faltan strike/dip/rake en event.ctl para el mecanismo focal.")
focal_mechanism = {
    "strike": ev["strike"],
    "dip": ev["dip"],
    "rake": ev["rake"],
    "magnitude": MAGNITUDE,
}

# Estaciones usadas en la inversion (se dibujan en ROJO; el resto en azul).
USED_STATIONS = {"A22F", "PB19", "PB15", "A20F", "A06F", "PB03", "PB09", "PB01"}

# SAC de aceleracion (generados por txt_to_sac.py).
SAC_ACC = os.path.join(BASE, "SAC", "ACC", "*.acc.sac")

# Ventana de la forma de onda a dibujar (s desde el origen).
window_length = 120.0


# ==============================================================================
# 1. Leer estaciones desde los SAC y definir la region
# ==============================================================================
stream = obspy.read(SAC_ACC)
stream.trim(
    starttime=ev["origin_time"],
    endtime=ev["origin_time"] + window_length,
    pad=True,
    fill_value=0,
)

stations = {}
for tr in stream:
    if tr.stats.channel.endswith("Z"):  # una traza (Z) por estacion
        stations[tr.stats.station] = {
            "lat": tr.stats.sac.stla,
            "lon": tr.stats.sac.stlo,
            "data": tr.data,
        }

lons = [s["lon"] for s in stations.values()] + [lon_event]
lats = [s["lat"] for s in stations.values()] + [lat_event]
region = [
    float(min(lons)) - 0.7,
    float(max(lons)) + 0.7,
    float(min(lats)) - 0.7,
    float(max(lats)) + 0.7,
]
proj = "M6i"

max_value = max(np.max(np.abs(s["data"])) for s in stations.values())
print(
    f"Region: {region} | {len(stations)} estaciones | max|acc| = {max_value:.3e} m/s2"
)


# ==============================================================================
# 2. Mapa base (relieve + costa)
# ==============================================================================
fig = pygmt.Figure()
grid = load_earth_relief(resolution="15s", region=region)
grid_land = grid.where(grid >= 0, other=np.nan)
grid_sea = grid.where(grid < 0, other=np.nan)

fig.grdimage(
    grid=grid_land,
    region=region,
    projection=proj,
    frame=True,
    shading=True,
    cmap="grayC",
    nan_transparent=True,
)
fig.grdimage(
    grid=grid_sea,
    region=region,
    projection=proj,
    frame=True,
    shading=False,
    cmap="ocean",
    nan_transparent=True,
)
fig.coast(
    region=region,
    projection=proj,
    shorelines="1p",
    resolution="f",
    borders=["1/1p,black"],
)
fig.basemap(
    region=region,
    projection=proj,
    frame=[
        f"WSne+t{ev['origin_time'].date} {MAGNITUDE} Mw "
        f"{depth_event:.0f} km - Calama, Chile"
    ],
)


# ==============================================================================
# 3. (Opcional) Contornos de la losa Slab2
# ==============================================================================
PLOT_SLAB2 = False  # True para dibujar los contornos de profundidad de la losa
try:
    if not PLOT_SLAB2:
        raise ImportError("Slab2 desactivado (PLOT_SLAB2=False)")
    import rockhound

    slab2 = rockhound.fetch_slab2(zone="south_america")
    slab2["longitude"] = slab2["longitude"] - 360
    slab2["depth"] = slab2["depth"] / 1000.0
    # Recortar el grid a la region del mapa antes de contornear.
    slab2_reg = slab2.depth.sel(
        longitude=slice(region[0], region[1]),
        latitude=slice(region[2], region[3]),
    )
    fig.grdcontour(
        grid=slab2_reg,
        region=region,
        projection=proj,
        levels=20,
        annotation="40+e+f8p+gwhite",
        pen="0.5p,gray",
    )
except ImportError:
    print("[map.py] rockhound no instalado: se omiten los contornos de Slab2.")
except Exception as exc:  # noqa: BLE001  (p.ej. fallo de descarga sin internet)
    print(f"[map.py] No se pudo cargar Slab2 ({exc}): se omiten los contornos.")


# ==============================================================================
# 4. Estaciones + formas de onda superpuestas
# ==============================================================================
total_length = 1.5  # largo (en grados de lon) que ocupa la forma de onda
amp_scale = 0.4  # amplitud maxima (en grados de lat)

for name, s in stations.items():
    lon, lat, data = s["lon"], s["lat"], s["data"]
    color = "red" if name in USED_STATIONS else "blue"
    norm = data / max_value * amp_scale
    t_x = np.arange(len(data)) / (len(data) - 1) * total_length
    fig.plot(x=t_x + lon, y=norm + lat, pen=f"0.5p,{color}")
    fig.plot(x=lon, y=lat, style="t0.4c", pen="black", fill=color)
    fig.text(
        x=lon,
        y=lat,
        text=name,
        font="8p,Helvetica-Bold,black",
        offset="0c/-0.3c",
        fill="white",
    )

# Escala de la forma de onda (esquina superior izquierda).
x0, y0 = region[0] + 0.3, region[3] - 0.8
fig.plot(x=[x0, x0], y=[y0, y0 + amp_scale], pen="2p,blue")
fig.text(
    x=x0 + 0.05,
    y=y0 + amp_scale,
    text=f"{100 * max_value / 9.8:.1f}% g",
    font="10p,Helvetica-Bold,white",
    justify="BL",
    fill="black",
)
fig.plot(x=[x0, x0 + total_length], y=[y0, y0], pen="2p,blue")
fig.text(
    x=x0 + total_length * 0.8,
    y=y0 - 0.05,
    text=f"{window_length:.0f} s",
    font="10p,Helvetica-Bold,white",
    justify="TL",
    fill="black",
)


# ==============================================================================
# 5. Mecanismo focal y epicentro
# ==============================================================================
fig.meca(
    spec=focal_mechanism,
    scale="1.5c",
    longitude=lon_event,
    latitude=lat_event,
    depth=depth_event,
    plot_longitude=lon_event - 1.0,
    plot_latitude=lat_event + 1.0,
    compressionfill="red",
    extensionfill="cornsilk",
    pen="0.5p,gray30,solid",
    offset="+p0.5p,black",
)
fig.plot(x=lon_event, y=lat_event, style="a0.8c", pen="black", fill="yellow")

# ==============================================================================
# 6. Inset global
# ==============================================================================
with fig.inset(position="jBR+w3.5c+o0.2c", margin=0, box="+p1p,gold"):
    fig.coast(
        region="g",
        projection=f"G{lon_event}/{lat_event}/?",
        land="gray",
        water="white",
        shorelines="0.2p,black",
        frame="g",
    )
    fig.plot(
        x=lon_event,
        y=lat_event,
        style="a0.4c",
        fill="red",
        pen="0.5p,black",
        projection=f"G{lon_event}/{lat_event}/?",
    )

# ==============================================================================
# 7. Subfigura a la derecha: zoom con circulo de 200 km
#    Solo epicentro, estaciones y circulo (sin formas de onda).
# ==============================================================================
RADIO_KM = 200.0
ZOOM_WIDTH_CM = 12.0
MAIN_WIDTH_CM = 6 * 2.54  # proj "M6i" -> 6 pulgadas

# Margen para que el circulo de 200 km entre completo (1 deg lat ~ 111 km).
margen = RADIO_KM / 111.0 + 0.4
zoom_region = [
    lon_event - margen,
    lon_event + margen,
    lat_event - margen,
    lat_event + margen,
]
zoom_proj = f"M{ZOOM_WIDTH_CM}c"


def mercator_height_cm(reg, width_cm):
    """Altura (cm) de un mapa Mercator dada su region y su ancho."""
    lon0, lon1, lat0, lat1 = reg
    merc = lambda p: np.degrees(np.log(np.tan(np.pi / 4 + np.radians(p) / 2)))
    scale = width_cm / (lon1 - lon0)  # cm por grado de longitud
    return scale * (merc(lat1) - merc(lat0))


# Posicion del panel: a la derecha del mapa principal y alineado por arriba.
# No usamos 'w'/'h' porque tras el inset global esas referencias cambian; en su
# lugar calculamos las dimensiones reales en cm.
main_h = mercator_height_cm(region, MAIN_WIDTH_CM)
zoom_h = mercator_height_cm(zoom_region, ZOOM_WIDTH_CM)
fig.shift_origin(xshift=f"{MAIN_WIDTH_CM + 1.5}c", yshift=f"{main_h - zoom_h}c")

fig.coast(
    region=zoom_region,
    projection=zoom_proj,
    land="gray80",
    water="lightblue",
    shorelines="0.5p,black",
    borders=["1/0.8p,black"],
    frame=[f"WSne+tZoom: radio {RADIO_KM:.0f} km", "xa1f0.5", "ya1f0.5"],
)

# Circulo geografico de radio 200 km (style 'E-' usa el DIAMETRO en km).
fig.plot(
    x=[lon_event],
    y=[lat_event],
    style="E-",
    size=[2 * RADIO_KM],
    pen="1.5p,red",
    projection=zoom_proj,
    region=zoom_region,
)

# Estaciones (solo el marcador + nombre, sin trazas).
for name, s in stations.items():
    fig.plot(
        x=s["lon"],
        y=s["lat"],
        style="t0.35c",
        fill="red" if name in USED_STATIONS else "blue",
        pen="0.4p,black",
        projection=zoom_proj,
        region=zoom_region,
    )
    fig.text(
        x=s["lon"],
        y=s["lat"],
        text=name,
        font="6p,Helvetica-Bold,black",
        offset="0c/-0.25c",
        fill="white",
        projection=zoom_proj,
        region=zoom_region,
    )

# Epicentro.
fig.plot(
    x=lon_event,
    y=lat_event,
    style="a0.7c",
    fill="yellow",
    pen="black",
    projection=zoom_proj,
    region=zoom_region,
)


# ==============================================================================
# 8. Panel de trazas de desplazamiento (record section), debajo del zoom.
#    Componente Z de las estaciones a < RADIO_KM; seleccionadas en rojo.
#    Eje X = tiempo desde el origen; linea roja punteada en t = 0.
# ==============================================================================
WF_WIDTH_CM = ZOOM_WIDTH_CM
WF_HEIGHT_CM = 18.0
TMIN, TMAX = -10.0, 130.0          # ventana de tiempo [s] respecto al origen

# Cargar desplazamiento (Z) de las estaciones cercanas, ordenadas por distancia.
disp = []
for f in glob.glob(os.path.join(BASE, "SAC", "DISP", "*.disp.sac")):
    tr = obspy.read(f)[0]
    if tr.stats.channel.endswith("Z") and tr.stats.sac.user0 < RADIO_KM:
        disp.append((tr.stats.sac.user0, tr.stats.station, tr))
disp.sort(key=lambda r: r[0])      # mas cercana primero
n = len(disp)

# Bajar el origen para quedar debajo del panel de zoom.
fig.shift_origin(yshift=f"-{WF_HEIGHT_CM + 1.5}c")

wf_region = [TMIN, TMAX, -1, n]
wf_proj = f"X{WF_WIDTH_CM}c/{WF_HEIGHT_CM}c"
fig.basemap(region=wf_region, projection=wf_proj,
            frame=["WSne+tDesplazamiento Z (rojo = usadas)",
                   "xa20f10+lTiempo desde el origen [s]", "ya100"])

# Linea vertical en el tiempo de origen.
fig.plot(x=[0, 0], y=[-1, n], pen="0.8p,red,dashed",
         region=wf_region, projection=wf_proj)

for i, (dist, name, tr) in enumerate(disp):
    offset = n - 1 - i                       # la mas cercana arriba
    dec = max(1, int(round(tr.stats.sampling_rate / 20.0)))  # ~20 Hz para dibujar
    t = tr.times() - (ev["origin_time"] - tr.stats.starttime)
    data = tr.data
    amax = np.max(np.abs(data)) or 1.0
    y = data / amax * 0.45 + offset          # normalizada, semi-espaciado 0.45
    color = "red" if name in USED_STATIONS else "black"
    fig.plot(x=t[::dec], y=y[::dec], pen=f"0.6p,{color}",
             region=wf_region, projection=wf_proj)
    fig.text(x=TMIN, y=offset, text=f"{name} ({dist:.0f} km)",
             font=f"7p,Helvetica-Bold,{color}", justify="LM",
             offset="0.1c/0.18c", region=wf_region, projection=wf_proj)


out = os.path.join(BASE, "SAC", "mapa_calama.png")
fig.savefig(out, dpi=300)
print(f"Mapa guardado en {out}")
fig.show()
