#!/usr/bin/env python3
"""
kde-map · Mapa del evento desde input.ctl + SAC/

Si pygmt está instalado, dibuja el mapa completo de 3 paneles (relieve + formas
de onda, zoom con círculo de 200 km y record section). Si no, cae al mapa simple
con matplotlib + obspy (estaciones + epicentro + beachball; costa si hay cartopy).

Estaciones del input.ctl = usadas (ROJO); el resto de las estaciones con SAC, azul.

    kde-map                       # ./input.ctl + ./SAC → ./SAC/mapa_calama.png
    kde-map ruta/al/proyecto --mag 6.9
"""
import argparse
import glob
from math import cos, radians
from pathlib import Path

from .cli import banner, section, info, ok, warn, err


def _sac_stations(sac_dir: Path):
    """{nombre: (lat, lon)} de la primera subcarpeta SAC con .sac (componente Z)."""
    from obspy import read
    for sub in ("ACC", "VEL", "DISP"):
        files = glob.glob(str(sac_dir / sub / "*.sac"))
        if files:
            out = {}
            for f in files:
                tr = read(f, headonly=True)[0]
                if tr.stats.channel.endswith("Z"):
                    out[tr.stats.station] = (tr.stats.sac.stla, tr.stats.sac.stlo)
            return out, sub
    return {}, None


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="kde-map", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("project_dir", type=Path, nargs="?", default=Path.cwd())
    p.add_argument("--input-ctl", type=Path, default=None)
    p.add_argument("--sac-dir", type=Path, default=None)
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--mag", type=float, default=6.9,
                   help="Magnitud (no está en el ctl): tamaño del beachball / título")
    p.add_argument("--simple", action="store_true",
                   help="Fuerza el mapa simple aunque pygmt esté instalado")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    import kdellipspy as kde

    banner("KDEllipsPy · Mapa del evento", "input.ctl + SAC → mapa_calama.png")

    ctl = args.input_ctl or args.project_dir / "input.ctl"
    sac = args.sac_dir or args.project_dir / "SAC"
    out = args.output or sac / "mapa_calama.png"

    section("Validando rutas")
    if not ctl.is_file():
        err(f"No existe input.ctl: {ctl}"); return 1
    ok(f"input.ctl : {ctl}")

    cfg = kde.ConfigParser(filepath=str(ctl))
    sp = cfg.source_position
    used = {s.name: (s.latitude, s.longitude) for s in cfg.stations.stations}

    section("Configuración")
    info(f"Evento     : {sp.event_name}  ({sp.latitude:.3f}, {sp.longitude:.3f}, {sp.depth:.0f} km)")
    info(f"Mecanismo  : strike {sp.strike:g} / dip {sp.dip:g} / rake {sp.rake:g}")
    info(f"Usadas ({len(used)}): {', '.join(used)}")

    have_pygmt = False
    if not args.simple:
        try:
            import pygmt  # noqa: F401
            have_pygmt = True
        except ImportError:
            warn("pygmt no instalado → mapa simple (instala pygmt para los 3 paneles)")

    out.parent.mkdir(parents=True, exist_ok=True)
    if have_pygmt:
        return _full_map(cfg, used, sac, out, args.mag)
    return _simple_map(cfg, used, sac, out, args.mag)


# --------------------------------------------------------------------------- #
#  Mapa simple (matplotlib + obspy; costa opcional con cartopy)
# --------------------------------------------------------------------------- #
def _simple_map(cfg, used, sac, out, mag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from obspy.imaging.beachball import beach

    sp = cfg.source_position
    sac_sta, sub = _sac_stations(sac) if sac.is_dir() else ({}, None)
    stations = {**sac_sta, **used}  # ctl manda en coordenadas
    if not stations:
        err("Sin estaciones (ni en input.ctl ni en SAC)"); return 1
    if sub:
        info(f"Contexto SAC/{sub}: {len(sac_sta)} estaciones")

    section("Dibujando (simple)")
    lons = [lo for _, lo in stations.values()] + [sp.longitude]
    lats = [la for la, _ in stations.values()] + [sp.latitude]
    pad = 0.6
    extent = [min(lons) - pad, max(lons) + pad, min(lats) - pad, max(lats) + pad]

    ax = None
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        fig = plt.figure(figsize=(8, 9))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke")
        ax.add_feature(cfeature.OCEAN, facecolor="lightsteelblue")
        ax.add_feature(cfeature.COASTLINE, edgecolor="gray", linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linestyle=":", edgecolor="gray")
        gl = ax.gridlines(draw_labels=True, alpha=0.3, linestyle="--")
        gl.top_labels = gl.right_labels = False
        tx = {"transform": ccrs.PlateCarree()}
    except ImportError:
        warn("cartopy no instalado: mapa sin costa (puntos sobre ejes lon/lat)")
        fig, ax = plt.subplots(figsize=(8, 9))
        ax.set_aspect(1.0 / cos(radians(sp.latitude)))
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud"); ax.grid(alpha=0.3)
        tx = {}

    for name, (la, lo) in stations.items():
        red = name in used
        ax.plot(lo, la, "^", ms=10, color="red" if red else "royalblue",
                mec="black", mew=0.6, zorder=5, **tx)
        ax.text(lo, la, f" {name}", fontsize=8, fontweight="bold",
                va="bottom", color="red" if red else "royalblue", zorder=6, **tx)

    bb = beach([sp.strike, sp.dip, sp.rake], xy=(sp.longitude, sp.latitude),
               width=0.5, linewidth=0.6, facecolor="red", zorder=8, axes=ax)
    ax.add_collection(bb)
    ax.plot(sp.longitude, sp.latitude, "*", ms=16, color="yellow",
            mec="black", zorder=10, **tx)

    date = str(sp.origin_time)[:10] if sp.origin_time else ""
    ax.set_title(f"{date}  Mw {mag:g}  {sp.depth:.0f} km — {sp.event_name}",
                 fontsize=12, fontweight="bold")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    ok(f"Mapa simple guardado → {out}")
    return 0


# --------------------------------------------------------------------------- #
#  Mapa completo de 3 paneles (pygmt). Porteado del map.py original, leyendo el
#  evento y las estaciones usadas desde input.ctl en vez de event.ctl/event.py.
# --------------------------------------------------------------------------- #
def _full_map(cfg, used, sac, out, mag):
    import numpy as np
    import obspy
    import pygmt
    from obspy import UTCDateTime
    from pygmt.datasets import load_earth_relief

    sp = cfg.source_position
    lon_event, lat_event, depth_event = sp.longitude, sp.latitude, sp.depth
    origin_time = UTCDateTime(sp.origin_time)
    used_names = set(used)
    focal = {"strike": sp.strike, "dip": sp.dip, "rake": sp.rake, "magnitude": mag}

    section("Dibujando (pygmt, 3 paneles)")
    window_length = 120.0

    # --- 1. Estaciones (ACC, componente Z) + región -----------------------
    stream = obspy.read(str(sac / "ACC" / "*.acc.sac"))
    stream.trim(origin_time, origin_time + window_length, pad=True, fill_value=0)
    stations = {}
    for tr in stream:
        if tr.stats.channel.endswith("Z"):
            stations[tr.stats.station] = {
                "lat": tr.stats.sac.stla, "lon": tr.stats.sac.stlo, "data": tr.data}

    lons = [s["lon"] for s in stations.values()] + [lon_event]
    lats = [s["lat"] for s in stations.values()] + [lat_event]
    region = [float(min(lons)) - 0.7, float(max(lons)) + 0.7,
              float(min(lats)) - 0.7, float(max(lats)) + 0.7]
    proj = "M6i"
    max_value = max(np.max(np.abs(s["data"])) for s in stations.values())
    info(f"Región {region} | {len(stations)} estaciones | max|acc|={max_value:.2e} m/s²")

    # --- 2. Base: relieve + costa -----------------------------------------
    fig = pygmt.Figure()
    grid = load_earth_relief(resolution="15s", region=region)
    fig.grdimage(grid=grid.where(grid >= 0, other=np.nan), region=region, projection=proj,
                 frame=True, shading=True, cmap="grayC", nan_transparent=True)
    fig.grdimage(grid=grid.where(grid < 0, other=np.nan), region=region, projection=proj,
                 frame=True, shading=False, cmap="ocean", nan_transparent=True)
    fig.coast(region=region, projection=proj, shorelines="1p", resolution="f",
              borders=["1/1p,black"])
    fig.basemap(region=region, projection=proj,
                frame=[f"WSne+t{origin_time.date} {mag} Mw {depth_event:.0f} km - Calama, Chile"])

    # --- 3. Estaciones + formas de onda -----------------------------------
    total_length, amp_scale = 1.5, 0.4
    for name, s in stations.items():
        color = "red" if name in used_names else "blue"
        norm = s["data"] / max_value * amp_scale
        t_x = np.arange(len(s["data"])) / (len(s["data"]) - 1) * total_length
        fig.plot(x=t_x + s["lon"], y=norm + s["lat"], pen=f"0.5p,{color}")
        fig.plot(x=s["lon"], y=s["lat"], style="t0.4c", pen="black", fill=color)
        fig.text(x=s["lon"], y=s["lat"], text=name, font="8p,Helvetica-Bold,black",
                 offset="0c/-0.3c", fill="white")

    x0, y0 = region[0] + 0.3, region[3] - 0.8
    fig.plot(x=[x0, x0], y=[y0, y0 + amp_scale], pen="2p,blue")
    fig.text(x=x0 + 0.05, y=y0 + amp_scale, text=f"{100 * max_value / 9.8:.1f}% g",
             font="10p,Helvetica-Bold,white", justify="BL", fill="black")
    fig.plot(x=[x0, x0 + total_length], y=[y0, y0], pen="2p,blue")
    fig.text(x=x0 + total_length * 0.8, y=y0 - 0.05, text=f"{window_length:.0f} s",
             font="10p,Helvetica-Bold,white", justify="TL", fill="black")

    # --- 4. Mecanismo focal + epicentro -----------------------------------
    fig.meca(spec=focal, scale="1.5c", longitude=lon_event, latitude=lat_event,
             depth=depth_event, plot_longitude=lon_event - 1.0, plot_latitude=lat_event + 1.0,
             compressionfill="red", extensionfill="cornsilk", pen="0.5p,gray30,solid",
             offset="+p0.5p,black")
    fig.plot(x=lon_event, y=lat_event, style="a0.8c", pen="black", fill="yellow")

    # --- 5. Inset global ---------------------------------------------------
    with fig.inset(position="jBR+w3.5c+o0.2c", margin=0, box="+p1p,gold"):
        fig.coast(region="g", projection=f"G{lon_event}/{lat_event}/?", land="gray",
                  water="white", shorelines="0.2p,black", frame="g")
        fig.plot(x=lon_event, y=lat_event, style="a0.4c", fill="red", pen="0.5p,black",
                 projection=f"G{lon_event}/{lat_event}/?")

    # --- 6. Panel zoom: círculo de 200 km ---------------------------------
    RADIO_KM, ZOOM_WIDTH_CM, MAIN_WIDTH_CM = 200.0, 12.0, 6 * 2.54
    margen = RADIO_KM / 111.0 + 0.4
    zoom_region = [lon_event - margen, lon_event + margen, lat_event - margen, lat_event + margen]
    zoom_proj = f"M{ZOOM_WIDTH_CM}c"

    def merc_h(reg, w):
        lon0, lon1, lat0, lat1 = reg
        m = lambda p: np.degrees(np.log(np.tan(np.pi / 4 + np.radians(p) / 2)))
        return w / (lon1 - lon0) * (m(lat1) - m(lat0))

    fig.shift_origin(xshift=f"{MAIN_WIDTH_CM + 1.5}c",
                     yshift=f"{merc_h(region, MAIN_WIDTH_CM) - merc_h(zoom_region, ZOOM_WIDTH_CM)}c")
    fig.coast(region=zoom_region, projection=zoom_proj, land="gray80", water="lightblue",
              shorelines="0.5p,black", borders=["1/0.8p,black"],
              frame=[f"WSne+tZoom: radio {RADIO_KM:.0f} km", "xa1f0.5", "ya1f0.5"])
    fig.plot(x=[lon_event], y=[lat_event], style="E-", size=[2 * RADIO_KM], pen="1.5p,red",
             projection=zoom_proj, region=zoom_region)
    for name, s in stations.items():
        fig.plot(x=s["lon"], y=s["lat"], style="t0.35c",
                 fill="red" if name in used_names else "blue", pen="0.4p,black",
                 projection=zoom_proj, region=zoom_region)
        fig.text(x=s["lon"], y=s["lat"], text=name, font="6p,Helvetica-Bold,black",
                 offset="0c/-0.25c", fill="white", projection=zoom_proj, region=zoom_region)
    fig.plot(x=lon_event, y=lat_event, style="a0.7c", fill="yellow", pen="black",
             projection=zoom_proj, region=zoom_region)

    # --- 7. Record section (DISP Z, estaciones a < 200 km) ----------------
    WF_W, WF_H, TMIN, TMAX = ZOOM_WIDTH_CM, 18.0, -10.0, 130.0
    disp = []
    for f in glob.glob(str(sac / "DISP" / "*.disp.sac")):
        tr = obspy.read(f)[0]
        if tr.stats.channel.endswith("Z") and tr.stats.sac.user0 < RADIO_KM:
            disp.append((tr.stats.sac.user0, tr.stats.station, tr))
    disp.sort(key=lambda r: r[0])
    n = len(disp)
    fig.shift_origin(yshift=f"-{WF_H + 1.5}c")
    wf_region, wf_proj = [TMIN, TMAX, -1, n], f"X{WF_W}c/{WF_H}c"
    fig.basemap(region=wf_region, projection=wf_proj,
                frame=["WSne+tDesplazamiento Z (rojo = usadas)",
                       "xa20f10+lTiempo desde el origen [s]", "ya100"])
    fig.plot(x=[0, 0], y=[-1, n], pen="0.8p,red,dashed", region=wf_region, projection=wf_proj)
    for i, (dist, name, tr) in enumerate(disp):
        offset = n - 1 - i
        dec = max(1, int(round(tr.stats.sampling_rate / 20.0)))
        t = tr.times() - (origin_time - tr.stats.starttime)
        amax = np.max(np.abs(tr.data)) or 1.0
        y = tr.data / amax * 0.45 + offset
        color = "red" if name in used_names else "black"
        fig.plot(x=t[::dec], y=y[::dec], pen=f"0.6p,{color}", region=wf_region, projection=wf_proj)
        fig.text(x=TMIN, y=offset, text=f"{name} ({dist:.0f} km)",
                 font=f"7p,Helvetica-Bold,{color}", justify="LM", offset="0.1c/0.18c",
                 region=wf_region, projection=wf_proj)

    fig.savefig(str(out), dpi=300)
    ok(f"Mapa completo (3 paneles) guardado → {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
