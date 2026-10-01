"""PyGMT version of the Navidad 2019 / Maule 2010 (Yue et al. 2014) map.
Writes maule_navidad_gmt.{png,pdf}."""

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pygmt  # noqa: E402
import xarray as xr  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

import maule_navidad as mn  # noqa: E402  (también fija DYN_CASE y rutas)

REGION = [-73.8, -71.4, -35.3, -33.7]
W = 13.0
PROJ = f"M{W}c"
PICH, KIN, DYN = "142/1/82", "31/79/191", "black"


def contour_paths(lon, lat, z, level):
    cs = plt.contour(lon, lat, z, levels=[level])
    segs = [p for p in cs.allsegs[0] if len(p) > 1]
    plt.close("all")
    return segs


def main():
    m = mn.read_yue()
    cfg = mn.om.kde.ConfigParser(str(mn.om.fd.CASE_DIR / "input.ctl"))
    mn.om.fd.setup_work_dir(cfg)
    fm, geom = mn.om.subfaults(cfg)
    fp, sp = cfg.fault_plane, cfg.source_position
    la = np.array([sf.lat for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    lo = np.array([sf.lon for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    ss, ds, _ = mn.om.kinematic_slip(cfg, fm)
    skin = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    ss, ds, _ = mn.om.dynamic_slip(cfg, geom, "B")
    sdyn = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    gj = mn.json.loads((mn.NAV / "faults" / "chaf_navidad.geojson").read_text())
    inter = mn.json.loads((mn.NAV / "faults" / "pichilemu_intersection.json").read_text())["line_in_mesh"]

    x = np.arange(REGION[0], REGION[1] + 1e-9, 0.01)
    y = np.arange(REGION[2], REGION[3] + 1e-9, 0.01)
    gx, gy = np.meshgrid(x, y)
    gz = np.clip(griddata(np.c_[m["lon"], m["lat"]], m["slip"], (gx, gy), method="cubic"), 0, None)
    grid = xr.DataArray(gz, coords={"lat": y, "lon": x}, dims=("lat", "lon"))

    fig = pygmt.Figure()
    pygmt.config(FONT_ANNOT_PRIMARY="9p", FONT_LABEL="10p", MAP_FRAME_TYPE="plain")
    relief = pygmt.datasets.load_earth_relief(resolution="15s", region=REGION)
    pygmt.makecpt(cmap="gray70,white", series=[-8000, 3000], continuous=True)
    fig.grdimage(grid=relief, region=REGION, projection=PROJ, cmap=True, shading="+a315+nt0.3")
    pygmt.makecpt(cmap="lajolla", reverse=True, series=[0, 18, 2], continuous=True)
    fig.grdimage(grid=grid, cmap=True, nan_transparent=True, transparency=40)
    fig.grdcontour(grid=grid, levels=2, annotation="4+f7p,Helvetica,darkred", pen="0.4p,darkred")
    fig.coast(shorelines="0.6p,black", frame=["WSne", "af"])
    fig.plot(x=m["lon"], y=m["lat"], style="x0.16c", pen="0.6p,gray30")
    for f in gj["features"]:
        if f["properties"].get("F_name") == "Pichilemu":
            fx, fy = np.array(f["geometry"]["coordinates"]).T
            fig.plot(x=fx, y=fy, pen=f"2.4p,{PICH}")
    fig.plot(x=[p["lon"] for p in inter], y=[p["lat"] for p in inter], pen=f"1.6p,{PICH},2_2")
    for g, col, pen in ((skin, KIN, "1.4p"), (sdyn, DYN, "1.2p")):
        for lev in (0.2, 0.6):
            for seg in contour_paths(lo, la, g, lev * g.max()):
                fig.plot(x=seg[:, 0], y=seg[:, 1], pen=f"{pen if lev > 0.5 else '0.7p'},{col}{',4_2' if col == DYN else ''}")
    for pen in ("1.4p,white", "0.9p,black"):
        fig.plot(x=[sp.longitude], y=[sp.latitude], style="a0.5c", pen=pen, fill="0/255/0" if "black" in pen else None)
    fig.plot(x=[mn.GCMT[1]], y=[mn.GCMT[0]], style="d0.3c", fill="purple", pen="0.6p,white")
    fig.basemap(map_scale="jBR+w30k+o0.4c/0.4c+f+lkm", box="+gwhite@30+p0.5p,black")
    fig.colorbar(position=f"JMR+o0.5c/0c+w8c/0.3c", frame=["xa4f2+lMaule 2010 slip (m)"])
    fig.legend(spec=io.StringIO(
        "N 2\n"
        "S 0.3c a 0.4c 0/255/0 0.9p,black 0.9c Navidad hypocentre (ISC)\n"
        "S 0.3c d 0.28c purple 0.6p,white 0.9c Navidad centroid (GCMT)\n"
        f"S 0.3c - 0.5c - 1.4p,{KIN} 0.9c Navidad kinematic slip (20, 60 %)\n"
        f"S 0.3c - 0.5c - 1.2p,{DYN},4_2 0.9c Navidad dynamic B slip (20, 60 %)\n"
        f"S 0.3c - 0.5c - 2.4p,{PICH} 0.9c Pichilemu fault (CHAF)\n"
        f"S 0.3c - 0.5c - 1.6p,{PICH},2_2 0.9c Pichilemu-interface intersection\n"
        "S 0.3c x 0.16c - 0.6p,gray30 0.9c Maule subfault nodes\n"
        "S 0.3c - 0.5c - 0.4p,darkred 0.9c Maule slip contours (2 m)\n"),
        position=f"JBC+jTC+o0/1.0c+w{W}c", box="+gwhite+p0.5p,black")
    fig.text(position="TL", offset="0.2c/-0.2c", justify="TL", font="8p,Helvetica,black", fill="white@20",
             text="Maule 2010: Yue et al. (2014), interpolated")
    for ext in ("png", "pdf"):
        fig.savefig(str(mn.HERE / f"maule_navidad_gmt.{ext}"), dpi=200)
    print("ok")


if __name__ == "__main__":
    main()
