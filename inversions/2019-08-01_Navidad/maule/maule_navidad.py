"""Navidad 2019 sobre el slip de Maule 2010 (Yue et al., 2014, JGR, doi:10.1002/2014JB011340;
FFM conjunto hr-GPS/telesísmico/GPS/InSAR/tsunami, subfallas de ~40x40 km).

Escribe maule_navidad.png (mapa + slip en función de la profundidad) y maule_navidad.json.
"""

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
os.environ.setdefault("DYN_CASE", "isc_tomo/b0.05-0.30_noM14L")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path[:0] = [str(NAV / "okada"), str(NAV / "dynamic")]
import cartopy.crs as ccrs  # noqa: E402
import cartopy.feature as cfeature  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

import okada_models as om  # noqa: E402

K = 111.19
COLS = "time dur lat lon depth strike dip rake m0 slip x y xinc yinc seg sub par".split()
GCMT = (-34.32, -72.24)


def read_yue(path=HERE / "yue2014_slip.txt"):
    rows = [l.split() for l in path.read_text().splitlines() if l[:1] in " 0123456789" and len(l.split()) == 17]
    return {c: np.array([float(r[i]) for r in rows]) for i, c in enumerate(COLS)}


def corners(m):
    """Esquinas (lon, lat) de cada subfalla: nodo al centro, xinc a lo largo del rumbo,
    yinc a lo largo del manteo (proyectado en horizontal con cos(dip))."""
    out = []
    for la, lo, st, dp, xi, yi in zip(m["lat"], m["lon"], m["strike"], m["dip"], m["xinc"], m["yinc"]):
        s = np.radians(st)
        u = np.array([np.sin(s), np.cos(s)]) * xi / 2
        d = np.array([np.cos(s), -np.sin(s)]) * yi / 2 * np.cos(np.radians(dp))
        kx = K * np.cos(np.radians(la))
        out.append([(lo + e / kx, la + n / K) for e, n in (-u - d, u - d, u + d, -u + d)])
    return out


def main():
    m = read_yue()
    polys = corners(m)
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, geom = om.subfaults(cfg)
    fp, sp = cfg.fault_plane, cfg.source_position
    la = np.array([sf.lat for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    lo = np.array([sf.lon for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    z = np.abs(np.array([sf.z_m for sf in geom.subfaults])).reshape(fp.ny, fp.nx) / 1e3
    ss, ds, _ = om.kinematic_slip(cfg, fm)
    skin = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    ss, ds, _ = om.dynamic_slip(cfg, geom, "B")
    sdyn = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    gj = json.loads((NAV / "faults" / "chaf_navidad.geojson").read_text())
    inter = json.loads((NAV / "faults" / "pichilemu_intersection.json").read_text())["line_in_mesh"]

    def which(lat, lon):
        k = next((i for i, p in enumerate(polys) if MPath(p).contains_point((lon, lat))), None)
        return None if k is None else dict(sub=int(m["sub"][k]), slip_m=float(m["slip"][k]), depth_km=float(m["depth"][k]),
                                           lat=float(m["lat"][k]), lon=float(m["lon"][k]))
    w = skin / skin.sum()
    kin_cen = (float((w * la).sum()), float((w * lo).sum()))
    info = dict(hypocentre_ISC=which(sp.latitude, sp.longitude), centroid_kin=which(*kin_cen), centroid_GCMT=which(*GCMT),
                navidad_peak_slip_kin=which(la.flat[skin.argmax()], lo.flat[skin.argmax()]),
                maule_max_slip_m=float(m["slip"].max()))
    # slip de Maule ponderado por el momento de Navidad (cinemático y dinámico)
    for name, s in (("kin", skin), ("dynB", sdyn)):
        sl = [which(a, b) for a, b in zip(la.flat, lo.flat)]
        mw = s.flatten() / s.sum()
        info[f"maule_slip_under_navidad_{name}"] = float(sum(wi * (x["slip_m"] if x else 0.0) for wi, x in zip(mw, sl)))
    (HERE / "maule_navidad.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))

    fig = plt.figure(figsize=(13, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])
    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    ax.set_extent([-74.0, -71.2, -35.4, -33.55])
    ax.add_feature(cfeature.LAND, facecolor="#efece4", zorder=0); ax.add_feature(cfeature.OCEAN, facecolor="#e6eef4", zorder=0)
    # Slip de Maule interpolado (nodos al centro de cada subfalla) en una grilla regular
    glon, glat = np.meshgrid(np.linspace(-74.6, -70.6, 400), np.linspace(-36.0, -33.3, 300))
    gs_ = griddata(np.c_[m["lon"], m["lat"]], m["slip"], (glon, glat), method="cubic")
    gs_ = np.clip(gs_, 0, None)  # ponytail: el cúbico sobreoscila bajo 0; se recorta
    pc = ax.pcolormesh(glon, glat, gs_, cmap="YlOrRd", vmin=0, vmax=18, alpha=0.6, shading="auto",
                       transform=ccrs.PlateCarree(), zorder=1)
    cs = ax.contour(glon, glat, gs_, levels=np.arange(2, 18, 2), colors="#7f0000", linewidths=0.5, alpha=0.7,
                    transform=ccrs.PlateCarree(), zorder=2)
    ax.clabel(cs, fmt="%d m", fontsize=6.5)
    ax.plot(m["lon"], m["lat"], "+", color="0.35", ms=5, mew=0.7, transform=ccrs.PlateCarree(), zorder=2)
    ax.plot([], [], "+", color="0.35", label="Maule subfault nodes (Yue et al. 2014)")
    fig.colorbar(pc, ax=ax, shrink=0.6, pad=0.02, label="Maule 2010 slip, Yue et al. (2014), interpolated (m)")
    ax.coastlines("10m", lw=0.9, zorder=3)
    gl = ax.gridlines(draw_labels=True, lw=0.3, alpha=0.5); gl.top_labels = gl.right_labels = False
    ax.contour(lo, la, skin, levels=[0.2 * skin.max(), 0.6 * skin.max()], colors="#1f4fbf", linewidths=[1.0, 2.0])
    ax.contour(lo, la, sdyn, levels=[0.2 * sdyn.max(), 0.6 * sdyn.max()], colors="k", linewidths=[0.8, 1.6], linestyles="--")
    ax.plot([], [], color="#1f4fbf", lw=1.6, label="Navidad 2019 kinematic slip (20, 60 % of max)")
    ax.plot([], [], color="k", ls="--", lw=1.3, label="Navidad 2019 dynamic B slip (20, 60 % of max)")
    for f in gj["features"]:
        if f["properties"].get("F_name") == "Pichilemu":
            x, y = np.array(f["geometry"]["coordinates"]).T
            ax.plot(x, y, color="#7a0177", lw=2.4, zorder=5)
    ax.plot([p["lon"] for p in inter], [p["lat"] for p in inter], ":", color="#7a0177", lw=2, zorder=5)
    ax.plot([], [], color="#7a0177", lw=2.4, label="Pichilemu fault (CHAF) and its interface intersection (dotted)")
    ax.plot(sp.longitude, sp.latitude, "*", color="lime", mec="k", ms=16, zorder=6, label="Navidad hypocentre (ISC)")
    ax.plot(GCMT[1], GCMT[0], "D", color="purple", mec="w", ms=8, zorder=6, label="Navidad centroid (GCMT)")
    ax.legend(loc="lower left", fontsize=7.2, framealpha=0.9)
    ax.set_title("Navidad 2019 on the Maule 2010 slip model", fontsize=11)

    # ---------- slip vs profundidad: filas de Maule cerca de Navidad y Navidad
    b = fig.add_subplot(gs[1])
    near = lambda x0: np.isclose(m["x"], x0, atol=2)
    for x0, col, lab in ((202, "#fd8d3c", "Maule row S of Navidad (~34.5–34.9°S)"),
                         (242, "#bd0026", "Maule row through Navidad (~34.1–34.5°S)"),
                         (282, "#feb24c", "Maule row N of Navidad (~33.8–34.2°S)")):
        k = near(x0)
        o = np.argsort(m["depth"][k])
        b.plot(m["slip"][k][o], m["depth"][k][o], "o-", color=col, lw=2 if x0 == 242 else 1.2, label=lab)
    b.plot(skin.max(axis=1) * 5, z.mean(axis=1), color="#1f4fbf", lw=2, label="Navidad kinematic slip ×5 (max per depth)")
    b.plot(sdyn.max(axis=1) * 5, z.mean(axis=1), color="k", ls="--", lw=1.6, label="Navidad dynamic B slip ×5")
    b.axhline(sp.depth, color="lime", lw=1.5, label="Navidad hypocentre depth (ISC)")
    b.set_ylim(60, 0); b.set_xlim(0, 18)
    b.set_xlabel("slip (m)"); b.set_ylabel("depth (km)")
    b.grid(alpha=0.25); b.legend(fontsize=7.2, loc="lower right", framealpha=0.9)
    b.set_title("Slip vs depth (Maule: 40×40 km subfaults)", fontsize=10)
    fig.savefig(HERE / "maule_navidad.png", dpi=130)


if __name__ == "__main__":
    main()
