"""Mapa: fallas CHAF, falla de Pichilemu, su intersección con el plano de Navidad
y el slip de los modelos (cinemático y dinámico B). Escribe pichilemu_map.png."""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
os.environ.setdefault("DYN_CASE", "isc_tomo/b0.05-0.30_noM14L")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path.insert(0, str(NAV / "okada"))
sys.path.insert(0, str(NAV / "dynamic"))
import cartopy.crs as ccrs  # noqa: E402
import cartopy.feature as cfeature  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import okada_models as om  # noqa: E402


def main():
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, geom = om.subfaults(cfg)
    fp, sp = cfg.fault_plane, cfg.source_position
    la = np.array([sf.lat for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    lo = np.array([sf.lon for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    ss, ds, _ = om.kinematic_slip(cfg, fm)
    skin = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    ss, ds, _ = om.dynamic_slip(cfg, geom, "B")
    sdyn = np.hypot(ss, ds).reshape(fp.ny, fp.nx)
    gj = json.loads((HERE / "chaf_navidad.geojson").read_text())
    inter = json.loads((HERE / "pichilemu_intersection.json").read_text())["line_in_mesh"]
    split = json.loads((HERE / "pichilemu_moment_split.json").read_text())

    fig = plt.figure(figsize=(8.5, 8.5))
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    ax.set_extent([-72.85, -71.55, -34.75, -33.85])
    ax.add_feature(cfeature.LAND, facecolor="#efece4")
    ax.add_feature(cfeature.OCEAN, facecolor="#dfeaf2")
    ax.coastlines("10m", lw=0.9)
    gl = ax.gridlines(draw_labels=True, lw=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    cf = ax.contourf(lo, la, sdyn, levels=np.linspace(0.1, sdyn.max(), 9), cmap="Oranges", alpha=0.85)
    fig.colorbar(cf, ax=ax, shrink=0.6, label="slip dinámico B (m)")
    ax.contour(lo, la, skin, levels=[0.05 * skin.max(), 0.5 * skin.max()], colors="royalblue", linewidths=[1.2, 2.2])
    ax.plot([], [], color="royalblue", lw=1.8, label="slip cinemático (5 % y 50 % del máx.)")
    for f in gj["features"]:
        x, y = np.array(f["geometry"]["coordinates"]).T
        pich = f["properties"].get("F_name") == "Pichilemu"
        ax.plot(x, y, color="crimson" if pich else "0.35", lw=2.6 if pich else 0.9, zorder=4)
    ax.plot([], [], color="crimson", lw=2.6, label="falla de Pichilemu (CHAF, 138°/55° SW, normal)")
    ax.plot([], [], color="0.35", lw=0.9, label="otras fallas CHAF")
    ax.plot([p["lon"] for p in inter], [p["lat"] for p in inter], "--", color="crimson", lw=2.2, zorder=5,
            label=f"intersección Pichilemu–Navidad ({inter[0]['depth_km']:.0f}–{inter[-1]['depth_km']:.0f} km)")
    c = [(0, 0), (0, fp.nx - 1), (fp.ny - 1, fp.nx - 1), (fp.ny - 1, 0), (0, 0)]
    ax.plot([lo[r, s] for r, s in c], [la[r, s] for r, s in c], "k:", lw=0.8)
    ax.plot(sp.longitude, sp.latitude, "*", color="lime", mec="k", ms=17, zorder=6, label="hipocentro ISC")
    for (name, (y, x)), mk in zip({"GCMT": (-34.32, -72.24), "NEIC Mww": (-34.2433, -72.3015)}.items(), ("D", "s")):
        ax.plot(x, y, mk, color="purple", mec="w", ms=9, zorder=6, label=f"centroide {name}")
    ax.text(-72.83, -34.73, f"Momento que cruza la intersección: cinemático {100 * split['kin']['frac_other_side']:.0f} %, "
            f"dinámico B {100 * split['B']['frac_other_side']:.0f} %", fontsize=8.5,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    ax.set_title("Navidad 2019: slip y falla de Pichilemu", fontsize=11)
    fig.savefig(HERE / "pichilemu_map.png", dpi=125, bbox_inches="tight")
    print("ok")


if __name__ == "__main__":
    main()
