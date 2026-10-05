"""Mapa del slip dinámico (mejor modelo) y del cinemático + dt0 sobre la malla de AXITRA.

Uso:  DYN_CASE=loc_mid python map_dynamic.py
Re-corre solo fd3d_TSN (~22 s) en tsn_work_report/. El slip se pasa a las
subfallas igual que en la inversión (bin_slip_rate_to_subfaults: columna 0 de
fd3d -> fila superior de AXITRA), así que el mapa muestra dónde está el slip
que genera los sintéticos. Escribe na_output_<caso>/report/map.png.
"""

import os
import sys

os.environ.setdefault("DYN_WORK", "tsn_work_report")

import cartopy.crs as ccrs  # noqa: E402
import cartopy.feature as cfeature  # noqa: E402
from pathlib import Path  # noqa: E402

import joblib  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import forward_dynamic as fd  # noqa: E402
import kdellipspy as kde  # noqa: E402
import na_dynamic as nd  # noqa: E402
from kdellipspy.core.geometry import TSNFaultGridSpec  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import TSNRunConfig, run_tsn_forward  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else fd.HERE / f"na_output_{fd.CASE_TAG}"
_MS = sorted((fd.CASE_DIR / "output_multiseed").glob("seed*.joblib"))
KIN = _MS or [fd.HERE.parent / f"{fd.CASE}_dt0" / "output" / "inversion_result.joblib"]
AGENCIES = {"GCMT": (-34.3200, -72.2400), "NEIC Mww": (-34.2433, -72.3015), "GFZ": (-34.2080, -72.1250)}


def grid_of(g, values, nx, ny):
    lat = np.array([sf.lat for sf in g.subfaults]).reshape(ny, nx)
    lon = np.array([sf.lon for sf in g.subfaults]).reshape(ny, nx)
    return lon, lat, np.asarray(values).reshape(ny, nx)


def main():
    cfg = kde.ConfigParser(str(fd.CASE_DIR / "input.ctl"))
    fd.setup_work_dir(cfg)
    fp, sp = cfg.fault_plane, cfg.source_position
    nx, ny = fp.nx, fp.ny
    fm = kde.AxitraForwardModel.from_config(cfg)
    g = fm.build_geometry()

    if (OUT / "na_result.joblib").exists():
        free = np.asarray(joblib.load(OUT / "na_result.joblib").best_model.model, float)
    else:  # corrida en curso: checkpoint en vivo
        vals = {l.split()[0]: float(l.split()[-1]) for l in open(OUT / "best_model_live.txt") if not l.startswith("#") and l.strip()}
        free = np.array([vals[n.split()[0]] for n in nd.NAMES])
    full = nd.to_full_model(free)
    from kdellipspy.core.geometry import tsn_hypocentre_coarse
    rate = run_tsn_forward(full, nx, ny, fd.tsn_grid(cfg), fd.tsn_run_cfg(), hypo=tsn_hypocentre_coarse(fp))
    slip_fine = np.hypot(rate["sliprateX"].sum(0), rate["sliprateZ"].sum(0)) * fd.DT  # (nxtT, nztT)
    if nd.M0_TARGET:  # modo M0 impuesto: slip real = slip de referencia x k
        mu = float(np.mean([sf.mu_pa for sf in g.subfaults]))
        slip_fine *= nd.M0_TARGET / (mu * slip_fine.sum() * fd.DH**2)
    fx, fz = fd.NXTT // nx, fd.NZTT // ny
    slip_sub = slip_fine.reshape(nx, fx, ny, fz).mean(axis=(1, 3))[:, ::-1]  # fd3d cuenta el dip desde abajo -> idip 0 = fila superior AXITRA
    lon, lat, sdyn = grid_of(g, slip_sub.T.ravel(), nx, ny)

    kin = min((joblib.load(k) for k in KIN), key=lambda r: r.best_model.misfit).best_model.model
    gk = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(kin[:7]), keep_all_sources=True)
    _, _, skin = grid_of(gk, [sf.slip_m for sf in gk.subfaults], nx, ny)

    corners = [(0, 0), (0, nx - 1), (ny - 1, nx - 1), (ny - 1, 0), (0, 0)]
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    st = cfg.stations.stations
    ax.set_extent([min(lon.min(), min(s.longitude for s in st)) - 0.15, max(lon.max(), max(s.longitude for s in st)) + 0.15,
                   min(lat.min(), min(s.latitude for s in st)) - 0.15, max(lat.max(), max(s.latitude for s in st)) + 0.15])
    ax.add_feature(cfeature.LAND, facecolor="#f2efe9")
    ax.add_feature(cfeature.OCEAN, facecolor="#dfeaf2")
    ax.coastlines("10m", lw=0.8)
    gl = ax.gridlines(draw_labels=True, lw=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False

    ax.plot([lon[r, c] for r, c in corners], [lat[r, c] for r, c in corners], "k--", lw=0.8, label="malla")
    ax.plot(lon[0, :], lat[0, :], "k-", lw=2.5, label="borde superior (somero)")
    cf = ax.contourf(lon, lat, sdyn, levels=np.linspace(0.1, sdyn.max(), 10), cmap="hot_r", alpha=0.85)
    fig.colorbar(cf, ax=ax, shrink=0.6, label="slip dinámico (m)")
    cs = ax.contour(lon, lat, skin, levels=[0.05 * skin.max(), 0.5 * skin.max()], colors="royalblue", linewidths=[1.2, 2.0])
    ax.plot([], [], color="royalblue", lw=1.5, label="cinemático+dt0 (5% y 50% del máx.)")
    ax.plot(sp.longitude, sp.latitude, "*", color="lime", mec="k", ms=18, label=f"hipocentro ({fd.CASE})")
    for (name, (la, lo)), mk in zip(AGENCIES.items(), ("D", "s", "^")):
        ax.plot(lo, la, mk, color="purple", mec="w", ms=9, label=f"centroide {name}")
    for s in st:
        ax.plot(s.longitude, s.latitude, "^", color="k", ms=8)
        ax.text(s.longitude + 0.02, s.latitude + 0.02, s.name, fontsize=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Navidad 2019 — dinámico ({OUT.name}) vs cinemático+dt0", fontsize=9)
    fig.savefig(OUT / "report" / "map.png", dpi=130, bbox_inches="tight")
    print("ok", OUT / "report" / "map.png")


if __name__ == "__main__":
    main()
