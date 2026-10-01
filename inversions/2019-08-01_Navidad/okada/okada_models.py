"""Deformación estática en superficie (Okada 1992, DC3D vía okada_wrapper) de los
modelos de slip de la configuración elegida (ISC + 1D regional + malla 32x40 km):
  - cinemático: mejor de las 5 semillas (isc_tomo/b0.05-0.30_noM14L)
  - dinámico A (M0 libre) y B (M0 del CSN, Mw 6.6)
Medio elástico homogéneo de Poisson (alpha = 2/3), cada subfalla un rectángulo.

Escribe okada_results.json y okada_<modelo>.png en esta carpeta.
Uso: python okada_models.py
"""

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
CASE = "isc_tomo/b0.05-0.30_noM14L"
os.environ["DYN_CASE"] = CASE
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path.insert(0, str(NAV / "dynamic"))

import cartopy.crs as ccrs  # noqa: E402
import cartopy.feature as cfeature  # noqa: E402
import joblib  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402
from okada_wrapper import dc3dwrapper  # noqa: E402

import forward_dynamic as fd  # noqa: E402
import kdellipspy as kde  # noqa: E402
import na_dynamic as nd  # noqa: E402
from kdellipspy.core.geometry import tsn_hypocentre_coarse  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import run_tsn_forward  # noqa: E402

ALPHA = 2.0 / 3.0  # (lambda+mu)/(lambda+2mu), Poisson solid
DYN_OUT = {"A": NAV / "dynamic" / "na_output_isc_tomo__b0.05-0.30_noM14L",
           "B": NAV / "dynamic" / "na_output_isc_tomo__b0.05-0.30_noM14L_M0"}
M0_B = 1.0e19


def subfaults(cfg):
    fm = kde.AxitraForwardModel.from_config(cfg)
    return fm, fm.build_geometry()


def kinematic_slip(cfg, fm):
    runs = [joblib.load(p) for p in sorted((fd.CASE_DIR / "output_multiseed").glob("seed*.joblib"))]
    best = min(runs, key=lambda r: r.best_model.misfit)
    g = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(best.best_model.model[:7]), keep_all_sources=True)
    slip = np.array([sf.slip_m for sf in g.subfaults])
    rake = np.radians(cfg.source_position.rake)
    return slip * np.cos(rake), slip * np.sin(rake), float(best.best_model.misfit)


def dynamic_slip(cfg, geom, which):
    """(strike-slip, dip-slip) por subfalla en el orden de la malla de AXITRA."""
    r = joblib.load(DYN_OUT[which] / "na_result.joblib")
    full = nd.to_full_model_named(r.param_names, np.asarray(r.best_model.model))  # A: 9 params, B: 8 (sin Te)
    fp = cfg.fault_plane
    rate = run_tsn_forward(full, fp.nx, fp.ny, fd.tsn_grid(cfg), fd.tsn_run_cfg(), hypo=tsn_hypocentre_coarse(fp))
    fx, fz = fd.NXTT // fp.nx, fd.NZTT // fp.ny

    def to_sub(rate_c):  # fino (nxtT, nztT), dip desde abajo -> (idip, istk) desde arriba
        s = rate_c.sum(axis=0) * fd.DT
        return s.reshape(fp.nx, fx, fp.ny, fz).mean(axis=(1, 3))[:, ::-1].T.ravel()

    ss, ds = to_sub(rate["sliprateX"]), to_sub(rate["sliprateZ"])
    if which == "B":  # mismo escalamiento que el forward (autosemejanza)
        mu = float(np.mean([sf.mu_pa for sf in geom.subfaults]))
        m0 = mu * np.hypot(rate["sliprateX"].sum(0), rate["sliprateZ"].sum(0)).sum() * fd.DT * fd.DH**2
        ss, ds = ss * M0_B / m0, ds * M0_B / m0
    return ss, ds, float(r.best_model.misfit)


def okada(geom, cfg, ss, ds, lat, lon):
    """Desplazamiento (E, N, Up) en metros en los puntos (lat, lon)."""
    sp, fp = cfg.source_position, cfg.fault_plane
    lat0, lon0 = sp.latitude, sp.longitude
    k = 111190.0
    E = (np.asarray(lon) - lon0) * k * np.cos(np.radians(lat0))
    N = (np.asarray(lat) - lat0) * k
    st = np.radians(sp.strike)
    dl, dw = fp.lx / fp.nx, fp.ly / fp.ny
    u = np.zeros((3,) + E.shape)
    keep = np.hypot(ss, ds) > 0.01 * np.hypot(ss, ds).max()
    for sf, s1, s2 in zip(np.array(geom.subfaults)[keep], ss[keep], ds[keep]):
        e0 = (sf.lon - lon0) * k * np.cos(np.radians(lat0))
        n0 = (sf.lat - lat0) * k
        # Marco DC3D: x a lo largo del strike, y 90° a la izquierda (lado somero).
        x = (E - e0) * np.sin(st) + (N - n0) * np.cos(st)
        y = -(E - e0) * np.cos(st) + (N - n0) * np.sin(st)
        for idx in np.ndindex(E.shape):
            _, uu, _ = dc3dwrapper(ALPHA, [x[idx], y[idx], 0.0], sf.z_m, sp.dip,
                                   [-dl / 2, dl / 2], [-dw / 2, dw / 2], [s1, s2, 0.0])
            u[0][idx] += uu[0] * np.sin(st) - uu[1] * np.cos(st)
            u[1][idx] += uu[0] * np.cos(st) + uu[1] * np.sin(st)
            u[2][idx] += uu[2]
    return u


def plot(name, title, lon, lat, u, cfg, geom, path):
    """Dos paneles: uh (magnitud horizontal + flechas de dirección) y uz (vertical)."""
    fp = cfg.fault_plane
    la = np.array([sf.lat for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    lo = np.array([sf.lon for sf in geom.subfaults]).reshape(fp.ny, fp.nx)
    uh = np.hypot(u[0], u[1]) * 1e3
    uz = u[2] * 1e3
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), subplot_kw={"projection": ccrs.PlateCarree()})
    vz = np.abs(uz).max()
    panels = [
        (axes[0], uh, "viridis", 0.0, uh.max(), "desplazamiento horizontal uh (mm)"),
        (axes[1], uz, "RdBu_r", -vz, vz, "desplazamiento vertical uz (mm)"),
    ]
    c = [(0, 0), (0, fp.nx - 1), (fp.ny - 1, fp.nx - 1), (fp.ny - 1, 0), (0, 0)]
    for ax, field, cmap, vmin, vmax, label in panels:
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])
        cf = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        fig.colorbar(cf, ax=ax, shrink=0.72, label=label, pad=0.02)
        ax.contour(lon, lat, field, levels=8, colors="k", linewidths=0.35, alpha=0.5)
        ax.coastlines("10m", lw=0.9)
        gl = ax.gridlines(draw_labels=True, lw=0.3, alpha=0.5)
        gl.top_labels = gl.right_labels = False
        gl.xlocator = mticker.FixedLocator([-73.0, -72.5, -72.0, -71.5, -71.0])
        ax.plot([lo[r, s_] for r, s_ in c], [la[r, s_] for r, s_ in c], "w--" if cmap == "viridis" else "k--", lw=0.8)
        ax.plot(lo[0, :], la[0, :], "w-" if cmap == "viridis" else "k-", lw=2.2)
        ax.plot(cfg.source_position.longitude, cfg.source_position.latitude, "*", color="lime", mec="k", ms=15)
        for st_ in cfg.stations.stations:
            ax.plot(st_.longitude, st_.latitude, "^", color="w" if cmap == "viridis" else "k", mec="k", ms=7)
            ax.text(st_.longitude + 0.02, st_.latitude + 0.02, st_.name, fontsize=8,
                    color="w" if cmap == "viridis" else "k")
    step = 4
    q = axes[0].quiver(lon[::step, ::step], lat[::step, ::step], u[0][::step, ::step] * 1e3, u[1][::step, ::step] * 1e3,
                       scale=uh.max() * 12, width=0.003, color="w")
    hmax = round(uh.max(), -1) or 10
    axes[0].quiverkey(q, 0.78, 0.04, hmax, f"{hmax:.0f} mm", labelpos="E", coordinates="axes", color="w",
                      labelcolor="w")
    axes[0].set_title(f"{title}\nuh · horizontal", fontsize=10)
    axes[1].set_title(f"{title}\nuz · vertical (+ alzamiento)", fontsize=10)
    fig.savefig(path, dpi=115, bbox_inches="tight")
    plt.close(fig)


def main():
    cfg = kde.ConfigParser(str(fd.CASE_DIR / "input.ctl"))
    fd.setup_work_dir(cfg)
    fm, geom = subfaults(cfg)
    lon1d = np.linspace(-73.1, -71.1, 61)
    lat1d = np.linspace(-35.0, -33.5, 51)
    lon, lat = np.meshgrid(lon1d, lat1d)
    sta = cfg.stations.stations
    models = {}
    for key, title, slipfun in [
        ("kin", "Cinemático (mejor de 5 semillas)", lambda: kinematic_slip(cfg, fm)),
        ("A", "Dinámico A (M0 libre)", lambda: dynamic_slip(cfg, geom, "A")),
        ("B", "Dinámico B (M0 CSN, Mw 6.6)", lambda: dynamic_slip(cfg, geom, "B")),
    ]:
        ss, ds, misfit = slipfun()
        u = okada(geom, cfg, ss, ds, lat, lon)
        us = okada(geom, cfg, ss, ds, np.array([s.latitude for s in sta]), np.array([s.longitude for s in sta]))
        mu = np.array([sf.mu_pa for sf in geom.subfaults])
        m0 = float((mu * np.hypot(ss, ds) * cfg.fault_plane.lx / cfg.fault_plane.nx * cfg.fault_plane.ly / cfg.fault_plane.ny).sum())
        plot(key, f"{title} · Okada", lon, lat, u, cfg, geom, HERE / f"okada_{key}.png")
        models[key] = {
            "title": title, "misfit": misfit, "mw": (np.log10(m0) - 9.1) / 1.5,
            "slip_max": float(np.hypot(ss, ds).max()),
            "uplift_max_mm": float(u[2].max() * 1e3), "subsidence_max_mm": float(u[2].min() * 1e3),
            "horiz_max_mm": float(np.hypot(u[0], u[1]).max() * 1e3),
            "stations": {s.name: [round(float(v) * 1e3, 1) for v in us[:, i]] for i, s in enumerate(sta)},
        }
        print(f"{key}: Mw {models[key]['mw']:.2f} | uz {models[key]['uplift_max_mm']:+.0f} / {models[key]['subsidence_max_mm']:+.0f} mm"
              f" | horiz max {models[key]['horiz_max_mm']:.0f} mm", flush=True)
    (HERE / "okada_results.json").write_text(json.dumps(models, ensure_ascii=False, indent=1))


def _check():
    """Falla inversa de prueba: el bloque colgante sube y se mueve hacia la fosa (lado somero, +y)."""
    _, u, _ = dc3dwrapper(ALPHA, [0.0, 5e3, 0.0], 10e3, 20.0, [-5e3, 5e3], [-5e3, 5e3], [0.0, 1.0, 0.0])
    assert u[2] > 0 and u[1] > 0
    _, u, _ = dc3dwrapper(ALPHA, [0.0, -20e3, 0.0], 10e3, 20.0, [-5e3, 5e3], [-5e3, 5e3], [0.0, 1.0, 0.0])
    assert u[2] < 0


if __name__ == "__main__":
    _check()
    main()
