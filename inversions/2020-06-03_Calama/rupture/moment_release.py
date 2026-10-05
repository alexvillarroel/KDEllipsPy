"""Liberación de momento en el espacio y el tiempo: dinámico B (M0 CSN) y mejor
cinemático (5 semillas), configuración ISC + 1D regional + 0.05-0.30 Hz sin M14L.

Mdot(x, z, t) = mu * v(x, z, t) * dA
  - dinámico: v = slip-rate de fd3d_TSN (celdas de 500 m, dt 0.015 s), escalado por
    k = M0_CSN / M0_simulado (autosemejanza), con el mismo mu medio del forward;
  - cinemático: cada subfalla libera M0_i con un triángulo de ancho T0 que empieza
    en su tiempo de ruptura t_i = distancia / Vr.

Escribe moment_release.png (STF, espacio-tiempo, tiempos de ruptura) y
snapshots_B.png (slip-rate del dinámico en el plano de falla).
"""

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
os.environ.setdefault("DYN_CASE", "dyn_cases/np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path[:0] = [str(NAV / "okada"), str(NAV / "dynamic")]

import joblib  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import okada_models as om  # noqa: E402
from kdellipspy.core.geometry import tsn_hypocentre_coarse  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import run_tsn_forward  # noqa: E402

M0_B = 1.0e19
DT_PLOT = 0.1  # s
T_MAX = 22.0
CMAP = "Oranges"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})


def dynamic_B(cfg, geom):
    fp = cfg.fault_plane
    r = joblib.load(om.DYN_OUT["B"] / "na_result.joblib")
    full = om.nd.to_full_model_named(r.param_names, np.asarray(r.best_model.model))
    rate = run_tsn_forward(full, fp.nx, fp.ny, om.fd.tsn_grid(cfg), om.fd.tsn_run_cfg(), hypo=tsn_hypocentre_coarse(fp))
    v = np.hypot(rate["sliprateX"], rate["sliprateZ"])[:, :, ::-1]  # (t, strike, dip desde el borde superior)
    dh, dt = om.fd.DH / 1e3, om.fd.DT
    mu = float(np.mean([sf.mu_pa for sf in geom.subfaults]))
    k = M0_B / (mu * v.sum() * dt * (dh * 1e3) ** 2)
    mdot = mu * v * k * (dh * 1e3) ** 2  # N·m/s por celda
    step = int(round(DT_PLOT / dt))
    nt = (mdot.shape[0] // step) * step
    mdot = mdot[:nt].reshape(nt // step, step, *mdot.shape[1:]).mean(axis=1)
    x = (np.arange(mdot.shape[1]) + 0.5) * dh
    w = (np.arange(mdot.shape[2]) + 0.5) * dh
    return mdot, x, w, np.arange(mdot.shape[0]) * DT_PLOT, k


def kinematic(cfg, fm):
    fp = cfg.fault_plane
    runs = [joblib.load(p) for p in sorted((om.fd.CASE_DIR / "output_multiseed").glob("seed*.joblib"))]
    best = min(runs, key=lambda r: r.best_model.misfit)
    g = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(best.best_model.model[:7]), keep_all_sources=True)
    m0 = np.array([sf.mu_pa * sf.area_m2 * sf.slip_m for sf in g.subfaults]).reshape(fp.ny, fp.nx).T  # (strike, dip)
    tr = np.array([sf.rupture_time_s for sf in g.subfaults]).reshape(fp.ny, fp.nx).T
    T0 = float(cfg.ellipse.t0)
    t = np.arange(0, T_MAX, DT_PLOT)
    tri = lambda s: np.clip(1 - np.abs(s - T0 / 2) / (T0 / 2), 0, None) * (2 / T0)  # área 1, ancho T0
    mdot = m0[None] * tri(t[:, None, None] - tr[None])
    dl, dw = fp.lx / fp.nx / 1e3, fp.ly / fp.ny / 1e3
    vr = float(best.best_model.model[6])
    return mdot, (np.arange(fp.nx) + 0.5) * dl, (np.arange(fp.ny) + 0.5) * dw, t, vr


def rupture_time(mdot, t):
    """Instante en que cada celda alcanza el 50 % de su momento (NaN si casi no desliza)."""
    cum = np.cumsum(mdot, axis=0)
    tot = cum[-1]
    t50 = t[np.argmax(cum >= 0.5 * tot[None], axis=0)]
    return np.where(tot > 0.05 * tot.max(), t50, np.nan)


def main():
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, geom = om.subfaults(cfg)
    fp = cfg.fault_plane
    hx, hy = fp.hx / 1e3, fp.hy / 1e3
    inter = json.loads((NAV / "faults" / "pichilemu_intersection.json").read_text())["line_in_mesh"]
    ix, iw = [p["strike_km"] for p in inter], [p["dip_km"] for p in inter]

    dB, xB, wB, tB, k = dynamic_B(cfg, geom)
    dK, xK, wK, tK, vrK = kinematic(cfg, fm)
    models = [("Dynamic B (CSN M0)", dB, xB, wB, tB), ("Kinematic (best of 5 seeds)", dK, xK, wK, tK)]

    fig = plt.figure(figsize=(12.5, 13.5), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1, 1, 1.25])
    ax = fig.add_subplot(gs[0, :])
    for (name, d, x, w, t), col in zip(models, ("black", "#4a3aa7")):
        stf = d.sum(axis=(1, 2))
        m0 = stf.sum() * DT_PLOT
        ax.plot(t, stf / 1e18, color=col, lw=2, label=f"{name}: M0 {m0:.2e} N·m, Mw {(np.log10(m0) - 9.1) / 1.5:.2f}")
        ax.fill_between(t, stf / 1e18, color=col, alpha=0.08)
    ax.set_xlim(0, T_MAX)
    ax.set_xlabel("Time since rupture onset (s)")
    ax.set_ylabel("Moment rate (10¹⁸ N·m/s)")
    ax.set_title("(a) Source time function")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.25)

    for col, (name, d, x, w, t) in enumerate(models):
        # (b) a lo largo del rumbo
        a = fig.add_subplot(gs[1, col])
        st = d.sum(axis=2) / 1e18
        im = a.pcolormesh(x - hx, t, st, cmap=CMAP, shading="auto")
        fig.colorbar(im, ax=a, label="10¹⁸ N·m/s")
        a.axvline(0, color="k", lw=0.6, ls=":")
        for vref, ls in ((1.7, "--"), (3.9, ":")):
            tt = np.linspace(0, T_MAX, 50)
            a.plot(vref * tt, tt, color="0.25", lw=0.8, ls=ls)
            a.plot(-vref * tt, tt, color="0.25", lw=0.8, ls=ls)
        a.set_xlim((x - hx).min(), (x - hx).max()); a.set_ylim(T_MAX, 0)
        a.set_xlabel("Along strike from hypocentre (km, + NNE)"); a.set_ylabel("Time (s)")
        a.set_title(f"({'bc'[col]}) {name}: strike–time")
        a.text(0.02, 0.04, "-- 1.7 km/s   ··· 3.9 km/s (≈Vs)", transform=a.transAxes, fontsize=7.5, color="0.25")
        # (d) a lo largo del manteo
        a = fig.add_subplot(gs[2, col])
        dp = d.sum(axis=1) / 1e18
        im = a.pcolormesh(w - hy, t, dp, cmap=CMAP, shading="auto")
        fig.colorbar(im, ax=a, label="10¹⁸ N·m/s")
        a.axvline(0, color="k", lw=0.6, ls=":")
        for vref, ls in ((1.7, "--"), (3.9, ":")):
            tt = np.linspace(0, T_MAX, 50)
            a.plot(vref * tt, tt, color="0.25", lw=0.8, ls=ls)
            a.plot(-vref * tt, tt, color="0.25", lw=0.8, ls=ls)
        a.set_xlim((w - hy).min(), (w - hy).max()); a.set_ylim(T_MAX, 0)
        a.set_xlabel("Along dip from hypocentre (km, + downdip)"); a.set_ylabel("Time (s)")
        a.set_title(f"({'de'[col]}) {name}: dip–time")
        # (f) tiempo de ruptura en el plano
        a = fig.add_subplot(gs[3, col])
        trp = rupture_time(d, t)
        slip_m0 = d.sum(axis=0) * DT_PLOT
        im = a.pcolormesh(x, w, (slip_m0 / slip_m0.max()).T, cmap=CMAP, shading="auto", vmin=0, vmax=1)
        fig.colorbar(im, ax=a, label="released moment (normalized)")
        cs = a.contour(x, w, trp.T, levels=np.arange(1, 20, 2), colors="k", linewidths=0.7)
        a.clabel(cs, fmt="%d s", fontsize=7)
        a.plot(ix, iw, color="#c4161c", lw=2, ls="--", label="Pichilemu fault intersection")
        a.plot(hx, hy, "*", color="lime", mec="k", ms=14, label="ISC hypocentre")
        a.set_xlim(0, fp.lx / 1e3); a.set_ylim(fp.ly / 1e3, 0); a.set_aspect("equal")
        a.set_xlabel("Along strike (km)"); a.set_ylabel("Along dip from top edge (km)")
        a.set_title(f"({'fg'[col]}) {name.split(' (')[0]}: moment and isochrones")
        a.legend(loc="lower right", fontsize=7.5, framealpha=0.85)
    fig.suptitle(f"Calama 2020: moment release (dynamic scaled by k={k:.2f}; kinematic Vr {vrK:.2f} km/s)", fontsize=11)
    fig.savefig(HERE / "moment_release.png", dpi=130)
    plt.close(fig)

    # Instantáneas del slip-rate del dinámico
    times = [1, 3, 5, 7, 9, 11, 13, 15]
    # Escala desde t = 2 s: los primeros instantes los domina el parche de nucleación
    # forzada (sobreesfuerzo artificial) y apagarían el resto del frente.
    vmax = dB[int(2 / DT_PLOT):].max() / 1e15
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 7.2), constrained_layout=True)
    for a, tt in zip(axes.ravel(), times):
        i = int(round(tt / DT_PLOT))
        im = a.pcolormesh(xB, wB, np.minimum(dB[i].T / 1e15, vmax), cmap=CMAP, shading="auto", vmin=0, vmax=vmax)
        a.plot(ix, iw, color="#c4161c", lw=1.6, ls="--")
        a.plot(hx, hy, "*", color="lime", mec="k", ms=11)
        a.set_xlim(0, fp.lx / 1e3); a.set_ylim(fp.ly / 1e3, 0); a.set_aspect("equal")
        a.set_title(f"t = {tt} s", fontsize=9)
        a.tick_params(labelsize=7)
    fig.colorbar(im, ax=axes, shrink=0.7, label="moment rate per cell (10¹⁵ N·m/s)")
    fig.supxlabel("Along strike (km)"); fig.supylabel("Along dip from top edge (km)")
    fig.suptitle("Dynamic B: moment rate on the fault plane (red: Pichilemu fault intersection; "
                 "colour scale saturated for t < 2 s by forced nucleation)", fontsize=10.5)
    fig.savefig(HERE / "snapshots_B.png", dpi=120)
    print("ok")


if __name__ == "__main__":
    main()
