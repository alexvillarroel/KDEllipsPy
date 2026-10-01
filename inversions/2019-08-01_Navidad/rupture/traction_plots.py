"""Tracción y slip-rate del dinámico B en el plano de falla: versión en Python de
PrintSnapshot.m / PrintTimeSeries.m de fd3d_TSN (examples/*), con el hipocentro y
la intersección con Pichilemu.

  - result/shearstressZ.res = tracción ABSOLUTA tz + T0Z (fd3d_theo.f90:424), Pa;
  - result/sliprateZ.res    = slip-rate en el manteo (m/s);
  - float32: nuestro build de gfortran no promueve a doble (los .m leen real*8 porque
    los ejemplos compilan con -r8 / -autodouble).
Modelo B en modo M0 impuesto: esfuerzos y slip-rate se escalan por k (autosemejanza).

Escribe traction_snapshots_B.png, traction_timeseries_B.png y stress_change_final_B.png.
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

import joblib  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import okada_models as om  # noqa: E402
from kdellipspy.core.geometry import tsn_hypocentre_coarse  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import read_tsn_fault_field, run_tsn_forward  # noqa: E402

M0_B = 1.0e19
PICH = "#c4161c"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5})


def load_B():
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, geom = om.subfaults(cfg)
    fp = cfg.fault_plane
    r = joblib.load(om.DYN_OUT["B"] / "na_result.joblib")
    full = om.nd.to_full_model_named(r.param_names, np.asarray(r.best_model.model))
    rate = run_tsn_forward(full, fp.nx, fp.ny, om.fd.tsn_grid(cfg), om.fd.tsn_run_cfg(), hypo=tsn_hypocentre_coarse(fp))
    nx, nz, dt, dh = om.fd.NXTT, om.fd.NZTT, om.fd.DT, om.fd.DH
    tau = read_tsn_fault_field(om.fd.WORK / "result" / "shearstressZ.res", nx, nz)
    v = np.hypot(rate["sliprateX"], rate["sliprateZ"])
    mu = float(np.mean([sf.mu_pa for sf in geom.subfaults]))
    k = M0_B / (mu * v.sum() * dt * dh**2)
    fr = np.loadtxt(om.fd.WORK / "result" / "friction.inp")  # T0X, T0Z, peak, Dc, ... (dip desde abajo, i rápido)
    t0 = fr[:, 1].reshape(nz, nx).T
    peak = fr[:, 2].reshape(nz, nx).T
    flip = lambda a: a[..., ::-1]  # dip desde el borde superior, como AXITRA
    return dict(cfg=cfg, k=k, dt=dt, dh=dh / 1e3, tau=flip(tau) * k, v=flip(v) * k,
                t0=flip(t0) * k, peak=flip(peak) * k, x=(np.arange(nx) + 0.5) * dh / 1e3, w=(np.arange(nz) + 0.5) * dh / 1e3)


def decorate(a, fp, ix, iw):
    a.plot(ix, iw, color=PICH, lw=1.6, ls="--")
    a.plot(fp.hx / 1e3, fp.hy / 1e3, "*", color="lime", mec="k", ms=11)
    a.set_xlim(0, fp.lx / 1e3); a.set_ylim(fp.ly / 1e3, 0); a.set_aspect("equal")
    a.tick_params(labelsize=7)


def main():
    B = load_B()
    fp = B["cfg"].fault_plane
    inter = json.loads((NAV / "faults" / "pichilemu_intersection.json").read_text())["line_in_mesh"]
    ix, iw = [p["strike_km"] for p in inter], [p["dip_km"] for p in inter]
    x, w, dt = B["x"], B["w"], B["dt"]
    dtau = (B["tau"] - B["t0"][None]) / 1e6  # MPa
    asp = B["t0"] > 0  # dentro de la aspereza (fuera: barrera, T0 = 0)

    # ---------------- instantáneas: Δτ y slip-rate
    times = [1, 3, 5, 7, 9, 11]
    lim = np.nanpercentile(np.abs(dtau[int(2 / dt):]), 99.5)
    vmax = B["v"][int(2 / dt):].max()
    fig, axes = plt.subplots(2, len(times), figsize=(14, 6.4), constrained_layout=True)
    for j, tt in enumerate(times):
        i = int(round(tt / dt))
        im1 = axes[0, j].pcolormesh(x, w, dtau[i].T, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto")
        im2 = axes[1, j].pcolormesh(x, w, np.minimum(B["v"][i].T, vmax), cmap="Oranges", vmin=0, vmax=vmax, shading="auto")
        for a in axes[:, j]:
            decorate(a, fp, ix, iw)
        axes[0, j].set_title(f"t = {tt} s")
    fig.colorbar(im1, ax=axes[0], shrink=0.85, label="Δτ = τ(t) − τ₀ (MPa)")
    fig.colorbar(im2, ax=axes[1], shrink=0.85, label="slip rate (m/s)")
    fig.supxlabel("Along strike (km)"); fig.supylabel("Along dip from top edge (km)")
    fig.suptitle("Dynamic B: shear stress change (top; red = loading, blue = drop) and slip rate (bottom). "
                 "Dashed red: Pichilemu fault intersection", fontsize=10)
    fig.savefig(HERE / "traction_snapshots_B.png", dpi=120)
    plt.close(fig)

    # ---------------- series de tiempo en 3 puntos
    hx, hy = fp.hx / 1e3, fp.hy / 1e3
    pts = {"hypocentre": (hx, hy), "NE of hypocentre (+6 km strike, +6 km dip)": (hx + 6, hy + 6),
           "SW, across Pichilemu (−3.5 km strike, +1 km dip)": (hx - 3.5, hy + 1)}
    t = np.arange(B["tau"].shape[0]) * dt
    fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True, constrained_layout=True)
    for j, (name, (px, pw)) in enumerate(pts.items()):
        i, kk = int(np.argmin(np.abs(x - px))), int(np.argmin(np.abs(w - pw)))
        a = axes[0, j]
        a.plot(t, B["tau"][:, i, kk] / 1e6, color="k", lw=1.4, label="τ(t)")
        a.axhline(B["peak"][i, kk] / 1e6, color=PICH, ls="--", lw=1, label=f"peak strength {B['peak'][i, kk] / 1e6:.1f} MPa")
        a.axhline(B["t0"][i, kk] / 1e6, color="0.5", ls=":", lw=1, label=f"τ₀ {B['t0'][i, kk] / 1e6:.1f} MPa")
        a.set_title(name); a.set_ylabel("traction (MPa)"); a.legend(fontsize=7, frameon=False)
        a.set_ylim(bottom=min(-0.5, a.get_ylim()[0]))
        b = axes[1, j]
        b.plot(t, B["v"][:, i, kk], color="#d94801", lw=1.4)
        b.set_ylabel("slip rate (m/s)"); b.set_xlabel("time (s)")
        slip = B["v"][:, i, kk].sum() * dt
        b.text(0.97, 0.9, f"final slip {slip:.2f} m", transform=b.transAxes, ha="right", fontsize=8)
    for a in axes.ravel():
        a.set_xlim(0, 20); a.grid(alpha=0.25)
    fig.suptitle("Dynamic B: traction and slip rate at three fault points", fontsize=11)
    fig.savefig(HERE / "traction_timeseries_B.png", dpi=120)
    plt.close(fig)

    # ---------------- cambio de esfuerzo final
    fin = dtau[-1]
    fig, a = plt.subplots(figsize=(6.4, 6.6), constrained_layout=True)
    lim = np.nanpercentile(np.abs(fin), 99.5)
    im = a.pcolormesh(x, w, fin.T, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto")
    a.contour(x, w, asp.T.astype(float), levels=[0.5], colors="k", linewidths=0.8)
    decorate(a, fp, ix, iw)
    fig.colorbar(im, ax=a, shrink=0.8, label="final Δτ (MPa): blue = drop, red = loading")
    a.set_xlabel("Along strike (km)"); a.set_ylabel("Along dip from top edge (km)")
    drop = -fin[asp & (B["v"].sum(0) * dt > 0.05 * (B["v"].sum(0) * dt).max())].mean()
    a.set_title(f"Dynamic B: static stress change (mean drop where slipping {drop:.1f} MPa)\n"
                "red ring: concentration at the imposed barrier edge (black); amplitude is mesh-dependent", fontsize=9)
    fig.savefig(HERE / "stress_change_final_B.png", dpi=130)
    print(f"ok: k={B['k']:.2f}, caída media {drop:.2f} MPa")


if __name__ == "__main__":
    main()
