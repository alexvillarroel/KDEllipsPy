"""Instantaneas de la ruptura, slip dinamico contra cinematico y perfiles.

Tres figuras, todas sobre el plano de falla:

  snapshots_<tag>.png  slip-rate en la malla fina a varios tiempos, con el
                       frente de ruptura (isocrona) encima. Muestra como crece
                       la ruptura, no solo el resultado final.
  slip_compare_<tag>.png  slip final del dinamico y del cinematico en la MISMA
                       malla de subfallas y la MISMA escala de color, con las
                       dos elipses dibujadas encima: la del dinamico (que es la
                       aspereza de prestress, o sea una condicion de borde) y
                       la del cinematico (que es la forma del slip impuesta).
                       Son objetos distintos y por eso se dibujan distinto.
  profiles_<tag>.png   cortes en rumbo y en manteo por el hipocentro.

Usa okada_models.kinematic_slip / dynamic_slip, que devuelven las dos cosas en
el orden de la malla de AXITRA (idip desde arriba, istk mas rapido), asi que la
comparacion es celda a celda.

Uso:  DYN_CASE=<caso dinamico> python slip_comparison.py [A|B]
"""

import os
import sys
from math import cos, pi, radians, sin
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

os.environ.setdefault("DYN_WORK", "tsn_work_cmp")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "okada"))

import kdellipspy as kde  # noqa: E402
import forward_dynamic as fd  # noqa: E402
import na_dynamic as nd  # noqa: E402
import okada_models as om  # noqa: E402
from kdellipspy.inversion.dynamic.tsn_bridge import read_tsn_fault_field  # noqa: E402

WHICH = (sys.argv[1] if len(sys.argv) > 1 else "A").upper()
TAG = f"{fd.CASE.split('/')[-1].split('_')[0]}_{WHICH}"
HERE = Path(__file__).resolve().parent


def dyn_ellipse(cfg):
    """Elipse de la aspereza del dinamico, en km sobre el plano (rumbo, manteo
    desde el borde SUPERIOR). fd3d cuenta el manteo desde el borde profundo, asi
    que la coordenada se invierte para dejarla en el marco de AXITRA."""
    r = joblib.load(om.DYN_OUT[WHICH] / "na_result.joblib")
    full = nd.to_full_model_named(r.param_names, np.asarray(r.best_model.model))
    fp = cfg.fault_plane
    dstk, ddip = fp.lx / fp.nx / 1e3, fp.ly / fp.ny / 1e3   # km por punto grueso
    a, b, xo, yo, phi = (float(v) for v in full[:5])
    xc = (xo - 0.5) * dstk
    yc = fp.ly / 1e3 - (yo - 0.5) * ddip                     # manteo desde arriba
    return xc, yc, a * dstk, b * ddip, phi, float(r.best_model.misfit)


def kin_ellipse(cfg):
    """Elipse del cinematico, misma convencion. Inversa de
    EllipticalSlipMapper._prepare (core/geometry.py)."""
    runs = [joblib.load(p) for p in sorted((om.KIN_CASE_DIR / "multiseed").glob("seed*.joblib"))]
    best = min(runs, key=lambda r: r.best_model.misfit)
    a1, a2, th, npf, tp = (float(v) for v in np.asarray(best.best_model.model)[:5])
    alpha, t = th * pi, tp * 2 * pi
    x01, y01 = a1 * npf * cos(t), a2 * npf * sin(t)
    fp = cfg.fault_plane
    xc = (x01 * cos(alpha) + y01 * sin(alpha) + fp.hx / 1e3)
    yc = (-x01 * sin(alpha) + y01 * cos(alpha) + fp.hy / 1e3)
    return xc, yc, a1, a2, alpha, float(best.best_model.misfit)


def draw_ellipse(ax, xc, yc, a, b, ang, color, ls, label):
    ax.add_patch(Ellipse((xc, yc), 2 * a, 2 * b, angle=-np.degrees(ang),
                         fill=False, ec=color, lw=1.8, ls=ls, label=label))


def main():
    cfg = kde.ConfigParser(str(fd.CASE_DIR / "input.ctl"))
    fp = cfg.fault_plane
    Lx, Ly = fp.lx / 1e3, fp.ly / 1e3
    hx, hy = fp.hx / 1e3, fp.hy / 1e3
    fd.setup_work_dir(cfg)   # crea tsn_work_cmp/ y compila fd3d si falta
    fm, geom = om.subfaults(cfg)

    ss_d, ds_d, mis_d = om.dynamic_slip(cfg, geom, WHICH)
    ss_k, ds_k, mis_k = om.kinematic_slip(cfg, fm)
    slip_d = np.hypot(ss_d, ds_d).reshape(fp.ny, fp.nx)
    slip_k = np.hypot(ss_k, ds_k).reshape(fp.ny, fp.nx)
    ext = [0, Lx, Ly, 0]            # origin arriba: fila 0 = borde superior
    vmax = max(slip_d.max(), slip_k.max())

    ed = dyn_ellipse(cfg)
    ek = kin_ellipse(cfg)
    print(f"dinamico  misfit {mis_d:.4f}  elipse centro ({ed[0]:.1f}, {ed[1]:.1f}) km  "
          f"semiejes {ed[2]:.1f} x {ed[3]:.1f} km  phi {np.degrees(ed[4]):.0f} deg  slip max {slip_d.max():.2f} m")
    print(f"cinematico misfit {mis_k:.4f}  elipse centro ({ek[0]:.1f}, {ek[1]:.1f}) km  "
          f"semiejes {ek[2]:.1f} x {ek[3]:.1f} km  phi {np.degrees(ek[4]):.0f} deg  slip max {slip_k.max():.2f} m")

    # ---------- 1. instantaneas ----------
    # dynamic_slip acaba de correr fd3d en fd.WORK, asi que el slip-rate de ESTE
    # modelo esta ahi: no hay que buscarlo en el directorio de otra corrida.
    srz = read_tsn_fault_field(fd.WORK / "result" / "sliprateZ.res", fd.NXTT, fd.NZTT)
    nt = srz.shape[0]
    mr = np.abs(srz).sum(axis=(1, 2))
    # El muestreo va sobre la FASE ACTIVA (ultimo instante con slip-rate
    # apreciable), no sobre el momento acumulado: el frente para a ~5 s pero
    # queda una cola lenta que se lleva el ultimo 1% hasta los ~11 s, y
    # muestrear hasta ahi dejaba media figura en negro.
    peak = np.abs(srz).max(axis=(1, 2))
    active = np.nonzero(peak > 0.02 * peak.max())[0]
    tend = (active[-1] if active.size else nt - 1) * fd.DT
    times = np.linspace(0.3, max(tend, 1.0), 8)
    cum = np.cumsum(np.abs(srz), axis=0) * fd.DT
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.4), sharex=True, sharey=True)
    vr = np.percentile(np.abs(srz), 99.9)
    for ax, t in zip(axes.ravel(), times):
        i = min(int(t / fd.DT), nt - 1)
        # fd3d: eje 1 = rumbo, eje 2 = manteo desde el borde PROFUNDO -> voltear
        ax.imshow(np.abs(srz[i])[:, ::-1].T, origin="upper", extent=ext,
                  cmap="magma", vmin=0, vmax=vr, aspect="equal")
        ax.contour(np.linspace(0, Lx, fd.NXTT), np.linspace(0, Ly, fd.NZTT),
                   cum[i][:, ::-1].T, levels=[0.05], colors="#5ad7ff", linewidths=1.1)
        ax.plot(hx, hy, "*", color="#5ad7ff", ms=11, mec="k", mew=.5)
        ax.set_title(f"t = {i * fd.DT:.1f} s", fontsize=10)
    for ax in axes[-1]:
        ax.set_xlabel("Along strike (km)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Down dip (km)")
    fig.suptitle(f"Calama 2020 · {TAG} · rupture snapshots (slip rate; cyan line = rupture front)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(HERE / f"snapshots_{TAG}.png", dpi=135)

    # ---------- 2. comparacion de slip ----------
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    for ax, s, ttl, mis in ((axes[0], slip_d, "Dynamic (fd3d_TSN)", mis_d),
                            (axes[1], slip_k, "Kinematic (elliptical patch)", mis_k)):
        im = ax.imshow(s, origin="upper", extent=ext, cmap="viridis", vmin=0, vmax=vmax, aspect="equal")
        draw_ellipse(ax, ed[0], ed[1], ed[2], ed[3], ed[4], "#ff6b35", "-", "Dynamic asperity")
        draw_ellipse(ax, ek[0], ek[1], ek[2], ek[3], ek[4], "w", "--", "Kinematic patch")
        ax.plot(hx, hy, "*", color="w", ms=13, mec="k", mew=.6)
        ax.set_title(f"{ttl}   misfit {mis:.4f}", fontsize=10)
        ax.set_xlabel("Along strike (km)")
    axes[0].set_ylabel("Down dip (km)")
    axes[0].legend(loc="upper left", fontsize=8, framealpha=.85)
    fig.colorbar(im, ax=axes[1], label="Slip (m)")
    d = slip_d - slip_k
    m = np.abs(d).max()
    im2 = axes[2].imshow(d, origin="upper", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m, aspect="equal")
    axes[2].plot(hx, hy, "*", color="k", ms=13)
    axes[2].set_title("Dynamic − kinematic", fontsize=10)
    axes[2].set_xlabel("Along strike (km)")
    fig.colorbar(im2, ax=axes[2], label="Δ slip (m)")
    fig.suptitle(f"Calama 2020 · {TAG} · final slip, same mesh and same colour scale", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(HERE / f"slip_compare_{TAG}.png", dpi=135)

    # ---------- 3. perfiles ----------
    istk = np.linspace(0, Lx, fp.nx)
    idip = np.linspace(0, Ly, fp.ny)
    jh = int(np.clip(round(hy / Ly * fp.ny), 0, fp.ny - 1))
    ih = int(np.clip(round(hx / Lx * fp.nx), 0, fp.nx - 1))
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    axes[0].plot(istk, slip_d[jh], color="#ff6b35", lw=2, label=f"Dynamic (max {slip_d.max():.2f} m)")
    axes[0].plot(istk, slip_k[jh], color="#1f6f7a", lw=2, ls="--", label=f"Kinematic (max {slip_k.max():.2f} m)")
    axes[0].axvline(hx, color="0.5", lw=1, ls=":")
    axes[0].set_xlabel("Along strike (km)"); axes[0].set_title("Profile along strike, through the hypocentre", fontsize=10)
    axes[1].plot(idip, slip_d[:, ih], color="#ff6b35", lw=2)
    axes[1].plot(idip, slip_k[:, ih], color="#1f6f7a", lw=2, ls="--")
    axes[1].axvline(hy, color="0.5", lw=1, ls=":")
    axes[1].set_xlabel("Down dip (km)"); axes[1].set_title("Profile down dip, through the hypocentre", fontsize=10)
    for ax in axes:
        ax.set_ylabel("Slip (m)"); ax.grid(alpha=.25)
    axes[0].legend(fontsize=9)
    fig.suptitle(f"Calama 2020 · {TAG} · slip profiles", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(HERE / f"profiles_{TAG}.png", dpi=135)
    print("figuras:", ", ".join(f"{n}_{TAG}.png" for n in ("snapshots", "slip_compare", "profiles")))


if __name__ == "__main__":
    main()
