"""Velocidad de ruptura del modelo dinámico, desde ruptime.res de fd3d_TSN.

En el dinámico Vr NO es un parámetro libre: emerge de la fricción. Esto la mide
sobre el modelo que ganó, que es la forma de contestar la pregunta del
supershear que el cinemático dejó abierta (ver kin_tests/VR_NO_RESUELTO.md: el
valle Vr-dt0 no tiene mínimo, así que la elipse no puede decir nada).

`ruptime.res` es TEXTO (un valor por celda de la malla fina, nxtT x nztT), no
binario: read_tsn_fault_field falla con él. 10000.0 es el centinela de "no
rompió".

Vr se reporta de dos formas, porque no son lo mismo:
  - media radial: mediana de |r| / t sobre las celdas rotas, con r medido desde
    el hipocentro. Es la velocidad media desde la nucleación.
  - gradiente local: 1 / |grad t|, la velocidad instantánea del frente. Es la
    que hay que comparar con beta, y la que mostraría supershear localizado.

Uso:  DYN_CASE=<caso> python rupture_velocity.py <tsn_work_dir>
"""

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import kdellipspy as kde  # noqa: E402
import forward_dynamic as fd  # noqa: E402
from kdellipspy.core.geometry import tsn_hypocentre_coarse  # noqa: E402

NOT_RUPTURED = 9999.0  # el centinela es 10000; cualquier cosa por encima no rompió


def load_ruptime(path, nxt, nzt):
    t = np.loadtxt(path, dtype=float)
    if t.size != nxt * nzt:
        raise ValueError(f"{path}: {t.size} valores, se esperaban nxt*nzt={nxt * nzt}")
    return t.reshape(nxt, nzt)


def beta_at_fault(cfg):
    """Vs en la capa que contiene el hipocentro."""
    z = cfg.source_position.depth
    layers = [(l.thickness / 1e3, l.vs / 1e3) for l in cfg.velocity_model.layers]
    return [v for top, v in layers if top <= z][-1]


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else fd.WORK
    cfg = kde.ConfigParser(str(fd.CASE_DIR / "input.ctl"))
    nxt, nzt, dh = fd.NXTT, fd.NZTT, fd.DH
    tr = load_ruptime(work / "result" / "ruptime.res", nxt, nzt)
    ruptured = tr < NOT_RUPTURED
    if not ruptured.any():
        raise SystemExit("ninguna celda rompió")

    # Hipocentro en la malla fina. tsn_hypocentre_coarse da la grilla gruesa
    # (1-based) y el manteo se cuenta desde el borde profundo, igual que fd3d.
    hx_c, hy_c = tsn_hypocentre_coarse(cfg.fault_plane)
    fx, fz = nxt / cfg.fault_plane.nx, nzt / cfg.fault_plane.ny
    hx, hz = (hx_c - 0.5) * fx, (hy_c - 0.5) * fz

    ix, iz = np.meshgrid(np.arange(nxt), np.arange(nzt), indexing="ij")
    r = np.hypot((ix - hx) * dh, (iz - hz) * dh) / 1e3  # km

    beta = beta_at_fault(cfg)
    m = ruptured & (tr > 0.05) & (r > 1.0)   # fuera del parche de nucleación
    v_rad = r[m] / tr[m]

    # Velocidad local del frente: 1/|grad t|. El gradiente se toma solo donde
    # las cuatro vecinas rompieron, para no cruzar el borde de la aspereza.
    t_fill = np.where(ruptured, tr, np.nan)
    gx, gz = np.gradient(t_fill, dh / 1e3)
    g = np.hypot(gx, gz)
    ok = np.isfinite(g) & (g > 1e-6) & m
    v_loc = 1.0 / g[ok]

    print(f"caso {fd.CASE}   work {work.name}")
    print(f"beta en el hipocentro        : {beta:.3f} km/s")
    print(f"celdas rotas                 : {ruptured.sum()} de {nxt * nzt} "
          f"({100 * ruptured.mean():.1f}%, {ruptured.sum() * (dh / 1e3) ** 2:.0f} km2)")
    print(f"duracion de la ruptura       : {tr[ruptured].max():.2f} s")
    for label, v in (("Vr media radial", v_rad), ("Vr local (1/|grad t|)", v_loc)):
        q = np.percentile(v, [5, 25, 50, 75, 95])
        print(f"{label:<29}: mediana {np.median(v):5.2f} km/s = {np.median(v) / beta:4.2f} beta"
              f"   p5-p95 {q[0]:.2f}-{q[4]:.2f}  ({q[0] / beta:.2f}-{q[4] / beta:.2f} beta)")
        print(f"{'':29}  fraccion supershear (>sqrt(2) beta = {np.sqrt(2) * beta:.2f}): "
              f"{100 * np.mean(v > np.sqrt(2) * beta):.1f}%")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ext = [0, nxt * dh / 1e3, 0, nzt * dh / 1e3]
    im = axes[0].imshow(np.where(ruptured, tr, np.nan).T, origin="lower", extent=ext, cmap="viridis")
    cs = axes[0].contour(np.linspace(0, ext[1], nxt), np.linspace(0, ext[3], nzt),
                         np.where(ruptured, tr, np.nan).T, levels=np.arange(0, 12, 1.0),
                         colors="w", linewidths=0.6)
    axes[0].clabel(cs, fmt="%.0f s", fontsize=7)
    axes[0].plot(hx * dh / 1e3, hz * dh / 1e3, "r*", ms=13)
    fig.colorbar(im, ax=axes[0], label="Rupture time (s)")
    axes[0].set_xlabel("Along strike (km)")
    axes[0].set_ylabel("Along dip (km, from bottom edge)")
    axes[0].set_title("Rupture time and isochrones")

    axes[1].hist(v_loc, bins=60, range=(0, 2.0 * beta), color="0.4")
    for frac, lab, st in ((0.9, r"0.9$\beta$", ":"), (1.0, r"$\beta$", "--"),
                          (np.sqrt(2), r"$\sqrt{2}\beta$", "--")):
        axes[1].axvline(frac * beta, color="r", ls=st, lw=1.2)
        axes[1].text(frac * beta, axes[1].get_ylim()[1] * 0.95, " " + lab, color="r", fontsize=8)
    axes[1].axvspan(beta, np.sqrt(2) * beta, color="r", alpha=0.10)
    axes[1].set_xlabel(r"Local rupture velocity $1/|\nabla t|$ (km s$^{-1}$)")
    axes[1].set_ylabel("Cells")
    axes[1].set_title("Rupture-velocity distribution")
    fig.suptitle(f"Calama 2020 · {fd.CASE.split('/')[-1][:24]} · dynamic rupture velocity", fontsize=10)
    fig.tight_layout()
    png = work.parent / f"rupture_velocity_{work.name}.png"
    fig.savefig(png, dpi=140)
    print(f"figura: {png}")


if __name__ == "__main__":
    main()
