"""Resumen del mejor modelo de una inversión dinámica (na_dynamic.py).

Uso:  DYN_CASE=loc_mid python report_dynamic.py [na_output_dir]
Toma el mejor modelo (na_result.joblib si existe; si no, best_model_live.txt),
re-corre fd3d_TSN en tsn_work_report/ (no toca corridas en curso) y escribe en
<na_output_dir>/report/: traces.png, slip.png, summary.txt.
"""

import os
import sys
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("DYN_WORK", "tsn_work_report")
import forward_dynamic as fd  # noqa: E402
import kdellipspy as kde  # noqa: E402
import na_dynamic as nd  # noqa: E402
from kdellipspy.core.geometry import TSNFaultGridSpec  # noqa: E402
from kdellipspy.inversion.dynamic import TSNRunConfig  # noqa: E402

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else fd.HERE / f"na_output_{fd.CASE_TAG}"


def best_free_model():
    res = OUT / "na_result.joblib"
    if res.exists():
        r = joblib.load(res)
        return np.asarray(r.best_model.model, float), float(r.best_model.misfit), "na_result.joblib"
    lines = [l.split() for l in open(OUT / "best_model_live.txt") if not l.startswith("#") and l.strip()]
    vals = {l[0]: float(l[-1]) for l in lines}
    misfit = float(open(OUT / "best_model_live.txt").read().split("misfit=")[1].split()[0])
    return np.array([vals[n.split()[0]] for n in nd.NAMES]), misfit, "best_model_live.txt (en curso)"


def main():
    rep = OUT / "report"
    rep.mkdir(exist_ok=True)
    cfg = kde.ConfigParser(str(fd.CASE_DIR / "input.ctl"))
    fd.setup_work_dir(cfg)
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(fd.CASE_DIR / "input.ctl"), data_dir=str(fd.CASE_DIR / "DATA"))
    inv = nd.FreeSubsetNA(
        config=cfg, observed_waveforms=observed, time_array=time_array,
        tsn_run_cfg=fd.tsn_run_cfg(),
        tsn_grid=fd.tsn_grid(cfg),
    )
    inv.dynamic_fm.m0_target = nd.M0_TARGET
    free, misfit_na, source = best_free_model()
    misfit, syn = inv._evaluate_model(free)
    full = nd.to_full_model(free)
    k = inv.dynamic_fm.last_m0_scale

    rate = inv.dynamic_fm._last_run[1]
    slip = np.hypot(rate["sliprateX"].sum(0), rate["sliprateZ"].sum(0)) * fd.DT * k
    mu = inv.dynamic_fm._mu_pa
    m0 = mu * slip.sum() * fd.DH**2
    mw = (np.log10(m0) - 9.1) / 1.5
    rupt = slip > 0.05 * slip.max()
    area = rupt.sum() * fd.DH**2
    r_eq = np.sqrt(area / np.pi)
    te, dc = full[5] * k, full[9] * k
    pts_km = fd.NXTT * fd.DH / 1e3 / (cfg.fault_plane.nx - 1)  # tamaño de un punto de la grilla gruesa

    lines = [
        f"Caso {fd.CASE}   modelo desde {source}",
        f"misfit NA {misfit_na:.4f}   (re-evaluado {misfit:.4f})",
        "",
        *[f"  {n:<12s} {v:10.4f}" for n, v in zip(nd.NAMES, free)],
        "",
        f"Semiejes           {full[0]*pts_km:.1f} x {full[1]*pts_km:.1f} km   phi {np.degrees(full[4]):.0f} deg",
        f"Te (stress drop)   {te:.2f} MPa" + (f"   (Te_ref {nd.TE_REF} x k={k:.3f}, M0 impuesto)" if nd.M0_TARGET else ""),
        f"Resistencia pico   {te*full[6]:.2f} MPa   (exceso {te*(full[6]-1):.2f} MPa)",
        f"Dc                 {dc:.2f} m",
        f"dt0                {free[nd.NAMES.index('dt0 (s)')]:+.2f} s",
        f"M0                 {m0:.3e} N·m   Mw {mw:.2f}",
        f"Slip max / medio   {slip.max():.2f} / {slip[rupt].mean():.2f} m",
        f"Area rota          {area/1e6:.0f} km2   (r_eq {r_eq/1e3:.1f} km)",
        f"Eshelby 7/16 M0/r^3 {7/16*m0/r_eq**3/1e6:.2f} MPa",
    ]
    (rep / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    names = [s.name for s in cfg.stations.stations]
    fig, ax = plt.subplots(len(names), 3, figsize=(10, 1.6 * len(names)), sharex=True)
    for i, n in enumerate(names):
        for c, comp in enumerate("XYZ"):
            ax[i, c].plot(time_array, observed[i, c], "k", lw=1)
            ax[i, c].plot(time_array, syn[i, c], "r", lw=1)
            ax[i, c].set_yticks([])
            if c == 0:
                ax[i, c].set_ylabel(n, rotation=0, ha="right")
            if i == 0:
                ax[i, c].set_title(comp)
    fig.suptitle(f"Dinámico {fd.CASE} — misfit {misfit:.3f}  (negro obs, rojo sint)", fontsize=9)
    fig.tight_layout()
    fig.savefig(rep / "traces.png", dpi=120)

    # fd3d_TSN cuenta las filas del dip desde la más profunda; AXITRA (y
    # bin_slip_rate_to_subfaults) desde la más somera: se invierte para graficar
    # con el borde superior arriba, igual que la malla de AXITRA.
    fig, axs = plt.subplots(figsize=(5.5, 5))
    ext = [0, fd.NXTT * fd.DH / 1e3, fd.NZTT * fd.DH / 1e3, 0]
    im = axs.imshow(slip[:, ::-1].T, extent=ext, cmap="hot_r", aspect="auto")
    fig.colorbar(im, label="slip (m)")
    fp = cfg.fault_plane
    axs.plot(fp.hx / 1e3, fp.hy / 1e3, "b*", ms=14, label="hipocentro")
    axs.legend(loc="lower right", fontsize=8)
    axs.set_xlabel("strike (km)")
    axs.set_ylabel("dip (km desde el borde superior, como en AXITRA)")
    axs.set_title(f"Mw {mw:.2f}  Te {te:.1f} MPa  slip max {slip.max():.1f} m", fontsize=9)
    fig.tight_layout()
    fig.savefig(rep / "slip.png", dpi=120)
    inv.clean()


if __name__ == "__main__":
    main()
