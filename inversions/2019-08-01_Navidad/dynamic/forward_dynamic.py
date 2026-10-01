"""Forward dinámico (fd3d_TSN) para Navidad 2019-08-01.

Compila fd3d_TSN (gfortran, DIPSLIP) en tsn_work/ si falta, escribe
inputfd3d.dat / inputinv.dat / crustal.dat desde ../$DYN_CASE/input.ctl (np1 por
defecto; p.ej. DYN_CASE=loc_mid), corre UN
modelo de 10 parámetros y compara sintéticos vs observados (misfit + figura).

Uso:  python forward_dynamic.py [a b xo yo phi Te cte1 cte2 r Dc]
      (a,b,xo,yo,r en puntos de la grilla gruesa Nx x Ny = 25x25 -> 2 km/pto)
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kdellipspy as kde
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicNAInversionModel, TSNRunConfig

HERE = Path(__file__).resolve().parent
CASE = os.environ.get("DYN_CASE", "np1")  # carpeta hermana con input.ctl y DATA/
CASE_DIR = HERE.parent / CASE
WORK = HERE / os.environ.get("DYN_WORK", "tsn_work")  # un directorio por corrida simultánea
FD3D_SRC = Path("/home/alex/fd3d_TSN/src")

# Grilla FD fina: 50 km x 50 km a dh=500 m -> 100 x 100 (divisible por 25x25).
# ponytail: dh=500 m resuelve la zona cohesiva con ~3-4 celdas para Dc~1 m;
# bajar a 250 m (200x200) si el slip-rate sale ruidoso.
DH, NXTT, NZTT, NYTT, NABC = 500.0, 100, 100, 30, 10
DT, NT = 0.015, 2667  # 40 s; fd3d_TSN exige CFL=Vpmax*dt/dh < 0.25 (Vp<=6.32 km/s en la falla)

# a, b, xo, yo, phi, Te[MPa], cte1, cte2, r, Dc[m]
DEFAULT_MODEL = [5.0, 3.0, 12.0, 12.0, 0.0, 15.0, 1.15, 1.1, 1.5, 1.0]


def setup_work_dir(cfg):
    WORK.mkdir(exist_ok=True)
    (WORK / "result").mkdir(exist_ok=True)
    if not (WORK / "fd3d_gnu_TSN").exists():
        srcs = ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90"]
        for f in srcs + ["inversion_com.f90"]:
            shutil.copy(FD3D_SRC / f, WORK / f)
        subprocess.run(["gfortran", "-DDIPSLIP", "-O3", "-march=native", "-cpp",
                        "-o", "fd3d_gnu_TSN", *srcs], cwd=WORK, check=True)

    dip = cfg.source_position.dip
    (WORK / "inputfd3d.dat").write_text(
        f"{NXTT} {NYTT} {NZTT}\n{DH}\n{NT}\n{DT}\n{dip}\n"
        f"{NABC} 7000. 25.\n0.3\n0\n0\n")
    nx, ny = cfg.fault_plane.nx, cfg.fault_plane.ny
    (WORK / "inputinv.dat").write_text(f"0\n{nx} {ny}\n")

    # input.ctl trae el tope de cada capa (m) en la columna "thickness".
    rows = "\n".join(f"{l.thickness/1e3:.3f} {l.vp/1e3:.3f} {l.vs/1e3:.3f} "
                     f"{l.rho/1e3:.3f} {l.qp:.1f} {l.qs:.1f}"
                     for l in cfg.velocity_model.layers)
    (WORK / "crustal.dat").write_text(
        f"Crustal model Navidad\nnumber of layers\n{len(cfg.velocity_model.layers)}\n"
        f"Parameters of the layers\ndepth(km) Vp Vs Rho Qp Qs\n{rows}\n")


def main():
    model = np.array([float(v) for v in sys.argv[1:]] or DEFAULT_MODEL, dtype=np.float32)
    assert model.size == 10, "se esperan 10 parámetros"

    cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
    setup_work_dir(cfg)

    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE_DIR / "input.ctl"), data_dir=str(CASE_DIR / "DATA"))

    inv = DynamicNAInversionModel(
        config=cfg,
        observed_waveforms=observed,
        time_array=time_array,
        tsn_run_cfg=TSNRunConfig(work_dir=WORK, nxtT=NXTT, nztT=NZTT, dt_s=DT),
        tsn_grid=TSNFaultGridSpec(dh=DH, dip_deg=cfg.source_position.dip, nztT=NZTT, nabc=NABC),
    )
    misfit = inv.objective_function(model)
    syn = inv.best_synthetics
    print(f"misfit = {misfit:.4f}")

    # Slip final y M0 desde la salida de fd3d_TSN (dip-slip -> sliprateZ).
    from kdellipspy.inversion.dynamic.tsn_bridge import read_tsn_fault_field
    srz = read_tsn_fault_field(WORK / "result" / "sliprateZ.res", NXTT, NZTT)
    slip = srz.sum(axis=0) * DT
    mu = inv.dynamic_fm._mu_pa
    m0 = mu * np.abs(slip).sum() * DH**2
    print(f"slip max = {np.abs(slip).max():.2f} m   M0 = {m0:.3e} Nm   Mw = {(np.log10(m0) - 9.1) / 1.5:.2f}")

    if syn is None:
        sys.exit("forward falló (ver log)")

    names = [s.name for s in cfg.stations.stations]
    fig, axes = plt.subplots(len(names), 3, figsize=(10, 1.6 * len(names)), sharex=True)
    for i, name in enumerate(names):
        for c, comp in enumerate("XYZ"):
            ax = axes[i, c]
            ax.plot(time_array, observed[i, c], "k", lw=1)
            ax.plot(time_array, syn[i, c], "r", lw=1)
            ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(name, rotation=0, ha="right")
            if i == 0:
                ax.set_title(comp)
    fig.suptitle(f"Navidad {CASE} dinámico — misfit {misfit:.3f}  model={model.tolist()}", fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "forward_dynamic_traces.png", dpi=120)

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(slip.T, origin="lower", extent=[0, NXTT * DH / 1e3, 0, NZTT * DH / 1e3], cmap="hot_r")
    fig.colorbar(im, label="slip (m)")
    ax.set_xlabel("strike (km)"); ax.set_ylabel("dip (km, desde abajo)")
    fig.savefig(HERE / "forward_dynamic_slip.png", dpi=120)
    inv.clean()


if __name__ == "__main__":
    main()
