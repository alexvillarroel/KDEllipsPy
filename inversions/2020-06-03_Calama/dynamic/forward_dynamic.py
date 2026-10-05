"""Forward dinámico (fd3d_TSN) para Calama 2020-06-03 (intraslab, falla NORMAL).

Compila fd3d_TSN (gfortran, DIPSLIP + FSPACE: falla enterrada, bordes absorbentes
arriba) en tsn_work/ si falta, escribe
inputfd3d.dat / inputinv.dat / crustal.dat desde ../$DYN_CASE/input.ctl (np1 por
defecto; p.ej. DYN_CASE=loc_mid), corre UN
modelo de 10 parámetros y compara sintéticos vs observados (misfit + figura).

Uso:  python forward_dynamic.py [a b xo yo phi Te cte1 cte2 r Dc]
      (a,b,xo,yo,r en puntos de la grilla gruesa Nx x Ny = 40x40 -> 1 km/pto)

Diferencias respecto de Navidad (ver docs/dynamic_rupture_guide.md, sección 5):
  - El mecanismo es NORMAL (np1 rake -85, np2 rake -101). fd3d_TSN solo recibe
    el buzamiento y la aspereza tiene prestress +Te, así que SIEMPRE desliza
    hacia +Z local (rake +90, inverso). DynamicForwardModel.slip_sign deriva el
    signo del rake del input.ctl y niega (Mx, Mz) = rake + 180. NO se fuerza el
    sentido con Te negativo: la lógica de barrera trata T0 < 0 como barrera.
  - CFL: a ~100-133 km de profundidad el modelo de Potin llega a Vp 8.415 km/s,
    con lo que dt = 0.015 s (Navidad) da CFL 0.2524 > 0.25. DT se calcula aquí
    desde el propio modelo de velocidades y se verifica con un assert.
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
from kdellipspy.core.geometry import TSNFaultGridSpec  # noqa: F401  (re-export para scripts)
from kdellipspy.inversion.dynamic import DynamicNAInversionModel, TSNRunConfig

HERE = Path(__file__).resolve().parent
CASE = os.environ.get("DYN_CASE", "np1")  # carpeta hermana con input.ctl y DATA/
# p.ej. DYN_CASE=kin_tests/np1_isola_potin_b0.04-0.15 (el caso cinemático elegido)
CASE_DIR = HERE.parent / CASE
CASE_TAG = CASE.replace("/", "__")  # nombre plano para directorios de salida
WORK = HERE / os.environ.get("DYN_WORK", "tsn_work")  # un directorio por corrida simultánea
FD3D_SRC = Path("/home/alex/fd3d_TSN/src")

# Grilla FD fina a dh=500 m sobre la falla del input.ctl (Lx x Ly), divisible
# por la malla de subfallas (nx x ny).
# ponytail: dh=500 m resuelve la zona cohesiva con ~3-4 celdas para Dc~1 m;
# bajar a 250 m si el slip-rate sale ruidoso.
DH, NYTT, NABC = 500.0, 30, 10
NFS = NABC  # -DFSPACE: fd3d_init.f90 usa nfs=nabc (sin superficie libre)
BINARY = "fd3d_gnu_TSN_fspace"
_FP = kde.ConfigParser(str(CASE_DIR / "input.ctl")).fault_plane
NXTT, NZTT = int(round(_FP.lx / DH)), int(round(_FP.ly / DH))
assert NXTT % _FP.nx == 0 and NZTT % _FP.ny == 0, "Lx/dh y Ly/dh deben ser múltiplos de nx, ny"
DURATION_S = 40.0
CFL_MAX = 0.25  # fd3d_TSN exige Vpmax*dt/dh < 0.25


def cfl_dt_max(cfg):
    """dt máximo admisible: el Vp más alto que ve fd3d_TSN (no el de la
    profundidad del hipocentro) sale de las mismas capas que escribe
    crustal_rows."""
    return CFL_MAX * DH / 1e3 / max(r[1] for r in crustal_rows(cfg))


DT = 0.0145  # CFL 0.244 con Vp 8.415 km/s (Potin et al. 2024 a ~112 km)
NT = int(round(DURATION_S / DT))

# a, b, xo, yo, phi, Te[MPa], cte1, cte2, r, Dc[m]
# a,b en ptos de la grilla gruesa (2 km/pto en dyn_cases): 3 ptos = 6 km de
# semieje, que es lo que da el cinemático (a1 5.9, a2 6.3 km). xo, yo van en
# el hipocentro -- tsn_hypocentre_coarse(_FP) = (10.5, 10.5) con 20x20.
DEFAULT_MODEL = [3.0, 3.0, 10.5, 10.5, 0.0, 15.0, 1.15, 1.1, 1.5, 1.0]


def tsn_grid(cfg):
    return TSNFaultGridSpec(dh=DH, dip_deg=cfg.source_position.dip, nztT=NZTT, nabc=NABC, nfs=NFS)


def tsn_run_cfg():
    return TSNRunConfig(work_dir=WORK, nxtT=NXTT, nztT=NZTT, dt_s=DT, binary=BINARY)


def fault_top_depth_km(cfg):
    """Profundidad del borde superior de la falla (Hy se mide desde ahí)."""
    sp, fp = cfg.source_position, cfg.fault_plane
    return sp.depth - fp.hy / 1e3 * np.sin(np.radians(sp.dip))


def crustal_rows(cfg):
    """Capas (tope km, Vp, Vs, rho, Qp, Qs) desplazadas para que la profundidad
    0 de fd3d_TSN (fila superior de la falla) sea el borde superior real."""
    ztop = fault_top_depth_km(cfg)
    layers = [(l.thickness / 1e3 - ztop, l.vp / 1e3, l.vs / 1e3, l.rho / 1e3, l.qp, l.qs)
              for l in cfg.velocity_model.layers]
    first = max(i for i, l in enumerate(layers) if l[0] <= 0.0)  # capa que contiene el borde superior
    return [(0.0,) + layers[first][1:]] + layers[first + 1:]


def setup_work_dir(cfg):
    dt_max = cfl_dt_max(cfg)
    print(f"[CFL] Vp max {max(r[1] for r in crustal_rows(cfg)):.3f} km/s  dh {DH:.0f} m  "
          f"dt {DT:.4f} s  ->  CFL {DT * max(r[1] for r in crustal_rows(cfg)) / (DH / 1e3):.4f} "
          f"(dt_max {dt_max:.5f} s, NT {NT} = {NT * DT:.1f} s)", flush=True)
    assert DT < dt_max, f"CFL violado: dt={DT} >= dt_max={dt_max:.5f} (bajar DT y subir NT)"
    WORK.mkdir(exist_ok=True)
    (WORK / "result").mkdir(exist_ok=True)
    if not (WORK / BINARY).exists():
        srcs = ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90"]
        for f in srcs + ["inversion_com.f90"]:
            shutil.copy(FD3D_SRC / f, WORK / f)
        subprocess.run(["gfortran", "-DDIPSLIP", "-DFSPACE", "-O3", "-march=native", "-cpp",
                        "-o", BINARY, *srcs], cwd=WORK, check=True)

    dip = cfg.source_position.dip
    (WORK / "inputfd3d.dat").write_text(
        f"{NXTT} {NYTT} {NZTT}\n{DH}\n{NT}\n{DT}\n{dip}\n"
        f"{NABC} 7000. 25.\n0.3\n0\n0\n")
    nx, ny = cfg.fault_plane.nx, cfg.fault_plane.ny
    (WORK / "inputinv.dat").write_text(f"0\n{nx} {ny}\n")

    # input.ctl trae el tope de cada capa (m) en la columna "thickness";
    # se desplaza al borde superior real de la falla (ver crustal_rows).
    layers = crustal_rows(cfg)
    rows = "\n".join(f"{t:.3f} {vp:.3f} {vs:.3f} {rho:.3f} {qp:.1f} {qs:.1f}" for t, vp, vs, rho, qp, qs in layers)
    (WORK / "crustal.dat").write_text(
        f"Crustal model Calama 2020 (z=0 en el borde superior de la falla, {fault_top_depth_km(cfg):.2f} km)\n"
        f"number of layers\n{len(layers)}\nParameters of the layers\ndepth(km) Vp Vs Rho Qp Qs\n{rows}\n")


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
        tsn_run_cfg=tsn_run_cfg(),
        tsn_grid=tsn_grid(cfg),
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
    fig.suptitle(f"Calama 2020 · {CASE} · dynamic forward — misfit {misfit:.3f}  rake sign {inv.dynamic_fm.slip_sign:+.0f}\nmodel={[round(float(v), 3) for v in model]}", fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "forward_dynamic_traces.png", dpi=120)

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(slip.T, origin="lower", extent=[0, NXTT * DH / 1e3, 0, NZTT * DH / 1e3], cmap="hot_r")
    fig.colorbar(im, label="slip (m)")
    ax.set_xlabel("Along strike (km)"); ax.set_ylabel("Along dip (km, from bottom edge)")
    ax.set_title("Final slip")
    fig.savefig(HERE / "forward_dynamic_slip.png", dpi=120)
    inv.clean()


if __name__ == "__main__":
    main()
