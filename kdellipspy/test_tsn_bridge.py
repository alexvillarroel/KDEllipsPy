"""End-to-end validation of the fd3d_TSN subprocess bridge: compiles a fresh
binary, runs a short forward, and checks the slip-rate reader against a raw
np.fromfile parse.

Skips if the fd3d_TSN clone or gfortran are unavailable.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic.tsn_bridge import TSNRunConfig, read_tsn_fault_field, run_tsn_forward

FD3D_TSN_SRC = Path("/home/alex/fd3d_TSN/src")
AMATRICE_EXAMPLE = Path("/home/alex/fd3d_TSN/examples/Amatrice2016")

needs_fd3d_tsn = pytest.mark.skipif(
    not FD3D_TSN_SRC.exists() or shutil.which("gfortran") is None,
    reason="fd3d_TSN clone or gfortran not available",
)


@pytest.fixture(scope="module")
def compiled_run_dir(tmp_path_factory):
    work_dir = tmp_path_factory.mktemp("fd3d_tsn_bridge_test")
    # inversion_com.f90 is pulled in via a textual INCLUDE (fd3d_init.f90),
    # not compiled directly, but still needs to be present in the cwd.
    for f in ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90", "inversion_com.f90"]:
        shutil.copy(FD3D_TSN_SRC / f, work_dir / f)
    shutil.copy(AMATRICE_EXAMPLE / "crustal.dat", work_dir / "crustal.dat")
    shutil.copy(AMATRICE_EXAMPLE / "inputinv.dat", work_dir / "inputinv.dat")

    binary = work_dir / "fd3d_gnu_TSN"
    subprocess.run(
        [
            "gfortran", "-DDIPSLIP", "-O2", "-cpp", "-o", str(binary),
            "fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90",
        ],
        cwd=work_dir,
        check=True,
        capture_output=True,
        text=True,
    )

    # Small/fast variant of the Amatrice2016 grid: 30 timesteps only.
    (work_dir / "inputfd3d.dat").write_text(
        "300 50 140       \\* Number of grid points in x,y,z direction (FD) *\\\n"
        "100.0             \\* Grid size in meters (FD) *\\\n"
        "30                \\* Number of time step (FD) *\\\n"
        "0.003             \\* Time step in seconds (FD) *\\\n"
        "45.               \\* Dip *\\\n"
        "10 6000. 25.      \\* nabc, vp, pml_fact *\\\n"
        "0.3               \\* damp_s *\\\n"
        "0\t\t  \\* number of seismic stations\n"
        "0\n"
    )
    return work_dir


@needs_fd3d_tsn
def test_run_tsn_forward_end_to_end(compiled_run_dir):
    nli, nwi = 13, 8
    grid = TSNFaultGridSpec(dh=100.0, dip_deg=45.0, nztT=140, nabc=10)
    model = np.array([4.0, 3.0, 6.0, 4.0, 0.3, 5.0, 1.05, 1.15, 1.5, 0.4], dtype=np.float32)
    run_cfg = TSNRunConfig(work_dir=compiled_run_dir, nxtT=300, nztT=140, dt_s=0.003)

    out = run_tsn_forward(model, nli, nwi, grid, run_cfg)

    assert set(out) == {"sliprateX", "sliprateZ"}
    for arr in out.values():
        assert arr.shape == (30, 300, 140)
        assert np.isfinite(arr).all()

    # Cross-check the reader against an independent raw parse.
    raw = np.fromfile(compiled_run_dir / "result" / "sliprateX.res", dtype=np.float32)
    assert raw.size == 30 * 300 * 140
    raw_reshaped = raw.reshape(30, 140, 300).transpose(0, 2, 1)
    assert np.array_equal(out["sliprateX"], raw_reshaped)


@needs_fd3d_tsn
def test_read_tsn_fault_field_rejects_bad_shape(compiled_run_dir, tmp_path):
    bogus = tmp_path / "bogus.res"
    np.zeros(100, dtype=np.float32).tofile(bogus)
    with pytest.raises(ValueError):
        read_tsn_fault_field(bogus, nxt=300, nzt=140)


@needs_fd3d_tsn
def test_rupture_stays_inside_ellipse(tmp_path):
    """Physics check of the params -> forwardmodel.dat -> fd3d_TSN chain: the
    asperity ruptures, the barrier around it does not. Regression for the
    legacy -1e8 Pa barrier, which fd3d_TSN (|traction| criterion) broke at
    t=0 over the whole fault.
    """
    for f in ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90", "inversion_com.f90"]:
        shutil.copy(FD3D_TSN_SRC / f, tmp_path / f)
    shutil.copy(AMATRICE_EXAMPLE / "crustal.dat", tmp_path / "crustal.dat")
    subprocess.run(
        ["gfortran", "-DDIPSLIP", "-O2", "-cpp", "-o", "fd3d_gnu_TSN",
         "fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90"],
        cwd=tmp_path, check=True, capture_output=True, text=True,
    )
    # 12 x 8 km fault, dh=200 m; coarse grid 13 x 9 -> 1 km between nodes.
    nxt, nzt, dh, nabc, dt, nt = 60, 40, 200.0, 10, 0.006, 700  # CFL = 6510*dt/dh ~ 0.2
    nli, nwi = 13, 9
    (tmp_path / "inputinv.dat").write_text(f"0\n{nli} {nwi}\n")
    (tmp_path / "inputfd3d.dat").write_text(
        f"{nxt} 20 {nzt}\n{dh}\n{nt}\n{dt}\n45.\n{nabc} 7000. 25.\n0.3\n0\n0\n")

    a, b, xo, yo = 3.0, 2.5, 7.0, 5.0  # coarse 1-based grid points (1 km each)
    model = np.array([a, b, xo, yo, 0.0, 10.0, 1.15, 1.1, 1.2, 0.3], dtype=np.float32)
    grid = TSNFaultGridSpec(dh=dh, dip_deg=45.0, nztT=nzt, nabc=nabc)
    out = run_tsn_forward(model, nli, nwi, grid, TSNRunConfig(work_dir=tmp_path, nxtT=nxt, nztT=nzt, dt_s=dt))
    slip = np.hypot(out["sliprateX"].sum(0), out["sliprateZ"].sum(0)) * dt  # (nxt, nzt)

    assert np.isfinite(slip).all()
    assert slip.max() > 0.1, "asperity did not rupture"

    # Fine node -> coarse 1-based coordinate (same mapping as inversion_modeltofd3d).
    x = np.arange(nxt)[:, None] * dh / (dh * nxt / (nli - 1)) + 1
    z = np.arange(nzt)[None, :] * dh / (dh * nzt / (nwi - 1)) + 1
    margin = 1.5  # bilinear interpolation smears the edge by ~1 coarse cell
    far_outside = ((x - xo) / (a + margin)) ** 2 + ((z - yo) / (b + margin)) ** 2 > 1
    assert far_outside.sum() > 0.3 * slip.size
    assert slip[far_outside].max() < 0.01 * slip.max(), "rupture escaped through the barrier"


@needs_fd3d_tsn
def test_slip_weakening_self_similarity(tmp_path):
    """Scaling prestress, strength and Dc by k scales slip by k with the same
    rupture timing (linear elastodynamics + linear slip-weakening). This is
    what lets the dynamic inversion impose an agency M0 by rescaling
    synthetics instead of re-running the solver (analogue of mt_strict).
    """
    for f in ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90", "inversion_com.f90"]:
        shutil.copy(FD3D_TSN_SRC / f, tmp_path / f)
    shutil.copy(AMATRICE_EXAMPLE / "crustal.dat", tmp_path / "crustal.dat")
    subprocess.run(
        ["gfortran", "-DDIPSLIP", "-O2", "-cpp", "-o", "fd3d_gnu_TSN",
         "fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90"],
        cwd=tmp_path, check=True, capture_output=True, text=True,
    )
    nxt, nzt, dh, nabc, dt, nt = 60, 40, 200.0, 10, 0.006, 700
    nli, nwi = 13, 9
    (tmp_path / "inputinv.dat").write_text(f"0\n{nli} {nwi}\n")
    (tmp_path / "inputfd3d.dat").write_text(
        f"{nxt} 20 {nzt}\n{dh}\n{nt}\n{dt}\n45.\n{nabc} 7000. 25.\n0.3\n0\n0\n")
    grid = TSNFaultGridSpec(dh=dh, dip_deg=45.0, nztT=nzt, nabc=nabc)
    run_cfg = TSNRunConfig(work_dir=tmp_path, nxtT=nxt, nztT=nzt, dt_s=dt)

    base = np.array([3.0, 2.5, 7.0, 5.0, 0.3, 4.0, 1.2, 1.1, 1.2, 0.4], dtype=np.float32)
    k = 2.5
    scaled = base.copy()
    scaled[5] *= k  # Te (strength = cte1*Te and nucleation = cte2*cte1*Te follow)
    scaled[9] *= k  # Dc
    rate1 = run_tsn_forward(base, nli, nwi, grid, run_cfg)["sliprateZ"]
    rate2 = run_tsn_forward(scaled, nli, nwi, grid, run_cfg)["sliprateZ"]

    assert np.abs(rate1).max() > 0.1
    assert np.abs(rate2 - k * rate1).max() < 1e-3 * np.abs(k * rate1).max()
