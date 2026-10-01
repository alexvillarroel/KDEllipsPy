"""Numerical validation of EllipticalStressMapper against the legacy Fortran.

Reference data: stressin.dat / peakin.dat written by the legacy fd3d_direct
run of 2026-07-13 (inversions/2026-05-25_Calama/Dynamic_inversion/Evento),
whose exact 10 input parameters are echoed in fd3d_direct.out.

Run with pytest or directly: python test_stress_mapper.py
"""

from pathlib import Path

import numpy as np
import pytest

from kdellipspy.core.geometry import (
    EllipticalStressMapper,
    TSNFaultGridSpec,
    build_tsn_dynamic_fields,
    write_fd3d_tsn_forwardmodel,
)

EVENTO = (
    Path(__file__).resolve().parents[1]
    / "inversions/2026-05-25_Calama/Dynamic_inversion/Evento"
)

# The 10 params printed by fd3d_direct.out (REAL*4 echo of the model that
# generated stressin.dat/peakin.dat). Grid: 160x160 (input.dat).
RMODEL = np.array(
    [
        57.8352814,   # a
        59.8246231,   # b
        101.124184,   # xo
        72.8244400,   # yo
        4.39626455,   # phi
        5.50441360,   # Te (MPa)
        1.05229867,   # cte1
        1.15203452,   # cte2
        7.65022230,   # r
        0.460266680,  # dmax (unused here)
    ],
    dtype=np.float32,
)
NXT = NYT = 160


needs_legacy_data = pytest.mark.skipif(
    not (EVENTO / "stressin.dat").exists(),
    reason=f"legacy reference data not found under {EVENTO}",
)


@needs_legacy_data
def test_prestress_matches_legacy_stressin():
    prestress, _ = EllipticalStressMapper(NXT, NYT).fields(RMODEL)
    # mkstress.f writes with j (dip) as the outer loop and i (strike) inner:
    # file order is [S(i,j) for j in 1..nyt for i in 1..nxt].
    ref = np.loadtxt(EVENTO / "stressin.dat", dtype=np.float64)
    assert ref.size == NXT * NYT
    ref_2d = ref.reshape(NYT, NXT).T  # -> [i, j]
    assert np.allclose(prestress, ref_2d, rtol=1e-5), (
        f"max abs diff = {np.max(np.abs(prestress - ref_2d)):.6g}, "
        f"n mismatched = {np.count_nonzero(~np.isclose(prestress, ref_2d, rtol=1e-5))}"
    )


@needs_legacy_data
def test_peak_matches_legacy_peakin():
    _, peak = EllipticalStressMapper(NXT, NYT).fields(RMODEL)
    # mkpeak.f writes with i (strike) as the outer loop and j (dip) inner:
    # file order is [P(i,j) for i in 1..nxt for j in 1..nyt].
    ref = np.loadtxt(EVENTO / "peakin.dat", dtype=np.float64)
    assert ref.size == NXT * NYT
    ref_2d = ref.reshape(NXT, NYT)  # -> [i, j]
    assert np.allclose(peak, ref_2d, rtol=1e-5), (
        f"max abs diff = {np.max(np.abs(peak - ref_2d)):.6g}"
    )


def _legacy_normstress_dipslip(k: int, dip_deg: float, dh: float, nzt: int, nfs: int) -> float:
    """Direct transcription of fd3d_init.f90:77-86 (DIPSLIP branch), indexed
    by the fine-grid k, used as ground truth for TSNFaultGridSpec."""
    dip_rad = np.radians(dip_deg)
    return max(1.0e5, 8520.0 * dh * float(nzt - nfs - k) * np.sin(dip_rad))


def test_normstress_matches_fortran_formula():
    # Amatrice2016 example: inputfd3d.dat -> nztT=140, dh=100, dip=45, nabc=10
    grid = TSNFaultGridSpec(dh=100.0, dip_deg=45.0, nztT=140, nabc=10)
    assert grid.nzt == 152  # nztT + nabc + nfs(=2)

    for k in [11, 50, 90, 130, 150]:  # k in [nabc+1, nzt-nfs]
        zs = grid.dh * (k - 1 - grid.nabc)
        got = float(grid.normstress_at_depth_m(np.array([zs]))[0])
        want = _legacy_normstress_dipslip(k, grid.dip_deg, grid.dh, grid.nzt, grid.nfs)
        assert got == pytest.approx(want, rel=1e-10), (k, got, want)


def test_build_tsn_dynamic_fields_shapes_and_conversion():
    nli, nwi = 25, 25
    grid = TSNFaultGridSpec(dh=100.0, dip_deg=45.0, nztT=140, nabc=10)
    model = np.array([8.0, 8.0, 12.0, 12.0, 0.3, 5.0, 1.05, 1.15, 3.0, 0.4], dtype=np.float32)

    t0, ts, dc = build_tsn_dynamic_fields(model, nli, nwi, grid)
    assert t0.shape == ts.shape == dc.shape == (nli, nwi)
    assert np.all(np.isfinite(t0)) and np.all(np.isfinite(ts)) and np.all(np.isfinite(dc))
    assert np.all(dc == np.float32(model[9]))

    # Inside the asperity TsI = peak_stress_Pa / normstress(depth); outside
    # (legacy negative prestress) it becomes a zero-prestress unbreakable
    # barrier, because fd3d_TSN compares |traction| against strength.
    t0_legacy, peak_pa = EllipticalStressMapper(nli, nwi).fields(model)
    normstress = grid.coarse_normstress_profile(nwi).astype(np.float32)
    inside = t0_legacy >= 0
    assert inside.any() and (~inside).any()
    assert np.allclose(ts[inside], (peak_pa / normstress[None, :])[inside], rtol=1e-6)
    assert np.all(t0[~inside] == 0) and np.all(t0 >= 0)
    assert np.all(ts[~inside] * normstress[None, :].repeat(nli, 0)[~inside] > 1e9)

    # forwardmodel.dat must round-trip through the writer without error.
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        write_fd3d_tsn_forwardmodel(Path(td) / "forwardmodel.dat", t0, ts, dc)


def test_forwardmodel_writer_roundtrip(tmp_path):
    rng = np.random.default_rng(0)
    t0 = rng.random((13, 8))
    ts = rng.random((13, 8))
    dc = rng.random((13, 8))
    out = tmp_path / "forwardmodel.dat"
    write_fd3d_tsn_forwardmodel(out, t0, ts, dc, header=(1.0, 2.0))
    # Read back exactly like fd3d_TSN: free-format, Fortran column-major order.
    vals = np.loadtxt(out).ravel()
    assert vals[0] == 1.0 and vals[1] == 2.0
    n = t0.size
    for arr, block in ((t0, vals[2 : 2 + n]), (ts, vals[2 + n : 2 + 2 * n]),
                       (dc, vals[2 + 2 * n :])):
        assert np.allclose(block.reshape(arr.shape, order="F"), arr, rtol=1e-4)


if __name__ == "__main__":
    test_prestress_matches_legacy_stressin()
    test_peak_matches_legacy_peakin()
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        test_forwardmodel_writer_roundtrip(Path(td))
    print("OK: EllipticalStressMapper matches legacy stressin.dat/peakin.dat")
