"""Validation of the dynamic slip-rate -> axitra synthetics bridge
(kdellipspy/inversion/dynamic/dynamic_convolution.py).

Two levels:
1. Pure-NumPy checks of the binning/rake-projection arithmetic (no axitra).
2. A real cross-check against the compiled axitra binary: for a separable
   (spatially-weighted x shared-time-shape) moment-rate field, summing one
   conv() call per subfault must exactly reproduce a single conv() call with
   all subfaults active and the same shared sfunc — the core assumption
   convolve_dynamic_sources relies on (axitra applies one sfunc per call,
   scaled per-source by hist's moment column).
"""

from pathlib import Path

import numpy as np
import pytest

from kdellipspy.inversion.dynamic.dynamic_convolution import (
    bin_slip_rate_to_subfaults,
    convolve_dynamic_sources,
    project_rake,
    resample_to_axitra_grid,
    axitra_aw,
    response_basis,
    synthetics_from_basis,
)

AXITRA_DIR = Path(__file__).resolve().parent / "axitra" / "src"
needs_axitra = pytest.mark.skipif(
    not (AXITRA_DIR / "axitra").exists(), reason="compiled axitra binary not found"
)


# ---------------------------------------------------------------------------
# Pure-NumPy arithmetic checks
# ---------------------------------------------------------------------------

def test_bin_slip_rate_to_subfaults_sums_moment():
    nt, nxt, nzt = 5, 4, 2
    # Uniform slip_x=1 m/s everywhere, slip_z=0.
    slip_x = np.ones((nt, nxt, nzt))
    slip_z = np.zeros((nt, nxt, nzt))
    dh, mu = 100.0, 3.0e10

    mrate_x, mrate_z = bin_slip_rate_to_subfaults(slip_x, slip_z, nx_sub=2, nz_sub=1, dh_fine_m=dh, mu_pa=mu)
    assert mrate_x.shape == (2, nt)
    # Each subfault covers (nxt/2)*(nzt/1) = 2*2 = 4 fine cells.
    expected = mu * dh**2 * 4  # sum(slip_x=1 over 4 cells) * mu * dA
    assert np.allclose(mrate_x, expected)
    assert np.allclose(mrate_z, 0.0)


def test_bin_slip_rate_to_subfaults_rejects_non_divisible_grid():
    slip = np.zeros((2, 5, 4))
    with pytest.raises(ValueError):
        bin_slip_rate_to_subfaults(slip, slip, nx_sub=2, nz_sub=1, dh_fine_m=100.0, mu_pa=1.0)


def test_project_rake_recovers_known_direction():
    # Subfault 0: pure X slip. Subfault 1: pure Z slip. Subfault 2: 45 deg.
    nt = 10
    shape = np.linspace(0, 1, nt)
    mrate_x = np.stack([shape, np.zeros(nt), shape])
    mrate_z = np.stack([np.zeros(nt), shape, shape])

    moment_rate, rake_deg = project_rake(mrate_x, mrate_z)
    assert rake_deg == pytest.approx([0.0, 90.0, 45.0], abs=1e-6)
    # Time integral of the projected scalar must equal |Mx_tot, Mz_tot|.
    integral = moment_rate.sum(axis=1)
    expected_mag = np.hypot(mrate_x.sum(axis=1), mrate_z.sum(axis=1))
    assert integral == pytest.approx(expected_mag, rel=1e-10)


def test_resample_to_axitra_grid_zero_pads_and_interpolates():
    moment_rate = np.array([[0.0, 1.0, 2.0, 3.0]])  # dt_fine=1.0 -> t=[0,1,2,3]
    out = resample_to_axitra_grid(moment_rate, dt_fine_s=1.0, npt_axitra=8, dt_axitra_s=0.5)
    assert out.shape == (1, 8)
    # t_axitra = [0,.5,1,1.5,2,2.5,3,3.5] -> linear interp of the ramp, then 0 past t=3.
    assert out[0] == pytest.approx([0, 0.5, 1, 1.5, 2, 2.5, 3, 0.0])


# ---------------------------------------------------------------------------
# Real cross-check against the compiled axitra binary
# ---------------------------------------------------------------------------

def _mock_config():
    from kdellipspy import ConfigParser

    # NOTE: ConfigParser.from_dict resolves keys via a fuzzy matcher that
    # requires the *long* descriptive names below (see config_parser.py's
    # `_get_param_value` / `FaultPlaneParams.from_dict` etc.) — short keys
    # like "nx"/"Lx" silently fall back to defaults (nx=ny=1!) instead of
    # raising, which is how test_forward_simple.py's mock config ended up
    # silently building a 1x1 mesh too.
    params = {
        "source_position": {
            "event_name": "TestEvent", "latitude": -33.45, "longitude": -70.66,
            "depth": 10.0, "strike": 20.0, "dip": 60.0, "rake": 90.0,
        },
        "fault_plane": {
            "Length along strike (Lx)": 20000.0,
            "Length along dip (Ly)": 20000.0,
            "Number of subfaults along strike (Nx)": 2,
            "Number of subfaults along dip (Ny)": 2,
            "Hypocenter position strike (Hx)": 10000.0,
            "Hypocenter position dip (Hy)": 10000.0,
        },
        "ellipse": {
            "Number of ellipses": 1, "Initial slip": 0, "Slip shape": 1,
            "Frequency 1 (Freq1)": 0.1, "Frequency 2 (Freq2)": 1.0, "Time shift (T0)": 2.0,
        },
        "observed_data": {
            "Time window start (t1)": -5.0, "Time window end (t2)": 30.0,
            "Number of points (Npts)": 256, "Delta / Time step": 0.1, "Units": 1,
        },
        # SI units, "thickness" = layer-top depth (m), same as input.ctl.
        # (km / g/cm3*1e-6 here silently gave all-zero Green's functions.)
        "velocity_model": [
            {"thickness": 0.0, "vp": 5000.0, "vs": 3000.0, "rho": 2500.0, "qp": 500.0, "qs": 200.0},
            {"thickness": 5000.0, "vp": 6500.0, "vs": 3800.0, "rho": 2800.0, "qp": 1000.0, "qs": 500.0},
        ],
        "moment_tensor": {"Moment Tensor Flag": 0},
        "stations": [
            {"name": "STA1", "latitude": -33.40, "longitude": -70.66, "height": 0.0, "use_n": True, "use_e": True, "use_z": True},
            {"name": "STA2", "latitude": -33.45, "longitude": -70.60, "height": 0.0, "use_n": True, "use_e": True, "use_z": True},
        ],
        "inversion_params": [],
    }
    return ConfigParser.from_dict(params)


@needs_axitra
def test_per_subfault_calls_reproduce_single_shared_sfunc_call():
    from kdellipspy import AxitraForwardModel

    cfg = _mock_config()
    fm = AxitraForwardModel.from_config(cfg, axitra_dir=str(AXITRA_DIR))
    geom = fm.build_geometry()
    # keep_all_sources=True -> all nx*ny=4 subfaults present regardless of slip.
    model_params = np.array([5.0, 3.0, 0.0, 0.5, 0.0, 0.001, 2.5])
    geom = fm.apply_ellipse_model_to_geometry(geom, model_params, keep_all_sources=True)
    nsub = geom.nsources
    assert nsub == 4

    ap = fm.build_axitra(geom, latlon=False)
    ap = fm.green(ap, quiet=True)

    from axitra import moment as axitra_moment  # available now: fm.green() added axpath to sys.path

    rng = np.random.default_rng(0)
    # axitra.py writes .hist and .sou with fmt="%.3g": pick weights (1 or 2 x 1e18 Nm)
    # and a 2-significant-digit shape so both sides are exact in 3 digits.
    weights = rng.integers(1, 3, size=nsub) * 1e18
    shared_shape = np.zeros(ap.npt)
    # A smooth causal-ish pulse, arbitrary but nonzero on a good chunk of the window.
    n_pulse = min(30, ap.npt)
    shared_shape[:n_pulse] = np.round(np.hanning(n_pulse), 2)
    strike, dip, rake = cfg.source_position.strike, cfg.source_position.dip, cfg.source_position.rake

    # Reference: ONE conv() call, all subfaults active with their own weight,
    # sharing the same sfunc.
    hist_ref = np.zeros((nsub, 8), dtype=np.float64)
    hist_ref[:, 0] = np.arange(1, nsub + 1)
    hist_ref[:, 1] = weights
    hist_ref[:, 2] = strike
    hist_ref[:, 3] = dip
    hist_ref[:, 4] = rake
    _, rx, ry, rz = axitra_moment.conv(ap, hist_ref, source_type=3, t0=0.0, unit=1, sfunc=shared_shape)
    rx, ry, rz = rx.copy(), ry.copy(), rz.copy()  # conv may return views of reused buffers

    # Candidate: per-subfault sfuncs (all equal to weight[i]*shared_shape),
    # summed via convolve_dynamic_sources (one conv() call per subfault).
    moment_rate_axitra = weights[:, None] * shared_shape[None, :]
    rake_arr = np.full(nsub, rake)
    _, cx, cy, cz = convolve_dynamic_sources(
        ap, moment_rate_axitra, strike_deg=strike, dip_deg=dip, rake_deg=rake_arr, unit=1
    )

    ref = np.concatenate([rx, ry, rz]).ravel()
    got = np.concatenate([cx, cy, cz]).ravel()
    assert np.abs(ref).max() > 1e-6, "reference synthetics are ~0: the check would be vacuous"
    assert np.allclose(got, ref, rtol=1e-6, atol=1e-4 * np.abs(ref).max())  # FFT round-off ~1e-5

    ap.clean()


@needs_axitra
def test_file_stf_is_moment_function_not_rate():
    """axitra's file source (type 3) is M(t), like its built-ins (type 7 =
    Heaviside, type 4 = integral of a triangle; see axitra/src/fsource.f90).
    Feeding the cumulative triangle must reproduce type 4; feeding the
    triangle itself (the rate) must not. This is what DynamicForwardModel
    relies on when it integrates the fd3d_TSN moment rate.
    """
    from kdellipspy import AxitraForwardModel

    cfg = _mock_config()
    cfg.fault_plane.nx = cfg.fault_plane.ny = 1
    fm = AxitraForwardModel.from_config(cfg, axitra_dir=str(AXITRA_DIR))
    ap = fm.green(fm.build_axitra(fm.build_geometry(), latlon=False), quiet=True)
    from axitra import moment as axitra_moment

    dt = ap.duration / ap.npt
    t = np.arange(ap.npt) * dt
    width = 20 * dt
    hist = np.array([[1, 1e18, 20.0, 60.0, 90.0, 0.0, 0.0, 0.0]])
    _, x4, y4, z4 = axitra_moment.conv(ap, hist, source_type=4, t0=width, unit=1)
    ref = np.concatenate([x4, y4, z4]).ravel().copy()
    assert np.abs(ref).max() > 1e-6

    rate = np.interp(t, [0, width / 2, width], [0, 2 / width, 0], right=0.0)  # area 1
    corr = {}
    for name, sfunc in (("rate", rate), ("cumulative", np.cumsum(rate) * dt)):
        _, x, y, z = axitra_moment.conv(ap, hist, source_type=3, t0=0.0, unit=1, sfunc=sfunc)
        corr[name] = np.corrcoef(np.concatenate([x, y, z]).ravel(), ref)[0, 1]
    ap.clean()

    assert corr["cumulative"] > 0.85, corr
    assert abs(corr["rate"]) < 0.5, corr


def test_resample_to_axitra_grid_can_hold_last_value():
    moment_fn = np.array([[0.0, 1.0, 2.0, 3.0]])
    out = resample_to_axitra_grid(moment_fn, dt_fine_s=1.0, npt_axitra=8, dt_axitra_s=0.5, right=None)
    assert out[0] == pytest.approx([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.0])


@needs_axitra
def test_response_basis_matches_per_subfault_conv():
    """Fast path (precomputed per-subfault axitra operator + numpy) must
    reproduce the reference path (one axitra conv() per subfault) for
    arbitrary per-subfault moment functions and rakes.
    """
    from kdellipspy import AxitraForwardModel

    cfg = _mock_config()
    fm = AxitraForwardModel.from_config(cfg, axitra_dir=str(AXITRA_DIR))
    ap = fm.green(fm.build_axitra(fm.build_geometry(), latlon=False), quiet=True)
    nsub, npt = ap.nsource, ap.npt
    strike, dip = cfg.source_position.strike, cfg.source_position.dip

    # Smooth ramps with different onsets/durations/final moments and rakes.
    t = np.arange(npt)
    rng = np.random.default_rng(0)
    # Values exact in the %.3g the reference path writes to .hist/.sou, so
    # any mismatch is the fast path's fault, not file rounding.
    m_tot = rng.integers(1, 3, nsub) * 1e18
    onset, dur = rng.integers(0, 15, nsub), rng.choice([5, 10, 20, 25], nsub)
    shape = np.clip((t[None, :] - onset[:, None]) / dur[:, None], 0, 1)
    rake = rng.integers(40, 140, nsub).astype(float)
    mx = m_tot[:, None] * shape * np.cos(np.radians(rake))[:, None]
    mz = m_tot[:, None] * shape * np.sin(np.radians(rake))[:, None]

    fast = synthetics_from_basis(response_basis(ap, strike, dip, unit=1), mx, mz, axitra_aw(ap))
    _, sx, sy, sz = convolve_dynamic_sources(
        ap, m_tot[:, None] * shape, strike_deg=strike, dip_deg=dip, rake_deg=rake, unit=1
    )
    slow = np.stack([sx, sy, sz], axis=1)
    ap.clean()

    assert np.abs(slow).max() > 1e-6
    assert np.abs(fast - slow).max() / np.abs(slow).max() < 1e-4


def test_resample_to_axitra_grid_delay_shifts_in_time():
    moment_fn = np.array([[0.0, 1.0, 2.0, 3.0]])
    out = resample_to_axitra_grid(moment_fn, dt_fine_s=1.0, npt_axitra=8, dt_axitra_s=1.0, right=None, delay_s=2.0)
    assert out[0] == pytest.approx([0, 0, 0, 1, 2, 3, 3, 3])
    early = resample_to_axitra_grid(moment_fn, dt_fine_s=1.0, npt_axitra=4, dt_axitra_s=1.0, right=None, delay_s=-1.0)
    assert early[0] == pytest.approx([1, 2, 3, 3])
