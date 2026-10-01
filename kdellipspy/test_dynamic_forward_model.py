"""End-to-end validation of Fase 4: DynamicForwardModel and
DynamicNAInversionModel — the full params -> fd3d_TSN -> slip_rate ->
axitra synthetics -> misfit pipeline, against the real compiled fd3d_TSN and
axitra binaries.

Does NOT run a full multi-iteration NA search (each dynamic forward
evaluation runs the real FD solver — too slow for a unit test); it validates
a single forward/objective_function call, which is what proves the wiring
into BaseInversionModel's misfit/NAResult/checkpointing machinery works.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from kdellipspy import ConfigParser, InversionParam
from kdellipspy.core.forward_model import AxitraForwardModel
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicForwardModel, DynamicNAInversionModel, TSNRunConfig
from kdellipspy.inversion.kinematic.model_na import NAConfig

FD3D_TSN_SRC = Path("/home/alex/fd3d_TSN/src")
AMATRICE_EXAMPLE = Path("/home/alex/fd3d_TSN/examples/Amatrice2016")
AXITRA_DIR = Path(__file__).resolve().parent / "axitra" / "src"

needs_both_solvers = pytest.mark.skipif(
    not FD3D_TSN_SRC.exists() or shutil.which("gfortran") is None or not (AXITRA_DIR / "axitra").exists(),
    reason="fd3d_TSN clone, gfortran or compiled axitra binary not available",
)

# Small, fast grid: 4x2 subfaults, FD grid divisible by that (40x20).
NX_SUB, NY_SUB = 4, 2
NXTT, NZTT = 40, 20
DH, DIP_DEG, NABC = 100.0, 45.0, 5
NT_FD, DT_FD = 20, 0.003


def _base_params():
    return {
        "source_position": {
            "event_name": "DynTest", "latitude": -33.45, "longitude": -70.66,
            "depth": 10.0, "strike": 20.0, "dip": DIP_DEG, "rake": 90.0,
        },
        "fault_plane": {
            "Length along strike (Lx)": NXTT * DH,
            "Length along dip (Ly)": NZTT * DH,
            "Number of subfaults along strike (Nx)": NX_SUB,
            "Number of subfaults along dip (Ny)": NY_SUB,
            "Hypocenter position strike (Hx)": NXTT * DH / 2,
            "Hypocenter position dip (Hy)": NZTT * DH / 2,
        },
        "ellipse": {
            "Number of ellipses": 1, "Initial slip": 0, "Slip shape": 1,
            "Frequency 1 (Freq1)": 0.05, "Frequency 2 (Freq2)": 2.0, "Time shift (T0)": 1.0,
        },
        "observed_data": {
            # Long enough window that the Butterworth filter (order 2*n_int,
            # zero-phase) has a valid padlen for filtfilt -- short windows
            # (<~150 samples here) silently produce NaN.
            "Time window start (t1)": -5.0, "Time window end (t2)": 30.0,
            "Number of points (Npts)": 512, "Delta / Time step": 0.1, "Units": 1,
        },
        # SI units, "thickness" = layer-top depth (m), same as input.ctl.
        "velocity_model": [
            {"thickness": 0.0, "vp": 5000.0, "vs": 3000.0, "rho": 2500.0, "qp": 500.0, "qs": 200.0},
            {"thickness": 5000.0, "vp": 6500.0, "vs": 3800.0, "rho": 2800.0, "qp": 1000.0, "qs": 500.0},
        ],
        "moment_tensor": {"Moment Tensor Flag": 0},
        "stations": [
            {"name": "STA1", "latitude": -33.40, "longitude": -70.66, "height": 0.0, "use_n": True, "use_e": True, "use_z": True},
            {"name": "STA2", "latitude": -33.45, "longitude": -70.60, "height": 0.0, "use_n": True, "use_e": True, "use_z": True},
        ],
        # Ranges bracketing the values used elsewhere in this file (model =
        # [8, 4, 12, 4, 0.3, 5, 1.05, 1.15, 2, 0.4]) -- just wide enough for
        # a tiny NA search to have somewhere to sample.
        "inversion_params": [
            InversionParam(name=n, min_val=lo, max_val=hi, flag=1)
            for n, lo, hi in zip(
                DynamicNAInversionModel._DYNAMIC_PARAM_NAMES,
                [6.0, 3.0, 9.0, 3.0, 0.1, 4.0, 0.9, 1.0, 1.5, 0.3],
                [10.0, 5.0, 15.0, 5.0, 0.5, 6.0, 1.2, 1.3, 2.5, 0.5],
            )
        ],
    }


@pytest.fixture(scope="module")
def solver_work_dir(tmp_path_factory):
    work_dir = tmp_path_factory.mktemp("dynamic_forward_model_test")
    for f in ["fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90", "inversion_com.f90"]:
        shutil.copy(FD3D_TSN_SRC / f, work_dir / f)
    shutil.copy(AMATRICE_EXAMPLE / "crustal.dat", work_dir / "crustal.dat")

    binary = work_dir / "fd3d_gnu_TSN"
    subprocess.run(
        ["gfortran", "-DDIPSLIP", "-O2", "-cpp", "-o", str(binary),
         "fd3d_init.f90", "fd3d_deriv.f90", "fd3d_theo.f90", "dynamicsolver.f90"],
        cwd=work_dir, check=True, capture_output=True, text=True,
    )
    (work_dir / "inputinv.dat").write_text(f"0\n{NX_SUB} {NY_SUB}\n")
    (work_dir / "inputfd3d.dat").write_text(
        f"{NXTT} 10 {NZTT}       \\* Number of grid points in x,y,z direction (FD) *\\\n"
        f"{DH}             \\* Grid size in meters (FD) *\\\n"
        f"{NT_FD}                \\* Number of time step (FD) *\\\n"
        f"{DT_FD}             \\* Time step in seconds (FD) *\\\n"
        f"{DIP_DEG}               \\* Dip *\\\n"
        f"{NABC} 6000. 25.      \\* nabc, vp, pml_fact *\\\n"
        "0.3               \\* damp_s *\\\n"
        "0\t\t  \\* number of seismic stations\n"
        "0\n"
    )
    return work_dir


def _build_config_and_dynfm(solver_work_dir):
    cfg = ConfigParser.from_dict(_base_params())
    fm_probe = AxitraForwardModel.from_config(cfg, axitra_dir=str(AXITRA_DIR))
    npt = fm_probe._estimate_npt()
    cfg.observed_data.npts = npt  # match axitra's own grid so misfit/filter shapes line up

    grid = TSNFaultGridSpec(dh=DH, dip_deg=DIP_DEG, nztT=NZTT, nabc=NABC)
    run_cfg = TSNRunConfig(work_dir=solver_work_dir, nxtT=NXTT, nztT=NZTT, dt_s=DT_FD)
    return cfg, grid, run_cfg, npt


@needs_both_solvers
def test_dynamic_forward_model_produces_finite_synthetics(solver_work_dir):
    cfg, grid, run_cfg, npt = _build_config_and_dynfm(solver_work_dir)
    dyn_fm = DynamicForwardModel(cfg, run_cfg, grid, axitra_dir=str(AXITRA_DIR))

    model = np.array([8.0, 4.0, 12.0, 4.0, 0.3, 5.0, 1.05, 1.15, 2.0, 0.4], dtype=np.float32)
    synthetics = dyn_fm.forward(model)

    nsta = len(cfg.stations.stations)
    assert synthetics.shape == (nsta, 3, npt)
    assert np.isfinite(synthetics).all()

    dyn_fm.clean()


@needs_both_solvers
def test_dynamic_na_inversion_model_objective_function(solver_work_dir):
    cfg, grid, run_cfg, npt = _build_config_and_dynfm(solver_work_dir)

    nsta = len(cfg.stations.stations)
    rng = np.random.default_rng(0)
    observed = rng.normal(scale=1e-6, size=(nsta, 3, npt))
    dt = float(cfg.observed_data.delta)
    time_array = np.arange(npt) * dt + float(cfg.observed_data.t1)

    inv_model = DynamicNAInversionModel(
        config=cfg,
        axitra_dir=str(AXITRA_DIR),
        observed_waveforms=observed,
        time_array=time_array,
        tsn_run_cfg=run_cfg,
        tsn_grid=grid,
    )
    assert inv_model.misfit_calc is not None, "azi_times auto-build failed; misfit_calc is None"

    model = np.array([8.0, 4.0, 12.0, 4.0, 0.3, 5.0, 1.05, 1.15, 2.0, 0.4], dtype=np.float32)
    misfit = inv_model.objective_function(model)

    assert np.isfinite(misfit)
    assert misfit != 1e10  # sentinel for a failed evaluation
    assert inv_model.best_synthetics is not None
    assert inv_model.best_synthetics.shape == (nsta, 3, npt)

    inv_model.clean()


@needs_both_solvers
def test_dynamic_na_search_runs_tiny_multi_iteration(solver_work_dir):
    """The one thing Fase 4/5 deliberately skipped: an actual multi-iteration
    run_na_search() loop (not just a single objective_function() call), on
    the smallest config that still exercises resampling (ni=2, ns=2, nr=1,
    n=1 -> 4 real fd3d_TSN forward evaluations total).
    """
    cfg, grid, run_cfg, npt = _build_config_and_dynfm(solver_work_dir)

    nsta = len(cfg.stations.stations)
    rng = np.random.default_rng(0)
    observed = rng.normal(scale=1e-6, size=(nsta, 3, npt))
    dt = float(cfg.observed_data.delta)
    time_array = np.arange(npt) * dt + float(cfg.observed_data.t1)

    inv_model = DynamicNAInversionModel(
        config=cfg,
        axitra_dir=str(AXITRA_DIR),
        observed_waveforms=observed,
        time_array=time_array,
        tsn_run_cfg=run_cfg,
        tsn_grid=grid,
    )

    result = inv_model.run_na_search(
        NAConfig(n_samples_initial=2, n_samples_iteration=2, n_iterations=1,
                 n_cells_resample=1, n_jobs=1, random_seed=0)
    )

    assert len(result.param_names) == 10
    assert np.isfinite(result.best_model.misfit)
    assert result.best_model.misfit != 1e10

    inv_model.clean()


@needs_both_solvers
def test_dynamic_forward_passes_cumulative_moment_to_axitra(monkeypatch, tmp_path):
    """fd3d_TSN gives slip RATE; axitra's file STF is M(t). With a known
    constant slip-rate on every fault cell for the first half of the run,
    the per-subfault function handed to axitra must ramp up and then hold
    at mu * slip * area (not drop back to zero like a rate would).
    """
    import kdellipspy.inversion.dynamic.forward_model_dynamic as fmd

    cfg = ConfigParser.from_dict(_base_params())
    grid = TSNFaultGridSpec(dh=DH, dip_deg=DIP_DEG, nztT=NZTT, nabc=NABC)
    dt, nt = 0.05, 100
    run_cfg = TSNRunConfig(work_dir=tmp_path, nxtT=NXTT, nztT=NZTT, dt_s=dt)
    dyn_fm = DynamicForwardModel(cfg, run_cfg, grid, axitra_dir=str(AXITRA_DIR))

    rate = np.zeros((nt, NXTT, NZTT), dtype=np.float32)
    rate[: nt // 2] = 1.0  # 1 m/s for 2.5 s -> 2.5 m final slip everywhere
    monkeypatch.setattr(fmd, "run_tsn_forward", lambda *a, **k: {"sliprateX": 0 * rate, "sliprateZ": rate})
    captured = {}

    def fake_synthetics(basis, moment_x, moment_z, aw):
        captured["sfunc"] = moment_z
        return np.zeros((len(cfg.stations.stations), 3, moment_z.shape[1]))

    monkeypatch.setattr(fmd, "synthetics_from_basis", fake_synthetics)
    monkeypatch.setattr(dyn_fm, "ensure_basis", lambda unit: None)
    dyn_fm.forward(np.zeros(10, dtype=np.float32))

    sfunc = captured["sfunc"]
    cells_per_sub = (NXTT // NX_SUB) * (NZTT // NY_SUB)
    m0_sub = dyn_fm._mu_pa * 2.5 * cells_per_sub * DH**2
    assert np.all(np.diff(sfunc, axis=1) >= -1e-6 * m0_sub), "not monotonic: looks like a rate"
    assert sfunc[:, -1] == pytest.approx(m0_sub, rel=0.02)  # holds final moment past the FD run
    dyn_fm.clean()


@needs_both_solvers
def test_dynamic_forward_m0_target_rescales_moment(monkeypatch, tmp_path):
    """With m0_target set, the moment handed to axitra totals m0_target."""
    import kdellipspy.inversion.dynamic.forward_model_dynamic as fmd

    cfg = ConfigParser.from_dict(_base_params())
    grid = TSNFaultGridSpec(dh=DH, dip_deg=DIP_DEG, nztT=NZTT, nabc=NABC)
    dt, nt = 0.05, 100
    dyn_fm = DynamicForwardModel(cfg, TSNRunConfig(work_dir=tmp_path, nxtT=NXTT, nztT=NZTT, dt_s=dt), grid,
                                 axitra_dir=str(AXITRA_DIR))
    rate = np.zeros((nt, NXTT, NZTT), dtype=np.float32)
    rate[: nt // 2] = 1.0
    monkeypatch.setattr(fmd, "run_tsn_forward", lambda *a, **k: {"sliprateX": 0 * rate, "sliprateZ": rate})
    captured = {}
    monkeypatch.setattr(fmd, "synthetics_from_basis",
                        lambda basis, mx, mz, aw: captured.setdefault("mz", mz) * 0 + np.zeros(1))
    monkeypatch.setattr(dyn_fm, "ensure_basis", lambda unit: None)

    dyn_fm.m0_target = 1.7e19
    dyn_fm.forward(np.zeros(10, dtype=np.float32))
    m0_free = dyn_fm._mu_pa * 2.5 * NXTT * NZTT * DH**2
    assert dyn_fm.last_m0_scale == pytest.approx(1.7e19 / m0_free, rel=1e-4)
    assert captured["mz"][:, -1].sum() == pytest.approx(1.7e19, rel=0.02)
    dyn_fm.clean()
