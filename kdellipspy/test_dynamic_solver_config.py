"""Fase 5: Section 10 (TSN dynamic solver config) parsing. Fast, no fd3d_TSN
or axitra binaries needed -- pure ConfigParser wiring check.
"""

from pathlib import Path

from kdellipspy import ConfigParser

MINIMAL_CTL = """\
# 1. Observed Data Parameters
Time window start (t1) : 0.0
Time window end (t2) : 10.0
Number of points (Npts) : 100
Delta / Time step : 0.1
Units : 1

# 10. TSN Dynamic Solver
Solver work directory : /tmp/tsn_run
FD grid points along strike (NXT) : 40
FD grid points along dip (NZT) : 20
FD grid spacing (dh) : 100.0
Absorbing boundary cells (nabc) : 5
FD time step (dt) : 0.003
Solver binary name : fd3d_gnu_TSN
"""


def test_section_10_parsed_from_file(tmp_path: Path):
    ctl_path = tmp_path / "input.ctl"
    ctl_path.write_text(MINIMAL_CTL)

    cfg = ConfigParser(filepath=str(ctl_path))

    assert cfg.dynamic_solver is not None
    assert cfg.dynamic_solver.work_dir == "/tmp/tsn_run"
    assert cfg.dynamic_solver.nxtT == 40
    assert cfg.dynamic_solver.nztT == 20
    assert cfg.dynamic_solver.dh == 100.0
    assert cfg.dynamic_solver.nabc == 5
    assert cfg.dynamic_solver.dt_s == 0.003
    assert cfg.dynamic_solver.binary == "fd3d_gnu_TSN"


def test_section_10_absent_by_default():
    cfg = ConfigParser.from_dict({})
    assert cfg.dynamic_solver is None


def test_section_10_from_dict():
    cfg = ConfigParser.from_dict({
        "dynamic_solver": {
            "Solver work directory": "/tmp/run2",
            "FD grid points along strike (NXT)": 8,
            "FD grid points along dip (NZT)": 4,
        }
    })
    assert cfg.dynamic_solver.work_dir == "/tmp/run2"
    assert cfg.dynamic_solver.nxtT == 8
    assert cfg.dynamic_solver.nztT == 4
    assert cfg.dynamic_solver.binary == "fd3d_gnu_TSN"  # default
