"""Subprocess bridge to a compiled fd3d_TSN binary (Plan B of the dynamic
integration plan): write ``forwardmodel.dat`` from the 10-param ellipse
model, run the solver, read the on-fault slip-rate back as NumPy arrays.

This is the deliberately "ugly" MVP path — file I/O per evaluation — used to
validate the physics end-to-end before investing in an f2py wrap of
``dynamicsolver.f90`` (a PROGRAM, not a SUBROUTINE: it reads
forwardmodel.dat/inputfd3d.dat/inputinv.dat from disk, so wrapping it
directly needs a solver-side I/O refactor first).
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np

from ...core.geometry import TSNFaultGridSpec, build_tsn_dynamic_fields, write_fd3d_tsn_forwardmodel


def read_tsn_fault_field(path: str | Path, nxt: int, nzt: int) -> np.ndarray:
    """Read a ``result/{sliprateX,sliprateZ,shearstressX,shearstressZ}.res``
    stream file into an array of shape ``(nt, nxt, nzt)``.

    fd3d_TSN writes one unformatted STREAM record per FD timestep, each a
    Fortran ``(nxt, nzt)`` slice (column-major, i.e. x varies fastest) —
    see ``fd3d_theo.f90:546-553``. No header, so ``nt`` is inferred from the
    file size.
    """
    path = Path(path)
    raw = np.fromfile(path, dtype=np.float32)
    per_record = nxt * nzt
    if raw.size % per_record != 0:
        raise ValueError(
            f"{path}: {raw.size} floats is not a multiple of nxt*nzt={per_record}"
        )
    nt = raw.size // per_record
    # Fortran column-major (x fastest) per record == C-order reshape (nzt, nxt).
    return raw.reshape(nt, nzt, nxt).transpose(0, 2, 1)  # -> (nt, nxt, nzt)


@dataclass
class TSNRunConfig:
    """Where to run the solver and how big its on-fault output grid is.
    (nxtT/nztT are the *fine* FD on-fault dimensions from ``inputfd3d.dat``
    — needed to parse the ``.res`` outputs, independent of the coarse
    ``nli``/``nwi`` the ellipse fields are generated on.)
    """

    work_dir: Path
    nxtT: int
    nztT: int
    dt_s: float
    binary: str = "fd3d_gnu_TSN"

    def binary_path(self) -> Path:
        b = Path(self.binary)
        return b if b.is_absolute() else self.work_dir / b


def run_tsn_forward(
    model: np.ndarray,
    nli: int,
    nwi: int,
    grid: TSNFaultGridSpec,
    run_cfg: TSNRunConfig,
    timeout_s: float = 3600.0,
) -> Dict[str, np.ndarray]:
    """Run one fd3d_TSN forward evaluation for a 10-param ellipse model.

    ``run_cfg.work_dir`` must already contain the solver binary plus
    ``inputfd3d.dat``, ``inputinv.dat`` and ``crustal.dat`` (the grid/medium
    setup is fixed across NA evaluations — only ``forwardmodel.dat``
    changes per model). Returns ``{"sliprateX", "sliprateZ"}``, each
    ``(nt, nxtT, nztT)`` float32.
    """
    work_dir = Path(run_cfg.work_dir)
    t0, ts, dc = build_tsn_dynamic_fields(model, nli, nwi, grid)
    write_fd3d_tsn_forwardmodel(work_dir / "forwardmodel.dat", t0, ts, dc)

    (work_dir / "result").mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [str(run_cfg.binary_path())],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"fd3d_TSN forward failed (exit {proc.returncode}):\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )

    result_dir = work_dir / "result"
    slip_x = read_tsn_fault_field(result_dir / "sliprateX.res", run_cfg.nxtT, run_cfg.nztT)
    slip_z = read_tsn_fault_field(result_dir / "sliprateZ.res", run_cfg.nxtT, run_cfg.nztT)
    return {"sliprateX": slip_x, "sliprateZ": slip_z}
