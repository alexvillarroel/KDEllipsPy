"""Dynamic-rupture forward model (Fase 4): orchestrates the full pipeline
from a 10-parameter elliptical model to synthetic seismograms —
``params -> forwardmodel.dat -> fd3d_TSN -> slip_rate -> subfault moment-rate
-> axitra synthetics``. Mirrors ``AxitraForwardModel`` (kinematic side), and
composes it directly for geometry/stations/Green's-function plumbing rather
than duplicating it.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ...core.config_parser import ConfigParser
from ...core.forward_model import AxitraForwardModel
from ...core.geometry import TSNFaultGridSpec, tsn_hypocentre_coarse
from .dynamic_convolution import (
    bin_slip_rate_to_subfaults,
    convolve_dynamic_sources,
    project_rake,
    resample_to_axitra_grid,
    axitra_aw,
    response_basis,
    synthetics_from_basis,
)
from .tsn_bridge import TSNRunConfig, run_tsn_forward


class DynamicForwardModel:
    """One forward evaluation: 10-param dynamic model -> synthetics (nsta, 3, npts).

    The subfault mesh (``cfg.fault_plane.nx`` x ``ny``) is used BOTH as the
    coarse grid fd3d_TSN's ``forwardmodel.dat`` is generated on (``nli=nx``,
    ``nwi=ny``) and as the axitra source grid the fine slip-rate is binned
    down to — one shared resolution knob instead of two independent ones.

    Green's functions are computed once (full fixed mesh, all subfaults
    always present — positions never change across models, only the dynamic
    slip does) and cached for reuse across evaluations.
    """

    def __init__(
        self,
        cfg: ConfigParser,
        tsn_run_cfg: TSNRunConfig,
        tsn_grid: TSNFaultGridSpec,
        axitra_dir: Optional[str] = None,
        axitra_aw: float = 0.5,
        axitra_ikmax: int = 100000,
    ):
        self.cfg = cfg
        self.fm = AxitraForwardModel.from_config(cfg, axitra_dir=axitra_dir)
        self.tsn_run_cfg = tsn_run_cfg
        self.tsn_grid = tsn_grid
        self.axitra_aw = float(axitra_aw)
        self.axitra_ikmax = int(axitra_ikmax)

        self.nx = int(cfg.fault_plane.nx)
        self.ny = int(cfg.fault_plane.ny)
        # Nucleation at the input.ctl hypocentre (fd3d coarse coords).
        self.hypo_coarse = tsn_hypocentre_coarse(cfg.fault_plane)
        # Fixed mesh, all nx*ny subfaults, positions only (no slip applied).
        self.base_geometry = self.fm.build_geometry()
        self._mu_pa = float(np.mean([sf.mu_pa for sf in self.base_geometry.subfaults]))

        self._ap = None
        self._basis = {}  # unit -> step-response basis
        # False -> legacy path: one axitra conv() per subfault per evaluation
        # (~200 s at 25x25); kept as the reference the fast path is tested against.
        self.use_basis = True
        # Origin-time correction (s, >0 = later) applied to the moment
        # functions: absorbs hypocentre/origin-time inconsistencies in the input.
        self.time_shift_s = 0.0
        self._last_run = (None, None)  # (model bytes, slip-rates): skip fd3d when only the shift changes
        # Imposed scalar moment (N·m), analogue of the kinematic mt_strict. Exact
        # thanks to slip-weakening self-similarity (stresses & Dc x k -> slip x k,
        # same timing): synthetics are rescaled by k = m0_target / M0 and the
        # model is equivalent to Te*k, Dc*k. ``last_m0_scale`` keeps that k.
        self.m0_target: Optional[float] = None
        self.last_m0_scale = 1.0

    def ensure_basis(self, unit: int) -> np.ndarray:
        """Compute (once per output unit) the per-subfault axitra operator."""
        if unit not in self._basis:
            ap = self.ensure_green()
            strike = float(self.cfg.source_position.strike)
            dip = float(self.cfg.source_position.dip)
            self._basis[unit] = response_basis(ap, strike, dip, unit=unit)
        return self._basis[unit]

    def ensure_green(self):
        """Compute (once) and cache Green's functions for the fixed mesh."""
        if self._ap is None:
            ap = self.fm.build_axitra(
                self.base_geometry,
                latlon=False,
                freesurface=True,
                aw=self.axitra_aw,
                ikmax=self.axitra_ikmax,
            )
            self._ap = self.fm.green(ap, quiet=True)
        return self._ap

    def forward(self, model: np.ndarray, unit: Optional[int] = None) -> np.ndarray:
        """Run one dynamic-rupture forward evaluation.

        Returns synthetics of shape ``(nsta, 3, npts)`` (X, Y, Z per axitra's
        ``moment.conv`` convention), sampled on axitra's own time grid
        (``ap.npt`` samples at ``ap.duration/ap.npt``) — the caller is
        responsible for matching this against ``cfg.observed_data`` sampling
        if they differ, same as the kinematic path.
        """
        ap = self.ensure_green()

        key = np.asarray(model, dtype=np.float32).tobytes()
        if self._last_run[0] == key:
            result = self._last_run[1]
        else:
            result = run_tsn_forward(model, self.nx, self.ny, self.tsn_grid, self.tsn_run_cfg, hypo=self.hypo_coarse)
            self._last_run = (key, result)
        slip_x, slip_z = result["sliprateX"], result["sliprateZ"]

        mrate_x, mrate_z = bin_slip_rate_to_subfaults(
            slip_x, slip_z, self.nx, self.ny, self.tsn_grid.dh, self._mu_pa
        )
        unit_val = int(self.cfg.observed_data.units) if unit is None else int(unit)
        dt_fd = self.tsn_run_cfg.dt_s
        dt_axitra = ap.duration / ap.npt

        # axitra's file STF is M(t), not dM/dt: integrate on the fine FD grid
        # (keeps the moment of sharp slip-rate peaks), then resample holding
        # the final moment past the FD run.
        def to_axitra_moment(rate):
            return resample_to_axitra_grid(
                np.cumsum(rate, axis=1) * dt_fd, dt_fd, ap.npt, dt_axitra, right=None, delay_s=self.time_shift_s
            )

        self.last_m0_scale = 1.0
        if self.m0_target is not None:
            m0 = float(np.hypot(mrate_x.sum(axis=1), mrate_z.sum(axis=1)).sum() * dt_fd)
            self.last_m0_scale = self.m0_target / m0 if m0 > 0 else 0.0
            mrate_x, mrate_z = mrate_x * self.last_m0_scale, mrate_z * self.last_m0_scale

        if self.use_basis:
            # Local X -> rake 0, local Z -> rake 90 (same convention as project_rake).
            basis = self.ensure_basis(unit_val)
            return synthetics_from_basis(basis, to_axitra_moment(mrate_x), to_axitra_moment(mrate_z), axitra_aw(ap))

        moment_rate, rake_deg = project_rake(mrate_x, mrate_z)
        strike = float(self.cfg.source_position.strike)
        dip = float(self.cfg.source_position.dip)
        _, sx, sy, sz = convolve_dynamic_sources(
            ap, to_axitra_moment(moment_rate), strike_deg=strike, dip_deg=dip, rake_deg=rake_deg, unit=unit_val
        )
        return np.transpose(np.array([sx, sy, sz]), (1, 0, 2))

    def clean(self) -> None:
        """Release the cached Green's-function axitra files."""
        if self._ap is not None:
            try:
                self._ap.clean()
            except Exception:
                pass
            self._ap = None
        self._basis = {}
