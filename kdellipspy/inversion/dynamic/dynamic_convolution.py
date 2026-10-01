"""Bridge from fd3d_TSN's spatially-distributed slip-rate field to axitra
synthetics (Fase 3 of the dynamic integration plan).

Why this can't just call ``moment.conv()`` once like the kinematic path:
axitra's ``convmPy.moment_conv`` applies ONE shared source-time-function
shape (``sfunc``/``source_type``) to every row of ``hist`` — each source only
gets its own scalar moment and delay, not its own distinct time history
(verified by reading ``axitra/src/convm.f90``: ``freqs=fsou(jf)`` is computed
once per frequency, outside the source loop). The dynamic slip-rate field has
a genuinely different moment-rate *shape* at every fault point (different
onset, duration and pulse shape, set by the physics of the spontaneous
rupture) — so it needs one ``conv()`` call per subfault, each activating only
that subfault (``disp=1``, everyone else ``disp=0``) with its own ``sfunc``,
summed by superposition (valid because the Green's functions in ``ap`` are
shared and convolution/summation is linear).

Moment convention (verified in ``convm.f90:cmoment``: ``if surf==0:
xmoment=disp``, matching ``FaultGeometry.to_axitra_hist()``): hist's moment
column is used directly as the scalar seismic moment. Here we set it to 1.0
and put the full physical moment-rate (Nm/s) into ``sfunc`` instead, since
that avoids a second, redundant amplitude scale living in two places.

ponytail: rake is FIXED per subfault (time-averaged slip direction), not
time-varying — axitra's per-source mechanism (strike/dip/rake) is static for
a given ``conv()`` call. A rotating rake during rupture would need the full
6-elementary-DC moment-tensor decomposition already used for MT sources on
the kinematic side (``GeometryBuilder._mt_basis_and_amplitudes``). Upgrade
path: reuse that decomposition here if a benchmark shows the fixed-rake
approximation matters.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np


def bin_slip_rate_to_subfaults(
    slip_x: np.ndarray,
    slip_z: np.ndarray,
    nx_sub: int,
    nz_sub: int,
    dh_fine_m: float,
    mu_pa: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Aggregate the fine FD slip-rate field into subfault moment-rate
    histories, by summing slip-rate over the fine cells in each subfault
    patch (moment is additive: ``dM/dt = mu * sum(dA_fine * slip_rate)``).

    Parameters
    ----------
    slip_x, slip_z : (nt, nxtT, nztT) float arrays, fd3d_TSN's local
        Cartesian slip-rate components (from
        :func:`kdellipspy.inversion.dynamic.tsn_bridge.read_tsn_fault_field`).
    nx_sub, nz_sub : subfault discretization (must evenly divide nxtT, nztT).
    dh_fine_m : fine FD grid spacing (metres).
    mu_pa : shear modulus at the fault depth (Pa), assumed uniform.

    Returns
    -------
    moment_rate_x, moment_rate_z : (nsub, nt) float64, Nm/s, subfault index
        ordered ``(idip-1)*nx_sub + istk`` (1-based, strike-fastest), matching
        :meth:`EllipticalSlipMapper._subfault_fault_plane_xy` / the standard
        ``GeometryBuilder`` subfault ordering.
    """
    nt, nxt, nzt = slip_x.shape
    if slip_z.shape != slip_x.shape:
        raise ValueError(f"slip_x/slip_z shape mismatch: {slip_x.shape} vs {slip_z.shape}")
    if nxt % nx_sub != 0 or nzt % nz_sub != 0:
        raise ValueError(
            f"fine grid ({nxt},{nzt}) must divide evenly by subfault grid ({nx_sub},{nz_sub})"
        )
    fx, fz = nxt // nx_sub, nzt // nz_sub
    dA = float(dh_fine_m) ** 2
    scale = float(mu_pa) * dA

    def _bin(field: np.ndarray) -> np.ndarray:
        # (nt, nx_sub, fx, nz_sub, fz) -> sum fine cells -> (nt, nx_sub, nz_sub)
        reshaped = field.reshape(nt, nx_sub, fx, nz_sub, fz)
        return reshaped.sum(axis=(2, 4)) * scale

    # fd3d_TSN writes dip rows from the DEEP edge up (k=nabc+1 .. nzt-nfs), the
    # axitra mesh numbers them from the SHALLOW edge (idip=1 on top): flip the
    # dip axis so each fd3d row lands on the subfault at its own depth.
    mrate_x = _bin(slip_x)[:, :, ::-1]  # (nt, nx_sub, nz_sub), nz_sub 0 = shallow
    mrate_z = _bin(slip_z)[:, :, ::-1]

    nsub = nx_sub * nz_sub
    # Flatten (nx_sub, nz_sub) -> subfault index, strike-fastest (istk varies first).
    mrate_x = mrate_x.transpose(0, 2, 1).reshape(nt, nsub).T  # -> (nsub, nt)
    mrate_z = mrate_z.transpose(0, 2, 1).reshape(nt, nsub).T
    return np.ascontiguousarray(mrate_x, dtype=np.float64), np.ascontiguousarray(mrate_z, dtype=np.float64)


def project_rake(
    moment_rate_x: np.ndarray, moment_rate_z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Collapse the (X, Z) moment-rate pair into one scalar moment-rate
    history per subfault plus a single representative rake angle.

    The rake is taken from the *time-integrated* slip direction
    (``atan2(sum(Mz), sum(Mx))``) — the direction of net slip, which is the
    standard notion of "average rake" and more robust than an instantaneous
    mean (early near-zero samples wouldn't dominate it). The scalar
    moment-rate is the projection of the (Mx, Mz) vector onto that fixed
    direction, chosen so its time integral equals the magnitude of the net
    moment vector (see docstring derivation in the module tests).

    Returns
    -------
    moment_rate : (nsub, nt) float64, Nm/s
    rake_deg    : (nsub,) float64, degrees (fd3d_TSN local X/Z convention:
        0 deg = pure local-X slip, 90 deg = pure local-Z slip — the caller
        is responsible for mapping this into axitra's strike/dip/rake frame
        if fd3d_TSN's local X/Z axes are not already strike/dip-aligned).
    """
    mx_tot = moment_rate_x.sum(axis=1)
    mz_tot = moment_rate_z.sum(axis=1)
    rake_rad = np.arctan2(mz_tot, mx_tot)
    moment_rate = moment_rate_x * np.cos(rake_rad)[:, None] + moment_rate_z * np.sin(rake_rad)[:, None]
    return moment_rate, np.degrees(rake_rad)


def resample_to_axitra_grid(
    moment_rate: np.ndarray, dt_fine_s: float, npt_axitra: int, dt_axitra_s: float,
    right: float | None = 0.0, delay_s: float = 0.0,
) -> np.ndarray:
    """Linearly resample ``(nsub, nt_fine)`` histories (sampled at the FD
    timestep) onto axitra's fixed ``npt_axitra``-sample grid
    (``dt_axitra_s = ap.duration / ap.npt``). Past the FD run's duration it
    pads with ``right`` (0 for a rate; ``None`` holds the last value, which
    is what a cumulative moment function needs). ``delay_s`` shifts the
    history in time (>0 later); a negative delay cuts the onset.
    """
    nsub, nt_fine = moment_rate.shape
    t_fine = np.arange(nt_fine) * dt_fine_s + delay_s
    t_axitra = np.arange(npt_axitra) * dt_axitra_s
    out = np.zeros((nsub, npt_axitra), dtype=np.float64)
    for i in range(nsub):
        out[i] = np.interp(t_axitra, t_fine, moment_rate[i], left=0.0, right=right)
    return out


def convolve_dynamic_sources(
    ap,
    moment_rate_axitra: np.ndarray,
    strike_deg,
    dip_deg,
    rake_deg: np.ndarray,
    unit: int = 1,
    delay_s: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sum one ``moment.conv()`` call per subfault (each activating only
    that subfault, ``disp=1.0``, ``sfunc=moment_rate_axitra[i]``) into the
    combined synthetics — the superposition described in the module
    docstring. ``ap`` must already have Green's functions computed
    (``moment.green(ap)``) for exactly ``ap.nsource == moment_rate_axitra.shape[0]``
    sources, in the same order as ``moment_rate_axitra``'s rows.

    NOTE: axitra's file source (type 3) is the source FUNCTION M(t) (step-like:
    type 7 Heaviside = 1/(i*omega), type 4 = integral of a triangle — see
    axitra/src/fsource.f90), NOT its rate. Pass the cumulative moment here.

    ``strike_deg``/``dip_deg`` may be a scalar (shared by all subfaults) or
    an ``(nsub,)`` array.
    """
    try:
        from axitra import moment
    except ImportError:
        import sys

        sys.path.append(str(ap.axpath))
        from axitra import moment

    nsub, npt = moment_rate_axitra.shape
    if nsub != ap.nsource:
        raise ValueError(f"moment_rate_axitra has {nsub} subfaults, ap.nsource={ap.nsource}")

    strike_arr = np.broadcast_to(np.asarray(strike_deg, dtype=float), (nsub,))
    dip_arr = np.broadcast_to(np.asarray(dip_deg, dtype=float), (nsub,))

    sismox = sismoy = sismoz = None
    for i in range(nsub):
        hist = np.zeros((nsub, 8), dtype=np.float64)
        hist[:, 0] = np.arange(1, nsub + 1)
        hist[i, 1] = 1.0
        hist[i, 2] = strike_arr[i]
        hist[i, 3] = dip_arr[i]
        hist[i, 4] = rake_deg[i]
        hist[i, 7] = delay_s

        result = moment.conv(ap, hist, source_type=3, t0=0.0, unit=unit, sfunc=moment_rate_axitra[i])
        if result is None:
            raise RuntimeError(f"moment.conv failed for subfault {i + 1}")
        time, sx, sy, sz = result
        if sismox is None:
            sismox, sismoy, sismoz = sx.copy(), sy.copy(), sz.copy()
        else:
            sismox += sx
            sismoy += sy
            sismoz += sz

    return time, sismox, sismoy, sismoz


def response_basis(ap, strike_deg: float, dip_deg: float, unit: int = 1) -> np.ndarray:
    """Per-subfault axitra operator for rake 0 and rake 90, as damped spectra
    ``(2, nsub, nsta, 3, npt//2+1)`` ready for :func:`synthetics_from_basis`.

    axitra's conv (convm.f90:243,256-258,340-341) computes
    ``u = circ(g, d*s) / d`` with ``d[n] = exp(-pi*aw*n/N)``: linear in the
    source, but CIRCULAR and only weakly damped (aw=0.5 -> d[N]=0.21), so it
    is not shift-invariant and a step/delay basis is wrong at the few-% level.
    The response to an impulse s=delta[0] is ``h = g`` (d[0]=1), so storing
    ``FFT(d*h)`` reproduces axitra exactly for any source. NOTE: the damping
    is the ``aw`` in axitra's ``<sid>.data`` namelist (axitra.py hard-codes
    ``aw=2.``, ignoring ``ap.aw``), read by :func:`axitra_aw`. A double couple's
    moment tensor is linear in the slip vector (rake r -> cos r*M(0) +
    sin r*M(90)), hence two rakes suffice. Costs 2*nsub conv calls once,
    replacing nsub conv calls PER forward evaluation.
    """
    try:
        from axitra import moment
    except ImportError:
        import sys

        sys.path.append(str(ap.axpath))
        from axitra import moment

    nsub, npt = ap.nsource, ap.npt
    impulse = np.zeros(npt)
    impulse[0] = 1.0
    damp = _axitra_damping(axitra_aw(ap), npt)
    basis = None
    for r_idx, rake in enumerate((0.0, 90.0)):
        for i in range(nsub):
            hist = np.zeros((nsub, 8), dtype=np.float64)
            hist[:, 0] = np.arange(1, nsub + 1)
            hist[i, 1:5] = [1.0, strike_deg, dip_deg, rake]
            result = moment.conv(ap, hist, source_type=3, t0=0.0, unit=unit, sfunc=impulse)
            if result is None:
                raise RuntimeError(f"moment.conv failed for subfault {i + 1}, rake {rake}")
            _, sx, sy, sz = result
            h = np.stack([sx, sy, sz], axis=1)  # (nsta, 3, npt)
            if basis is None:
                basis = np.zeros((2, nsub) + h.shape[:2] + (npt // 2 + 1,), dtype=np.complex128)
            basis[r_idx, i] = np.fft.rfft(h * damp, axis=-1)
    return basis


def axitra_aw(ap) -> float:
    """Damping ``aw`` axitra actually uses: the one in ``<sid>.data``."""
    import re

    text = Path(f"{ap.sid}.data").read_text()
    match = re.search(r"aw\s*=\s*([-+0-9.eEdD]+)", text)
    if match is None:
        raise ValueError(f"no aw= in {ap.sid}.data")
    return float(match.group(1).lower().replace("d", "e").rstrip("."))


def _axitra_damping(aw: float, npt: int) -> np.ndarray:
    return np.exp(-np.pi * float(aw) * np.arange(npt) / npt)


def synthetics_from_basis(basis: np.ndarray, moment_x: np.ndarray, moment_z: np.ndarray, aw: float) -> np.ndarray:
    """Synthetics ``(nsta, 3, npt)`` for cumulative moment functions
    ``moment_x/z`` ``(nsub, npt)`` (Nm, axitra time grid), identical to what
    one axitra conv() per subfault would give (minus its %.3g file rounding).
    """
    npt = moment_x.shape[-1]
    damp = _axitra_damping(aw, npt)
    src = np.fft.rfft(np.stack([moment_x, moment_z]) * damp, axis=-1)  # (2, nsub, nf)
    return np.fft.irfft(np.einsum("rsf,rsjcf->jcf", src, basis), npt, axis=-1) / damp
