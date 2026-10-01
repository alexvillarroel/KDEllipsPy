"""Fast kinematic synthetics: axitra's conv reproduced in numpy from a
per-source impulse basis computed once.

Per forward evaluation the kinematic model only changes each fixed source's
scalar moment and rupture delay (the mesh, medium, stations and every
source's strike/dip/rake stay fixed — also in MT mode, where each elementary
source keeps its mechanism). axitra's conv (convmPy.F90) is, per frequency
bin, ``u(w) = sum_i M0_i * fsource(w_c) * exp(-i*w_c*delay_i) * G_i(w)`` at the
complex frequency ``w_c = 2*pi*f - i*pi*aw/T``, then IFFT and undamping by
``exp(+pi*aw*t/T)``. Its response to a unit Dirac source (fsource = 1, no
delay) is ``h_i = G_i`` undamped, so storing ``rfft(d*h_i)`` with
``d = exp(-pi*aw*t/T)`` reproduces axitra exactly (the same identity as
dynamic_convolution.response_basis). One conv per source once, then each model
is a numpy sum — replaces one conv call (and its .hist file) per evaluation.
"""

from __future__ import annotations

import numpy as np

from .dynamic.dynamic_convolution import _axitra_damping, axitra_aw

SUPPORTED_SOURCE_TYPES = (0, 4)  # Dirac, integral of a triangle (fsource.f90)


def impulse_basis(ap, hist_template: np.ndarray, unit: int) -> np.ndarray:
    """Damped spectra ``(nsrc, nsta, 3, npt//2+1)`` of each source's response
    to a unit-moment Dirac source, with the mechanism of its ``hist`` row."""
    try:
        from axitra import moment
    except ImportError:
        import sys

        sys.path.append(str(ap.axpath))
        from axitra import moment

    nsrc, npt = hist_template.shape[0], ap.npt
    damp = _axitra_damping(axitra_aw(ap), npt)
    basis = None
    for i in range(nsrc):
        hist = hist_template.copy()
        hist[:, 1] = 0.0
        hist[:, 7] = 0.0
        hist[i, 1] = 1.0
        result = moment.conv(ap, hist, source_type=0, t0=0.0, unit=unit)
        if result is None:
            raise RuntimeError(f"moment.conv failed for source {i + 1}")
        _, sx, sy, sz = result
        h = np.stack([sx, sy, sz], axis=1)
        if basis is None:
            basis = np.zeros((nsrc,) + h.shape[:2] + (npt // 2 + 1,), dtype=np.complex128)
        basis[i] = np.fft.rfft(h * damp, axis=-1)
    return basis


def fsource_spectrum(source_type: int, t0: float, omega_c: np.ndarray) -> np.ndarray:
    """axitra's fsource (fsource.f90) at complex angular frequencies."""
    if source_type == 0:
        return np.ones_like(omega_c)
    if source_type == 4:
        iw = 1j * omega_c
        ramp = (1.0 - np.exp(-iw * t0)) / (iw * t0)
        half = iw * t0 / 2.0
        finite = (np.exp(half) - np.exp(-half)) / half / 2.0
        return ramp * finite / iw
    raise ValueError(f"source_type {source_type} not supported by the fast path")


def kinematic_synthetics(basis: np.ndarray, moments: np.ndarray, delays: np.ndarray,
                         source_type: int, t0: float, aw: float, duration: float) -> np.ndarray:
    """Synthetics ``(nsta, 3, npt)`` for per-source moments (N·m) and delays (s)."""
    npt = 2 * (basis.shape[-1] - 1)
    omega_c = 2.0 * np.pi * np.arange(basis.shape[-1]) / duration - 1j * np.pi * aw / duration
    src = np.asarray(moments)[:, None] * np.exp(-1j * np.outer(delays, omega_c))  # (nsrc, nf)
    spec = fsource_spectrum(source_type, t0, omega_c) * np.einsum("sf,sjcf->jcf", src, basis)
    return np.fft.irfft(spec, npt, axis=-1) / _axitra_damping(aw, npt)
