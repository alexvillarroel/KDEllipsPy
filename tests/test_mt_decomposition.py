"""
Regression tests for the moment-tensor -> elementary-mechanism decomposition
(core.geometry).

Background: a previous implementation hard-coded the 6 basis tensors in a frame
that did not match the GCMT/USE (RTP) convention of the input moment tensor, so
each MT component was routed to the wrong elementary mechanism (rt/rp/tp slots
cyclically permuted; dip-slip diagonal terms with wrong signs). The radiation
pattern came out polarity-inverted. These tests pin the convention so the bug
cannot silently return.

The key invariant: the decomposition basis must be the *actual* moment tensors
of the 6 elementary mechanisms (b_str/b_dip/b_rak). Therefore decomposing a
single mechanism's own moment tensor must return a unit weight in that
mechanism's slot and (almost) zero elsewhere.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from kdellipspy.core.geometry import (
    _MT_BASIS_DIP,
    _MT_BASIS_RAK,
    _MT_BASIS_STR,
    _elementary_mt_matrix_rtp,
    _mt_signed_coeffs,
    _sdr_to_mt_rtp,
)


def _mt(mrr, mtt, mpp, mrt, mrp, mtp, exponent=0.0):
    """Minimal stand-in for the config MomentTensor object."""
    return SimpleNamespace(
        mrr=mrr, mtt=mtt, mpp=mpp, mrt=mrt, mrp=mrp, mtp=mtp, exponent=exponent
    )


def test_basis_is_invertible():
    B = _elementary_mt_matrix_rtp()
    assert B.shape == (6, 6)
    assert abs(np.linalg.det(B)) > 1e-6


@pytest.mark.parametrize("k", range(5))
def test_mechanism_decomposes_to_its_own_slot(k):
    """Feeding elementary mechanism k's own MT must isolate slot k.

    This is exactly what the convention bug broke: the off-diagonal slots were
    permuted and the diagonal ones sign-flipped, so a mechanism decomposed into
    the *wrong* slots.
    """
    rr, tt, pp, rt, rp, tp = _sdr_to_mt_rtp(
        _MT_BASIS_STR[k], _MT_BASIS_DIP[k], _MT_BASIS_RAK[k]
    )
    coeffs = _mt_signed_coeffs(_mt(rr, tt, pp, rt, rp, tp))
    expected = np.zeros(6)
    expected[k] = 1.0
    np.testing.assert_allclose(coeffs, expected, atol=1e-9)


def test_isotropic_maps_to_isotropic_slot():
    """A pure explosion (identity MT) must land entirely in the isotropic slot 5."""
    coeffs = _mt_signed_coeffs(_mt(1.0, 1.0, 1.0, 0.0, 0.0, 0.0))
    expected = np.zeros(6)
    expected[5] = 1.0
    np.testing.assert_allclose(coeffs, expected, atol=1e-9)


def test_decomposition_reconstructs_arbitrary_tensor():
    """B @ coeffs must reproduce the input tensor (RTP layout)."""
    mt = _mt(-2.486, 7.932, 0.083, 0.038, -0.739, 0.473, exponent=18.0)
    coeffs = _mt_signed_coeffs(mt)
    recovered = _elementary_mt_matrix_rtp() @ coeffs
    scale = 10.0 ** mt.exponent
    expected = np.array([mt.mrr, mt.mtt, mt.mpp, mt.mrt, mt.mrp, mt.mtp]) * scale
    np.testing.assert_allclose(recovered, expected, rtol=1e-9, atol=1.0)


def test_sdr_to_mt_recovers_plane_via_obspy():
    """Round-trip check: the MT built by _sdr_to_mt_rtp, fed to obspy's
    independent mt2plane (RTP convention), must recover the input nodal plane
    (or its auxiliary)."""
    obspy_bb = pytest.importorskip("obspy.imaging.beachball")
    strike, dip, rake = 72.0, 47.0, -111.0
    rr, tt, pp, rt, rp, tp = _sdr_to_mt_rtp(strike, dip, rake)
    mt = obspy_bb.MomentTensor(rr, tt, pp, rt, rp, tp, 0)
    np1 = obspy_bb.mt2plane(mt)
    np2 = obspy_bb.aux_plane(np1.strike, np1.dip, np1.rake)

    def close(p):
        ds = min(abs(p[0] - strike), 360 - abs(p[0] - strike))
        dr = min(abs(p[2] - rake), 360 - abs(p[2] - rake))
        return ds < 2 and abs(p[1] - dip) < 2 and dr < 2

    planes = [(np1.strike, np1.dip, np1.rake), (np2[0], np2[1], np2[2])]
    assert any(close(p) for p in planes), f"input plane not recovered: {planes}"
