"""Fast kinematic path (inversion/fast_axitra.py) must reproduce one axitra
conv() call with per-source moments, rupture delays and the triangle STF."""

import numpy as np
import pytest

from kdellipspy.inversion.fast_axitra import impulse_basis, kinematic_synthetics
from kdellipspy.inversion.dynamic.dynamic_convolution import axitra_aw
from kdellipspy.test_dynamic_convolution import AXITRA_DIR, _mock_config, needs_axitra


@needs_axitra
@pytest.mark.parametrize("source_type", [0, 4])
def test_fast_kinematic_matches_axitra_conv(source_type):
    from kdellipspy import AxitraForwardModel

    cfg = _mock_config()
    fm = AxitraForwardModel.from_config(cfg, axitra_dir=str(AXITRA_DIR))
    ap = fm.green(fm.build_axitra(fm.build_geometry(), latlon=False), quiet=True)
    from axitra import moment

    nsrc = ap.nsource
    rng = np.random.default_rng(0)
    hist = np.zeros((nsrc, 8))
    hist[:, 0] = np.arange(1, nsrc + 1)
    # Exact in the %.3g / %.3f that axitra.py writes to .hist.
    hist[:, 1] = rng.integers(1, 3, nsrc) * 1e18
    hist[:, 2:5] = [20.0, 60.0, 90.0]
    hist[:, 7] = rng.integers(0, 40, nsrc) * 0.125  # rupture delays (s)
    t0 = 2.0

    _, sx, sy, sz = moment.conv(ap, hist, source_type=source_type, t0=t0, unit=1)
    slow = np.stack([sx, sy, sz], axis=1)
    fast = kinematic_synthetics(impulse_basis(ap, hist, unit=1), hist[:, 1], hist[:, 7],
                                source_type, t0, axitra_aw(ap), float(ap.duration))
    ap.clean()

    assert np.abs(slow).max() > 1e-6
    assert np.abs(fast - slow).max() / np.abs(slow).max() < 1e-4
