"""NA-search driver for the dynamic-rupture forward model (Fase 4, second
half): wires :class:`DynamicForwardModel` into the SAME
``BaseInversionModel``/``NAInversionModel`` machinery the kinematic side
uses (misfit, ``NAResult``, checkpointing, ``run_na_search``), by overriding
only ``_evaluate_model`` — the one method that is genuinely kinematic-specific
(everything else in ``NAInversionModel`` only calls ``objective_function``,
which stays untouched).
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import numpy as np

from ..kinematic.model_na import NAInversionModel
from ...core.geometry import TSNFaultGridSpec
from .forward_model_dynamic import DynamicForwardModel
from .tsn_bridge import TSNRunConfig

logger = logging.getLogger(__name__)


class DynamicNAInversionModel(NAInversionModel):
    """Same NA search as :class:`NAInversionModel`, but the forward model is
    fd3d_TSN's dynamic rupture solver instead of the prescribed-STF ellipse.

    Extra constructor args beyond ``NAInversionModel``: ``tsn_run_cfg``
    (solver work directory + FD grid dims) and ``tsn_grid`` (physical FD
    geometry for the friction-coefficient conversion). See
    :mod:`kdellipspy.inversion.dynamic.forward_model_dynamic`.

    ``self.param_ranges``/``param_names`` must describe the 10 dynamic
    params (a,b,xo,yo,phi,Te,cte1,cte2,r,dmax — see
    ``EllipticalStressMapper``), in COARSE (nx,ny) grid-point units — the
    caller's ``cfg.inversion_params`` section is responsible for that (Fase
    5, dynamic input parser, still pending).
    """

    _DYNAMIC_PARAM_NAMES = [
        "a (pts)", "b (pts)", "xo (pts)", "yo (pts)", "phi (rad)",
        "Te (MPa)", "cte1", "cte2", "r (pts)", "Dc (m)",
    ]

    def __init__(self, *args, tsn_run_cfg: TSNRunConfig, tsn_grid: TSNFaultGridSpec, **kwargs):
        super().__init__(*args, **kwargs)
        self.dynamic_fm = DynamicForwardModel(
            self.cfg,
            tsn_run_cfg,
            tsn_grid,
            axitra_dir=str(self.fm.axitra_dir),
            axitra_aw=self.axitra_aw,
            axitra_ikmax=self.axitra_ikmax,
        )
        self.param_names = list(self._DYNAMIC_PARAM_NAMES)

    # ------------------------------------------------------------------
    def _evaluate_model(self, model: np.ndarray) -> Tuple[float, Optional[np.ndarray]]:
        """Dynamic-rupture counterpart of ``BaseInversionModel._evaluate_model``:
        same misfit/bandpass/exception-handling contract, different forward.
        """
        try:
            if self.misfit_calc is None:
                return 1e10, None

            synthetics = self.dynamic_fm.forward(model)

            from kdellipspy.core.signal_utils import bandpass_filter_waveforms

            n_int = 3 - int(self.cfg.observed_data.units)
            synthetics = bandpass_filter_waveforms(
                synthetics,
                self.time_array,
                freq1=float(self.cfg.ellipse.freq1),
                freq2=float(self.cfg.ellipse.freq2),
                corners=2 * n_int,
                zerophase=bool(getattr(self.cfg.ellipse, "zerophase", True)),
            )

            misfit = float(self.misfit_calc.l2_misfit(synthetics, use_full_signal=self.use_full_signal))
            return misfit, synthetics

        except Exception as exc:
            logger.error("Dynamic model evaluation failed: %s", exc)
            return 1e10, None

    def clean(self) -> None:
        self.dynamic_fm.clean()
