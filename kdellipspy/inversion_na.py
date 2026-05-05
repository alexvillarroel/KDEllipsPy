"""
Neighbourhood Algorithm (NA) inversion driver.
(Motor de inversión con el Algoritmo Neighbourhood — NA.)

This module contains ONLY the NA search implementation.
For MCMC inversion, see ``inversion_mcmc.py``.

Exports
-------
- NAConfig          : NA search hyperparameters
- NAInversionModel  : Full inversion driver; call ``run_na_search()``

Dependencies
------------
  inversion_base  (NAModel, NAResult, BaseInversionModel)
  neighpy         (pip install neighpy)
  config_parser, forward_model  (from this package)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

# Handle both relative imports (package) and direct imports (notebooks/scripts)
try:
    from .inversion_base import (
        BaseInversionModel,
        NAModel,
        NAResult,
    )
except ImportError:
    from inversion_base import (
        BaseInversionModel,
        NAModel,
        NAResult,
    )


# ---------------------------------------------------------------------------
# NAConfig
# ---------------------------------------------------------------------------

@dataclass
class NAConfig:
    """Hyperparameters for the Neighbourhood Algorithm search.
    (Hiperparámetros para la búsqueda con el Algoritmo Neighbourhood.)

    Attributes
    ----------
    n_samples_initial    : Number of random samples drawn in the first iteration (ni)
    n_samples_iteration  : Samples drawn around the best Voronoi cells per iteration (ns)
    n_iterations         : Number of NA iterations after the initial random stage (n)
    n_cells_resample     : Number of best Voronoi cells to resample from (nr)
    random_seed          : Optional seed for reproducibility
    keep_axitra_files    : Keep temporary AXITRA files (useful for debugging)
    """

    n_samples_initial: int = 30
    n_samples_iteration: int = 30
    n_iterations: int = 10
    n_cells_resample: int = 7
    random_seed: Optional[int] = None
    keep_axitra_files: bool = False


# ---------------------------------------------------------------------------
# NAInversionModel
# ---------------------------------------------------------------------------

class NAInversionModel(BaseInversionModel):
    """
    Kinematic Inversion Model using the Neighbourhood Algorithm (NA).
    (Modelo de Inversión Cinemática utilizando el Algoritmo Neighbourhood — NA.)

    Integrates observed data, arrival times, and event configuration to evaluate
    kinematic rupture models. Communicates with AXITRA to simulate synthetic
    seismograms and compute misfit against real data.
    (Integra los datos observados, tiempos de llegada y configuración del evento
    para evaluar distintos modelos de ruptura cinemática. Se comunica con AXITRA
    para simular sismogramas sintéticos y calcular el desajuste con los datos reales.)

    Parameters
    ----------
    input_ctl_path     : Path to input.ctl configuration file
    axitra_dir         : Path to AXITRA2024 binary directory (optional)
    observed_waveforms : 3-component seismograms, shape (nsta, 3, npts)
    time_array         : Time vector, shape (npts,)
    azi_times_array    : P/S arrival time table, shape (nsta, 3)

    Attributes
    ----------
    best_synthetics : np.ndarray | None
        Synthetic seismograms of the best model found so far (nsta, 3, npts).
        Set to None until the first objective function evaluation.
    param_names     : List[str]  — Human-readable parameter labels
    param_ranges    : np.ndarray, shape (n_params, 2) — [min, max] per parameter

    Example
    -------
    >>> model = NAInversionModel(
    ...     "path/to/run/input.ctl",
    ...     axitra_dir="path/to/AXITRA2024",
    ...     observed_waveforms=obs,
    ...     time_array=t,
    ...     azi_times_array=azi,
    ... )
    >>> result = model.run_na_search()
    >>> print(result.best_model.misfit)
    """

    # NAInversionModel inherits __init__ directly from BaseInversionModel —
    # no additional state is required for the NA algorithm.

    # ------------------------------------------------------------------
    def run_na_search(self, na_config: Optional[NAConfig] = None) -> NAResult:
        """Run the Neighbourhood Algorithm search and return all sampled models.
        (Ejecuta la búsqueda con el Algoritmo Neighbourhood y retorna todos los modelos.)

        Parameters
        ----------
        na_config : NAConfig, optional
            Search hyperparameters.  If ``None``, values are read from
            ``self.cfg.inversion_process`` (input.ctl).

        Returns
        -------
        NAResult
            Container with every evaluated model, the best model, and JSON/CSV
            export helpers.

        Raises
        ------
        ImportError
            If ``neighpy`` is not installed.
        """
        if na_config is None:
            na_config = NAConfig(
                n_samples_initial=self.cfg.inversion_process.ss1,
                n_samples_iteration=self.cfg.inversion_process.ss_other,
                n_iterations=self.cfg.inversion_process.num_iterations,
                n_cells_resample=self.cfg.inversion_process.cells_resample,
            )

        # ------------------------------------------------------------------
        # Workaround for neighpy sampling bug when ns % nr != 0.
        # In neighpy.search.NASearcher, internal allocation uses ns // nr, but
        # the best Voronoi cell walks ns % nr extra steps, which can overflow.
        # We silently round ns up to the nearest multiple of nr.
        # ------------------------------------------------------------------
        nr = max(1, int(na_config.n_cells_resample))
        ns_requested = int(na_config.n_samples_iteration)
        ns_effective = max(nr, ns_requested)
        if ns_effective % nr != 0:
            ns_effective = ((ns_effective // nr) + 1) * nr
            print(
                f"[NA] Adjusting n_samples_iteration from {ns_requested} to {ns_effective} "
                f"to satisfy neighpy constraint (must be a multiple of n_cells_resample={nr}).",
                flush=True,
            )

        if na_config.random_seed is not None:
            np.random.seed(na_config.random_seed)

        # Store a (possibly adjusted) copy so objective_function can read n0/ns
        self._mcmc_cfg_runtime = None
        self._mcmc_step_index = None
        self._na_cfg_runtime = (
            NAConfig(
                n_samples_initial=int(na_config.n_samples_initial),
                n_samples_iteration=int(ns_effective),
                n_iterations=int(na_config.n_iterations),
                n_cells_resample=int(nr),
                random_seed=na_config.random_seed,
                keep_axitra_files=bool(na_config.keep_axitra_files),
            )
            if ns_effective != int(na_config.n_samples_iteration)
            else na_config
        )
        self._eval_count = 0
        self._best_misfit_seen = float("inf")
        self._axitra_id_counter = 0

        try:
            from neighpy import NASearcher
        except Exception as exc:
            raise ImportError(
                "neighpy is required for NA search. Install with: pip install neighpy"
            ) from exc

        bounds = tuple((float(lo), float(hi)) for lo, hi in self.param_ranges)

        searcher = NASearcher(
            self.objective_function,
            ns=ns_effective,
            nr=nr,
            ni=int(na_config.n_samples_initial),
            n=int(na_config.n_iterations),
            bounds=bounds,
        )

        expected_models = (
            int(na_config.n_samples_initial)
            + ns_effective * int(na_config.n_iterations)
        )
        print(
            f"[NA] Starting search: ni={int(na_config.n_samples_initial)}, "
            f"ns={ns_effective}, n={int(na_config.n_iterations)}, nr={nr} "
            f"-> expected evaluations={expected_models}",
            flush=True,
        )
        searcher.run()

        samples = np.asarray(searcher.samples)
        objectives = np.asarray(searcher.objectives)

        # Assign approximate iteration labels for compatibility with existing notebooks.
        n0 = int(na_config.n_samples_initial)
        models: List[NAModel] = []
        for i, (sample, obj) in enumerate(zip(samples, objectives)):
            iteration = 0 if i < n0 else 1 + ((i - n0) // max(1, ns_effective))
            models.append(
                NAModel(
                    model=np.asarray(sample, dtype=float),
                    misfit=float(obj),
                    iteration=iteration,
                )
            )

        self._na_cfg_runtime = None
        self._mcmc_step_index = None

        return NAResult(
            models,
            param_names=self.param_names,
            extra_metadata={"algorithm": "NA", "neighpy": True},
        )


__all__ = [
    "NAConfig",
    "NAInversionModel",
    # Re-export shared types so callers can import everything from one place
    "NAModel",
    "NAResult",
]
