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
    from ..base import (
        BaseInversionModel,
        NAModel,
        NAResult,
    )
except ImportError:
    try:
        from kdellipspy.inversion.base import (
            BaseInversionModel,
            NAModel,
            NAResult,
        )
    except ImportError:
        from base import (
            BaseInversionModel,
            NAModel,
            NAResult,
        )


# ---------------------------------------------------------------------------
# NAConfig
# ---------------------------------------------------------------------------

try:
    from joblib import Parallel, delayed
except ImportError:
    Parallel = None

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
    n_jobs               : Number of parallel workers for objective function evaluation (-1 for all)
    random_seed          : Optional seed for reproducibility
    keep_axitra_files    : Keep temporary axitra files (useful for debugging)
    """

    n_samples_initial: int = 30
    n_samples_iteration: int = 30
    n_iterations: int = 10
    n_cells_resample: int = 7
    n_jobs: int = 1
    random_seed: Optional[int] = None
    keep_axitra_files: bool = False


# ---------------------------------------------------------------------------
# ParallelNASearcher
# ---------------------------------------------------------------------------

try:
    from neighpy import NASearcher
    
    class ParallelNASearcher(NASearcher):
        """
        Subclass of NASearcher that parallelizes the objective function evaluation.
        (Subclase de NASearcher que paraleliza la evaluación de la función objetivo.)
        """
        def __init__(self, *args, n_jobs: int = 1, inversion_model: Optional[BaseInversionModel] = None, **kwargs):
            super().__init__(*args, **kwargs)
            self.n_jobs = n_jobs
            self.inversion_model = inversion_model

        def _update_ensemble(self, new_samples: np.ndarray):
            n = new_samples.shape[0]
            self.samples[self.np : self.np + n] = new_samples
            
            if self.n_jobs == 1 or Parallel is None or self.inversion_model is None:
                # Sequential fallback
                for i in range(n):
                    self.objectives[self.np + i] = self.objective(new_samples[i])
            else:
                # Parallel evaluation
                results = Parallel(n_jobs=self.n_jobs)(
                    delayed(self.inversion_model._evaluate_model)(new_samples[i]) 
                    for i in range(n)
                )
                
                for i, (misfit, synthetics) in enumerate(results):
                    self.inversion_model._eval_count += 1
                    if misfit < self.inversion_model._best_misfit_seen:
                        self.inversion_model._best_misfit_seen = misfit
                        if synthetics is not None:
                            self.inversion_model.best_synthetics = synthetics.copy()
                    
                    self.objectives[self.np + i] = misfit

                    # Logging in parallel mode
                    iter_est = 0
                    if self.inversion_model._na_cfg_runtime is not None:
                        n0 = int(self.inversion_model._na_cfg_runtime.n_samples_initial)
                        ns = max(1, int(self.inversion_model._na_cfg_runtime.n_samples_iteration))
                        if self.inversion_model._eval_count > n0:
                            iter_est = 1 + ((self.inversion_model._eval_count - n0 - 1) // ns)
                    print(
                        f"[NA-Parallel] iter={iter_est:05d} eval={self.inversion_model._eval_count:05d} "
                        f"misfit={misfit:.6e} best={self.inversion_model._best_misfit_seen:.6e}",
                        flush=True,
                    )

            self.np += n
except ImportError:
    NASearcher = object
    ParallelNASearcher = object


# ---------------------------------------------------------------------------
# NAInversionModel
# ---------------------------------------------------------------------------

class NAInversionModel(BaseInversionModel):
    """
    Kinematic Inversion Model using the Neighbourhood Algorithm (NA).
    (Modelo de Inversión Cinemática utilizando el Algoritmo Neighbourhood — NA.)

    Integrates observed data, arrival times, and event configuration to evaluate
    kinematic rupture models. Communicates with axitra to simulate synthetic
    seismograms and compute misfit against real data.
    (Integra los datos observados, tiempos de llegada y configuración del evento
    para evaluar distintos modelos de ruptura cinemática. Se comunica con axitra
    para simular sismogramas sintéticos y calcular el desajuste con los datos reales.)

    Parameters
    ----------
    input_ctl_path     : Path to input.ctl configuration file (optional if config is provided)
    axitra_dir         : Path to axitra binary directory (optional)
    observed_waveforms : 3-component seismograms, shape (nsta, 3, npts)
    time_array         : Time vector, shape (npts,)
    azi_times_array    : P/S arrival time table, shape (nsta, 3)
    config             : ConfigParser object (optional if input_ctl_path is provided)

    Attributes
    ----------
    best_synthetics : np.ndarray | None
        Synthetic seismograms of the best model found so far (nsta, 3, npts).
        Set to None until the first objective function evaluation.
    param_names     : List[str]  — Human-readable parameter labels
    param_ranges    : np.ndarray, shape (n_params, 2) — [min, max] per parameter

    Example
    -------
    >>> # Using a config file path:
    >>> model = NAInversionModel(
    ...     "path/to/run/input.ctl",
    ...     axitra_dir="path/to/axitra",
    ...     observed_waveforms=obs,
    ...     time_array=t,
    ...     azi_times_array=azi,
    ... )
    >>> # Or using a ConfigParser object:
    >>> cfg = ConfigParser("path/to/input.ctl")
    >>> model = NAInversionModel(
    ...     config=cfg,
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
            # Safely get n_jobs from config, defaulting to 1 if the field is missing
            n_jobs_cfg = 1
            if self.cfg.inversion_process is not None:
                n_jobs_cfg = getattr(self.cfg.inversion_process, "n_jobs", 1)

            na_config = NAConfig(
                n_samples_initial=self.cfg.inversion_process.ss1,
                n_samples_iteration=self.cfg.inversion_process.ss_other,
                n_iterations=self.cfg.inversion_process.num_iterations,
                n_cells_resample=self.cfg.inversion_process.cells_resample,
                n_jobs=n_jobs_cfg,
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
                n_jobs=int(na_config.n_jobs),
                random_seed=na_config.random_seed,
                keep_axitra_files=bool(na_config.keep_axitra_files),
            )
            if ns_effective != int(na_config.n_samples_iteration)
            else na_config
        )
        self._eval_count = 0
        self._best_misfit_seen = float("inf")
        self._axitra_id_counter = 0
        # Reset any stale Green's-function cache from a previous search.
        self.clear_green_cache()

        if NASearcher is object:
            raise ImportError(
                "neighpy is required for NA search. Install with: pip install neighpy"
            )

        bounds = tuple((float(lo), float(hi)) for lo, hi in self.param_ranges)

        # Use ParallelNASearcher if n_jobs != 1
        if na_config.n_jobs != 1 and Parallel is not None:
            searcher_cls = ParallelNASearcher
            extra_kwargs = {"n_jobs": na_config.n_jobs, "inversion_model": self}
        else:
            searcher_cls = NASearcher
            extra_kwargs = {}

        searcher = searcher_cls(
            self.objective_function,
            ns=ns_effective,
            nr=nr,
            ni=int(na_config.n_samples_initial),
            n=int(na_config.n_iterations),
            bounds=bounds,
            **extra_kwargs
        )

        expected_models = (
            int(na_config.n_samples_initial)
            + ns_effective * int(na_config.n_iterations)
        )
        print(
            f"[NA] Starting search: ni={int(na_config.n_samples_initial)}, "
            f"ns={ns_effective}, n={int(na_config.n_iterations)}, nr={nr} "
            f"-> expected evaluations={expected_models} (n_jobs={na_config.n_jobs})",
            flush=True,
        )
        searcher.run(parallel=(na_config.n_jobs != 1))

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
        # Clean up the cached Green's-function axitra files.
        self.clear_green_cache()

        return NAResult(
            models,
            param_names=self.param_names,
            extra_metadata={"algorithm": "NA", "neighpy": True},
            best_synthetics=self.best_synthetics,
            observed=self.observed_waveforms,
            time=self.time_array,
            config=self.cfg,
            azi_times_array=self.azi_times_array,
        )


__all__ = [
    "NAConfig",
    "NAInversionModel",
    # Re-export shared types so callers can import everything from one place
    "NAModel",
    "NAResult",
]
