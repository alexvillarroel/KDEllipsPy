"""
MCMC inversion driver using PyMC + ArviZ.
(Motor de inversión MCMC usando PyMC + ArviZ.)

This module contains ONLY the MCMC search implementation.
For Neighbourhood Algorithm inversion, see ``inversion_na.py``.

Exports
-------
- MCMCConfig          : MCMC sampling hyperparameters
- MCMCInversionModel  : Full inversion driver; call ``run_mcmc_search()``

Dependencies
------------
  inversion_base  (NAModel, NAResult, BaseInversionModel)
  pymc            (pip install pymc)   — includes arviz and pytensor
  config_parser, forward_model  (from this package)
"""

from __future__ import annotations

import json
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
# MCMCConfig
# ---------------------------------------------------------------------------

@dataclass
class MCMCConfig:
    """Hyperparameters for the MCMC search with PyMC + ArviZ.
    (Hiperparámetros para la búsqueda MCMC con PyMC + ArviZ.)

    PyMC runs a Metropolis sampler with a uniform prior on the parameter box.
    ArviZ is used for posterior summary statistics.

    Attributes
    ----------
    total_steps     : Total number of MCMC steps per chain (tune + draws)
    burn_in         : Number of tuning steps (excluded from posterior)
    proposal_scale  : Metropolis proposal standard deviation as a fraction of
                      each parameter's range, e.g. 0.08 → 8 % of width per param
    thin            : Thinning factor applied to the posterior draws
    random_seed     : Optional seed for reproducibility
    keep_axitra_files : Keep temporary AXITRA files (useful for debugging)
    chains          : Number of parallel MCMC chains.
                      Warning: chains > 1 may fail when AXITRA uses shared temp files.
    """

    total_steps: int = 500
    burn_in: int = 0
    proposal_scale: float = 0.08
    thin: int = 1
    random_seed: Optional[int] = None
    keep_axitra_files: bool = False
    chains: int = 1


# ---------------------------------------------------------------------------
# MCMCInversionModel
# ---------------------------------------------------------------------------

class MCMCInversionModel(BaseInversionModel):
    """
    Kinematic Inversion Model using MCMC (Metropolis, PyMC + ArviZ).
    (Modelo de Inversión Cinemática usando MCMC — PyMC + ArviZ.)

    Integrates observed data, arrival times, and event configuration to sample
    the posterior distribution of kinematic rupture model parameters.
    (Integra los datos observados, tiempos de llegada y configuración del evento
    para muestrear la distribución posterior de los parámetros del modelo cinemático.)

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

    Example
    -------
    >>> model = MCMCInversionModel(
    ...     "path/to/run/input.ctl",
    ...     axitra_dir="path/to/AXITRA2024",
    ...     observed_waveforms=obs,
    ...     time_array=t,
    ...     azi_times_array=azi,
    ... )
    >>> result = model.run_mcmc_search()
    >>> print(result.best_model.misfit)
    """

    # MCMCInversionModel inherits __init__ from BaseInversionModel.

    # ------------------------------------------------------------------
    def run_mcmc_search(self, mc_config: Optional[MCMCConfig] = None) -> NAResult:
        """Run the PyMC Metropolis sampler and return posterior models.
        (Ejecuta el muestreador Metropolis de PyMC y retorna los modelos posteriores.)

        Parameters
        ----------
        mc_config : MCMCConfig, optional
            MCMC hyperparameters.  If ``None``, values are read from
            ``self.cfg.inversion_process`` (input.ctl).

        Returns
        -------
        NAResult
            Container with posterior samples, the best model (minimum misfit),
            ArviZ summary in ``extra_metadata``, and JSON/CSV export helpers.

        Raises
        ------
        ImportError
            If ``pymc`` (or ``arviz`` / ``pytensor``) is not installed.
        ValueError
            If ``burn_in >= total_steps``.
        """
        if mc_config is None:
            ip = self.cfg.inversion_process
            mc_config = MCMCConfig(
                total_steps=int(ip.mcmc_total_steps),
                burn_in=int(ip.mcmc_burn_in),
                proposal_scale=float(ip.mcmc_proposal_scale),
                thin=int(ip.mcmc_thin),
                random_seed=None,
                keep_axitra_files=False,
                chains=int(ip.mcmc_chains),
            )

        if int(mc_config.total_steps) < 1:
            raise ValueError("MCMCConfig.total_steps must be >= 1")
        if int(mc_config.burn_in) >= int(mc_config.total_steps):
            raise ValueError("MCMCConfig.burn_in must be < total_steps")

        return self._run_mcmc_pymc(mc_config)

    # ------------------------------------------------------------------
    def _run_mcmc_pymc(self, mc_config: MCMCConfig) -> NAResult:
        """Internal: PyMC Metropolis over uniform prior; ArviZ for posterior summary.
        (Interno: Metropolis de PyMC sobre prior uniforme; ArviZ para el resumen posterior.)
        """
        try:
            import arviz as az
            import pymc as pm
            import pytensor.tensor as pt
            from pytensor.graph.basic import Apply
            from pytensor.graph.op import Op
        except ImportError as exc:
            raise ImportError(
                "MCMC inversion requires pymc (includes arviz and pytensor). "
                "Install with: pip install pymc"
            ) from exc

        thin = max(1, int(mc_config.thin))
        tune = int(mc_config.burn_in)
        total = int(mc_config.total_steps)
        draws = max(1, (total - tune) // thin)
        chains = max(1, int(mc_config.chains))

        if tune >= total:
            raise ValueError("MCMCConfig.burn_in must be < total_steps")
        if chains > 1:
            print(
                "[MCMC/PyMC] Warning: chains > 1 may fail with AXITRA (shared temp files). "
                "If errors occur, set 'MCMC chains : 1' in input.ctl.",
                flush=True,
            )

        low = self.param_ranges[:, 0].astype(np.float64)
        high = self.param_ranges[:, 1].astype(np.float64)
        width = np.maximum(high - low, 1e-30)
        prop_sd = float(mc_config.proposal_scale) * width
        S_prop = np.diag(prop_sd ** 2)

        # Keep a reference for the inner Op (avoids 'self' cell closure issues)
        inv = self

        # ------------------------------------------------------------------
        # PyTensor Op that wraps objective_function as a black-box log-likelihood.
        # Gradient is zero → PyMC will use Metropolis (no NUTS).
        # ------------------------------------------------------------------
        class _NegLogLike(Op):
            """log p(data|theta) ∝ -misfit(theta); no gradient (Metropolis).
            (log p(datos|theta) ∝ -misfit(theta); sin gradiente — usa Metropolis.)
            """

            def make_node(self, theta):
                theta = pt.as_tensor_variable(theta)
                return Apply(self, [theta], [pt.dscalar()])

            def perform(self, node, inputs, outputs):
                (theta,) = inputs
                inv._pymc_inner = True
                try:
                    m = inv.objective_function(
                        np.asarray(theta, dtype=np.float64).ravel()
                    )
                    outputs[0][0] = np.asarray(-float(m), dtype=np.float64)
                finally:
                    inv._pymc_inner = False

            def grad(self, inputs, output_grads):
                (theta,) = inputs
                return [pt.zeros_like(theta)]

        neg_op = _NegLogLike()
        coords = {"par": np.asarray(self.param_names, dtype=str)}

        # Reset runtime counters
        self._na_cfg_runtime = None
        self._mcmc_cfg_runtime = mc_config
        self._eval_count = 0
        self._best_misfit_seen = float("inf")
        self._axitra_id_counter = 0
        self._mcmc_step_index = None
        self._pymc_inner = False

        print(
            f"[MCMC/PyMC] tune={tune} draws={draws} thin(post)={thin} chains={chains} "
            f"(~{total} steps per chain aligned with MCMCConfig.total_steps)",
            flush=True,
        )

        # ------------------------------------------------------------------
        # Build and sample the PyMC model
        # ------------------------------------------------------------------
        with pm.Model(coords=coords):
            theta = pm.Uniform("theta", lower=low, upper=high, dims="par")
            pm.Potential("loglike", neg_op(theta))

            try:
                step = pm.Metropolis([theta], S=S_prop)
            except TypeError:
                # Older PyMC versions may not accept the S argument
                step = pm.Metropolis([theta])

            idata = pm.sample(
                draws=draws,
                tune=tune,
                step=step,
                chains=chains,
                cores=1,
                random_seed=mc_config.random_seed,
                return_inferencedata=True,
                compute_convergence_checks=bool(chains > 1),
                progressbar=True,
            )

        # Post-sampling thinning (if thin > 1)
        if thin > 1:
            try:
                idata = idata.sel(draw=slice(None, None, thin))
            except Exception:
                idata = idata.isel(draw=slice(None, None, thin))

        # ------------------------------------------------------------------
        # ArviZ summary (mean, sd, hdi, etc.)
        # ------------------------------------------------------------------
        summary_df = az.summary(idata, var_names=["theta"], kind="stats")
        arviz_records = summary_df.reset_index().rename(columns={"index": "var"})
        try:
            js = arviz_records.to_json(
                orient="records", double_precision=15, default_handler=str
            )
        except TypeError:
            js = arviz_records.to_json(orient="records", double_precision=15)
        arviz_summary = json.loads(js)

        # ------------------------------------------------------------------
        # Re-evaluate forward model for every posterior sample to align misfits.
        # (Re-evalúa el modelo forward en cada muestra posterior para alinear misfits.)
        # PyMC does not expose exact log-likelihood values for black-box Ops.
        # ------------------------------------------------------------------
        post = idata.posterior["theta"]
        models: List[NAModel] = []
        idx = 0
        print(
            "[MCMC/PyMC] Re-evaluating forward model on posterior samples to align misfit "
            "in export. (Re-evaluando el forward en muestras posteriores para exportar misfit.)",
            flush=True,
        )
        self._mcmc_step_index = None
        self._pymc_inner = False

        for c in range(int(post.sizes["chain"])):
            for d in range(int(post.sizes["draw"])):
                th = post.isel(chain=c, draw=d).values.astype(np.float64).ravel()
                mf = float(self.objective_function(th))
                models.append(NAModel(model=th.copy(), misfit=mf, iteration=idx))
                idx += 1

        # Clean up runtime state
        self._mcmc_step_index = None
        self._mcmc_cfg_runtime = None

        return NAResult(
            models,
            param_names=self.param_names,
            extra_metadata={
                "algorithm": "MCMC",
                "mcmc_total_steps": total,
                "mcmc_burn_in": tune,
                "mcmc_thin": thin,
                "mcmc_proposal_scale": float(mc_config.proposal_scale),
                "mcmc_chains": chains,
                "mcmc_stored_samples": len(models),
                "arviz_summary": arviz_summary,
            },
        )


__all__ = [
    "MCMCConfig",
    "MCMCInversionModel",
    # Re-export shared types so callers can import everything from one place
    "NAModel",
    "NAResult",
]
