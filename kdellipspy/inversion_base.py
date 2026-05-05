"""
Base module for kinematic inversion: shared dataclasses, misfit logic and abstract model.
(Módulo base para inversión cinemática: dataclasses compartidos, lógica de misfit y modelo base.)

Exports
-------
- NAModel          : Single sampled model + misfit
- MisfitCalculator : L2 waveform misfit with P/S time windows
- NAResult         : Container for all sampled models with export helpers
- BaseInversionModel : Abstract base class for NA and MCMC inversion drivers

Dependencies (solo stdlib + numpy)
-----------------------------------
  numpy, pathlib, csv, json, datetime, logging, os, time, copy
"""

from __future__ import annotations

import csv
import json
import logging
import os
import time
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

# Handle both relative imports (package) and direct imports (notebooks/scripts)
try:
    from .config_parser import ConfigParser
    from .forward_model import AxitraForwardModel
except ImportError:
    from config_parser import ConfigParser
    from forward_model import AxitraForwardModel

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class NAModel:
    """Single sampled model and its objective value.
    (Modelo muestreado individual y su valor objetivo.)

    Attributes
    ----------
    model     : Parameter vector (np.ndarray)
    misfit    : Objective/misfit value (float)
    iteration : Iteration index at which this model was sampled (int)
    """

    model: np.ndarray
    misfit: float
    iteration: int


# ---------------------------------------------------------------------------
# MisfitCalculator
# ---------------------------------------------------------------------------

class MisfitCalculator:
    """L2 waveform misfit calculator with P/S arrival-time windows.
    (Calculador de desajuste L2 de formas de onda con ventanas de llegada P/S.)

    Parameters
    ----------
    observed_waveforms : np.ndarray, shape (nsta, 3, npts)
    time_array         : np.ndarray, shape (npts,)
    azi_times_path     : Path to ASCII azi_times.txt (3 columns: azi, tP, tS)
    azi_times_array    : Pre-loaded azi_times as np.ndarray (takes priority)
    time_window_s      : Duration (s) of the P and S analysis windows
    """

    def __init__(
        self,
        observed_waveforms: np.ndarray,
        time_array: np.ndarray,
        azi_times_path: Optional[Path] = None,
        azi_times_array: Optional[np.ndarray] = None,
        time_window_s: float = 20.0,
    ):
        self.observed = observed_waveforms
        self.time = time_array
        self.time_window_s = float(time_window_s)

        if azi_times_array is not None:
            arr = np.asarray(azi_times_array, dtype=float)
        elif azi_times_path is not None:
            if not azi_times_path.exists():
                raise FileNotFoundError(f"Missing required azi_times file: {azi_times_path}")
            arr = np.loadtxt(str(azi_times_path), dtype=float)
        else:
            raise ValueError("Provide either azi_times_array (in-memory) or azi_times_path.")

        arr = np.atleast_2d(arr)
        if arr.shape[1] < 3:
            src_desc = "azi_times_array" if azi_times_array is not None else str(azi_times_path)
            raise ValueError(
                f"Invalid azi_times format in {src_desc}: expected 3 columns (azi, tP, tS)"
            )

        self.azi = arr[:, 0]
        self.tp = arr[:, 1]
        self.ts = arr[:, 2]

    # ------------------------------------------------------------------
    def diagnostics_summary(self, synthetic: np.ndarray, max_stations: int = 3) -> str:
        """Compact per-station diagnostics string for one model evaluation.
        (Resumen diagnóstico compacto por estación para una evaluación de modelo.)
        """
        if synthetic.shape != self.observed.shape:
            return (
                f"[MISFIT DIAG] shape mismatch: synthetic {synthetic.shape} "
                f"vs observed {self.observed.shape}"
            )

        nsta, _, npts = self.observed.shape
        dt = float(self.time[1] - self.time[0]) if len(self.time) > 1 else 1.0
        sampling = max(1, int(np.rint(1.0 / dt)))
        win = max(1, int(np.rint(self.time_window_s * sampling)))

        obs_rms_global = float(np.sqrt(np.mean(self.observed ** 2)))
        syn_rms_global = float(np.sqrt(np.mean(synthetic ** 2)))
        global_ratio = syn_rms_global / max(obs_rms_global, 1e-30)

        obs_energy = 0.0
        syn_energy = 0.0
        lines = [
            (
                "[MISFIT DIAG] "
                f"rms_global(obs)={obs_rms_global:.3e} "
                f"rms_global(syn)={syn_rms_global:.3e} "
                f"syn/obs={global_ratio:.3e}"
            )
        ]

        nshow = min(max_stations, nsta)
        for j in range(nsta):
            az = float(self.azi[j])
            start_p = int(np.rint(float(self.tp[j]))) * sampling
            start_s = int(np.rint(float(self.ts[j]))) * sampling
            kp0 = max(0, start_p - 1)
            kp1 = min(npts - 1, start_p + win - 1)
            ks0 = max(0, start_s - 1)
            ks1 = min(npts - 1, start_s + win - 1)

            r_obs_rms = r_syn_rms = z_obs_rms = z_syn_rms = 0.0
            if kp1 >= kp0:
                x_obs = self.observed[j, 0, kp0:kp1 + 1]
                y_obs = self.observed[j, 1, kp0:kp1 + 1]
                z_obs = self.observed[j, 2, kp0:kp1 + 1]
                x_syn = synthetic[j, 0, kp0:kp1 + 1]
                y_syn = synthetic[j, 1, kp0:kp1 + 1]
                z_syn = synthetic[j, 2, kp0:kp1 + 1]
                r_obs = x_obs * np.cos(az) + y_obs * np.sin(az)
                r_syn = x_syn * np.cos(az) + y_syn * np.sin(az)
                r_obs_rms = float(np.sqrt(np.mean(r_obs ** 2)))
                r_syn_rms = float(np.sqrt(np.mean(r_syn ** 2)))
                z_obs_rms = float(np.sqrt(np.mean(z_obs ** 2)))
                z_syn_rms = float(np.sqrt(np.mean(z_syn ** 2)))
                obs_energy += float(np.sum(r_obs ** 2) + np.sum(z_obs ** 2))
                syn_energy += float(np.sum(r_syn ** 2) + np.sum(z_syn ** 2))

            t_obs_rms = t_syn_rms = 0.0
            if ks1 >= ks0:
                x_obs = self.observed[j, 0, ks0:ks1 + 1]
                y_obs = self.observed[j, 1, ks0:ks1 + 1]
                x_syn = synthetic[j, 0, ks0:ks1 + 1]
                y_syn = synthetic[j, 1, ks0:ks1 + 1]
                t_obs = y_obs * np.cos(az) - x_obs * np.sin(az)
                t_syn = y_syn * np.cos(az) - x_syn * np.sin(az)
                t_obs_rms = float(np.sqrt(np.mean(t_obs ** 2)))
                t_syn_rms = float(np.sqrt(np.mean(t_syn ** 2)))
                obs_energy += float(np.sum(t_obs ** 2))
                syn_energy += float(np.sum(t_syn ** 2))

            if j < nshow:
                lines.append(
                    (
                        f"[MISFIT DIAG] sta={j+1:02d} "
                        f"P(R): obs={r_obs_rms:.3e} syn={r_syn_rms:.3e} | "
                        f"P(Z): obs={z_obs_rms:.3e} syn={z_syn_rms:.3e} | "
                        f"S(T): obs={t_obs_rms:.3e} syn={t_syn_rms:.3e}"
                    )
                )

        window_ratio = syn_energy / max(obs_energy, 1e-30)
        lines.append(
            (
                "[MISFIT DIAG] "
                f"window_energy(obs)={obs_energy:.3e} "
                f"window_energy(syn)={syn_energy:.3e} "
                f"syn/obs={window_ratio:.3e}"
            )
        )
        return "\n".join(lines)

    # ------------------------------------------------------------------
    def l2_misfit(self, synthetic: np.ndarray) -> float:
        """Compute normalised L2 misfit in P (radial+vertical) and S (transverse) windows.
        (Calcula el desajuste L2 normalizado en ventanas P (radial+vertical) y S (transversal).)
        """
        if synthetic.shape != self.observed.shape:
            raise ValueError(
                f"Shape mismatch: synthetic {synthetic.shape} vs observed {self.observed.shape}"
            )

        nsta, ncomp, npts = self.observed.shape
        if ncomp != 3:
            raise ValueError("Expected 3 components per station (x,y,z).")
        if len(self.azi) != nsta:
            raise ValueError(
                f"azi_times rows ({len(self.azi)}) do not match number of stations ({nsta})"
            )

        dt = float(self.time[1] - self.time[0]) if len(self.time) > 1 else 1.0
        sampling = max(1, int(np.rint(1.0 / dt)))
        win = max(1, int(np.rint(self.time_window_s * sampling)))

        num = 0.0
        den = 0.0

        for j in range(nsta):
            az = float(self.azi[j])
            start_p = int(np.rint(float(self.tp[j]))) * sampling
            start_s = int(np.rint(float(self.ts[j]))) * sampling
            kp0 = max(0, start_p - 1)
            kp1 = min(npts - 1, start_p + win - 1)
            ks0 = max(0, start_s - 1)
            ks1 = min(npts - 1, start_s + win - 1)

            # P window: radial + vertical
            if kp1 >= kp0:
                x_obs = self.observed[j, 0, kp0:kp1 + 1]
                y_obs = self.observed[j, 1, kp0:kp1 + 1]
                z_obs = self.observed[j, 2, kp0:kp1 + 1]
                x_syn = synthetic[j, 0, kp0:kp1 + 1]
                y_syn = synthetic[j, 1, kp0:kp1 + 1]
                z_syn = synthetic[j, 2, kp0:kp1 + 1]
                r_obs = x_obs * np.cos(az) + y_obs * np.sin(az)
                r_syn = x_syn * np.cos(az) + y_syn * np.sin(az)
                num += float(np.sum((r_obs - r_syn) ** 2))
                den += float(np.sum(r_obs ** 2))
                num += float(np.sum((z_obs - z_syn) ** 2))
                den += float(np.sum(z_obs ** 2))

            # S window: transverse
            if ks1 >= ks0:
                x_obs = self.observed[j, 0, ks0:ks1 + 1]
                y_obs = self.observed[j, 1, ks0:ks1 + 1]
                x_syn = synthetic[j, 0, ks0:ks1 + 1]
                y_syn = synthetic[j, 1, ks0:ks1 + 1]
                t_obs = y_obs * np.cos(az) - x_obs * np.sin(az)
                t_syn = y_syn * np.cos(az) - x_syn * np.sin(az)
                num += float(np.sum((t_obs - t_syn) ** 2))
                den += float(np.sum(t_obs ** 2))

        return num / den if den > 0.0 else num


# ---------------------------------------------------------------------------
# NAResult
# ---------------------------------------------------------------------------

class NAResult:
    """Container for all sampled models with export helpers.
    (Contenedor de todos los modelos muestreados con exportadores.)

    Parameters
    ----------
    all_models     : List of NAModel instances
    param_names    : Parameter name labels (optional, defaults to 7-param names)
    extra_metadata : Algorithm-specific metadata written to JSON export
    """

    _DEFAULT_PARAM_NAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]

    def __init__(
        self,
        all_models: List[NAModel],
        param_names: Optional[Sequence[str]] = None,
        extra_metadata: Optional[Dict[str, Any]] = None,
    ):
        self.all_models = all_models
        self.best_model = min(all_models, key=lambda m: m.misfit) if all_models else None
        self.param_names = (
            list(param_names) if param_names is not None else list(self._DEFAULT_PARAM_NAMES)
        )
        self.extra_metadata = dict(extra_metadata) if extra_metadata else {}

    def export_results(self, filepath: Path) -> None:
        """Export all models + metadata as JSON.
        (Exporta todos los modelos y metadatos como JSON.)
        """
        meta: Dict[str, Any] = {
            "timestamp": datetime.now().isoformat(),
            "n_models": len(self.all_models),
            "best_misfit": self.best_model.misfit if self.best_model else None,
            "param_names": list(self.param_names),
        }
        meta.update(self.extra_metadata)
        payload = {
            "metadata": meta,
            "models": [
                {
                    "model": m.model.tolist(),
                    "misfit": float(m.misfit),
                    "iteration": int(m.iteration),
                }
                for m in self.all_models
            ],
        }
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with filepath.open("w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    def export_csv(self, filepath: Path) -> None:
        """Export all models as CSV (iteration, misfit, param columns).
        (Exporta todos los modelos como CSV con columnas de iteración, misfit y parámetros.)
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = ["iteration", "misfit", *self.param_names]
        with filepath.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for model in self.all_models:
                row: Dict[str, Any] = {
                    "iteration": int(model.iteration),
                    "misfit": float(model.misfit),
                }
                for name, value in zip(self.param_names, model.model):
                    row[name] = float(value)
                writer.writerow(row)


# ---------------------------------------------------------------------------
# BaseInversionModel
# ---------------------------------------------------------------------------

class BaseInversionModel:
    """Abstract base class for kinematic inversion drivers (NA and MCMC).
    (Clase base abstracta para los motores de inversión cinemática NA y MCMC.)

    Sub-classes must call ``super().__init__(...)`` and then implement either
    ``run_na_search`` or ``run_mcmc_search``.

    Parameters
    ----------
    input_ctl_path     : Path to input.ctl configuration file
    axitra_dir         : Path to AXITRA2024 binary directory (optional)
    observed_waveforms : Observed 3-component seismograms, shape (nsta, 3, npts)
    time_array         : Time vector, shape (npts,)
    azi_times_array    : Pre-computed P/S arrival time table, shape (nsta, 3)

    Attributes
    ----------
    best_synthetics : np.ndarray | None  — Updated whenever a new best misfit is found
    """

    def __init__(
        self,
        input_ctl_path: str,
        axitra_dir: Optional[str] = None,
        observed_waveforms: Optional[np.ndarray] = None,
        time_array: Optional[np.ndarray] = None,
        azi_times_array: Optional[np.ndarray] = None,
    ):
        self.input_ctl_path = Path(input_ctl_path).resolve()
        self.cfg = ConfigParser(str(self.input_ctl_path))
        self.fm = AxitraForwardModel(str(self.input_ctl_path), axitra_dir=axitra_dir)
        self.base_geometry = self.fm.build_geometry()

        self.observed_waveforms = observed_waveforms
        self.time_array = time_array
        self.misfit_calc: Optional[MisfitCalculator] = None

        if observed_waveforms is not None and time_array is not None:
            if azi_times_array is not None:
                self.misfit_calc = MisfitCalculator(
                    observed_waveforms,
                    time_array,
                    azi_times_array=azi_times_array,
                    time_window_s=20.0,
                )
            else:
                azi_times_path = self.input_ctl_path.parent / "Event" / "azi_times.txt"
                self.misfit_calc = MisfitCalculator(
                    observed_waveforms,
                    time_array,
                    azi_times_path=azi_times_path,
                    time_window_s=20.0,
                )

        self.param_names: List[str] = [
            "a1 (km)",
            "a2 (km)",
            "theta (x pi)",
            "np (frac)",
            "tp (x 2pi)",
            "dmax (m)",
            "vr (km/s)",
        ]

        n_params = len(self.cfg.inversion_params.parameters)
        self.param_ranges = np.zeros((n_params, 2), dtype=np.float64)
        for i, param in enumerate(self.cfg.inversion_params.parameters):
            self.param_ranges[i, 0] = param.min_val
            self.param_ranges[i, 1] = param.max_val

        # Runtime state — typed as Any so subclasses can store their own config objects
        # without creating circular imports between NA and MCMC modules.
        self._na_cfg_runtime: Optional[Any] = None
        self._mcmc_cfg_runtime: Optional[Any] = None
        self._mcmc_step_index: Optional[int] = None
        self._pymc_inner: bool = False
        self._eval_count: int = 0
        self._best_misfit_seen: float = np.inf
        self._axitra_id_counter: int = 0

        # Best synthetic seismograms (nsta, 3, npts) — updated inside objective_function
        self.best_synthetics: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    def _next_axitra_id(self) -> int:
        """Generate a unique-enough AXITRA ID to avoid temporary file collisions.
        (Genera un ID único para AXITRA para evitar colisiones de archivos temporales.)
        """
        self._axitra_id_counter += 1
        ns_stamp = int(time.time_ns() % 1_000_000_000)
        pid_part = (os.getpid() % 10_000) * 100_000
        ctr_part = self._axitra_id_counter % 100_000
        return ns_stamp + pid_part + ctr_part

    # ------------------------------------------------------------------
    def _build_geometry_from_parameters(self, model: np.ndarray):
        """Build fault geometry with ellipse slip from a model parameter vector.
        (Construye la geometría de falla con deslizamiento elíptico desde un vector de parámetros.)

        Reuses the precomputed invariant mesh stored in ``self.base_geometry``
        and applies model-dependent fields on a deep copy.
        """
        geom = deepcopy(self.base_geometry)
        return self.fm.apply_ellipse_model_to_geometry(geometry=geom, model=model)

    # ------------------------------------------------------------------
    def objective_function(self, model: np.ndarray) -> float:
        """Evaluate the forward model and return L2 misfit for one parameter vector.
        (Evalúa el modelo forward y retorna el desajuste L2 para un vector de parámetros.)

        Side-effects
        ------------
        - Updates ``self._best_misfit_seen`` and ``self.best_synthetics`` when improved.
        - Cleans AXITRA temporary files unless the active config sets ``keep_axitra_files=True``.
        """
        self._eval_count += 1

        # Determine logging tag and iteration estimate
        if self._mcmc_step_index is not None:
            log_tag = "MCMC"
            iter_est = int(self._mcmc_step_index)
        elif self._pymc_inner:
            log_tag = "PYMC"
            iter_est = self._eval_count
        else:
            log_tag = "NA"
            iter_est = 0
            if self._na_cfg_runtime is not None:
                n0 = int(self._na_cfg_runtime.n_samples_initial)
                ns = max(1, int(self._na_cfg_runtime.n_samples_iteration))
                if self._eval_count > n0:
                    iter_est = 1 + ((self._eval_count - n0 - 1) // ns)

        ap = None
        try:
            if self.misfit_calc is None:
                return 1e10

            geom = self._build_geometry_from_parameters(model)

            ap = self.fm.build_axitra(geom, latlon=False, freesurface=True)
            ap = self.fm.green(ap, quiet=True)
            _, sx, sy, sz = self.fm.conv(
                ap, geom, source_type=1, t0=float(self.cfg.ellipse.t0), quiet=True
            )

            synthetics = np.array([sx, sy, sz])
            synthetics = np.transpose(synthetics, (1, 2, 0))
            synthetics = np.transpose(synthetics, (0, 2, 1))

            # Filter synthetics to the same frequency band as observed data
            try:
                from signal_utils import bandpass_filter_waveforms
            except ImportError:
                from .signal_utils import bandpass_filter_waveforms

            synthetics = bandpass_filter_waveforms(
                synthetics,
                self.time_array,
                freq1=float(self.cfg.ellipse.freq1),
                freq2=float(self.cfg.ellipse.freq2),
            )

            misfit = float(self.misfit_calc.l2_misfit(synthetics))
            if misfit < self._best_misfit_seen:
                self._best_misfit_seen = misfit
                self.best_synthetics = synthetics.copy()

            # Optional slip/moment scaling diagnostics
            diag = getattr(geom, "slip_scale_diagnostics", None)
            diag_s = ""
            if isinstance(diag, dict):
                m0_t = diag.get("m0_target_nm")
                m0_t_s = f"{float(m0_t):.6e}" if m0_t is not None else "n/a"
                diag_s = (
                    f" mt_mode={diag.get('mt_scaling_mode')} "
                    f"dmax_req={float(diag.get('dmax_requested_m', 0.0)):.6e} "
                    f"dmax_eff={float(diag.get('dmax_effective_m', 0.0)):.6e} "
                    f"M0_target={m0_t_s} "
                    f"M0_L1={float(diag.get('m0_L1_sum_abs_moments_nm', 0.0)):.6e}"
                )
                warn = diag.get("strict_scale_warning")
                if warn:
                    diag_s += f" WARN={warn}"

            print(
                f"[{log_tag}] iter={iter_est:05d} eval={self._eval_count:05d} "
                f"misfit={misfit:.6e} best={self._best_misfit_seen:.6e}"
                f"{diag_s}",
                flush=True,
            )
            return misfit

        except Exception as exc:
            print(
                f"[{log_tag}] iter={iter_est:05d} eval={self._eval_count:05d} "
                f"misfit=1.000000e+10 best={self._best_misfit_seen:.6e}",
                flush=True,
            )
            logger.error("Model evaluation failed: %s", exc)
            return 1e10

        finally:
            keep = False
            if self._na_cfg_runtime is not None:
                keep = bool(self._na_cfg_runtime.keep_axitra_files)
            elif self._mcmc_cfg_runtime is not None:
                keep = bool(self._mcmc_cfg_runtime.keep_axitra_files)
            if ap is not None and not keep:
                try:
                    ap.clean()
                except Exception:
                    pass


__all__ = [
    "NAModel",
    "MisfitCalculator",
    "NAResult",
    "BaseInversionModel",
]
