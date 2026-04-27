"""
NA inversion driver using neighpy as search backend.

This module replaces the local prototype search implementation and keeps a
compatible API for notebook usage:

- NAConfig
- NAInversionModel.run_na_search(...)
- result.best_model / result.all_models / result.export_results(...)
"""

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence
import csv
from copy import deepcopy
import json
import logging
import os
import time

import numpy as np

# Handle both relative imports (when used as package) and direct imports (from notebooks)
try:
	from .config_parser import ConfigParser
	from .forward_model import AxitraForwardModel
except ImportError:
	from config_parser import ConfigParser
	from forward_model import AxitraForwardModel

logger = logging.getLogger(__name__)


@dataclass
class NAConfig:
	"""Search configuration compatible with existing code."""

	n_samples_initial: int = 30
	n_samples_iteration: int = 30
	n_iterations: int = 10
	n_cells_resample: int = 7
	random_seed: Optional[int] = None
	keep_axitra_files: bool = False


@dataclass
class NAModel:
	"""Single sampled model and its objective value."""

	model: np.ndarray
	misfit: float
	iteration: int


class MisfitCalculator:
	"""Waveform misfit calculator."""

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
				raise FileNotFoundError(f"Missing required old-misfit file: {azi_times_path}")
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

	def diagnostics_summary(self, synthetic: np.ndarray, max_stations: int = 3) -> str:
		"""Return compact diagnostics to explain misfit behavior for one model."""
		if synthetic.shape != self.observed.shape:
			return (
				f"[MISFIT DIAG] shape mismatch: synthetic {synthetic.shape} "
				f"vs observed {self.observed.shape}"
			)

		nsta, _, npts = self.observed.shape
		if len(self.time) > 1:
			dt = float(self.time[1] - self.time[0])
		else:
			dt = 1.0
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
			end_p = start_p + win
			start_s = int(np.rint(float(self.ts[j]))) * sampling
			end_s = start_s + win

			kp0 = max(0, start_p - 1)
			kp1 = min(npts - 1, end_p - 1)
			ks0 = max(0, start_s - 1)
			ks1 = min(npts - 1, end_s - 1)

			# P window: radial + vertical
			r_obs_rms = 0.0
			r_syn_rms = 0.0
			z_obs_rms = 0.0
			z_syn_rms = 0.0
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

			# S window: transverse
			t_obs_rms = 0.0
			t_syn_rms = 0.0
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

	def l2_misfit(self, synthetic: np.ndarray) -> float:
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

		if len(self.time) > 1:
			dt = float(self.time[1] - self.time[0])
		else:
			dt = 1.0
		sampling = max(1, int(np.rint(1.0 / dt)))
		win = max(1, int(np.rint(self.time_window_s * sampling)))

		num = 0.0
		den = 0.0

		for j in range(nsta):
			az = float(self.azi[j])

			# Old Fortran behavior:
			# start_data = nint(tP(j)) * sampling
			# start_data = nint(tS(j)) * sampling
			start_p = int(np.rint(float(self.tp[j]))) * sampling
			end_p = start_p + win
			start_s = int(np.rint(float(self.ts[j]))) * sampling
			end_s = start_s + win

			# Convert Fortran 1-based inclusive indices to Python slices.
			kp0 = max(0, start_p - 1)
			kp1 = min(npts - 1, end_p - 1)
			ks0 = max(0, start_s - 1)
			ks1 = min(npts - 1, end_s - 1)

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

		if den > 0.0:
			return num / den
		return num


class NAResult:
	"""Compatibility wrapper for neighpy outputs."""

	def __init__(self, all_models: List[NAModel], param_names: Optional[Sequence[str]] = None):
		self.all_models = all_models
		self.best_model = min(all_models, key=lambda m: m.misfit) if all_models else None
		self.param_names = list(param_names) if param_names is not None else [
			"a1",
			"a2",
			"theta",
			"np",
			"tp",
			"dmax",
			"vr",
		]

	def export_results(self, filepath: Path) -> None:
		payload = {
			"metadata": {
				"timestamp": datetime.now().isoformat(),
				"n_models": len(self.all_models),
				"best_misfit": self.best_model.misfit if self.best_model else None,
				"param_names": list(self.param_names),
			},
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
		filepath = Path(filepath)
		filepath.parent.mkdir(parents=True, exist_ok=True)

		fieldnames = ["iteration", "misfit", *self.param_names]
		with filepath.open("w", newline="", encoding="utf-8") as f:
			writer = csv.DictWriter(f, fieldnames=fieldnames)
			writer.writeheader()
			for model in self.all_models:
				row = {
					"iteration": int(model.iteration),
					"misfit": float(model.misfit),
				}
				for name, value in zip(self.param_names, model.model):
					row[name] = float(value)
				writer.writerow(row)


class NAInversionModel:
	"""
	Kinematic Inversion Model using the Neighbourhood Algorithm (NA).
	(Modelo de Inversión Cinemática utilizando el algoritmo Neighbourhood Algorithm.)

	Integrates observed data, arrival times, and event configuration to evaluate 
	different kinematic rupture models. It communicates with AXITRA to simulate 
	synthetic seismograms and calculate the misfit against real data.
	(Integra los datos observados, tiempos de llegada y configuración del evento 
	para evaluar distintos modelos de ruptura. Se comunica con AXITRA para simular 
	sismogramas sintéticos y calcular el error respecto a los datos reales.)

	Attributes:
		input_ctl_path (str): Path to the configuration file.
		                      (Ruta al archivo de configuración.)
		axitra_dir (str): Path to the AXITRA simulation binaries directory.
		                  (Ruta al directorio de binarios de AXITRA.)
		observed_waveforms (np.ndarray): Processed and filtered real seismograms.
		                                 (Sismogramas reales procesados y filtrados.)
		time_array (np.ndarray): Time vector corresponding to the waveforms.
		                         (Vector de tiempo correspondiente a las formas de onda.)
		azi_times_array (np.ndarray): Theoretical arrival times (P/S) per station.
		                              (Tiempos de llegada teóricos P/S por estación.)
		param_names (list): Names of the parameters being inverted.
		                    (Nombres de los parámetros que se están invirtiendo.)
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

		self.param_names = [
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

		self._na_cfg_runtime: Optional[NAConfig] = None
		self._eval_count: int = 0
		self._best_misfit_seen: float = np.inf
		self._axitra_id_counter: int = 0

	def _next_axitra_id(self) -> int:
		"""Generate a large unique-enough AXITRA ID to avoid file collisions."""
		self._axitra_id_counter += 1
		ns_stamp = int(time.time_ns() % 1_000_000_000)
		pid_part = (os.getpid() % 10_000) * 100_000
		ctr_part = self._axitra_id_counter % 100_000
		return ns_stamp + pid_part + ctr_part

	def _build_geometry_from_parameters(self, model: np.ndarray):
		"""Build fault geometry with ellipse-to-subfault slip mapping.

		Reuses a precomputed invariant mesh and applies model-dependent fields.
		"""
		geom = deepcopy(self.base_geometry)
		return self.fm.apply_ellipse_model_to_geometry(geometry=geom, model=model)

	def objective_function(self, model: np.ndarray) -> float:
		self._eval_count += 1
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

			ap = self.fm.build_axitra(
				geom,
				latlon=False,
				freesurface=True,
			)

			ap = self.fm.green(ap, quiet=True)

			_, sx, sy, sz = self.fm.conv(ap, geom, source_type=1, t0=float(self.cfg.ellipse.t0), quiet=True)

			synthetics = np.array([sx, sy, sz])
			synthetics = np.transpose(synthetics, (1, 2, 0))
			synthetics = np.transpose(synthetics, (0, 2, 1))
			misfit = float(self.misfit_calc.l2_misfit(synthetics))
			if misfit < self._best_misfit_seen:
				self._best_misfit_seen = misfit
			print(
				f"[NA] iter={iter_est:03d} eval={self._eval_count:05d} "
				f"misfit={misfit:.6e} best={self._best_misfit_seen:.6e}"
				,
				flush=True,
			)
			return misfit
		except Exception as exc:
			print(
				f"[NA] iter={iter_est:03d} eval={self._eval_count:05d} "
				f"misfit=1.000000e+10 best={self._best_misfit_seen:.6e}"
				,
				flush=True,
			)
			logger.error("Model evaluation failed: %s", exc)
			return 1e10
		finally:
			if ap is not None and (self._na_cfg_runtime is None or not self._na_cfg_runtime.keep_axitra_files):
				try:
					ap.clean()
				except Exception:
					pass

	def run_na_search(self, na_config: Optional[NAConfig] = None) -> NAResult:
		if na_config is None:
			na_config = NAConfig(
				n_samples_initial=self.cfg.inversion_process.ss1,
				n_samples_iteration=self.cfg.inversion_process.ss_other,
				n_iterations=self.cfg.inversion_process.num_iterations,
				n_cells_resample=self.cfg.inversion_process.cells_resample,
			)

		# Workaround for neighpy sampling bug when ns % nr != 0.
		# In neighpy.search.NASearcher, internal allocation uses ns // nr, but the
		# best Voronoi cell walks ns % nr extra steps, which can overflow.
		nr = max(1, int(na_config.n_cells_resample))
		ns_requested = int(na_config.n_samples_iteration)
		ns_effective = max(nr, ns_requested)
		if ns_effective % nr != 0:
			ns_effective = ((ns_effective // nr) + 1) * nr
			print(
				f"[NA] Adjusting n_samples_iteration from {ns_requested} to {ns_effective} "
				f"to satisfy neighpy constraint (multiple of n_cells_resample={nr})."
			)

		if na_config.random_seed is not None:
			np.random.seed(na_config.random_seed)

		self._na_cfg_runtime = na_config
		if ns_effective != int(na_config.n_samples_iteration):
			self._na_cfg_runtime = NAConfig(
				n_samples_initial=int(na_config.n_samples_initial),
				n_samples_iteration=int(ns_effective),
				n_iterations=int(na_config.n_iterations),
				n_cells_resample=int(nr),
				random_seed=na_config.random_seed,
				keep_axitra_files=bool(na_config.keep_axitra_files),
			)
		self._eval_count = 0
		self._best_misfit_seen = np.inf
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
		expected_models = int(na_config.n_samples_initial) + int(ns_effective) * int(na_config.n_iterations)
		print(
			f"[NA] Starting search: ni={int(na_config.n_samples_initial)}, ns={int(ns_effective)}, "
			f"n={int(na_config.n_iterations)}, nr={int(nr)} -> expected evaluations={expected_models}",
			flush=True,
		)
		searcher.run()

		samples = np.asarray(searcher.samples)
		objectives = np.asarray(searcher.objectives)

		# Approximate iteration labels for compatibility with existing notebooks.
		n0 = int(na_config.n_samples_initial)
		niter = ns_effective
		models: List[NAModel] = []
		for i, (sample, obj) in enumerate(zip(samples, objectives)):
			if i < n0:
				iteration = 0
			else:
				iteration = 1 + ((i - n0) // max(1, niter))
			models.append(NAModel(model=np.asarray(sample, dtype=float), misfit=float(obj), iteration=iteration))

		return NAResult(models, param_names=self.param_names)
