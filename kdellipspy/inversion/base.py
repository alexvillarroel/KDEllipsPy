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
    from ..core.config_parser import ConfigParser
    from ..core.forward_model import AxitraForwardModel
except ImportError:
    try:
        from kdellipspy.core.config_parser import ConfigParser
        from kdellipspy.core.forward_model import AxitraForwardModel
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
        station_flags: Optional[np.ndarray] = None,
    ):
        self.observed = observed_waveforms
        self.time = time_array
        self.time_window_s = float(time_window_s)
        
        nsta = self.observed.shape[0]
        if station_flags is not None:
            self.station_flags = np.asarray(station_flags, dtype=bool)
        else:
            self.station_flags = np.ones((nsta, 3), dtype=bool)

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
            use_n, use_e, use_z = self.station_flags[j]
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
                
                # Apply N/Z flags to Radial/Vertical
                if use_n:
                    r_obs = x_obs * np.cos(az) + y_obs * np.sin(az)
                    r_syn = x_syn * np.cos(az) + y_syn * np.sin(az)
                    r_obs_rms = float(np.sqrt(np.mean(r_obs ** 2)))
                    r_syn_rms = float(np.sqrt(np.mean(r_syn ** 2)))
                    obs_energy += float(np.sum(r_obs ** 2))
                    syn_energy += float(np.sum(r_syn ** 2))
                if use_z:
                    z_obs_rms = float(np.sqrt(np.mean(z_obs ** 2)))
                    z_syn_rms = float(np.sqrt(np.mean(z_syn ** 2)))
                    obs_energy += float(np.sum(z_obs ** 2))
                    syn_energy += float(np.sum(z_syn ** 2))

            t_obs_rms = t_syn_rms = 0.0
            if ks1 >= ks0:
                x_obs = self.observed[j, 0, ks0:ks1 + 1]
                y_obs = self.observed[j, 1, ks0:ks1 + 1]
                x_syn = synthetic[j, 0, ks0:ks1 + 1]
                y_syn = synthetic[j, 1, ks0:ks1 + 1]
                
                # Apply E flag to Transverse
                if use_e:
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
    def l2_misfit(self, synthetic: np.ndarray, use_full_signal: bool = False) -> float:
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
        if use_full_signal == False:
            for j in range(nsta):
                use_n, use_e, use_z = self.station_flags[j]
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
                    
                    if use_n:  # Map N flag to Radial
                        r_obs = x_obs * np.cos(az) + y_obs * np.sin(az)
                        r_syn = x_syn * np.cos(az) + y_syn * np.sin(az)
                        num += float(np.sum((r_obs - r_syn) ** 2))
                        den += float(np.sum(r_obs ** 2))
                    if use_z:
                        num += float(np.sum((z_obs - z_syn) ** 2))
                        den += float(np.sum(z_obs ** 2))

                # S window: transverse
                if ks1 >= ks0:
                    x_obs = self.observed[j, 0, ks0:ks1 + 1]
                    y_obs = self.observed[j, 1, ks0:ks1 + 1]
                    x_syn = synthetic[j, 0, ks0:ks1 + 1]
                    y_syn = synthetic[j, 1, ks0:ks1 + 1]
                    
                    if use_e:  # Map E flag to Transverse
                        t_obs = y_obs * np.cos(az) - x_obs * np.sin(az)
                        t_syn = y_syn * np.cos(az) - x_syn * np.sin(az)
                        num += float(np.sum((t_obs - t_syn) ** 2))
                        den += float(np.sum(t_obs ** 2))
        else:
            for j in range(nsta):
                use_n, use_e, use_z = self.station_flags[j]
                az = float(self.azi[j])
                x_obs = self.observed[j, 0, :]
                y_obs = self.observed[j, 1, :]
                z_obs = self.observed[j, 2, :]
                x_syn = synthetic[j, 0, :]
                y_syn = synthetic[j, 1, :]
                z_syn = synthetic[j, 2, :]
                
                # Rotate the full signal
                r_obs = x_obs * np.cos(az) + y_obs * np.sin(az)
                r_syn = x_syn * np.cos(az) + y_syn * np.sin(az)
                t_obs = y_obs * np.cos(az) - x_obs * np.sin(az)
                t_syn = y_syn * np.cos(az) - x_syn * np.sin(az)
                
                if use_n:  # Map N flag to Radial
                    num += float(np.sum((r_obs - r_syn)**2))
                    den += float(np.sum(r_obs**2))
                if use_e:  # Map E flag to Transverse
                    num += float(np.sum((t_obs - t_syn)**2))
                    den += float(np.sum(t_obs**2))
                if use_z:
                    num += float(np.sum((z_obs - z_syn)**2))
                    den += float(np.sum(z_obs**2))

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
    best_synthetics : Synthetic waveforms of the best model (optional)
    observed       : Observed waveforms used for inversion (optional)
    time           : Time array used for inversion (optional)
    config         : ConfigParser instance (optional)
    """

    _DEFAULT_PARAM_NAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]

    def __init__(
        self,
        all_models: List[NAModel],
        param_names: Optional[Sequence[str]] = None,
        extra_metadata: Optional[Dict[str, Any]] = None,
        best_synthetics: Optional[np.ndarray] = None,
        observed: Optional[np.ndarray] = None,
        time: Optional[np.ndarray] = None,
        config: Optional[ConfigParser] = None,
        azi_times_array: Optional[np.ndarray] = None,
    ):
        self.all_models = all_models
        self.best_model = min(all_models, key=lambda m: m.misfit) if all_models else None
        self.param_names = (
            list(param_names) if param_names is not None else list(self._DEFAULT_PARAM_NAMES)
        )
        self.extra_metadata = dict(extra_metadata) if extra_metadata else {}

        # Extended fields for full persistence
        self.best_synthetics = best_synthetics
        self.observed = observed
        self.time = time
        self.config = config
        # azi_times table (nsta,3) = [azimuth_rad, tP_s, tS_s]; enables rotated
        # R/T/Z waveform plots with phase windows.
        self.azi_times_array = azi_times_array

    def save(self, filepath: str | Path) -> None:
        """
        Save the entire NAResult object to a file for later reloading.
        Uses joblib if available (efficient for numpy arrays), otherwise uses pickle.
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            import joblib
            joblib.dump(self, filepath)
        except ImportError:
            import pickle
            with open(filepath, 'wb') as f:
                pickle.dump(self, f)
        
        print(f"Inversion results saved to {filepath}")

    @classmethod
    def load(cls, filepath: str | Path) -> 'NAResult':
        """
        Load a NAResult object from a saved file.
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Result file not found: {filepath}")
            
        try:
            import joblib
            obj = joblib.load(filepath)
        except ImportError:
            import pickle
            with open(filepath, 'rb') as f:
                obj = pickle.load(f)
        
        if not isinstance(obj, cls):
            raise TypeError(f"Loaded object is not a {cls.__name__}, got {type(obj)}")
            
        return obj

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

    def plot(self, show: bool = True, save_path: Optional[str] = None) -> Tuple[Any, Any]:
        """
        Grafica el resumen de resultados de la búsqueda NA.
        """
        from ..core.plotting import plot_na_results
        
        # Convert all_models to list of dicts for the plotting function
        rows = []
        for model in self.all_models:
            row = {"iteration": float(model.iteration), "misfit": float(model.misfit)}
            for name, value in zip(self.param_names, model.model):
                row[name] = float(value)
            rows.append(row)
            
        return plot_na_results(rows, self.param_names, show=show, save_path=save_path)

    def plot_convergence(self, show: bool = True, save_path: Optional[str] = None) -> Tuple[Any, Any]:
        """
        Grafica la convergencia detallada por cada parámetro.
        """
        from ..core.plotting import plot_parameter_convergence
        
        rows = []
        misfits = []
        iterations = []
        for model in self.all_models:
            row = {"iteration": float(model.iteration), "misfit": float(model.misfit)}
            for name, value in zip(self.param_names, model.model):
                row[name] = float(value)
            rows.append(row)
            misfits.append(model.misfit)
            iterations.append(model.iteration)
            
        return plot_parameter_convergence(
            rows, self.param_names, np.array(misfits), np.array(iterations),
            show=show, save_path=save_path
        )

    def plot_fit(self, show: bool = True, save_path: Optional[str] = None,
                 rotate: bool = True, mark_windows: bool = True) -> Tuple[Any, Any]:
        """
        Grafica el mejor ajuste de formas de onda encontrado.
        Requiere que best_synthetics, observed y time estén presentes en el objeto.

        Por defecto rota a R/T/Z y sombrea la ventana del misfit (P→R,Z ; S→T),
        consistente con cómo se calcula el misfit. Si no hay datos de azimut/
        tiempos disponibles, cae a las componentes N/E/Z sin ventana.
        """
        if self.best_synthetics is None or self.observed is None or self.time is None:
            print("Datos insuficientes para graficar el ajuste (best_synthetics, observed o time son None).")
            return None, None

        from ..core.plotting import plot_waveform_fit

        station_names = None
        station_flags = None
        if self.config is not None and getattr(self.config, "stations", None) is not None:
            station_names = [s.name for s in self.config.stations.stations]
            station_flags = np.array(
                [[s.use_n, s.use_e, s.use_z] for s in self.config.stations.stations], dtype=bool
            )
        else:
            nsta = self.observed.shape[0]
            station_names = [f"STA{i+1}" for i in range(nsta)]

        # azi_times: usar el guardado; si falta (joblibs antiguos), recalcular.
        azt = getattr(self, "azi_times_array", None)
        if azt is None and self.config is not None:
            try:
                from ..core.signal_utils import build_azi_times_array
                azt = build_azi_times_array(config=self.config)
            except Exception:
                azt = None
        azimuths = tp_s = ts_s = None
        if azt is not None:
            azt = np.asarray(azt, dtype=float)
            azimuths, tp_s, ts_s = azt[:, 0], azt[:, 1], azt[:, 2]

        # Ventana del misfit (0.0 = señal completa -> no se marca ventana).
        window_s = None
        try:
            tw = float(self.config.inversion_process.misfit_time_window)
            window_s = tw if tw > 0.0 else None
        except Exception:
            window_s = None

        return plot_waveform_fit(
            observed=self.observed,
            synthetic=self.best_synthetics,
            time=self.time,
            station_names=station_names,
            misfit=self.best_model.misfit if self.best_model else None,
            show=show,
            save_path=save_path,
            azimuths=azimuths,
            tp_s=tp_s,
            ts_s=ts_s,
            window_s=window_s,
            station_flags=station_flags,
            rotate=rotate and azimuths is not None,
            mark_windows=mark_windows and azimuths is not None,
        )

    def plot_ellipse(self, show: bool = True, save_path: Optional[str] = None, title: Optional[str] = None) -> Tuple[Any, Any]:
        """
        Grafica la distribución de slip de la elipse para el mejor misfit.
        Requiere que config y best_model estén presentes en el objeto.
        """
        if self.config is None or self.best_model is None:
            print("Datos insuficientes para graficar la elipse (config o best_model es None).")
            return None, None

        from ..core.forward_model import AxitraForwardModel
        from ..core.plotting import plot_slip_distribution

        fm = AxitraForwardModel.from_config(self.config)
        geom = fm.build_geometry_with_ellipse_slip(self.best_model.model)
        plot_title = title or "2D Slip Distribution (Best misfit)"
        return plot_slip_distribution(
            geom,
            title=plot_title,
            show=show,
            save_path=save_path,
        )

    def plot_elipse(self, show: bool = True, save_path: Optional[str] = None, title: Optional[str] = None) -> Tuple[Any, Any]:
        """Alias por compatibilidad con el nombre en español."""
        return self.plot_ellipse(show=show, save_path=save_path, title=title)

    def plot_azimuthal(self, show: bool = True, save_path: Optional[str] = None) -> Tuple[Any, Any]:
        """Diagrama polar (radar) de la cobertura azimutal de las estaciones
        respecto al epicentro, con el mayor hueco azimutal sombreado."""
        if self.config is None or getattr(self.config, "stations", None) is None:
            print("Sin config/estaciones: no se puede graficar la cobertura azimutal.")
            return None, None
        from ..core.plotting import plot_azimuthal_coverage
        sts = self.config.stations.stations
        sp = self.config.source_position
        return plot_azimuthal_coverage(
            [s.latitude for s in sts], [s.longitude for s in sts],
            [s.name for s in sts], sp.latitude, sp.longitude,
            show=show, save_path=save_path,
        )

    def plot_ellipse_map(self, show: bool = True, save_path: Optional[str] = None,
                         pad: float = 0.25) -> Tuple[Any, Any]:
        """Mapa cartopy de la elipse de slip proyectada a la superficie (footprint
        coloreado por slip + borde + hipocentro + estaciones)."""
        if self.config is None or self.best_model is None:
            print("Sin config/best_model: no se puede mapear la elipse.")
            return None, None
        from ..core.forward_model import AxitraForwardModel
        from ..core.plotting import plot_ellipse_map as _pem
        fm = AxitraForwardModel.from_config(self.config)
        geom = fm.apply_ellipse_model_to_geometry(
            fm.build_geometry(), self.best_model.model, keep_all_sources=True)
        src = geom.to_axitra_sources(latlon=True)        # [idx, lat, lon, z]
        slip = np.array([sp.displacement for sp in geom.source_points], dtype=float)
        sp_cfg = self.config.source_position
        sts = self.config.stations.stations if self.config.stations else []
        return _pem(src[:, 1], src[:, 2], slip, sp_cfg.latitude, sp_cfg.longitude,
                    [s.latitude for s in sts], [s.longitude for s in sts],
                    [s.name for s in sts], pad=pad, show=show, save_path=save_path)

    def _misfit_breakdown_arrays(self, window_s="auto"):
        """Devuelve (num, den, window_s) por estación×(R,T,Z) del mejor modelo,
        con la misma lógica del misfit (P→R,Z desde tP; S→T desde tS; o señal
        completa si window_s es None/0)."""
        obs = np.asarray(self.observed, float)
        syn = np.asarray(self.best_synthetics, float)
        time = np.asarray(self.time, float)
        cfg = self.config
        flags = np.array([[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations], bool)
        azt = getattr(self, "azi_times_array", None)
        if azt is None:
            from ..core.signal_utils import build_azi_times_array
            azt = build_azi_times_array(config=cfg)
        azt = np.asarray(azt, float)
        azi, tP, tS = azt[:, 0], azt[:, 1], azt[:, 2]
        if window_s == "auto":
            try:
                tw = float(cfg.inversion_process.misfit_time_window)
                window_s = tw if tw > 0 else None
            except Exception:
                window_s = None

        nsta, _, npts = obs.shape
        num = np.zeros((nsta, 3)); den = np.zeros((nsta, 3))
        use_win = window_s is not None and window_s > 0
        if use_win:
            dt = float(time[1] - time[0]); samp = max(1, int(np.rint(1.0 / dt)))
            win = max(1, int(np.rint(window_s * samp)))
        for j in range(nsta):
            az = float(azi[j]); un, ue, uz = flags[j]
            if use_win:
                sp = int(np.rint(tP[j])) * samp; ss = int(np.rint(tS[j])) * samp
                kp0, kp1 = max(0, sp - 1), min(npts - 1, sp + win - 1)
                ks0, ks1 = max(0, ss - 1), min(npts - 1, ss + win - 1)
            else:
                kp0, kp1 = 0, npts - 1; ks0, ks1 = 0, npts - 1
            if kp1 >= kp0:
                xo, yo, zo = obs[j, 0, kp0:kp1+1], obs[j, 1, kp0:kp1+1], obs[j, 2, kp0:kp1+1]
                xs, ys, zs = syn[j, 0, kp0:kp1+1], syn[j, 1, kp0:kp1+1], syn[j, 2, kp0:kp1+1]
                if un:
                    ro = xo*np.cos(az)+yo*np.sin(az); rs = xs*np.cos(az)+ys*np.sin(az)
                    num[j, 0] = np.sum((ro-rs)**2); den[j, 0] = np.sum(ro**2)
                if uz:
                    num[j, 2] = np.sum((zo-zs)**2); den[j, 2] = np.sum(zo**2)
            if ks1 >= ks0:
                xo, yo = obs[j, 0, ks0:ks1+1], obs[j, 1, ks0:ks1+1]
                xs, ys = syn[j, 0, ks0:ks1+1], syn[j, 1, ks0:ks1+1]
                if ue:
                    to = yo*np.cos(az)-xo*np.sin(az); ts = ys*np.cos(az)-xs*np.sin(az)
                    num[j, 1] = np.sum((to-ts)**2); den[j, 1] = np.sum(to**2)
        return num, den, window_s

    def plot_misfit_breakdown(self, show: bool = True, save_path: Optional[str] = None,
                              window_s="auto") -> Tuple[Any, Any]:
        """Visualiza la contribución al misfit por estación/componente (heatmaps)."""
        if self.best_synthetics is None or self.observed is None or self.config is None:
            print("Datos insuficientes para el desglose de misfit.")
            return None, None
        num, den, ws = self._misfit_breakdown_arrays(window_s)
        E = num.sum() / den.sum() if den.sum() > 0 else float("nan")
        names = [s.name for s in self.config.stations.stations]
        mode = f"ventana {ws:.0f}s" if ws else "señal completa"
        from ..core.plotting import plot_misfit_contribution
        return plot_misfit_contribution(num, den, names, E, mode,
                                        show=show, save_path=save_path)

    def plot_yolo(self, save_path: str | Path = "dashboard.pdf",
                  show: bool = False, dpi: int = 200) -> Optional[Path]:
        """🎲 YOLO: arma un DASHBOARD de una sola página con todos los paneles
        encajados (GridSpec):

            ┌───────────────── título (evento · misfit · Mw) ─────────────────┐
            │  mapa elipse (cartopy)          │  ajuste R/T/Z (7 est, ventanas)│
            │  azimutal (radar) │ heatmap mf  │  (panel alto)                  │
            │  convergencia de parámetros     │  appraisal (corner)            │
            └─────────────────────────────────────────────────────────────────┘

        Cada panel se renderiza con su método y se compone como imagen, así se
        reaprovecha todo el código existente (cartopy, corner, etc.). El appraisal
        se corre automáticamente si no estaba."""
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        # Corre el appraisal si no estaba. El costo del resampleo escala con
        # n_resample * n_modelos, así que se baja n_resample para ensembles grandes
        # (con 52k modelos, 20000 resamples tardaría muchísimo).
        if getattr(self, "appraisal_samples", None) is None:
            try:
                n_models = len(self.all_models)
                n_res = 20000 if n_models <= 6000 else max(4000, int(20000 * 6000 / n_models))
                self.run_appraisal(n_resample=n_res, verbose=False)
            except Exception as exc:  # noqa: BLE001
                print(f"  [plot_yolo] appraisal no disponible: {exc}")

        def _autocrop(buf):
            """Recorta los márgenes blancos uniformes del panel para que llene mejor."""
            rgb = buf[:, :, :3]
            mask = np.any(rgb < 248, axis=2)
            if not mask.any():
                return buf
            rows = np.where(np.any(mask, axis=1))[0]
            cols = np.where(np.any(mask, axis=0))[0]
            pad = 4  # margen mínimo para no cortar al ras
            r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad, buf.shape[0] - 1)
            c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad, buf.shape[1] - 1)
            return buf[r0:r1 + 1, c0:c1 + 1]

        def _render(maker):
            """Genera una figura componente y la devuelve como array RGBA recortado."""
            try:
                fig, _ = maker()
                if fig is None:
                    return None
                fig.set_dpi(150)
                fig.canvas.draw()
                buf = np.asarray(fig.canvas.buffer_rgba()).copy()
                plt.close(fig)
                return _autocrop(buf)
            except Exception as exc:  # noqa: BLE001
                print(f"  [plot_yolo] panel falló: {exc}")
                return None

        imgs = {
            "map":  _render(lambda: self.plot_ellipse_map(show=False, pad=2.0)),
            "fit":  _render(lambda: self.plot_fit(show=False)),
            "azi":  _render(lambda: self.plot_azimuthal(show=False)),
            "heat": _render(lambda: self.plot_misfit_breakdown(show=False)),
            "conv": _render(lambda: self.plot_convergence(show=False)),
            "appr": (_render(lambda: self.plot_appraisal(show=False))
                     if getattr(self, "appraisal_samples", None) is not None else None),
        }

        # --- Componer el dashboard --------------------------------------------
        fig = plt.figure(figsize=(19, 16))
        gs = GridSpec(3, 3, figure=fig, width_ratios=[1, 1, 1.7],
                      height_ratios=[1.1, 1.0, 1.0], hspace=0.04, wspace=0.03,
                      left=0.01, right=0.99, top=0.955, bottom=0.01)
        axA = fig.add_subplot(gs[0, 0:2])   # mapa elipse
        axC = fig.add_subplot(gs[1, 0])     # azimutal
        axD = fig.add_subplot(gs[1, 1])     # heatmap misfit
        axB = fig.add_subplot(gs[0:2, 2])   # ajuste R/T/Z (alto)
        axE = fig.add_subplot(gs[2, 0:2])   # convergencia
        axF = fig.add_subplot(gs[2, 2])     # appraisal

        for ax, key in [(axA, "map"), (axB, "fit"), (axC, "azi"),
                        (axD, "heat"), (axE, "conv"), (axF, "appr")]:
            ax.axis("off")
            if imgs[key] is not None:
                ax.imshow(imgs[key], aspect="auto")

        # Título: evento, misfit, Mw.
        title = "Resultados de la inversión"
        try:
            sp = self.config.source_position
            mw_txt = ""
            try:
                from ..core.forward_model import AxitraForwardModel
                fm = AxitraForwardModel.from_config(self.config)
                _m0, mw = fm.estimate_total_moment_and_mw(self.best_model.model)
                mw_txt = f"  ·  Mw {mw:.2f}"
            except Exception:
                pass
            mf = self.best_model.misfit if self.best_model else float("nan")
            title = (f"{getattr(sp, 'event_name', 'Evento')}  ·  "
                     f"{sp.latitude:.2f}, {sp.longitude:.2f}  ·  {sp.depth:.0f} km  ·  "
                     f"misfit {mf:.4f}{mw_txt}")
        except Exception:
            pass
        fig.suptitle(title, fontsize=16, fontweight="bold", y=0.985)

        fig.savefig(str(save_path), dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)
        n_ok = sum(1 for v in imgs.values() if v is not None)
        print(f"🎲 plot_yolo dashboard: {n_ok} paneles → {save_path}")
        return save_path

    # ------------------------------------------------------------------
    # Uncertainty appraisal (NA second stage, Sambridge 1999 Part II)
    # ------------------------------------------------------------------
    def _get_bounds(self) -> np.ndarray:
        """Return (n_params, 2) [min, max] per parameter.

        Taken from the stored config (the inversion bounds); if no config is
        available, falls back to the min/max spanned by the sampled models.
        (Devuelve los límites por parámetro desde el config; si no hay config,
        usa el rango cubierto por los modelos muestreados.)
        """
        cfg = self.config
        if cfg is not None and getattr(cfg, "inversion_params", None) is not None:
            params = cfg.inversion_params.parameters
            return np.array([[p.min_val, p.max_val] for p in params], dtype=float)
        M = np.array([m.model for m in self.all_models], dtype=float)
        return np.column_stack([M.min(axis=0), M.max(axis=0)])

    def run_appraisal(
        self,
        n_resample: int = 20000,
        n_walkers: int = 1,
        temperature: Optional[float] = None,
        bounds: Optional[np.ndarray] = None,
        seed: Optional[int] = None,
        verbose: bool = True,
        save: bool = True,
    ) -> Optional[np.ndarray]:
        """Run the NA appraisal stage to approximate the posterior.
        (Ejecuta la etapa de 'appraisal' del NA para aproximar la posterior.)

        Reuses the models already evaluated during the NA search (``all_models``)
        and resamples their Voronoi cells with ``neighpy.NAAppraiser`` — it does
        NOT evaluate the forward model again. Populates ``appraisal_samples``,
        ``appraisal_mean`` and ``appraisal_covariance``.

        Parameters
        ----------
        n_resample  : length of the resampling random walk (more = smoother PDFs).
        n_walkers   : parallel walkers (>1 uses neighpy's multiprocessing).
        temperature : misfit-to-log-posterior scale, ``log_ppd = -misfit / T``.
            If ``None``, ``T = 2 * best_misfit`` (auto-scale). SMALLER T → narrower
            posterior (more trust in the data).
            NOTE: the misfit here is the *normalized L2* misfit, not a noise-calibrated
            chi-square, so ``temperature`` implicitly sets the assumed noise level.
            Tune it (or pass an explicit value) if you need calibrated uncertainties.
        bounds      : (n_params, 2) optional; derived from the config if ``None``.
        seed, verbose, save : forwarded to ``NAAppraiser`` / ``NAAppraiser.run``.

        Returns
        -------
        np.ndarray | None
            Posterior ensemble, shape (n_samples, n_params), or ``None`` if
            ``save=False``.
        """
        try:
            from neighpy import NAAppraiser
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "neighpy is required for the appraisal stage. Install with: "
                "pip install neighpy"
            ) from exc

        if not self.all_models:
            raise ValueError("No sampled models available for appraisal.")

        ensemble = np.array([m.model for m in self.all_models], dtype=float)
        misfits = np.array([m.misfit for m in self.all_models], dtype=float)

        # Drop non-finite misfits (failed forward evaluations) so log_ppd is finite.
        finite = np.isfinite(misfits)
        if not finite.all():
            ensemble = ensemble[finite]
            misfits = misfits[finite]
        if ensemble.shape[0] < 2:
            raise ValueError("Not enough finite models for appraisal.")

        best = float(np.min(misfits))
        T = float(temperature) if temperature is not None else max(2.0 * best, 1e-9)
        log_ppd = -misfits / T  # unnormalized log posterior (higher = better fit)

        bnds = self._get_bounds() if bounds is None else np.asarray(bounds, dtype=float)
        bounds_tuple = tuple((float(lo), float(hi)) for lo, hi in bnds)

        appraiser = NAAppraiser(
            n_resample=int(n_resample),
            n_walkers=int(n_walkers),
            initial_ensemble=ensemble,
            log_ppd=log_ppd,
            bounds=bounds_tuple,
            verbose=verbose,
            seed=seed,
        )
        appraiser.run(save=save)

        self.appraisal_samples = getattr(appraiser, "samples", None)
        self.appraisal_mean = getattr(appraiser, "mean", None)
        self.appraisal_covariance = getattr(appraiser, "covariance", None)
        self.appraisal_temperature = T
        return self.appraisal_samples

    def plot_appraisal(
        self,
        samples: Optional[np.ndarray] = None,
        show: bool = True,
        save_path: Optional[str] = None,
        bins: int = 40,
        title: Optional[str] = None,
        **run_kwargs: Any,
    ) -> Tuple[Any, Any]:
        """Corner plot of the posterior uncertainty (1D marginals + 2D trade-offs).
        (Gráfico 'corner' de la incertidumbre: marginales 1D y trade-offs 2D.)

        If ``samples`` is not given and the appraisal has not been run yet, it is
        executed via :meth:`run_appraisal` (extra kwargs are forwarded there, e.g.
        ``n_resample``, ``temperature``, ``seed``).
        """
        if samples is None:
            samples = getattr(self, "appraisal_samples", None)
        if samples is None:
            samples = self.run_appraisal(**run_kwargs)
        if samples is None:
            print("Appraisal did not store samples (save=False?); nothing to plot.")
            return None, None

        from ..core.plotting import plot_uncertainty_corner

        truths = self.best_model.model if self.best_model is not None else None
        mean = getattr(self, "appraisal_mean", None)
        bounds = self._get_bounds()
        if title is None:
            T = getattr(self, "appraisal_temperature", None)
            title = "Posterior uncertainty — NA appraisal"
            if T is not None:
                title += f"  (T={T:.3g}, N={len(samples)})"
        return plot_uncertainty_corner(
            samples,
            self.param_names,
            bounds=bounds,
            truths=truths,
            mean=mean,
            bins=bins,
            show=show,
            save_path=save_path,
            title=title,
        )


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
    axitra_dir         : Path to axitra binary directory (optional)
    observed_waveforms : Observed 3-component seismograms, shape (nsta, 3, npts)
    time_array         : Time vector, shape (npts,)
    azi_times_array    : Pre-computed P/S arrival time table, shape (nsta, 3)

    Attributes
    ----------
    best_synthetics : np.ndarray | None  — Updated whenever a new best misfit is found
    """

    def __init__(
        self,
        input_ctl_path: Optional[str | Path] = None,
        axitra_dir: Optional[str | Path] = None,
        observed_waveforms: Optional[np.ndarray] = None,
        time_array: Optional[np.ndarray] = None,
        azi_times_array: Optional[np.ndarray] = None,
        axitra_aw: float = 0.5,
        axitra_ikmax: int = 100000,
        config: Optional[ConfigParser] = None,
    ):
        if config is not None:
            self.cfg = config
            if hasattr(config, "filepath") and config.filepath not in ("<manual>", "<from_dict>"):
                self.input_ctl_path = Path(config.filepath).resolve()
            else:
                self.input_ctl_path = Path.cwd() / "input.ctl" # Fallback dummy path
            
            self.fm = AxitraForwardModel.from_config(self.cfg, axitra_dir=axitra_dir)
        elif input_ctl_path is not None:
            self.input_ctl_path = Path(input_ctl_path).resolve()
            self.cfg = ConfigParser(str(self.input_ctl_path))
            self.fm = AxitraForwardModel(str(self.input_ctl_path), axitra_dir=axitra_dir)
        else:
            raise ValueError("Either input_ctl_path or config must be provided.")

        self.base_geometry = self.fm.build_geometry()

        self.observed_waveforms = observed_waveforms
        self.time_array = time_array
        self.misfit_calc: Optional[MisfitCalculator] = None
        self.azi_times_array: Optional[np.ndarray] = None
        self.use_full_signal: bool = False
        self.axitra_aw = float(axitra_aw)
        self.axitra_ikmax = int(axitra_ikmax)

        if observed_waveforms is not None and time_array is not None:
            # Validate observed waveform array consistency
            if observed_waveforms.ndim != 3:
                raise ValueError(
                    f"observed_waveforms must be 3D (nsta, 3, npts), got shape {observed_waveforms.shape}"
                )
            nsta_obs, ncomp, npts = observed_waveforms.shape
            if ncomp != 3:
                raise ValueError(
                    f"observed_waveforms must have 3 components (N,E,Z), got {ncomp}"
                )
            nsta_cfg = len(self.cfg.stations.stations) if self.cfg.stations is not None else 0
            if nsta_obs != nsta_cfg:
                raise ValueError(
                    f"Mismatch: observed_waveforms has {nsta_obs} stations "
                    f"but input.ctl defines {nsta_cfg} stations. "
                    f"Ensure you load data for exactly the stations in Section 8."
                )
            if len(time_array) != npts:
                raise ValueError(
                    f"time_array length ({len(time_array)}) must match npts ({npts})"
                )

            # Prepare station flags from ConfigParser
            station_flags = None
            if hasattr(self.cfg, "stations") and self.cfg.stations is not None:
                station_flags = np.array(
                    [[s.use_n, s.use_e, s.use_z] for s in self.cfg.stations.stations],
                    dtype=bool,
                )

            time_window_s = 20.0
            if (
                hasattr(self.cfg, "inversion_process")
                and self.cfg.inversion_process is not None
            ):
                tw = self.cfg.inversion_process.misfit_time_window
                if tw > 0.0:
                    time_window_s = tw
                elif tw == 0.0:
                    self.use_full_signal = True

            # If no azi_times_array is provided, try to calculate it directly
            if azi_times_array is None:
                try:
                    from kdellipspy.core.signal_utils import build_azi_times_array
                    azi_times_array = build_azi_times_array(config=self.cfg)
                    logger.info("Automatically calculated azimuth and travel times for misfit calculation.")
                except Exception as e:
                    logger.warning(f"Could not calculate azi_times_array automatically: {e}")
                    # Fallback to file search if automatic calculation fails
                    if input_ctl_path is not None or (hasattr(config, "filepath") and config.filepath not in ("<manual>", "<from_dict>")):
                        azi_times_path = self.input_ctl_path.parent / "Event" / "azi_times.txt"
                        if azi_times_path.exists():
                            azi_times_array = np.loadtxt(str(azi_times_path), dtype=float)
                            logger.info(f"Loaded azi_times from fallback file: {azi_times_path}")

            if azi_times_array is not None:
                self.azi_times_array = np.asarray(azi_times_array, dtype=float)
                self.misfit_calc = MisfitCalculator(
                    observed_waveforms,
                    time_array,
                    azi_times_array=azi_times_array,
                    time_window_s=time_window_s,
                    station_flags=station_flags,
                )
            else:
                 logger.warning("No azimuth/travel time data available. Misfit calculation will likely fail.")

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
        self._best_model_vec: Optional[np.ndarray] = None
        self._axitra_id_counter: int = 0
        # Si se fija, escribe el mejor modelo actual a este .txt cada vez que mejora
        # (checkpoint en vivo durante la búsqueda).
        self.checkpoint_path: Optional[str | Path] = None

        # Green's-function cache: when ``use_green_cache`` is True the Green's
        # functions for the FULL fixed subfault mesh are computed once and reused
        # in every model evaluation (only the cheap ``conv`` runs per model).
        # Valid because Green's functions depend only on source/station positions,
        # the velocity model and frequencies — not on slip, rupture time or
        # mechanism (which ``conv`` applies from the per-model source history).
        self.use_green_cache: bool = False
        self._green_cache_ap: Optional[Any] = None

        # Best synthetic seismograms (nsta, 3, npts) — updated inside objective_function
        self.best_synthetics: Optional[np.ndarray] = None

    @classmethod
    def from_config(
        cls,
        config: ConfigParser,
        axitra_dir: Optional[str | Path] = None,
        observed_waveforms: Optional[np.ndarray] = None,
        time_array: Optional[np.ndarray] = None,
        azi_times_array: Optional[np.ndarray] = None,
        **kwargs
    ) -> BaseInversionModel:
        """Create an inversion model directly from a ConfigParser object."""
        return cls(
            config=config,
            axitra_dir=axitra_dir,
            observed_waveforms=observed_waveforms,
            time_array=time_array,
            azi_times_array=azi_times_array,
            **kwargs
        )

    @classmethod
    def from_params(
        cls,
        params: Dict[str, Any],
        axitra_dir: Optional[str | Path] = None,
        observed_waveforms: Optional[np.ndarray] = None,
        time_array: Optional[np.ndarray] = None,
        azi_times_array: Optional[np.ndarray] = None,
        **kwargs
    ) -> BaseInversionModel:
        """Create an inversion model from a dictionary of parameters."""
        config = ConfigParser.from_dict(params)
        return cls.from_config(
            config=config,
            axitra_dir=axitra_dir,
            observed_waveforms=observed_waveforms,
            time_array=time_array,
            azi_times_array=azi_times_array,
            **kwargs
        )

    # ------------------------------------------------------------------
    def _next_axitra_id(self) -> int:
        """Generate a unique-enough axitra ID to avoid temporary file collisions.
        (Genera un ID único para axitra para evitar colisiones de archivos temporales.)
        """
        self._axitra_id_counter += 1
        ns_stamp = int(time.time_ns() % 1_000_000_000)
        pid_part = (os.getpid() % 10_000) * 100_000
        ctr_part = self._axitra_id_counter % 100_000
        return ns_stamp + pid_part + ctr_part

    # ------------------------------------------------------------------
    def _build_geometry_from_parameters(self, model: np.ndarray, keep_all_sources: bool = False):
        """Build fault geometry with ellipse slip from a model parameter vector.
        (Construye la geometría de falla con deslizamiento elíptico desde un vector de parámetros.)

        Reuses the precomputed invariant mesh stored in ``self.base_geometry``
        and applies model-dependent fields on a deep copy.

        ``keep_all_sources=True`` keeps the full source set (zero slip outside the
        ellipse) so source positions/indices stay constant across models, enabling
        the cached-Green's-function path.
        """
        geom = deepcopy(self.base_geometry)
        return self.fm.apply_ellipse_model_to_geometry(
            geometry=geom, model=model, keep_all_sources=keep_all_sources
        )

    # ------------------------------------------------------------------
    def _ensure_green_cache(self):
        """Compute (once) and return the cached Green's functions for the full mesh.
        (Calcula una vez y retorna las funciones de Green cacheadas de la malla completa.)

        The source set is the full fixed mesh (``keep_all_sources=True``), so the
        same Green's functions serve every model — only ``conv`` runs per model.
        The returned axitra instance is NOT cleaned during the run.
        """
        if self._green_cache_ap is None:
            # The ellipse parameters do not change WHICH sources exist when
            # keep_all_sources=True (always the full mesh), so any valid model works.
            dummy = np.mean(self.param_ranges, axis=1)
            geom_full = self._build_geometry_from_parameters(dummy, keep_all_sources=True)
            ap = self.fm.build_axitra(
                geom_full,
                latlon=False,
                freesurface=True,
                aw=self.axitra_aw,
                ikmax=self.axitra_ikmax,
            )
            ap = self.fm.green(ap, quiet=True)
            self._green_cache_ap = ap
        return self._green_cache_ap

    # ------------------------------------------------------------------
    def clear_green_cache(self) -> None:
        """Clean up the cached Green's-function axitra files (call at end of run)."""
        if self._green_cache_ap is not None:
            try:
                self._green_cache_ap.clean()
            except Exception:
                pass
            self._green_cache_ap = None

    # ------------------------------------------------------------------
    def _evaluate_model(self, model: np.ndarray) -> Tuple[float, Optional[np.ndarray]]:
        """
        Internal method to evaluate one model vector.
        (Método interno para evaluar un solo vector de modelo.)
        
        Returns
        -------
        misfit : float
        synthetics : np.ndarray or None
        """
        ap = None
        try:
            if self.misfit_calc is None:
                return 1e10, None

            if self.use_green_cache:
                # Cached path: full mesh, Green's functions computed once, reused.
                # Only the cheap conv runs per model (zero-moment sources outside
                # the ellipse contribute nothing, so the result is identical).
                geom = self._build_geometry_from_parameters(model, keep_all_sources=True)
                ap = self._ensure_green_cache()
            else:
                # Standard path: rebuild geometry + recompute Green's functions per model.
                geom = self._build_geometry_from_parameters(model)
                ap = self.fm.build_axitra(
                    geom,
                    latlon=False,
                    freesurface=True,
                    aw=self.axitra_aw,
                    ikmax=self.axitra_ikmax,
                )
                ap = self.fm.green(ap, quiet=True)

            source_type = int(getattr(self.cfg.ellipse, "source_type", 4))
            _, sx, sy, sz = self.fm.conv(
                ap, geom, source_type=source_type, t0=float(self.cfg.ellipse.t0), quiet=True
            )

            synthetics = np.array([sx, sy, sz])
            synthetics = np.transpose(synthetics, (1, 2, 0))
            synthetics = np.transpose(synthetics, (0, 2, 1))

            # Filter synthetics to the same frequency band as observed data
            from kdellipspy.core.signal_utils import bandpass_filter_waveforms

            synthetics = bandpass_filter_waveforms(
                synthetics,
                self.time_array,
                freq1=float(self.cfg.ellipse.freq1),
                freq2=float(self.cfg.ellipse.freq2),
                zerophase=bool(getattr(self.cfg.ellipse, "zerophase", True)),
            )

            misfit = float(self.misfit_calc.l2_misfit(synthetics, use_full_signal=self.use_full_signal))
            return misfit, synthetics

        except Exception as exc:
            logger.error("Model evaluation failed: %s", exc)
            return 1e10, None

        finally:
            keep = False
            if self._na_cfg_runtime is not None:
                keep = bool(self._na_cfg_runtime.keep_axitra_files)
            elif self._mcmc_cfg_runtime is not None:
                keep = bool(self._mcmc_cfg_runtime.keep_axitra_files)
            # Never clean the cached Green's-function instance here — it is reused
            # across evaluations and cleaned once via clear_green_cache().
            if ap is not None and ap is not self._green_cache_ap and not keep:
                try:
                    ap.clean()
                except Exception:
                    pass

    # ------------------------------------------------------------------
    def objective_function(self, model: np.ndarray) -> float:
        """Evaluate the forward model and return L2 misfit for one parameter vector.
        (Evalúa el modelo forward y retorna el desajuste L2 para un vector de parámetros.)

        Side-effects
        ------------
        - Updates ``self._best_misfit_seen`` and ``self.best_synthetics`` when improved.
        - Cleans axitra temporary files unless the active config sets ``keep_axitra_files=True``.
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

        misfit, synthetics = self._evaluate_model(model)

        if misfit < self._best_misfit_seen:
            self._best_misfit_seen = misfit
            self._best_model_vec = np.asarray(model, dtype=float).copy()
            if synthetics is not None:
                self.best_synthetics = synthetics.copy()
            if self.checkpoint_path is not None:
                self._write_checkpoint(iter_est)

        print(
            f"[{log_tag}] iter={iter_est:05d} eval={self._eval_count:05d} "
            f"misfit={misfit:.6e} best={self._best_misfit_seen:.6e}",
            flush=True,
        )

        return misfit

    def _write_checkpoint(self, iter_est: int) -> None:
        """Vuelca el mejor modelo actual a ``self.checkpoint_path`` (escritura atómica)."""
        try:
            path = Path(self.checkpoint_path)
            lines = [
                "# Checkpoint del mejor modelo (se actualiza al mejorar el misfit)",
                f"# eval={self._eval_count}  iter={iter_est}  "
                f"misfit={self._best_misfit_seen:.6f}",
            ]
            for name, val in zip(self.param_names, self._best_model_vec):
                lines.append(f"{name:<14s} {val:12.4f}")
            tmp = path.with_suffix(path.suffix + ".tmp")
            tmp.write_text("\n".join(lines) + "\n")
            tmp.replace(path)   # rename atómico: nunca se lee a medias
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"No se pudo escribir el checkpoint: {exc}")

    def plot_fit(self, show: bool = True, save_path: Optional[str] = None) -> Tuple[Any, Any]:
        """
        Grafica el mejor ajuste de formas de onda encontrado hasta ahora.
        """
        if self.best_synthetics is None:
            print("No se han generado sintéticos aún. Corre la inversión primero.")
            return None, None
            
        from ..core.plotting import plot_waveform_fit
        
        station_names = [s.name for s in self.cfg.stations.stations]
        return plot_waveform_fit(
            observed=self.observed_waveforms,
            synthetic=self.best_synthetics,
            time=self.time_array,
            station_names=station_names,
            misfit=self._best_misfit_seen,
            show=show,
            save_path=save_path
        )


__all__ = [
    "NAModel",
    "MisfitCalculator",
    "NAResult",
    "BaseInversionModel",
]
