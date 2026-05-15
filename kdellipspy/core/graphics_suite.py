"""Unified graphics module for kinematic inversion outputs.

This module centralizes plotting tasks previously split across Graphics/\*.py.
It provides both a Python API and a CLI entrypoint.
"""

from __future__ import annotations

from dataclasses import dataclass
from csv import DictReader
import json
from pathlib import Path
from typing import Any, Iterable, Tuple
import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

@dataclass
class GraphicsConfig:
    base_dir: Path
    units_param: int = 1  # 1: displacement(cm), 2: velocity(cm/s)
    show: bool = True
    dpi: int = 300


class GraphicsSuite:
    """
    Visualization suite for the Kinematic Inversion results.
    (Suite de visualización para los resultados de la Inversión Cinemática.)

    Generates and saves convergence plots, waveform fits, and visualizations 
    of the spatial distribution of the rupture.
    (Genera y guarda gráficos de convergencia, ajuste de formas de onda, y 
    visualizaciones de la distribución espacial de la ruptura.)

    Attributes:
        base_dir (Path): Current run directory where plots will be saved.
                         (Directorio de ejecución actual donde se guardarán los plots.)
        show (bool): If True, displays interactive plots in addition to saving them.
                     (Si es True, muestra los gráficos interactivamente además de guardarlos.)
        output_dir (Path): Output folder within base_dir for exporting images.
                           (Carpeta 'output' dentro de base_dir para exportar imágenes.)
    """
    def __init__(self, base_dir: str | Path = "..", units_param: int = 1, show: bool = True, cfg: Any = None):
        self.base_dir = Path(base_dir).resolve()
        self.cfg = GraphicsConfig(base_dir=self.base_dir, units_param=units_param, show=show)
        self.inversion_cfg = cfg

        self.event_dir = self.base_dir / "Event" / "kine_files"
        self.fig_dir = self.base_dir / "Figures"
        self.fig_dir.mkdir(parents=True, exist_ok=True)

    def _read_dt_npts(self) -> Tuple[float, int]:
        kine_param_inc = self.base_dir / "src-covm_inkm" / "kine_param.inc"
        dt = None
        npts = None
        with kine_param_inc.open("r", encoding="utf-8") as f:
            for row in f:
                cols = row.split(",")
                if "parameter" in cols[0] and "dt=" in cols[0]:
                    dt = float(cols[0].split("=")[1])
                    npts = int(cols[4].strip().split("=")[1].replace(")", ""))
                    break
        if dt is None or npts is None:
            raise ValueError(f"Could not parse dt/npts from {kine_param_inc}")
        return dt, npts

    def _read_station_names(self) -> list[str]:
        stationn = self.base_dir / "Stations" / "stationn"
        with stationn.open("r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]

    def _read_seismograms(self):
        real_x = 100.0 * np.loadtxt(self.event_dir / "real_disp_x")
        best_x = 100.0 * np.loadtxt(self.event_dir / "best_disp_x")
        real_y = 100.0 * np.loadtxt(self.event_dir / "real_disp_y")
        best_y = 100.0 * np.loadtxt(self.event_dir / "best_disp_y")
        real_z = 100.0 * np.loadtxt(self.event_dir / "real_disp_z")
        best_z = 100.0 * np.loadtxt(self.event_dir / "best_disp_z")
        return real_x, best_x, real_y, best_y, real_z, best_z

    def _units_label(self) -> str:
        if self.cfg.units_param == 1:
            return "Displacement (cm)"
        if self.cfg.units_param == 2:
            return "Velocity (cm/s)"
        return ""

    @staticmethod
    def _round_to_1(x: float) -> float:
        if x == 0:
            return 0.0
        return round(x, -int(math.floor(math.log10(abs(x)))))

    def _normalize_na_rows(self, source: Any) -> tuple[list[dict[str, float]], list[str]]:
        if hasattr(source, "all_models"):
            rows = []
            param_names = list(getattr(source, "param_names", []))
            for model in source.all_models:
                row = {"iteration": float(model.iteration), "misfit": float(model.misfit)}
                for name, value in zip(param_names, model.model):
                    row[name] = float(value)
                rows.append(row)
            return rows, param_names

        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(path)

        if path.suffix.lower() == ".json":
            with path.open("r", encoding="utf-8") as f:
                payload = json.load(f)
            models = payload.get("models", [])
            param_names = list(payload.get("metadata", {}).get("param_names", []))
            rows = []
            for item in models:
                row = {"iteration": float(item.get("iteration", 0)), "misfit": float(item.get("misfit", 0.0))}
                for key, value in item.items():
                    if key not in {"model", "misfit", "iteration"}:
                        row[key] = float(value)
                if not param_names and "model" in item:
                    param_names = [f"p{i+1}" for i in range(len(item["model"]))]
                if "model" in item:
                    for name, value in zip(param_names, item["model"]):
                        row[name] = float(value)
                rows.append(row)
            return rows, param_names

        if path.suffix.lower() == ".csv":
            with path.open("r", encoding="utf-8") as f:
                reader = DictReader(f)
                rows = []
                fieldnames = reader.fieldnames or []
                param_names = [name for name in fieldnames if name not in {"iteration", "misfit"}]
                for item in reader:
                    row = {"iteration": float(item.get("iteration", 0)), "misfit": float(item.get("misfit", 0.0))}
                    for name in param_names:
                        if item.get(name, "") != "":
                            row[name] = float(item[name])
                    rows.append(row)
                return rows, param_names

        raise ValueError(f"Unsupported NA result source: {source}")

    def plot_na_results(self, source: Any) -> None:
        rows, param_names = self._normalize_na_rows(source)
        if not rows:
            raise ValueError("No NA results to plot")

        misfits = np.array([row["misfit"] for row in rows], dtype=float)
        iterations = np.array([row["iteration"] for row in rows], dtype=float)

        fig, axes = plt.subplots(2, 2, figsize=(13, 9))

        axes[0, 0].semilogy(np.arange(len(misfits)), misfits, color="tab:red", lw=1.2)
        axes[0, 0].set_title("Misfit evolution")
        axes[0, 0].set_xlabel("Model index")
        axes[0, 0].set_ylabel("Misfit")
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].scatter(iterations, misfits, s=14, alpha=0.7, color="tab:blue")
        axes[0, 1].set_title("Misfit by iteration")
        axes[0, 1].set_xlabel("Iteration")
        axes[0, 1].set_ylabel("Misfit")
        axes[0, 1].set_yscale("log")
        axes[0, 1].grid(True, alpha=0.3)

        best_idx = int(np.argmin(misfits))
        best_row = rows[best_idx]
        if param_names:
            axes[1, 0].bar(param_names, [best_row.get(name, np.nan) for name in param_names], color="tab:green")
            axes[1, 0].tick_params(axis="x", rotation=35)
            axes[1, 0].set_title("Best model parameters")
            axes[1, 0].set_ylabel("Value")
            axes[1, 0].grid(True, axis="y", alpha=0.3)
        else:
            axes[1, 0].axis("off")

        if param_names:
            p0 = param_names[0]
            p1 = param_names[1] if len(param_names) > 1 else param_names[0]
            axes[1, 1].scatter(
                [row.get(p0, np.nan) for row in rows],
                [row.get(p1, np.nan) for row in rows],
                c=misfits,
                cmap="viridis",
                s=20,
            )
            axes[1, 1].set_xlabel(p0)
            axes[1, 1].set_ylabel(p1)
            axes[1, 1].set_title("Parameter cloud colored by misfit")
            axes[1, 1].grid(True, alpha=0.3)
            cbar = fig.colorbar(axes[1, 1].collections[0], ax=axes[1, 1])
            cbar.set_label("Misfit")
        else:
            axes[1, 1].axis("off")

        fig.tight_layout()
        fig.savefig(self.fig_dir / "NA_results_summary.png", dpi=self.cfg.dpi)
        if self.cfg.show:
            plt.show()
        else:
            plt.close(fig)

        # Convergencia por parámetro
        self.plot_parameter_convergence(rows, param_names, misfits, iterations)

    def plot_parameter_convergence(
        self,
        rows: list[dict],
        param_names: list[str],
        misfits: np.ndarray,
        iterations: np.ndarray,
    ) -> None:
        """Una figura con 7 subplots (uno por parámetro): valor vs índice de modelo,
        coloreado por misfit, con la línea del mejor modelo superpuesta.
        """
        if not param_names:
            return

        n_params = len(param_names)
        # Distribuir en filas de 4 columnas máximo → para 7 params: 2 filas (4+3)
        ncols = 4
        nrows = math.ceil(n_params / ncols)

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * 3.5, nrows * 3.0),
            squeeze=False,
        )

        model_indices = np.arange(len(rows))
        best_idx = int(np.argmin(misfits))

        # Normalizar misfit para colormap
        vmin, vmax = float(np.nanmin(misfits)), float(np.nanmax(misfits))
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.get_cmap("plasma_r")

        for k, name in enumerate(param_names):
            row_i, col_i = divmod(k, ncols)
            ax = axes[row_i][col_i]

            values = np.array([row.get(name, np.nan) for row in rows], dtype=float)
            best_val = float(rows[best_idx].get(name, np.nan))

            sc = ax.scatter(
                model_indices,
                values,
                c=misfits,
                cmap=cmap,
                norm=norm,
                s=12,
                alpha=0.75,
                zorder=2,
            )
            # Línea horizontal del mejor modelo
            ax.axhline(best_val, color="crimson", lw=1.2, ls="--", zorder=3,
                       label=f"best = {best_val:.3g}")

            # Marcar separación de iteraciones con líneas verticales tenues
            iter_changes = np.where(np.diff(iterations.astype(int)) > 0)[0] + 1
            for xc in iter_changes:
                ax.axvline(xc, color="gray", lw=0.5, alpha=0.4, zorder=1)

            ax.set_title(name, fontsize=9, pad=3)
            ax.set_xlabel("Model index", fontsize=7)
            ax.set_ylabel("Value", fontsize=7)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=6, loc="upper right", framealpha=0.6)
            ax.grid(True, alpha=0.2)

        # Ocultar subplots sobrantes
        for k in range(n_params, nrows * ncols):
            row_i, col_i = divmod(k, ncols)
            axes[row_i][col_i].axis("off")

        # Colorbar global
        fig.subplots_adjust(right=0.88, hspace=0.55, wspace=0.4)
        cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Misfit", fontsize=9)

        fig.suptitle("Parameter convergence — NA search", fontsize=11, y=1.01)

        out_path = self.fig_dir / "NA_parameter_convergence.png"
        fig.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Gráfico guardado: {out_path}", flush=True)
        if self.cfg.show:
            plt.show()
        else:
            plt.close(fig)

    def load_na_results(self, source: Any) -> tuple[list[dict[str, float]], list[str]]:
        return self._normalize_na_rows(source)

    def plot_waveform_fit(
        self,
        observed: np.ndarray,
        synthetic: np.ndarray,
        time_array: np.ndarray,
        station_names: list[str],
        misfit: float | None = None,
    ) -> None:
        """
        Grafica las formas de onda observadas vs sintéticas.
        
        Parameters
        ----------
        observed : np.ndarray, shape (nsta, 3, npts)
        synthetic : np.ndarray, shape (nsta, 3, npts)
        time_array : np.ndarray, shape (npts,)
        station_names : list[str], length nsta
        misfit : float, opcional
        """
        nsta = observed.shape[0]
        fig, axs = plt.subplots(nsta, 3, figsize=(10, 1.5 * nsta), squeeze=False, sharex=True)
        
        comps = ['North', 'East', 'Z']
        colors_obs = ['black', 'black', 'black']
        colors_syn = ['red', 'red', 'red']

        title = "Mejor ajuste de formas de onda"
        if misfit is not None:
            title += f" (Misfit: {misfit:.4f})"
        fig.suptitle(title, fontsize=14, y=1.02)

        for i in range(nsta):
            for j in range(3):
                ax = axs[i, j]
                
                # Observado
                ax.plot(time_array, observed[i, j, :], color=colors_obs[j], label='Observed' if i==0 else None, linewidth=1.5)
                # Sintético
                ax.plot(time_array, synthetic[i, j, :], color=colors_syn[j], label='Synthetic' if i==0 else None, linewidth=1.2, linestyle='--')
                
                if i == 0:
                    ax.set_title(comps[j])
                
                if j == 0:
                    ax.set_ylabel(station_names[i], fontweight='bold')
                
                if i == nsta - 1:
                    ax.set_xlabel("Time (s)")
                
                ax.grid(True, alpha=0.3)
                ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelleft=False)

        fig.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
        plt.tight_layout()
        
        out_path = self.fig_dir / "Best_Waveform_Fit.png"
        fig.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Gráfico guardado: {out_path}", flush=True)
        if self.cfg.show:
            plt.show()
        else:
            plt.close(fig)

    def plot_record_section(
        self,
        waveforms: np.ndarray,
        time: np.ndarray,
        station_names: List[str],
        comp_name: str = "Z",
        scale: float = 1.0,
    ) -> None:
        """
        Plot waveforms as a record section (amplitude vs distance from hypocenter).
        (Grafica formas de onda como una sección de registro (amplitud vs distancia al hipocentro).)
        """
        # Calculate distances
        cfg = self.inversion_cfg
        if cfg is None:
            raise ValueError("plot_record_section requires an inversion config; pass cfg=... when creating GraphicsSuite")

        ev_lat = cfg.source_position.latitude
        ev_lon = cfg.source_position.longitude
        
        distances = []
        for s in cfg.stations.stations:
            dist = np.sqrt((s.latitude - ev_lat)**2 + (s.longitude - ev_lon)**2) * 111.0 # km approx
            distances.append(dist)
        
        distances = np.array(distances)
        nsta = waveforms.shape[0]

        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Normalize for plotting
        max_amp = np.max(np.abs(waveforms))
        norm_factor = (np.max(distances) - np.min(distances)) / (nsta * 2) if nsta > 1 else 1.0
        if max_amp == 0: max_amp = 1.0

        for i in range(nsta):
            # Scale trace for visibility in the section
            trace = waveforms[i, :] / max_amp * norm_factor * scale
            ax.plot(time, trace + distances[i], color='black', lw=0.7)
            ax.text(time[-1], distances[i], f" {station_names[i]} ({distances[i]:.1f} km)", 
                    va='center', fontsize=9)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Distance (km)")
        ax.set_title(f"Record Section - Component {comp_name}")
        ax.grid(True, alpha=0.3)
        
        out_path = self.fig_dir / f"Record_Section_{comp_name}.png"
        fig.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Record Section saved: {out_path}")
        plt.show()

    def plot_seismograms_all(self) -> None:
        dt, npts = self._read_dt_npts()
        stations = self._read_station_names()
        nsta = len(stations)
        real_x, best_x, real_y, best_y, real_z, best_z = self._read_seismograms()
        time = np.arange(npts) * dt

        fig, ax = plt.subplots(nsta, 3, figsize=(9, 10), squeeze=False)
        for j in range(nsta):
            i0 = npts * j
            i1 = npts * (j + 1)
            maxv = np.max(
                np.abs(
                    np.concatenate(
                        [
                            real_x[i0:i1],
                            best_x[i0:i1],
                            real_y[i0:i1],
                            best_y[i0:i1],
                            real_z[i0:i1],
                            best_z[i0:i1],
                        ]
                    )
                )
            )
            limit = self._round_to_1(maxv * 1.5)

            ax[j, 0].plot(time, real_x[i0:i1], "b")
            ax[j, 0].plot(time, best_x[i0:i1], "r")
            ax[j, 1].plot(time, real_y[i0:i1], "b")
            ax[j, 1].plot(time, best_y[i0:i1], "r")
            ax[j, 2].plot(time, real_z[i0:i1], "b")
            ax[j, 2].plot(time, best_z[i0:i1], "r")

            for k in range(3):
                ax[j, k].set_xlim(0, npts * dt)
                ax[j, k].set_ylim(-limit, limit)
                ax[j, k].set_yticks([-limit, -limit / 2, 0, limit / 2, limit])
                ax[j, k].text(0.03, 0.9, stations[j], transform=ax[j, k].transAxes)

        fig.suptitle("All Seismograms")
        fig.text(0.22, 0.95, "N-S", ha="center")
        fig.text(0.55, 0.95, "E-W", ha="center")
        fig.text(0.88, 0.95, "Z", ha="center")
        fig.text(0.5, 0.02, "Time [s]", ha="center")
        fig.text(0.02, 0.5, self._units_label(), va="center", rotation="vertical")
        fig.subplots_adjust(wspace=0.53, hspace=0.55, left=0.12, right=0.98, bottom=0.06, top=0.94)
        fig.savefig(self.fig_dir / "Seismog_all.png", dpi=self.cfg.dpi)

    def plot_seismograms_detailed(self) -> None:
        dt, npts = self._read_dt_npts()
        stations = self._read_station_names()
        nsta = len(stations)
        nrow_plot = nsta // 2 if nsta % 2 == 0 else (nsta // 2 + 1)
        real_x, best_x, real_y, best_y, real_z, best_z = self._read_seismograms()
        time = np.arange(npts) * dt

        components = [
            (real_x, best_x, "North component", "Seismog_north.png"),
            (real_y, best_y, "East component", "Seismog_east.png"),
            (real_z, best_z, "Vertical component", "Seismog_vertical.png"),
        ]

        for real_c, best_c, title, outname in components:
            fig = plt.figure(figsize=(7, 9))
            for j in range(nsta):
                i0 = npts * j
                i1 = npts * (j + 1)
                ax = plt.subplot(nrow_plot, 2, j + 1)
                ax.plot(time, real_c[i0:i1], "b")
                ax.plot(time, best_c[i0:i1], "r")
                ax.set_xlim(0, npts * dt)
                maxv = max(np.max(np.abs(real_c[i0:i1])), np.max(np.abs(best_c[i0:i1])))
                limit = self._round_to_1(maxv * 1.5)
                ax.set_ylim(-limit, limit)
                ax.set_yticks([-limit, -limit / 2, 0, limit / 2, limit])
                ax.text(0.03, 0.9, stations[j], transform=ax.transAxes)

            fig.text(0.5, 0.98, title, ha="center")
            fig.text(0.5, 0.01, "Time (s)", ha="center")
            fig.text(0.02, 0.5, self._units_label(), va="center", rotation="vertical")
            fig.subplots_adjust(wspace=0.35, hspace=0.50, bottom=0.06, top=0.95, left=0.15, right=0.96)
            fig.savefig(self.fig_dir / outname, dpi=self.cfg.dpi)

    def plot_parameters(self) -> None:
        param_file = self.event_dir / "kine_param"
        with param_file.open("r", encoding="utf-8") as f:
            first = f.readline().split()
        nparams = int(first[0])

        models = np.loadtxt(self.event_dir / "models.dat")
        nmodels = int(models.size / nparams)
        models = models.reshape((nmodels, nparams))
        misfits = np.loadtxt(self.event_dir / "misfits.dat")

        res_min = misfits.min()

        fig = plt.figure(figsize=(8, 7))
        ax1 = plt.subplot(2, 1, 1)
        ax1.plot(models[:, 0], "ro", label="Axis 1")
        ax1.plot(models[:, 1], "bo", label="Axis 2")
        ax1.set_title("Ellipse axes")
        ax1.set_xlabel("Models")
        ax1.set_ylabel("Length (km)")
        ax1.legend(prop={"size": 10})
        ax1.grid(True)

        ax2 = plt.subplot(2, 1, 2)
        x0 = models[:, 0] * models[:, 3] * np.cos(2 * math.pi * models[:, 4])
        y0 = models[:, 1] * models[:, 3] * np.sin(2 * math.pi * models[:, 4])
        ax2.plot(x0, "bo", label="x")
        ax2.plot(y0, "ro", label="y")
        ax2.legend(prop={"size": 10})
        ax2.set_title("Ellipse center with respect to hypocenter")
        ax2.set_xlabel("Models")
        ax2.set_ylabel("Location (km)")
        ax2.grid(True)
        plt.subplots_adjust(hspace=0.50)
        fig.savefig(self.fig_dir / "Convergence_axes_center.png", dpi=self.cfg.dpi)

        fig = plt.figure()
        plt.plot(models[:, 2] * 180.0, "ro")
        plt.title("Ellipse rotation angle")
        plt.xlabel("Models")
        plt.ylabel("Angle (degrees)")
        plt.grid(True)
        fig.savefig(self.fig_dir / "Convergence_angle.png", dpi=self.cfg.dpi)

        fig = plt.figure()
        plt.plot(models[:, 5], "ro")
        plt.title("Maximum slip")
        plt.xlabel("Models")
        plt.ylabel("Maximum slip (m)")
        plt.grid(True)
        fig.savefig(self.fig_dir / "Convergence_Dmax.png", dpi=self.cfg.dpi)

        fig = plt.figure()
        plt.plot(models[:, 6], "ro")
        plt.title("Rupture velocity")
        plt.xlabel("Models")
        plt.ylabel("Rupture velocity (km/s)")
        plt.grid(True)
        fig.savefig(self.fig_dir / "Convergence_Vr.png", dpi=self.cfg.dpi)

        fig = plt.figure()
        plt.plot(misfits, "ro")
        plt.ylim(0, 5)
        plt.title("Global convergence")
        plt.xlabel("Models")
        plt.ylabel("Misfit")
        plt.grid(True)
        text_str = "\n".join((f"Min. misfit = {res_min:0.3f}", f"Last misfit = {misfits[-1]:0.3f}"))
        plt.text(
            0.98,
            0.97,
            text_str,
            horizontalalignment="right",
            verticalalignment="top",
            transform=plt.gca().transAxes,
            bbox=dict(facecolor="w", alpha=0.8, edgecolor="k"),
        )
        fig.savefig(self.fig_dir / "Convergence_misfit.png", dpi=self.cfg.dpi)

    def plot_source_ellipse(self, geometry: 'FaultGeometry', title: str = "2D Slip Distribution") -> None:
        """
        Visualización 2D de la distribución de slip interpolada.
        """
        # 1. Obtener el mapa de slip máximo por subfalla
        slip_map = {}
        for sp in geometry.source_points:
            sf_idx = sp.subfault_index
            slip_map[sf_idx] = max(slip_map.get(sf_idx, 0.0), sp.displacement)
            
        # 2. Extraer los datos como puntos discretos (en km)
        x_pts = np.array([sf.x_m for sf in geometry.subfaults]) / 1000.0
        y_pts = np.array([sf.y_m for sf in geometry.subfaults]) / 1000.0
        z_pts = np.array([slip_map.get(sf.index, 0.0) for sf in geometry.subfaults])

        # 3. Malla regular para la interpolación (100x100 para un gradiente suave)
        resolucion = 100
        xi = np.linspace(x_pts.min(), x_pts.max(), resolucion)
        yi = np.linspace(y_pts.min(), y_pts.max(), resolucion)
        X, Y = np.meshgrid(xi, yi)

        # 4. Interpolar (Suavizado de la superficie)
        Z = griddata((x_pts, y_pts), z_pts, (X, Y), method='cubic', fill_value=0.0)
        Z = np.clip(Z, 0, None) # Evitar valores negativos por la interpolación
        max_slip = np.max(Z)

        # 5. Configurar la figura 2D
        fig, ax = plt.subplots(figsize=(10, 6))

        # 6. Dibujar el mapa de calor (niveles rellenos)
        # Puedes cambiar 'gnuplot' por 'magma' o 'inferno' si prefieres
        mapa_calor = ax.contourf(X, Y, Z, levels=50, cmap='gnuplot')

        # 7. Dibujar el contorno verde (Aspereza al 15% del máximo)
        umbral_verde = max_slip * 0.15
        if max_slip > 0:
            ax.contour(X, Y, Z, levels=[umbral_verde], colors='lime', 
                    linestyles='dashed', linewidths=1.5)

        # 8. Dibujar el Hipocentro
        ax.scatter(0, 0, marker='*', s=350, color='yellow', 
                edgecolors='black', linewidth=1.2, label='Hypocenter', zorder=10)

        # 9. Estilos y ejes
        cbar = fig.colorbar(mapa_calor, ax=ax)
        cbar.set_label("Slip (m)", fontsize=10)
        
        ax.set_xlabel("Along Strike (km)", fontsize=11)
        ax.set_ylabel("Along Dip (km)", fontsize=11)
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        
        # Invertir eje Y para la convención de profundidad
        ax.invert_yaxis()
        
        # Fundamental en 2D: Asegurar que la proporción X/Y sea real (1:1)
        # para que la falla no se vea estirada artificialmente
        ax.set_aspect('equal') 
        
        ax.legend(loc='upper right', frameon=True, shadow=True)

        # 10. Guardar la figura
        out_path = self.fig_dir / "Source_Slip_2D_Smooth.png"
        plt.tight_layout()
        plt.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Gráfico 2D de slip suavizado guardado: {out_path}", flush=True)
        
        if self.cfg.show:
            plt.show()
        else:
            plt.close(fig)

    def plot_stations_map(self, cfg: 'ConfigParser') -> None:
        """
        Plot station locations and their active components on a map.
        (Grafica la ubicación de las estaciones y sus componentes activas en un mapa.)
        """
        stations = cfg.stations.stations
        if not stations:
            print("No stations to plot.")
            return

        lats = [s.latitude for s in stations]
        lons = [s.longitude for s in stations]
        names = [s.name for s in stations]
        
        # Source position
        src_lat = cfg.source_position.latitude
        src_lon = cfg.source_position.longitude

        plt.figure(figsize=(10, 8))
        
        # Plot stations
        for i, s in enumerate(stations):
            # Base marker for station
            plt.scatter(s.longitude, s.latitude, marker='^', s=100, color='blue', label='Station' if i==0 else "")
            plt.text(s.longitude, s.latitude, f"  {s.name}", verticalalignment='bottom', fontsize=9)
            
            # Indicators for components
            comps = []
            if s.use_n: comps.append('N')
            if s.use_e: comps.append('E')
            if s.use_z: comps.append('Z')
            
            comp_str = f"({','.join(comps)})"
            plt.text(s.longitude, s.latitude, f"  {comp_str}", verticalalignment='top', fontsize=7, color='gray')

        # Plot source (Hypocenter)
        plt.scatter(src_lon, src_lat, marker='*', s=200, color='red', label='Hypocenter')
        
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"Station Distribution - {cfg.source_position.event_name}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        out_path = self.fig_dir / "Station_Distribution.png"
        plt.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Gráfico de estaciones guardado: {out_path}", flush=True)
        if self.cfg.show:
            plt.show()
        else:
            plt.close()

    def plot_synthetic_components(self, time: np.ndarray, synthetic: np.ndarray, station_names: list[str]) -> None:
        """
        Plot synthetic seismograms for all stations and components.
        (Grafica sismogramas sintéticos para todas las estaciones y componentes.)
        
        Parameters
        ----------
        time : np.ndarray, shape (npts,)
        synthetic : np.ndarray, shape (nsta, 3, npts)
        station_names : list[str]
        """
        nsta = synthetic.shape[0]
        npts = synthetic.shape[2]
        
        fig, axs = plt.subplots(nsta, 3, figsize=(12, 2 * nsta), squeeze=False, sharex=True)
        comps = ['North (x)', 'East (y)', 'Vertical (z)']
        colors = ['red', 'blue', 'green']
        
        for i in range(nsta):
            for j in range(3):
                ax = axs[i, j]
                ax.plot(time, synthetic[i, j, :], color=colors[j], linewidth=1.0)
                
                if i == 0:
                    ax.set_title(comps[j])
                if j == 0:
                    ax.set_ylabel(f"{station_names[i]}\nAmp", fontweight='bold')
                if i == nsta - 1:
                    ax.set_xlabel("Time (s)")
                
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        out_path = self.fig_dir / "Synthetic_Seismograms.png"
        fig.savefig(out_path, dpi=self.cfg.dpi, bbox_inches="tight")
        print(f"✓ Gráfico de componentes guardado: {out_path}", flush=True)
        if self.cfg.show:
            plt.show()
        else:
            plt.close(fig)

    def plot_all(self) -> None:
        self.plot_seismograms_all()
        self.plot_seismograms_detailed()
        self.plot_parameters()
        self.plot_source()
        if self.cfg.show:
            plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified plotting suite for kinematic inversion")
    parser.add_argument("task", choices=["all", "seismog-all", "seismog-detailed", "parameters", "source"])
    parser.add_argument("--base-dir", default="..", help="Kinematic_inversion base directory")
    parser.add_argument("--units", type=int, default=1, help="1 displacement (cm), 2 velocity (cm/s)")
    parser.add_argument("--no-show", action="store_true", help="Save figures without opening windows")
    args = parser.parse_args()

    suite = GraphicsSuite(base_dir=args.base_dir, units_param=args.units, show=not args.no_show)

    if args.task == "all":
        suite.plot_all()
    elif args.task == "seismog-all":
        suite.plot_seismograms_all()
    elif args.task == "seismog-detailed":
        suite.plot_seismograms_detailed()
    elif args.task == "parameters":
        suite.plot_parameters()
    elif args.task == "source":
        suite.plot_source()

    if suite.cfg.show:
        plt.show()


if __name__ == "__main__":
    main()
