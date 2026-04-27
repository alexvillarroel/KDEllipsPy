"""Unified graphics module for kinematic inversion outputs.

This module centralizes plotting tasks previously split across Graphics/*.py.
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
    def __init__(self, base_dir: str | Path = "..", units_param: int = 1, show: bool = True):
        self.base_dir = Path(base_dir).resolve()
        self.cfg = GraphicsConfig(base_dir=self.base_dir, units_param=units_param, show=show)

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

    def load_na_results(self, source: Any) -> tuple[list[dict[str, float]], list[str]]:
        return self._normalize_na_rows(source)

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

    def plot_source(self) -> None:
        kine_in = self.event_dir / "kine.in"
        with kine_in.open("r", encoding="utf-8") as f:
            for _ in range(5):
                f.readline()
            nsource = int(f.readline().split()[0])
            xhypo, yhypo, nstr, ndip = f.readline().split()[0:4]
            nstr = int(nstr)
            ndip = int(ndip)

        faille_in = self.base_dir / "Faille" / "faille.in"
        with faille_in.open("r", encoding="utf-8") as f:
            xsize, ysize = f.readline().split()[0:2]
        nx1 = int(xsize) / 1000
        ny1 = int(ysize) / 1000

        hist = np.loadtxt(self.base_dir / "Event" / "axi.hist").reshape((nsource, 8))
        dx = float(hist[1, 5]) / 4.0
        dy = float(hist[1, 6]) / 4.0
        xhypo = float(xhypo) / dx
        yhypo = float(yhypo) / dy

        slip = np.loadtxt(self.event_dir / "slip_fin.dat")
        tr = np.loadtxt(self.event_dir / "tr_fin.dat")
        nx = nstr * 4
        ny = ndip * 4
        slip = slip.reshape((ny, nx))
        tr = tr.reshape((ny, nx))

        fig = plt.figure(figsize=(8, 7))
        plt.imshow(slip, vmin=0, vmax=float(np.max(slip)), interpolation="bilinear", cmap="binary")
        plt.plot([xhypo], [yhypo], "wo", ms=7)
        plt.xticks(np.arange(2) * nx, ("0", nx1), fontsize=12)
        plt.yticks(np.arange(2) * ny, (ny1, ""), fontsize=12)
        plt.xlim(0, nx)
        cbar = plt.colorbar()
        cbar.ax.set_title("(m)")
        plt.title("Slip distribution")
        plt.xlabel("Along strike (km)")
        plt.ylabel("Along dip (km)")
        plt.subplots_adjust(left=0.08, right=1.00, bottom=0.08, top=0.92)
        fig.savefig(self.fig_dir / "Source_slip_distribution.png", dpi=self.cfg.dpi)

        fig = plt.figure(figsize=(8, 7))
        plt.imshow(tr, vmin=0, vmax=float(np.max(tr)), interpolation="bilinear", cmap="binary")
        plt.plot([xhypo], [yhypo], "wo", ms=7)
        plt.xticks(np.arange(2) * nx, ("0", nx1), fontsize=12)
        plt.yticks(np.arange(2) * ny, (ny1, ""), fontsize=12)
        plt.xlim(0, nx)
        cbar = plt.colorbar()
        cbar.ax.set_title("(s)")
        plt.title("Rupture time")
        plt.xlabel("Along strike (km)")
        plt.ylabel("Along dip (km)")
        plt.subplots_adjust(left=0.08, right=1.00, bottom=0.08, top=0.92)
        fig.savefig(self.fig_dir / "Source_rupture_time.png", dpi=self.cfg.dpi)

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
