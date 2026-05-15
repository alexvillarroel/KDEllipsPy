from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config_parser import ConfigParser
from src.forward_model import AxitraForwardModel
from src.inversion_na import MisfitCalculator
from src.signal_utils import build_azi_times_array, load_and_filter_observed_data


def _parse_model_values(values: str) -> np.ndarray:
    parts = [p.strip() for p in values.split(",") if p.strip()]
    if len(parts) != 7:
        raise ValueError("--model-values must contain 7 comma-separated numbers")
    return np.asarray([float(p) for p in parts], dtype=float)


def _midpoint_model(cfg: ConfigParser) -> np.ndarray:
    vals = []
    for p in cfg.inversion_params.parameters:
        vals.append(0.5 * (float(p.min_val) + float(p.max_val)))
    if len(vals) != 7:
        raise ValueError(f"Expected 7 inversion parameters, got {len(vals)}")
    return np.asarray(vals, dtype=float)


def _best_model_from_csv(csv_path: Path) -> np.ndarray:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing NA CSV file: {csv_path}")

    keys = [
        "a1 (km)",
        "a2 (km)",
        "theta (x pi)",
        "np (frac)",
        "tp (x 2pi)",
        "dmax (m)",
        "vr (km/s)",
    ]

    best_misfit = np.inf
    best_row = None
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            misfit = float(row["misfit"])
            if misfit < best_misfit:
                best_misfit = misfit
                best_row = row

    if best_row is None:
        raise ValueError(f"No rows found in {csv_path}")
    return np.asarray([float(best_row[k]) for k in keys], dtype=float)


def _active_subfault_ratio(geom, eps: float = 1e-14) -> float:
    by_sf = {}
    for sp in geom.source_points:
        idx = int(sp.subfault_index)
        by_sf[idx] = max(by_sf.get(idx, 0.0), abs(float(sp.amplitude)))
    if not by_sf:
        return 0.0
    active = sum(1 for amp in by_sf.values() if amp > eps)
    return float(active) / float(len(by_sf))


def _select_station_indices(station_names: List[str], selected: Optional[str], max_stations: int) -> List[int]:
    if selected:
        wanted = {s.strip().upper() for s in selected.split(",") if s.strip()}
        idx = [i for i, name in enumerate(station_names) if name.upper() in wanted]
    else:
        idx = list(range(len(station_names)))

    if max_stations > 0:
        idx = idx[:max_stations]
    return idx


def _plot_obs_vs_syn(
    out_dir: Path,
    time: np.ndarray,
    observed: np.ndarray,
    synthetic: np.ndarray,
    station_names: List[str],
    indices: List[int],
    show: bool,
) -> None:
    comp_names = ["X", "Y", "Z"]

    nsel = len(indices)
    if nsel == 0:
        raise ValueError("No stations selected for plotting")

    ncols = 2
    nrows = (nsel + ncols - 1) // ncols

    for comp_idx, comp in enumerate(comp_names):
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 2.6 * nrows), squeeze=False)
        for k, ista in enumerate(indices):
            r = k // ncols
            c = k % ncols
            ax = axes[r][c]

            obs = observed[ista, comp_idx, :]
            syn = synthetic[ista, comp_idx, :]
            peak = float(max(np.max(np.abs(obs)), np.max(np.abs(syn)), 1e-12))

            ax.plot(time, obs, color="tab:blue", lw=1.25, label="Observed")
            ax.plot(time, syn, color="tab:red", lw=1.25, label="Synthetic")
            ax.set_title(station_names[ista])
            ax.set_xlim(float(time[0]), float(time[-1]))
            ax.set_ylim(-1.2 * peak, 1.2 * peak)
            ax.grid(True, alpha=0.3)

            if k == 0:
                ax.legend(loc="upper right", fontsize=8)

        for k in range(nsel, nrows * ncols):
            r = k // ncols
            c = k % ncols
            axes[r][c].axis("off")

        fig.suptitle(f"Forward check - component {comp}")
        fig.supxlabel("Time (s)")
        fig.supylabel("Amplitude")
        fig.tight_layout()
        out_file = out_dir / f"forward_check_{comp}.png"
        fig.savefig(out_file, dpi=180)
        print(f"Saved: {out_file}")

    if show:
        plt.show()
    else:
        plt.close("all")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a single forward model and plot Observed vs Synthetic waveforms"
    )
    parser.add_argument("--input-ctl", default=str(ROOT / "input.ctl"))
    parser.add_argument("--axitra-dir", default=str(ROOT / "AXITRA2024"))
    parser.add_argument("--data-dir", default=str(ROOT / "DATA"))
    parser.add_argument("--output-dir", default=str(ROOT / "Figures"))
    parser.add_argument("--model-source", choices=["midpoint", "best_csv"], default="midpoint")
    parser.add_argument("--best-csv", default=str(ROOT / "output" / "na_results.csv"))
    parser.add_argument(
        "--model-values",
        default=None,
        help="7 comma-separated values: a1,a2,theta,np,tp,dmax,vr (overrides model-source)",
    )
    parser.add_argument("--stations", default=None, help="Comma-separated station names")
    parser.add_argument("--max-stations", type=int, default=10)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    input_ctl = Path(args.input_ctl).resolve()
    axitra_dir = Path(args.axitra_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = ConfigParser(str(input_ctl))
    station_names = [s.name for s in cfg.stations.stations]

    observed, time = load_and_filter_observed_data(
        input_ctl_path=input_ctl,
        data_dir=data_dir,
        prefer_raw=False,
    )

    if args.model_values is not None:
        model = _parse_model_values(args.model_values)
        model_origin = "model_values"
    elif args.model_source == "best_csv":
        model = _best_model_from_csv(Path(args.best_csv).resolve())
        model_origin = "best_csv"
    else:
        model = _midpoint_model(cfg)
        model_origin = "midpoint"

    fm = AxitraForwardModel(str(input_ctl), axitra_dir=str(axitra_dir))
    geom = fm.build_geometry_with_ellipse_slip(model)
    active_ratio = _active_subfault_ratio(geom)

    ap = None
    try:
        ap = fm.build_axitra(geom, latlon=False, freesurface=True)
        ap = fm.green(ap, quiet=True)
        _, sx, sy, sz = fm.conv(ap, geom, source_type=1, t0=float(cfg.ellipse.t0), quiet=True)
    finally:
        if ap is not None:
            try:
                ap.clean()
            except Exception:
                pass

    synthetic = np.array([sx, sy, sz])
    synthetic = np.transpose(synthetic, (1, 2, 0))
    synthetic = np.transpose(synthetic, (0, 2, 1))

    if synthetic.shape != observed.shape:
        raise ValueError(
            f"Observed shape {observed.shape} != synthetic shape {synthetic.shape}. "
            "Revisa Npts/Stations/units en input.ctl y archivos DATA."
        )

    azi_times = build_azi_times_array(input_ctl_path=input_ctl)
    misfit_calc = MisfitCalculator(observed, time, azi_times_array=azi_times, time_window_s=20.0)
    misfit = misfit_calc.l2_misfit(synthetic)

    print("=" * 72)
    print("FORWARD CHECK")
    print("=" * 72)
    print(f"model_origin: {model_origin}")
    print(f"model       : {model.tolist()}")
    print(f"active_subfault_ratio: {active_ratio:.3f}")
    print(f"misfit (num/den): {misfit:.6e}")

    indices = _select_station_indices(station_names, args.stations, args.max_stations)
    print(f"stations plotted: {len(indices)} / {len(station_names)}")

    _plot_obs_vs_syn(
        out_dir=out_dir,
        time=time,
        observed=observed,
        synthetic=synthetic,
        station_names=station_names,
        indices=indices,
        show=bool(args.show),
    )


if __name__ == "__main__":
    main()
