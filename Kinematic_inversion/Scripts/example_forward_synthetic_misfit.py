from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config_parser import ConfigParser
from src.forward_model import AxitraForwardModel
from src.inversion_na import MisfitCalculator
from src.signal_utils import build_azi_times_array


def _midpoint_model(cfg: ConfigParser) -> np.ndarray:
    vals = []
    for p in cfg.inversion_params.parameters:
        vals.append(0.5 * (float(p.min_val) + float(p.max_val)))
    if len(vals) != 7:
        raise ValueError(f"Expected 7 inversion parameters in input.ctl, got {len(vals)}")
    return np.asarray(vals, dtype=float)


def _stack_components(sx: np.ndarray, sy: np.ndarray, sz: np.ndarray) -> np.ndarray:
    arr = np.array([sx, sy, sz], dtype=float)
    arr = np.transpose(arr, (1, 2, 0))
    arr = np.transpose(arr, (0, 2, 1))
    return arr


def _add_noise(clean_waveforms: np.ndarray, noise_ratio: float, rng: np.random.Generator) -> tuple[np.ndarray, float]:
    if noise_ratio < 0.0:
        raise ValueError("noise_ratio must be >= 0")

    rms = float(np.sqrt(np.mean(clean_waveforms**2)))
    noise_std = max(rms * float(noise_ratio), 1e-14)
    noisy = clean_waveforms + noise_std * rng.standard_normal(clean_waveforms.shape)
    return noisy, noise_std


def _plot_comparison(
    output_dir: Path,
    time: np.ndarray,
    clean: np.ndarray,
    noisy: np.ndarray,
    station_names: list[str],
    max_stations: int,
    show: bool,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    nsta = clean.shape[0]
    nsel = min(int(max_stations), int(nsta)) if int(max_stations) > 0 else int(nsta)
    if nsel <= 0:
        return None

    comp_names = ["X", "Y", "Z"]
    fig, axes = plt.subplots(nsel, 3, figsize=(13, 2.6 * nsel), squeeze=False)

    for i in range(nsel):
        sta_name = station_names[i] if i < len(station_names) else f"STA_{i+1:02d}"
        for icomp, comp_name in enumerate(comp_names):
            ax = axes[i][icomp]
            y_clean = clean[i, icomp, :]
            y_noisy = noisy[i, icomp, :]
            peak = float(max(np.max(np.abs(y_clean)), np.max(np.abs(y_noisy)), 1e-12))

            ax.plot(time, y_clean, color="tab:blue", lw=1.2, label="Forward (clean)")
            ax.plot(time, y_noisy, color="tab:red", lw=1.0, alpha=0.85, label="Observed (clean+noise)")
            ax.set_xlim(float(time[0]), float(time[-1]))
            ax.set_ylim(-1.2 * peak, 1.2 * peak)
            ax.grid(True, alpha=0.3)

            if i == 0:
                ax.set_title(f"{comp_name}")
            if icomp == 0:
                ax.set_ylabel(sta_name)
            if i == 0 and icomp == 0:
                ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("Synthetic Forward Check: clean vs noisy")
    fig.supxlabel("Time (s)")
    fig.supylabel("Amplitude")
    fig.tight_layout()

    out_path = output_dir / "synthetic_forward_comparison.png"
    fig.savefig(out_path, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Example: run a synthetic forward with midpoint model from input.ctl, "
            "build noisy observed traces, and compute misfit."
        )
    )
    parser.add_argument("--input-ctl", default=str(ROOT / "input.ctl"))
    parser.add_argument("--axitra-dir", default=str(ROOT / "AXITRA2024"))
    parser.add_argument("--output-dir", default=str(ROOT / "output" / "synthetic_forward_example"))
    parser.add_argument("--noise-ratio", type=float, default=0.02, help="Noise std as fraction of global RMS")
    parser.add_argument("--seed", type=int, default=20260420)
    parser.add_argument("--max-stations-plot", type=int, default=6)
    parser.add_argument("--show", action="store_true", help="Show plots interactively")
    args = parser.parse_args()

    input_ctl = Path(args.input_ctl).resolve()
    axitra_dir = Path(args.axitra_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = ConfigParser(str(input_ctl))
    station_names = [s.name for s in cfg.stations.stations]
    midpoint_model = _midpoint_model(cfg)

    fm = AxitraForwardModel(str(input_ctl), axitra_dir=str(axitra_dir))
    geom = fm.build_geometry_with_ellipse_slip(midpoint_model)

    ap = None
    try:
        ap = fm.build_axitra(geom, latlon=False, freesurface=True)
        ap = fm.green(ap, quiet=True)
        time, sx, sy, sz = fm.conv(ap, geom, source_type=1, t0=float(cfg.ellipse.t0), quiet=True)
    finally:
        if ap is not None:
            try:
                ap.clean()
            except Exception:
                pass

    time = np.asarray(time, dtype=float)
    clean = _stack_components(sx, sy, sz)

    rng = np.random.default_rng(int(args.seed))
    observed_noisy, noise_std = _add_noise(clean, noise_ratio=float(args.noise_ratio), rng=rng)

    azi_times = build_azi_times_array(input_ctl_path=input_ctl)
    misfit_calc = MisfitCalculator(
        observed_waveforms=observed_noisy,
        time_array=time,
        azi_times_array=azi_times,
        time_window_s=20.0,
    )
    misfit = float(misfit_calc.l2_misfit(clean))
    diag_text = misfit_calc.diagnostics_summary(clean, max_stations=3)

    npz_path = output_dir / "synthetic_forward_example.npz"
    np.savez_compressed(
        npz_path,
        time=time,
        clean=clean,
        observed_noisy=observed_noisy,
        model=midpoint_model,
        azi_times=azi_times,
        misfit=np.array([misfit], dtype=float),
    )

    txt_path = output_dir / "synthetic_forward_summary.txt"
    with txt_path.open("w", encoding="utf-8") as f:
        f.write("Synthetic forward example summary\n")
        f.write("=" * 72 + "\n")
        f.write(f"Event: {cfg.source_position.event_name}\n")
        f.write(f"Stations: {len(station_names)}\n")
        f.write(f"Noise ratio: {float(args.noise_ratio):.4f}\n")
        f.write(f"Noise std (global): {noise_std:.6e}\n")
        f.write(f"Misfit (clean vs noisy-observed): {misfit:.6e}\n")
        f.write("Midpoint model (a1,a2,theta,np,tp,dmax,vr):\n")
        f.write("  " + ", ".join(f"{v:.6f}" for v in midpoint_model) + "\n")
        f.write("\nMisfit diagnostics:\n")
        f.write(diag_text + "\n")

    fig_path = _plot_comparison(
        output_dir=output_dir,
        time=time,
        clean=clean,
        noisy=observed_noisy,
        station_names=station_names,
        max_stations=int(args.max_stations_plot),
        show=bool(args.show),
    )

    print("=" * 72)
    print("SYNTHETIC FORWARD EXAMPLE")
    print("=" * 72)
    print(f"input.ctl: {input_ctl}")
    print(f"axitra_dir: {axitra_dir}")
    print(f"n_stations: {clean.shape[0]}")
    print(f"n_samples: {clean.shape[2]}")
    print(f"noise_ratio: {float(args.noise_ratio):.4f}")
    print(f"noise_std: {noise_std:.6e}")
    print(f"misfit: {misfit:.6e}")
    print(f"saved: {npz_path}")
    print(f"saved: {txt_path}")
    if fig_path is not None:
        print(f"saved: {fig_path}")
    else:
        print("matplotlib not available: comparison figure not generated")


if __name__ == "__main__":
    main()
