from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config_parser import ConfigParser
from src.geometry import EllipseDiagnostics


def midpoint_model(cfg: ConfigParser) -> np.ndarray:
    vals = [0.5 * (float(p.min_val) + float(p.max_val)) for p in cfg.inversion_params.parameters]
    return np.asarray(vals, dtype=float)


def best_model_from_csv(path: Path) -> np.ndarray:
    keys = [
        "a1 (km)",
        "a2 (km)",
        "theta (x pi)",
        "np (frac)",
        "tp (x 2pi)",
        "dmax (m)",
        "vr (km/s)",
    ]
    best_row = None
    best_misfit = np.inf
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            m = float(row["misfit"])
            if m < best_misfit:
                best_misfit = m
                best_row = row

    if best_row is None:
        raise ValueError(f"No rows found in {path}")

    return np.asarray([float(best_row[k]) for k in keys], dtype=float)


def parse_model_values(values: str) -> np.ndarray:
    parts = [p.strip() for p in values.split(",") if p.strip()]
    if len(parts) != 7:
        raise ValueError("--model-values requires 7 comma-separated values")
    return np.asarray([float(p) for p in parts], dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot ellipse/slip factors used by forward model")
    parser.add_argument("--input-ctl", default=str(ROOT / "input.ctl"))
    parser.add_argument("--model-source", choices=["midpoint", "best_csv", "manual"], default="best_csv")
    parser.add_argument("--best-csv", default=str(ROOT / "output" / "na_results.csv"))
    parser.add_argument("--model-values", default=None)
    parser.add_argument("--out", default=str(ROOT / "Figures" / "forward_ellipse_diagnostics.png"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    input_ctl = Path(args.input_ctl).resolve()
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    cfg = ConfigParser(str(input_ctl))
    if args.model_values:
        model = parse_model_values(args.model_values)
        source = "manual"
    elif args.model_source == "midpoint":
        model = midpoint_model(cfg)
        source = "midpoint"
    else:
        model = best_model_from_csv(Path(args.best_csv).resolve())
        source = "best_csv"

    diag = EllipseDiagnostics(cfg)
    result = diag.evaluate(model)

    print("=" * 72)
    print("ELLIPSE FORWARD DIAGNOSTICS")
    print("=" * 72)
    print(f"model_source: {source}")
    print(f"model: {model.tolist()}")
    print(f"active_subfault_ratio: {result.active_subfault_ratio:.3f}")
    print(f"ellipse center (km): ({result.xe_m/1000.0:.3f}, {result.ye_m/1000.0:.3f})")
    print(f"a1, a2 (km): ({result.a1_m/1000.0:.3f}, {result.a2_m/1000.0:.3f})")

    diag.plot(result, save_path=str(out), show=bool(args.show))
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
