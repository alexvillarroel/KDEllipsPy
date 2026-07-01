"""kde-dashboard · genera el dashboard (plot_yolo) desde un inversion_result.joblib.

    kde-dashboard                       # ./output/inversion_result.joblib -> ./output/figures/dashboard.pdf
    kde-dashboard ruta/al/result.joblib
    kde-dashboard result.joblib -o salida.pdf
"""
import argparse
from pathlib import Path

from kdellipspy.inversion.base import NAResult


def main(argv=None):
    p = argparse.ArgumentParser(prog="kde-dashboard", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("joblib", nargs="?", type=Path,
                   default=Path("output/inversion_result.joblib"))
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="PDF de salida (default <joblib>/../figures/dashboard.pdf)")
    args = p.parse_args(argv)

    if not args.joblib.exists():
        p.error(f"no existe el joblib: {args.joblib}")
    out = args.output or args.joblib.parent / "figures" / "dashboard.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    NAResult.load(str(args.joblib)).plot_yolo(save_path=str(out), show=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
