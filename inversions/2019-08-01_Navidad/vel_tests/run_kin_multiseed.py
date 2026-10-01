"""Cinemático + dt0 con varias semillas NA sobre UNA base de axitra (fast path).

Uso:  KIN_CASE=loc_mid python run_kin_multiseed.py      (caso relativo a ../)
Escribe ../<caso>/output_multiseed/: seed<k>.joblib, summary.txt
"""

import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "loc_mid_dt0"))
import kdellipspy as kde  # noqa: E402
from kdellipspy.inversion.kinematic.model_na import NAConfig  # noqa: E402
from run_kin_dt0 import DT0_RANGE, KinematicDt0NA  # noqa: E402

KIN_CASE = os.environ.get("KIN_CASE", "loc_mid")
CASE = HERE.parent / KIN_CASE
SEEDS = [int(s) for s in os.environ.get("NA_SEEDS", "0,1,2,3,4").split(",")]
NA = dict(n_samples_initial=400, n_samples_iteration=100, n_iterations=30, n_cells_resample=20, n_jobs=1)


class MultiSeedNA(KinematicDt0NA):
    """run_na_search limpia el caché de Green al empezar (y con él la base del
    fast path): aquí se conserva entre semillas y se limpia una vez al final."""

    def clear_green_cache(self):
        if getattr(self, "_keep_cache", False):
            return
        super().clear_green_cache()


def main():
    cfg = kde.ConfigParser(filepath=str(CASE / "input.ctl"))
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE / "input.ctl"), data_dir=str(CASE / "DATA"))
    inv = MultiSeedNA(config=cfg, observed_waveforms=observed, time_array=time_array,
                      azi_times_array=kde.build_azi_times_array(config=cfg, model_name="iasp91"))
    inv.use_green_cache = True
    inv.param_ranges = np.vstack([inv.param_ranges, DT0_RANGE])
    inv.param_names = list(inv.param_names) + ["dt0 (s)"]
    out = CASE / "output_multiseed"
    out.mkdir(exist_ok=True)

    inv._evaluate_model(np.mean(inv.param_ranges, axis=1))  # construye caché de Green + base
    inv._keep_cache = True
    rows = []
    for seed in SEEDS:
        inv.checkpoint_path = out / f"seed{seed}_best_live.txt"
        result = inv.run_na_search(NAConfig(random_seed=seed, **NA))
        result.save(str(out / f"seed{seed}.joblib"))
        rows.append((result.best_model.misfit, seed, np.asarray(result.best_model.model)))
        print(f"[seed {seed}] best misfit {result.best_model.misfit:.4f}", flush=True)
    inv._keep_cache = False
    inv.clear_green_cache()

    ms = np.array([r[0] for r in rows])
    models = np.array([r[2] for r in rows])
    lines = [f"{KIN_CASE}: {len(SEEDS)} semillas x {NA['n_samples_initial'] + NA['n_samples_iteration'] * NA['n_iterations']} evals",
             f"misfit  best {ms.min():.4f}  mean {ms.mean():.4f}  std {ms.std():.4f}  (por semilla: {', '.join(f'{m:.4f}' for m in ms)})",
             f"{'param':<22s} {'mejor':>9s} {'media':>9s} {'std':>8s}"]
    best = models[np.argmin(ms)]
    for i, name in enumerate(inv.param_names):
        lines.append(f"{name:<22s} {best[i]:9.3f} {models[:, i].mean():9.3f} {models[:, i].std():8.3f}")
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
