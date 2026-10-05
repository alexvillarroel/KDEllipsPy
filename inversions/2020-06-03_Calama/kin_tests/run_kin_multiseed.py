"""Cinematico + dt0 con 5 semillas NA sobre UNA base de axitra (camino rapido).

Un ranking con una sola semilla engana: el NA es estocastico y la dispersion
entre semillas suele ser del orden de las diferencias entre configuraciones.
Por eso cada caso se corre con NA_SEEDS semillas y se reporta best/media/std.

Uso:  KIN_CASE=np1_isola_potin_b0.04-0.15 python run_kin_multiseed.py
Escribe <caso>/multiseed/: seed<k>.joblib, summary.txt, station_misfit.txt
"""

import os
import sys
from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde
from kdellipspy.inversion.kinematic.model_na import NAConfig

HERE = Path(__file__).resolve().parent
CASE = HERE / os.environ["KIN_CASE"]
SEEDS = [int(s) for s in os.environ.get("NA_SEEDS", "0,1,2,3,4").split(",")]
NA = dict(n_samples_initial=400, n_samples_iteration=100, n_iterations=30,
          n_cells_resample=20, n_jobs=1)
DT0_RANGE = (-5.0, 5.0)  # s, >0 = mas tarde


class KinematicDt0NA(kde.NAInversionModel):
    """8o parametro dt0: corrige la hora de origen sumandose al tiempo de
    ruptura de todas las subfallas (axitra lo aplica como exp(-i*w*delay),
    desplazamiento exacto). run_na_search limpia el cache de Green al empezar,
    y con el la base del camino rapido: aqui se conserva entre semillas."""

    def _build_geometry_from_parameters(self, model, keep_all_sources=False):
        geom = super()._build_geometry_from_parameters(model[:7], keep_all_sources=keep_all_sources)
        for sp in geom.source_points:
            sp.rupture_time_s += float(model[7])
        return geom

    def clear_green_cache(self):
        if getattr(self, "_keep_cache", False):
            return
        super().clear_green_cache()


def station_misfit(result, cfg):
    """Misfit por estacion x componente (R, T, Z) del mejor modelo."""
    num, den, win = result._misfit_breakdown_arrays()
    names = [s.name for s in cfg.stations.stations]
    lines = [f"misfit por estacion (ventana {win} s; R/Z desde tP, T desde tS)",
             f"{'sta':<7s}{'R':>9s}{'T':>9s}{'Z':>9s}{'total':>9s}{'% del num':>11s}"]
    tot_num = num.sum()
    for j, n in enumerate(names):
        per = [num[j, c] / den[j, c] if den[j, c] > 0 else np.nan for c in range(3)]
        tot = num[j].sum() / den[j].sum() if den[j].sum() > 0 else np.nan
        lines.append(f"{n:<7s}{per[0]:9.3f}{per[1]:9.3f}{per[2]:9.3f}{tot:9.3f}"
                     f"{100 * num[j].sum() / tot_num:11.1f}")
    lines.append(f"{'TOTAL':<7s}{'':>27s}{num.sum() / den.sum():9.3f}")
    return "\n".join(lines) + "\n"


def main():
    cfg = kde.ConfigParser(filepath=str(CASE / "input.ctl"))
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE / "input.ctl"), data_dir=str(CASE / "DATA"))
    inv = KinematicDt0NA(config=cfg, observed_waveforms=observed, time_array=time_array,
                         azi_times_array=kde.build_azi_times_array(config=cfg, model_name="iasp91"))
    inv.use_green_cache = True
    inv.param_ranges = np.vstack([inv.param_ranges, DT0_RANGE])
    inv.param_names = list(inv.param_names) + ["dt0 (s)"]
    out = CASE / "multiseed"
    out.mkdir(exist_ok=True)

    inv._evaluate_model(np.mean(inv.param_ranges, axis=1))  # construye cache de Green + base
    inv._keep_cache = True
    rows, best_result = [], None
    for seed in SEEDS:
        inv.checkpoint_path = out / f"seed{seed}_best_live.txt"
        result = inv.run_na_search(NAConfig(random_seed=seed, **NA))
        result.save(str(out / f"seed{seed}.joblib"))
        rows.append((result.best_model.misfit, seed, np.asarray(result.best_model.model)))
        if best_result is None or result.best_model.misfit < best_result.best_model.misfit:
            best_result = result
        print(f"[seed {seed}] best misfit {result.best_model.misfit:.4f}", flush=True)
    inv._keep_cache = False
    inv.clear_green_cache()

    ms = np.array([r[0] for r in rows])
    models = np.array([r[2] for r in rows])
    evals = NA["n_samples_initial"] + NA["n_samples_iteration"] * NA["n_iterations"]
    lines = [f"{os.environ['KIN_CASE']}: {len(SEEDS)} semillas x {evals} evals",
             f"misfit  best {ms.min():.4f}  mean {ms.mean():.4f}  std {ms.std():.4f}"
             f"  (por semilla: {', '.join(f'{m:.4f}' for m in ms)})",
             f"{'param':<22s} {'mejor':>9s} {'media':>9s} {'std':>8s}"]
    best = models[np.argmin(ms)]
    for i, name in enumerate(inv.param_names):
        lines.append(f"{name:<22s} {best[i]:9.3f} {models[:, i].mean():9.3f} {models[:, i].std():8.3f}")
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)

    try:
        (out / "station_misfit.txt").write_text(station_misfit(best_result, cfg))
        print("\n" + (out / "station_misfit.txt").read_text(), flush=True)
    except Exception as exc:  # noqa: BLE001
        print(f"[aviso] no se pudo escribir station_misfit.txt: {exc}", file=sys.stderr, flush=True)
    joblib.dump({"case": os.environ["KIN_CASE"], "misfits": ms, "models": models,
                 "param_names": list(inv.param_names)}, out / "multiseed.joblib")


if __name__ == "__main__":
    main()
