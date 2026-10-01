"""Inversión cinemática + corrección de hora de origen dt0 (caso: KIN_CASE, loc_mid por defecto).

Igual que `kdellipspy ../loc_mid` (mismo input.ctl, DATA, NA, caché de Green),
con un 8º parámetro dt0 (s, >0 = más tarde) sumado al tiempo de ruptura de
todas las subfallas (axitra lo aplica como exp(-i*w*delay): desplazamiento
exacto). loc_mid usa el hipocentro medio CSN-ISC con la hora de origen CSN;
el dinámico midió ~+3 s de desfase.
"""

import os
from pathlib import Path

import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
KIN_CASE = os.environ.get("KIN_CASE", "loc_mid")
CASE = HERE.parent / KIN_CASE
DT0_RANGE = (-5.0, 5.0)


class KinematicDt0NA(kde.NAInversionModel):
    def _build_geometry_from_parameters(self, model, keep_all_sources=False):
        geom = super()._build_geometry_from_parameters(model[:7], keep_all_sources=keep_all_sources)
        for sp in geom.source_points:
            sp.rupture_time_s += float(model[7])
        return geom


def main():
    cfg = kde.ConfigParser(filepath=str(CASE / "input.ctl"))
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE / "input.ctl"), data_dir=str(CASE / "DATA"))
    inv = KinematicDt0NA(
        config=cfg, observed_waveforms=observed, time_array=time_array,
        azi_times_array=kde.build_azi_times_array(config=cfg, model_name="iasp91"),
    )
    inv.use_green_cache = True
    inv.param_ranges = np.vstack([inv.param_ranges, DT0_RANGE])
    inv.param_names = list(inv.param_names) + ["dt0 (s)"]
    out = HERE / "output" if KIN_CASE == "loc_mid" else CASE / "output_dt0"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    inv.checkpoint_path = out / "best_model_live.txt"

    # Referencia: mejor modelo original de loc_mid (dt0=0 -> 0.2538 con su propio
    # modelo de velocidades; con otro medio solo informa).
    ref = np.array([6.5028, 9.2309, 0.9997, 0.9093, 0.7115, 1.5450, 1.3605, 0.0])
    print(f"[preflight] {KIN_CASE}: mejor original de loc_mid, dt0=0: misfit {inv._evaluate_model(ref)[0]:.4f}", flush=True)

    result = inv.run_na_search()
    result.save(str(out / "inversion_result.joblib"))
    print("best misfit", result.best_model.misfit, flush=True)
    for name, val in zip(result.param_names, result.best_model.model):
        print(f"  {name:<28s} {val:10.4f}", flush=True)
    for method, fname in [("plot_fit", "waveform_fit.png"), ("plot_convergence", "parameter_convergence.png")]:
        try:
            getattr(result, method)(show=False, save_path=str(out / "figures" / fname))
        except Exception as exc:  # noqa: BLE001
            print(f"no se pudo generar {fname}: {exc}")


if __name__ == "__main__":
    main()
