"""barrido_fuente.py — sensibilidad del misfit a (source_type, T0).

Fija la mejor geometría (del joblib) y reevalúa el misfit para varias funciones
fuente de Axitra y tiempos de rise T0. Rápido: las Green se cachean (solo corre
conv por evaluación). NO re-invierte; es un sondeo de sensibilidad.

    conda run -n kdellipspy python codigos/barrido_fuente.py
"""
from pathlib import Path

import kdellipspy as kde
from kdellipspy.inversion.base import NAResult

BASE = Path(__file__).resolve().parent.parent

cfg = kde.ConfigParser(filepath=str(BASE / "input.ctl"))
observed, time_array = kde.load_and_filter_observed_data(
    input_ctl_path=str(BASE / "input.ctl"), data_dir=str(BASE / "DATA"),
    freq1=cfg.ellipse.freq1, freq2=cfg.ellipse.freq2)
azi = kde.build_azi_times_array(config=cfg, model_name="iasp91")
model = kde.NAInversionModel(config=cfg, observed_waveforms=observed,
                             time_array=time_array, azi_times_array=azi)
model.use_green_cache = True

best = NAResult.load(str(BASE / "output/inversion_result.joblib")).best_model.model

SRC = {1: "Ricker", 2: "smooth step", 4: "triangulo", 5: "rampa", 7: "true step"}
T0S = [1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
cur_st, cur_t0 = int(cfg.ellipse.source_type), float(cfg.ellipse.t0)
print(f"geometría fija (mejor modelo): a1={best[0]:.2f} a2={best[1]:.2f} "
      f"Vr={best[6]:.2f}  ·  actual: src={cur_st} T0={cur_t0:g}\n")
print(f"{'source':>13s} | " + "  ".join(f"T0={t:<4g}" for t in T0S))
print("-" * 70)
best_cell = (1e9, None)
for st, name in SRC.items():
    cfg.ellipse.source_type = st
    row = []
    for t0 in T0S:
        cfg.ellipse.t0 = t0
        mf, _ = model._evaluate_model(best)
        row.append(mf)
        if mf < best_cell[0]:
            best_cell = (mf, (name, t0))
    print(f"{name:>13s} | " + "  ".join(f"{m:6.4f}" for m in row))
print("-" * 70)
print(f"mejor: {best_cell[1][0]} con T0={best_cell[1][1]:g}  ->  misfit {best_cell[0]:.4f}")
