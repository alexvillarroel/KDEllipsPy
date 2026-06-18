#!/usr/bin/env python3
"""Appraisal del resultado (run limpio tiempo 07:35:34, dt=1/128)."""
import numpy as np, joblib
from pathlib import Path
import kdellipspy as kde
B = Path(__file__).resolve().parent.parent
r = joblib.load(B / "output" / "inversion_result.joblib")
b = np.asarray(r.best_model.model, float)
print("misfit =", round(r.best_model.misfit, 4))
print("model  =", dict(zip(["a1", "a2", "theta", "np", "tp", "dmax", "vr"], np.round(b, 3))))
try:
    r.run_appraisal(n_resample=20000, verbose=False, save=False)
    r.plot_appraisal(save_path=str(B / "output" / "appraisal_corner_t34.png"), show=False)
    print("appraisal -> output/appraisal_corner_t34.png")
except Exception as e:
    print("appraisal fallo:", e)
