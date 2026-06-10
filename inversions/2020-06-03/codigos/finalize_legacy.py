#!/usr/bin/env python3
"""Cierra el estudio legacy-equiv con las corridas ya hechas: appraisal de
rep_seed11 + informe/dispersion con los recs disponibles."""
import sys, json
from pathlib import Path
import joblib
sys.path.insert(0, str(Path(__file__).resolve().parent))
import robustez_legacy as R

RESDIR = R.RESDIR
records = [json.loads(p.read_text()) for p in sorted(RESDIR.glob("rec_*.json"))]
RANGES = [(4, 13), (4, 13), (0, 2), (0, 1), (0, 1), (1, 6), (0.5, 3.5)]  # config legacy
print(f"Recs disponibles: {[r['label'] for r in records]}")

# Appraisal del de referencia (reusa el ensamble, no re-evalua)
try:
    ref = joblib.load(RESDIR / "res_rep_seed11.joblib")
    print("Corriendo appraisal de rep_seed11 ...")
    ref.run_appraisal(n_resample=20000, verbose=False, save=False)
    ref.plot_appraisal(save_path=str(R.OUT / "appraisal_corner.png"), show=False)
    print("  -> appraisal_corner.png")
except Exception as e:
    print("appraisal fallo:", e)

R.aggregate(records, RANGES)
print("LISTO")
