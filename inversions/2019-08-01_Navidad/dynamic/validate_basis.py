"""Valida el camino rápido (operador de axitra precalculado por subfalla) contra el lento
(una conv() de axitra por subfalla) con el mejor modelo del NA en NP1."""
import time
import joblib
import numpy as np
import kdellipspy as kde
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicForwardModel, TSNRunConfig
from forward_dynamic import DH, DT, NABC, CASE_DIR, NXTT, NZTT, WORK, setup_work_dir
from na_dynamic import to_full_model

cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
setup_work_dir(cfg)
fm = DynamicForwardModel(cfg, TSNRunConfig(work_dir=WORK, nxtT=NXTT, nztT=NZTT, dt_s=DT),
                         TSNFaultGridSpec(dh=DH, dip_deg=cfg.source_position.dip, nztT=NZTT, nabc=NABC))
model = to_full_model(joblib.load("na_output/na_result.joblib").best_model.model)

fm.use_basis = False
t0 = time.time(); slow = fm.forward(model); print(f"slow eval: {time.time()-t0:.0f} s", flush=True)
fm.use_basis = True
t0 = time.time(); fm.ensure_basis(int(cfg.observed_data.units)); print(f"basis (once): {time.time()-t0:.0f} s", flush=True)
t0 = time.time(); fast = fm.forward(model); print(f"fast eval: {time.time()-t0:.0f} s", flush=True)
rel = np.abs(fast - slow).max() / np.abs(slow).max()
cc = np.corrcoef(fast.ravel(), slow.ravel())[0, 1]
print(f"max rel err {rel:.2e}   corr {cc:.6f}   |slow| {np.abs(slow).max():.3e}")
np.savez("validate_basis.npz", slow=slow, fast=fast)
fm.clean()
