"""Diagnóstico: ¿el desajuste es polaridad o desfase? (Te=3 MPa, 14x8 km)."""
import numpy as np
import kdellipspy as kde
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicNAInversionModel, TSNRunConfig
from forward_dynamic import DEFAULT_MODEL, DH, DT, HERE, NABC, CASE_DIR, NXTT, NZTT, WORK

cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
obs, t = kde.load_and_filter_observed_data(input_ctl_path=str(CASE_DIR / "input.ctl"), data_dir=str(CASE_DIR / "DATA"))
inv = DynamicNAInversionModel(config=cfg, observed_waveforms=obs, time_array=t,
    tsn_run_cfg=TSNRunConfig(work_dir=WORK, nxtT=NXTT, nztT=NZTT, dt_s=DT),
    tsn_grid=TSNFaultGridSpec(dh=DH, dip_deg=cfg.source_position.dip, nztT=NZTT, nabc=NABC))
m = np.array(DEFAULT_MODEL, dtype=np.float32); m[0], m[1], m[5] = 7, 4, 3
mis, syn = inv._evaluate_model(m)
np.savez(HERE / "diag_phase.npz", obs=obs, syn=syn, t=t)
l2 = lambda s: ((obs - s) ** 2).sum() / (obs ** 2).sum()
print(f"misfit pipeline {mis:.3f} | L2 +syn {l2(syn):.3f} | L2 -syn {l2(-syn):.3f}")
for lag in range(-10, 11, 2):
    print(f"lag {lag:+d} s: L2 {l2(np.roll(syn, lag, axis=-1)):.3f}")
for i, s in enumerate(cfg.stations.stations):
    cc = [np.corrcoef(obs[i, c], syn[i, c])[0, 1] for c in range(3)]
    print(s.name, " ".join(f"{c:+.2f}" for c in cc))
inv.clean()
