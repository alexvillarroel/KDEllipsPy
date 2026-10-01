"""Diagnóstico: mismo modelo (mejor NA de NP1, misfit 0.652 por camino lento),
camino rápido, en np1 y loc_mid. Directorio de trabajo propio."""
import os, sys, time
import numpy as np
import kdellipspy as kde
import forward_dynamic as fd
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicNAInversionModel, TSNRunConfig

fd.WORK = fd.HERE / "tsn_work_diag"
MODEL = np.array([4.879, 3.571, 12.204, 10.127, 0.941, 3.029, 1.15, 1.1, 1.5, 1.0], np.float32)
for case in sys.argv[1:]:
    d = fd.HERE.parent / case
    cfg = kde.ConfigParser(str(d / "input.ctl")); fd.setup_work_dir(cfg)
    obs, t = kde.load_and_filter_observed_data(input_ctl_path=str(d / "input.ctl"), data_dir=str(d / "DATA"))
    inv = DynamicNAInversionModel(config=cfg, observed_waveforms=obs, time_array=t,
        tsn_run_cfg=TSNRunConfig(work_dir=fd.WORK, nxtT=fd.NXTT, nztT=fd.NZTT, dt_s=fd.DT),
        tsn_grid=TSNFaultGridSpec(dh=fd.DH, dip_deg=cfg.source_position.dip, nztT=fd.NZTT, nabc=fd.NABC))
    t0 = time.time(); mis, syn = inv._evaluate_model(MODEL)
    amp = np.ptp(syn, -1).sum() / np.ptp(obs, -1).sum()
    cc = np.corrcoef(syn.ravel(), obs.ravel())[0, 1]
    print(f"{case}: misfit {mis:.4f}  amp {amp:.3f}  corr {cc:+.3f}  ({time.time()-t0:.0f} s)", flush=True)
    np.save(fd.HERE / f"diag_syn_{case}.npy", syn)
    inv.clean()
