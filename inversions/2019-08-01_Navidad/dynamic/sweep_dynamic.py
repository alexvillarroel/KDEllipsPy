"""Barrido Te x tamaño de elipse del forward dinámico (Navidad NP1).

Reusa un solo DynamicNAInversionModel -> Green de AXITRA calculadas una vez.
Escribe sweep_results.csv (se va llenando) y sweep_best_traces.png.
"""

import itertools
import time

import numpy as np

import kdellipspy as kde
from kdellipspy.core.geometry import TSNFaultGridSpec
from kdellipspy.inversion.dynamic import DynamicNAInversionModel, TSNRunConfig
from kdellipspy.inversion.dynamic.tsn_bridge import read_tsn_fault_field

from forward_dynamic import DEFAULT_MODEL, DH, DT, HERE, NABC, CASE_DIR, NXTT, NZTT, WORK, setup_work_dir

TE_MPA = [3.0, 6.0, 10.0, 15.0]
SIZES = [(3.0, 2.0), (5.0, 3.0), (7.0, 4.0)]  # semiejes (a, b) en ptos de 2 km


def main():
    cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
    setup_work_dir(cfg)
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE_DIR / "input.ctl"), data_dir=str(CASE_DIR / "DATA"))
    inv = DynamicNAInversionModel(
        config=cfg, observed_waveforms=observed, time_array=time_array,
        tsn_run_cfg=TSNRunConfig(work_dir=WORK, nxtT=NXTT, nztT=NZTT, dt_s=DT),
        tsn_grid=TSNFaultGridSpec(dh=DH, dip_deg=cfg.source_position.dip, nztT=NZTT, nabc=NABC),
    )

    out = HERE / "sweep_results.csv"
    out.write_text("Te_MPa,a_km,b_km,misfit,slip_max_m,Mw,amp_ratio,seconds\n")
    best = (np.inf, None, None)
    for te, (a, b) in itertools.product(TE_MPA, SIZES):
        model = np.array(DEFAULT_MODEL, dtype=np.float32)
        model[0], model[1], model[5] = a, b, te
        t0 = time.time()
        misfit, syn = inv._evaluate_model(model)
        srz = read_tsn_fault_field(WORK / "result" / "sliprateZ.res", NXTT, NZTT)
        slip = srz.sum(axis=0) * DT
        m0 = inv.dynamic_fm._mu_pa * np.abs(slip).sum() * DH**2
        mw = (np.log10(m0) - 9.1) / 1.5 if m0 > 0 else float("nan")
        # sintético/observado en amplitud pico-a-pico global: >1 = sintético grande
        amp = np.ptp(syn, axis=-1).sum() / np.ptp(observed, axis=-1).sum() if syn is not None else float("nan")
        line = f"{te},{2*a:.0f},{2*b:.0f},{misfit:.4f},{np.abs(slip).max():.2f},{mw:.2f},{amp:.2f},{time.time()-t0:.0f}"
        print(line, flush=True)
        with out.open("a") as f:
            f.write(line + "\n")
        if misfit < best[0]:
            best = (misfit, model.copy(), syn)

    misfit, model, syn = best
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = [s.name for s in cfg.stations.stations]
    fig, axes = plt.subplots(len(names), 3, figsize=(10, 1.6 * len(names)), sharex=True)
    for i, name in enumerate(names):
        for c, comp in enumerate("XYZ"):
            ax = axes[i, c]
            ax.plot(time_array, observed[i, c], "k", lw=1)
            ax.plot(time_array, syn[i, c], "r", lw=1)
            ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(name, rotation=0, ha="right")
            if i == 0:
                ax.set_title(comp)
    fig.suptitle(f"Mejor del barrido: Te={model[5]:g} MPa, ejes {2*model[0]:g}x{2*model[1]:g} km, misfit {misfit:.3f}", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "sweep_best_traces.png", dpi=120)
    inv.clean()


if __name__ == "__main__":
    main()
