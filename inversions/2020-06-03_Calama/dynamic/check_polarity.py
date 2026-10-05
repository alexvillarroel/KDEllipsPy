"""Verificación de polaridad del dinámico en una falla NORMAL (guía, sección 5.1).

fd3d_TSN no recibe el rake: solo el buzamiento. Dentro de la aspereza el
prestress es +Te, así que con DIPSLIP el deslizamiento va SIEMPRE hacia +Z
local, es decir rake +90 (inverso). Para Calama 2020 (np1 rake -85, np2 rake
-101) eso daría sintéticos con la polaridad invertida y una inversión
completa tirada a la basura.

Este script corre UN forward de fd3d_TSN y arma los sintéticos con las dos
convenciones (el mismo slip, solo cambia el signo del vector de momento):

    slip_sign = +1  ->  lo que hacía el código de Navidad (rake +90, inverso)
    slip_sign = -1  ->  rake + 180 = -90, el sentido normal

y las compara con los observados: misfit, rake efectivo que sale de
project_rake, y correlación por estación en la ventana P.

NO se usa el signo del primer arribo P: los dos sintéticos son negativos
exactos uno del otro, así que los conteos de coincidencia suman siempre el
número de estaciones y la métrica no puede discriminar (en np1 dio 4/9 vs 5/9
mientras la cc daba +0.67 vs +0.87). La evidencia es el rake efectivo y la cc.

Cada signo se evalúa en SU mejor dt0 (barrido), no en dt0 = 0: con un desfase
de hora de origen el misfit satura en ~1 y la comparación no distingue nada
(es el mismo motivo por el que el preflight de na_dynamic.py barre dt0).

Uso:  DYN_CASE=dyn_cases/np1 python check_polarity.py
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kdellipspy as kde
from kdellipspy.inversion.dynamic import DynamicNAInversionModel
from kdellipspy.inversion.dynamic.dynamic_convolution import (
    bin_slip_rate_to_subfaults, project_rake,
)

from forward_dynamic import (CASE, CASE_DIR, DEFAULT_MODEL, DH, HERE, NXTT, NZTT,
                             setup_work_dir, tsn_grid, tsn_run_cfg)

DT0_SCAN = np.arange(-3.0, 10.1, 0.5)  # s, corrección de hora de origen


def main():
    cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
    setup_work_dir(cfg)
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE_DIR / "input.ctl"), data_dir=str(CASE_DIR / "DATA"))
    azt = np.asarray(kde.build_azi_times_array(config=cfg, model_name="iasp91"), float)

    inv = DynamicNAInversionModel(config=cfg, observed_waveforms=observed, time_array=time_array,
                                  tsn_run_cfg=tsn_run_cfg(), tsn_grid=tsn_grid(cfg))
    rake_ctl = float(cfg.source_position.rake)
    print(f"[caso] {CASE}  strike {cfg.source_position.strike:g} / dip {cfg.source_position.dip:g} "
          f"/ rake {rake_ctl:g}   slip_sign automático = {inv.dynamic_fm.slip_sign:+.0f}", flush=True)

    model = np.array(DEFAULT_MODEL, dtype=np.float32)
    # _evaluate_model devuelve (misfit, sintéticos) de ESA evaluación;
    # inv.best_synthetics es el mejor hasta ahora y no sirve para comparar.
    # fd3d se cachea por modelo, así que todo el barrido cuesta un solo run.
    out = {}
    for sign in (+1.0, -1.0):
        inv.dynamic_fm.slip_sign = sign
        best = None
        for dt0 in DT0_SCAN:
            inv.dynamic_fm.time_shift_s = float(dt0)
            misfit, syn = inv._evaluate_model(model)
            if syn is not None and (best is None or misfit < best[0]):
                best = (misfit, np.asarray(syn, float), float(dt0))
        assert best is not None, f"ningún forward válido con slip_sign {sign:+.0f}"
        out[sign] = best
        print(f"[slip_sign {sign:+.0f}] mejor misfit = {best[0]:.4f} en dt0 = {best[2]:+.1f} s", flush=True)

    # Rake efectivo que sale del slip de fd3d, con y sin el flip.
    from kdellipspy.inversion.dynamic.tsn_bridge import read_tsn_fault_field
    work = Path(tsn_run_cfg().work_dir) / "result"
    srx = read_tsn_fault_field(work / "sliprateX.res", NXTT, NZTT)
    srz = read_tsn_fault_field(work / "sliprateZ.res", NXTT, NZTT)
    mx, mz = bin_slip_rate_to_subfaults(srx, srz, cfg.fault_plane.nx, cfg.fault_plane.ny,
                                        DH, inv.dynamic_fm._mu_pa)
    w = np.hypot(mx.sum(axis=1), mz.sum(axis=1))  # pondera por momento de la subfalla
    active = w > 0.01 * w.max()
    rakes = {}
    for sign in (+1.0, -1.0):
        _, rk = project_rake(sign * mx, sign * mz)
        rakes[sign] = float(np.average(rk[active], weights=w[active]))
        print(f"[slip_sign {sign:+.0f}] rake efectivo (medio pesado por M0) = {rakes[sign]:+7.1f}°  "
              f"| rake del input.ctl = {rake_ctl:+.1f}°", flush=True)

    # Polaridad del primer arribo P en la vertical y correlación, por estación.
    names = [s.name for s in cfg.stations.stations]
    dt = float(time_array[1] - time_array[0])
    print(f"\n{'sta':<7s}{'tP(s)':>7s}{'cc+1':>8s}{'cc-1':>8s}{'mejor':>7s}", flush=True)
    rows = []
    for j, n in enumerate(names):
        tp = float(azt[j, 1])
        i0, i1 = max(0, int(tp / dt)), min(len(time_array), int((tp + 25.0) / dt))
        cc = {s: float(np.corrcoef(observed[j, 2, i0:i1], out[s][1][j, 2, i0:i1])[0, 1])
              for s in (+1.0, -1.0)}
        rows.append((n, cc[+1.0], cc[-1.0]))
        print(f"{n:<7s}{tp:7.1f}{cc[+1.0]:8.3f}{cc[-1.0]:8.3f}"
              f"{('-1' if cc[-1.0] > cc[+1.0] else '+1'):>7s}", flush=True)

    ccm = {s: float(np.nanmean([r[1 if s > 0 else 2] for r in rows])) for s in (+1.0, -1.0)}
    nwin = sum(1 for r in rows if r[2] > r[1])
    print(f"\ncc media (ventana P, 25 s): +1 -> {ccm[+1.0]:+.3f}   -1 -> {ccm[-1.0]:+.3f}"
          f"   (el flip gana en {nwin}/{len(rows)} estaciones)", flush=True)
    winner = -1.0 if (out[-1.0][0] < out[+1.0][0]) else +1.0
    print(f"\n==> gana slip_sign {winner:+.0f} "
          f"(misfit {out[winner][0]:.4f} en dt0 {out[winner][2]:+.1f} s  vs  "
          f"{out[-winner][0]:.4f} en dt0 {out[-winner][2]:+.1f} s); "
          f"automático = {np.sign(np.sin(np.radians(rake_ctl))):+.0f}", flush=True)

    fig, axes = plt.subplots(len(names), 3, figsize=(11, 1.5 * len(names)), sharex=True)
    for i, n in enumerate(names):
        for c, comp in enumerate("NEZ"):
            ax = axes[i, c]
            ax.plot(time_array, observed[i, c], "k", lw=1.2, label="Observed" if i == 0 and c == 0 else None)
            ax.plot(time_array, out[+1.0][1][i, c], color="tab:orange", lw=1, ls="--",
                    label=f"fd3d as-is (rake +90, reverse), dt0 {out[+1.0][2]:+.1f} s"
                    if i == 0 and c == 0 else None)
            ax.plot(time_array, out[-1.0][1][i, c], color="tab:red", lw=1,
                    label=f"flipped (rake +180, normal), dt0 {out[-1.0][2]:+.1f} s"
                    if i == 0 and c == 0 else None)
            ax.axvline(float(azt[i, 1]), color="0.6", lw=0.8)
            ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(n, rotation=0, ha="right")
            if i == 0:
                ax.set_title(comp)
            if i == len(names) - 1:
                ax.set_xlabel("Time (s)")
    fig.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"Calama 2020 · {CASE} (rake {rake_ctl:g}°, normal) · dynamic slip-direction check\n"
                 f"misfit: as-is {out[+1.0][0]:.3f} / flipped {out[-1.0][0]:.3f}   "
                 f"effective rake: {rakes[+1.0]:+.0f}° / {rakes[-1.0]:+.0f}°", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = HERE / f"check_polarity_{CASE.replace('/', '__')}.png"
    fig.savefig(png, dpi=130)
    print(f"\nfigura: {png}", flush=True)
    inv.clean()


if __name__ == "__main__":
    main()
