#!/usr/bin/env python3
"""
misfit_breakdown.py
===================
Desglose del misfit L2  E = Sum(obs-syn)^2 / Sum(obs^2)  por estacion y
componente (rotado a R/T/Z). Funciona en dos modos:
  - señal completa (window_s=None / 0): usa toda la traza.
  - ventana (window_s>0): P -> R,Z desde tP ; S -> T desde tS  (ventana de
    window_s segundos), igual que la rama use_full_signal=False del misfit.

Expone ``breakdown_from_result(result, window_s="auto")`` reutilizable por
run_inversion.py; el modo se toma de la ventana del input.ctl si es "auto".

Uso (CLI):
    python codigos/misfit_breakdown.py
"""

from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
JOBLIB = BASE / "output" / "inversion_result.joblib"


def _num_den(obs, syn, time, azt, flags, window_s):
    """Devuelve (num, den) shape (nsta,3) en componentes rotadas R,T,Z."""
    nsta, _, npts = obs.shape
    azi, tP, tS = azt[:, 0], azt[:, 1], azt[:, 2]
    num = np.zeros((nsta, 3)); den = np.zeros((nsta, 3))
    use_win = window_s is not None and window_s > 0
    if use_win:
        dt = float(time[1] - time[0])
        sampling = max(1, int(np.rint(1.0 / dt)))
        win = max(1, int(np.rint(window_s * sampling)))

    for j in range(nsta):
        az = float(azi[j])
        un, ue, uz = flags[j]
        # indices de ventana (o señal completa)
        if use_win:
            sp = int(np.rint(tP[j])) * sampling
            ss = int(np.rint(tS[j])) * sampling
            kp0, kp1 = max(0, sp - 1), min(npts - 1, sp + win - 1)
            ks0, ks1 = max(0, ss - 1), min(npts - 1, ss + win - 1)
        else:
            kp0, kp1 = 0, npts - 1
            ks0, ks1 = 0, npts - 1
        # P window -> R, Z
        if kp1 >= kp0:
            xo, yo, zo = obs[j, 0, kp0:kp1+1], obs[j, 1, kp0:kp1+1], obs[j, 2, kp0:kp1+1]
            xs, ys, zs = syn[j, 0, kp0:kp1+1], syn[j, 1, kp0:kp1+1], syn[j, 2, kp0:kp1+1]
            if un:
                ro = xo*np.cos(az)+yo*np.sin(az); rs = xs*np.cos(az)+ys*np.sin(az)
                num[j, 0] = np.sum((ro-rs)**2); den[j, 0] = np.sum(ro**2)
            if uz:
                num[j, 2] = np.sum((zo-zs)**2); den[j, 2] = np.sum(zo**2)
        # S window -> T
        if ks1 >= ks0:
            xo, yo = obs[j, 0, ks0:ks1+1], obs[j, 1, ks0:ks1+1]
            xs, ys = syn[j, 0, ks0:ks1+1], syn[j, 1, ks0:ks1+1]
            if ue:
                to = yo*np.cos(az)-xo*np.sin(az); ts = ys*np.cos(az)-xs*np.sin(az)
                num[j, 1] = np.sum((to-ts)**2); den[j, 1] = np.sum(to**2)
    return num, den


def breakdown_from_result(result, window_s="auto"):
    """Devuelve un string con la tabla de desglose del misfit del mejor modelo."""
    obs = np.asarray(result.observed, float)
    syn = np.asarray(result.best_synthetics, float)
    time = np.asarray(result.time, float)
    cfg = result.config
    names = [s.name for s in cfg.stations.stations]
    flags = np.array([[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations], bool)
    azt = getattr(result, "azi_times_array", None)
    if azt is None:
        azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")
    azt = np.asarray(azt, float)

    if window_s == "auto":
        try:
            tw = float(cfg.inversion_process.misfit_time_window)
            window_s = tw if tw > 0 else None
        except Exception:
            window_s = None

    num, den = _num_den(obs, syn, time, azt, flags, window_s)
    num_tot, den_tot = num.sum(), den.sum()
    E = num_tot / den_tot if den_tot > 0 else float("nan")
    comps = ["R", "T", "Z"]
    mode = f"ventana {window_s:.0f}s (P->R,Z ; S->T)" if window_s else "señal completa"

    L = []
    L.append(f"Desglose de misfit ({mode}) — E = {E:.6f}")
    L.append("=" * 78)
    L.append(f"{'Est':<6}{'Cmp':<4}{'%resid':>9}{'aporte_E':>10}{'mf_local':>10}{'peso_obs':>10}")
    L.append("-" * 78)
    for j in np.argsort(num.sum(axis=1))[::-1]:
        for c in range(3):
            pres = 100*num[j, c]/num_tot if num_tot else 0
            ap = num[j, c]/den_tot if den_tot else 0
            mfl = num[j, c]/den[j, c] if den[j, c] > 0 else 0.0
            peso = 100*den[j, c]/den_tot if den_tot else 0
            L.append(f"{names[j]:<6}{comps[c]:<4}{pres:>8.1f}%{ap:>10.4f}{mfl:>10.3f}{peso:>9.1f}%")
        L.append(f"{'  -> '+names[j]+' total':<28}{100*num[j].sum()/num_tot:>6.1f}% del residuo"
                 f"  (peso obs {100*den[j].sum()/den_tot:.1f}%)")
        L.append("-" * 78)
    L.append("Por componente:")
    for c in range(3):
        L.append(f"  {comps[c]}: {100*num[:,c].sum()/num_tot:5.1f}% del residuo   "
                 f"(peso obs {100*den[:,c].sum()/den_tot:4.1f}%)")
    return "\n".join(L)


def main():
    result = joblib.load(JOBLIB)
    print(breakdown_from_result(result, window_s="auto"))


if __name__ == "__main__":
    main()
