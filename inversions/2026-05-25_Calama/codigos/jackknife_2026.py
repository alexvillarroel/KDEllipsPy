#!/usr/bin/env python3
"""
jackknife_2026.py
=================
Jackknife de estaciones/componentes para la inversion cinematica de Calama 2026.
Corre una matriz de inversiones NA (en memoria, sin tocar input.ctl/DATA) variando
que estaciones/componentes estan activas, y rankea por:

    score = misfit + LAMBDA * |Mw - MW_REF|        (MW_REF = 6.9, USGS)

Matriz (build_matrix):
  * all_in                : las 14 estaciones activas (R,T,Z).
  * jack_no_<EST>         : deja fuera una estacion a la vez (flags 0 0 0).
  * jack_<EST>_no{R,T,Z}  : quita una componente en estaciones sospechosas.

Adaptado de inversions/2020-06-03_Calama/codigos/robustez_2020.py. Reusa la MISMA
ruta de carga que run_inversion.py (kde.load_and_filter_observed_data): los datos
observados se cargan UNA vez (banda y estaciones fijas; solo cambian las flags).

Salida en jackknife/:  results/res_<label>.joblib + rec_<label>.json,
ranking.csv, INFORME_jackknife.md, dispersion.png.

Uso:
    python codigos/jackknife_2026.py --check          # valida la matriz (sin NA)
    python codigos/jackknife_2026.py --smoke          # end-to-end minimo (~1-2 min)
    python codigos/jackknife_2026.py --n-jobs 4       # corrida completa (de noche)
"""
import sys, json, time, argparse, csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import joblib

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "codigos"))
import kdellipspy as kde
try:
    from kdellipspy import NAConfig
except Exception:
    from kdellipspy.inversion.kinematic.model_na import NAConfig

CTL = BASE / "input.ctl"
DATA = BASE / "DATA"
OUT = BASE / "jackknife"; OUT.mkdir(exist_ok=True)
RESDIR = OUT / "results"; RESDIR.mkdir(exist_ok=True)

PNAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]
REF_BAND = (0.06, 0.15)
REF_SEED = 11
MW_REF = 6.9      # solo referencia para reportar el desvio; NO penaliza
LAMBDA = 0.0      # 0 => el objetivo es solo el misfit (Mw = validacion a posteriori)

# 14 estaciones (orden = input.ctl). Se lee del ctl para no duplicar.
STATIONS = [s.name for s in kde.ConfigParser(filepath=str(CTL)).stations.stations]

# Drops de componente curados en estaciones sospechosas (editable).
# A22F = reemplazo del outlier AF01 (azimut este); PB20 = deconvuelta a mano;
# A20F = muy cercana (15.9 km, campo cercano).
DROP_IDX = {"noR": 0, "noT": 1, "noZ": 2}   # use_n=R, use_e=T, use_z=Z
SUSPECT_DROPS = [(st, d) for st in ("A22F", "PB20", "A20F")
                 for d in ("noR", "noT", "noZ")]


def make_na(smoke, seed, n_jobs):
    if smoke:
        return NAConfig(n_samples_initial=20, n_samples_iteration=10, n_iterations=2,
                        n_cells_resample=2, n_jobs=n_jobs, random_seed=seed)
    return NAConfig(n_samples_initial=500, n_samples_iteration=50, n_iterations=50,
                    n_cells_resample=15, n_jobs=n_jobs, random_seed=seed)


def sel_all():
    return {s: (1, 1, 1) for s in STATIONS}


def build_matrix():
    """Lista de (label, sel) — sel = {estacion: (use_n,use_e,use_z)}."""
    configs = [("all_in", sel_all())]
    for st in STATIONS:                       # leave-one-out por estacion
        sel = sel_all(); sel[st] = (0, 0, 0)
        configs.append((f"jack_no_{st}", sel))
    for st, drop in SUSPECT_DROPS:            # drop de una componente
        sel = sel_all()
        flags = list(sel[st]); flags[DROP_IDX[drop]] = 0
        sel[st] = tuple(flags)
        configs.append((f"jack_{st}_{drop}", sel))
    return configs


def run_one(label, sel, observed, time_array, azt, seed, n_jobs, smoke):
    t0 = time.time()
    cfg = kde.ConfigParser(filepath=str(CTL))
    cfg.ellipse.freq1, cfg.ellipse.freq2 = REF_BAND
    for s in cfg.stations.stations:
        un, ue, uz = sel.get(s.name, (1, 1, 1))
        s.use_n, s.use_e, s.use_z = bool(un), bool(ue), bool(uz)
    model = kde.NAInversionModel(config=cfg, observed_waveforms=observed,
                                 time_array=time_array, azi_times_array=azt)
    model.use_green_cache = True
    res = model.run_na_search(make_na(smoke, seed, n_jobs))
    b = np.asarray(res.best_model.model, float)
    m0, mw = model.fm.estimate_total_moment_and_mw(b)
    misfit = float(res.best_model.misfit)
    score = misfit + LAMBDA * abs(mw - MW_REF)   # LAMBDA=0 por defecto => score==misfit
    try:
        model.clear_green_cache()
    except Exception:
        pass
    rec = dict(label=label, seed=int(seed), misfit=misfit, mw=float(mw),
               m0=float(m0), score=float(score),
               sel={k: list(v) for k, v in sel.items()},
               model={n: float(x) for n, x in zip(PNAMES, b)},
               secs=round(time.time() - t0, 1))
    joblib.dump(res, RESDIR / f"res_{label}.joblib")
    (RESDIR / f"rec_{label}.json").write_text(json.dumps(rec, indent=2))
    print(f"[OK] {label:18s} misfit={misfit:.4f} Mw={mw:.2f} "
          f"(dMw={mw-MW_REF:+.2f}) ({rec['secs']:.0f}s)", flush=True)
    return rec


def aggregate(records):
    # Objetivo = solo misfit. Mw se reporta como validacion (desvio vs MW_REF).
    rows = sorted(records, key=lambda r: r["misfit"])
    # ranking.csv
    with (OUT / "ranking.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "misfit", "Mw", "dMw_vs_ref", *PNAMES])
        for r in rows:
            w.writerow([r["label"], f"{r['misfit']:.4f}", f"{r['mw']:.2f}",
                        f"{r['mw']-MW_REF:+.2f}",
                        *[f"{r['model'][p]:.3f}" for p in PNAMES]])
    # informe markdown
    base = next((r for r in records if r["label"] == "all_in"), None)
    L = ["# Jackknife Calama 2026 — ranking por misfit (Mw = validacion)\n",
         f"Objetivo = solo misfit. Mw_ref={MW_REF} solo para reportar desvio. "
         f"Corridas: {len(records)}.\n",
         "| # | corrida | misfit | Mw | ΔMw vs ref | Δmisfit vs all_in |",
         "|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        dm = "" if base is None else f"{r['misfit']-base['misfit']:+.4f}"
        L.append(f"| {i} | {r['label']} | {r['misfit']:.4f} | {r['mw']:.2f} | "
                 f"{r['mw']-MW_REF:+.2f} | {dm} |")
    L += ["\n**Lectura:** el mejor candidato es el de menor misfit. Si quitar una "
          "estacion (`jack_no_<EST>`) dispara `vr` o mueve el parche respecto a "
          "`all_in`, esa estacion era el ancla fisica (leverage). El Mw se valida "
          f"vs {MW_REF} pero NO se fuerza.\n"]
    (OUT / "INFORME_jackknife.md").write_text("\n".join(L))

    # dispersion de parametros (normalizada a su rango observado)
    fig, ax = plt.subplots(figsize=(11, 5))
    data = [np.array([r["model"][p] for r in records]) for p in PNAMES]
    norm = [(d - d.min()) / (d.max() - d.min()) if d.max() > d.min() else d*0
            for d in data]
    ax.boxplot(norm, tick_labels=PNAMES, showmeans=True)
    ax.set_ylabel("valor normalizado al rango observado [0,1]")
    ax.set_title(f"Jackknife 2026: dispersion de parametros ({len(records)} corridas)")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(OUT / "dispersion.png", dpi=150)
    print("\n".join(L[:24]))
    print(f"\n-> {OUT}/ranking.csv  +  INFORME_jackknife.md  +  dispersion.png")


def main():
    global LAMBDA, MW_REF
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="NA minimo, 2 configs")
    ap.add_argument("--check", action="store_true", help="solo valida la matriz")
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--lambda", dest="lam", type=float, default=LAMBDA)
    ap.add_argument("--mw-ref", type=float, default=MW_REF)
    args = ap.parse_args()
    LAMBDA, MW_REF = args.lam, args.mw_ref

    if args.check:
        _check(); return

    configs = build_matrix()
    if args.smoke:
        configs = [configs[0], configs[1]]   # all_in + primer jackknife
    print(f"=== JACKKNIFE 2026 — {len(configs)} corridas "
          f"(smoke={args.smoke}, n_jobs={args.n_jobs}) ===", flush=True)

    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CTL), data_dir=str(DATA),
        prefer_raw=False, freq1=REF_BAND[0], freq2=REF_BAND[1])
    cfg0 = kde.ConfigParser(filepath=str(CTL))
    azt = np.asarray(kde.build_azi_times_array(config=cfg0, model_name="iasp91"), float)
    print(f"observed {observed.shape} | azi_times {azt.shape}", flush=True)

    records = []
    for label, sel in configs:
        records.append(run_one(label, sel, observed, time_array, azt,
                               REF_SEED, args.n_jobs, args.smoke))
    aggregate(records)
    print("\n=== JACKKNIFE COMPLETO ===", flush=True)


def _check():
    """Valida la logica de la matriz sin correr NA."""
    cfgs = dict(build_matrix())
    assert cfgs["all_in"] == sel_all(), "all_in deberia tener todo (1,1,1)"
    # leave-one-out: la estacion fuera en (0,0,0), el resto intacto
    for st in STATIONS:
        sel = cfgs[f"jack_no_{st}"]
        assert sel[st] == (0, 0, 0), f"{st} deberia estar fuera"
        assert all(sel[o] == (1, 1, 1) for o in STATIONS if o != st)
    # drop de componente: solo esa flag a 0
    sel = cfgs["jack_A22F_noT"]
    assert sel["A22F"] == (1, 0, 1), "noT deberia apagar solo use_e (T)"
    # objetivo = solo misfit: con LAMBDA=0, score == misfit (Mw no penaliza)
    assert LAMBDA == 0.0, "LAMBDA por defecto debe ser 0 (Mw es validacion, no objetivo)"
    n = 1 + len(STATIONS) + len(SUSPECT_DROPS)
    assert len(cfgs) == n, f"esperaba {n} configs, hay {len(cfgs)}"
    print(f"check OK: {len(cfgs)} configs (1 all_in + {len(STATIONS)} jackknife "
          f"+ {len(SUSPECT_DROPS)} drops); objetivo=misfit, Mw=validacion")


if __name__ == "__main__":
    main()
