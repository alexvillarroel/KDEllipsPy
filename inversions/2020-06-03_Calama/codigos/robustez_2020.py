#!/usr/bin/env python3
"""
robustez_2020.py
================
Estudio de ROBUSTEZ de la inversion cinematica del 2020. Corre una matriz de
inversiones (en memoria, sin tocar input.ctl/real_disp.py) variando:

  * REPETIBILIDAD: misma config, 5 semillas distintas  -> ¿el mejor es estable?
  * JACKKNIFE: deja fuera una estacion a la vez         -> ¿quien lo controla?
  * BANDAS: 0.02-0.1, 0.04-0.12, 0.06-0.15 Hz           -> ¿persiste en frecuencia?
  * APPRAISAL bayesiano del caso de referencia          -> marginales + correlaciones

Para cada corrida guarda el NAResult (joblib) y agrega un informe de robustez:
dispersion por parametro (mediana + rango intercuartil), figuras y un .md/.json.

Reusa la convencion de datos del proyecto (real_disp.procesar_estacion: SAC/DISP,
filtro causal, ventana 0-128 s, 256 muestras, orden N,E,Z) filtrando por banda.

Uso:
    python codigos/robustez_2020.py            # matriz completa (~5-6 h)
    python codigos/robustez_2020.py --smoke    # prueba rapida (NA minimo, 1-2 min)
"""
import sys, json, time, argparse
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import joblib

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "codigos"))
import kdellipspy as kde
import real_disp as rd                      # reuso de procesar_estacion / SELECTED / DISP_DIR
from event import load_event
try:
    from kdellipspy import NAConfig
except Exception:
    from kdellipspy.inversion.kinematic.model_na import NAConfig

OUT = BASE / "robustez"; OUT.mkdir(exist_ok=True)
RESDIR = OUT / "results"; RESDIR.mkdir(exist_ok=True)

PNAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]
STATIONS = list(rd.SELECTED)                  # orden de estaciones (= input.ctl)
EV = load_event(str(BASE / "event.ctl"))

# Seleccion de componentes de referencia (flags: use_n=R, use_e=T, use_z=Z)
REFSEL = {s: (1, 1, 1) for s in STATIONS}
REFSEL["A08F"] = (0, 0, 0)   # fuera
REFSEL["A22F"] = (0, 1, 1)   # sin R
REFSEL["A05F"] = (1, 1, 0)   # sin Z
REFSEL["A23F"] = (1, 0, 1)   # sin T

REF_BAND = (0.06, 0.15)
REF_SEED = 11


def load_observed(band):
    """Carga las 9 estaciones (SAC/DISP) filtradas a la banda dada -> (9,3,npts) N,E,Z."""
    rd.F1, rd.F2 = float(band[0]), float(band[1])
    obs = []
    for s in STATIONS:
        r = rd.procesar_estacion(s, EV, rd.DISP_DIR, "disp")
        if r is None:
            raise RuntimeError(f"falta {s}")
        n, e, z, la, lo = r
        obs.append([n, e, z])
    return np.asarray(obs, float)


def make_na(smoke, seed):
    if smoke:
        return NAConfig(n_samples_initial=20, n_samples_iteration=10, n_iterations=2,
                        n_cells_resample=2, n_jobs=1, random_seed=seed)
    return NAConfig(n_samples_initial=500, n_samples_iteration=50, n_iterations=50,
                    n_cells_resample=10, n_jobs=1, random_seed=seed)


def run_one(label, band, sel, seed, smoke=False):
    t0 = time.time()
    cfg = kde.ConfigParser(filepath=str(BASE / "input.ctl"))
    cfg.ellipse.freq1, cfg.ellipse.freq2 = float(band[0]), float(band[1])
    for s in cfg.stations.stations:
        un, ue, uz = sel.get(s.name, (1, 1, 1))
        s.use_n, s.use_e, s.use_z = bool(un), bool(ue), bool(uz)
    obs = load_observed(band)
    t = np.arange(rd.NPTS) * rd.DT
    azt = np.asarray(kde.build_azi_times_array(config=cfg, model_name="iasp91"), float)
    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs,
                                 time_array=t, azi_times_array=azt)
    model.use_green_cache = True
    res = model.run_na_search(make_na(smoke, seed))
    try:
        model.clear_green_cache()             # libera el axi_*.res de esta corrida
    except Exception:
        pass
    b = np.asarray(res.best_model.model, float)
    rec = dict(label=label, band=list(band), seed=int(seed),
               sel={k: list(v) for k, v in sel.items()},
               misfit=float(res.best_model.misfit),
               model={n: float(x) for n, x in zip(PNAMES, b)},
               secs=round(time.time() - t0, 1))
    joblib.dump(res, RESDIR / f"res_{label}.joblib")
    (RESDIR / f"rec_{label}.json").write_text(json.dumps(rec, indent=2))
    ranges = list(getattr(model, "param_ranges", []))
    print(f"[OK] {label:22s} misfit={rec['misfit']:.4f}  ({rec['secs']:.0f}s)  "
          f"{ {n: round(v,2) for n,v in rec['model'].items()} }", flush=True)
    return rec, res, ranges


def build_matrix():
    configs = []
    # 1) repetibilidad
    for sd in [11, 22, 33, 44, 55]:
        configs.append((f"rep_seed{sd}", REF_BAND, REFSEL, sd))
    # 2) jackknife (deja fuera cada estacion activa, sobre REFSEL)
    for st in STATIONS:
        if REFSEL[st] == (0, 0, 0):
            continue
        sel = dict(REFSEL); sel[st] = (0, 0, 0)
        configs.append((f"jack_no_{st}", REF_BAND, sel, REF_SEED))
    # 3) bandas (0.06-0.15 ya esta en rep_seed11)
    for band in [(0.02, 0.10), (0.04, 0.12)]:
        configs.append((f"band_{band[0]:.2f}-{band[1]:.2f}", band, REFSEL, REF_SEED))
    return configs


def aggregate(records, ranges):
    """Tabla + figura de dispersion por parametro + informe."""
    import numpy as np
    blocks = {"rep": [], "jack": [], "band": []}
    for r in records:
        k = "rep" if r["label"].startswith("rep") else "jack" if r["label"].startswith("jack") else "band"
        blocks[k].append(r)
    allv = {p: np.array([r["model"][p] for r in records]) for p in PNAMES}

    # --- figura: boxplot normalizado por rango ---
    fig, ax = plt.subplots(figsize=(11, 5))
    norm = []
    for i, p in enumerate(PNAMES):
        lo, hi = (ranges[i] if i < len(ranges) else (allv[p].min(), allv[p].max()))
        norm.append((allv[p] - lo) / (hi - lo) if hi > lo else allv[p] * 0)
    ax.boxplot(norm, labels=PNAMES, showmeans=True)
    colors = {"rep": "tab:blue", "jack": "tab:orange", "band": "tab:green"}
    for i, p in enumerate(PNAMES):
        lo, hi = ranges[i]
        for k, rs in blocks.items():
            for r in rs:
                y = (r["model"][p] - lo) / (hi - lo)
                ax.plot(i + 1 + np.random.uniform(-0.12, 0.12), y, "o", ms=4,
                        color=colors[k], alpha=0.7)
    ax.set_ylim(-0.05, 1.05); ax.set_ylabel("valor normalizado al rango [0,1]")
    ax.set_title("Robustez 2020: dispersion de parametros (azul=semillas, naranja=jackknife, verde=bandas)")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(OUT / "robustez_dispersion.png", dpi=150)

    # --- informe markdown ---
    L = ["# Informe de robustez — inversion cinematica 2020\n",
         f"Corridas: {len(records)} (repetibilidad {len(blocks['rep'])}, "
         f"jackknife {len(blocks['jack'])}, bandas {len(blocks['band'])}).\n",
         "## Dispersion por parametro (todas las corridas)\n",
         "| Param | mediana | min | max | IQR | rango | ¿robusto? |",
         "|---|---|---|---|---|---|---|"]
    for i, p in enumerate(PNAMES):
        v = allv[p]; lo, hi = ranges[i]
        q1, q3 = np.percentile(v, [25, 75]); iqr = q3 - q1
        frac = iqr / (hi - lo) if hi > lo else 0
        rob = "si" if frac < 0.10 else ("medio" if frac < 0.25 else "NO")
        L.append(f"| {p} | {np.median(v):.3f} | {v.min():.3f} | {v.max():.3f} | "
                 f"{iqr:.3f} | [{lo:g},{hi:g}] | {rob} ({100*frac:.0f}%) |")
    L += ["\n## Misfit por corrida\n", "| corrida | banda | misfit |", "|---|---|---|"]
    for r in sorted(records, key=lambda r: r["misfit"]):
        L.append(f"| {r['label']} | {r['band'][0]}-{r['band'][1]} | {r['misfit']:.4f} |")
    L += ["\n**Lectura:** un parametro es robusto si su IQR es <10% del rango y "
          "coincide entre semillas, jackknife y bandas. Parametros con IQR alto o "
          "que dependen de quitar una estacion (jackknife) NO estan bien restringidos.\n"]
    (OUT / "INFORME_robustez.md").write_text("\n".join(L))
    print("\n".join(L[:30]))
    print(f"\n-> {OUT}/INFORME_robustez.md  +  robustez_dispersion.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    configs = build_matrix()
    if args.smoke:
        configs = configs[:1]   # una sola corrida minima
    print(f"=== ROBUSTEZ 2020 — {len(configs)} corridas (smoke={args.smoke}) ===", flush=True)
    records, ranges = [], []
    for label, band, sel, seed in configs:
        rec, res, rg = run_one(label, band, sel, seed, smoke=args.smoke)
        records.append(rec)
        if rg: ranges = rg

    # Appraisal del caso de referencia (rep_seed11) — reusa su ensamble
    if not args.smoke:
        try:
            ref = joblib.load(RESDIR / f"res_rep_seed{REF_SEED}.joblib")
            print("\n=== Appraisal bayesiano (referencia rep_seed%d) ===" % REF_SEED, flush=True)
            ref.run_appraisal(n_resample=20000, verbose=False, save=False)
            try:
                ref.plot_appraisal(save_path=str(OUT / "appraisal_corner.png"), show=False)
            except Exception as e:
                print("plot_appraisal fallo:", e)
        except Exception as e:
            print("appraisal fallo:", e)

    if ranges:
        aggregate(records, ranges)
    print("\n=== ROBUSTEZ COMPLETA ===", flush=True)


if __name__ == "__main__":
    main()
