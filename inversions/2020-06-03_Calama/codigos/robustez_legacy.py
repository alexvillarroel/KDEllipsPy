#!/usr/bin/env python3
"""
robustez_legacy.py
==================
Reproduce el escenario del LEGACY FORTRAN (Kinematic_inversion) en kdellipspy y le
aplica el MISMO estudio de robustez que robustez_2020.py, para responder: con
velocidad + banda 0.02-0.1 + las 8 estaciones del legacy (donde el misfit daba
~0.2), ¿la GEOMETRIA queda restringida? (misfit bajo != geometria unica).

Escenario legacy:
  * datos = VELOCIDAD (los real_disp_* del legacy son velocidad, mal nombrados;
    verificado por correlacion con SAC/VEL). Se usan tal cual (ya 0.02-0.1, 512pts).
  * 8 estaciones: A05C AC01 AC02 PB05 PB06 PB08 PB19 T15A
  * plano 32x32 km, 30x30 subfallas, banda 0.02-0.1, modelo de velocidades identico.

Matriz: 5 semillas + 8 jackknife + 2 sub-bandas (0.03-0.08, 0.04-0.10; dentro de
0.02-0.1) + appraisal. Salidas en legacy_equiv/robustez/.

Uso:
    python codigos/robustez_legacy.py            # matriz completa
    python codigos/robustez_legacy.py --smoke    # prueba rapida
"""
import sys, json, time, argparse
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import joblib

BASE = Path(__file__).resolve().parent.parent
LEGDATA = BASE / "legacy" / "Kinematic_inversion" / "DATA"
EQUIV = BASE / "legacy_equiv"; EQUIV.mkdir(exist_ok=True)
OUT = EQUIV / "robustez"; OUT.mkdir(exist_ok=True)
RESDIR = OUT / "results"; RESDIR.mkdir(exist_ok=True)
sys.path.insert(0, str(BASE / "codigos"))
import kdellipspy as kde
try:
    from kdellipspy import NAConfig
except Exception:
    from kdellipspy.inversion.kinematic.model_na import NAConfig

PNAMES = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]
# 8 estaciones del legacy (orden de los real_disp_*): (nombre, lat, lon)
STA = [("A05C", -27.3612, -70.3390), ("AC01", -26.1479, -70.5987),
       ("AC02", -26.8355, -69.1291), ("PB05", -22.8528, -70.2023),
       ("PB06", -22.7058, -69.5717), ("PB08", -20.1411, -69.1535),
       ("PB19", -23.9048, -69.2900), ("T15A", -20.2393, -70.0540)]
NP, DT = 512, 0.25
REF_BAND = (0.02, 0.10)
REF_SEED = 11

CONFIG_TXT = """\
################################################################################
#             Control file for KINEMATIC INVERSION (Ellipse Patch)  LEGACY-EQUIV
################################################################################
#===============================================:===============================
# 1. Observed Data Parameters                   :     Values                   |
#===============================================:===============================
 Time window start (t1)                         :     0.0
 Time window end (t2)                           :   128.0
 Number of points (Npts)                        :     512
 Delta / Time step                              :    0.25
 Units                        (1:disp, 2:vel)   :       2

#===============================================:===============================
# 2. Source & Focal Mechanism                   :     Values                   |
#===============================================:===============================
 Event Name                                     :     Calama 2020 Intraplate
 Origin Time (UTC)                              :     2020-06-03T07:35:28.000Z
 Latitude                                       :   -23.2470
 Longitude                                      :   -68.530
 Depth                                     (km) :      123.400
 Strike                                  (deg.) :      333.0
 Dip                                     (deg.) :      60.0
 Rake                                    (deg.) :    -91.0

#===============================================:===============================
# 3. Fault Plane Parameters                     :     Values                   |
#===============================================:===============================
 Length along strike (Lx)                  (m)  : 32000.0
 Length along dip (Ly)                     (m)  : 32000.0
 Hypocenter position strike (Hx)           (m)  : 16000.0
 Hypocenter position dip (Hy)              (m)  : 16000.0
 Number of subfaults along strike (Nx)          :       30
 Number of subfaults along dip (Ny)             :       30

#===============================================:===============================
# 4. Ellipse Parameters & Frequency Band        :     Values                   |
#===============================================:===============================
 Number of ellipses                             :       1
 Initial slip                 (0:no, 1:yes)     :       0
 Slip shape            (0:cst, 1:gauss, 2:ell)  :       1
 Frequency 1 (Freq1)                      (Hz)  :    0.02
 Frequency 2 (Freq2)                      (Hz)  :    0.1
 Time shift (T0)                           (s)  :     3.0
 Source type (axitra)                           :     4
 Zerophase filter (0:causal, 1:acausal)         :     0

#===============================================:===============================
# 5. Parameters to Invert (Min, Max, Flag)      :   Min      Max     Flag      |
#===============================================:===============================
 Number of parameters                           :       7
 Param 1: Length of axis 1                 (km) :     4.0    13.0      1
 Param 2: Length of axis 2                 (km) :     4.0    13.0      1
 Param 3: Rotation angle                 (x pi) :     0.0    2.0       1
 Param 4: Position of the center np             :     0.0    1.0       1
 Param 5: Position of the center tp     (x 2pi) :     0.0    1.0       1
 Param 6: Maximum slip (Dmax)               (m) :     1.0    6.0       1
 Param 7: Rupture velocity (Vr)          (km/s) :     0.5    3.5       1

#===============================================:===============================
# 6. Inversion Process Parameters               :     Values                   |
#===============================================:===============================
 Algorithm type                 (0:NA, 1:MC)    :      0
 NA: Number of iterations                           :      50
 NA: Sample size for first iteration (SS1)          :      500
 NA: Sample size for other iterations               :      50
 NA: Cells to resample                              :      10
 Misfit time window     (s, 0.0=Full signal)    :      20.0
 MCMC: total steps                               :      500
 MCMC: burn-in                                   :      100
 MCMC: proposal scale                            :      0.05
 MCMC: thinning                                  :      1
 MCMC: chains                                    :      1

#===============================================:===============================
# 7. Optional: Moment Tensor (Full MT)          :     Values                   |
#===============================================:===============================
 Moment Tensor Flag           (0:no, 1:yes)     :       0
 MT Scaling Mode (mt_strict/mt_factored)      : mt_factored
 Mrr                                            :   0.0
 Mtt                                            :   0.0
 Mpp                                            :   0.0
 Mrt                                            :   0.0
 Mrp                                            :   0.0
 Mtp                                            :   0.0
 Exponent (iexp)                                :    18.0

#===============================================:===============================
# 8. Station Parameters (Lat, Lon, Elev, Name, N, E, Z)                         |
#===============================================:===============================
 Number of stations                             :       8
__STATIONS__

#===============================================:===============================
# 9. Velocity Model 1D                              : Thickness Vp Vs Rho Qp Qs |
#===============================================:===============================
 Number of layers                               :       16
 0.000 5210.000 2990.000 2500.000  600.000  400.000
 2500.000 5370.000 3090.000 2500.000  600.000  400.000
 4500.000 5550.000 3190.000 2500.000  600.000  400.000
 6500.000 5720.000 3290.000 2600.000  600.000  400.000
 8500.000 5890.000 3390.000 2600.000  600.000  400.000
 10500.000 5980.000 3440.000 2600.000  600.000  400.000
 15000.000 6800.000 3750.000 2800.000  600.000  400.000
 20000.000 6810.000 3880.000 2800.000  600.000  400.000
 25000.000 6950.000 3940.000 3000.000  600.000  400.000
 30000.000 6980.000 4050.000 3000.000  600.000  400.000
 35000.000 7110.000 4110.000 3100.000  600.000  400.000
 40000.000 7410.000 4180.000 3300.000  600.000  400.000
 45000.000 7690.000 4300.000 3300.000  600.000  400.000
 50000.000 8050.000 4390.000 3300.000  600.000  400.000
 60000.000 8480.000 4730.000 3400.000  600.000  400.000
 70000.000 8480.000 4780.000 3400.000  600.000  400.000
"""


def write_config():
    rows = "\n".join(f"{la:.4f} {lo:.4f} 0.0000 {nm}   1 1 1" for nm, la, lo in STA)
    (EQUIV / "input.ctl").write_text(CONFIG_TXT.replace("__STATIONS__", rows))
    return EQUIV / "input.ctl"


def load_observed(band):
    """Data de VELOCIDAD del legacy (8,3,512) N,E,Z; re-filtra causal si sub-banda."""
    raw = {c: np.loadtxt(LEGDATA / f"real_disp_{c}").reshape(len(STA), NP) for c in "xyz"}
    obs = np.stack([raw["x"], raw["y"], raw["z"]], axis=1)
    if (round(band[0], 3), round(band[1], 3)) != REF_BAND:
        t = np.arange(NP) * DT
        obs = kde.bandpass_filter_waveforms(obs, t, float(band[0]), float(band[1]),
                                            corners=4, zerophase=False)
    return obs


def make_na(smoke, seed):
    if smoke:
        return NAConfig(n_samples_initial=20, n_samples_iteration=10, n_iterations=2,
                        n_cells_resample=2, n_jobs=1, random_seed=seed)
    return NAConfig(n_samples_initial=500, n_samples_iteration=50, n_iterations=50,
                    n_cells_resample=10, n_jobs=1, random_seed=seed)


def run_one(label, band, drop, seed, smoke=False):
    t0 = time.time()
    cfg = kde.ConfigParser(filepath=str(EQUIV / "input.ctl"))
    cfg.ellipse.freq1, cfg.ellipse.freq2 = float(band[0]), float(band[1])
    for s in cfg.stations.stations:
        on = (s.name != drop)
        s.use_n = s.use_e = s.use_z = on
    obs = load_observed(band)
    t = np.arange(NP) * DT
    azt = np.asarray(kde.build_azi_times_array(config=cfg, model_name="iasp91"), float)
    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs,
                                 time_array=t, azi_times_array=azt)
    model.use_green_cache = True
    res = model.run_na_search(make_na(smoke, seed))
    try:
        model.clear_green_cache()
    except Exception:
        pass
    b = np.asarray(res.best_model.model, float)
    rec = dict(label=label, band=list(band), drop=drop, seed=int(seed),
               misfit=float(res.best_model.misfit),
               model={n: float(x) for n, x in zip(PNAMES, b)},
               secs=round(time.time() - t0, 1))
    joblib.dump(res, RESDIR / f"res_{label}.joblib")
    (RESDIR / f"rec_{label}.json").write_text(json.dumps(rec, indent=2))
    print(f"[OK] {label:18s} misfit={rec['misfit']:.4f} ({rec['secs']:.0f}s) "
          f"{ {n: round(v,2) for n,v in rec['model'].items()} }", flush=True)
    return rec, list(getattr(model, "param_ranges", []))


def build_matrix():
    cfgs = []
    for sd in [11, 22, 33, 44, 55]:
        cfgs.append((f"rep_seed{sd}", REF_BAND, None, sd))
    for nm, _, _ in STA:
        cfgs.append((f"jack_no_{nm}", REF_BAND, nm, REF_SEED))
    for band in [(0.03, 0.08), (0.04, 0.10)]:
        cfgs.append((f"band_{band[0]:.2f}-{band[1]:.2f}", band, None, REF_SEED))
    return cfgs


def aggregate(records, ranges):
    blocks = {"rep": [], "jack": [], "band": []}
    for r in records:
        k = "rep" if r["label"].startswith("rep") else "jack" if r["label"].startswith("jack") else "band"
        blocks[k].append(r)
    allv = {p: np.array([r["model"][p] for r in records]) for p in PNAMES}
    fig, ax = plt.subplots(figsize=(11, 5))
    norm = [(allv[p] - ranges[i][0]) / (ranges[i][1] - ranges[i][0]) for i, p in enumerate(PNAMES)]
    ax.boxplot(norm, labels=PNAMES, showmeans=True)
    colors = {"rep": "tab:blue", "jack": "tab:orange", "band": "tab:green"}
    for i, p in enumerate(PNAMES):
        lo, hi = ranges[i]
        for k, rs in blocks.items():
            for r in rs:
                ax.plot(i + 1 + np.random.uniform(-0.12, 0.12), (r["model"][p] - lo) / (hi - lo),
                        "o", ms=4, color=colors[k], alpha=0.7)
    ax.set_ylim(-0.05, 1.05); ax.set_ylabel("valor normalizado al rango [0,1]")
    ax.set_title("Robustez LEGACY-EQUIV (vel, 0.02-0.1, 8 est): azul=semillas, naranja=jackknife, verde=bandas")
    ax.grid(alpha=0.3, axis="y"); fig.tight_layout()
    fig.savefig(OUT / "robustez_dispersion.png", dpi=150)
    L = ["# Robustez — escenario LEGACY-EQUIV (velocidad, 0.02-0.1 Hz, 8 estaciones)\n",
         f"Corridas: {len(records)} (rep {len(blocks['rep'])}, jack {len(blocks['jack'])}, band {len(blocks['band'])}).\n",
         "| Param | mediana | min | max | IQR | rango | ¿robusto? |", "|---|---|---|---|---|---|---|"]
    for i, p in enumerate(PNAMES):
        v = allv[p]; lo, hi = ranges[i]; q1, q3 = np.percentile(v, [25, 75]); iqr = q3 - q1
        frac = iqr / (hi - lo); rob = "si" if frac < 0.10 else ("medio" if frac < 0.25 else "NO")
        L.append(f"| {p} | {np.median(v):.3f} | {v.min():.3f} | {v.max():.3f} | {iqr:.3f} | [{lo:g},{hi:g}] | {rob} ({100*frac:.0f}%) |")
    L += ["\n| corrida | banda | misfit |", "|---|---|---|"]
    for r in sorted(records, key=lambda r: r["misfit"]):
        L.append(f"| {r['label']} | {r['band'][0]}-{r['band'][1]} | {r['misfit']:.4f} |")
    (OUT / "INFORME_robustez_legacy.md").write_text("\n".join(L))
    print("\n".join(L[:14]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    write_config()
    cfgs = build_matrix()
    if args.smoke:
        cfgs = cfgs[:1]
    print(f"=== ROBUSTEZ LEGACY-EQUIV — {len(cfgs)} corridas (smoke={args.smoke}) ===", flush=True)
    records, ranges = [], []
    for label, band, drop, seed in cfgs:
        rec, rg = run_one(label, band, drop, seed, smoke=args.smoke)
        records.append(rec)
        if rg: ranges = rg
    if not args.smoke:
        try:
            ref = joblib.load(RESDIR / f"res_rep_seed{REF_SEED}.joblib")
            print("=== Appraisal (rep_seed%d) ===" % REF_SEED, flush=True)
            ref.run_appraisal(n_resample=20000, verbose=False, save=False)
            ref.plot_appraisal(save_path=str(OUT / "appraisal_corner.png"), show=False)
        except Exception as e:
            print("appraisal fallo:", e)
    if ranges:
        aggregate(records, ranges)
    print("=== LEGACY-EQUIV COMPLETO ===", flush=True)


if __name__ == "__main__":
    main()
