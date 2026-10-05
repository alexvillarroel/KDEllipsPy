"""NA-Bayes (Sambridge 1999b) para el cinemático robusto: ISC + 1D regional,
0.05-0.30 Hz, sin M14L, 5 semillas x 3400 modelos ya evaluados (sin forwards nuevos).

Verosimilitud: misfit = sum r^2 / sum d^2 (señal completa, 128 s a 1 Hz, 21 trazas).
Gaussiana con el mejor modelo al nivel del ruido:
    log L = -N_eff * misfit / (2 * misfit_min)
    N_eff = 2 * B * T * n_trazas = 2 * 0.25 Hz * 128 s * 21 = 1344
Se corre también N_eff/4 (intervalos ~2x más anchos) para mostrar la sensibilidad.
El default de run_appraisal (T = 2*misfit_min) equivale a N_eff = 1: casi el prior.

Además de los 8 parámetros, cada muestra se transforma a magnitudes físicas
(Mw, centroide, desplazamiento centroide-hipocentro, dimensiones), que no tienen
las simetrías de la elipse (a1<->a2 con theta+1/2).

Escribe posterior_N<n>.npz, summary_N<n>.json, corner_N<n>.png, derived_N<n>.png.
"""

import json
import os
import sys
import time
from pathlib import Path

import joblib
import numpy as np

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
CASE = NAV / "dyn_cases" / "np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F"
os.environ.setdefault("DYN_CASE", "dyn_cases/np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path[:0] = [str(NAV / "okada"), str(NAV / "dynamic"), str(NAV / "kin_tests")]

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from neighpy import NAAppraiser  # noqa: E402

import okada_models as om  # noqa: E402
from kdellipspy.core.plotting import plot_uncertainty_corner  # noqa: E402
from run_kin_dt0 import DT0_RANGE  # noqa: E402

N_FULL = 2 * 0.25 * 128 * 21
N_RESAMPLE = int(os.environ.get("N_RESAMPLE", 20000))
N_CHAINS = 4  # ponytail: una cadena por proceso (CHAIN=k); el multiprocessing de neighpy se colgó con 4 walkers
N_DERIVED = 4000
K = 111.19


def prior_bounds(cfg):
    lines = [l for l in (CASE / "input.ctl").read_text().splitlines() if l.strip().startswith("Param ")]
    b = [[float(v) for v in l.split(":")[-1].split()[:2]] for l in lines]
    return np.array(b + [list(DT0_RANGE)])


def derived(cfg, fm, models):
    """Mw, centroide y dimensiones (momento-ponderados) para cada modelo."""
    sp = cfg.source_position
    s = np.radians(sp.strike)
    kx = K * np.cos(np.radians(sp.latitude))
    out = []
    for m in models:
        g = fm.apply_ellipse_model_to_geometry(fm.build_geometry(), np.asarray(m[:7]), keep_all_sources=True)
        sf = g.subfaults
        m0 = np.array([f.mu_pa * f.area_m2 * f.slip_m for f in sf])
        if m0.sum() <= 0:
            out.append([np.nan] * 9); continue
        w = m0 / m0.sum()
        e = (np.array([f.lon for f in sf]) - sp.longitude) * kx
        n = (np.array([f.lat for f in sf]) - sp.latitude) * K
        z = np.abs(np.array([f.z_m for f in sf])) / 1e3
        along = e * np.sin(s) + n * np.cos(s)                      # + NNE
        down = np.hypot(e * np.cos(s) - n * np.sin(s), z - sp.depth) * np.sign(z - sp.depth)  # + más profundo
        ca, cd = (w * along).sum(), (w * down).sum()
        out.append([(np.log10(m0.sum()) - 9.1) / 1.5,
                    sp.latitude + (w * n).sum() / K, sp.longitude + (w * e).sum() / kx, (w * z).sum(),
                    ca, cd,
                    4 * np.sqrt((w * (along - ca) ** 2).sum()), 4 * np.sqrt((w * (down - cd) ** 2).sum()),
                    float(np.degrees(np.arctan2(cd, ca)))])
    return np.array(out)


DERIVED_NAMES = ["Mw", "centroid lat", "centroid lon", "centroid depth (km)",
                 "centroid along strike (km, +NNE)", "centroid along dip (km, +downdip)",
                 "length along strike (km, 4σ)", "width along dip (km, 4σ)", "centroid direction (°, 0=NNE, 90=downdip)"]


def stats(x):
    x = x[np.isfinite(x)]
    p = np.percentile(x, [2.5, 16, 50, 84, 97.5])
    return dict(mean=float(x.mean()), std=float(x.std()), p2_5=float(p[0]), p16=float(p[1]), median=float(p[2]),
                p84=float(p[3]), p97_5=float(p[4]))


def main():
    cfg = om.kde.ConfigParser(str(CASE / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, _ = om.subfaults(cfg)
    rs = [joblib.load(p) for p in sorted((NAV / "kin_tests" / "np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F" / "multiseed").glob("seed*.joblib"))]
    names = list(rs[0].param_names)
    E = np.vstack([[m.model for m in r.all_models] for r in rs])
    M = np.concatenate([[m.misfit for m in r.all_models] for r in rs])
    ok = np.isfinite(M); E, M = E[ok], M[ok]
    best = E[np.argmin(M)]
    bounds = prior_bounds(cfg)
    print(f"ensemble {E.shape}, misfit min {M.min():.4f}", flush=True)

    chain = os.environ.get("CHAIN")
    for n_eff in (N_FULL, N_FULL / 4):
        tag = f"N{int(round(n_eff))}"
        if chain is not None:  # modo cadena: solo remuestrea y guarda
            t = time.time()
            a = NAAppraiser(n_resample=N_RESAMPLE, n_walkers=1, initial_ensemble=E,
                            log_ppd=-n_eff * M / (2 * M.min()), bounds=tuple(map(tuple, bounds)), verbose=False,
                            seed=100 + int(chain))
            a.run(save=True)
            np.save(HERE / f"chain_{tag}_{chain}.npy", np.asarray(a.samples))
            print(f"[{tag} chain {chain}] {N_RESAMPLE} pasos en {time.time() - t:.0f} s", flush=True)
            continue
        S = np.vstack([np.load(HERE / f"chain_{tag}_{k}.npy") for k in range(N_CHAINS)])
        print(f"[{tag}] {len(S)} muestras de {N_CHAINS} cadenas", flush=True)
        rng = np.random.default_rng(0)
        D = derived(cfg, fm, S[rng.choice(len(S), size=min(N_DERIVED, len(S)), replace=False)])
        D_best = derived(cfg, fm, [best])[0]
        np.savez(HERE / f"posterior_{tag}.npz", samples=S, derived=D, names=names, derived_names=DERIVED_NAMES)
        summ = dict(n_eff=n_eff, n_samples=len(S), n_ensemble=len(E), misfit_min=float(M.min()),
                    params={nm: stats(S[:, i]) | {"best": float(best[i])} for i, nm in enumerate(names)},
                    derived={nm: stats(D[:, i]) | {"best": float(D_best[i])} for i, nm in enumerate(DERIVED_NAMES)},
                    p_downdip=float(np.mean((D[:, 8] > 45) & (D[:, 8] < 135))),
                    p_centroid_NNE=float(np.mean(D[:, 4] > 0)))
        (HERE / f"summary_{tag}.json").write_text(json.dumps(summ, indent=1))
        plot_uncertainty_corner(S, names, bounds=bounds, truths=best, mean=S.mean(0), show=False,
                                save_path=HERE / f"corner_{tag}.png", dpi=110,
                                title=f"NA-Bayes posterior, kinematic 0.05–0.30 Hz (N_eff={n_eff:.0f}, {len(S)} samples)")
        plt.close("all")
        sel = [0, 3, 4, 5, 6, 7, 8]
        fig, axes = plt.subplots(2, 4, figsize=(14, 6.4), constrained_layout=True)
        for ax, i in zip(axes.ravel(), sel):
            x = D[:, i][np.isfinite(D[:, i])]
            ax.hist(x, bins=40, color="0.55", edgecolor="white")
            st = summ["derived"][DERIVED_NAMES[i]]
            for q, ls in (("p2_5", ":"), ("p16", "--"), ("median", "-"), ("p84", "--"), ("p97_5", ":")):
                ax.axvline(st[q], color="k", ls=ls, lw=0.9)
            ax.axvline(D_best[i], color="#c4161c", lw=1.4)
            ax.set_title(f"{DERIVED_NAMES[i]}\n{st['median']:.2f} [{st['p2_5']:.2f}, {st['p97_5']:.2f}]", fontsize=8.5)
            ax.set_yticks([])
        ax = axes.ravel()[-1]
        for i, nm in ((6, "vr (km/s)"), (7, "dt0 (s)")):
            pass
        ax.hist(S[:, 6], bins=40, color="#1f4fbf", alpha=0.7, label="Vr (km/s)")
        st = summ["params"]["vr (km/s)"]
        ax.set_title(f"Vr (km/s)\n{st['median']:.2f} [{st['p2_5']:.2f}, {st['p97_5']:.2f}]", fontsize=8.5)
        ax.axvline(best[6], color="#c4161c", lw=1.4); ax.set_yticks([])
        fig.suptitle(f"NA-Bayes derived quantities (N_eff={n_eff:.0f}): median [95 % interval]; "
                     f"red = best model; dashed = 68 %, dotted = 95 %.  P(centroid downdip) = {summ['p_downdip']:.2f}, "
                     f"P(centroid NNE of hypocentre) = {summ['p_centroid_NNE']:.2f}", fontsize=9.5)
        fig.savefig(HERE / f"derived_{tag}.png", dpi=120)
        plt.close(fig)
        print(json.dumps({k: summ[k] for k in ("p_downdip", "p_centroid_NNE")}), flush=True)
    print("ok", flush=True)


if __name__ == "__main__":
    main()
