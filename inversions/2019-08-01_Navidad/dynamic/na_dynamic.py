"""Inversión NA dinámica (Navidad): geometría de la aspereza + fricción.

Libres : a, b, u, v, phi, Te, cte1, Dc
         a, b = semiejes en ptos de la grilla 25x25 (~2.08 km)
         u, v = posición del hipocentro dentro de la elipse, en unidades de
                (a-r, b-r): |u|,|v| <= 0.7 garantiza que la ruptura nuclea (el
                parche de nucleación está fijo en el pto 12,12)
         cte1 = resistencia de pico / Te;  Dc = distancia crítica (m)
         dt0  = corrección de hora de origen (s, >0 = más tarde): loc_* usa
                hipocentro movido con la hora de origen CSN -> desfase de ~3 s
Fijos  : cte2=1.1, r=1.5
Caso   : DYN_CASE (np1 por defecto), p.ej.  DYN_CASE=loc_mid python na_dynamic.py
M0     : DYN_M0=<N·m> impone el momento escalar (análogo a mt_strict): Te queda
         fijo en TE_REF y sale de los libres; por autosemejanza del
         slip-weakening el modelo equivale a Te*k, Dc*k con k = M0/M0_simulado.
Serie (n_jobs=1) dentro de cada corrida; DYN_WORK separa corridas simultáneas.
"""

import os

import joblib
import numpy as np

import kdellipspy as kde
from kdellipspy.core.geometry import tsn_hypocentre_coarse
from kdellipspy.inversion.dynamic import DynamicNAInversionModel
from kdellipspy.inversion.kinematic.model_na import NAConfig

from forward_dynamic import _FP, CASE, CASE_DIR, DEFAULT_MODEL, HERE, setup_work_dir, tsn_grid, tsn_run_cfg

M0_TARGET = float(os.environ["DYN_M0"]) if os.environ.get("DYN_M0") else None
TE_REF = 3.0  # MPa, Te de referencia en modo M0 impuesto

ALL_NAMES = ["a (pts)", "b (pts)", "u", "v", "phi (rad)", "Te (MPa)", "cte1", "Dc (m)", "dt0 (s)"]
ALL_RANGES = np.array([
    [2.0, 8.0],    # a    (semieje, ptos ~2.1 km; > r)
    [2.0, 8.0],    # b    (semieje, ptos ~2.1 km; > r)
    [-0.7, 0.7],   # u
    [-0.7, 0.7],   # v
    [0.0, np.pi],  # phi  (rad)
    [1.0, 8.0],    # Te   (MPa)
    [1.05, 1.8],   # cte1 (pico = cte1*Te)
    [0.4, 2.5],    # Dc   (m); >=0.4 -> zona cohesiva >= ~2 celdas de 500 m
    [-3.0, 8.0],   # dt0  (s): cinemático ISC da ~+3.9 s; la nucleación dinámica suma retardo
])
FREE = [i for i, n in enumerate(ALL_NAMES) if not (M0_TARGET and n.startswith("Te"))]
NAMES = [ALL_NAMES[i] for i in FREE]
RANGES = ALL_RANGES[FREE]
# ~15-20 s/eval con la malla 32x40 km -> 1500 evals ~ 7 h
NA = NAConfig(n_samples_initial=300, n_samples_iteration=60, n_iterations=20,
              n_cells_resample=12, n_jobs=1, random_seed=0)
# Referencia para el chequeo previo: mejor modelo dinámico de loc_mid, en
# parámetros relativos al hipocentro (válido para cualquier malla); se prueba
# también con v invertido porque la convención del dip cambió (fd3d -> axitra).
REFERENCE_FREE = {"a": 6.52, "b": 4.59, "u": 0.095, "v": 0.556, "phi": 1.442, "Te": 3.0, "cte1": 1.34, "Dc": 1.12}
PREFLIGHT_MAX_MISFIT = float(os.environ.get("DYN_PREFLIGHT_MAX", "0.95"))
HYPO_X, HYPO_Y = tsn_hypocentre_coarse(_FP)  # nucleación en el hipocentro del input.ctl


def _as_dict(m_free):
    d = dict(zip((ALL_NAMES[i].split()[0] for i in FREE), m_free))
    d.setdefault("Te", TE_REF)
    return d


def to_full_model(m_free):
    """Libres -> vector de 10 del solver, con el centro (xo, yo) tal que el
    hipocentro cae en (u*(a-r), v*(b-r)) del marco de la elipse."""
    d = _as_dict(m_free)
    a, b, u, v, phi, te, cte1, dc = (d[k] for k in ("a", "b", "u", "v", "phi", "Te", "cte1", "Dc"))
    r = DEFAULT_MODEL[8]
    xh, yh = u * (a - r), v * (b - r)
    dx = xh * np.cos(phi) - yh * np.sin(phi)
    dy = xh * np.sin(phi) + yh * np.cos(phi)
    full = np.array(DEFAULT_MODEL, dtype=np.float32)
    full[[0, 1, 2, 3, 4, 5, 6, 9]] = [a, b, HYPO_X - dx, HYPO_Y - dy, phi, te, cte1, dc]
    return full


class FreeSubsetNA(DynamicNAInversionModel):
    """NA ve solo los parámetros libres (neighpy normaliza por max-min, así
    que un rango de ancho cero para 'fijar' un parámetro no sirve). Se
    convierte en _evaluate_model para que log/checkpoint muestren los libres."""

    def _evaluate_model(self, model):
        self.dynamic_fm.time_shift_s = float(_as_dict(model)["dt0"])
        return super()._evaluate_model(to_full_model(model))


def preflight(inv):
    """Referencia (y su espejo en v) con dt0 = -2..8 s; un fd3d por variante
    gracias al caché del forward. Aborta si nada ajusta: mejor perder minutos
    que una noche."""
    rows = []
    for vsign in (1.0, -1.0):
        ref = dict(REFERENCE_FREE, v=vsign * REFERENCE_FREE["v"])
        if M0_TARGET:
            ref["Te"] = TE_REF
        for dt0 in np.arange(-2.0, 8.5, 1.0):
            free = np.array([dict(ref, dt0=dt0)[n.split()[0]] for n in NAMES])
            rows.append((FreeSubsetNA._evaluate_model(inv, free)[0], dt0, vsign))
            print(f"[preflight] v*{vsign:+.0f} dt0={dt0:+.0f} s  misfit={rows[-1][0]:.4f}", flush=True)
    best = min(rows)
    if best[0] > PREFLIGHT_MAX_MISFIT:
        raise SystemExit(f"[preflight] ABORT: best reference misfit {best[0]:.3f} > {PREFLIGHT_MAX_MISFIT}")
    print(f"[preflight] OK: best dt0={best[1]:+.0f} s (v*{best[2]:+.0f}) misfit={best[0]:.4f}", flush=True)


def main():
    cfg = kde.ConfigParser(str(CASE_DIR / "input.ctl"))
    setup_work_dir(cfg)
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE_DIR / "input.ctl"), data_dir=str(CASE_DIR / "DATA"))
    inv = FreeSubsetNA(
        config=cfg, observed_waveforms=observed, time_array=time_array,
        tsn_run_cfg=tsn_run_cfg(),
        tsn_grid=tsn_grid(cfg),
    )
    inv.param_ranges = RANGES
    inv.param_names = list(NAMES)
    inv.dynamic_fm.m0_target = M0_TARGET
    out = HERE / (f"na_output_{CASE}" + ("_M0" if M0_TARGET else ""))
    out.mkdir(exist_ok=True)
    inv.checkpoint_path = out / "best_model_live.txt"

    # Operador de axitra por subfalla: una vez (~7 min), antes del NA.
    inv.dynamic_fm.ensure_basis(int(cfg.observed_data.units))
    print("[setup] axitra basis ready", flush=True)
    preflight(inv)

    result = inv.run_na_search(NA)
    joblib.dump(result, out / "na_result.joblib")
    np.save(out / "best_synthetics.npy", inv.best_synthetics)
    np.save(out / "best_full_model.npy", to_full_model(result.best_model.model))
    print("best misfit", result.best_model.misfit, flush=True)
    if M0_TARGET:
        inv._evaluate_model(result.best_model.model)  # recomputa k del mejor modelo
        k = inv.dynamic_fm.last_m0_scale
        full = to_full_model(result.best_model.model)
        print(f"M0 impuesto {M0_TARGET:.3e} N·m: k={k:.3f} -> Te real {TE_REF*k:.2f} MPa, "
              f"pico {TE_REF*k*full[6]:.2f} MPa, Dc real {full[9]*k:.2f} m", flush=True)
    inv.clean()


def _check():
    """El hipocentro debe caer dentro de la elipse reducida (nuclea)."""
    from kdellipspy.core.geometry import EllipticalStressMapper
    rng = np.random.default_rng(1)
    for _ in range(500):
        m = to_full_model(rng.uniform(RANGES[:, 0], RANGES[:, 1]))  # dt0 no afecta la nucleación
        pre, peak = EllipticalStressMapper(_FP.nx, _FP.ny, hypo=(HYPO_X, HYPO_Y)).fields(m)
        i, j = int(round(HYPO_X)) - 1, int(round(HYPO_Y)) - 1
        assert pre[i, j] > peak[i, j], m  # nucleation patch above strength


if __name__ == "__main__":
    main()
