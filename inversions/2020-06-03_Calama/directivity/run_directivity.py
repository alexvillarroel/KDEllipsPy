"""Test de directividad: inversión cinemática con el centro de la elipse obligado
a quedar en un sector respecto del hipocentro, en coordenadas del plano de falla.

Caso de Calama 2020: np1 (335/63/-85), hipocentro ISOLA, modelo de Potin,
0.04-0.15 Hz, sin A09F ni la transversal de A08F.

NO se espera que esto RESUELVA la dirección, sino que la ACOTE: el gap azimutal
es de 146 grados (vacío de 82 a 228) y, con la fuente a 115.6 km y estaciones a
49-156 km, todos los rayos salen en un cono de take-off de 23 a 53 grados.
Además el cinemático ya mostró que Vr no está resuelto (valle Vr-dt0 sin mínimo,
ver kin_tests/VR_NO_RESUELTO.md), y directividad y duración son la misma
información. Lo esperable es que varios sectores salgan indistinguibles; eso
también es un resultado.

El parámetro tp (ángulo en el marco de la elipse) se reemplaza por psi: dirección
del desplazamiento centro - hipocentro en el plano de falla (0° = +rumbo, NNE;
90° = +manteo, hacia abajo). Para cada (a1, a2, alpha) se calcula el tp que da
esa dirección, así NA explora solo el sector pedido sin rechazar modelos.

Uso:  SECTOR=nne|downdip|ssw|updip|centred python run_directivity.py
Escribe <SECTOR>/seed<k>.joblib y <SECTOR>/summary.txt
"""

import os
import sys
from math import atan2, cos, pi, radians, sin
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
sys.path[:0] = [str(NAV / "kin_tests")]
os.environ.setdefault("KIN_CASE", "kin_tests/np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F")
import kdellipspy as kde  # noqa: E402
from kdellipspy.inversion.kinematic.model_na import NAConfig  # noqa: E402
from run_kin_multiseed import DT0_RANGE, NA, KinematicDt0NA  # noqa: E402

CASE = NAV / os.environ["KIN_CASE"]
# (psi_min, psi_max) en grados; np mínimo para que el desplazamiento sea real
SECTORS = {"nne": (-45, 45), "downdip": (45, 135), "ssw": (135, 225), "updip": (225, 315), "centred": None}
SECTOR = os.environ.get("SECTOR", "nne")
NP_DIRECTIVE = (0.4, 1.0)
NP_CENTRED = (0.0, 0.2)
SEEDS = [int(s) for s in os.environ.get("NA_SEEDS", "0,1,2,3,4,5,6,7").split(",")]


def tp_for_direction(psi_deg, a1, a2, alpha_frac):
    """tp (x 2pi) tal que el centro de la elipse quede en la dirección psi desde el
    hipocentro (inversa de EllipticalSlipMapper._prepare)."""
    alpha = alpha_frac * pi
    dx, dy = cos(radians(psi_deg)), sin(radians(psi_deg))
    x01 = dx * cos(alpha) - dy * sin(alpha)
    y01 = dx * sin(alpha) + dy * cos(alpha)
    return (atan2(y01 / a2, x01 / a1) / (2 * pi)) % 1.0


class DirectivityNA(KinematicDt0NA):
    def _build_geometry_from_parameters(self, model, keep_all_sources=False):
        m = np.array(model, dtype=float)
        if SECTORS[SECTOR] is not None:
            m[4] = tp_for_direction(m[4], m[0], m[1], m[2])  # m[4] llega como psi (grados)
        return super()._build_geometry_from_parameters(m, keep_all_sources=keep_all_sources)


def check(inv):
    """psi = 0/180 debe mover el centroide de slip al N/S y psi = 90/270 a mayor/menor
    profundidad (verifica la convención de signos de la inversa de tp)."""
    def centroid(psi):
        m = np.array([4.0, 8.0, 0.3, 0.8, psi, 2.0, 2.0, 0.0])
        sf = inv._build_geometry_from_parameters(m, keep_all_sources=True).subfaults
        w = np.array([f.slip_m for f in sf])
        return (w * [f.lat for f in sf]).sum() / w.sum(), (w * np.abs([f.z_m for f in sf])).sum() / w.sum()
    (n, _), (s, _), (_, dd), (_, du) = (centroid(p) for p in (0.0, 180.0, 90.0, 270.0))
    assert n > s and dd > du, (n, s, dd, du)
    print(f"[check] lat NNE {n:.3f} > SSW {s:.3f};  prof. downdip {dd/1e3:.1f} > updip {du/1e3:.1f} km", flush=True)


def main():
    cfg = kde.ConfigParser(filepath=str(CASE / "input.ctl"))
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(CASE / "input.ctl"), data_dir=str(CASE / "DATA"))
    inv = DirectivityNA(config=cfg, observed_waveforms=observed, time_array=time_array,
                        azi_times_array=kde.build_azi_times_array(config=cfg, model_name="iasp91"))
    inv.use_green_cache = True
    pr = np.vstack([inv.param_ranges, DT0_RANGE]).astype(float)
    names = list(inv.param_names) + ["dt0 (s)"]
    if SECTORS[SECTOR] is None:
        pr[3] = NP_CENTRED
    else:
        pr[3] = NP_DIRECTIVE
        pr[4] = SECTORS[SECTOR]
        names[4] = "psi (deg, 0=NNE, 90=downdip)"
    inv.param_ranges, inv.param_names = pr, names
    if SECTOR == "nne":
        check(inv)

    out = HERE / SECTOR
    out.mkdir(exist_ok=True)
    inv._evaluate_model(np.mean(inv.param_ranges, axis=1))
    inv._keep_cache = True
    rows = []
    for seed in SEEDS:
        inv.checkpoint_path = out / f"seed{seed}_best_live.txt"
        r = inv.run_na_search(NAConfig(random_seed=seed, **NA))
        r.save(str(out / f"seed{seed}.joblib"))
        rows.append((r.best_model.misfit, np.asarray(r.best_model.model)))
        print(f"[{SECTOR} seed {seed}] best misfit {r.best_model.misfit:.4f}", flush=True)
    ms = np.array([x[0] for x in rows])
    models = np.array([x[1] for x in rows])
    lines = [f"{SECTOR}: {len(SEEDS)} semillas",
             f"misfit best {ms.min():.4f} mean {ms.mean():.4f} std {ms.std():.4f} ({', '.join(f'{v:.4f}' for v in ms)})"]
    best = models[np.argmin(ms)]
    lines += [f"{n:<32s} {best[i]:9.3f} {models[:, i].mean():9.3f} {models[:, i].std():8.3f}" for i, n in enumerate(names)]
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
