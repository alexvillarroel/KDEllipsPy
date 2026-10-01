"""Robustez: fracción de momento que cruza el plano de Pichilemu en los N mejores
modelos cinemáticos de las 5 semillas (y en todos los modelos, por comparación)."""
import json
import sys
from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
CASE = NAV / "isc_tomo" / "b0.05-0.30_noM14L"
sys.path.insert(0, str(HERE))
from intersection import K, frame  # noqa: E402

cfg = kde.ConfigParser(str(CASE / "input.ctl"))
fm = kde.AxitraForwardModel.from_config(cfg)
base = fm.build_geometry()
sp = cfg.source_position
en = lambda lo, la: np.array([(lo - sp.longitude) * K * np.cos(np.radians(sp.latitude)), (la - sp.latitude) * K, 0.0])
gj = json.loads((HERE / "chaf_navidad.geojson").read_text())
pich = next(f for f in gj["features"] if f["properties"].get("F_name") == "Pichilemu" and f["properties"].get("FT_name") == "Principal")
q = np.mean([en(lo, la) for lo, la in pich["geometry"]["coordinates"]], axis=0)
_, _, n2 = frame(float(pich["properties"]["strike"]), float(pich["properties"]["dip"]))
P = np.array([en(sf.lon, sf.lat) + np.array([0, 0, -sf.z_m / 1e3]) for sf in base.subfaults])
hyp = np.sign((np.array([0, 0, -sp.depth]) - q) @ n2)
cross = np.sign((P - q) @ n2) != hyp
mu = np.array([sf.mu_pa for sf in base.subfaults])

models = [m for p in sorted((CASE / "output_multiseed").glob("seed*.joblib")) for m in joblib.load(p).all_models]
models.sort(key=lambda m: m.misfit)
from copy import deepcopy


def frac(m):
    g = fm.apply_ellipse_model_to_geometry(deepcopy(base), np.asarray(m.model[:7]), keep_all_sources=True)
    s = mu * np.array([sf.slip_m for sf in g.subfaults])
    return s[cross].sum() / s.sum() if s.sum() > 0 else np.nan


res = {}
for n in (10, 100, 500):
    f = np.array([frac(m) for m in models[:n]])
    res[f"top{n}"] = {"misfit_max": models[n - 1].misfit, "median": float(np.nanmedian(f)), "p90": float(np.nanpercentile(f, 90)),
                      "share_lt_10pct": float(np.mean(f < 0.10))}
    print(f"top {n:4d} (misfit ≤ {models[n - 1].misfit:.3f}): cruza mediana {100 * np.nanmedian(f):.1f} %, p90 {100 * np.nanpercentile(f, 90):.1f} %, "
          f"modelos con <10 % cruzando: {100 * np.mean(f < 0.10):.0f} %")
rng = np.random.default_rng(0)
f = np.array([frac(m) for m in rng.choice(models, 500, replace=False)])
res["random500"] = {"median": float(np.nanmedian(f)), "share_lt_10pct": float(np.mean(f < 0.10))}
print(f"500 al azar (todo el ensamble): cruza mediana {100 * np.nanmedian(f):.1f} %, modelos con <10 %: {100 * np.mean(f < 0.10):.0f} %")
(HERE / "pichilemu_ensemble_crossing.json").write_text(json.dumps(res, indent=1))

# Referencia sin datos: modelos uniformes dentro de los rangos de la inversión (dt0 no afecta).
lo_hi = np.array([[p.min_val, p.max_val] for p in cfg.inversion_params.parameters])
U = rng.uniform(lo_hi[:, 0], lo_hi[:, 1], size=(2000, len(lo_hi)))


class M:  # mismo interfaz que NAModel para frac()
    def __init__(self, v):
        self.model = v


f = np.array([frac(M(u)) for u in U])
f = f[np.isfinite(f)]
res["prior_uniform"] = {"n": int(f.size), "median": float(np.median(f)), "share_lt_10pct": float(np.mean(f < 0.10)),
                        "share_gt_30pct": float(np.mean(f > 0.30))}
print(f"{f.size} uniformes en los rangos (sin datos): cruza mediana {100 * np.median(f):.1f} %, <10 %: {100 * np.mean(f < 0.10):.0f} %, >30 %: {100 * np.mean(f > 0.30):.0f} %")
(HERE / "pichilemu_ensemble_crossing.json").write_text(json.dumps(res, indent=1))
