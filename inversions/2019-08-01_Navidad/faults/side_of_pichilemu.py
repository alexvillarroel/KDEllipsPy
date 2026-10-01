"""Fracción del momento de cada modelo a cada lado del plano de Pichilemu
(cinemático: mejor de 5 semillas; dinámicos A y B). Usa el slip por subfalla de
okada/okada_models.py y el plano de faults/intersection.py.
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
os.environ.setdefault("DYN_CASE", "isc_tomo/b0.05-0.30_noM14L")
os.environ.setdefault("DYN_WORK", "tsn_work_okada")
sys.path.insert(0, str(NAV / "okada"))
sys.path.insert(0, str(NAV / "dynamic"))
import numpy as np  # noqa: E402

import okada_models as om  # noqa: E402
from intersection import K, frame  # noqa: E402


def main():
    cfg = om.kde.ConfigParser(str(om.fd.CASE_DIR / "input.ctl"))
    om.fd.setup_work_dir(cfg)
    fm, geom = om.subfaults(cfg)
    sp = cfg.source_position
    lat0, lon0 = sp.latitude, sp.longitude
    gj = json.loads((HERE / "chaf_navidad.geojson").read_text())
    pich = next(f for f in gj["features"] if f["properties"].get("F_name") == "Pichilemu" and f["properties"].get("FT_name") == "Principal")
    en = lambda lo, la: np.array([(lo - lon0) * K * np.cos(np.radians(lat0)), (la - lat0) * K, 0.0])
    q = np.mean([en(lo, la) for lo, la in pich["geometry"]["coordinates"]], axis=0)
    _, _, n2 = frame(float(pich["properties"]["strike"]), float(pich["properties"]["dip"]))
    P = np.array([en(sf.lon, sf.lat) + np.array([0, 0, -sf.z_m / 1e3]) for sf in geom.subfaults])
    side = np.sign((P - q) @ n2)  # mismo signo que el hipocentro = lado del hipocentro
    hyp_side = np.sign((np.array([0, 0, -sp.depth]) - q) @ n2)
    mu = np.array([sf.mu_pa for sf in geom.subfaults])
    res = {}
    for key, fn in [("kin", lambda: om.kinematic_slip(cfg, fm)), ("A", lambda: om.dynamic_slip(cfg, geom, "A")),
                    ("B", lambda: om.dynamic_slip(cfg, geom, "B"))]:
        ss, ds, _ = fn()
        m = mu * np.hypot(ss, ds)
        frac = float(m[side == hyp_side].sum() / m.sum())
        res[key] = {"frac_hypocentre_side": frac, "frac_other_side": 1 - frac}
        print(f"{key}: {100 * frac:.0f} % del momento al lado del hipocentro, {100 * (1 - frac):.0f} % cruzando Pichilemu")
    res["hypocentre_side_is"] = "NE (bloque bajo el plano de Pichilemu)" if hyp_side > 0 else "SW"
    (HERE / "pichilemu_moment_split.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
