"""Casos de prueba hipocentro x modelo de velocidades (cinemático + dt0).

Para cada hipocentro base (../loc_*) escribe vel_tests/<base>_<vel>/input.ctl
(= input.ctl base con la sección 9 reemplazada) y DATA -> ../../<base>/DATA:
  orig   : modelo original (sin cambios)
  tomo   : promedio de lentitud de la tomografía CSN 2019 en la caja falla+estaciones
  hybrid : original sobre HYBRID_Z_KM, tomografía debajo
Correr con:  KIN_CASE=vel_tests/<caso> python ../loc_mid_dt0/run_kin_dt0.py
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "loc_mid_tomo"))
import kdellipspy as kde  # noqa: E402
from make_tomo_model import BLOCK, C111, ZIP, ZMAX_KM, brocher_rho, load  # noqa: E402

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
HYBRID_Z_KM = 12.0
CASES = {"loc_csn": ("orig", "tomo", "hybrid"), "loc_mid": ("hybrid",), "loc_isc": ("orig",), "loc_neic": ("orig",)}


def tomo_rows(cfg, zf):
    geo = kde.AxitraForwardModel.from_config(cfg).build_geometry()
    lats = [sf.lat for sf in geo.subfaults] + [s.latitude for s in cfg.stations.stations]
    lons = [sf.lon for sf in geo.subfaults] + [s.longitude for s in cfg.stations.stations]
    la, lo = np.meshgrid(np.linspace(min(lats), max(lats), 40), np.linspace(min(lons), max(lons), 40))
    la, lo = la.ravel(), lo.ravel()
    prof = {}
    for ph in ("P", "S"):
        s, (x0, y0, z0, d, lat0, lon0) = load(zf, ph)
        ix = np.rint(((lo - lon0) * C111 * np.cos(np.radians(la)) - x0) / d).astype(int)
        iy = np.rint(((la - lat0) * C111 - y0) / d).astype(int)
        prof[ph] = s[ix, iy, :].mean(axis=0)
        z = z0 + d * np.arange(s.shape[2])
    rows = []
    for top in np.arange(0.0, ZMAX_KM, 4.0):
        k = int(np.argmin(np.abs(z - top)))
        vp, vs = 1.0 / prof["P"][k:k + 2].mean(), 1.0 / prof["S"][k:k + 2].mean()
        q = 300.0 if top < 50 else 1000.0
        rows.append((top * 1e3, vp * 1e3, vs * 1e3, brocher_rho(vp), q, q))
    return rows


def write_case(base, vel, rows):
    out = HERE / f"{base}_{vel}"
    out.mkdir(exist_ok=True)
    text = (NAV / base / "input.ctl").read_text()
    if rows is not None:
        body = "\n".join(f" {t:10.1f} {vp:14.2f} {vs:12.2f} {rho:10.1f} {qp:11.1f} {qs:9.1f}" for t, vp, vs, rho, qp, qs in rows)
        text, n = re.subn(r"( Number of layers\s*:\s*)\d+\n(?:\s*[-0-9.eE]+(?:\s+[-0-9.eE]+){5}\s*\n)+",
                          lambda m: f"{m.group(1)}{len(rows)}\n{body}\n", text)
        assert n == 1
    (out / "input.ctl").write_text(text)
    if not (out / "DATA").exists():
        (out / "DATA").symlink_to(NAV / base / "DATA")
    print(f"{out.name}: {len(rows) if rows else 'orig'} capas")


def main():
    import zipfile
    with zipfile.ZipFile(ZIP) as zf:
        for base, vels in CASES.items():
            cfg = kde.ConfigParser(str(NAV / base / "input.ctl"))
            orig = [(l.thickness, l.vp, l.vs, l.rho, l.qp, l.qs) for l in cfg.velocity_model.layers]
            tomo = tomo_rows(cfg, zf) + [r for r in orig if r[0] >= ZMAX_KM * 1e3]
            for vel in vels:
                if vel == "orig":
                    write_case(base, vel, None)
                elif vel == "tomo":
                    write_case(base, vel, tomo)
                else:
                    write_case(base, vel, [r for r in orig if r[0] < HYBRID_Z_KM * 1e3] +
                               [r for r in tomo if r[0] >= HYBRID_Z_KM * 1e3])


if __name__ == "__main__":
    main()
