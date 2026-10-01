"""Modelo 1D de loc_mid desde la tomografía CSN 2019 (Potin et al., 2024, SRL).

Promedia la LENTITUD (conserva tiempos de viaje) en la caja que cubre la falla
de loc_mid y las 8 estaciones, bloque 30_38 (NonLinLoc, celdas de 4 km).
Capas de 4 km hasta 60 km: cada capa usa el promedio de los nodos en su tope
y su base. Bajo 60 km y Q: los del modelo original de loc_mid. Densidad:
Nafe-Drake de Brocher (2005) a partir de Vp. Escribe input.ctl (= loc_mid con
la sección 9 reemplazada) y enlaza DATA -> ../loc_mid/DATA.
"""

import io
import re
import zipfile
from pathlib import Path

import numpy as np

import kdellipspy as kde

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "loc_mid"
ZIP = Path("/home/alex/Drive/Publicacion_Calama_2026/Source_code/Tomography/"
           "chilean_tomography_csn_version_2019-publication_release.zip")
BLOCK = "30_38"
C111 = 10000.0 / 90.0  # NonLinLoc SIMPLE transform (km/deg)
ZMAX_KM = 60.0


def load(zf, phase):
    stem = f"chilean_tomography_csn_version_2019-publication_release/{BLOCK}/CHILE_{BLOCK}-s4_Simple.{phase}.mod"
    h = zf.read(stem + ".hdr").decode().split()
    nx, ny, nz = map(int, h[:3])
    x0, y0, z0, d = map(float, h[3:7])
    lat0, lon0 = float(h[h.index("latOrig") + 1]), float(h[h.index("longOrig") + 1])
    v = np.frombuffer(zf.read(stem + ".buf"), dtype=np.float32).reshape(nx, ny, nz)
    return v / d, (x0, y0, z0, d, lat0, lon0)  # slowness s/km


def brocher_rho(vp_kms):
    v = vp_kms
    return 1000.0 * (1.6612 * v - 0.4721 * v**2 + 0.0671 * v**3 - 0.0043 * v**4 + 0.000106 * v**5)


def main():
    cfg = kde.ConfigParser(str(BASE / "input.ctl"))
    geo = kde.AxitraForwardModel.from_config(cfg).build_geometry()
    lats = [sf.lat for sf in geo.subfaults] + [s.latitude for s in cfg.stations.stations]
    lons = [sf.lon for sf in geo.subfaults] + [s.longitude for s in cfg.stations.stations]
    la, lo = np.meshgrid(np.linspace(min(lats), max(lats), 40), np.linspace(min(lons), max(lons), 40))
    la, lo = la.ravel(), lo.ravel()

    prof = {}
    with zipfile.ZipFile(ZIP) as zf:
        for ph in ("P", "S"):
            s, (x0, y0, z0, d, lat0, lon0) = load(zf, ph)
            ix = np.rint(((lo - lon0) * C111 * np.cos(np.radians(la)) - x0) / d).astype(int)
            iy = np.rint(((la - lat0) * C111 - y0) / d).astype(int)
            prof[ph] = s[ix, iy, :].mean(axis=0)
            z = z0 + d * np.arange(s.shape[2])

    rows = []
    for top in np.arange(0.0, ZMAX_KM, 4.0):
        k = int(np.argmin(np.abs(z - top)))
        vp = 1.0 / prof["P"][k:k + 2].mean()
        vs = 1.0 / prof["S"][k:k + 2].mean()
        q = 300.0 if top < 50 else 1000.0
        rows.append((top * 1e3, vp * 1e3, vs * 1e3, brocher_rho(vp), q, q))
    rows += [(l.thickness, l.vp, l.vs, l.rho, l.qp, l.qs) for l in cfg.velocity_model.layers if l.thickness >= ZMAX_KM * 1e3]

    text = (BASE / "input.ctl").read_text()
    body = "\n".join(f" {t:10.1f} {vp:14.2f} {vs:12.2f} {rho:10.1f} {qp:11.1f} {qs:9.1f}" for t, vp, vs, rho, qp, qs in rows)
    new, n = re.subn(r"( Number of layers\s*:\s*)\d+\n(?:\s*[-0-9.eE]+(?:\s+[-0-9.eE]+){5}\s*\n)+",
                     lambda m: f"{m.group(1)}{len(rows)}\n{body}\n", text)
    assert n == 1, "no se encontró la sección 9 del input.ctl"
    (HERE / "input.ctl").write_text(new)
    if not (HERE / "DATA").exists():
        (HERE / "DATA").symlink_to(BASE / "DATA")
    print(body)


if __name__ == "__main__":
    main()
