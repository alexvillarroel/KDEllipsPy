"""Casos ISC + 1D regional (tomografía CSN 2019) + malla desplazada, por banda
de frecuencia y con/sin M14L. Cada caso: input.ctl, SAC/codigos -> ../../loc_isc,
y DATA reprocesado con kde-prep en su banda (acc cruda -> desplazamiento).

Uso: python make_band_cases.py
"""

import re
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAV = HERE.parent
BASE = NAV / "loc_isc"
sys.path.insert(0, str(NAV / "vel_tests"))
sys.path.insert(0, str(NAV / "loc_mid_tomo"))
import kdellipspy as kde  # noqa: E402
from make_cases import ZIP, ZMAX_KM, tomo_rows  # noqa: E402

# Malla: hipocentro arriba (Hy=10 km) porque el slip se va hacia abajo en el
# dip; ±16 km en strike. 2 km por subfalla; Lx/dh y Ly/dh múltiplos de nx, ny.
FAULT = {"Length along strike (Lx)": 32000.0, "Length along dip (Ly)": 40000.0,
         "Hypocenter position strike (Hx)": 16000.0, "Hypocenter position dip (Hy)": 10000.0,
         "Number of subfaults along strike (Nx)": 16, "Number of subfaults along dip (Ny)": 20}
BANDS = [(0.06, 0.15), (0.05, 0.20), (0.05, 0.25), (0.05, 0.30)]
DROP = "M14L"
KDE_PREP = Path(sys.executable).parent / "kde-prep"


def set_value(text, label, value):
    pat = re.compile(r"^(\s*" + re.escape(label) + r"[^:\n]*:\s*)\S+", re.M)
    text, n = pat.subn(lambda m: f"{m.group(1)}{value}", text)
    assert n == 1, label
    return text


def set_velocity(text, rows):
    body = "\n".join(f" {t:10.1f} {vp:14.2f} {vs:12.2f} {rho:10.1f} {qp:11.1f} {qs:9.1f}" for t, vp, vs, rho, qp, qs in rows)
    text, n = re.subn(r"( Number of layers\s*:\s*)\d+\n(?:\s*[-0-9.eE]+(?:\s+[-0-9.eE]+){5}\s*\n)+",
                      lambda m: f"{m.group(1)}{len(rows)}\n{body}\n", text)
    assert n == 1
    return text


def drop_station(text, name):
    lines = text.splitlines(keepends=True)
    lines = [l for l in lines if not re.match(rf"^\s*\S+\s+\S+\s+\S+\s+{name}\s", l)]
    text = "".join(lines)
    n = int(re.search(r"Number of stations\s*:\s*(\d+)", text).group(1))
    return set_value(text, "Number of stations", n - 1)


def main():
    text = (BASE / "input.ctl").read_text()
    for label, value in FAULT.items():
        text = set_value(text, label, value)
    tmp = HERE / "_geom_input.ctl"
    tmp.write_text(text)
    cfg = kde.ConfigParser(str(tmp))  # caja falla+estaciones con la malla nueva
    orig = [(l.thickness, l.vp, l.vs, l.rho, l.qp, l.qs) for l in cfg.velocity_model.layers]
    with zipfile.ZipFile(ZIP) as zf:
        rows = tomo_rows(cfg, zf) + [r for r in orig if r[0] >= ZMAX_KM * 1e3]
    tmp.unlink()
    text = set_velocity(text, rows)

    for f1, f2 in BANDS:
        for stations in ("all", f"no{DROP}"):
            case = HERE / f"b{f1:.2f}-{f2:.2f}_{stations}"
            case.mkdir(exist_ok=True)
            t = set_value(set_value(text, "Frequency 1 (Freq1)", f"{f1:.6f}"), "Frequency 2 (Freq2)", f"{f2:.6f}")
            if stations != "all":
                t = drop_station(t, DROP)
            (case / "input.ctl").write_text(t)
            for link in ("SAC", "codigos"):
                if not (case / link).exists():
                    (case / link).symlink_to(BASE / link)
            r = subprocess.run([str(KDE_PREP), str(case)], capture_output=True, text=True)
            ok = (case / "DATA" / "real_disp_z").exists()
            print(f"{case.name}: kde-prep {'ok' if r.returncode == 0 and ok else 'FALLÓ'}"
                  + ("" if r.returncode == 0 else f"\n{r.stdout[-800:]}\n{r.stderr[-800:]}"))


if __name__ == "__main__":
    main()
