"""Reúne los resultados de Navidad 2019 en dashboard/data.json (+ figuras en base64)
para la página del dashboard."""

import base64
import json
import re
from pathlib import Path

import joblib
import numpy as np

NAV = Path(__file__).resolve().parent.parent
DYN = NAV / "dynamic"
OUT = Path(__file__).resolve().parent


def evolution(path, label, kind, note):
    r = joblib.load(path)
    mis = np.array([m.misfit for m in r.all_models], float)
    mis = np.where(mis > 1e5, np.nan, mis)  # 1e10 = evaluación fallida
    best = np.fmin.accumulate(np.nan_to_num(mis, nan=np.inf))
    return {"label": label, "kind": kind, "note": note, "misfit": np.round(mis, 4).tolist(),
            "best": np.round(best, 4).tolist(), "final": float(np.nanmin(mis)),
            "iter": [int(m.iteration) for m in r.all_models]}


def params_evolution(path):
    r = joblib.load(path)
    models = np.array([m.model for m in r.all_models], float)
    return {"names": list(r.param_names), "models": np.round(models, 3).tolist(),
            "misfit": np.round([m.misfit for m in r.all_models], 4).tolist()}


def png(path):
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def read_summary(path):
    vals = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r"\s*([^:]+?)\s{2,}(.+)$", line)
        if m:
            vals[m.group(1).strip()] = m.group(2).strip()
    return vals


def main():
    case = NAV / "isc_tomo" / "b0.05-0.30_noM14L"
    dyn_A = DYN / "na_output_isc_tomo__b0.05-0.30_noM14L"
    dyn_B = DYN / "na_output_isc_tomo__b0.05-0.30_noM14L_M0"
    data = {}

    # 1. Evolución del NA
    data["evolution"] = [
        evolution(DYN / "na_output/na_result.joblib", "Dinámico NP1 (CSN, 1ª versión)", "dyn", "60 evals, hipocentro CSN, malla 50x50"),
        evolution(DYN / "na_output_loc_mid_sin_dt0_fallida/na_result.joblib", "Dinámico loc_mid sin dt0 (fallida)", "dyn",
                  "1500 evals; desfase de hora de origen -> nada correlaciona"),
        evolution(DYN / "na_output_loc_mid/na_result.joblib", "Dinámico loc_mid + dt0", "dyn", "1000 evals, malla 50x50"),
        evolution(dyn_A / "na_result.joblib", "Dinámico A: ISC, M0 libre", "dyn", "0.05-0.30 Hz, sin M14L, malla nueva"),
        evolution(dyn_B / "na_result.joblib", "Dinámico B: ISC, M0 CSN (Mw 6.6)", "dyn", "0.05-0.30 Hz, sin M14L, malla nueva"),
    ]
    data["kin_seeds"] = [evolution(p, f"Semilla {p.stem[4:]}", "kin", "")
                         for p in sorted((case / "output_multiseed").glob("seed*.joblib"))]

    # 2. Parámetros a lo largo de la búsqueda (dinámico B y A)
    data["params_B"] = params_evolution(dyn_B / "na_result.joblib")
    data["params_A"] = params_evolution(dyn_A / "na_result.joblib")

    # 3. Matriz de bandas (cinemático multisemilla)
    bands = []
    blocks = (NAV / "isc_tomo" / "analysis.txt").read_text().split("== ")[1:]
    for b in blocks:
        name = b.splitlines()[0].strip()
        num = lambda pat: [float(x) for x in re.search(pat, b).groups()]
        mis, mstd, mbest = num(r"misfit ([\d.]+) ± ([\d.]+) \(mejor ([\d.]+)\)")
        dt, dts = num(r"dt0 ([+-][\d.]+) ± ([\d.]+)")
        vr, vrs = num(r"Vr ([\d.]+) ± ([\d.]+)")
        mw, = num(r"Mw ([\d.]+)")
        dg, = num(r"a GCMT ([\d.]+) km")
        sta = dict((k, float(v)) for k, v in re.findall(r"(\w{4}) ([\d.]+)", b.split("misfit por estación:")[1]))
        f1, f2, st = re.match(r"b([\d.]+)-([\d.]+)_(\w+)", name).groups()
        bands.append({"band": f"{f1}–{f2} Hz", "stations": "todas" if st == "all" else "sin M14L", "misfit": mis,
                      "misfit_std": mstd, "best": mbest, "dt0": dt, "dt0_std": dts, "vr": vr, "vr_std": vrs,
                      "mw": mw, "d_gcmt": dg, "per_station": sta})
    data["bands"] = bands

    # 4. Hipocentros (multisemilla, 1D original, malla 50x50) + centroides
    hyp = []
    cent = {l.split()[0]: l for l in (NAV / "vel_tests" / "centroids.txt").read_text().splitlines()}
    for name, label in [("loc_isc", "ISC"), ("loc_neic", "NEIC"), ("loc_mid", "Punto medio CSN–ISC"), ("loc_csn", "CSN")]:
        s = (NAV / name / "output_multiseed" / "summary.txt").read_text()
        mean, std, best = (float(x) for x in re.search(r"best ([\d.]+)\s+mean ([\d.]+)\s+std ([\d.]+)", s).groups()[::-1][::-1])
        best, mean, std = (float(x) for x in re.search(r"best ([\d.]+)\s+mean ([\d.]+)\s+std ([\d.]+)", s).groups())
        rows = {l.split()[0]: l.split() for l in s.splitlines() if l.startswith(("dt0", "vr"))}
        dg = float(re.search(r"GCMT:\s+([\d.]+) km", cent[name]).group(1))
        hyp.append({"label": label, "misfit": mean, "misfit_std": std, "best": best,
                    "dt0": float(rows["dt0"][3]), "dt0_std": float(rows["dt0"][4]),
                    "vr": float(rows["vr"][3]), "vr_std": float(rows["vr"][4]), "d_gcmt": dg})
    data["hypocentres"] = hyp

    # 5. Modelos de velocidad (1 semilla, malla 50x50; primera batería)
    vel = []
    for lab, path, hypo, model in [
        ("loc_mid", NAV / "loc_mid_dt0/run.log", "Punto medio", "Original"),
        ("loc_mid_tomo", NAV / "loc_mid_tomo/run_dt0.log", "Punto medio", "Regional (tomografía)"),
        ("loc_mid_hybrid", NAV / "vel_tests/loc_mid_hybrid/run_dt0.log", "Punto medio", "Híbrido"),
        ("loc_csn_orig", NAV / "vel_tests/loc_csn_orig/run_dt0.log", "CSN", "Original"),
        ("loc_csn_tomo", NAV / "vel_tests/loc_csn_tomo/run_dt0.log", "CSN", "Regional (tomografía)"),
        ("loc_csn_hybrid", NAV / "vel_tests/loc_csn_hybrid/run_dt0.log", "CSN", "Híbrido"),
    ]:
        m = re.search(r"best misfit ([\d.]+)", path.read_text())
        vel.append({"hypo": hypo, "model": model, "misfit": float(m.group(1))})
    data["velocity"] = vel

    # 6. Modelos dinámicos (resúmenes)
    data["dynamic"] = {
        "A": read_summary(dyn_A / "report/summary.txt"),
        "B": read_summary(dyn_B / "report/summary.txt"),
        "loc_mid": read_summary(DYN / "na_output_loc_mid/report/summary.txt"),
    }

    # 7. Figuras
    data["figs"] = {
        "map_B": png(dyn_B / "report/map.png"), "map_A": png(dyn_A / "report/map.png"),
        "traces_B": png(dyn_B / "report/traces.png"), "traces_A": png(dyn_A / "report/traces.png"),
        "slip_B": png(dyn_B / "report/slip.png"), "slip_A": png(dyn_A / "report/slip.png"),
        "map_old": png(DYN / "na_output_loc_mid/report/map.png"),
    }
    # 8. Deformación estática (Okada)
    ok = NAV / "okada"
    if (ok / "okada_results.json").exists():
        data["okada"] = json.loads((ok / "okada_results.json").read_text())
        for k in ("kin", "A", "B"):
            data["figs"][f"okada_{k}"] = png(ok / f"okada_{k}.png")
    (OUT / "data.json").write_text(json.dumps(data, ensure_ascii=False))
    print("ok", round((OUT / "data.json").stat().st_size / 1e6, 2), "MB")


if __name__ == "__main__":
    main()
