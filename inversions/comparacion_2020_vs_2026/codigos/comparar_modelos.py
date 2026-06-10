#!/usr/bin/env python3
"""
comparar_modelos.py
===================
Compara los mejores modelos de las inversiones cinematicas de los terremotos de
Calama 2020-06-03 y 2026-05-25. Carga los dos `inversion_result.joblib` ya
existentes (NO re-invierte nada) e imprime una tabla lado a lado con:

  - metadatos del evento (lat, lon, profundidad, mecanismo strike/dip/rake)
  - parametros del mejor modelo de la elipse (a1, a2, theta, np, tp, dmax, vr)
  - cantidades derivadas (angulo en grados, area de la elipse, M0, Mw)
  - condiciones de la inversion (banda de frecuencia, plano de falla, N estaciones,
    misfit, iteracion del mejor modelo)

AVISO: las dos inversiones NO usaron la misma banda de frecuencias ni el mismo
numero de estaciones, asi que el *misfit* NO es comparable directamente entre
eventos. Los parametros de la elipse (geometria, slip, Vr) y el M0/Mw si lo son
y son el foco de la comparacion. Se imprime la banda de cada uno para tenerlo
presente.

Uso:
    python codigos/comparar_modelos.py
"""

from pathlib import Path

import joblib
import numpy as np

import kdellipspy as kde  # noqa: F401  (asegura que los objetos del joblib resuelvan)
from kdellipspy.core.forward_model import AxitraForwardModel

# ---------------------------------------------------------------------------- #
BASE = Path(__file__).resolve().parent.parent          # comparacion_2020_vs_2026/
INVERSIONS = BASE.parent                                # inversions/
OUTPUT = BASE / "output"

EVENTOS = {
    "Calama 2020": INVERSIONS / "2020-06-03",
    "Calama 2026": INVERSIONS / "2026-05-25_Calama",
}

PARAM_NAMES = ["a1 (km)", "a2 (km)", "theta (x pi)", "np (frac)",
               "tp (x 2pi)", "dmax (m)", "vr (km/s)"]


def load_event_ctl(event_dir):
    """Lee event.ctl -> dict con lat, lon, depth, strike, dip, rake, utc."""
    line = None
    with open(event_dir / "event.ctl") as fh:
        for raw in fh:
            raw = raw.strip()
            if raw and not raw.startswith("#"):
                line = raw
                break
    if line is None:
        return {}
    parts = line.split()
    keys = ["lat", "lon", "depth", "strike", "dip", "rake", "utc"]
    out = {}
    for k, v in zip(keys, parts):
        out[k] = v if k == "utc" else float(v)
    return out


def collect(name, event_dir):
    """Carga el joblib + event.ctl de un evento y arma un dict de metricas."""
    jb = event_dir / "output" / "inversion_result.joblib"
    result = joblib.load(jb)
    model = np.asarray(result.best_model.model, dtype=float)
    cfg = result.config
    ev = load_event_ctl(event_dir)

    # cantidades derivadas
    a1, a2, theta, npf, tp, dmax, vr = model
    theta_deg = theta * 180.0          # theta esta en unidades de pi rad
    area_km2 = np.pi * a1 * a2          # area de la elipse (semiejes en km)

    # M0 / Mw desde el forward (no necesita datos observados)
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.build_geometry_with_ellipse_slip(model)
    m0, mw = fm.estimate_total_moment_and_mw(model, geometry=geom)

    return {
        "name": name,
        "event": ev,
        "model": model,
        "theta_deg": theta_deg,
        "area_km2": area_km2,
        "m0": m0,
        "mw": mw,
        "misfit": result.best_model.misfit,
        "iteration": result.best_model.iteration,
        "freq1": cfg.ellipse.freq1,
        "freq2": cfg.ellipse.freq2,
        "lx_km": cfg.fault_plane.lx / 1000.0,
        "nx": cfg.fault_plane.nx,
        "nsta": len(cfg.stations.stations),
        "npts": result.observed.shape[-1] if result.observed is not None else None,
    }


def _row(label, va, vb, fmt="{}"):
    """Formatea una fila de la tabla de 3 columnas."""
    sa = fmt.format(va) if va is not None else "-"
    sb = fmt.format(vb) if vb is not None else "-"
    return f"  {label:<26s} | {sa:>18s} | {sb:>18s}"


def main():
    data = {name: collect(name, d) for name, d in EVENTOS.items()}
    a = data["Calama 2020"]
    b = data["Calama 2026"]

    lines = []
    lines.append("=" * 70)
    lines.append("  COMPARACION DE INVERSIONES CINEMATICAS  (mejor modelo)")
    lines.append("=" * 70)
    lines.append(f"  {'':<26s} | {a['name']:>18s} | {b['name']:>18s}")
    lines.append("  " + "-" * 66)

    # --- Evento ----------------------------------------------------------- #
    lines.append("  [ EVENTO ]")
    ea, eb = a["event"], b["event"]
    lines.append(_row("Fecha UTC", ea.get("utc"), eb.get("utc")))
    lines.append(_row("Latitud", ea.get("lat"), eb.get("lat"), "{:.3f}"))
    lines.append(_row("Longitud", ea.get("lon"), eb.get("lon"), "{:.3f}"))
    lines.append(_row("Profundidad (km)", ea.get("depth"), eb.get("depth"), "{:.1f}"))
    lines.append(_row("Strike / Dip / Rake",
                      f"{ea.get('strike'):.0f}/{ea.get('dip'):.0f}/{ea.get('rake'):.0f}",
                      f"{eb.get('strike'):.0f}/{eb.get('dip'):.0f}/{eb.get('rake'):.0f}"))

    # --- Parametros del mejor modelo -------------------------------------- #
    lines.append("  " + "-" * 66)
    lines.append("  [ PARAMETROS ELIPSE (mejor modelo) ]")
    for i, pn in enumerate(PARAM_NAMES):
        lines.append(_row(pn, a["model"][i], b["model"][i], "{:.4f}"))

    # --- Derivados -------------------------------------------------------- #
    lines.append("  " + "-" * 66)
    lines.append("  [ DERIVADOS ]")
    lines.append(_row("Angulo theta (grados)", a["theta_deg"], b["theta_deg"], "{:.1f}"))
    lines.append(_row("Area elipse (km^2)", a["area_km2"], b["area_km2"], "{:.1f}"))
    lines.append(_row("M0 (N.m)", a["m0"], b["m0"], "{:.3e}"))
    lines.append(_row("Mw", a["mw"], b["mw"], "{:.2f}"))

    # --- Condiciones de inversion ----------------------------------------- #
    lines.append("  " + "-" * 66)
    lines.append("  [ CONDICIONES DE INVERSION ]")
    lines.append(_row("Banda freq (Hz)",
                      f"{a['freq1']:.3f}-{a['freq2']:.3f}",
                      f"{b['freq1']:.3f}-{b['freq2']:.3f}"))
    lines.append(_row("Plano Lx=Ly (km)", a["lx_km"], b["lx_km"], "{:.0f}"))
    lines.append(_row("Subfallas Nx=Ny", a["nx"], b["nx"], "{:d}"))
    lines.append(_row("N estaciones", a["nsta"], b["nsta"], "{:d}"))
    lines.append(_row("N puntos por traza", a["npts"], b["npts"], "{:d}"))
    lines.append(_row("Misfit (NO comparable)", a["misfit"], b["misfit"], "{:.4f}"))
    lines.append(_row("Iteracion mejor modelo", a["iteration"], b["iteration"], "{:d}"))
    lines.append("=" * 70)
    lines.append("  NOTA: el misfit y N estaciones difieren entre eventos; el misfit")
    lines.append("  NO es comparable directo. Geometria, slip, Vr, M0 y Mw si lo son.")
    lines.append("=" * 70)

    text = "\n".join(lines)
    print(text)

    OUTPUT.mkdir(exist_ok=True)
    out_txt = OUTPUT / "comparacion_modelos.txt"
    out_txt.write_text(text + "\n")
    print(f"\n  -> tabla guardada en {out_txt}")


if __name__ == "__main__":
    main()
