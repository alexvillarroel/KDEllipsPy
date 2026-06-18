#!/usr/bin/env python3
"""
test_plano_grande.py
====================
Prueba (NO re-inversion): toma los parametros que convergieron en la corrida de
noche (output/inversion_result.joblib) y, dejandolos FIJOS, agranda el plano de
falla a 64 km (Nx=Ny=60, hipocentro centrado Hx=Hy=32 km). Hace UN solo forward y
grafica la elipse completa + el ajuste de formas de onda.

Objetivo: ver la elipse SIN recortar y confirmar que el misfit/sinteticas no
cambian respecto a la corrida original (misfit 0.4679). Si el misfit se mantiene
(~1%), la parte que quedaba cortada tenia slip despreciable (cola gaussiana) y la
solucion es confiable; el recorte era solo un artefacto de geometria/plano chico.

No toca input.ctl ni real_disp.py: la geometria grande vive solo en memoria.
Las figuras llevan sufijo _plano64 para no pisar las canonicas de output/.

Uso:
    python codigos/test_plano_grande.py
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import joblib
import numpy as np

import kdellipspy as kde
from kdellipspy.core.forward_model import AxitraForwardModel
from kdellipspy.core.plotting import plot_waveform_fit

# ---------------------------------------------------------------------------- #
BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
DATA_DIR = BASE / "DATA"
OUTPUT = BASE / "output"
JOBLIB = OUTPUT / "inversion_result.joblib"

# Plano grande para la prueba (half-plano 32 km > los ~31 km que alcanza el
# borde lejano de la elipse con np=1 y a2~15.7 -> garantiza no cortar).
LX = LY = 64000.0
HX = HY = 32000.0
NX = NY = 60

MISFIT_REF = 0.4679  # mejor misfit de la corrida de noche (plano 32 km)


def main():
    print(f"Cargando best_model de {JOBLIB} ...")
    result = joblib.load(JOBLIB)
    best7 = np.asarray(result.best_model.model, dtype=float)
    names = ["a1", "a2", "theta", "np", "tp", "dmax", "vr"]
    print("  best_model:")
    for n, v in zip(names, best7):
        print(f"    {n:<6s} = {v:.4f}")

    # --- Config con el plano agrandado (en memoria) ------------------------- #
    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    cfg.fault_plane.lx = LX
    cfg.fault_plane.ly = LY
    cfg.fault_plane.hx = HX
    cfg.fault_plane.hy = HY
    cfg.fault_plane.nx = NX
    cfg.fault_plane.ny = NY
    dstk = LX / NX
    print(f"\nPlano agrandado: Lx=Ly={LX/1000:.0f} km, Hx=Hy={HX/1000:.0f} km, "
          f"Nx=Ny={NX} (subfalla ~{dstk/1000:.2f} km)")

    # --- Datos observados + tiempos P/S (igual que run_inversion.py) -------- #
    observed, time_array = kde.load_and_filter_observed_data(
        input_ctl_path=str(INPUT_CTL),
        data_dir=str(DATA_DIR),
        prefer_raw=False,
        freq1=cfg.ellipse.freq1,
        freq2=cfg.ellipse.freq2,
    )
    azi_times = kde.build_azi_times_array(config=cfg, model_name="iasp91")
    print(f"observed shape = {observed.shape}  | azi_times = {azi_times.shape}")

    # --- Un solo forward con el modelo fijo --------------------------------- #
    model = kde.NAInversionModel(
        config=cfg,
        observed_waveforms=observed,
        time_array=time_array,
        azi_times_array=azi_times,
    )
    model.use_green_cache = True
    print("\nCorriendo UN forward con el plano de 64 km (calcula green una vez)...")
    misfit, synth = model._evaluate_model(best7)
    print(f"\n  misfit (plano 64 km) = {misfit:.6f}")
    print(f"  misfit (referencia)  = {MISFIT_REF:.6f}")
    dpc = 100.0 * (misfit - MISFIT_REF) / MISFIT_REF
    print(f"  diferencia           = {dpc:+.2f} %  "
          f"({'OK, recorte despreciable' if abs(dpc) < 1.5 else 'revisar'})")

    # --- Figuras (sufijo _plano64) ------------------------------------------ #
    fm = AxitraForwardModel.from_config(cfg)
    geom = fm.build_geometry_with_ellipse_slip(best7)
    m0, mw = fm.estimate_total_moment_and_mw(best7, geometry=geom)
    print(f"\n  M0 = {m0:.3e} N.m   Mw = {mw:.2f}")

    fig_ell = OUTPUT / "ellipse_fit_plano64.png"
    geom.plot(title="2D Slip Distribution (plano 64 km, params fijos)",
              show=False, save_path=str(fig_ell))
    print(f"  figura elipse  -> {fig_ell}")

    station_names = [s.name for s in cfg.stations.stations]
    fig_wf = OUTPUT / "waveform_fit_plano64.png"
    plot_waveform_fit(observed, synth, time_array, station_names,
                      misfit=misfit, show=False, save_path=str(fig_wf))
    print(f"  figura formas  -> {fig_wf}")
    print("\nListo.")


if __name__ == "__main__":
    main()
