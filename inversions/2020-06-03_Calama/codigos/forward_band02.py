#!/usr/bin/env python3
"""
forward_band02.py
=================
Forward del MEJOR modelo de la inversion (output/inversion_result.joblib) pero en
banda 0.06-0.2 Hz (en vez de 0.06-0.15 de la inversion), para ver si se captura
mas señal. NO re-invierte: usa los 7 parametros convergidos.

- Re-procesa los observados desde SAC/VEL a 0.06-0.2 Hz (reusa
  real_disp.procesar_estacion con F2 sobreescrito) — los archivos planos estaban
  a 0.15, por eso hay que re-filtrar desde la velocidad.
- Corre el forward; el sintetico se filtra a 0.06-0.2 (cfg.ellipse.freq2).
- Grafica el ajuste rotado R/T/Z con ventanas y reporta el misfit en esa banda.

Uso:
    python codigos/forward_band02.py
"""

from pathlib import Path
import numpy as np
import joblib

import kdellipspy as kde
from kdellipspy.core.plotting import plot_waveform_fit
from event import load_event
import real_disp                      # reusar el procesado de observados

BASE = Path(__file__).resolve().parent.parent
JOBLIB = BASE / "output" / "inversion_result.joblib"
CTL_FILE = BASE / "event.ctl"
OUT = BASE / "output" / "waveform_fit_006-02.png"

F2_NEW = 0.2


def main():
    ev = load_event(str(CTL_FILE))
    r = joblib.load(JOBLIB)
    best7 = np.asarray(r.best_model.model, float)
    cfg = r.config

    # --- Observados re-procesados a 0.06-0.2 Hz (mismo pipeline, F2 distinto) --
    real_disp.F2 = F2_NEW
    print(f"Re-procesando observados a {real_disp.F1}-{real_disp.F2} Hz "
          f"(causal={real_disp.CAUSAL})...")
    nx, ny, nz, names_proc = [], [], [], []
    for st in real_disp.SELECTED:
        res = real_disp.procesar_estacion(st, ev)
        if res is None:
            print(f"  [!] {st} omitida"); continue
        n, e, z, _, _ = res
        nx.append(n); ny.append(e); nz.append(z); names_proc.append(st)
    obs = np.stack([np.array(nx), np.array(ny), np.array(nz)], axis=1)  # (nsta,3,npts)
    time = np.arange(obs.shape[2]) * real_disp.DT
    print(f"observed {obs.shape}  estaciones {names_proc}")

    # Coherencia de orden con el config.
    cfg_names = [s.name for s in cfg.stations.stations]
    assert names_proc == cfg_names, f"orden distinto: {names_proc} vs {cfg_names}"

    # --- Forward en banda 0.06-0.2 -------------------------------------------
    cfg.ellipse.freq2 = F2_NEW
    azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")
    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs,
                                 time_array=time, azi_times_array=azt)
    model.use_green_cache = True
    print("Forward (Green hasta Nyquist 1 Hz, sintetico filtrado a 0.06-0.2)...")
    misfit, synth = model._evaluate_model(best7)
    print(f"\n  misfit (0.06-0.2 Hz, ventana 20s) = {misfit:.4f}")
    print(f"  (referencia inversion 0.06-0.15 Hz = {r.best_model.misfit:.4f})")

    # --- Plot rotado R/T/Z + ventana -----------------------------------------
    flags = np.array([[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations], bool)
    plot_waveform_fit(obs, synth, time, cfg_names, misfit=misfit, show=False,
                      save_path=str(OUT),
                      azimuths=azt[:, 0], tp_s=azt[:, 1], ts_s=azt[:, 2],
                      window_s=20.0, station_flags=flags,
                      rotate=True, mark_windows=True)
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
