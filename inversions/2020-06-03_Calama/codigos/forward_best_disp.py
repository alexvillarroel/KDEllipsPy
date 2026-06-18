#!/usr/bin/env python3
"""
forward_best_disp.py
====================
Forward de la MEJOR SOLUCION de la corrida previa (que fue en VELOCIDAD) pero
evaluada en DESPLAZAMIENTO (units=1) contra los observados real_disp.
Reporta el misfit en desplazamiento y grafica el ajuste rotado R/T/Z + ventana.
"""
from pathlib import Path
import numpy as np
import kdellipspy as kde
from kdellipspy.core.plotting import plot_waveform_fit

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
OUT = BASE / "output" / "waveform_fit_best_disp.png"

BEST = np.array([9.8454, 4.5649, 0.4062, 0.6339, 0.9116, 3.4557, 1.5239])  # mejor sol. (vel run)


def main():
    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))   # units=1, real_disp, causal, ventana 20s
    print(f"units={cfg.observed_data.units} (1=disp) | source_type={cfg.ellipse.source_type} | "
          f"banda {cfg.ellipse.freq1}-{cfg.ellipse.freq2} | ventana {cfg.inversion_process.misfit_time_window}s")
    obs, t = kde.load_and_filter_observed_data(
        input_ctl_path=str(INPUT_CTL), data_dir=str(BASE / "DATA"), prefer_raw=False,
        freq1=cfg.ellipse.freq1, freq2=cfg.ellipse.freq2)
    azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")

    model = kde.NAInversionModel(config=cfg, observed_waveforms=obs, time_array=t, azi_times_array=azt)
    model.use_green_cache = False
    misfit, synth = model._evaluate_model(BEST)
    print(f"\n  misfit (DESPLAZAMIENTO, ventana 20s) = {misfit:.4f}")

    names = [s.name for s in cfg.stations.stations]
    flags = np.array([[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations], bool)
    plot_waveform_fit(obs, synth, t, names, misfit=misfit, show=False, save_path=str(OUT),
                      azimuths=azt[:, 0], tp_s=azt[:, 1], ts_s=azt[:, 2], window_s=20.0,
                      station_flags=flags, rotate=True, mark_windows=True)
    print(f"figura -> {OUT}")


if __name__ == "__main__":
    main()
