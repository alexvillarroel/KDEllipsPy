#!/usr/bin/env python3
"""
calc_azi_times_dynamic.py
=========================
Calcula azi_times.txt para la inversion DINAMICA (legacy FD3D) usando obspy,
reutilizando la MISMA convencion del lado cinematico (kde.build_azi_times_array:
TauPyModel iasp91, azimut en radianes, tP/tS con shift -1 s).

Columnas (formato que lee calcmisfit.f):  azi[rad]  tP[s]  tS[s]
Una linea por estacion, EN EL ORDEN del input.ctl / real_disp.

Salida: Dynamic_inversion/Evento/azi_times.txt
Uso:    python codigos/calc_azi_times_dynamic.py
"""
from pathlib import Path
import numpy as np
import kdellipspy as kde

BASE = Path(__file__).resolve().parent.parent
INPUT_CTL = BASE / "input.ctl"
OUT = BASE / "Dynamic_inversion" / "Evento" / "azi_times.txt"


def main():
    cfg = kde.ConfigParser(filepath=str(INPUT_CTL))
    azt = kde.build_azi_times_array(config=cfg, model_name="iasp91")  # [azi_rad, tP, tS]
    names = [s.name for s in cfg.stations.stations]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        for (azi, tp, ts), name in zip(azt, names):
            f.write(f"{azi:.4f} {tp:.1f} {ts:.1f}\n")

    print(f"azi_times.txt ({len(names)} estaciones) -> {OUT}")
    for (azi, tp, ts), name in zip(azt, names):
        print(f"  {name:5s}  azi={np.degrees(azi):6.1f} deg  tP={tp:6.1f}s  tS={ts:6.1f}s")


if __name__ == "__main__":
    main()
