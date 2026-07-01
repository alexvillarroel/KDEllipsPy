"""
regen_data.py  —  regenera DATA/real_vel_{x,y,z} para ESTA banda
================================================================
Versión corregida de real_disp.py específica de la carpeta de la banda:
- Lee estaciones (EN ORDEN), banda y ventana del propio ``input.ctl`` de la
  banda (vía kdellipspy.ConfigParser), no de constantes hardcodeadas.
- Procesa los SAC de velocidad del evento (../SAC/VEL): detrend, bandpass
  causal, trim a la ventana, interpola a npts.
- Escribe ../DATA/real_vel_{x,y,z} en el orden del input.ctl (x=N, y=E, z=Z).
- NO reescribe input.ctl (ya está correcto).

Uso:
    python codigos/regen_data.py
"""

import glob
import sys
from pathlib import Path

import numpy as np
import obspy

BAND = Path(__file__).resolve().parent.parent  # carpeta de la banda
EVENT = BAND.parent  # raíz del evento
sys.path.insert(0, str(EVENT / "codigos"))  # para importar event.py

import kdellipspy as kde  # noqa: E402

VEL_DIR = EVENT / "SAC" / "VEL"
DATA_DIR = BAND / "DATA"
CTL = BAND / "input.ctl"
EVENT_CTL = EVENT / "event.ctl"


def procesar(station, ev, f1, f2, causal, t1, t2, npts, dt):
    """3 componentes (N,E,Z) de una estación → arrays a npts, o None si falta."""
    out = {}
    for comp in ("N", "E", "Z"):
        matches = glob.glob(str(VEL_DIR / f"{station}.*{comp}.vel.sac"))
        if not matches:
            print(f"  [!] {station}: falta {comp}; se omite la estación.")
            return None
        tr = obspy.read(matches[0])[0]
        tr.detrend("linear")
        tr.detrend("demean")
        tr.filter("bandpass", freqmin=f1, freqmax=f2, corners=4, zerophase=not causal)
        a = ev["origin_time"] + t1
        b = ev["origin_time"] + t2
        tr.trim(starttime=a, endtime=b, pad=True, fill_value=0.0)
        tr.interpolate(sampling_rate=1.0 / dt, starttime=a, npts=npts)
        out[comp] = tr.data[:npts].astype(np.float64)
    return out["N"], out["E"], out["Z"]


def main():
    cfg = kde.ConfigParser(filepath=str(CTL))
    stations = [s.name for s in cfg.stations.stations]  # EN ORDEN del input.ctl
    f1, f2 = cfg.ellipse.freq1, cfg.ellipse.freq2
    causal = not cfg.ellipse.zerophase
    od = cfg.observed_data
    npts, dt, t1, t2 = od.npts, od.delta, od.t1, od.t2
    ev = load_event(str(EVENT_CTL))

    print(
        f"Banda {f1:g}-{f2:g} Hz (causal={causal}) | npts={npts} dt={dt:g} "
        f"| ventana [{t1:g},{t2:g}] s"
    )
    print(f"Estaciones ({len(stations)}, orden input.ctl): {', '.join(stations)}")

    nx, ny, nz, usadas = [], [], [], []
    for st in stations:
        r = procesar(st, ev, f1, f2, causal, t1, t2, npts, dt)
        if r is None:
            continue
        n, e, z = r
        nx.append(n)
        ny.append(e)
        nz.append(z)
        usadas.append(st)
        print(f"  {st:6s} OK")

    assert len(usadas) == len(stations), (
        f"Faltan estaciones en SAC: pedidas {len(stations)}, generadas {len(usadas)}"
    )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for arr, suf in ((nx, "x"), (ny, "y"), (nz, "z")):
        path = DATA_DIR / f"real_vel_{suf}"
        np.savetxt(path, np.concatenate(arr), fmt="%.8e")
        rows = sum(len(a) for a in arr)
        assert rows == npts * len(usadas), (
            f"{path}: {rows} filas, esperaba {npts * len(usadas)}"
        )
        print(f"  -> {path}  ({rows} filas = {npts}×{len(usadas)})")

    print("Listo. DATA/real_vel_* regenerado; input.ctl intacto.")


if __name__ == "__main__":
    main()
