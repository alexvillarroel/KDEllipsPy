"""
real_disp.py  (adaptado para 2020-06-03)
========================================
Lee los SAC de banda ancha generados por ``txt_to_sac.py`` (SAC/VEL y SAC/DISP),
aplica el filtro fino de inversion (0.02-0.1 Hz causal), recorta a la ventana del
evento y remuestrea, y escribe los archivos planos:
    DATA/real_vel_{x,y,z}    (velocidad,     m/s)   <- desde SAC/VEL
    DATA/real_disp_{x,y,z}   (desplazamiento, m)    <- desde SAC/DISP
para las estaciones SELECTED, EN ESE ORDEN (x=N, y=E, z=Z).

Ademas actualiza, en input.ctl, el bloque de estaciones (seccion 8) y el campo
Units, dejando todo consistente con el orden de los archivos planos.

Uso:
    python codigos/real_disp.py
"""
import os
import glob

import numpy as np
import obspy

from event import load_event

# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTL_FILE = os.path.join(BASE, "event.ctl")
VEL_DIR = os.path.join(BASE, "SAC", "VEL")
DISP_DIR = os.path.join(BASE, "SAC", "DISP")
DATA_DIR = os.path.join(BASE, "DATA")
INPUT_CTL = os.path.join(BASE, "input.ctl")

# Set azimutal aprobado (orden por azimut). El orden es el de los archivos planos
# y el del bloque de estaciones del input.ctl.
# Orden por DISTANCIA epicentral (cerca -> lejos): A22F 50, A08F 54, A05F 90,
# A18F 97, A23F 102, PB19 107, PB06 122, A09F 128, A06F 153 km
SELECTED = ["A22F", "A08F", "A05F", "A18F", "A23F", "PB19", "PB06", "A09F", "A06F"]

# Ventana / muestreo (calama2020).
T1, T2 = 0.0, 128.0
NPTS, DT = 128, 1.0
# Banda de inversion + fase.
F1, F2 = 0.1, 0.3
CAUSAL = True          # True -> filtro causal (zerophase=False), consistente con input.ctl


def procesar_estacion(station, ev, src_dir, kind):
    """Procesa las 3 componentes de una estacion desde SAC (kind='vel'|'disp').
    Devuelve (n, e, z, lat, lon) a NPTS muestras, o None si falta algo."""
    comps, coords = {}, None
    for comp in ("N", "E", "Z"):
        matches = glob.glob(os.path.join(src_dir, f"{station}.*{comp}.{kind}.sac"))
        if not matches:
            print(f"  [!] {station}: falta {comp} ({kind}); se omite.")
            return None
        tr = obspy.read(matches[0])[0]
        tr.detrend("linear")
        tr.detrend("demean")
        tr.filter("bandpass", freqmin=F1, freqmax=F2, corners=4, zerophase=not CAUSAL)
        t1 = ev["origin_time"] + T1
        t2 = ev["origin_time"] + T2
        tr.trim(starttime=t1, endtime=t2, pad=True, fill_value=0.0)
        tr.interpolate(sampling_rate=1.0 / DT, starttime=t1, npts=NPTS)
        comps[comp] = tr.data[:NPTS].astype(np.float64)
        coords = (float(tr.stats.sac.stla), float(tr.stats.sac.stlo))
    return comps["N"], comps["E"], comps["Z"], coords[0], coords[1]


def generar(ev, kind, src_dir, prefix):
    nx, ny, nz, est = [], [], [], []
    print(f"\n{prefix}  (banda {F1}-{F2} Hz, causal={CAUSAL}):")
    for station in SELECTED:
        r = procesar_estacion(station, ev, src_dir, kind)
        if r is None:
            continue
        n, e, z, la, lo = r
        nx.append(n); ny.append(e); nz.append(z); est.append((station, la, lo))
        print(f"  {station:6s} lat {la:8.3f} lon {lo:8.3f} OK")
    os.makedirs(DATA_DIR, exist_ok=True)
    np.savetxt(os.path.join(DATA_DIR, f"{prefix}_x"), np.concatenate(nx), fmt="%.8e")
    np.savetxt(os.path.join(DATA_DIR, f"{prefix}_y"), np.concatenate(ny), fmt="%.8e")
    np.savetxt(os.path.join(DATA_DIR, f"{prefix}_z"), np.concatenate(nz), fmt="%.8e")
    print(f"  -> {len(est)} est x {NPTS} = {len(est)*NPTS} filas en {prefix}_*")
    return est


def actualizar_input_ctl(est, units=1):
    """Reemplaza el bloque de estaciones (sec. 8) y el campo Units en input.ctl."""
    with open(INPUT_CTL) as f:
        lines = f.readlines()
    out, i = [], 0
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if s.startswith("Units"):
            out.append(f" Units                        (1:disp, 2:vel)   :       {units}\n")
            i += 1
            continue
        if s.startswith("Number of stations"):
            out.append(f" Number of stations                             :       {len(est)}\n")
            for name, la, lo in est:
                out.append(f"{la:.4f} {lo:.4f} 0.0000 {name:<6s} 1 1 1\n")
            # saltar las lineas de estaciones viejas (hasta blanco o seccion '#')
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].lstrip().startswith("#"):
                i += 1
            continue
        out.append(ln)
        i += 1
    with open(INPUT_CTL, "w") as f:
        f.writelines(out)
    print(f"\ninput.ctl: {len(est)} estaciones, units={units} (desplazamiento)")


def main():
    ev = load_event(CTL_FILE)
    print(f"Evento: {ev['origin_time']} | {ev['lat']},{ev['lon']} {ev['depth']} km")
    generar(ev, "vel", VEL_DIR, "real_vel")
    est = generar(ev, "disp", DISP_DIR, "real_disp")
    actualizar_input_ctl(est, units=1)
    print("\nListo. real_vel_* y real_disp_* en DATA/, input.ctl actualizado.")


if __name__ == "__main__":
    main()
