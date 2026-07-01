"""
real_disp.py
============
Prepara los datos para la inversion cinematica con KDEllipsPy (estilo del
``real_disp.py`` de copiapo) y escribe el ``input.ctl``.

Para las estaciones seleccionadas (subconjunto azimutal a < 200 km):
  1. Lee la velocidad (SAC/VEL).
  2. detrend (linear+demean) -> integra a desplazamiento (cumtrapz).
  3. Filtro pasabanda CAUSAL 0.06-0.15 Hz (zerophase=False), igual que el
     que usara KDEllipsPy en los sinteticos (input.ctl: Zerophase filter = 0).
  4. Recorta a la ventana [origen, origen+128 s] e interpola a dt=0.25 s
     (512 muestras), formato ISOLA/KDEllipsPy.
  5. Concatena en archivos de una columna:
        DATA/real_disp_x  (componente N)
        DATA/real_disp_y  (componente E)
        DATA/real_disp_z  (componente Z)
     con npts*nestaciones filas, en el MISMO orden que el bloque de estaciones
     del input.ctl.

Tambien escribe ``input.ctl`` con:
  - fuente y mecanismo focal desde event.ctl (solo doble cupla -> MT flag 0),
  - banda 0.06-0.15 Hz, filtro causal,
  - parametros de inversion (plano de falla, elipse, NA) y modelo de
    velocidades 1D tomados de inversions/calama2020.

Uso:
    python codigos/real_disp.py
"""

import os

import numpy as np
import obspy

from event import load_event

# ==============================================================================
# CONFIGURACION
# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTL_FILE = os.path.join(BASE, "event.ctl")
VEL_DIR = os.path.join(BASE, "SAC", "VEL")
DATA_DIR = os.path.join(BASE, "DATA")
INPUT_CTL = os.path.join(BASE, "input.ctl")

# Estaciones de la inversion (subconjunto azimutal), ordenadas por distancia.
# AF01 quitada (2026-06-03): outlier del azimut este, lag sistematico +10 s
# por imprecision del modelo 1D hacia la cordillera (no es ruido). Ver memoria
# calama-2026-inversion-diag. Backup de su corrida en output_AF01_backup/.
SELECTED = ["A20F", "PB09", "A06F", "PB03", "PB15", "PB01", "PB19"]

# Ventana / muestreo (formato KDEllipsPy, igual que calama2020).
T1 = 0.0
T2 = 128.0
NPTS = 256
DT = 0.5                       # s  -> Nyquist 1 Hz (banda 0.06-0.15 Hz, margen de sobra)
# Banda de inversion y fase del filtro.
F1, F2 = 0.06, 0.15            # Hz
CAUSAL = True                  # True -> filtro causal (zerophase=False)

EVENT_NAME = "Event Calama 2026 INTRAPLATE"


# ==============================================================================
# 1. PROCESAR TRAZAS -> DESPLAZAMIENTO
# ==============================================================================
def procesar_estacion(station, ev):
    """Devuelve (n, e, z, lat, lon) procesados a NPTS muestras, o None."""
    comps = {}
    coords = None
    for comp in ("N", "E", "Z"):
        import glob
        matches = glob.glob(os.path.join(VEL_DIR, f"{station}.*{comp}.vel.sac"))
        if not matches:
            print(f"  [!] {station}: falta componente {comp}, se omite la estacion.")
            return None
        tr = obspy.read(matches[0])[0]

        # detrend -> integra (vel->disp) -> pasabanda causal
        tr.detrend("linear")
        tr.detrend("demean")
        tr.integrate(method="cumtrapz")
        tr.filter("bandpass", freqmin=F1, freqmax=F2, corners=4,
                  zerophase=not CAUSAL)

        # ventana + remuestreo exacto a NPTS @ DT
        t1 = ev["origin_time"] + T1
        t2 = ev["origin_time"] + T2
        tr.trim(starttime=t1, endtime=t2, pad=True, fill_value=0.0)
        tr.interpolate(sampling_rate=1.0 / DT, starttime=t1, npts=NPTS)

        comps[comp] = tr.data[:NPTS].astype(np.float64)
        coords = (tr.stats.sac.stla, tr.stats.sac.stlo)
    return comps["N"], comps["E"], comps["Z"], coords[0], coords[1]


def generar_datos(ev):
    os.makedirs(DATA_DIR, exist_ok=True)
    nx, ny, nz, estaciones = [], [], [], []
    print("Procesando estaciones (vel -> disp, causal %.2f-%.2f Hz):"
          % (F1, F2))
    for station in SELECTED:
        res = procesar_estacion(station, ev)
        if res is None:
            continue
        n, e, z, lat, lon = res
        nx.append(n)
        ny.append(e)
        nz.append(z)
        estaciones.append((station, lat, lon))
        print(f"  {station:6s}  lat {lat:8.3f}  lon {lon:8.3f}  OK")

    # Guardar archivos de una columna (npts*nestaciones).
    np.savetxt(os.path.join(DATA_DIR, "real_disp_x"), np.concatenate(nx), fmt="%.8e")
    np.savetxt(os.path.join(DATA_DIR, "real_disp_y"), np.concatenate(ny), fmt="%.8e")
    np.savetxt(os.path.join(DATA_DIR, "real_disp_z"), np.concatenate(nz), fmt="%.8e")
    print(f"\n  -> {len(estaciones)} estaciones x {NPTS} muestras "
          f"= {len(estaciones) * NPTS} filas por archivo en DATA/")
    return estaciones


# ==============================================================================
# 2. ESCRIBIR input.ctl
# ==============================================================================
# Modelo de velocidades 1D (inversions/calama2020): Thickness Vp Vs Rho Qp Qs
VELOCITY_MODEL = """ 0.000 5210.000 2990.000 2500.000  600.000  400.000
 2500.000 5370.000 3090.000 2500.000  600.000  400.000
 4500.000 5550.000 3190.000 2500.000  600.000  400.000
 6500.000 5720.000 3290.000 2600.000  600.000  400.000
 8500.000 5890.000 3390.000 2600.000  600.000  400.000
 10500.000 5980.000 3440.000 2600.000  600.000  400.000
 15000.000 6800.000 3750.000 2800.000  600.000  400.000
 20000.000 6810.000 3880.000 2800.000  600.000  400.000
 25000.000 6950.000 3940.000 3000.000  600.000  400.000
 30000.000 6980.000 4050.000 3000.000  600.000  400.000
 35000.000 7110.000 4110.000 3100.000  600.000  400.000
 40000.000 7410.000 4180.000 3300.000  600.000  400.000
 45000.000 7690.000 4300.000 3300.000  600.000  400.000
 50000.000 8050.000 4390.000 3300.000  600.000  400.000
 60000.000 8480.000 4730.000 3400.000  600.000  400.000
 70000.000 8480.000 4780.000 3400.000  600.000  400.000"""
N_LAYERS = 16


def escribir_input_ctl(ev, estaciones):
    zerophase_flag = 0 if CAUSAL else 1
    sta_lines = "\n".join(
        f"{lat:.4f} {lon:.4f} 0.0000 {name:<6s} 1 1 1"
        for name, lat, lon in estaciones
    )
    origin = ev["origin_time"]

    txt = f"""################################################################################
#                                                                              #
#             Control file for KINEMATIC INVERSION (Ellipse Patch)             #
#                                  v2025                                       #
#                                                                              #
################################################################################
#===============================================:===============================
# 1. Observed Data Parameters                   :     Values                   |
#===============================================:===============================
 Time window start (t1)                         :     {T1:.6f}
 Time window end (t2)                           :     {T2:.6f}
 Number of points (Npts)                        :     {NPTS}
 Delta / Time step                              :     {DT:.6f}
 Units                        (1:disp, 2:vel)   :     1

#===============================================:===============================
# 2. Source Position & Focal Mechanism          :     Values                   |
#===============================================:===============================
 Event Name                                     :     {EVENT_NAME}
 Origin Time (UTC)                              :     {origin.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]}Z
 Latitude                                       :     {ev['lat']:.6f}
 Longitude                                      :     {ev['lon']:.6f}
 Depth                                     (km) :     {ev['depth']:.6f}
 Strike                                  (deg.) :     {ev['strike']:.6f}
 Dip                                     (deg.) :     {ev['dip']:.6f}
 Rake                                    (deg.) :     {ev['rake']:.6f}

#===============================================:===============================
# 3. Fault Plane Parameters                     :     Values                   |
#===============================================:===============================
 Length along strike (Lx)                  (m)  :     50000.000000
 Length along dip (Ly)                     (m)  :     50000.000000
 Hypocenter position strike (Hx)           (m)  :     25000.000000
 Hypocenter position dip (Hy)              (m)  :     25000.000000
 Number of subfaults along strike (Nx)          :     25
 Number of subfaults along dip (Ny)             :     25

#===============================================:===============================
# 4. Ellipse Parameters & Frequency Band        :     Values                   |
#===============================================:===============================
 Number of ellipses                             :     1
 Initial slip                 (0:no, 1:yes)     :     0
 Slip shape            (0:cst, 1:gauss, 2:ell)  :     1
 Frequency 1 (Freq1)                      (Hz)  :     {F1:.6f}
 Frequency 2 (Freq2)                      (Hz)  :     {F2:.6f}
 Time shift (T0)                           (s)  :     3.000000
 Source type (axitra)                           :     4
 Zerophase filter (0:causal, 1:acausal)         :     {zerophase_flag}

#===============================================:===============================
# 5. Parameters to Invert (Min, Max, Flag)      :   Min      Max     Flag      |
#===============================================:===============================
 Number of parameters                           :     7
 Param 1 : Length of axis 1                (km) :     3.000000   20.000000   1
 Param 2 : Length of axis 2                (km) :     3.000000   20.000000   1
 Param 3 : Rotation angle                (x pi) :     0.000000   1.000000    1
 Param 4 : Position of the center np            :     0.000000   1.000000    1
 Param 5 : Position of the center tp    (x 2pi) :     0.000000   1.000000    1
 Param 6 : Maximum slip (Dmax)              (m) :     0.500000   4.000000    1
 Param 7 : Rupture velocity (Vr)         (km/s) :     1.000000   4.500000    1

#===============================================:===============================
# 6. Inversion Process Parameters               :     Values                   |
#===============================================:===============================
 Algorithm type                 (0:NA, 1:MC)    :     0
 NA: Number of iterations                       :     100
 NA: Sample size for first iteration (SS1)      :     2000
 NA: Sample size for other iterations           :     500
 NA: Cells to resample                          :     50
 Misfit time window     (s, 0.0=Full signal)    :     20.000
 MCMC: total steps                              :     500
 MCMC: burn-in                                  :     100
 MCMC: proposal scale                           :     0.050000
 MCMC: thinning                                 :     1
 MCMC: chains                                   :     1

#===============================================:===============================
# 7. Optional: Moment Tensor (Full MT)          :     Values                   |
#===============================================:===============================
 Moment Tensor Flag           (0:no, 1:yes)     :     0
 MT Scaling Mode (mt_strict/mt_factored)        :     mt_factored
 Mrr                                            :     0.000000
 Mtt                                            :     0.000000
 Mpp                                            :     0.000000
 Mrt                                            :     0.000000
 Mrp                                            :     0.000000
 Mtp                                            :     0.000000
 Exponent (iexp)                                :     18.0

#===============================================:===============================
# 8. Station Parameters (Lat, Lon, Elev, Name, R, T, Z) (0:no, 1:yes)          |
#===============================================:===============================
 Number of stations                             :     {len(estaciones)}
{sta_lines}

#===============================================:===============================
# 9. Velocity Model 1D                              : Thickness Vp Vs Rho Qp Qs |
#===============================================:===============================
 Number of layers                               :     {N_LAYERS}
{VELOCITY_MODEL}
"""
    with open(INPUT_CTL, "w") as f:
        f.write(txt)
    print(f"\ninput.ctl escrito en {INPUT_CTL}")
    print(f"  Banda {F1}-{F2} Hz | filtro {'CAUSAL' if CAUSAL else 'acausal'} "
          f"(Zerophase filter={zerophase_flag}) | solo doble cupla (MT flag 0)")


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    ev = load_event(CTL_FILE)
    print(f"Evento: {ev['origin_time']} | {ev['lat']}, {ev['lon']}, "
          f"{ev['depth']} km | strike/dip/rake = "
          f"{ev['strike']}/{ev['dip']}/{ev['rake']}\n")
    estaciones = generar_datos(ev)
    escribir_input_ctl(ev, estaciones)
    print("\nListo. Datos en DATA/ e input.ctl en la raiz del evento.")


if __name__ == "__main__":
    main()
