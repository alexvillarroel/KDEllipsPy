"""
real_disp.py
============
Prepara los datos para la inversion cinematica con KDEllipsPy (estilo del
``real_disp.py`` de copiapo) y escribe el ``input.ctl``.

Para las estaciones seleccionadas (subconjunto azimutal):
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
  - fuente y mecanismo focal desde event.ctl + MT completo de ISOLA (MT flag 1),
  - banda 0.06-0.15 Hz, filtro causal,
  - parametros de inversion (plano de falla, elipse, NA) y modelo de
    velocidades 1D del evento de Copiapo 2025.

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

# Estaciones de la inversion (subconjunto azimutal recomendado), ordenadas por
# azimut para cubrir lo mejor posible la red (el cuadrante O-NO es oceano).
#   AC01 340 | A24F 2.5 | A28F 27.7 | A15C 52.6 | A17C 104 | A06C 176
#   A18C 206 | A30C 221 | A16C 61 (88.8 km, agregada como 9a; refuerza NE)
SELECTED = ["AC01", "A24F", "A28F", "A15C", "A17C", "A06C", "A18C", "A30C", "A16C"]

# Ventana / muestreo (formato KDEllipsPy).
T1 = 0.0
T2 = 128.0
NPTS = 512
DT = 0.25                      # s  -> 4 Hz
# Banda de inversion y fase del filtro.
F1, F2 = 0.06, 0.15            # Hz
CAUSAL = True                  # True -> filtro causal (zerophase=False)

EVENT_NAME = "Event Copiapo 2025 INTRAPLATE"


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
# Modelo de velocidades 1D del evento de Copiapo 2025 (Thickness Vp Vs Rho Qp Qs;
# la primera columna es profundidad al techo de cada capa, en m).
VELOCITY_MODEL = """        0.0        4700.00      2750.00     2562.0       300.0     300.0
     5000.0        5730.00      3300.00     2692.0       300.0     300.0
    10000.0        6320.00      3620.00     2817.0       300.0     300.0
    15000.0        6680.00      3830.00     2946.0       300.0     300.0
    20000.0        6900.00      3950.00     3026.0       300.0     300.0
    25000.0        7020.00      4020.00     3069.0       300.0     300.0
    30000.0        7100.00      4060.00     3098.0       300.0     300.0
    35000.0        7250.00      4140.00     3152.0       300.0     300.0
    40000.0        7480.00      4260.00     3235.0       300.0     300.0
    45000.0        7790.00      4430.00     3346.0       300.0     300.0
    50000.0        8070.00      4570.00     3447.0      1000.0    1000.0
    55000.0        8230.00      4640.00     3505.0      1000.0    1000.0
    60000.0        8360.00      4680.00     3551.0      1000.0    1000.0
    65000.0        8480.00      4730.00     3595.0      1000.0    1000.0
    70000.0        8550.00      4770.00     3620.0      1000.0    1000.0
    75000.0        8570.00      4780.00     3627.0      1000.0    1000.0
    80000.0        8550.00      4760.00     3620.0      1000.0    1000.0
    85000.0        8500.00      4740.00     3602.0      1000.0    1000.0
    90000.0        8440.00      4720.00     3580.0      1000.0    1000.0
    95000.0        8380.00      4690.00     3559.0      1000.0    1000.0"""
N_LAYERS = 20


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
#                                  v2025                                        #
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
 Length along strike (Lx)                  (m)  :     32000.000000
 Length along dip (Ly)                     (m)  :     32000.000000
 Hypocenter position strike (Hx)           (m)  :     16000.000000
 Hypocenter position dip (Hy)              (m)  :     16000.000000
 Number of subfaults along strike (Nx)          :     30
 Number of subfaults along dip (Ny)             :     30

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
 Param 1 : Length of axis 1                (km) :     3.000000   10.000000   1
 Param 2 : Length of axis 2                (km) :     3.000000   10.000000   1
 Param 3 : Rotation angle                (x pi) :     0.000000   2.000000    1
 Param 4 : Position of the center np            :     0.000000   1.000000    1
 Param 5 : Position of the center tp    (x 2pi) :     0.000000   1.000000    1
 Param 6 : Maximum slip (Dmax)              (m) :     0.300000   4.000000    1
 Param 7 : Rupture velocity (Vr)         (km/s) :     1.800000   3.500000    1

#===============================================:===============================
# 6. Inversion Process Parameters               :     Values                   |
#===============================================:===============================
 Algorithm type                 (0:NA, 1:MC)    :     0
 NA: Number of iterations                       :     8
 NA: Sample size for first iteration (SS1)      :     600
 NA: Sample size for other iterations           :     60
 NA: Cells to resample                          :     8
 Misfit time window     (s, 0.0=Full signal)    :     0.000
 MCMC: total steps                              :     500
 MCMC: burn-in                                  :     100
 MCMC: proposal scale                           :     0.050000
 MCMC: thinning                                 :     1
 MCMC: chains                                   :     1

#===============================================:===============================
# 7. Optional: Moment Tensor (Full MT)          :     Values                   |
#===============================================:===============================
 Moment Tensor Flag           (0:no, 1:yes)     :     1
 MT Scaling Mode (mt_strict/mt_factored)        :     mt_factored
 Mrr                                            :     -2.486000e+00
 Mtt                                            :     7.932000e+00
 Mpp                                            :     8.300000e-02
 Mrt                                            :     3.800000e-02
 Mrp                                            :     -7.390000e-01
 Mtp                                            :     4.730000e-01
 Exponent (iexp)                                :     18.0

#===============================================:===============================
# 8. Station Parameters (Lat, Lon, Elev, Name, N, E, Z) (0:no, 1:yes)          |
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
          f"(Zerophase filter={zerophase_flag}) | MT completo activado (flag 1, ISOLA)")


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
