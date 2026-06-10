# Códigos — Evento Copiapó 2025-06-06

Procesamiento de los registros de aceleración del CSN para el sismo de Copiapó
(2025-06-06 17:15:06 UTC, lat −26.639, lon −70.404, prof 75 km, intraplaca,
strike/dip/rake = 71/46/−110). Misma estructura que el evento de Calama 2026.

## Estructura (módulos seccionados)

| Archivo | Qué hace |
|---|---|
| `event.py` | Lee `event.ctl` → dict con lat/lon/prof/strike/dip/rake/tiempo de origen. |
| `io_sismica.py` | Lee los `.txt` del CSN (`registros_sismicos/<EST>/*.txt`) → `obspy.Stream`, parseando la cabecera comentada (estación, canal, Fs, npts, lat/lon, unidades). |
| `integracion.py` | Port a Python de los MATLAB `integraf.m`, `time.m` y `pssa3.m`. Integra aceleración → velocidad → desplazamiento con filtrado Butterworth. |
| `txt_to_sac.py` | **Script principal.** Orquesta todo y escribe los SAC (acc/vel/desp) con geometría completa. |
| `real_disp.py` | Prepara la inversión KDEllipsPy: integra vel→desp, filtra **causal 0.06–0.15 Hz**, genera `DATA/real_disp_{x,y,z}` (8 estaciones) y escribe `input.ctl` (doble cupla, modelo 1D y parámetros de Copiapó). |
| `run_inversion.py` | CLI que corre la búsqueda NA con KDEllipsPy y guarda resultados/figuras en `output/`. |
| `map.py` | Mapa PyGMT del evento (estaciones + formas de onda + mecanismo focal + zoom con círculo de 200 km). |
| `plot_check.py` | Figura acc/vel/desp (grilla 3×3 N/E/Z) por estación a <200 km. Salida en `SAC/checks/`. |
| `plot_secciones.py` | 3 PDFs (acc/vel/desp), filas = estaciones <200 km, columnas = N/E/Z, eje X alineado al origen. |

## Flujo

```
registros_sismicos/<EST>/*.txt
        │  txt_to_sac.py  (io_sismica + integracion.pssa3)
        ▼
   SAC/ACC/   SAC/VEL/   SAC/DISP/
        │  real_disp.py  (vel→desp, causal 0.06–0.15 Hz)
        ▼
   DATA/real_disp_{x,y,z}  +  input.ctl
        │  run_inversion.py  (Neighbourhood Algorithm)
        ▼
   output/  (inversion_result.joblib + figuras)
```

## Uso

```bash
cd codigos
python txt_to_sac.py        # 1) .txt -> SAC (acc, vel, desp)
python real_disp.py         # 2) prepara DATA/ + input.ctl
python plot_check.py        # (opcional) figuras de revisión por estación
python plot_secciones.py    # (opcional) secciones acc/vel/desp
python run_inversion.py     # 3) corre la inversión NA
python map.py               # (opcional) mapa — requiere pygmt
```

## Estaciones de la inversión (cobertura azimutal)

8 estaciones elegidas para el mejor spread azimutal posible (el cuadrante O–NO
es océano y queda vacío):

| Estación | Azimut | Epic (km) | Sector |
|---|---|---|---|
| AC01 | 340° | 58 | NO |
| A24F | 2.5° | 336 | N |
| A28F | 27.7° | 191 | NNE |
| A15C | 52.6° | 45 | NE |
| A17C | 104° | 51 | E |
| A06C | 176° | 79 | S |
| A18C | 206° | 107 | SSO |
| A30C | 221° | 63 | SO |

## Notas

- Los `.txt` ya vienen en m/s² (no se remueve respuesta instrumental).
- Canales `HN*` a 200 Hz, `HL*`/otros a 100 Hz (acelerómetros).
- El "Tiempo de Origen" en la cabecera de cada `.txt` es el tiempo de la primera
  muestra de la traza; el origen real del sismo está en `event.ctl`.
- Modelo de velocidades 1D: 20 capas (mismo del proyecto antiguo `copiapo-06-06-2025`).
- `pygmt` se corre con `/home/alex/.conda/envs/kdellipspy/bin/python`.
