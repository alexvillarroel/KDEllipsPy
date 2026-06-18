# Códigos — Evento Calama 2026-05-25

Procesamiento de los registros de aceleración del CSN para el sismo de Calama
(2026-05-25 21:52:18 UTC, lat −22.38, lon −68.76, prof 114 km, Mw≈6.9).

## Estructura (módulos seccionados)

| Archivo | Qué hace |
|---|---|
| `event.py` | Lee `event.ctl` → dict con lat/lon/prof/strike/rake/dip/tiempo de origen. |
| `io_sismica.py` | Lee los `.txt` del CSN (`registros_sismicos/<EST>/*.txt`) → `obspy.Stream`, parseando la cabecera comentada (estación, canal, Fs, npts, lat/lon, unidades). |
| `integracion.py` | Port a Python de los MATLAB `integraf.m`, `time.m` y `pssa3.m`. Integra aceleración → velocidad → desplazamiento con filtrado Butterworth (filtfilt). |
| `txt_to_sac.py` | **Script principal.** Orquesta todo y escribe los SAC. |
| `map.py` | Mapa PyGMT del evento (estaciones + formas de onda + mecanismo focal + panel de zoom con círculo de 200 km). |
| `plot_check.py` | Figura acc/vel/desp (grilla 3×3 N/E/Z) por estación a <200 km del epicentro, con línea vertical en el tiempo de origen. Salida en `SAC/checks/`. |
| `plot_secciones.py` | 3 PDFs (uno por modo: acc/vel/desp). Filas = estaciones <200 km (por distancia), columnas = N/E/Z, eje X alineado al origen. Salida `SAC/checks/seccion_*.pdf`. |
| `real_disp.py` | Prepara la inversión KDEllipsPy: integra vel→desp, filtra causal 0.06–0.15 Hz, genera `DATA/real_disp_{x,y,z}` (8 estaciones) y escribe `input.ctl` (doble cupla, modelo y parámetros de calama2020). |

## Flujo

```
registros_sismicos/<EST>/*.txt
        │  io_sismica.read_all()
        ▼
   obspy.Stream (aceleración, m/s²)
        │  integracion.pssa3()   (acc → vel → desp)
        ▼
   SAC/ACC/   SAC/VEL/   SAC/DISP/
```

## Uso

```bash
# 1) Convertir los .txt a SAC (acc, vel, desp). Se ejecuta desde codigos/
cd codigos
python txt_to_sac.py

# 2) (opcional) Mapa — requiere pygmt (y rockhound para los contornos de slab)
python map.py
```

## Salida SAC

Se generan **9 archivos por estación** (3 componentes × 3 tipos), nombrados
`{ESTACION}.{CANAL}.{tipo}.sac` en:

- `SAC/ACC/`  → aceleración corregida (m/s²), `idep=8`
- `SAC/VEL/`  → velocidad (m/s), `idep=7`
- `SAC/DISP/` → desplazamiento (m), `idep=6`

Cada SAC lleva en su cabecera: `kstnm`, `kcmpnm`, `stla`, `stlo`,
`evla`, `evlo`, `evdp`, `dist` (**distancia hipocentral en km**),
`gcarc`/`az`/`baz` (distancia epicentral en grados, azimut, back-azimut) y
`user0` (distancia epicentral en km).

## Parámetros de filtrado

En `txt_to_sac.py`:

```python
FILT_F1 = 0.02   # Hz, esquina baja (estabiliza la integración)
FILT_F2 = 40.0   # Hz, esquina alta (se limita a 0.45·Fs por traza)
FILT_ORDER = 4   # orden Butterworth
ZEROPHASE = True # filtfilt (acausal, fase cero)
```

Esto entrega productos de **banda ancha**. Para la inversión de momento tensor
(KDEllipsPy) normalmente se vuelve a filtrar a la banda fina (p.ej. 0.04–0.2 Hz)
sobre el desplazamiento ya integrado.

## Notas

- Los `.txt` ya vienen en unidades físicas (m/s²): **no** se remueve respuesta
  instrumental.
- Canales `HN*` muestreados a 200 Hz, `HL*` a 100 Hz (ambos acelerómetros).
- El "Tiempo de Origen" de la cabecera de cada `.txt` es el tiempo de la primera
  muestra de la traza, no el origen del sismo; se usa como `starttime`.
- `pygmt`/`rockhound` **no** están instalados en el entorno actual
  (`base`/`kdellipspy`); instálalos para correr `map.py`.
```
