# NP1 — plano nodal (Calama 2020-06-03, mecanismo ISOLA actualizado 2026-08-06)

Copia autónoma de la carpeta padre, con el mecanismo focal ISOLA actualizado
(Alex, MATLAB, `ISOLA_MT/invert/MTsol.txt`), que **reemplaza** la corrida
anterior (ver `../np1_old_mecanismo_v1/`, respaldada, NO borrada).

| Plane | Strike | Dip | Rake |
|---|---|---|---|
| NP1 (esta corrida) | 335° | 63° | -85° |
| NP2 (`../np2`) | 143° | 27° | -101° |
| NP1 anterior (`../np1_old_mecanismo_v1`) | 160° | 31° | -83° (mecanismo previo, superado) |

Hipocentro ISOLA: lat −23.2536, lon −68.4991, depth 115.557 km (antes
−23.253586/−68.479575/115.557 — cambio menor, ~2 km en longitud).

Se mantuvo la banda/ventana YA establecida en el `input.ctl` anterior
(128 pts, dt=1.0 s, 0.04–0.15 Hz, causal) — se corrigió `codigos/real_disp.py`
para que filtre los observados en esa misma banda (antes decía 0.1-0.3 Hz,
desajustado respecto al input.ctl). `DATA/` regenerada con esa banda; el resto
de `input.ctl` (fault plane, rangos de la elipse, parámetros NA, modelo de
velocidad, 9 estaciones) se dejó intacto, solo se tocó Lat/Lon/Depth/Strike/
Dip/Rake. `output/` queda vacía, a la espera de correr la inversión.

## Corrección 2026-10-05 (rama `calama-dynamic`)

`event.ctl` traía 143/27/−101 (los ángulos de NP2) y lo leen las figuras
(`mapa_azimutal.py`, `map_elipse_pygmt.py`, `plot_secciones.py`). Corregido a
335/63/−85. `input.ctl` ya estaba bien y es el que leyó la inversión
(confirmado en `output/run.log`: «strike 335 / dip 63 / rake -85»).

Resultado de la corrida de una semilla: misfit 0.1576 (eval 6048, iter 30),
a1 5.93 km, a2 6.34 km, theta 0.753π, Dmax 3.73 m, Vr 2.30 km/s.
