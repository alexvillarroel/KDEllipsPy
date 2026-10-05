# NP2 — plano nodal (Calama 2020-06-03, mecanismo ISOLA actualizado 2026-08-06)

Copia autónoma de la carpeta padre, con el mecanismo focal ISOLA actualizado
(Alex, MATLAB, `ISOLA_MT/invert/MTsol.txt`), que **reemplaza** la corrida
anterior (ver `../np1_old_mecanismo_v1/`, respaldada, NO borrada).

| Plane | Strike | Dip | Rake |
|---|---|---|---|
| NP1 (`../np1`) | 335° | 63° | -85° |
| NP2 (esta corrida) | 143° | 27° | -101° |
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

`input.ctl` había quedado sobrescrito con los ángulos de NP1 (335/63/−85),
pese a que la inversión sí corrió con NP2 — `output/run.log` dice
«strike 143 / dip 27 / rake -101». Restaurado a 143/27/−101 para que un
re-run reproduzca lo que está en `output/`. `np1/` y `np2/` están sin trackear
en git (solo se versionan `input.ctl`, `event.ctl` y `README.md`), así que la
geometría correcta solo existía en el log.

Resultado de la corrida de una semilla: misfit 0.1812 (eval 5942, iter 29),
a1 6.03 km, a2 5.99 km, theta 0.201π, Dmax 3.89 m, Vr 1.61 km/s.
