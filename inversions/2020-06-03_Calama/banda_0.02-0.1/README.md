# Banda 0.02–0.1 Hz — única corrida disponible

Corrida canónica para la banda 0.02–0.1 Hz del terremoto de Calama 2020-06-03.
Proviene de la antigua `output_prev_002-01` (renombrada). **Es la única corrida
que existe en esta banda** (las carpetas con nombre `band002-01` en realidad eran
0.06–0.15 Hz; ver `../banda_0.06-0.15/` y `../archive/`).

Carpeta autónoma:

```
banda_0.02-0.1/
├── input.ctl   ← config exacto de esta corrida (banda 0.02–0.1, plano 32 km)
├── codigos/    ← scripts (copia de ../codigos)
└── output/     ← resultados (inversion_result.joblib, figuras, logs)
└── Herreraetal/  Replica de la carpeta Herreraetal, misma configuracion
```

| Campo | Valor |
|---|---|
| Banda | 0.02–0.1 Hz, filtro **causal** |
| Plano | 32 × 32 km |
| Estaciones | 9 |
| Evaluaciones NA | 2480 |
| **Misfit** | **0.6194** |

Mejor modelo (`output/inversion_result.joblib`):

| a1 (km) | a2 (km) | θ (×π) | np | tp (×2π) | dmax (m) | vr (km/s) |
|---|---|---|---|---|---|---|
| 6.42 | 5.87 | 1.72 | 1.00 | 0.258 | 1.27 | 1.64 |

> Nota: plano de solo 32 km y `np=1` (centro contra el borde). Si se quiere
> comparar de igual a igual con la banda 0.06–0.15 (plano 50 km, 3000 evals),
> conviene re-correr esta banda con la misma geometría.
