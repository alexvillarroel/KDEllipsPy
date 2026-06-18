# Banda 0.06–0.15 Hz — corrida de producción (pipeline nuevo)

Corrida canónica para la banda 0.06–0.15 Hz del terremoto de Calama 2020-06-03.
Proviene de la antigua `output_vel_t34_006-015` (renombrada). Las demás corridas
de esta banda (disp, plano64, pruebas, etc.) están en `../archive/`.

Carpeta autónoma:

```
banda_0.06-0.15/
├── input.ctl   ← config exacto de esta corrida (banda 0.06–0.15, plano 50 km)
├── codigos/    ← scripts (copia de ../codigos)
└── output/     ← resultados (inversion_result.joblib, figuras, logs)
```

| Campo | Valor |
|---|---|
| Banda | 0.06–0.15 Hz, filtro **causal** |
| Datos | formas de onda en **velocidad**, ventana misfit 20 s (P→R,Z ; S→T) |
| Plano | 50 × 50 km |
| Estaciones | 9 |
| Evaluaciones NA | 3000 |
| **Misfit** | **0.2063** |

Mejor modelo (`output/inversion_result.joblib`):

| a1 (km) | a2 (km) | θ (×π) | np | tp (×2π) | dmax (m) | vr (km/s) |
|---|---|---|---|---|---|---|
| 5.24 | 10.3 | 0.111 | 0.997 | 0.554 | 2.98 | 4.27 |

Es la corrida que tu `resumen_trabajo.tex` reporta como producción (misfit ~0.204).
