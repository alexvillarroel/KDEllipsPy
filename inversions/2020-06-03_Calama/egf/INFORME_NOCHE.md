# Robustez inversion EGF — Calama 2020 (0.1-0.3 Hz)

| run | misfit | a1 | a2 | theta | np | tp | dmax | vr |
|---|---|---|---|---|---|---|---|---|
| principal (seed11, 7 est) | 0.882 | 3.46 | 6.07 | 0.85 | 0.82 | 0.04 | 2.92 | 1.59 |
| seed22 | 0.870 | 3.00 | 5.08 | 0.49 | 1.00 | 0.49 | 4.00 | 2.07 |
| seed33 | 0.869 | 3.03 | 4.99 | 0.50 | 1.00 | 0.48 | 4.00 | 1.98 |
| sin AF01 | 0.955 | 4.15 | 3.45 | 0.44 | 0.99 | 0.51 | 3.74 | 3.03 |
| sin PB19 | 0.874 | 3.11 | 5.95 | 0.90 | 0.99 | 0.04 | 3.34 | 1.63 |
| EGF M4.6 (3 est) | 0.941 | 9.87 | 5.94 | 0.10 | 0.80 | 0.18 | 0.76 | 1.79 |

## Dispersion entre runs M4.9 (semillas + jackknife)

| param | spread | % rango | veredicto |
|---|---|---|---|
| a1 | 1.15 | 10% | ROBUSTO |
| a2 | 2.62 | 24% | MEDIO |
| theta | 0.45 | 45% | NO |
| np | 0.18 | 18% | MEDIO |
| tp | 0.47 | 47% | NO |
| dmax | 1.08 | 31% | MEDIO |
| vr | 1.44 | 48% | NO |

EGF cruzada (M4.6, 3 estaciones) como consistencia independiente — no entra al spread.

## Interpretacion (noche 2026-06-12)

### Lo que la EGF SI logro (primera vez en todo el estudio 2020)

1. **El paisaje discrimina**: modelos aleatorios dan misfit mediana 1.12
   (p90 2.1) vs best 0.87-0.88; solo 0.2% de los aleatorios se acerca al
   optimo. Con AXITRA a 0.06-0.15 Hz esto NUNCA paso (marginales planas).
2. **TAMANO ROBUSTO**: a1 = 3.0-4.2 km en TODOS los runs M4.9 (10% del rango);
   a2 = 3.5-6.1 km. Aspereza COMPACTA (~3x6 km de semiejes, area ~50-70 km2)
   — mucho menor que el plano de 50 km. Consistente con alto stress drop
   tipico de eventos intraslab profundos. Incluso la EGF cruzada M4.6
   (3 estaciones) da a2 = 5.9 km.
3. **Vr lenta**: 1.6-2.1 km/s en 4/5 runs (~0.35-0.45 Vs).
4. **Directividad real en los datos**: lags S azimutales (norte +5/+6 s,
   oeste +3.4 s, sur -1 s).

### Lo que queda bimodal (no resuelto)

Dos cuencas con Dmisfit = 0.012 (marginal vs piso EGF ~0.87):
- **Cuenca A** (seed11, sinPB19; tambien np/tp de la EGF M4.6): theta ~0.85-0.90 pi,
  np ~0.8, tp ~0.04, vr ~1.6. Solucion interior (sin rieles). Coincide
  notablemente con el best AXITRA causal 0.02-0.1 (theta 0.91, np 0.80,
  tp 0.14, vr 1.61) — dos Greens INDEPENDIENTES.
- **Cuenca B** (seed22/33; misfit 0.870, levemente mejor): theta ~0.5,
  tp ~0.49, pero con np = 1.00 y dmax = 4.00 RIELANDO (sospechosa; pide
  salirse de los limites).

Verificacion cruzada de cuencas dentro de cada ensamble: cada semilla solo
muestreo bien su cuenca (A: 0.882 vs B mal muestreada 0.923 en seed11;
B: 0.870 vs A mal muestreada 0.994 en seed22) -> para zanjar hace falta un
run con limites ampliados (dmax hasta ~8 m) y/o NA mas exploratorio
(mas semillas / ns mayor), o muestreo dirigido de ambas cuencas.

### Caveats

- **Piso de misfit ~0.87** fijado por la calidad de la EGF (CC_S ~0.7
  -> el misfit absoluto NO es comparable con los runs AXITRA).
- **M0 subestimado ~2x** (Mw 6.6 vs 6.8 catalogo): el escalado optimo
  subpredice amplitud cuando hay desfase residual (clasico en EGF); ademas
  M0_egf viene de Mw de catalogo (incierto facilmente 1.5-2x). Trade-off
  directo con dmax (presiona el riel 4.0 m); NO afecta la geometria.
- Appraisal NA no confiable aqui (posterior picuda en 7D + caminata Gibbs
  no mezcla; IQR=50% espurio a toda T). Diagnostico valido: ensamble
  directo (misfit_vs_param.png) + robustez empirica entre runs.

### Recomendacion para manana

1. Re-correr con dmax hasta 8 m (desriela B y resuelve el trade-off M0).
2. 3-4 semillas mas para mapear la frecuencia de cada cuenca.
3. Si A persiste como atractor sin rieles y B sigue rielando -> reportar A
   (que ademas coincide con AXITRA independiente).

## Test dmax<=8 m (3 semillas, desempate de cuencas)

| run | misfit | a1 | a2 | theta | np | tp | dmax | vr |
|---|---|---|---|---|---|---|---|---|
| dmax8_s11 | 0.8804 | 4.00 | 5.92 | 0.44 | 1.00 | 0.50 | 2.88 | 3.00 |
| dmax8_s22 | **0.8665** | 2.58 | 4.44 | 0.48 | 0.96 | 0.50 | 6.31 | 1.79 |
| dmax8_s44 | 0.8901 | 7.75 | 4.44 | 0.05 | 0.89 | 0.21 | 1.86 | 2.94 |

- Con dmax libre, la **cuenca B gana** (2/3 runs; mejor global 0.8665 con
  dmax=6.31 -> M0 ya del orden del catalogo). El riel de dmax=4 si era
  vinculante para B.
- **np >= 0.8 en los 9 runs de la noche** -> hipocentro en el BORDE de la
  aspereza: rasgo robusto adicional (nucleacion al borde, comun en rupturas
  reales).
- s44 encuentra una tercera variante (theta 0.05 con a1<->a2 intercambiados:
  familia parcialmente degenerada de la elipse).
- **Vr deja de ser estable** con limites amplios (1.8-3.0) -> retirarla de
  los robustos.

## VEREDICTO FINAL DE LA NOCHE

ROBUSTO (9/9 runs):
  1. Aspereza COMPACTA: a2 = 4.4-6.1 km, a1 mayormente 2.6-4.2 km
     (un outlier 7.8 por degeneracion a1<->a2). Area ~40-80 km2 en un plano
     de 2500 km2 -> alto stress drop, consistente con intraslab profundo.
  2. Hipocentro en el borde de la aspereza (np = 0.8-1.0 SIEMPRE).
  3. El paisaje EGF discrimina (aleatorios 1.12 vs best 0.87) — la banda
     0.1-0.3 con Green empirica SI contiene informacion de fuente finita.

NO RESUELTO al piso de calidad EGF (~0.87, CC~0.7):
  - Orientacion theta y posicion angular tp (familia multi-modal,
    Dmisfit < 0.02 entre configuraciones).
  - Vr (1.8-3.0).
  - M0 dentro de factor ~2 (incertidumbre M0_egf + desfase residual).

Para mejorar el piso: EGF de mas calidad (apilar 2-3 replicas validas como
EGF compuesta; deconvolucion de la STF de la EGF; o banda 0.15-0.4 si el SNR
del M4.9 lo permite).
