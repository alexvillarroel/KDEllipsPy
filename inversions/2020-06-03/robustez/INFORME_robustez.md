# Informe de robustez — inversion cinematica 2020

Corridas: 15 (repetibilidad 5, jackknife 8, bandas 2).

## Dispersion por parametro (todas las corridas)

| Param | mediana | min | max | IQR | rango | ¿robusto? |
|---|---|---|---|---|---|---|
| a1 | 6.410 | 5.022 | 12.249 | 4.264 | [4,13] | NO (47%) |
| a2 | 7.483 | 4.117 | 12.661 | 5.571 | [4,13] | NO (62%) |
| theta | 0.498 | 0.081 | 0.901 | 0.372 | [0,1] | NO (37%) |
| np | 0.999 | 0.109 | 1.000 | 0.010 | [0,1] | si (1%) |
| tp | 0.454 | 0.048 | 0.961 | 0.259 | [0,1] | NO (26%) |
| dmax | 4.573 | 1.012 | 5.914 | 1.748 | [1,6] | NO (35%) |
| vr | 0.979 | 0.576 | 1.275 | 0.169 | [0.5,3.5] | si (6%) |

## Misfit por corrida

| corrida | banda | misfit |
|---|---|---|
| jack_no_A18F | 0.06-0.15 | 0.3685 |
| rep_seed44 | 0.06-0.15 | 0.3764 |
| jack_no_A22F | 0.06-0.15 | 0.4155 |
| rep_seed22 | 0.06-0.15 | 0.4155 |
| jack_no_A09F | 0.06-0.15 | 0.4162 |
| jack_no_A23F | 0.06-0.15 | 0.4162 |
| jack_no_PB19 | 0.06-0.15 | 0.4233 |
| rep_seed33 | 0.06-0.15 | 0.4293 |
| jack_no_PB06 | 0.06-0.15 | 0.4355 |
| rep_seed11 | 0.06-0.15 | 0.4389 |
| rep_seed55 | 0.06-0.15 | 0.4455 |
| jack_no_A06F | 0.06-0.15 | 0.4557 |
| band_0.02-0.10 | 0.02-0.1 | 0.4754 |
| jack_no_A05F | 0.06-0.15 | 0.5214 |
| band_0.04-0.12 | 0.04-0.12 | 0.6763 |

**Lectura:** un parametro es robusto si su IQR es <10% del rango y coincide entre semillas, jackknife y bandas. Parametros con IQR alto o que dependen de quitar una estacion (jackknife) NO estan bien restringidos.
