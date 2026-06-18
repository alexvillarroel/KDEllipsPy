# Inversion EGF — Calama 2020 (banda 0.15-0.4 Hz)

EGF: 20200616T175441 M4.89 (M0_egf=2.43e+16 Nm, de Mw catalogo).
Estaciones (8): PB01, PB02, PB03, PB05, PB06, PB09, PB19, AF01

**Mejor misfit: 0.6208**  (M0=5.50e+19 Nm, Mw=7.13)

| param | valor |
|---|---|
| a1 (km) | 9.781 |
| a2 (km) | 8.380 |
| theta (x pi) | 0.387 |
| np (frac) | 0.253 |
| tp (x 2pi) | 0.976 |
| dmax (m) | 4.489 |
| vr (km/s) | 1.549 |
Lags estaticos (max +-4.0 s) del mejor modelo; ventanas P excluidas: PB02, PB09

| sta | lag P (s) | lag S (s) |
|---|---|---|
| PB01 | -2.50 | +4.00 |
| PB02 | -- | +2.00 |
| PB03 | +0.25 | -1.25 |
| PB05 | -1.00 | +4.00 |
| PB06 | +3.00 | +4.00 |
| PB09 | -- | -1.25 |
| PB19 | -1.00 | -2.00 |
| AF01 | -2.00 | +4.00 |


Notas: EGF alineada solo por dt de viaje S (taup); los lags de directividad quedan como senal. M0_egf incierto -> trade-off directo con dmax (geometria no afectada). P (R,Z) con CC moderada; la senal dominante es S (T).

![ajuste](ajuste_ondas.png)

![appraisal](appraisal_corner_egf.png)
