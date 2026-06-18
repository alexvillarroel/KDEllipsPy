# Inversion EGF — Calama 2020 (banda 0.15-0.4 Hz)

EGF: 20200616T175441 M4.89 (M0_egf=2.43e+16 Nm, de Mw catalogo).
Estaciones (7): PB01, PB03, PB05, PB06, PB09, PB19, AF01

**Mejor misfit: 0.6013**  (M0=9.48e+19 Nm, Mw=7.28)

| param | valor |
|---|---|
| a1 (km) | 6.601 |
| a2 (km) | 12.039 |
| theta (x pi) | 0.711 |
| np (frac) | 0.759 |
| tp (x 2pi) | 0.444 |
| dmax (m) | 7.687 |
| vr (km/s) | 1.866 |
Lags estaticos (max +-4.0 s) del mejor modelo; ventanas P excluidas: PB02, PB09

| sta | lag P (s) | lag S (s) |
|---|---|---|
| PB01 | -2.50 | -1.00 |
| PB03 | -0.25 | +3.25 |
| PB05 | -1.75 | -2.25 |
| PB06 | +2.00 | +3.50 |
| PB09 | -- | -2.25 |
| PB19 | +2.25 | -0.75 |
| AF01 | -4.00 | +2.00 |


Notas: EGF alineada solo por dt de viaje S (taup); los lags de directividad quedan como senal. M0_egf incierto -> trade-off directo con dmax (geometria no afectada). P (R,Z) con CC moderada; la senal dominante es S (T).

![ajuste](ajuste_ondas.png)

![appraisal](appraisal_corner_egf.png)
