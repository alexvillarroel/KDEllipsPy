# Inversion EGF — Calama 2020 (banda 0.15-0.4 Hz)

EGF: 20200616T175441 M4.89 (M0_egf=2.43e+16 Nm, de Mw catalogo).
Estaciones (8): PB01, PB02, PB03, PB05, PB06, PB09, PB19, AF01

**Mejor misfit: 0.7618**  (M0=1.32e+19 Nm, Mw=6.71)

| param | valor |
|---|---|
| a1 (km) | 5.832 |
| a2 (km) | 4.239 |
| theta (x pi) | 0.080 |
| np (frac) | 0.231 |
| tp (x 2pi) | 0.112 |
| dmax (m) | 3.416 |
| vr (km/s) | 2.396 |
Lags estaticos (max +-4.0 s) del mejor modelo; ventanas P excluidas: PB02, PB09

| sta | lag P (s) | lag S (s) |
|---|---|---|
| PB01 | -0.50 | +1.75 |
| PB02 | -- | -2.75 |
| PB03 | +2.00 | +0.50 |
| PB05 | +1.25 | +1.50 |
| PB06 | +2.50 | +1.75 |
| PB09 | -- | +1.00 |
| PB19 | +1.50 | +0.75 |
| AF01 | +1.25 | +1.25 |


Notas: EGF alineada solo por dt de viaje S (taup); los lags de directividad quedan como senal. M0_egf incierto -> trade-off directo con dmax (geometria no afectada). P (R,Z) con CC moderada; la senal dominante es S (T).

![ajuste](ajuste_ondas.png)
