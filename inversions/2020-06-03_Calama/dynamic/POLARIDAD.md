# Polaridad del dinámico en Calama 2020 — resultado

Un forward de fd3d por plano, los dos signos del vector de momento, cada signo
en su mejor dt0 (barrido -3..+10 s).

| caso | rake ctl | signo | misfit | dt0 | rake efectivo | cc media |
|---|---|---|---|---|---|---|
| np1 (335/63) | -85  | +1 (tal cual) | 0.688 | -2.0 s | +90.0 | +0.667 |
| np1 (335/63) | -85  | **-1 (flip)** | **0.403** | **+3.0 s** | **-90.0** | **+0.870** (9/9) |
| np2 (143/27) | -101 | +1 (tal cual) | 0.721 | -2.0 s | +90.0 | +0.589 |
| np2 (143/27) | -101 | **-1 (flip)** | **0.458** | **+3.0 s** | **-90.0** | **+0.725** (8/9) |

- Los dos planos, con dip 63 y 27, piden el mismo flip y el mismo dt0 (+3.0 s),
  que coincide con el `Time shift (T0) = 3.0` del input.ctl del cinemático. El
  camino de la base dinámica usa t0 = 0 y no lo lee: los dos caminos llegan al
  mismo desfase por separado.
- fd3d solo produce deslizamiento puro en el buzamiento (rake +-90): el prestress
  de la aspereza no le impone componente de rumbo. Por eso el rake efectivo da
  -90 en los dos planos, a 5 grados del objetivo en np1 y a 11 en np2. Al
  comparar misfits entre planos, np2 arrastra esa penalización de geometría.
- El signo del primer arribo P no se usa: los dos sintéticos son negativos
  exactos uno del otro, los conteos suman siempre 9 y no discrimina.
