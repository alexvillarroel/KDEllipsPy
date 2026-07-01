# Forward legacy vs Python — mismo mejor modelo (A05C)

**Pregunta:** ¿el desajuste de A05C viene de la STF, el time shift o la función de Green?
**Método:** forwardeé el **mismo mejor modelo** de la corrida Python (Herreraetal,
8 est, misfit 0.178) a través del **código Fortran legacy** (`kine_direct`), reusando
las Green's `axi.res` existentes (geometría idéntica: plano 32 km, 30×30, hipocentro
16/16, mismas 8 estaciones en el mismo orden). Comparé sintético vs sintético.

Modelo: `a1=5.21, a2=10.91, θ=1.921, np=0.736, tp=0.667, dmax=2.5, vr=1.73`

## Resultado (correlación legacy vs Python, por estación)

| Estación | dist | corr0 (N/E/Z) aprox | veredicto |
|---|---|---|---|
| PB05, PB06, PB08, PB19 | cercanas | 0.6–0.85 | **coinciden** (offset ~+1 s) |
| T15A | media | 0.3–0.83 | mayormente coincide |
| AC01, AC02 | lejanas sur | 0.26–0.61, lags −12…−20 s | divergen |
| **A05C** | **más lejana (~450 km)** | **0.04–0.37** | **divergen por completo** |

Ver `comparacion_legacy_vs_python.png`: PB06 se superpone; A05C no.

## Conclusión

- **No es la STF, ni el time shift, ni el momento/unidades.** Esos afectarían a
  TODAS las estaciones por igual, y las cercanas coinciden bien.
- **Sí es la función de Green / forward en el campo lejano.** Los dos códigos
  producen la misma respuesta a distancia normal pero **divergen en las 3
  estaciones lejanas del sur (A05C, AC01, AC02)**. En A05C el Python coloca la
  energía en ~55–70 s y el legacy en ~100–128 s.

## Causa probable (a confirmar)

Diferencia que escala con la distancia → candidatos:
1. **Muestreo / nfreq:** Python dt=1.0 → nfreq=64; legacy dt=0.25 → nfreq=128.
   A distancia, la Green tiene más fases/cola; con nfreq=64 puede quedar
   sub-resuelta en tiempo. (Aunque la banda 0.02–0.1 está cubierta espectralmente.)
2. **Convención de origen / reducción temporal por distancia** (el desfase grande
   solo en lejanas apunta a esto).
3. Integración de número de onda (xl periodicidad ≈ 1.2e6 vs 1e6; ikmax igual).

## CONFIRMADO (test 2): era la discretización de la Green

Re-corrí el forward **Python con los parámetros del legacy** (nfreq=128, tl=128,
xl=1e6) — ver `match_legacy_params.png`. Resultado:

| | A05C corr (antes→ahora) | todas las estaciones |
|---|---|---|
| Python dt=1.0 (nfreq=64) | 0.04–0.37 (diverge) | lejanas mal |
| **Python params legacy (nfreq=128)** | **0.80–0.91** (coincide) | **0.75–0.91, lag uniforme −2.5 s** |

→ Al igualar la discretización, **A05C converge al legacy**. La divergencia era
el **muestreo grueso** de producción (dt=1.0 → nfreq=64), que **sub-resuelve la
Green de las estaciones lejanas**. Queda solo un desfase **uniforme de ~2.5 s**
(convención de origen temporal entre los dos códigos; afecta a todas por igual,
benigno para el misfit interno Python-vs-observado).

## CAUSA RAÍZ DEFINITIVA (test 3): bug en el auto-`xl`

NO era el muestreo (dt). Aislando variables: cambiar dt no movió el misfit;
cambiar `xl` (periodicidad espacial / muestreo de número de onda) sí. Convergencia
del misfit vs observados (mejor modelo, dt=1.0):

| xl (m) | misfit |
|---|---|
| 0.5e6 | 0.351 |
| **0.673e6 (auto, producción)** | **0.178** |
| 0.8e6 | 0.103 |
| 0.9e6 | 0.095 |
| **≥1.0e6 (hasta 16e6)** | **0.093** (convergido) |

**Bug:** en `kdellipspy/axitra/src/axitra.py`, el auto-`xl` calcula
`vmax = max(model[:,2], model[:,3])`. El `model_array` es
`[depth, Vp, Vs, rho, Qp, Qs]`, así que col2,col3 = **(Vs, rho)** → vmax=4780,
**se salta Vp (col1=8480)**. Da `xl = 1.1·128·4780 = 0.673e6` (corto) en vez de
`1.1·128·8480 = 1.19e6` (convergido). La periodicidad la gobierna la onda más
rápida (Vp), no Vs.

→ A `xl=0.673e6` la integral de número de onda está **sub-muestreada para las
estaciones lejanas** (A05C/AC01/AC02), inflando el misfit a más del doble.

## Impacto y acción

- **TODAS las inversiones que usan el auto-`xl`** (2020, 2026, Copiapó) tienen el
  misfit inflado ~2× y modelos posiblemente subóptimos (la NA optimizó contra un
  forward sub-resuelto).
- **Fix:** corregir el auto-`xl` para usar Vp (col1), o `max` sobre Vp y Vs
  `max(model[:,1], model[:,2])`. Alternativa inmediata: pasar `xl` explícito ≥1.5e6.
- **Re-correr** las inversiones con el `xl` corregido (el misfit cae ~a la mitad).

_Tests: `/tmp/legacy_fwd/`, `/tmp/pyfwd_work/` (borrables). No se modificó el legacy._
