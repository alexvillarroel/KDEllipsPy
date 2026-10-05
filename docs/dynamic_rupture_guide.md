# Ruptura dinámica con KDEllipsPy + fd3d_TSN: guía para un evento nuevo

Escrita a partir de Navidad 2019-08-01 (interplaca, Mw 6.6), para aplicar el
mismo flujo a los sismos **intraplaca (intraslab) del norte**. Lee primero la
sección 5: varias cosas validadas en Navidad **no** están validadas para fallas
normales profundas.

- Repo: `/home/alex/KDEllipsPy`, rama `navidad-dynamic` (GitHub `alexvillarroel/KDEllipsPy`).
- Python: `/home/alex/.conda/envs/kdellipspy/bin/python` (pygmt, neighpy, okada_wrapper).
- fd3d_TSN (Premus/Gallovič): `/home/alex/fd3d_TSN/src`. Se compila solo desde `forward_dynamic.py`.
- Caso de referencia completo: `inversions/2019-08-01_Navidad/` (el `dynamic/` es la plantilla).
- Tests: `kdellipspy/test_{dynamic_convolution,dynamic_forward_model,tsn_bridge,stress_mapper,fast_axitra,dynamic_solver_config}.py`
  y `tests/`. Correr con `rtk pytest` antes y después de tocar código del paquete.

## 1. Cómo funciona el acoplamiento

```
modelo elíptico (a, b, u, v, phi, Te, cte1, Dc, dt0)
  -> EllipticalStressMapper / build_tsn_dynamic_fields   (kdellipspy/core/geometry.py)
     prestress T0, resistencia de pico, Dc en la grilla gruesa (2 km)
  -> fd3d_TSN (slip-weakening, DIPSLIP, FSPACE)          (tsn_bridge.run_tsn_forward)
     slip-rate en la falla, celdas de 500 m, dt 0.015 s
  -> bin_slip_rate_to_subfaults + project_rake           (dynamic_convolution.py)
     momento acumulado M(t) por subfalla (NO la tasa)
  -> base de respuestas impulsivas de AXITRA, rake 0/90  (response_basis / synthetics_from_basis)
     convolución analítica en frecuencia compleja
  -> sintéticos -> misfit L2 normalizado -> NA (Neighbourhood Algorithm)
```

- La nucleación es un parche circular sobreesforzado **en el hipocentro**
  (`tsn_hypocentre_coarse`), y la ruptura no puede salir de la elipse: fuera de
  ella hay una barrera (T0 = 0, pico = 1e10 Pa). Nucleación y detención están
  **impuestas**: el dinámico no puede demostrar dónde nace o se detiene.
- Modo **B (M0 fijo)**: por autosemejanza del slip-weakening (Te, Dc × k ⇒ slip × k)
  se fija M0 al de la agencia, se elimina Te (`TE_REF = 3`) y el parámetro libre
  pasa a ser la forma. En Navidad, B ajustó mejor que A (M0 libre).
- **Convolución rápida** (por defecto): se calcula una vez la base de Green por
  subfalla y después cada forward son milisegundos (~25 ms el cinemático). Esto
  es lo que hace viables muchas semillas, tests de sectores y NA-Bayes.

## 2. Convenciones y errores ya corregidos (no reintroducirlos)

| Tema | Regla | Qué pasaba |
|---|---|---|
| Barrera | fuera de la elipse: T0 = 0, pico = `TSN_BARRIER_PEAK_PA` (1e10) | con T0 = −1e8, fd3d compara \|tracción\| con la resistencia y la barrera rompía |
| STF de AXITRA | el tipo 3 de archivo es **M(t)**, no la tasa; integrar | sintéticos con forma y amplitud erróneas |
| Amortiguamiento | `aw = 2` está fijo en `axitra.py` (`.data`); leerlo con `_axitra_damping` | la base rápida no reproducía AXITRA |
| Eje de manteo | fd3d cuenta las filas **desde el borde profundo**; `bin_slip_rate_to_subfaults` invierte `[:, :, ::-1]` | slip espejado en el manteo |
| Hipocentro | en la grilla gruesa, 1-based, manteo desde el borde profundo (`tsn_hypocentre_coarse`) | nucleación en el lugar equivocado |
| Hora de origen | invertir **dt0** (s); en Navidad salió +3.9 s con el hipocentro ISC | una noche entera perdida por un desfase de ~3 s |
| Reconstruir modelos | **siempre** `to_full_model_named(result.param_names, values)` (`dynamic/na_dynamic.py`) | el modelo B (8 parámetros) se leyó con el mapeo de 9: parámetros corridos y análisis falsos |
| Velocidades | `crustal_rows` desplaza el modelo 1D para que z = 0 de fd3d sea el **borde superior de la falla** | — |
| Superficie libre | build `-DDIPSLIP -DFSPACE`, `nfs = nabc` (falla enterrada) | — |
| CFL | fd3d exige Vp·dt/dh < 0.25 (con dh = 500 m y dt = 0.015 s ⇒ Vp ≤ 8.33 km/s) | inestabilidad numérica |

## 3. Flujo para un evento nuevo

0. **Datos observados: siempre `kde-prep`** (`python -m kdellipspy.prep_data <caso>`),
   que integra desde la aceleración cruda (`SAC/ACC/`) con la cadena de ISOLA y la banda
   del propio `input.ctl`. El `codigos/real_disp.py` de cada evento **no** sirve: lee
   `SAC/DISP/` ya integrado, integra dos veces en el tiempo y filtra solo al final, lo que
   deja energía antes de la llegada P. `kde-prep` necesita `integracion.py` en
   `<caso>/codigos` o `<caso>/../codigos`.
1. **Cinemático robusto primero** (de él salen la geometría, la banda, dt0 y Vr de referencia).
   - Probar hipocentros (ISC, CSN, NEIC), modelo de velocidades regional y bandas,
     con **5 semillas NA por configuración** (`vel_tests/run_kin_multiseed.py`).
     Un ranking con una sola semilla engaña.
   - Bandas: subir la frecuencia máxima mientras Vr y dt0 se estabilicen. En Navidad,
     0.06–0.15 Hz no resolvía Vr (±0.96) y 0.05–0.30 Hz sí (±0.06), con más misfit.
   - Revisar el misfit por estación y sacar las que no se ajustan en ninguna banda (en Navidad, M14L).
   - Malla: Lx/dh y Ly/dh deben ser múltiplos de nx, ny (lo verifica `forward_dynamic.py`).
2. **Caso dinámico**: copiar `inversions/2019-08-01_Navidad/dynamic/` al evento y
   apuntar `DYN_CASE` al caso cinemático elegido. Ajustar `DT`/`NT` (CFL y duración).
   Usar un `DYN_WORK` distinto para cada corrida simultánea.
3. **Preflight obligatorio** (`na_dynamic.py` ya lo trae: `REFERENCE_FREE`, `DYN_PREFLIGHT_MAX`):
   un forward con un modelo razonable en **el caso exacto** que se va a invertir.
   Revisar el misfit, el Mw y la figura de formas de onda, y abortar si el misfit es malo.
4. **NA dinámico**: A (M0 libre) y B (`DYN_M0` = M0 de la agencia), varias semillas.
   En Navidad: NA 300 + 20 × 60; rangos a, b ∈ [2, 8] puntos gruesos; dt0 ∈ [−3, 8] s.
5. **Análisis** (todo en `inversions/2019-08-01_Navidad/`, se puede copiar):
   - `dynamic/report_dynamic.py`, `dynamic/map_dynamic.py`: ajuste, slip y métricas.
   - `rupture/moment_release.py`, `traction_plots.py`, `movie_B.py`: función fuente,
     diagramas espacio–tiempo, isócronas, tracción, película mp4.
   - `okada/okada_models.py`: deformación estática (figuras uh y uz por separado).
   - `directivity/run_directivity.py`: NA con el centroide obligado a un sector (resuelve la directividad).
   - `egf/astf.py`: funciones fuente aparentes con funciones de Green empíricas (necesita una réplica con un mecanismo similar).
   - `uncertainty/na_bayes.py`: posterior NA-Bayes **calibrada** (sección 4).
6. Dashboard (`dashboard/build_data.py` + `template.html`) si se quiere compartir.

## 4. Incertidumbres (NA-Bayes)

- `run_appraisal()` (`kdellipspy/inversion/base.py`) usa por defecto T = 2·misfit_min,
  lo que equivale a N = 1 dato independiente: un posterior casi igual al prior. **No usar el default.**
- Verosimilitud calibrada: log L = −N·misfit / (2·misfit_min), con
  N = 2 × ancho de banda × duración de la ventana × número de trazas. Correr también N/4 como sensibilidad.
- **No usar `n_walkers > 1` de neighpy**: se colgó 5 h sin CPU. Correr una cadena
  por proceso (`CHAIN=k`, `n_walkers=1`) y combinar (así está `na_bayes.py`).
- Reportar **magnitudes derivadas** (Mw, centroide, desplazamiento del centroide
  respecto del hipocentro en rumbo y manteo, dimensiones): los ejes de la elipse
  tienen simetrías (a1↔a2 con θ+½) y sus marginales son multimodales.
- Revisar si algún parámetro queda pegado al borde del prior (en Navidad, dt0 en +5 s): si pasa, ampliar el rango.

## 5. Diferencias para sismos intraslab (VERIFICAR antes de invertir)

1. **Sentido del deslizamiento (lo más crítico). RESUELTO (Calama 2020, 2026-10-05).**
   Confirmado: fd3d_TSN solo recibe el buzamiento y dentro de la aspereza el prestress
   es +Te, así que con DIPSLIP el slip va **siempre** hacia +Z local → rake +90 (inverso).
   Para un mecanismo normal hay que invertir el vector de momento completo (rake + 180).
   Lo hace `DynamicForwardModel.slip_sign`: se deriva del rake del `input.ctl`
   (`sin(rake) < 0 → −1`) y niega `(Mx, Mz)` después de `bin_slip_rate_to_subfaults`.
   En el camino rápido (`use_basis`, el de por defecto) **no** interviene `project_rake`:
   los sintéticos son `base(rake 0)·Mx + base(rake 90)·Mz`, lineal, así que negar ambos
   es exactamente rake + 180. Navidad (rake +100) mantiene +1 y queda bit-idéntico.
   Test: `test_normal_rake_flips_moment_handed_to_axitra` (+90, +100, −85, −101).
   Verificación en datos (`inversions/2020-06-03_Calama/dynamic/check_polarity.py`,
   np1, un forward de fd3d, cada signo en su mejor dt0):

   | | misfit | dt0 | rake efectivo | cc media (ventana P) |
   |---|---|---|---|---|
   | fd3d tal cual (+1) | 0.688 | −2.0 s | +90.0° | +0.667 |
   | con flip (−1) | **0.403** | +3.0 s | **−90.0°** | **+0.870** |

   La cc es mayor con el flip en las 9 estaciones. **No** usar el signo del primer
   arribo P como métrica: los dos sintéticos son negativos exactos uno del otro, los
   conteos suman siempre el número de estaciones y no discrimina (dio 4/9 vs 5/9).
   **No** invertir el sentido poniendo Te negativo: la lógica de barrera trata T0 < 0
   como barrera (T0 = 0) y mataría la aspereza.
2. **CFL. VERIFICADO: sí se viola en Calama 2020.** No manda el Vp del hipocentro
   sino el **máximo** que ve fd3d. Con el modelo de Potin et al. (2024) la falla de
   np1 va de 97.7 a 133.4 km y la de np2 de 106.5 a 124.6 km, y ambas incluyen la
   capa de 112.5 km con **Vp = 8.415 km/s** → dt = 0.015 s da CFL **0.2524 > 0.25**.
   Se usa **dt = 0.0145 s** (CFL 0.2440) y NT = 2759 para los mismos 40 s.
   `forward_dynamic.cfl_dt_max(cfg)` lo recalcula desde `crustal_rows` del caso y lo
   verifica con un `assert` en `setup_work_dir`, para que no se cuele en silencio al
   cambiar de hipocentro o de modelo 1D.
3. **Esfuerzo normal.** `normstress` crece con la distancia desde el borde superior
   de la falla, no con la profundidad absoluta. Para slip-weakening con resistencia
   de pico explícita (cte1·Te) no debería importar, pero hay que tenerlo presente si se interpreta μ.
4. **Rangos.** Las caídas de esfuerzo intraslab suelen ser mayores: ampliar Te/cte1
   y revisar Dc. Comprobar la caída de esfuerzo final con Eshelby y la tracción (`traction_plots.py`).
5. **Plano nodal.** Hay ambigüedad entre np1 y np2 (el caso Calama 2020 ya tiene
   `np1/` y `np2/`): invertir ambos con el mismo protocolo y comparar misfits entre semillas.
   Ojo: `np1/` y `np2/` **no están trackeados** en git (solo `input.ctl`, `event.ctl` y
   `README.md`), así que un `input.ctl` sobrescrito se pierde sin dejar rastro. Pasó:
   el de np2 tenía la geometría de np1 y la correcta solo sobrevivía en `output/run.log`.
6. **Velocidades.** Comprobar que el modelo 1D cubra la profundidad de la falla y
   que `crustal_rows` encuentre la capa que contiene el borde superior.
7. **Directividad.** La geometría de estaciones decide qué se resuelve (en Navidad,
   todas al E: el manteo se resolvió y el rumbo no). Hacer el test de sectores antes de interpretar.

## 6. Reglas de trabajo con el usuario

- Prefijar los comandos de shell con `rtk` (ver `CLAUDE.md`).
- Trabajar en una rama y hacer commit y push de cada paso, con
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` (el usuario quiere poder volver atrás).
- **No** hacer commit en `/home/alex/Projects/Paper_Calama_2026` sin preguntar: tiene cambios propios pendientes.
- No cambiar la cuenta activa de `gh` sin preguntar.
- Figuras con textos **en inglés**; la conversación, en español.
- Corridas largas: preflight, **una sola cosa cambiada por vez**, y verificar con
  `ps -o pcpu` a los 2 minutos que el proceso de verdad avanza.
- Reportar con honestidad lo que no se resuelve (por ejemplo, la directividad en el rumbo)
  y corregir explícitamente las conclusiones anteriores cuando un test las contradice.
