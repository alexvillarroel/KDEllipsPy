# EGF — Calama 2020 (Mw 6.8, z=123 km)

Objetivo: sustituir la Green 1D (AXITRA) por una **Green empirica** (replica
co-localizada) para poder trabajar en la banda **0.1-0.3 Hz**, donde el modelo
1D ya no ajusta (misfit ~0.67) pero donde lambda ~ L permitiria por fin
restringir la geometria de la ruptura.

## Pipeline (codigos/)

1. `buscar_egf.py`   — ranking de candidatos (USGS+GFZ) -> `candidatos_ranking.csv`.
2. `descargar_egf.py [n]` — baja ondas + respuestas (GFZ/CX y EarthScope/C1)
   para los top-n -> `<tag>/RAW/`.
3. `validar_egf.py [tag]` — SNR + correlacion cruzada mainshock<->EGF en
   0.1-0.3 Hz -> `INFORME_egf.md` + `validacion_<tag>.png`.
   Criterio: CC_S mediana >= 0.6 valida; < 0.4 probar otro candidato.
4. `egf_forward.py`  — operador shift-and-sum (una EGF + retardos de ruptura y
   directividad). Tests de cordura incluidos (`python codigos/egf_forward.py`).

## Candidato top (ranking 2026-06-12)

- **2020-10-25T21:32:42 M4.6**, d3D = 5.5 km del hipocentro, dM = 2.2.
  Ventana espectral valida: fc(M6.8) ~ 0.05 Hz < banda 0.1-0.3 < fc(M4.6) ~ 0.6 Hz.
- Respaldo: 2020-06-07 M4.0 (d3D 3.6 km, dM 2.8 — SNR mas justo).

## Proximos pasos (Fase 2, manana)

- Si la validacion da CC alta: conectar `EGFForward.predict()` con la
  geometria elipse de kdellipspy (`build_geometry_with_ellipse_slip` ->
  `to_axitra_sources(latlon=True)` + momentos por subfalla) y armar el misfit
  en 0.1-0.3 con las estaciones validadas.
- M0_egf: estimar de Mw del catalogo (M0 = 10^(1.5*Mw + 9.05)); si hay MT
  regional publicado, usarlo.
- Limitacion: la inversion EGF solo usa las estaciones con EGF valida
  (PB05/PB06/PB08/PB19 + las C1 que aparezcan) — menos cobertura azimutal que
  las 8-9 del setup AXITRA; evaluar mezcla EGF (forma) + AXITRA (amplitud LF).
