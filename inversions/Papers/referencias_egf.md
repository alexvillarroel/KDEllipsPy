# Referencias — Funciones de Green Empíricas (EGF)

Lista curada para entender el método EGF (de fundacional → método/límites → razón espectral).
Los DOIs marcados ✓ están verificados; el resto convendría confirmarlos al descargar.

## Fundacionales (de dónde viene la idea)

- **Hartzell, S. H. (1978).** *Earthquake aftershocks as Green's functions.*
  Geophysical Research Letters, 5(1), 1–4. ✓ doi:[10.1029/GL005i001p00001](https://doi.org/10.1029/GL005i001p00001)
  → Paper original: usar un sismo pequeño co-localizado como Green's function empírica.

- **Mueller, C. S. (1985).** *Source pulse enhancement by deconvolution of an empirical Green's function.*
  Geophysical Research Letters, 12(1), 33–36. ✓ doi:[10.1029/GL012i001p00033](https://doi.org/10.1029/GL012i001p00033)
  → La deconvolución EGF para recuperar la función fuente (STF). Lo más directo a lo que haces.
  (PDF abierto en USGS: https://pubs.usgs.gov)

## Método y sus límites (lo más útil en la práctica)

- **Hough, S. E. (1997).** *Empirical Green's function analysis: Taking the next step.*
  Journal of Geophysical Research: Solid Earth, 102(B3), 5369–5384. ✓ doi:[10.1029/96JB03488](https://doi.org/10.1029/96JB03488)
  → Extensión del método + inversión de atenuación común; condiciones que debe cumplir el EGF.
  (PDF abierto en USGS: https://pubs.usgs.gov/publication/70019415)

- **Abercrombie, R. E. (2015).** *Investigating uncertainties in empirical Green's function
  analysis of earthquake source parameters.* JGR: Solid Earth, 120(6), 4263–4277.
  ✓ doi:[10.1002/2015JB011984](https://doi.org/10.1002/2015JB011984)
  → IMPRESCINDIBLE: cuantifica incertidumbres y errores típicos. Leer antes de confiar en un Δσ/STF.

## Razón espectral / stress drop (si vas por escalamiento)

- **Prieto, G. A., Shearer, P. M., Vernon, F. L., & Kilb, D. (2004).** *Earthquake source
  scaling and self-similarity estimation from stacking P and S spectra.*
  JGR: Solid Earth, 109, B08310. ✓ doi:[10.1029/2004JB003084](https://doi.org/10.1029/2004JB003084)
  → Razón espectral con EGF (stacking) para estimar escalamiento y self-similarity.

- **Imanishi, K., & Ellsworth, W. L. (2006).** *Source scaling relationships of microearthquakes
  using a Green's function method.* En *Earthquakes: Radiated Energy and the Physics of Faulting*,
  AGU Geophysical Monograph 170, 81–90. (verificar DOI del capítulo)
  → Método de razón espectral con EGF aplicado a microsismos.

---
Ruta sugerida de lectura para tu caso (evento profundo intraplaca, recuperar STF/duración):
**Mueller 1985** (qué/cómo) → **Abercrombie 2015** (cuándo NO creértelo).
