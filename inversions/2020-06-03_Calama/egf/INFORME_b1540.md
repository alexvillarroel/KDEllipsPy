# Batería EGF banda 0.15–0.4 Hz + lags estáticos (2026-06-12 tarde)

Motivación: (1) espectros señal/ruido (espectros_snr.pdf) mostraron el pico del
microsismo justo en 0.1–0.3 Hz y SNR de la EGF 3× mejor en 0.15–0.4; (2) CC
mainshock↔EGF validada en la banda nueva (mediana 0.72, 8/8 estaciones, PB02
reincorporada; cc_bandas.csv); (3) sospecha de que la "aspereza compacta" de la
batería nocturna era artefacto de errores de alineación (~2–4 s) — se agregó al
misfit un lag estático por estación y ventana (P y S independientes, ±4 s) y se
excluyeron las ventanas P con CC<0.4 (PB02, PB09).

## Resultados (best model por run; 3000 evals NA c/u)

| run | banda | lags | misfit | a1 | a2 | θ(π) | np | tp | dmax | vr | área πa1a2 | Mw |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| principal | 0.15–0.4 | ±4 | 0.629 | 8.0 | 9.8 | 0.30 | 1.00 | 0.47 | 6.8 | 4.1 | 246 | 7.24 |
| s22 | 0.15–0.4 | ±4 | 0.556 | 8.5 | 12.3 | 0.45 | 0.82 | 0.07 | 7.4 | 1.8 | 330 | 7.35 |
| s33 | 0.15–0.4 | ±4 | 0.621 | 9.8 | 8.4 | 0.39 | 0.25 | 0.98 | 4.5 | 1.6 | 257 | 7.13 |
| sinPB02 | 0.15–0.4 | ±4 | 0.601 | 6.6 | 12.0 | 0.71 | 0.76 | 0.44 | 7.7 | 1.9 | 249 | 7.28 |
| **nolag** (control) | 0.15–0.4 | no | 0.767 | 7.7 | 3.8 | 0.23 | 0.86 | 0.70 | 7.3 | 1.5 | **92** | 6.97 |
| **b0103+lag** (control) | 0.1–0.3 | ±4 | 0.543 | 5.2 | 11.6 | 0.71 | 0.95 | 0.09 | 7.7 | 1.6 | 189 | 7.21 |

Contraste con aleatorios (mediana de los 500 iniciales → best):
con lags 0.91 → 0.54–0.63 (Δ≈0.3); nolag 1.44 → 0.77 (Δ=0.68, el mayor).
La banda nueva por sí sola baja el piso 0.87 → 0.77.

## Conclusión principal: la compacidad ERA artefacto de timing

Control limpio (misma banda, estaciones y exclusiones; solo lags on/off):
área 92 km² sin lags vs 246–330 km² con lags. La banda vieja con lags también
crece (189 km²). Los errores de alineación taup/localización (~2–4 s)
penalizaban las fuentes extendidas y favorecían artificialmente la aspereza
compacta (~40–80 km²) de la batería nocturna.

## Robusto ahora (runs con lags)

- **Fuente extendida**: a1 5.2–9.8, a2 8.4–12.3 km; área ~190–330 km².
- **dmax ~7 m** (rielaría sin el tope 8) en 5/6.
- Paisaje sigue discriminando (Δ≈0.3 sobre aleatorios).

## NO robusto / caveats

- θ, tp, vr (1.6–4.1) dispersos entre semillas; np alto en 5/6 (s33: 0.25).
- **Lags S rielando en +4.0 s** en varias estaciones (s22/s33: 4 c/u); la
  validación mostró residuos hasta +6.2 s → ±4 s es vinculante. Próxima
  iteración: lagmax 6 s, o alinear por lag CC medido y permitir ±2.
- **M0 sobreestimado** (Mw 7.0–7.35 vs 6.8): factor 3–6 en M0, conjugado con
  la incertidumbre de M0_egf (de Mw 4.89 de catálogo). Calibrar M0_egf por
  cociente espectral de plateaus a baja frecuencia mainshock/EGF antes de
  interpretar dmax/M0 en absoluto.
- Con lags, parte de la información de directividad se va a los lags: la
  geometría queda constreñida por forma/duración, no por moveout absoluto.

## Próximos pasos sugeridos

1. lagmax 6 s (o alineación por CC) + semillas, ver si θ/tp/vr se ordenan.
2. Calibrar M0_egf con cociente espectral LF (quita el factor de Mw).
3. EGF compuesta (M4.9+M4.6) para bajar el piso de ruido.

Archivos: inversion_b1540[_s22|_s33|_sinPB02|_nolag], inversion_b0103_lag,
bateria_b1540.log. Configs vía CLI de invertir_egf.py
(--banda/--csv/--lagmax/--sinP/--dmax).
