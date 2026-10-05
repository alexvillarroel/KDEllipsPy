# Vr no está resuelto por el cinemático (y por qué)

`vr_dt0_tradeoff.py`, caso `np1_isola_potin_b0.04-0.15_n20_wide_noA08FEA09F`,
con los otros seis parámetros fijos en el mejor modelo:

- **No hay mínimo en Vr.** El misfit baja monótonamente hasta el borde del
  barrido: el "mínimo global" queda en Vr = 7.50 km/s (1.57 beta) porque ahí se
  cortó, no porque haya un óptimo.
- El valle es la curva `t = L/Vr + dt0 = cte`: ruptura más rápida, hora de
  origen más tardía. El contorno de misfit 0.125 va de ~3.2 km/s a más de 7.5
  km/s, un factor 2.3 dentro del 5% del mejor ajuste.
- Imponer sub-Rayleigh (Vr <= 0.92 beta = 4.39 km/s) cuesta **+3.3%**
  (0.1228 contra 0.1188).
- Llevado al límite, el modelo prefiere Vr -> infinito con dt0 -> ~2 s, es
  decir una **fuente instantánea**: entre 0.04 y 0.15 Hz no se ve el efecto de
  fuente finita de una aspereza de este tamaño.

## Consecuencias

1. **No se puede afirmar supershear desde el cinemático.** Que el rango
   [5.0, 7.5] dé un misfit 3-4% mejor que [1.0, 4.5] no es evidencia: es el
   valle. Ver el commit del test de supershear.
2. **La dispersión entre semillas NO es la incertidumbre de Vr.** Las 5
   semillas dan std 0.13-0.30 km/s, pero eso mide convergencia del NA dentro
   del valle, no el ancho del valle. Reportar Vr como no resuelto, con el rango
   del valle.
3. El tope original de 4.0 km/s del input.ctl sí estaba mal puesto (0.84 beta,
   por debajo de la hipótesis sub-Rayleigh estándar de 0.9 beta = 4.29 km/s),
   pero corregirlo no resuelve el parámetro: solo mueve dónde se pega.
4. La pregunta corresponde al **dinámico**, donde Vr no es libre sino que sale
   del slip-weakening. fd3d puede producir supershear si la razón S es
   suficientemente baja, así que el test no está sesgado de entrada.

## Números de contexto

beta = 4.770 km/s en el hipocentro (Potin et al. 2024, capa de 112.5 km).
0.9 beta = 4.29 | 0.92 beta (Rayleigh) = 4.39 | sqrt(2) beta = 6.75 km/s.

Con semiejes de ~5-6.5 km, la ruptura recorre 6.5-13 km y la diferencia de
tiempo entre 0.9 beta y 1.26 beta es de 0.43 a 0.86 s: un 6-13% del periodo
más corto de la banda (6.7 s a 0.15 Hz). El valor absoluto de dt0 varía entre
configuraciones de -0.11 a +2.83 s, mucho más que eso.
