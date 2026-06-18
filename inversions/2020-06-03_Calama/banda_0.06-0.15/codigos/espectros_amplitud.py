"""
espectros_amplitud.py
=====================
Espectro de amplitud de DESPLAZAMIENTO de todas las estaciones del evento
Calama 2020-06-03, para visualizar el *plateau* de baja frecuencia (proporcional
al momento sismico M0) y la *frecuencia esquina* (fc).

Los registros crudos (registros_20200603/<EST>/*.txt) son ACELERACION a 200 Hz.
El espectro de desplazamiento se obtiene en el dominio de la frecuencia como

        |U(f)| = |A(f)| / (2*pi*f)^2

evitando la deriva numerica de integrar dos veces en el tiempo. Se combina el
aporte de las 3 componentes como suma vectorial:

        |U|(f) = sqrt(|U_N|^2 + |U_E|^2 + |U_Z|^2)

Las estaciones USADAS en la inversion se dibujan en ROJO; el resto en gris.
Salida: output/espectros_amplitud.pdf

Uso:
    python codigos/espectros_amplitud.py
"""
import os
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from obspy import UTCDateTime

from io_sismica import read_all

# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TXT_ROOT = os.path.join(BASE, "registros_20200603")
OUT_PDF = os.path.join(BASE, "output", "espectros_amplitud.pdf")

# Hora de origen corregida del evento (07:35:34, no 07:35:28).
ORIGIN = UTCDateTime("2020-06-03T07:35:34.000")

# Estaciones efectivamente usadas en la inversion (real_disp.py SELECTED).
USADAS = ["A22F", "A08F", "A05F", "A18F", "A23F", "PB19", "PB06", "A09F", "A06F"]

# Ventana de analisis relativa al origen (incluye P y S del campo regional).
WIN_PRE, WIN_POST = 5.0, 120.0     # s antes/despues del origen
TAPER = 0.05                       # cosine taper (fraccion)

# Banda de la inversion (se sombrea como referencia).
BAND = (0.1, 0.3)                  # Hz

# Rango de frecuencias a mostrar.
FMIN_PLOT, FMAX_PLOT = 0.02, 20.0  # Hz


def espectro_desplazamiento(tr):
    """Devuelve (freqs, |U(f)|) del desplazamiento a partir de una traza de
    aceleracion, calculado en el dominio de la frecuencia."""
    tr = tr.copy()
    tr.detrend("linear")
    tr.detrend("demean")
    # Ventana alrededor del origen.
    tr.trim(ORIGIN - WIN_PRE, ORIGIN + WIN_POST, pad=True, fill_value=0.0)
    tr.taper(max_percentage=TAPER, type="cosine")

    fs = tr.stats.sampling_rate
    data = tr.data.astype(np.float64)
    n = len(data)
    acc = np.fft.rfft(data)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    # |U| = |A| / (2 pi f)^2 ; evitar division por 0 en f=0.
    w = 2.0 * math.pi * freqs
    disp = np.zeros_like(np.abs(acc))
    nz = freqs > 0
    disp[nz] = np.abs(acc[nz]) / (w[nz] ** 2)
    return freqs, disp


def main():
    st = read_all(TXT_ROOT)
    # Agrupar trazas por estacion.
    estaciones = {}
    for tr in st:
        estaciones.setdefault(tr.stats.station, []).append(tr)
    nombres = sorted(estaciones.keys())
    print(f"Estaciones leidas: {len(nombres)}")

    # Espectro vectorial por estacion.
    espectros = {}
    for est in nombres:
        comps = estaciones[est]
        f_ref, suma = None, None
        for tr in comps:
            f, u = espectro_desplazamiento(tr)
            if f_ref is None:
                f_ref = f
                suma = np.zeros_like(u)
            if len(u) == len(suma):
                suma += u ** 2
        espectros[est] = (f_ref, np.sqrt(suma))

    # Grilla de subplots.
    n = len(nombres)
    ncols = 7
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.4, nrows * 2.1),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, est in zip(axes, nombres):
        f, u = espectros[est]
        mask = (f >= FMIN_PLOT) & (f <= FMAX_PLOT)
        usada = est in USADAS
        color = "tab:red" if usada else "0.45"
        lw = 1.4 if usada else 0.9
        ax.loglog(f[mask], u[mask], color=color, lw=lw)
        # Banda de inversion.
        ax.axvspan(BAND[0], BAND[1], color="tab:red" if usada else "0.7",
                   alpha=0.12 if usada else 0.10, lw=0)
        ax.set_title(est, fontsize=9,
                     color="tab:red" if usada else "0.2",
                     fontweight="bold" if usada else "normal")
        ax.set_xlim(FMIN_PLOT, FMAX_PLOT)
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=6)

    # Apagar ejes sobrantes.
    for ax in axes[len(nombres):]:
        ax.axis("off")

    # Etiquetas comunes.
    fig.supxlabel("Frecuencia (Hz)", fontsize=11)
    fig.supylabel("Amplitud espectral de desplazamiento  |U(f)|  (m·s)", fontsize=11)
    fig.suptitle(
        "Espectros de amplitud (desplazamiento) — Calama 2020-06-03\n"
        "rojo = estaciones usadas en la inversion ; banda sombreada = "
        f"{BAND[0]}–{BAND[1]} Hz",
        fontsize=12)

    fig.tight_layout(rect=[0.02, 0.02, 1, 0.96])

    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    with PdfPages(OUT_PDF) as pdf:
        pdf.savefig(fig, dpi=200)
    plt.close(fig)
    print(f"-> {OUT_PDF}")
    print(f"   {len(nombres)} estaciones, {len(USADAS)} en rojo (usadas).")


if __name__ == "__main__":
    main()
