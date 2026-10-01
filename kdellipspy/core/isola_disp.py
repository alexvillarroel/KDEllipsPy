"""acc -> desplazamiento por la cadena de ISOLA, leyendo SAC directamente.

Nada de esto reimplementa el procesamiento: pssa3 (integracion.py del proyecto)
reproduce los .vel.sac con error 0.0000%, y filter_play es el Fortran de ISOLA
envuelto con f2py, verificado bit a bit contra ISOLA_MT/invert/*fil.dat
(cc 1.0000, error 3e-4 %, desfase 0 s en 18 trazas).

Sustituye a _preprocess_trace para units=1. Aquella integraba dos veces en el
tiempo y filtraba una sola vez al final, lo que dejaba una deriva cuadratica:
la energia anterior a la llegada P alcanzaba el 367% de la senal (PB19 N, 2022)
frente al 0-1% de ISOLA.
"""
import numpy as np
import obspy

from .signal_utils import integrate_waveforms

from .. import isola_src  # noqa: F401  (asegura que el paquete existe)
from ..isola_src import fplay

KEYDIS_DISP = 1         # filter_play: salida en desplazamiento (integra)
KEYDIS_VEL = 0          # salida en velocidad (solo filtra)
KEYFIL_BP = 0           # 0 = pasa banda
F1_ANCHO = 0.04         # esquina baja de acc_to_disp.py; NO es la banda de inversion
LOWPASS_RAW = 2.0       # raw_isola.py
TAPER_RAW = 0.15        # raw_isola.py
RATE_FILTRO = 10.0      # Hz: filtrar/integrar a 0.1 s, como ISOLA


def esquina_alta(sampling_rate, objetivo=40.0, tope=20.0, frac_nyquist=0.8):
    """Esquina alta utilizable segun el muestreo (port de acc_to_disp.py)."""
    ny = sampling_rate / 2.0
    if sampling_rate >= 100.0:
        return min(objetivo, frac_nyquist * ny)
    return min(tope, frac_nyquist * ny)


def sac_a_disp(sac_acc, origen, freq1, freq2, npts, delta, pssa3,
               keydis=KEYDIS_DISP, zerophase=False):
    """SAC de aceleracion -> desplazamiento (keydis=1) o velocidad (keydis=0).

    Controla la deriva de la doble integracion filtrando ANCHO y en FASE CERO
    entre pasos -- ancho para no distorsionar la banda util, fase cero para no
    introducir retardo -- y aplica el filtro ESTRECHO de inversion al final,
    con el mismo scipy/obspy y la misma causalidad que reciben los sinteticos.

    Esa simetria es el requisito: si el dato pasa por filter_play (XAPIIR) y el
    sintetico por scipy, aparece un desfase de 5-8 s entre ambos y el misfit se
    queda en ~1 aunque las trazas sean correctas. Por eso filter_play NO se usa
    aqui; queda disponible para reproducir ISOLA bit a bit (ver _selftest).
    """
    tr = obspy.read(str(sac_acc))[0]
    fs = tr.stats.sampling_rate
    t = tr.copy()
    t.detrend("linear"); t.detrend("demean")
    t.taper(max_percentage=0.05, type="cosine")
    ancho = dict(freqmin=F1_ANCHO, freqmax=esquina_alta(fs), corners=4, zerophase=True)
    t.filter("bandpass", **ancho)
    pasos = 1 if keydis == KEYDIS_VEL else 2
    for k in range(pasos):
        t.data = integrate_waveforms(t.data, float(t.stats.delta), steps=1)
        t.detrend("linear")
        if k < pasos - 1:                      # el ultimo filtro es el estrecho
            t.filter("bandpass", **ancho)
    # Filtro de inversion: mismo que los sinteticos (banda y causalidad del ctl)
    t.filter("bandpass", freqmin=freq1, freqmax=freq2, corners=4,
             zerophase=bool(zerophase))
    # margen de una muestra: con el limite justo, interpolate se sale por
    # redondeo y lanza 'No extrapolation can be performed'
    t.trim(origen, origen + npts * delta, pad=True, fill_value=0.0)
    t.interpolate(sampling_rate=1.0 / delta, starttime=origen, npts=npts)
    return t.data


def _selftest():
    """filter_play debe reproducir un fil.dat de ISOLA desde su raw.dat.

    No usa SAC a proposito: aisla el unico trozo que este modulo no controla.
    Si esto falla, el .so esta compilado contra otro Python o cambio ISOLA.
    """
    from pathlib import Path
    inv = Path("/home/alex/Projects/Paper_Calama_2026/Source_code/inversions"
               "/2022-07-27_Calama/ISOLA_MT/invert")
    if not (inv / "PB19raw.dat").exists():
        print("selftest omitido: no estan los ficheros de referencia")
        return
    raw = np.loadtxt(inv / "PB19raw.dat")
    ref = np.loadtxt(inv / "PB19fil.dat")[:, 1]
    b = np.asfortranarray(raw[:, 1], dtype=np.float32)
    fplay.filter_play(KEYFIL_BP, KEYDIS_DISP, 0.04, 0.09, raw[1, 0] - raw[0, 0], b)
    err = np.abs(b - ref).max() / np.abs(ref).max()
    assert err < 1e-3, f"filter_play no reproduce ISOLA: error {100*err:.4f}%"
    print(f"ok  filter_play reproduce ISOLA (error {100*err:.5f}%)")


if __name__ == "__main__":
    _selftest()
