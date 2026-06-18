"""
integracion.py
==============
Port a Python de los codigos MATLAB de integracion del grupo (R. Boroschek):

    Filtro/integraf.m   -> integraf()   (integracion trapezoidal via filtro)
    Filtro/time.m       -> time_vector()
    Filtro/pssa3.m      -> pssa3()       (aceleracion -> velocidad -> desplazamiento
                                          con filtrado Butterworth)

La idea es exactamente la misma que en MATLAB:
  * ``integraf`` integra con el metodo del trapecio implementado como un filtro
    IIR  b=[1 1]/2 , a=[1 -1] , escalado por dt, y se le resta el primer valor.
  * ``pssa3`` aplica un Butterworth (pasa-banda / pasa-alto / pasa-bajo) y va
    integrando: acc -> vel -> desp, re-filtrando despues de cada integracion.
    Con ``nophase=True`` se usa filtfilt (acausal, fase cero), igual que el
    ``filtfilt`` de MATLAB.

Diferencia respecto al MATLAB: aqui ``f1``/``f2`` se pasan como ``None`` (en vez
del string 'NO') para indicar "sin esa esquina".
"""

import numpy as np
from scipy.signal import butter, filtfilt, lfilter, detrend


def time_vector(n, Fs):
    """Vector de tiempo de longitud ``n`` muestreado a ``Fs`` (port de time.m)."""
    return np.arange(n) / Fs


def integraf(d, Fs):
    """Integracion trapezoidal por filtrado (port fiel de integraf.m).

    Parameters
    ----------
    d  : array_like  senal a integrar
    Fs : float       frecuencia de muestreo [Hz]

    Returns
    -------
    ix : ndarray  integral de ``d``, con ix[0] = 0.
    """
    d = np.asarray(d, dtype=np.float64).ravel()
    dt = 1.0 / Fs
    b = np.array([1.0, 1.0]) / 2.0
    a = np.array([1.0, -1.0])
    ix = lfilter(b, a, d) * dt
    ix = ix - ix[0]
    return ix


def _design_butter(Fs, f1, f2, N):
    """Disena el Butterworth segun que esquinas esten definidas.

    f1=None, f2=val -> pasa-bajos en f2
    f1=val, f2=None -> pasa-altos en f1
    f1=val, f2=val  -> pasa-banda [f1, f2]
    """
    nyq = Fs / 2.0
    if f1 is None and f2 is not None:
        f2 = min(f2, 0.999 * nyq)
        return butter(N, f2 / nyq, btype="low")
    if f1 is not None and f2 is None:
        return butter(N, f1 / nyq, btype="high")
    if f1 is not None and f2 is not None:
        f2 = min(f2, 0.999 * nyq)
        return butter(N, [f1 / nyq, f2 / nyq], btype="band")
    raise ValueError("Al menos una de f1/f2 debe estar definida para filtrar.")


def pssa3(acc, Fs, f1=0.167, f2=25.0, N=4, nophase=True, izero=True):
    """Aceleracion -> (acc corregida, velocidad, desplazamiento).

    Port de pssa3.m. Procesa la senal aplicando un Butterworth y dos
    integraciones trapezoidales, re-filtrando tras cada integracion.

    Parameters
    ----------
    acc     : array_like  serie de aceleracion.
    Fs      : float        frecuencia de muestreo [Hz].
    f1, f2  : float | None esquinas del filtro [Hz]. Usa ``None`` para omitir
                           una de ellas (pasa-alto o pasa-bajo). Si ambas son
                           ``None`` no se filtra (solo detrend + integracion).
    N       : int          orden del Butterworth (def 4).
    nophase : bool         True -> filtfilt (fase cero, acausal); False -> filter.
    izero   : bool         True -> agrega ceros al inicio/fin antes de filtrar
                           (reduce efectos de borde) y luego los recorta.

    Returns
    -------
    ac, v, d : ndarray  aceleracion corregida, velocidad y desplazamiento,
                        todos con la misma longitud que ``acc``.
    """
    acc = np.asarray(acc, dtype=np.float64).ravel()
    m = acc.size

    filtrar = not (f1 is None and f2 is None)

    def _filt(b, a, x):
        return filtfilt(b, a, x) if nophase else lfilter(b, a, x)

    if filtrar:
        b, a = _design_butter(Fs, f1, f2, N)
        nfi = len(b) * 10  # padding de ceros (igual que en pssa3.m)

        if izero:
            ac = np.concatenate([np.zeros(nfi), detrend(acc), np.zeros(nfi)])
        else:
            ac = detrend(acc)

        ac = _filt(b, a, ac)

        v = integraf(ac, Fs)
        v = detrend(v)
        v = _filt(b, a, v)

        d = integraf(v, Fs)
        d = _filt(b, a, d)

        if izero:
            ac = ac[nfi:nfi + m]
            v = v[nfi:nfi + m]
            d = d[nfi:nfi + m]
    else:
        ac = detrend(acc)
        v = integraf(ac, Fs)
        d = integraf(v, Fs)

    return ac, v, d


if __name__ == "__main__":
    # Prueba rapida: integrar un seno -> deberia dar -cos escalado.
    Fs = 200.0
    t = time_vector(2000, Fs)
    sig = np.sin(2 * np.pi * 1.0 * t)
    ac, v, d = pssa3(sig, Fs, f1=0.1, f2=40.0)
    print("len acc/vel/disp:", ac.size, v.size, d.size)
    print("max |vel|:", np.max(np.abs(v)), " max |disp|:", np.max(np.abs(d)))
