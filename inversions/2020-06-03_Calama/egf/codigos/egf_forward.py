"""
egf_forward.py
==============
Operador forward con Funcion de Green Empirica (EGF) para una fuente finita:

    synth_i(t) = sum_k  (m_k / M0_egf) * EGF_i( t - tau_k,i )

con   tau_k,i = t_r(k) + [ r(k,i) - r(hipo,i) ] / c

  t_r(k)   : tiempo de ruptura de la subfalla k = dist(hipo->k)/Vr
  r(k,i)   : distancia subfalla k -> estacion i (termino de directividad)
  c        : velocidad de fase efectiva (~Vs en la fuente para ondas S)

Una sola EGF por estacion + corrimientos: valido porque la ruptura
(~10-20 km) << distancia fuente-estacion (>130 km), de modo que todas las
subfallas comparten esencialmente el mismo trayecto; solo cambia el retardo.

Supuestos / limites de validez:
  - EGF co-localizada y mecanismo similar (validar antes con validar_egf.py).
  - Banda de uso por encima de la esquina del mainshock y por debajo de la
    esquina de la EGF (aqui 0.1-0.3 Hz: f_c(M6.8)~0.05, f_c(M4.6)~0.6 Hz),
    de modo que la EGF actua como delta efectiva.
  - m_k / M0_egf reparte el momento; no se deconvoluciona la STF de la EGF.

El corrimiento se aplica en el dominio de la frecuencia (retardo fraccional
exacto). Sin dependencias de kdellipspy: recibe posiciones/momentos de las
subfallas como arrays, asi se puede alimentar desde
AxitraForwardModel.build_geometry_with_ellipse_slip(...).to_axitra_sources().

Uso (test de cordura):
    python codigos/egf_forward.py
"""
import numpy as np

R_TIERRA = 6371.0


def geo2km(lat, lon, lat0, lon0):
    """Conversion local lat/lon -> km (plano tangente; valido a estas escalas)."""
    x = (np.asarray(lon) - lon0) * np.cos(np.deg2rad(lat0)) * np.pi / 180.0 * R_TIERRA
    y = (np.asarray(lat) - lat0) * np.pi / 180.0 * R_TIERRA
    return x, y


def shift_frac(data, dt, tau):
    """Retardo fraccional exacto via FFT (tau>0 retrasa la senal)."""
    n = len(data)
    spec = np.fft.rfft(data)
    f = np.fft.rfftfreq(n, d=dt)
    return np.fft.irfft(spec * np.exp(-2j * np.pi * f * tau), n=n)


class EGFForward:
    """Forward de fuente finita por suma de EGFs corridas en tiempo.

    Parameters
    ----------
    egf       : ndarray (n_sta, 3, npts)  EGF por estacion/componente.
    dt        : float   paso temporal de la EGF [s].
    m0_egf    : float   momento sismico de la EGF [Nm].
    sta_lat, sta_lon : (n_sta,) coordenadas de las estaciones.
    hypo      : (lat, lon, depth_km) hipocentro del mainshock (= ancla EGF).
    c_phase   : float   velocidad de fase efectiva [km/s] (def. 4.6 ~ Vs fuente).
    """

    def __init__(self, egf, dt, m0_egf, sta_lat, sta_lon, hypo, c_phase=4.6):
        self.egf = np.asarray(egf, float)
        self.dt = float(dt)
        self.m0_egf = float(m0_egf)
        self.hypo = hypo
        self.c = float(c_phase)
        lat0, lon0, z0 = hypo
        sx, sy = geo2km(np.asarray(sta_lat), np.asarray(sta_lon), lat0, lon0)
        # Posicion 3D de las estaciones rel. al hipocentro (superficie z=0).
        self.sta_xyz = np.column_stack([sx, sy, np.full_like(sx, -z0)])

    def predict(self, sub_lat, sub_lon, sub_z_km, sub_moment, vr_kms,
                t_rup=None):
        """Sintetico (n_sta, 3, npts) para una distribucion de subfallas.

        sub_lat/lon/z_km : posiciones de las K subfallas (z positivo abajo).
        sub_moment       : momentos m_k [Nm] (solo entran las con m_k>0).
        vr_kms           : velocidad de ruptura [km/s].
        t_rup            : (K,) tiempos de ruptura precalculados [s] (p.ej. de
                           FaultGeometry.to_axitra_hist()); si None se calculan
                           como dist(hipo->k)/vr_kms.
        """
        lat0, lon0, z0 = self.hypo
        kx, ky = geo2km(np.asarray(sub_lat), np.asarray(sub_lon), lat0, lon0)
        kz = np.asarray(sub_z_km, float) - z0          # rel. al hipocentro
        m = np.asarray(sub_moment, float)
        act = m > 0
        if not np.any(act):
            raise ValueError("Ninguna subfalla con momento > 0.")
        if t_rup is not None:
            t_rup = np.asarray(t_rup, float)[act]
        kx, ky, kz, m = kx[act], ky[act], kz[act], m[act]
        sub_xyz = np.column_stack([kx, ky, kz])

        # Tiempo de ruptura de cada subfalla.
        if t_rup is not None:
            t_r = t_rup
        else:
            t_r = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) / float(vr_kms)

        n_sta, _, npts = self.egf.shape
        synth = np.zeros_like(self.egf)
        for i in range(n_sta):
            # Directividad: diferencia de caminos subfalla->estacion.
            r_k = np.linalg.norm(sub_xyz - self.sta_xyz[i], axis=1)
            r_0 = np.linalg.norm(self.sta_xyz[i])
            tau = t_r + (r_k - r_0) / self.c
            w = m / self.m0_egf
            for c in range(3):
                base = self.egf[i, c]
                spec = np.fft.rfft(base)
                f = np.fft.rfftfreq(npts, d=self.dt)
                # Suma vectorizada de retardos: sum_k w_k exp(-2pi i f tau_k)
                phase = np.exp(-2j * np.pi * np.outer(tau, f))
                synth[i, c] = np.fft.irfft((w[:, None] * phase).sum(0) * spec,
                                           n=npts)
        return synth


# ---------------------------------------------------------------------------
# Tests de cordura
# ---------------------------------------------------------------------------
def _ricker(t, t0, f0):
    a = (np.pi * f0 * (t - t0)) ** 2
    return (1 - 2 * a) * np.exp(-a)


def _test():
    dt, npts = 0.05, 2048
    t = np.arange(npts) * dt
    hypo = (-23.247, -68.530, 123.4)
    sta_lat = np.array([-22.85, -23.90])
    sta_lon = np.array([-70.20, -69.29])
    egf = np.zeros((2, 3, npts))
    for i in range(2):
        for c in range(3):
            egf[i, c] = _ricker(t, 30.0 + 3 * i, 0.2) * (1 + 0.2 * c)
    m0_egf = 1.0e16

    fwd = EGFForward(egf, dt, m0_egf, sta_lat, sta_lon, hypo)

    # 1) Limite puntual: toda la energia en el hipocentro -> EGF escalada.
    ratio = 2000.0
    syn = fwd.predict([hypo[0]], [hypo[1]], [hypo[2]], [ratio * m0_egf], vr_kms=3.0)
    err = np.max(np.abs(syn - ratio * egf)) / np.max(np.abs(ratio * egf))
    assert err < 1e-10, f"limite puntual fallo: err={err:.2e}"
    print(f"[OK] limite puntual: synth = {ratio:.0f} x EGF (err {err:.1e})")

    # 2) Dos subfallas separadas 10 km en profundidad, Vr=3 -> retardo ~3.3 s
    syn2 = fwd.predict([hypo[0], hypo[0]], [hypo[1], hypo[1]],
                       [hypo[2], hypo[2] + 10.0],
                       [ratio * m0_egf, ratio * m0_egf], vr_kms=3.0)
    # el pico de la suma debe estar entre el pico simple y el retardado
    p0 = t[np.argmax(egf[0, 0])]
    p2 = t[np.argmax(syn2[0, 0])]
    assert p0 <= p2 <= p0 + 4.0, f"retardo fuera de rango: {p0} vs {p2}"
    print(f"[OK] directividad/retardo: pico {p0:.2f}s -> {p2:.2f}s con 2a "
          "subfalla a 10 km / Vr=3")

    # 3) Conservacion de momento: amplitud espectral a f->0 escala con sum(m).
    s1 = np.abs(np.fft.rfft(syn[0, 0]))[1]
    s2 = np.abs(np.fft.rfft(syn2[0, 0]))[1]
    assert abs(s2 / s1 - 2.0) < 0.05, f"momento no conservado: {s2/s1:.3f}"
    print(f"[OK] momento total: razon espectral LF = {s2/s1:.3f} (~2 esperado)")

    print("\nTodos los tests de cordura pasaron.")


if __name__ == "__main__":
    _test()
