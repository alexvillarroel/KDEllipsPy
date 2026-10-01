"""Funciones fuente aparentes (ASTF) de Navidad 2019 con función de Green empírica.

EGF: réplica 2019-08-01 20:01:28 (CSN -34.31, -72.42, 21 km; GCMT Mw 5.9,
358/13/87 ~ sismo principal 19/14/106). Registros de aceleración del CSN (evtdb)
en las estaciones comunes -> desplazamiento, 0.03-0.4 Hz, 10 Hz.

ASTF por estación y componente: deconvolución por Landweber proyectado
(ASTF >= 0, soporte [-5, 30] s respecto del origen). Cada registro se alinea en
la hora de origen de SU evento, así que el centroide temporal de la ASTF mide el
centroide del sismo respecto de la ubicación de la EGF:

    tc(sta) = tc0 - s_i · (r_c - r_EGF),   s_i = sin(i)/Vs * (dirección a la estación)

y la duración aparente (Haskell unilateral): T(az) = T0 - L·s_i·cos(az - phi).
Escribe astf.png, astf_directivity.png y astf_results.json.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from obspy import Stream, Trace, UTCDateTime  # noqa: E402
from obspy.geodetics import gps2dist_azimuth  # noqa: E402

HERE = Path(__file__).resolve().parent
MAIN = dict(t0=UTCDateTime("2019-08-01T18:28:03"), lat=-34.28, lon=-72.51, dep=13.0)   # CSN
EGF = dict(t0=UTCDateTime("2019-08-01T20:01:28"), lat=-34.31, lon=-72.42, dep=21.0)    # CSN
GCMT_RATIO = 1.86e26 / 8.53e24
FMIN, FMAX, FS = 0.03, 0.4, 10.0
WIN = 100.0                 # s de registro desde el origen
SUPPORT = (-3.0, 22.0)      # soporte de la ASTF (s)
VS_SRC = 3.9                # km/s cerca de la fuente
STRIKE = 19.0               # GCMT del principal


def read_txt(path):
    hdr = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            k, _, v = line.lstrip("# ").partition(":")
            hdr[k.strip()] = v.strip()
    data = np.loadtxt(path, comments="#")
    sta, cha = hdr["Estacion"].split()[0], hdr["Estacion"].split()[-1]
    lat, lon = [float(x) for x in hdr["Latitud"].replace("Longitud:", "").split()]
    tr = Trace(data, header=dict(station=sta, channel=cha, starttime=UTCDateTime(hdr["Tiempo de Origen"]),
                                 sampling_rate=float(hdr["Tasa de muestreo"].split()[0])))
    tr.stats.coordinates = (lat, lon)
    return tr


def to_disp(tr, t0):
    tr = tr.copy()
    tr.detrend("demean"); tr.detrend("linear"); tr.taper(0.05)
    tr.filter("highpass", freq=FMIN, corners=4, zerophase=False)
    tr.integrate(); tr.filter("highpass", freq=FMIN, corners=4, zerophase=False)
    tr.integrate(); tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=False)
    tr.resample(FS) if tr.stats.sampling_rate % FS else tr.decimate(int(tr.stats.sampling_rate / FS), no_filter=True)
    tr.trim(t0, t0 + WIN, pad=True, fill_value=0.0)
    return tr.data[: int(WIN * FS)]


def landweber(d, g, n_iter=300):  # ponytail: parada temprana = regularización
    """d = astf * g (convolución lineal, astf soportado en SUPPORT). Devuelve (t, astf, VR)."""
    n0, n1 = int(SUPPORT[0] * FS), int(SUPPORT[1] * FS)
    nfft = 1 << int(np.ceil(np.log2(2 * len(d) + n1 - n0)))
    G = np.fft.rfft(g, nfft)
    D = np.fft.rfft(d, nfft)
    tau = 1.0 / np.max(np.abs(G)) ** 2
    lag = np.arange(nfft); lag[lag > nfft // 2] -= nfft  # lags circulares (negativos al final)
    mask = (lag >= n0) & (lag < n1)
    f = np.zeros(nfft)
    for _ in range(n_iter):
        F = np.fft.rfft(f)
        r = D - G * F
        r = np.fft.rfft(np.fft.irfft(r, nfft)[: len(d)] * 1.0, nfft)  # residuo solo en la ventana de datos
        f = f + tau * np.fft.irfft(np.conj(G) * r, nfft)
        f[~mask] = 0.0
        np.maximum(f, 0.0, out=f)
    pred = np.fft.irfft(G * np.fft.rfft(f), nfft)[: len(d)]
    vr = 1 - np.sum((d - pred) ** 2) / np.sum(d ** 2)
    idx = np.r_[np.arange(n0, 0) % nfft, np.arange(0, n1)]
    return np.arange(n0, n1) / FS, f[idx], vr, pred


def moments(t, a):
    m0 = a.sum() / FS
    tc = (t * a).sum() / FS / m0
    sd = np.sqrt(((t - tc) ** 2 * a).sum() / FS / m0)
    c = np.cumsum(a) / a.sum()
    return m0, tc, 2 * sd, t[np.searchsorted(c, 0.05)], t[np.searchsorted(c, 0.95)]


def main():
    rows = []
    for d in sorted((HERE / "egf_2001").iterdir()):
        sta = d.name
        if not (HERE / "main" / sta).exists():
            continue
        for comp in ("HNZ", "HNN", "HNE"):
            fe = list(d.glob(f"*{comp}.txt")); fm = list((HERE / "main" / sta).glob(f"*{comp}.txt"))
            if not fe or not fm:
                continue
            te, tm = read_txt(fe[0]), read_txt(fm[0])
            if tm.stats.starttime > MAIN["t0"] or te.stats.starttime > EGF["t0"]:
                continue  # registro empieza después del origen
            g, u = to_disp(te, EGF["t0"]), to_disp(tm, MAIN["t0"])
            t, a, vr, pred = landweber(u, g)
            lat, lon = tm.stats.coordinates
            dist, az, _ = gps2dist_azimuth(EGF["lat"], EGF["lon"], lat, lon)
            rows.append(dict(sta=sta, comp=comp, lat=lat, lon=lon, dist_km=dist / 1e3, az=az, vr=float(vr),
                             t=t, a=a, u=u, pred=pred, mom=moments(t, a)))
            print(f"{sta} {comp} dist {dist/1e3:5.1f} az {az:5.1f}  VR {vr:.2f}  ratio {rows[-1]['mom'][0]:.1f} "
                  f"tc {rows[-1]['mom'][1]:.1f} T {rows[-1]['mom'][2]:.1f}", flush=True)

    good = [r for r in rows if r["vr"] > 0.6]
    stas = sorted({r["sta"] for r in good})
    st_rows = []
    for s in stas:
        rr = [r for r in good if r["sta"] == s]
        dist, az = rr[0]["dist_km"], rr[0]["az"]
        sin_i = dist / np.hypot(dist, EGF["dep"])
        st_rows.append(dict(sta=s, az=az, dist=dist, p=sin_i / VS_SRC, n=len(rr),
                            tc=float(np.mean([r["mom"][1] for r in rr])), T=float(np.mean([r["mom"][2] for r in rr])),
                            ratio=float(np.mean([r["mom"][0] for r in rr]))))
    az = np.radians([r["az"] for r in st_rows]); p = np.array([r["p"] for r in st_rows])
    tc = np.array([r["tc"] for r in st_rows]); T = np.array([r["T"] for r in st_rows])
    # tc = tc0 - p (dE sin az + dN cos az)
    A = np.c_[np.ones_like(tc), -p * np.sin(az), -p * np.cos(az)]
    sol_c, *_ = np.linalg.lstsq(A, tc, rcond=None)
    res_c = tc - A @ sol_c
    sol_T, *_ = np.linalg.lstsq(A, T, rcond=None)  # T = T0 - p (LE sin az + LN cos az)
    # jackknife
    jk = []
    for k in range(len(tc)):
        m = np.arange(len(tc)) != k
        jk.append(np.linalg.lstsq(A[m], tc[m], rcond=None)[0])
    jk = np.array(jk)
    se = np.sqrt((len(tc) - 1) / len(tc) * ((jk - jk.mean(0)) ** 2).sum(0))
    dE, dN = sol_c[1:]
    kx = 111.19 * np.cos(np.radians(EGF["lat"]))
    cen_lat, cen_lon = EGF["lat"] + dN / 111.19, EGF["lon"] + dE / kx
    # respecto del hipocentro CSN del principal
    hE, hN = (cen_lon - MAIN["lon"]) * kx, (cen_lat - MAIN["lat"]) * 111.19
    s = np.radians(STRIKE)
    along, across = hE * np.sin(s) + hN * np.cos(s), hE * np.cos(s) - hN * np.sin(s)
    LE, LN = sol_T[1:]
    out = dict(stations=[{k: v for k, v in r.items()} for r in st_rows], n_traces=len(good), n_total=len(rows),
               centroid_from_egf_km=dict(E=float(dE), N=float(dN), se_E=float(se[1]), se_N=float(se[2])),
               centroid_latlon=[float(cen_lat), float(cen_lon)], tc0=float(sol_c[0]), rms_tc=float(np.sqrt(np.mean(res_c ** 2))),
               centroid_from_csn_hypo_km=dict(E=float(hE), N=float(hN), along_strike_NNE=float(along), across_downdip=float(across)),
               duration_fit=dict(T0=float(sol_T[0]), LE=float(LE), LN=float(LN), dir_deg=float(np.degrees(np.arctan2(LE, LN)) % 360),
                                 L_km=float(np.hypot(LE, LN))),
               gcmt_moment_ratio=GCMT_RATIO)
    (HERE / "astf_results.json").write_text(json.dumps(out, indent=1))

    # ---------- figura ASTF por estación
    n = len(rows); ncol = 6
    fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(14, 2.0 * np.ceil(n / ncol)), sharex=True, constrained_layout=True)
    for a_, r in zip(axes.ravel(), sorted(rows, key=lambda r: (r["az"], r["comp"]))):
        a_.fill_between(r["t"], r["a"], color="k" if r["vr"] > 0.6 else "0.7", alpha=0.6)
        a_.set_title(f"{r['sta']} {r['comp'][-1]} az {r['az']:.0f}° VR {r['vr']:.2f}", fontsize=7.5)
        a_.axvline(r["mom"][1], color="#c4161c", lw=0.8)
        a_.tick_params(labelsize=6); a_.set_yticks([])
    for a_ in axes.ravel()[n:]:
        a_.axis("off")
    fig.supxlabel("time since mainshock origin (s)")
    fig.suptitle("Navidad 2019: apparent source time functions (EGF: 2019-08-01 20:01 Mw 5.9); grey = VR < 0.6; red = centroid time", fontsize=10)
    fig.savefig(HERE / "astf.png", dpi=120); plt.close(fig)

    # ---------- figura directividad
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    azd = np.degrees(az); gg = np.linspace(0, 360, 361); gr = np.radians(gg); pm = p.mean()
    ax[0].plot(azd, tc, "o", color="k")
    for r, x, y in zip(st_rows, azd, tc):
        ax[0].annotate(r["sta"], (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[0].plot(gg, sol_c[0] - pm * (dE * np.sin(gr) + dN * np.cos(gr)), color="#c4161c",
               label=f"fit: centroid {dE:+.1f}±{se[1]:.1f} km E, {dN:+.1f}±{se[2]:.1f} km N of EGF")
    ax[0].set_xlim(azd.min() - 20, azd.max() + 20); ax[0].set_xlabel("station azimuth (°)"); ax[0].set_ylabel("ASTF centroid time (s)")
    ax[0].axvline(STRIKE, color="0.5", ls=":"); ax[0].axvline(STRIKE + 180, color="0.5", ls=":"); ax[0].axvline(STRIKE + 90, color="0.5", ls="--")
    ax[0].legend(fontsize=7.5, frameon=False); ax[0].set_title("Centroid time vs azimuth")
    ax[1].plot(azd, T, "o", color="k")
    ax[1].plot(gg, sol_T[0] - pm * (LE * np.sin(gr) + LN * np.cos(gr)), color="#c4161c",
               label=f"fit: directivity toward {out['duration_fit']['dir_deg']:.0f}°")
    ax[1].set_xlim(azd.min() - 20, azd.max() + 20); ax[1].set_xlabel("station azimuth (°)"); ax[1].set_ylabel("apparent duration 2σ (s)")
    ax[1].legend(fontsize=7.5, frameon=False); ax[1].set_title("Apparent duration vs azimuth (dotted: strike, dashed: downdip)")
    fig.savefig(HERE / "astf_directivity.png", dpi=130)
    print(json.dumps({k: v for k, v in out.items() if k != "stations"}, indent=1))


if __name__ == "__main__":
    main()
