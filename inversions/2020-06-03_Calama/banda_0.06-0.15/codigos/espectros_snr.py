"""
espectros_snr.py
================
Espectros de SENAL vs RUIDO (signal & noise spectra) para control de calidad
de la banda de inversion EGF (0.1-0.3 Hz).

Pagina 1: mainshock 2020-06-03 (todas las estaciones de registros_20200603).
Pagina 2: EGF 2020-06-16 M4.9 (mseed descargados, respuesta removida).

Metodologia ISOLA (Sokos & Zahradnik): la ventana de ruido (NTW) tiene el
MISMO largo que la de senal (STW) y termina donde esta comienza; el SNR(f)
es el promedio de las curvas S/N por componente, y el valor unico (SNR
average) es el promedio de SNR(f) sobre la banda de inversion.

Por estacion:
  - STW = [tP-2, tP+90] s   (P + S + coda)
  - NTW = mismo largo, inmediatamente antes de la STW. Los txt del CSN
    parten pocos s antes del origen: si no alcanza el largo completo se usa
    lo disponible (>= 8 s, espectros normalizados por duracion) y el SNR se
    marca con '*'; si ni eso hay -> "sin pre-P".
  - SNR(f) se guarda en output/snr_curvas.csv y se grafica (pagina EGF).

Las amplitudes se convierten a desplazamiento (|A|/(2pi f)^2 o |V|/(2pi f))
para que la figura sea comparable con espectros_amplitud.pdf; el SNR es
invariante a esa conversion.

Rojo = estaciones usadas en la inversion EGF.
Salida: output/espectros_snr.pdf
Uso:    python codigos/espectros_snr.py
"""
import os
import sys
import glob
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import obspy
from obspy import UTCDateTime, read_inventory
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel
from scipy.signal.windows import tukey

from io_sismica import read_all

# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TXT_ROOT = os.path.join(BASE, "registros_20200603")
EGF_RAW = os.path.join(BASE, "egf", "20200616T175441", "RAW")
OUT_PDF = os.path.join(BASE, "output", "espectros_snr.pdf")

MAIN = {"time": UTCDateTime("2020-06-03T07:35:34.000"),
        "lat": -23.247, "lon": -68.530, "depth_km": 123.4}
EGF = {"time": UTCDateTime("2020-06-16T17:54:41.000"),
       "lat": -23.307, "lon": -68.578, "depth_km": 108.0}

# Estaciones usadas en la inversion EGF (estaciones_validadas.csv, valida=1).
USADAS = ["PB01", "PB03", "PB05", "PB06", "PB09", "PB19", "AF01"]

BAND = (0.1, 0.3)                  # banda de inversion EGF (sombreada)
WIN_SIG = (-2.0, 90.0)             # STW en s respecto a tP
NOISE_MIN = 8.0                    # s minimos de ruido para calcular SNR
FMIN_PLOT, FMAX_PLOT = 0.02, 20.0  # Hz
SNR_OK, SNR_MARG = 3.0, 1.5        # umbrales para colorear el valor

model = TauPyModel(model="iasp91")


def t_p(dist_deg, depth_km):
    arr = model.get_travel_times(source_depth_in_km=depth_km,
                                 distance_in_degree=dist_deg,
                                 phase_list=["p", "P"])
    return arr[0].time if arr else None


def densidad_espectral(data, fs):
    """|X(f)| * dt / sqrt(T): amplitud por Hz, comparable entre ventanas de
    distinta duracion."""
    d = np.asarray(data, float)
    d = d - d.mean()
    d = d * tukey(len(d), alpha=0.2)
    n = len(d)
    A = np.abs(np.fft.rfft(d)) / fs
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    return f, A / math.sqrt(n / fs)


def suavizar_log(f, A, nf=240, half_width=1.15):
    """Promedio en ventanas multiplicativas sobre grilla log comun."""
    fg = np.logspace(math.log10(FMIN_PLOT), math.log10(FMAX_PLOT), nf)
    out = np.full(nf, np.nan)
    for i, fc in enumerate(fg):
        m = (f >= fc / half_width) & (f <= fc * half_width)
        if m.any():
            out[i] = A[m].mean()
    return fg, out


def a_desplazamiento(fg, A, unidad):
    """acel -> despl (/(2pi f)^2) o vel -> despl (/(2pi f))."""
    w = 2.0 * math.pi * fg
    return A / w ** 2 if unidad == "acc" else A / w


def espectro_estacion(comps, fs, t1_sig, t2_sig, t0_data, unidad):
    """comps: lista de Trace. ISOLA: NTW de igual largo que la STW, terminando
    en su inicio; si la traza no alcanza se usa lo disponible (flag corto).

    Devuelve (fg, S(f), N(f), snr_curve, corto): S y N espectros vectoriales
    en desplazamiento (solo display); snr_curve = promedio por componente de
    S_c(f)/N_c(f), invariante a la conversion."""
    L = t2_sig - t1_sig
    t2_noi = t1_sig
    t1_noi = max(t1_sig - L, t0_data + 1.0)
    if t2_noi - t1_noi < NOISE_MIN:
        t1_noi = None                          # sin pre-evento utilizable
    corto = t1_noi is not None and (t2_noi - t1_noi) < 0.95 * L

    sig2 = noi2 = snr_sum = None
    n_comp = 0
    fg = None
    for tr in comps:
        Ss = Ns = None
        for (t1, t2, dest) in [(t1_sig, t2_sig, "s"), (t1_noi, t2_noi, "n")]:
            if t1 is None:
                continue
            w = tr.copy().trim(t1, t2, nearest_sample=True)
            if len(w) < 16:
                continue
            w.detrend("linear")
            f, A = densidad_espectral(w.data, fs)
            fg, As = suavizar_log(f, A)
            if dest == "s":
                Ss = As
            else:
                Ns = As
        if Ss is None:
            continue
        sig2 = Ss ** 2 if sig2 is None else sig2 + Ss ** 2
        if Ns is not None:
            noi2 = Ns ** 2 if noi2 is None else noi2 + Ns ** 2
            r = np.where(Ns > 0, Ss / np.maximum(Ns, 1e-30), np.nan)
            snr_sum = r if snr_sum is None else snr_sum + r
            n_comp += 1
    if sig2 is None:
        return None
    S = np.sqrt(sig2)
    N = np.sqrt(noi2) if noi2 is not None else None
    snr_curve = snr_sum / n_comp if n_comp else None
    return fg, a_desplazamiento(fg, S, unidad), \
        (a_desplazamiento(fg, N, unidad) if N is not None else None), \
        snr_curve, corto


def snr_promedio(fg, snr_curve, band):
    """SNR average (ISOLA): promedio de la curva SNR(f) sobre la banda."""
    if snr_curve is None:
        return np.nan
    m = (fg >= band[0]) & (fg <= band[1]) & np.isfinite(snr_curve)
    return float(snr_curve[m].mean()) if m.any() else np.nan


# ------------------------------------------------------------------ mainshock
def datos_mainshock():
    st = read_all(TXT_ROOT)
    grupos = {}
    for tr in st:
        grupos.setdefault(tr.stats.station, []).append(tr)

    res = {}
    for est in sorted(grupos):
        comps = grupos[est]
        la = comps[0].stats.coordinates["latitude"]
        lo = comps[0].stats.coordinates["longitude"]
        d_m, _, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], la, lo)
        tp = t_p(kilometers2degrees(d_m / 1000.0), MAIN["depth_km"])
        t_pabs = MAIN["time"] + tp
        fs = comps[0].stats.sampling_rate
        t0 = max(tr.stats.starttime for tr in comps)
        r = espectro_estacion(comps, fs, t_pabs + WIN_SIG[0],
                              t_pabs + WIN_SIG[1], t0, "acc")
        if r is not None:
            res[est] = r
    return res


# ------------------------------------------------------------------ EGF
def datos_egf():
    inv = read_inventory(os.path.join(EGF_RAW, "inventory.xml"))
    res = {}
    for fn in sorted(glob.glob(os.path.join(EGF_RAW, "*.mseed"))):
        est = os.path.basename(fn).split(".")[1]
        st = obspy.read(fn)
        for pref in ("HH", "HL", "HN", "BH"):
            sel = st.select(channel=f"{pref}?")
            if len(sel) >= 3:
                st = sel
                break
        st.detrend("linear")
        st.taper(max_percentage=0.05, type="cosine")
        try:
            st.remove_response(inventory=inv, output="VEL",
                               pre_filt=(0.01, 0.02, 15.0, 18.0))
        except Exception as exc:                                  # noqa: BLE001
            print(f"  [egf] {est}: sin respuesta ({exc})")
            continue
        coo = inv.get_coordinates(f"{st[0].stats.network}.{est}."
                                  f"{st[0].stats.location}.{st[0].stats.channel}",
                                  st[0].stats.starttime)
        d_m, _, _ = gps2dist_azimuth(EGF["lat"], EGF["lon"],
                                     coo["latitude"], coo["longitude"])
        tp = t_p(kilometers2degrees(d_m / 1000.0), EGF["depth_km"])
        t_pabs = EGF["time"] + tp
        fs = st[0].stats.sampling_rate
        t0 = max(tr.stats.starttime for tr in st)
        r = espectro_estacion(list(st), fs, t_pabs + WIN_SIG[0],
                              t_pabs + WIN_SIG[1], t0, "vel")
        if r is not None:
            res[est] = r
    return res


# ------------------------------------------------------------------ figura
def pagina(res, titulo, ncols, pdf, ylab, curva_snr=False):
    nombres = sorted(res)
    nrows = math.ceil(len(nombres) / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 2.4, nrows * 2.1), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, est in zip(axes, nombres):
        fg, S, N, snr_c, corto = res[est]
        snr = snr_promedio(fg, snr_c, BAND)
        usada = est in USADAS
        c_sig = "tab:red" if usada else "0.35"
        ax.loglog(fg, S, color=c_sig, lw=1.3 if usada else 0.9)
        if N is not None:
            ax.loglog(fg, N, color="tab:blue", lw=0.8, ls="--", alpha=0.8)
        ax.axvspan(BAND[0], BAND[1], color="tab:red" if usada else "0.7",
                   alpha=0.12 if usada else 0.10, lw=0)
        if curva_snr and snr_c is not None:
            ax2 = ax.twinx()
            ax2.loglog(fg, snr_c, color="forestgreen", lw=0.8, alpha=0.9)
            ax2.axhline(SNR_OK, color="forestgreen", lw=0.5, ls=":")
            ax2.set_ylim(0.1, 1e3)
            ax2.tick_params(labelsize=5, colors="forestgreen")
        if np.isfinite(snr):
            c = ("forestgreen" if snr >= SNR_OK
                 else "darkorange" if snr >= SNR_MARG else "crimson")
            ax.text(0.97, 0.95, f"SNR={snr:.1f}" + ("*" if corto else ""),
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=7, color=c, fontweight="bold", zorder=5)
        else:
            ax.text(0.97, 0.95, "sin pre-P", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6.5, color="0.5")
        ax.set_title(est, fontsize=9, color="tab:red" if usada else "0.2",
                     fontweight="bold" if usada else "normal")
        ax.set_xlim(FMIN_PLOT, FMAX_PLOT)
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=6)
    for ax in axes[len(nombres):]:
        ax.axis("off")
    fig.supxlabel("Frecuencia (Hz)", fontsize=11)
    fig.supylabel(ylab, fontsize=11)
    fig.suptitle(titulo, fontsize=12)
    fig.tight_layout(rect=[0.02, 0.02, 1, 0.95])
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def main():
    print("Mainshock 2020-06-03 ...")
    rm = datos_mainshock()
    print(f"  {len(rm)} estaciones")
    print("EGF 2020-06-16 M4.9 ...")
    re_ = datos_egf()
    print(f"  {len(re_)} estaciones")

    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
    with PdfPages(OUT_PDF) as pdf:
        pagina(rm, "Senal vs ruido — MAINSHOCK Calama 2020-06-03 (Mw 6.8)\n"
               "negro/rojo = senal (STW: P..coda), azul punteado = ruido "
               "(NTW de igual largo; '*' = NTW recortada por inicio del "
               f"registro) ; SNR = promedio de SNR(f) en {BAND[0]}-{BAND[1]} "
               "Hz (rojo = usadas en inversion EGF)",
               ncols=7, pdf=pdf,
               ylab="Amplitud espectral de desplazamiento (m·s)")
        pagina(re_, "Senal vs ruido — EGF 2020-06-16 M4.9\n"
               "linea = senal (STW), azul punteado = ruido (NTW igual largo), "
               "verde = curva SNR(f) (eje derecho, punteado = umbral 3) ; "
               f"SNR = promedio en {BAND[0]}-{BAND[1]} Hz",
               ncols=4, pdf=pdf, curva_snr=True,
               ylab="Amplitud espectral de desplazamiento (m·s)")
    print(f"-> {OUT_PDF}")

    # curvas SNR(f) a archivo (formato largo), como pide la metodologia
    csvfn = os.path.join(BASE, "output", "snr_curvas.csv")
    with open(csvfn, "w") as f:
        f.write("dataset,sta,freq_hz,snr\n")
        for lbl, res in [("main", rm), ("egf", re_)]:
            for est in sorted(res):
                fg, _, _, snr_c, _ = res[est]
                if snr_c is None:
                    continue
                for fv, sv in zip(fg, snr_c):
                    if np.isfinite(sv):
                        f.write(f"{lbl},{est},{fv:.4f},{sv:.3f}\n")
    print(f"-> {csvfn}")

    print(f"\nSNR average (ISOLA) en {BAND} y en (0.15, 0.4):")
    for lbl, res in [("MAIN", rm), ("EGF ", re_)]:
        for est in sorted(res):
            fg, _, _, snr_c, corto = res[est]
            s1 = snr_promedio(fg, snr_c, BAND)
            s2 = snr_promedio(fg, snr_c, (0.15, 0.4))
            tagu = "*" if est in USADAS else " "
            tc = "(NTW corta)" if corto else ""
            if np.isfinite(s1):
                print(f"  {lbl} {est:5s}{tagu} {s1:7.1f} {s2:7.1f}  {tc}")
            else:
                print(f"  {lbl} {est:5s}{tagu}     n/a     n/a  {tc}")


if __name__ == "__main__":
    main()
