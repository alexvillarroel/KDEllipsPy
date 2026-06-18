"""
validar_egf.py
==============
Valida un candidato EGF contra el mainshock de Calama 2020 en la banda
objetivo 0.1-0.3 Hz:

  1. SNR del EGF por traza (ventana P+S vs ruido pre-evento).
  2. Correlacion cruzada normalizada mainshock <-> EGF por estacion/componente
     (Z en ventana P; transversal en ventana S). CC alta => comparten
     trayectoria + sitio + patron de radiacion => EGF valida.

Mainshock: SAC locales de aceleracion (calama2020/DATA/RAW, fisica m/s^2),
integrados a velocidad. EGF: mseed descargado, respuesta removida a VEL.
Ambos con el MISMO bandpass causal 0.1-0.3 y remuestreo comun a 20 Hz.

Salidas: egf/INFORME_egf.md, egf/validacion_<tag>.png
Uso:     python codigos/validar_egf.py [tag]    (default: top-1 del ranking)
"""
import os
import sys
import csv
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import obspy
from obspy import UTCDateTime, read_inventory
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))      # .../egf
RANKING = os.path.join(BASE, "candidatos_ranking.csv")
MAIN_RAW = os.path.normpath(os.path.join(BASE, "..", "..", "calama2020", "DATA", "RAW"))

MAIN = {"time": UTCDateTime("2020-06-03T07:35:34.000"),
        "lat": -23.247, "lon": -68.530, "depth_km": 123.4}

# Coordenadas de las estaciones comunes (de calama2020/input.ctl).
STA_COORDS = {
    "PB05": (-22.8528, -70.2023),
    "PB06": (-22.7058, -69.5717),
    "PB08": (-20.1411, -69.1535),
    "PB19": (-23.9048, -69.2900),
}

F1, F2 = 0.1, 0.3        # banda objetivo (acordada)
FS = 20.0                # Hz, remuestreo comun
WIN_P, WIN_S = 25.0, 40.0  # s, ventanas de analisis
MAX_LAG = 15.0           # s, busqueda de CC

model = TauPyModel(model="iasp91")


def tt(dist_deg, depth_km):
    """(tP, tS) en segundos."""
    tp = ts = None
    for arr in model.get_travel_times(source_depth_in_km=depth_km,
                                      distance_in_degree=dist_deg,
                                      phase_list=["p", "P", "s", "S"]):
        if arr.name.lower() == "p" and tp is None:
            tp = arr.time
        if arr.name.lower() == "s" and ts is None:
            ts = arr.time
    return tp, ts


def preparar(tr, t0, integrar):
    """Detrend, (integra), filtra causal y remuestrea; recorta a t0-60..t0+360."""
    tr = tr.copy()
    tr.detrend("linear")
    tr.detrend("demean")
    tr.taper(max_percentage=0.05, type="cosine")
    if integrar:
        tr.integrate(method="cumtrapz")
        tr.detrend("linear")
    tr.filter("bandpass", freqmin=F1, freqmax=F2, corners=4, zerophase=False)
    # padding extra para que interpolate no necesite extrapolar en los bordes
    tr.trim(t0 - 62.0, t0 + 362.0, pad=True, fill_value=0.0)
    n_out = int((360.0 + 60.0) * FS)
    tr.interpolate(sampling_rate=FS, starttime=t0 - 60.0, npts=n_out)
    return tr


def rotar(n, e, az_deg):
    """N,E -> R,T con azimut evento->estacion."""
    az = np.deg2rad(az_deg)
    r = n * np.cos(az) + e * np.sin(az)
    t = e * np.cos(az) - n * np.sin(az)
    return r, t


def cc_max(a, b, fs, max_lag_s):
    """CC normalizada maxima y su lag, en ventana ya recortada (mismo largo)."""
    n = min(len(a), len(b))
    a, b = a[:n] - a[:n].mean(), b[:n] - b[:n].mean()
    den = np.sqrt(np.sum(a ** 2) * np.sum(b ** 2))
    if den == 0:
        return 0.0, 0.0
    full = np.correlate(a, b, mode="full") / den
    lags = (np.arange(len(full)) - (n - 1)) / fs
    sel = np.abs(lags) <= max_lag_s
    k = np.argmax(full[sel])
    return float(full[sel][k]), float(lags[sel][k])


def cargar_egf(tag, ev, station):
    raw = os.path.join(BASE, tag, "RAW")
    fn = None
    for cand in glob.glob(os.path.join(raw, f"*.{station}.mseed")):
        fn = cand
        break
    if fn is None:
        return None
    st = obspy.read(fn)
    inv = read_inventory(os.path.join(raw, "inventory.xml"))
    # Prioriza canales HH (banda ancha); si no, HL/HN.
    for pref in ("HH", "HL", "HN", "BH"):
        sel = st.select(channel=f"{pref}?")
        if len(sel) >= 3:
            st = sel
            break
    try:
        st.remove_response(inventory=inv, output="VEL",
                           pre_filt=(0.02, 0.05, 8.0, 9.0))
    except Exception:
        return None
    comps = {}
    for tr in st:
        comps[tr.stats.channel[-1].replace("1", "N").replace("2", "E")] = \
            preparar(tr, ev["time"], integrar=False)
    if not all(c in comps for c in "NEZ"):
        return None
    return comps


def cargar_main(station):
    comps = {}
    for c in "NEZ":
        m = glob.glob(os.path.join(MAIN_RAW, f"*{station}*H?{c}*.sac"))
        if not m:
            return None
        tr = obspy.read(m[0])[0]
        tr.data = np.asarray(tr.data, float)
        comps[c] = preparar(tr, MAIN["time"], integrar=True)   # accel -> vel
    return comps


def main():
    with open(RANKING) as f:
        rows = list(csv.DictReader(f))
    tag = sys.argv[1] if len(sys.argv) > 1 else \
        UTCDateTime(rows[0]["time"]).strftime("%Y%m%dT%H%M%S")
    row = next(r for r in rows
               if UTCDateTime(r["time"]).strftime("%Y%m%dT%H%M%S") == tag)
    ev = {"time": UTCDateTime(row["time"]), "lat": float(row["lat"]),
          "lon": float(row["lon"]), "depth_km": float(row["depth_km"])}
    print(f"Validando EGF {tag} (M{row['mag']}, d3D={row['dist3d_km']} km), "
          f"banda {F1}-{F2} Hz")

    resultados = []
    figs = []
    for sta, (la, lo) in STA_COORDS.items():
        egf = cargar_egf(tag, ev, sta)
        mn = cargar_main(sta)
        if egf is None or mn is None:
            print(f"  {sta}: faltan datos (egf={egf is not None}, main={mn is not None})")
            continue

        d_m, az, _ = gps2dist_azimuth(MAIN["lat"], MAIN["lon"], la, lo)
        deg = kilometers2degrees(d_m / 1000.0)
        tp_m, ts_m = tt(deg, MAIN["depth_km"])
        d_m2, az2, _ = gps2dist_azimuth(ev["lat"], ev["lon"], la, lo)
        deg2 = kilometers2degrees(d_m2 / 1000.0)
        tp_e, ts_e = tt(deg2, ev["depth_km"])

        t_axis = np.arange(len(mn["Z"].data)) / FS - 60.0   # rel. al origen
        def win(tr, t_on, dur):
            i0 = int((t_on + 60.0) * FS)
            return tr.data[max(0, i0):i0 + int(dur * FS)]

        # SNR del EGF (Z, ventana P) y (T, ventana S)
        rn_e, te_e = rotar(egf["N"].data, egf["E"].data, az2)
        noise = egf["Z"].data[:int(50 * FS)]
        snr_p = np.sqrt(np.mean(win(egf["Z"], tp_e, WIN_P) ** 2)) / \
            max(np.sqrt(np.mean(noise ** 2)), 1e-20)
        # CC: Z en P, T en S (cada uno con su tP/tS para alinear fases)
        rn_m, te_m = rotar(mn["N"].data, mn["E"].data, az)
        zp_m = win(mn["Z"], tp_m, WIN_P); zp_e = win(egf["Z"], tp_e, WIN_P)
        ccP, lagP = cc_max(zp_m, zp_e, FS, MAX_LAG)
        i0m = int((ts_m + 60.0 - 5.0) * FS); i0e = int((ts_e + 60.0 - 5.0) * FS)
        ts_w = int(WIN_S * FS)
        ccS, lagS = cc_max(te_m[i0m:i0m + ts_w], te_e[i0e:i0e + ts_w], FS, MAX_LAG)

        resultados.append({"sta": sta, "dist_km": round(d_m / 1000, 1),
                           "snrP_egf": round(float(snr_p), 1),
                           "ccP_Z": round(ccP, 2), "lagP": round(lagP, 1),
                           "ccS_T": round(ccS, 2), "lagS": round(lagS, 1)})
        figs.append((sta, t_axis, mn["Z"].data, egf["Z"].data, tp_m, tp_e,
                     te_m, te_e, ts_m, ts_e))
        print(f"  {sta}: SNR_P(egf)={snr_p:.1f}  CC_P(Z)={ccP:.2f} ({lagP:+.1f}s)"
              f"  CC_S(T)={ccS:.2f} ({lagS:+.1f}s)")

    if not resultados:
        print("Sin estaciones validables."); return

    # ---- figura ----
    nsta = len(figs)
    fig, axes = plt.subplots(nsta, 2, figsize=(13, 2.4 * nsta), sharex=True)
    axes = np.atleast_2d(axes)
    for i, (sta, t, zm, ze, tpm, tpe, tm, te_, tsm, tse) in enumerate(figs):
        for j, (m_, e_, ton_m, ton_e, lbl) in enumerate(
                [(zm, ze, tpm, tpe, "Z (P)"), (tm, te_, tsm, tse, "T (S)")]):
            ax = axes[i, j]
            nm = max(np.max(np.abs(m_)), 1e-20)
            ne = max(np.max(np.abs(e_)), 1e-20)
            ax.plot(t, m_ / nm, "k", lw=0.9, label="mainshock")
            # alinear EGF por la diferencia de tiempos de viaje
            ax.plot(t + (ton_m - ton_e), e_ / ne, "r", lw=0.9, alpha=0.8,
                    label="EGF (alineada)")
            ax.axvline(ton_m, color="b", ls=":", lw=0.8)
            ax.set_ylabel(f"{sta}\n{lbl}", fontsize=8)
            ax.set_xlim(ton_m - 20, ton_m + 80)
            ax.tick_params(labelsize=7)
            if i == 0:
                ax.legend(fontsize=7, loc="upper right")
    axes[-1, 0].set_xlabel("t desde origen (s)")
    axes[-1, 1].set_xlabel("t desde origen (s)")
    fig.suptitle(f"Validacion EGF {tag} vs mainshock — banda {F1}-{F2} Hz "
                 "(amplitudes normalizadas)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    figpath = os.path.join(BASE, f"validacion_{tag}.png")
    fig.savefig(figpath, dpi=150)
    print(f"figura -> {figpath}")

    # ---- informe ----
    cc_med = np.median([r["ccS_T"] for r in resultados])
    veredicto = ("VALIDA (CC_S mediana >= 0.6)" if cc_med >= 0.6 else
                 "MARGINAL (0.4 <= CC_S < 0.6)" if cc_med >= 0.4 else
                 "DEBIL (CC_S < 0.4) — probar otro candidato")
    inf = os.path.join(BASE, "INFORME_egf.md")
    with open(inf, "a") as f:
        f.write(f"\n## Candidato {tag}  (M{row['mag']}, dM={row['dM']}, "
                f"d3D={row['dist3d_km']} km)\n\n")
        f.write(f"Banda {F1}-{F2} Hz. Mainshock local (accel->vel) vs EGF "
                f"(respuesta removida a VEL), mismo filtro causal.\n\n")
        f.write("| sta | dist (km) | SNR_P egf | CC_P (Z) | lag | CC_S (T) | lag |\n")
        f.write("|-----|-----------|-----------|----------|-----|----------|-----|\n")
        for r in resultados:
            f.write(f"| {r['sta']} | {r['dist_km']} | {r['snrP_egf']} | "
                    f"{r['ccP_Z']} | {r['lagP']} | {r['ccS_T']} | {r['lagS']} |\n")
        f.write(f"\n**CC_S mediana = {cc_med:.2f} -> {veredicto}**\n")
        f.write(f"\n![validacion](validacion_{tag}.png)\n")
    print(f"informe -> {inf}\nVeredicto: {veredicto}")


if __name__ == "__main__":
    main()
