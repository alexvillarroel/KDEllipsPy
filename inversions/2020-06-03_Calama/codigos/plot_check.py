"""
plot_check.py
=============
Genera una figura de revision (estilo ``check_A02F_HNZ.png``) por cada estacion
ubicada a menos de ``DIST_MAX_KM`` del epicentro (distancia EPICENTRAL).

Cada figura tiene una grilla 3x3:
    filas    = aceleracion / velocidad / desplazamiento
    columnas = componentes N / E / Z
y una linea vertical roja en el TIEMPO DE ORIGEN del evento.

Lee los SAC generados por ``txt_to_sac.py`` (SAC/ACC, SAC/VEL, SAC/DISP) y los
guarda en ``SAC/checks/check_<ESTACION>.png``.

Uso:
    python codigos/plot_check.py
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import obspy  # noqa: E402

from event import load_event  # noqa: E402

# ==============================================================================
# CONFIGURACION
# ==============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CTL_FILE = os.path.join(BASE, "event.ctl")
SAC_DIR = {"ACC": os.path.join(BASE, "SAC", "ACC"),
           "VEL": os.path.join(BASE, "SAC", "VEL"),
           "DISP": os.path.join(BASE, "SAC", "DISP")}
OUT_DIR = os.path.join(BASE, "SAC", "checks")

DIST_MAX_KM = 200.0            # umbral de distancia epicentral
USE_DISTANCE = "epicentral"    # "epicentral" (sac.user0) o "hypocentral" (sac.dist)

KIND_UNITS = {"ACC": "m/s²", "VEL": "m/s", "DISP": "m"}
COMPS = ["N", "E", "Z"]


def listar_estaciones(ev):
    """Devuelve {estacion: dist_km} para las que estan a < DIST_MAX_KM."""
    sel = {}
    for f in glob.glob(os.path.join(SAC_DIR["DISP"], "*.sac")):
        tr = obspy.read(f)[0]
        if not tr.stats.channel.endswith("Z"):
            continue
        dist = tr.stats.sac.user0 if USE_DISTANCE == "epicentral" else tr.stats.sac.dist
        if dist < DIST_MAX_KM:
            sel[tr.stats.station] = dist
    return dict(sorted(sel.items(), key=lambda kv: kv[1]))


def cargar(station, kind, comp):
    """Carga el Trace de una estacion/tipo/componente (o None si no existe)."""
    matches = glob.glob(os.path.join(SAC_DIR[kind], f"{station}.*{comp}.*.sac"))
    return obspy.read(matches[0])[0] if matches else None


def figura_estacion(station, dist, ev):
    fig, axes = plt.subplots(3, 3, figsize=(13, 8), sharex=True)
    t_event = ev["origin_time"]                      # tiempo de origen absoluto

    for i, kind in enumerate(["ACC", "VEL", "DISP"]):
        for j, comp in enumerate(COMPS):
            ax = axes[i, j]
            tr = cargar(station, kind, comp)
            if tr is not None:
                # Eje de tiempo relativo al inicio de la traza.
                ax.plot(tr.times(), tr.data, "k", lw=0.6)
                # Linea vertical en el tiempo de origen del evento.
                t0 = t_event - tr.stats.starttime    # s desde inicio de traza
                ax.axvline(t0, color="red", lw=1.2, ls="--",
                           label="Origen" if (i == 0 and j == 0) else None)
            else:
                ax.text(0.5, 0.5, f"sin {comp}", ha="center", va="center",
                        transform=ax.transAxes, color="gray")
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(f"Componente {comp}")
            if j == 0:
                ax.set_ylabel(f"{kind}\n[{KIND_UNITS[kind]}]")
    for ax in axes[-1, :]:
        ax.set_xlabel("Tiempo [s] desde inicio de traza")

    fig.suptitle(f"{station} — Calama 2026-05-25  |  dist epicentral = {dist:.0f} km",
                 fontsize=13)
    # Leyenda con la linea de origen (una sola).
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, f"check_{station}.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def main():
    ev = load_event(CTL_FILE)
    estaciones = listar_estaciones(ev)
    print(f"{len(estaciones)} estaciones a < {DIST_MAX_KM:.0f} km "
          f"({USE_DISTANCE}). Generando figuras...\n")
    for station, dist in estaciones.items():
        out = figura_estacion(station, dist, ev)
        print(f"  {station:6s} {dist:6.1f} km -> {os.path.relpath(out, BASE)}")
    print(f"\nListo. Figuras en {os.path.relpath(OUT_DIR, BASE)}/")


if __name__ == "__main__":
    main()
