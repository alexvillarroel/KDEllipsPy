"""
plot_secciones.py
=================
Genera UN PDF por modo (aceleracion, velocidad, desplazamiento). Cada PDF
contiene una grilla:

    filas    = estaciones a < DIST_MAX_KM del epicentro (ordenadas por distancia)
    columnas = componentes N / E / Z

El eje X esta alineado al TIEMPO DE ORIGEN del evento (t = 0 s en el origen),
de modo que la linea vertical roja queda en 0 para todas las estaciones y se
pueden comparar las llegadas.

Salida:
    SAC/checks/seccion_acc.pdf
    SAC/checks/seccion_vel.pdf
    SAC/checks/seccion_disp.pdf

Uso:
    python codigos/plot_secciones.py
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
KIND_NOMBRE = {"ACC": "Aceleración", "VEL": "Velocidad", "DISP": "Desplazamiento"}
COMPS = ["N", "E", "Z"]


def listar_estaciones(ev):
    """Devuelve [(estacion, dist_km), ...] de las que estan a < DIST_MAX_KM,
    ordenadas por distancia ascendente."""
    sel = {}
    for f in glob.glob(os.path.join(SAC_DIR["DISP"], "*.sac")):
        tr = obspy.read(f)[0]
        if not tr.stats.channel.endswith("Z"):
            continue
        dist = tr.stats.sac.user0 if USE_DISTANCE == "epicentral" else tr.stats.sac.dist
        if dist < DIST_MAX_KM:
            sel[tr.stats.station] = dist
    return sorted(sel.items(), key=lambda kv: kv[1])


def cargar(station, kind, comp):
    matches = glob.glob(os.path.join(SAC_DIR[kind], f"{station}.*{comp}.*.sac"))
    return obspy.read(matches[0])[0] if matches else None


def pdf_modo(kind, estaciones, ev):
    n = len(estaciones)
    fig, axes = plt.subplots(n, 3, figsize=(11, 1.6 * n), squeeze=False)

    for i, (station, dist) in enumerate(estaciones):
        for j, comp in enumerate(COMPS):
            ax = axes[i, j]
            tr = cargar(station, kind, comp)
            if tr is not None:
                # Tiempo relativo al ORIGEN del evento (origen = 0 s).
                t = tr.times() - (ev["origin_time"] - tr.stats.starttime)
                ax.plot(t, tr.data, "k", lw=0.5)
                ax.axvline(0, color="red", lw=1.0, ls="--")
            else:
                ax.text(0.5, 0.5, f"sin {comp}", ha="center", va="center",
                        transform=ax.transAxes, color="gray")
            ax.grid(alpha=0.3)
            ax.tick_params(labelsize=7)
            if i == 0:
                ax.set_title(f"Componente {comp}", fontsize=11)
            if j == 0:
                ax.set_ylabel(f"{station}\n{dist:.0f} km", fontsize=8)
            if i == n - 1:
                ax.set_xlabel("Tiempo desde el origen [s]", fontsize=9)

    fig.suptitle(f"{KIND_NOMBRE[kind]} ({KIND_UNITS[kind]}) — Calama 2026-05-25  "
                 f"|  estaciones a < {DIST_MAX_KM:.0f} km  "
                 f"(línea roja = tiempo de origen)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.99])

    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, f"seccion_{kind.lower()}.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out


def main():
    ev = load_event(CTL_FILE)
    estaciones = listar_estaciones(ev)
    print(f"{len(estaciones)} estaciones a < {DIST_MAX_KM:.0f} km ({USE_DISTANCE}).")
    print("Estaciones:", ", ".join(f"{s}({d:.0f})" for s, d in estaciones), "\n")
    for kind in ["ACC", "VEL", "DISP"]:
        out = pdf_modo(kind, estaciones, ev)
        print(f"  {KIND_NOMBRE[kind]:14s} -> {os.path.relpath(out, BASE)}")
    print(f"\nListo. PDFs en {os.path.relpath(OUT_DIR, BASE)}/")


if __name__ == "__main__":
    main()
