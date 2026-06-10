#!/usr/bin/env python3
"""
comparar_modelo_velocidad.py
============================
Compara en PROFUNDIDAD el modelo de velocidades 1D usado en la inversion (seccion
9 del input.ctl, identico en 2020 y 2026) contra el modelo de referencia iasp91
(via obspy.taup). Esto importa porque los tiempos tP/tS de las VENTANAS del misfit
salen de iasp91 (build_azi_times_array), pero las sinteticas de AXITRA usan el
modelo 1D del input.ctl: si difieren a la profundidad de la fuente (~114-123 km),
las ventanas quedan mal alineadas.

Dibuja: Vp(z) y Vs(z) de ambos modelos + panel de diferencia (modelo - iasp91),
marcando las profundidades de los dos eventos.

Uso:
    python codigos/comparar_modelo_velocidad.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel

BASE = Path(__file__).resolve().parent.parent          # comparacion_2020_vs_2026/
INVERSIONS = BASE.parent
INPUT_CTL = INVERSIONS / "2020-06-03" / "input.ctl"     # modelo identico al 2026
OUT = BASE / "output" / "modelo_velocidad_vs_iasp91.png"

ZMAX = 150.0   # km (cubre la fuente ~114-123 km)
EVENTOS = {"Calama 2020": 123.4, "Calama 2026": 114.0}   # profundidades (km)


def leer_modelo_input_ctl(path):
    """Lee la seccion 9 (Velocity Model): devuelve depth[km], vp[km/s], vs[km/s].
    La 1a columna es la profundidad del TOPE de cada capa (m); Vp/Vs en m/s."""
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines) and "Number of layers" not in lines[i]:
        i += 1
    n = int(lines[i].split(":")[-1])
    depth, vp, vs = [], [], []
    i += 1
    while i < len(lines) and len(depth) < n:
        s = lines[i].strip()
        if s and not s.startswith("#"):
            c = s.split()
            depth.append(float(c[0]) / 1000.0)
            vp.append(float(c[1]) / 1000.0)
            vs.append(float(c[2]) / 1000.0)
        i += 1
    return np.array(depth), np.array(vp), np.array(vs)


def perfil_step(depth_top, val, zmax):
    """Convierte (tope de capa, valor) en arrays escalonados para graficar; la
    ultima capa se extiende como semiespacio hasta zmax."""
    zz, vv = [], []
    for k in range(len(depth_top)):
        z0 = depth_top[k]
        z1 = depth_top[k + 1] if k + 1 < len(depth_top) else zmax
        zz += [z0, z1]; vv += [val[k], val[k]]
    return np.array(zz), np.array(vv)


def perfil_iasp91(zmax, fase):
    """Perfil iasp91 (km/s) hasta zmax, usando los topes/bases de cada capa."""
    vm = TauPyModel("iasp91").model.s_mod.v_mod
    ly = vm.layers
    key0 = "top_p_velocity" if fase == "p" else "top_s_velocity"
    key1 = "bot_p_velocity" if fase == "p" else "bot_s_velocity"
    zz, vv = [], []
    for L in ly:
        if L["top_depth"] > zmax:
            break
        zz += [L["top_depth"], min(L["bot_depth"], zmax)]
        vv += [L[key0], L[key1]]
    return np.array(zz), np.array(vv)


def iasp91_en(z, zp_i, vp_i):
    """Velocidad iasp91 a profundidad z (km) por interpolacion del perfil."""
    return float(np.interp(z, zp_i, vp_i))


def main():
    zt, vp, vs = leer_modelo_input_ctl(INPUT_CTL)

    zp_m, vp_m = perfil_step(zt, vp, ZMAX)
    zs_m, vs_m = perfil_step(zt, vs, ZMAX)
    zp_i, vp_i = perfil_iasp91(ZMAX, "p")
    zs_i, vs_i = perfil_iasp91(ZMAX, "s")

    fig, ax = plt.subplots(1, 3, figsize=(13, 8), sharey=True)

    # --- Vp ---
    ax[0].plot(vp_m, zp_m, "-", color="tab:blue", lw=2, label="Modelo input.ctl")
    ax[0].plot(vp_i, zp_i, "-", color="tab:red", lw=1.5, label="iasp91")
    ax[0].set_xlabel("Vp [km/s]"); ax[0].set_ylabel("Profundidad [km]")
    ax[0].set_title("Vp")
    ax[0].legend(loc="lower left", fontsize=9)

    # --- Vs ---
    ax[1].plot(vs_m, zs_m, "-", color="tab:blue", lw=2, label="Modelo input.ctl")
    ax[1].plot(vs_i, zs_i, "-", color="tab:red", lw=1.5, label="iasp91")
    ax[1].set_xlabel("Vs [km/s]"); ax[1].set_title("Vs")

    # --- Diferencia (modelo - iasp91) en una grilla fina ---
    zz = np.linspace(0, ZMAX, 600)
    def step_at(z, depth_top, val):
        idx = np.searchsorted(depth_top, z, side="right") - 1
        return val[np.clip(idx, 0, len(val) - 1)]
    dvp = step_at(zz, zt, vp) - np.interp(zz, zp_i, vp_i)
    dvs = step_at(zz, zt, vs) - np.interp(zz, zs_i, vs_i)
    ax[2].plot(dvp, zz, color="tab:blue", lw=2, label="ΔVp")
    ax[2].plot(dvs, zz, color="tab:green", lw=2, label="ΔVs")
    ax[2].axvline(0, color="0.5", lw=0.8)
    ax[2].set_xlabel("Δ (modelo − iasp91) [km/s]")
    ax[2].set_title("Diferencia")
    ax[2].legend(loc="lower left", fontsize=9)

    # marcar profundidades de los eventos en todos los paneles
    for a in ax:
        for name, z in EVENTOS.items():
            a.axhline(z, color="k", ls=":", lw=1.0)
        a.grid(alpha=0.3)
    for name, z in EVENTOS.items():
        ax[2].text(ax[2].get_xlim()[1], z, f" {name} ({z:.0f} km)",
                   va="center", ha="right", fontsize=7,
                   bbox=dict(fc="white", ec="none", alpha=0.6))

    ax[0].set_ylim(ZMAX, 0)   # profundidad creciente hacia abajo
    fig.suptitle("Modelo de velocidades de la inversion vs. iasp91", y=0.98, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"figura -> {OUT}")

    # tabla a las profundidades de los eventos
    print(f"\n{'z [km]':>8}{'Vp_mod':>9}{'Vp_iasp':>9}{'ΔVp':>8}{'Vs_mod':>9}{'Vs_iasp':>9}{'ΔVs':>8}")
    for name, z in EVENTOS.items():
        vpm = step_at(np.array([z]), zt, vp)[0]; vsm = step_at(np.array([z]), zt, vs)[0]
        vpi = iasp91_en(z, zp_i, vp_i); vsi = iasp91_en(z, zs_i, vs_i)
        print(f"{z:>8.1f}{vpm:>9.2f}{vpi:>9.2f}{vpm-vpi:>+8.2f}"
              f"{vsm:>9.2f}{vsi:>9.2f}{vsm-vsi:>+8.2f}  <- {name}")


if __name__ == "__main__":
    main()
