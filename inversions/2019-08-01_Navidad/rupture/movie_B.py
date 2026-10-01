"""Película mp4 del dinámico B: Δτ, slip-rate y slip acumulado en el plano de falla,
más la función fuente con un cursor de tiempo. Reutiliza load_B de traction_plots.

Escribe rupture_B.mp4 (0.1 s de ruptura por cuadro, 15 fps -> ~1.5x más lento que real).
"""

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter

from traction_plots import HERE, NAV, decorate, load_B

T_MAX, DT_FRAME, FPS = 18.0, 0.1, 15


def main():
    B = load_B()
    fp = B["cfg"].fault_plane
    inter = json.loads((NAV / "faults" / "pichilemu_intersection.json").read_text())["line_in_mesh"]
    ix, iw = [p["strike_km"] for p in inter], [p["dip_km"] for p in inter]
    x, w, dt = B["x"], B["w"], B["dt"]
    dtau = (B["tau"] - B["t0"][None]) / 1e6
    v = B["v"]
    slip = np.cumsum(v, axis=0) * dt
    stf = v.sum(axis=(1, 2))  # proporcional a la tasa de momento (mu y dA constantes)
    t = np.arange(v.shape[0]) * dt
    i2 = int(2 / dt)  # escalas desde t = 2 s: la nucleación forzada saturaría todo
    lim, vmax, smax = np.percentile(np.abs(dtau[i2:]), 99.5), v[i2:].max(), slip[-1].max()

    fig = plt.figure(figsize=(13, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[3, 1])
    axs = [fig.add_subplot(gs[0, j]) for j in range(3)]
    fields = [(dtau, "RdBu_r", -lim, lim, "Δτ = τ(t) − τ₀ (MPa)"),
              (v, "Oranges", 0, vmax, "slip-rate (m/s)"),
              (slip, "viridis", 0, smax, "slip acumulado (m)")]
    ims = []
    for a, (f, cm, lo, hi, lab) in zip(axs, fields):
        im = a.pcolormesh(x, w, np.clip(f[0].T, lo, hi), cmap=cm, vmin=lo, vmax=hi, shading="auto")
        fig.colorbar(im, ax=a, orientation="horizontal", shrink=0.85, pad=0.02, label=lab)
        decorate(a, fp, ix, iw)
        a.set_xlabel("rumbo (km)")
        ims.append(im)
    axs[0].set_ylabel("manteo desde el borde superior (km)")
    ax = fig.add_subplot(gs[1, :])
    ax.plot(t, stf / stf.max(), color="k", lw=1.5)
    ax.fill_between(t, stf / stf.max(), color="k", alpha=0.08)
    cur = ax.axvline(0, color="#c4161c", lw=1.5)
    ax.set_xlim(0, T_MAX); ax.set_ylim(0, 1.05); ax.set_yticks([])
    ax.set_xlabel("tiempo (s)"); ax.set_ylabel("tasa de momento")
    title = fig.suptitle("")

    step = int(round(DT_FRAME / dt))
    writer = FFMpegWriter(fps=FPS, bitrate=4000, metadata={"title": "Navidad 2019 dinámico B"})
    with writer.saving(fig, str(HERE / "rupture_B.mp4"), dpi=110):
        for i in range(0, int(T_MAX / dt), step):
            for im, (f, _, lo, hi, _) in zip(ims, fields):
                im.set_array(np.clip(f[i].T, lo, hi).ravel())
            cur.set_xdata([t[i]])
            title.set_text(f"Navidad 2019 – ruptura dinámica B   t = {t[i]:5.1f} s   "
                           "(rojo discontinuo: intersección con Pichilemu; escalas saturadas en t < 2 s)")
            writer.grab_frame()
    print("ok", HERE / "rupture_B.mp4")


if __name__ == "__main__":
    main()
