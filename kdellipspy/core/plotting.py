"""Unified plotting functions for kdellipspy.
"""

from __future__ import annotations

from typing import Any, List, Optional, Tuple, TYPE_CHECKING
import math
import json
from pathlib import Path
from csv import DictReader

import numpy as np

if TYPE_CHECKING:
    from .config_parser import ConfigParser
    from .geometry import FaultGeometry

def plot_waveform_fit(
    observed: np.ndarray,
    synthetic: np.ndarray,
    time: np.ndarray,
    station_names: List[str],
    misfit: Optional[float] = None,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300,
    azimuths: Optional[np.ndarray] = None,
    tp_s: Optional[np.ndarray] = None,
    ts_s: Optional[np.ndarray] = None,
    window_s: Optional[float] = None,
    station_flags: Optional[np.ndarray] = None,
    rotate: bool = False,
    mark_windows: bool = False,
) -> Tuple[Any, Any]:
    """
    Plot observed vs synthetic waveforms.

    Parameters
    ----------
    rotate : if True and ``azimuths`` (rad, one per station) is given, rotate the
        N/E components to Radial/Transverse and label columns R / T / Z, matching
        how the L2 misfit is computed (P -> R+Z, S -> T).
    mark_windows : if True and ``tp_s``, ``ts_s`` (s) and ``window_s`` (>0) are
        given, shade the misfit window per component: P window [tP, tP+window_s]
        on R and Z, S window [tS, tS+window_s] on T.
    station_flags : (nsta, 3) booleans (use_n, use_e, use_z). Components that are
        NOT used in the misfit are drawn translucent and tagged "no usada".
    """
    import matplotlib.pyplot as plt

    obs = np.asarray(observed, dtype=float)
    syn = np.asarray(synthetic, dtype=float)
    nsta = obs.shape[0]

    do_rotate = bool(rotate) and azimuths is not None
    if do_rotate:
        az = np.asarray(azimuths, dtype=float)
        c = np.cos(az)[:, None]
        s = np.sin(az)[:, None]
        R_o = obs[:, 0] * c + obs[:, 1] * s
        T_o = obs[:, 1] * c - obs[:, 0] * s
        R_s = syn[:, 0] * c + syn[:, 1] * s
        T_s = syn[:, 1] * c - syn[:, 0] * s
        plot_o = [R_o, T_o, obs[:, 2]]
        plot_s = [R_s, T_s, syn[:, 2]]
        comps = ['R (P)', 'T (S)', 'Z (P)']
    else:
        plot_o = [obs[:, 0], obs[:, 1], obs[:, 2]]
        plot_s = [syn[:, 0], syn[:, 1], syn[:, 2]]
        comps = ['North', 'East', 'Z']

    do_win = (
        bool(mark_windows) and tp_s is not None and ts_s is not None
        and window_s is not None and float(window_s) > 0.0
    )
    if do_win:
        tp = np.asarray(tp_s, dtype=float)
        ts = np.asarray(ts_s, dtype=float)
        win_start = [tp, ts, tp]   # R -> P, T -> S, Z -> P

    flags = np.asarray(station_flags, dtype=bool) if station_flags is not None else None

    fig, axs = plt.subplots(nsta, 3, figsize=(11, 1.6 * nsta), squeeze=False, sharex=True)

    title = "Waveform Fit (R/T/Z)" if do_rotate else "Waveform Fit"
    if do_win:
        title += f" · ventana {float(window_s):.0f}s (P→R,Z ; S→T)"
    if misfit is not None:
        title += f"  (Misfit: {misfit:.4f})"
    fig.suptitle(title, fontsize=14, y=1.01)

    for i in range(nsta):
        for j in range(3):
            ax = axs[i, j]
            used = True if flags is None else bool(flags[i, j])
            a_obs = 1.0 if used else 0.30
            a_syn = 1.0 if used else 0.40
            ax.plot(time, plot_o[j][i], color='black', linewidth=1.4, alpha=a_obs,
                    label='Observed' if (i == 0 and j == 0) else None)
            ax.plot(time, plot_s[j][i], color='red', linewidth=1.8, linestyle='--',
                    alpha=a_syn, label='Synthetic' if (i == 0 and j == 0) else None)

            if do_win:
                t0 = float(win_start[j][i])
                ax.axvspan(t0, t0 + float(window_s), color='orange', alpha=0.15, lw=0)
                ax.axvline(t0, color='orange', lw=0.7, alpha=0.6)
            if not used:
                ax.text(0.97, 0.92, 'no usada', transform=ax.transAxes, ha='right',
                        va='top', fontsize=7, color='gray', style='italic')

            if i == 0:
                ax.set_title(comps[j])
            if j == 0:
                ax.set_ylabel(station_names[i], fontweight='bold')
            if i == nsta - 1:
                ax.set_xlabel("Time (s)")

            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelleft=False)

    fig.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axs

def plot_azimuthal_coverage(
    station_lats,
    station_lons,
    station_names,
    ev_lat: float,
    ev_lon: float,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 200,
) -> Tuple[Any, Any]:
    """Diagrama polar (radar) de la cobertura azimutal de las estaciones respecto
    al epicentro. Ángulo = azimut evento→estación (N arriba, horario), radio =
    distancia epicentral (km). Sombrea el mayor hueco azimutal (gap)."""
    import matplotlib.pyplot as plt
    from obspy.geodetics import gps2dist_azimuth

    az_deg, dist_km = [], []
    for la, lo in zip(station_lats, station_lons):
        d_m, az, _baz = gps2dist_azimuth(ev_lat, ev_lon, la, lo)
        az_deg.append(az)
        dist_km.append(d_m / 1000.0)
    az_deg = np.array(az_deg)
    dist_km = np.array(dist_km)

    a_sorted = np.sort(az_deg)
    gaps = np.diff(np.concatenate([a_sorted, [a_sorted[0] + 360.0]]))
    gi = int(np.argmax(gaps))
    gap_max = float(gaps[gi]); gap_lo = a_sorted[gi]; gap_hi = gap_lo + gap_max

    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    rmax = dist_km.max() * 1.15
    ax.set_ylim(0, rmax)
    ax.fill_between(np.deg2rad(np.linspace(gap_lo, gap_hi, 50)), 0, rmax,
                    color="red", alpha=0.12, label=f"Gap máx = {gap_max:.0f}°")
    for az, d, nm in zip(az_deg, dist_km, station_names):
        t = np.deg2rad(az)
        ax.plot([t, t], [0, d], color="0.5", lw=1.0, zorder=2)
        ax.plot(t, d, marker="^", ms=12, color="dodgerblue", markeredgecolor="k", zorder=3)
        ax.text(t, d + rmax * 0.05, f"{nm}\n{az:.0f}° · {d:.0f}km", ha="center",
                va="center", fontsize=8, zorder=4)
    ax.plot(0, 0, marker="*", ms=14, color="yellow", markeredgecolor="k", zorder=5)
    ax.set_rlabel_position(135)
    ax.set_title(f"Cobertura azimutal — {len(az_deg)} estaciones\n"
                 f"gap máximo {gap_max:.0f}° ({gap_lo:.0f}°–{gap_hi:.0f}°)",
                 fontsize=12, pad=20)
    ax.legend(loc="lower right", bbox_to_anchor=(1.15, -0.05), fontsize=9)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def plot_ellipse_map(
    lats, lons, slip,
    ev_lat: float, ev_lon: float,
    station_lats=None, station_lons=None, station_names=None,
    pad: float = 0.25,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 200,
) -> Tuple[Any, Any]:
    """Mapa cartopy de la elipse de slip proyectada a la superficie: footprint
    (subfallas coloreadas por slip) + contorno del borde a 0.15*max_slip,
    hipocentro/epicentro y estaciones dentro del recuadro."""
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    lats = np.asarray(lats, float); lons = np.asarray(lons, float)
    slip = np.asarray(slip, float)
    max_slip = float(slip.max()) if slip.size else 0.0

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=proj)
    if max_slip > 0:
        m = slip > 0.05 * max_slip
        west, east = lons[m].min() - pad, lons[m].max() + pad
        south, north = lats[m].min() - pad, lats[m].max() + pad
    else:
        west, east = ev_lon - pad, ev_lon + pad
        south, north = ev_lat - pad, ev_lat + pad
    ax.set_extent([west, east, south, north], crs=proj)

    for feat, kw in [(cfeature.LAND, dict(facecolor="0.93")),
                     (cfeature.BORDERS, dict(edgecolor="0.4", linewidth=0.8)),
                     (cfeature.COASTLINE, dict(edgecolor="0.3", linewidth=0.8))]:
        try:
            ax.add_feature(feat, **kw)
        except Exception:
            pass

    if max_slip > 0:
        levels = np.linspace(0.15 * max_slip, max_slip, 12)
        tcf = ax.tricontourf(lons, lats, slip, levels=levels, cmap="hot_r",
                             alpha=0.85, transform=proj, extend="max")
        fig.colorbar(tcf, ax=ax, shrink=0.6, pad=0.02).set_label("Slip (m)")
        ax.tricontour(lons, lats, slip, levels=[0.15 * max_slip], colors="lime",
                      linewidths=2.0, transform=proj)

    ax.plot(ev_lon, ev_lat, marker="*", ms=11, color="yellow", markeredgecolor="k",
            markeredgewidth=0.8, transform=proj, zorder=6, label="Epicentro")
    if station_lats is not None:
        for la, lo, nm in zip(station_lats, station_lons, station_names):
            if west <= lo <= east and south <= la <= north:
                ax.plot(lo, la, marker="^", ms=10, color="dodgerblue",
                        markeredgecolor="k", transform=proj, zorder=6)
                ax.text(lo, la + (north - south) * 0.02, nm, fontsize=8,
                        ha="center", transform=proj, zorder=7)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = gl.right_labels = False
    ax.set_title("Elipse de slip (proyectada a superficie)", fontsize=12)
    ax.legend(loc="upper right", fontsize=8)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def plot_misfit_contribution(
    num: np.ndarray,
    den: np.ndarray,
    station_names: List[str],
    total_E: float,
    mode_label: str = "",
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 200,
) -> Tuple[Any, Any]:
    """Visualiza la 'tabla' de contribución al misfit como dos heatmaps
    estaciones×(R,T,Z): (izq) %% del residuo total que aporta cada componente;
    (der) misfit local num/den (calidad de ajuste). Anotados."""
    import matplotlib.pyplot as plt

    num = np.asarray(num, float); den = np.asarray(den, float)
    num_tot = num.sum()
    pct = 100.0 * num / num_tot if num_tot > 0 else num * 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        mfl = np.where(den > 0, num / den, 0.0)
    # Ordenar estaciones por contribución total (descendente).
    order = np.argsort(num.sum(axis=1))[::-1]
    pct, mfl = pct[order], mfl[order]
    names = [station_names[i] for i in order]
    comps = ["R", "T", "Z"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 0.5 * len(names) + 2.5))
    for ax, data, ttl, cmap, vmax in [
        (axes[0], pct, "% del residuo total", "magma_r", None),
        (axes[1], mfl, "misfit local (num/den)", "RdYlGn_r", 1.0),
    ]:
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=0,
                       vmax=(np.nanmax(data) if vmax is None else vmax))
        ax.set_xticks(range(3)); ax.set_xticklabels(comps)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
        for i in range(len(names)):
            for j in range(3):
                ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center",
                        fontsize=8, color="0.15")
        ax.set_title(ttl, fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    fig.suptitle(f"Contribución al misfit por estación/componente "
                 f"({mode_label})  ·  E = {total_E:.4f}", fontsize=13, y=1.02)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes


def plot_synthetic_components(
    time: np.ndarray,
    synthetic: np.ndarray,
    station_names: List[str],
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot purely synthetic seismograms.
    """
    import matplotlib.pyplot as plt
    
    nsta = synthetic.shape[0]
    fig, axs = plt.subplots(nsta, 3, figsize=(12, 2 * nsta), squeeze=False, sharex=True)
    comps = ['North (x)', 'East (y)', 'Vertical (z)']
    colors = ['red', 'blue', 'green']
    
    for i in range(nsta):
        for j in range(3):
            ax = axs[i, j]
            ax.plot(time, synthetic[i, j, :], color=colors[j], linewidth=1.0)
            
            if i == 0:
                ax.set_title(comps[j])
            if j == 0:
                ax.set_ylabel(f"{station_names[i]}\nAmp", fontweight='bold')
            if i == nsta - 1:
                ax.set_xlabel("Time (s)")
            
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axs

def plot_record_section(
    waveforms: np.ndarray,
    time: np.ndarray,
    station_names: List[str],
    distances: np.ndarray,
    comp_name: str = "Z",
    scale: float = 1.0,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot waveforms as a record section.
    """
    import matplotlib.pyplot as plt
    
    nsta = waveforms.shape[0]
    fig, ax = plt.subplots(figsize=(10, 8))
    
    max_amp = np.max(np.abs(waveforms))
    norm_factor = (np.max(distances) - np.min(distances)) / (nsta * 2) if nsta > 1 else 1.0
    if max_amp == 0: max_amp = 1.0

    for i in range(nsta):
        trace = waveforms[i, :] / max_amp * norm_factor * scale
        ax.plot(time, trace + distances[i], color='black', lw=0.7)
        ax.text(time[-1], distances[i], f" {station_names[i]} ({distances[i]:.1f} km)", 
                va='center', fontsize=9)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Distance (km)")
    ax.set_title(f"Record Section - Component {comp_name}")
    ax.grid(True, alpha=0.3)
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax

def plot_slip_distribution(
    geometry: FaultGeometry,
    title: str = "2D Slip Distribution",
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot 2D slip distribution on the fault-plane subfault grid (no interpolation).
    """
    import matplotlib.pyplot as plt

    # Slip por subfalla: maximo de los source points que caen en ella.
    slip_map = {}
    for sp in geometry.source_points:
        sf_idx = sp.subfault_index
        slip_map[sf_idx] = max(slip_map.get(sf_idx, 0.0), sp.displacement)

    nx, ny = int(geometry.nx), int(geometry.ny)
    dstk = geometry.length_strike_m / nx / 1000.0   # km por subfalla (strike)
    ddip = geometry.length_dip_m / ny / 1000.0       # km por subfalla (dip)
    hx = geometry.hypo_strike_m / 1000.0
    hy = geometry.hypo_dip_m / 1000.0

    # Malla rectangular real en coords de falla (along-strike, along-dip),
    # relativa al hipocentro -> hipocentro en (0,0). Cada celda = una subfalla.
    Z = np.zeros((ny, nx))
    for sf in geometry.subfaults:
        istk = (int(sf.index) - 1) % nx
        idip = (int(sf.index) - 1) // nx
        Z[idip, istk] = max(0.0, slip_map.get(sf.index, 0.0))

    strike = (np.arange(nx) + 0.5) * dstk - hx
    dip = (np.arange(ny) + 0.5) * ddip - hy
    X, Y = np.meshgrid(strike, dip)

    max_slip = float(np.max(Z))
    vmax = max(1.0, np.ceil(max_slip * 10.0) / 10.0)   # colorbar fijo 0 -> 1+

    fig, ax = plt.subplots(figsize=(10, 6))
    mapa_calor = ax.pcolormesh(X, Y, Z, cmap='gnuplot', shading='nearest',
                               vmin=0.0, vmax=vmax)

    # Contorno del borde de la elipse al 15% del slip maximo.
    if max_slip > 0:
        ax.contour(X, Y, Z, levels=[0.15 * max_slip], colors='lime',
                   linestyles='dashed', linewidths=1.5)

    ax.scatter(0, 0, marker='*', s=350, color='yellow', edgecolors='black', linewidth=1.2, label='Hypocenter', zorder=10)

    cbar = fig.colorbar(mapa_calor, ax=ax)
    cbar.set_label("Slip (m)", fontsize=10)

    ax.set_xlabel("Along Strike (km)", fontsize=11)
    ax.set_ylabel("Along Dip (km)", fontsize=11)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.legend(loc='upper right', frameon=True, shadow=True)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax

def plot_stations_map(
    cfg: ConfigParser,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300,
    use_cartopy: bool = True
) -> Tuple[Any, Any]:
    """
    Plot station locations on a map using Cartopy for geographic features if available.
    """
    import matplotlib.pyplot as plt
    
    stations = cfg.stations.stations
    src_lat = cfg.source_position.latitude
    src_lon = cfg.source_position.longitude

    lons = [s.longitude for s in stations] + [src_lon]
    lats = [s.latitude for s in stations] + [src_lat]
    
    # Calculate map extent with some padding
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)
    pad = 0.5
    extent = [lon_min - pad, lon_max + pad, lat_min - pad, lat_max + pad]

    ax = None
    fig = None
    cartopy_active = False

    if use_cartopy:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
            
            fig = plt.figure(figsize=(10, 8))
            projection = ccrs.PlateCarree()
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_extent(extent, crs=projection)

            # Add geographical features
            ax.add_feature(cfeature.LAND, facecolor='whitesmoke')
            ax.add_feature(cfeature.OCEAN, facecolor='lightsteelblue')
            ax.add_feature(cfeature.COASTLINE, edgecolor='gray', linewidth=0.8)
            ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
            ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='silver')

            gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                              alpha=0.3, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            
            cartopy_active = True
        except ImportError:
            print("Cartopy not found. Falling back to basic matplotlib plot.")

    if not cartopy_active:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.3)

    # Plot stations
    for i, s in enumerate(stations):
        ax.scatter(s.longitude, s.latitude, marker='^', s=120, color='blue', 
                   edgecolor='black', label='Station' if i==0 else "", zorder=5,
                   transform=ccrs.PlateCarree() if cartopy_active else None)
        ax.text(s.longitude, s.latitude, f"  {s.name}", verticalalignment='bottom', 
                fontsize=10, fontweight='bold', zorder=6,
                transform=ccrs.PlateCarree() if cartopy_active else None)
        
        comps = []
        if getattr(s, 'use_n', True): comps.append('N')
        if getattr(s, 'use_e', True): comps.append('E')
        if getattr(s, 'use_z', True): comps.append('Z')
        comp_str = f"({','.join(comps)})"
        ax.text(s.longitude, s.latitude, f"  {comp_str}", verticalalignment='top', 
                fontsize=8, color='dimgray', zorder=6,
                transform=ccrs.PlateCarree() if cartopy_active else None)

    # Plot hypocenter
    ax.scatter(src_lon, src_lat, marker='*', s=300, color='red', edgecolor='black', 
               label='Hypocenter', zorder=10,
               transform=ccrs.PlateCarree() if cartopy_active else None)
    
    ax.set_title(f"Station Distribution - {cfg.source_position.event_name}", 
                 fontsize=14, pad=20, fontweight='bold')
    ax.legend(loc='lower right', frameon=True, shadow=True)
    
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax

def plot_na_results(
    source_data: list[dict[str, float]],
    param_names: list[str],
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot summary of NA results.
    """
    import matplotlib.pyplot as plt
    
    misfits = np.array([row["misfit"] for row in source_data], dtype=float)
    iterations = np.array([row["iteration"] for row in source_data], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    axes[0, 0].semilogy(np.arange(len(misfits)), misfits, color="tab:red", lw=1.2)
    axes[0, 0].set_title("Misfit evolution")
    axes[0, 0].set_xlabel("Model index")
    axes[0, 0].set_ylabel("Misfit")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].scatter(iterations, misfits, s=14, alpha=0.7, color="tab:blue")
    axes[0, 1].set_title("Misfit by iteration")
    axes[0, 1].set_xlabel("Iteration")
    axes[0, 1].set_ylabel("Misfit")
    axes[0, 1].set_yscale("log")
    axes[0, 1].grid(True, alpha=0.3)

    best_idx = int(np.argmin(misfits))
    best_row = source_data[best_idx]
    if param_names:
        axes[1, 0].bar(param_names, [best_row.get(name, np.nan) for name in param_names], color="tab:green")
        axes[1, 0].tick_params(axis="x", rotation=35)
        axes[1, 0].set_title("Best model parameters")
        axes[1, 0].set_ylabel("Value")
        axes[1, 0].grid(True, axis="y", alpha=0.3)
    else:
        axes[1, 0].axis("off")

    if param_names:
        p0 = param_names[0]
        p1 = param_names[1] if len(param_names) > 1 else param_names[0]
        sc = axes[1, 1].scatter(
            [row.get(p0, np.nan) for row in source_data],
            [row.get(p1, np.nan) for row in source_data],
            c=misfits,
            cmap="viridis",
            vmin=0.0, vmax=1.0,          # misfit fijo 0 -> 1 (>1 se recorta al tope)
            s=20,
        )
        axes[1, 1].set_xlabel(p0)
        axes[1, 1].set_ylabel(p1)
        axes[1, 1].set_title("Parameter cloud colored by misfit")
        axes[1, 1].grid(True, alpha=0.3)
        cbar = fig.colorbar(sc, ax=axes[1, 1], extend="max")
        cbar.set_label("Misfit")
    else:
        axes[1, 1].axis("off")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes

def plot_parameter_convergence(
    rows: list[dict],
    param_names: list[str],
    misfits: np.ndarray,
    iterations: np.ndarray,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot convergence for each parameter.
    """
    import matplotlib.pyplot as plt
    
    if not param_names:
        return None, None

    n_params = len(param_names)
    ncols = 4
    nrows = math.ceil(n_params / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.0), squeeze=False)

    model_indices = np.arange(len(rows))
    best_idx = int(np.argmin(misfits))

    norm = plt.Normalize(vmin=0.0, vmax=1.0)   # misfit fijo 0 -> 1 (>1 se recorta al tope)
    cmap = plt.cm.get_cmap("plasma_r")

    for k, name in enumerate(param_names):
        row_i, col_i = divmod(k, ncols)
        ax = axes[row_i][col_i]

        values = np.array([row.get(name, np.nan) for row in rows], dtype=float)
        best_val = float(rows[best_idx].get(name, np.nan))

        ax.scatter(model_indices, values, c=misfits, cmap=cmap, norm=norm, s=12, alpha=0.75, zorder=2)
        ax.axhline(best_val, color="crimson", lw=1.2, ls="--", zorder=3, label=f"best = {best_val:.3g}")

        iter_changes = np.where(np.diff(iterations.astype(int)) > 0)[0] + 1
        for xc in iter_changes:
            ax.axvline(xc, color="gray", lw=0.5, alpha=0.4, zorder=1)

        ax.set_title(name, fontsize=9, pad=3)
        ax.set_xlabel("Model index", fontsize=7)
        ax.set_ylabel("Value", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper right", framealpha=0.6)
        ax.grid(True, alpha=0.2)

    for k in range(n_params, nrows * ncols):
        row_i, col_i = divmod(k, ncols)
        axes[row_i][col_i].axis("off")

    fig.subplots_adjust(right=0.88, hspace=0.55, wspace=0.4)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, extend="max")
    cbar.set_label("Misfit", fontsize=9)

    fig.suptitle("Parameter convergence — NA search", fontsize=11, y=1.01)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes


def plot_uncertainty_corner(
    samples: np.ndarray,
    param_names: List[str],
    bounds: Optional[np.ndarray] = None,
    truths: Optional[np.ndarray] = None,
    mean: Optional[np.ndarray] = None,
    bins: int = 40,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    dpi: int = 300,
    title: Optional[str] = None,
) -> Tuple[Any, Any]:
    """Corner plot of the NA appraisal posterior (1D marginals + 2D trade-offs).
    (Gráfico 'corner' de la posterior del appraisal NA: marginales 1D en la diagonal
    y densidades 2D —trade-offs entre parámetros— en el triángulo inferior.)

    Parameters
    ----------
    samples     : np.ndarray, shape (n_samples, n_params)
        Posterior ensemble, e.g. ``NAAppraiser.samples``.
    param_names : list of str — axis labels (one per parameter).
    bounds      : np.ndarray, shape (n_params, 2), optional
        [min, max] per parameter, used to fix the axis ranges.
    truths      : np.ndarray, shape (n_params,), optional
        Reference values marked in red (e.g. the best-fit model).
    mean        : np.ndarray, shape (n_params,), optional
        Posterior mean, marked with a green cross on the 2D panels.
    bins        : int — histogram bins for 1D and 2D panels.

    Returns
    -------
    (fig, axes)
    """
    import matplotlib.pyplot as plt

    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 2:
        raise ValueError("samples must be 2D (n_samples, n_params).")
    nd = samples.shape[1]
    q = np.percentile(samples, [16, 50, 84], axis=0)  # (3, nd)

    fig, axes = plt.subplots(nd, nd, figsize=(2.3 * nd, 2.3 * nd), squeeze=False)
    for i in range(nd):
        for j in range(nd):
            ax = axes[i, j]
            if j > i:                      # upper triangle: hide
                ax.axis("off")
                continue
            if i == j:                     # diagonal: 1D marginal
                ax.hist(samples[:, i], bins=bins, density=True,
                        color="steelblue", histtype="stepfilled", alpha=0.75)
                lo, med, hi = q[:, i]
                for val, style in ((lo, ":"), (med, "--"), (hi, ":")):
                    ax.axvline(val, color="k", lw=0.8, ls=style)
                if truths is not None:
                    ax.axvline(truths[i], color="red", lw=1.6)
                ax.set_yticks([])
                ax.set_title(f"{param_names[i]}\n"
                             f"{med:.2f}  (+{hi-med:.2f} / -{med-lo:.2f})",
                             fontsize=8)
                if bounds is not None:
                    ax.set_xlim(bounds[i, 0], bounds[i, 1])
            else:                          # lower triangle: 2D density
                rng = None
                if bounds is not None:
                    rng = [bounds[j], bounds[i]]
                ax.hist2d(samples[:, j], samples[:, i], bins=bins,
                          range=rng, cmap="Blues")
                if truths is not None:
                    ax.plot(truths[j], truths[i], marker="*", color="red",
                            ms=11, mec="k", mew=0.5)
                if mean is not None:
                    ax.plot(mean[j], mean[i], marker="x", color="green", ms=7)
                if bounds is not None:
                    ax.set_xlim(bounds[j, 0], bounds[j, 1])
                    ax.set_ylim(bounds[i, 0], bounds[i, 1])

            # outer labels only
            if i == nd - 1:
                ax.set_xlabel(param_names[j], fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0 and i != 0:
                ax.set_ylabel(param_names[i], fontsize=8)
            elif i != j:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=6)

    if title:
        fig.suptitle(title, fontsize=11, y=1.0)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes
