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
    dpi: int = 300
) -> Tuple[Any, Any]:
    """
    Plot observed vs synthetic waveforms.
    """
    import matplotlib.pyplot as plt
    
    nsta = observed.shape[0]
    fig, axs = plt.subplots(nsta, 3, figsize=(10, 1.5 * nsta), squeeze=False, sharex=True)
    
    comps = ['North', 'East', 'Z']
    colors_obs = ['black', 'black', 'black']
    colors_syn = ['red', 'red', 'red']

    title = "Waveform Fit"
    if misfit is not None:
        title += f" (Misfit: {misfit:.4f})"
    fig.suptitle(title, fontsize=14, y=1.02)

    for i in range(nsta):
        for j in range(3):
            ax = axs[i, j]
            ax.plot(time, observed[i, j, :], color=colors_obs[j], label='Observed' if i==0 else None, linewidth=1.5)
            ax.plot(time, synthetic[i, j, :], color=colors_syn[j], label='Synthetic' if i==0 else None, linewidth=1.2, linestyle='--')
            
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
    Plot 2D slip distribution interpolated on fault plane.
    """
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    
    slip_map = {}
    for sp in geometry.source_points:
        sf_idx = sp.subfault_index
        slip_map[sf_idx] = max(slip_map.get(sf_idx, 0.0), sp.displacement)
        
    x_pts = np.array([sf.x_m for sf in geometry.subfaults]) / 1000.0
    y_pts = np.array([sf.y_m for sf in geometry.subfaults]) / 1000.0
    z_pts = np.array([slip_map.get(sf.index, 0.0) for sf in geometry.subfaults])

    resolucion = 100
    xi = np.linspace(x_pts.min(), x_pts.max(), resolucion)
    yi = np.linspace(y_pts.min(), y_pts.max(), resolucion)
    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x_pts, y_pts), z_pts, (X, Y), method='cubic', fill_value=0.0)
    Z = np.clip(Z, 0, None)
    max_slip = np.max(Z)

    fig, ax = plt.subplots(figsize=(10, 6))
    mapa_calor = ax.contourf(X, Y, Z, levels=50, cmap='gnuplot')

    umbral_verde = max_slip * 0.15
    if max_slip > 0:
        ax.contour(X, Y, Z, levels=[umbral_verde], colors='lime', linestyles='dashed', linewidths=1.5)

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
        axes[1, 1].scatter(
            [row.get(p0, np.nan) for row in source_data],
            [row.get(p1, np.nan) for row in source_data],
            c=misfits,
            cmap="viridis",
            s=20,
        )
        axes[1, 1].set_xlabel(p0)
        axes[1, 1].set_ylabel(p1)
        axes[1, 1].set_title("Parameter cloud colored by misfit")
        axes[1, 1].grid(True, alpha=0.3)
        cbar = fig.colorbar(axes[1, 1].collections[0], ax=axes[1, 1])
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

    vmin, vmax = float(np.nanmin(misfits)), float(np.nanmax(misfits))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
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
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Misfit", fontsize=9)

    fig.suptitle("Parameter convergence — NA search", fontsize=11, y=1.01)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes
