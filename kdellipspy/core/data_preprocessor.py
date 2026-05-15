import numpy as np
from pathlib import Path
import logging
from typing import Any, List, Optional, Tuple, Dict
from .config_parser import ConfigParser
from .signal_utils import bandpass_filter_waveforms, integrate_waveforms, _load_from_raw
from pathlib import Path

logger = logging.getLogger(__name__)

class DataPreprocessor:
    """
    Utility to prepare raw seismic data for AXITRA inversion.
    (Utilidad para preparar datos sísmicos crudos para la inversión con AXITRA.)
    """

    def __init__(self, cfg: ConfigParser):
        self.cfg = cfg
        self.nsta_total = len(cfg.stations.stations)
        self.npts_cfg = cfg.observed_data.npts
        self.delta_cfg = cfg.observed_data.delta

    def plot_record_section(
        self,
        raw_dir: Path,
        t_start: Optional[Any] = None,
        t_end: Optional[Any] = None,
        freqmin: float = 0.05,
        freqmax: float = 0.5,
        scale: float = 2.0,
        components: List[str] = ['Z'],
        station_names: Optional[List[str]] = None,
    ):
        """
        Visualizes waveforms from SAC files in a record section, similar to the provided notebook.
        (Visualiza formas de onda de archivos SAC en una sección de registro, similar al notebook proporcionado.)
        """
        try:
            from obspy import read, UTCDateTime
            import matplotlib.pyplot as plt
            from obspy.geodetics import gps2dist_azimuth
        except ImportError:
            logger.error("ObsPy and Matplotlib are required for plot_record_section.")
            return

        raw_dir = Path(raw_dir)
        pattern = "*.SAC"
        st = read(str(raw_dir / pattern))
        
        if station_names:
            st = st.select(station=",".join(station_names))

        evlat = self.cfg.source_position.latitude
        evlon = self.cfg.source_position.longitude

        for comp in components:
            st_comp = st.select(component=comp).copy()
            if len(st_comp) == 0:
                logger.warning(f"No traces found for component {comp}")
                continue

            # Filter and trim
            st_comp.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
            
            if t_start is not None and t_end is not None:
                st_comp.trim(starttime=t_start, endtime=t_end, fill_value=0, pad=True)

            for tr in st_comp:
                try:
                    stlat, stlon = tr.stats.sac.stla, tr.stats.sac.stlo
                    tr.stats.distance = gps2dist_azimuth(evlat, evlon, stlat, stlon)[0]
                except AttributeError:
                    logger.warning(f"Metadata missing for {tr.id}, skipping distance calculation.")
                    tr.stats.distance = 0

            fig = plt.figure(figsize=(10, 12))
            st_comp.plot(
                type='section',
                orientation='horizontal',
                scale=scale,
                fig=fig,
                ev_coord=(evlat, evlon)
            )
            ax = fig.axes[0]
            ax.set_title(f'Component {comp} - Filtered {freqmin}-{freqmax} Hz')
            
            for tr in st_comp:
                ax.text(
                    0, tr.stats.distance/1000.0, tr.stats.station,
                    transform=ax.get_yaxis_transform(),
                    ha='right', va='center', fontsize=8
                )
            
            plt.show()

    def process_raw_files(
        self,
        raw_dir: Path,
        output_dir: Path,
        freq1: float,
        freq2: float,
        t_start: Any = 0.0,
        t_end: Optional[Any] = None,
        data_start_time: Optional[Any] = None,
        units: int = 2,  # 1: disp, 2: vel
        station_indices: Optional[List[int]] = None,
        station_names: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Loads concatenated raw velocity files, filters, trims with zero-padding,
        and saves Axitra-ready files. Supports UTCDateTime or float for timing.
        
        If station_names is provided, it filters the stations by name.
        (Si se proporciona station_names, filtra las estaciones por nombre.)
        """
        try:
            from obspy import UTCDateTime
        except ImportError:
            UTCDateTime = None

        raw_dir = Path(raw_dir)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # If RAW waveform files exist (SAC/MSEED), prefer RAW mode and load them.
        # Check both the provided raw_dir and a possible raw_dir/'RAW' subfolder.
        raw_candidates = [raw_dir, raw_dir / "RAW"]
        has_raw = False
        selected_raw_dir = None
        for cand in raw_candidates:
            if any(cand.glob("*.SAC")) or any(cand.glob("*.sac")) or any(cand.glob("*.mseed")) or any(cand.glob("*.MSEED")):
                has_raw = True
                selected_raw_dir = cand
                break

        if has_raw and selected_raw_dir is not None:
            # _load_from_raw returns (observed, time) with shape (n_stations, 3, npts)
            observed, time = _load_from_raw(selected_raw_dir, self.cfg, freq1, freq2)

            # Map components to axis 1: we follow the convention returned by _load_from_raw
            # which uses component order [E, N, Z] by internal mapping. We'll export files
            # as real_vel_x -> E, real_vel_y -> N, real_vel_z -> Z for compatibility.
            comp_map = { 'x': 0, 'y': 1, 'z': 2 }
            processed_data = {}
            for comp in ['x', 'y', 'z']:
                arr = observed[:, comp_map[comp], :]
                if units == 1:
                    arr = integrate_waveforms(arr[:, np.newaxis, :], self.delta_cfg).squeeze(1)
                    prefix = 'real_disp'
                else:
                    prefix = 'real_vel'

                if station_indices is not None:
                    arr = arr[station_indices, :]

                out_name = f"{prefix}_{comp}"
                np.savetxt(output_dir / out_name, arr.flatten(), fmt='%.8e')
                processed_data[comp] = arr
                print(f"✓ Saved: {output_dir / out_name} (npts={arr.shape[1]})")

            return processed_data

        # Convert UTCDateTime to relative seconds if needed
        def to_rel(t, ref):
            if UTCDateTime and isinstance(t, UTCDateTime):
                if ref is None:
                    # If no reference provided, try to use origin time from config
                    if self.cfg.source_position.origin_time:
                        ref = UTCDateTime(self.cfg.source_position.origin_time)
                    else:
                        raise ValueError("data_start_time (UTCDateTime) is required if t_start/t_end are UTCDateTime and no Origin Time is in config.")
                return float(t - ref)
            return float(t)

        t_start_rel = to_rel(t_start, data_start_time)
        
        if t_end is None:
            t_end_rel = t_start_rel + (self.npts_cfg * self.delta_cfg)
        else:
            t_end_rel = to_rel(t_end, data_start_time)

        # Handle station filtering
        if station_names is not None:
            all_sta_names = [s.name for s in self.cfg.stations.stations]
            station_indices = [all_sta_names.index(name) for name in station_names if name in all_sta_names]

        components = ['x', 'y', 'z']
        processed_data = {}

        print(f"\n--- Data Preprocessor ---")
        print(f"(!) Warning: Input files must be in Velocity (m/s).")
        print(f"Target: npts={self.npts_cfg}, delta={self.delta_cfg}s")
        print(f"Window: rel [{t_start_rel:.2f}, {t_end_rel:.2f}]s")
        
        for comp in components:
            file_name = f"real_vel_{comp}"
            file_path = raw_dir / file_name
            if not file_path.exists():
                logger.warning(f"File {file_path} not found. Skipping.")
                continue

            raw_array = np.loadtxt(file_path)
            npts_raw_total = len(raw_array) // self.nsta_total
            data_raw = raw_array.reshape(self.nsta_total, npts_raw_total)

            # Extract and pad
            extracted = np.zeros((self.nsta_total, self.npts_cfg))
            
            for i in range(self.nsta_total):
                time_raw = np.arange(npts_raw_total) * self.delta_cfg
                
                # Selection logic with zero padding
                # Find indices in raw that fall within [t_start_rel, t_end_rel]
                mask = (time_raw >= t_start_rel) & (time_raw < t_end_rel)
                idx_in_raw = np.where(mask)[0]
                
                if len(idx_in_raw) > 0:
                    # How many points to copy?
                    pts_to_copy = min(len(idx_in_raw), self.npts_cfg)
                    extracted[i, :pts_to_copy] = data_raw[i, idx_in_raw[:pts_to_copy]]

            # Filter
            time_target = np.arange(self.npts_cfg) * self.delta_cfg
            data_filtered = bandpass_filter_waveforms(
                extracted[:, np.newaxis, :],
                time_target,
                freq1=freq1,
                freq2=freq2
            ).squeeze(1)

            if units == 1:
                data_filtered = integrate_waveforms(
                    data_filtered[:, np.newaxis, :], 
                    self.delta_cfg
                ).squeeze(1)
                prefix = "real_disp"
            else:
                prefix = "real_vel"

            if station_indices is not None:
                data_filtered = data_filtered[station_indices, :]

            out_name = f"{prefix}_{comp}"
            np.savetxt(output_dir / out_name, data_filtered.flatten(), fmt='%.8e')
            processed_data[comp] = data_filtered
            print(f"✓ Saved: {output_dir / out_name} (npts={data_filtered.shape[1]})")

        return processed_data
