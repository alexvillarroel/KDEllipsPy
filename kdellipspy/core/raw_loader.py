"""
Integration layer between SAC processor and kdellipspy inversion workflows.

Provides convenient functions to load and preprocess SAC/MiniSEED data
from DATA/RAW/ for use in kinematic and moment tensor inversion.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional, Tuple
import numpy as np
from .config_parser import ConfigParser
from .sac_processor import ProcessingConfig, WaveformExtractor
from obspy import UTCDateTime


def load_observed_from_raw(
    input_ctl_path: str | Path = "input.ctl",
    data_dir: str | Path = "DATA",
    config: Optional[ProcessingConfig] = None,
    time_window_offset: Tuple[float, float] = (0.0, 128.0),
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load observed seismic waveforms from DATA/RAW/ SAC/MiniSEED files.

    Uses configuration from input.ctl and applies robust preprocessing:
    - Polynomial baseline correction
    - Zero-phase bandpass filtering with padding
    - Multi-step integration for ACC→VEL or ACC→DISP conversion
    - Instrument response removal (if inventory available)
    - Uniform resampling and time windowing

    Args:
        input_ctl_path: Path to input.ctl configuration file.
        data_dir: Base data directory (will look for DATA/RAW/). Defaults to "DATA".
        config: ProcessingConfig for signal processing. If None, uses defaults.
        time_window_offset: (t_start, t_end) relative to event origin (seconds).

    Returns:
        observed: (n_stations, 3, npts) array of [N, E, Z] components.
        time: (npts,) time axis starting from 0.0.

    Raises:
        FileNotFoundError: If DATA/RAW/ or input.ctl not found.
        ValueError: If stations or waveform components are missing.
    """
    input_ctl_path = Path(input_ctl_path).resolve()
    cfg = ConfigParser(str(input_ctl_path))

    # Resolve data directory
    data_dir = Path(data_dir)
    if not data_dir.is_absolute():
        data_dir = (input_ctl_path.parent / data_dir).resolve()

    raw_dir = data_dir / "RAW"

    # Initialize extractor
    if config is None:
        config = ProcessingConfig(
            freq_min=cfg.ellipse.freq1,
            freq_max=cfg.ellipse.freq2,
        )
    else:
        # Override frequency band from config if explicitly provided
        if config.freq_min != 0.04 or config.freq_max != 0.2:  # defaults
            pass
        else:
            config.freq_min = cfg.ellipse.freq1
            config.freq_max = cfg.ellipse.freq2

    extractor = WaveformExtractor(raw_dir, config=config)

    # Determine target output mode
    target_mode = "DISP" if cfg.observed_data.units == 1 else "VEL"

    # Parse event origin time
    event_origin_date = getattr(cfg.source_position, "event_origin_date", "") or ""
    event_origin_time = getattr(cfg.source_position, "event_origin_time", "") or ""

    # Validate that all required stations are available in waveform files
    available_stations = set(extractor.extract_stations())
    required_stations = set(s.name.strip().upper() for s in cfg.stations.stations)
    missing_stations = required_stations - available_stations
    if missing_stations:
        raise ValueError(
            f"Missing SAC/MiniSEED files for stations: {sorted(missing_stations)}. "
            f"Found: {sorted(available_stations)}"
        )

    if event_origin_date and event_origin_time:
        event_time = UTCDateTime(f"{event_origin_date}T{event_origin_time}")
    else:
        # Fallback: use earliest trace start time
        earliest = min(tr.stats.starttime for tr in extractor.stream)
        event_time = earliest + cfg.observed_data.t1

    # Extract station list from configuration
    station_names = [s.name.strip().upper() for s in cfg.stations.stations]
    npts = int(cfg.observed_data.npts)
    delta = float(cfg.observed_data.delta)

    # Preprocess each station and build output array
    time = cfg.observed_data.t1 + np.arange(npts, dtype=float) * delta
    observed = np.zeros((len(station_names), 3, npts), dtype=float)

    comp_idx = {"N": 0, "E": 1, "Z": 2}

    for ista, sta_name in enumerate(station_names):
        try:
            waveforms = extractor.preprocess_station(
                sta_name,
                event_time,
                target_mode=target_mode,
                time_window=time_window_offset,
            )

            for comp in ("N", "E", "Z"):
                data = waveforms[comp]

                # Resample and pad/trim to exact npts
                if len(data) < npts:
                    padded = np.zeros(npts, dtype=float)
                    padded[: len(data)] = data
                    data = padded
                elif len(data) > npts:
                    data = data[:npts]

                observed[ista, comp_idx[comp], :] = data

        except Exception as e:
            raise ValueError(f"Failed to preprocess station {sta_name}: {e}") from e

    # Final validation: check shape consistency
    if observed.shape[0] != len(station_names):
        raise RuntimeError(
            f"Internal consistency error: observed shape[0]={observed.shape[0]} "
            f"does not match station_names length={len(station_names)}"
        )
    if observed.shape[1] != 3:
        raise RuntimeError(
            f"Expected 3 components (N, E, Z), got {observed.shape[1]}"
        )
    if observed.shape[2] != npts:
        raise RuntimeError(
            f"Expected {npts} time points, got {observed.shape[2]}"
        )

    return observed, time
