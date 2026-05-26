"""
SAC/MiniSEED waveform processor for seismic inversion.

Implements robust signal processing following Boore (2001) and Boore & Bommer (2005)
standards for seismic data preprocessing:
- Zero-phase bandpass filtering with padding
- Polynomial baseline correction
- Multi-step integration with proper detrending
- Station and component extraction from metadata

This module bridges ObsPy with inversion workflows in kdellipspy.
"""

from __future__ import annotations
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
from obspy import Stream, Trace, UTCDateTime, read, read_inventory


# ==============================================================================
# CONFIGURATION & CONSTANTS
# ==============================================================================

class ProcessingConfig:
    """Configuration for waveform processing pipeline."""

    def __init__(
        self,
        freq_min: float = 0.04,
        freq_max: float = 0.2,
        filter_corners: int = 4,
        zerophase: bool = True,
        zeropad_duration: float = 100.0,
        taper_percentage: float = 0.05,
        taper_type: str = "cosine",
        polynomial_degree_vel: int = 1,
        polynomial_degree_disp: int = 2,
    ):
        """
        Args:
            freq_min: Minimum frequency (Hz) for bandpass filter.
            freq_max: Maximum frequency (Hz) for bandpass filter.
            filter_corners: Number of filter corners (order).
            zerophase: Apply forward-backward filtering for zero phase distortion.
            zeropad_duration: Padding duration (seconds) to avoid edge effects.
            taper_percentage: Taper percentage (0-1) applied before integration.
            taper_type: Taper type ('cosine', 'hamming', 'blackman', etc.).
            polynomial_degree_vel: Polynomial degree for baseline correction after ACC→VEL.
            polynomial_degree_disp: Polynomial degree for baseline correction after VEL→DISP.
        """
        self.freq_min = float(freq_min)
        self.freq_max = float(freq_max)
        self.filter_corners = int(filter_corners)
        self.zerophase = bool(zerophase)
        self.zeropad_duration = float(zeropad_duration)
        self.taper_percentage = float(taper_percentage)
        self.taper_type = str(taper_type)
        self.polynomial_degree_vel = int(polynomial_degree_vel)
        self.polynomial_degree_disp = int(polynomial_degree_disp)

        # Validate
        if not (0 < self.freq_min < self.freq_max):
            raise ValueError(f"Invalid frequency band: {self.freq_min} < {self.freq_max}")
        if not (0 <= self.taper_percentage <= 1):
            raise ValueError(f"Taper percentage must be in [0, 1], got {self.taper_percentage}")


# ==============================================================================
# SIGNAL PROCESSING FUNCTIONS
# ==============================================================================


def polynomial_baseline_correction(tr: Trace, degree: int = 1) -> Trace:
    """
    Apply polynomial baseline correction to waveform.

    References:
    - Boore (2001): Effect of baseline corrections on displacements and response spectra.
    - Graizer (2010): Baseline correction of strong-motion records.

    Typical usage:
    - degree=1 (linear): After ACC→VEL integration.
    - degree=2 (quadratic): After VEL→DISP integration (preserves static offset).

    Args:
        tr: ObsPy Trace.
        degree: Polynomial degree for fit.

    Returns:
        Modified Trace (operates in-place, but returns for chaining).
    """
    times = np.arange(tr.stats.npts) * tr.stats.delta
    coeffs = np.polyfit(times, tr.data, deg=degree)
    polynomial = np.polyval(coeffs, times)
    tr.data = tr.data - polynomial
    return tr


def acausal_bandpass(
    tr: Trace,
    freq_min: float,
    freq_max: float,
    corners: int = 4,
    zerophase: bool = True,
    zeropad_duration: float = 100.0,
) -> Trace:
    """
    Apply zero-phase bandpass filter with manual padding (Boore 2001, §3.2).

    Manual padding avoids ObsPy's automatic tapering, which can degrade
    low-frequency content. Critical for periods > 5 s.

    Args:
        tr: ObsPy Trace.
        freq_min: Low-frequency cutoff (Hz).
        freq_max: High-frequency cutoff (Hz).
        corners: Filter order.
        zerophase: Apply forward-backward filtering.
        zeropad_duration: Padding duration (seconds).

    Returns:
        Filtered Trace (modifies in-place).
    """
    dt = tr.stats.delta
    npad = int(zeropad_duration / dt)

    # Detrending and tapering BEFORE padding
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(max_percentage=0.05, type="cosine")

    # Manual padding with zeros
    padded_data = np.concatenate(
        [np.zeros(npad), tr.data, np.zeros(npad)]
    )
    tr_pad = tr.copy()
    tr_pad.data = padded_data
    tr_pad.stats.starttime = tr.stats.starttime - npad * dt

    # Apply filter
    tr_pad.filter(
        "bandpass",
        freqmin=freq_min,
        freqmax=freq_max,
        corners=corners,
        zerophase=zerophase,
    )

    # Extract original window
    tr_pad.trim(starttime=tr.stats.starttime, endtime=tr.stats.endtime)
    tr.data = tr_pad.data

    return tr


def process_trace(
    tr: Trace,
    event_time: UTCDateTime,
    target_mode: str = "VEL",
    config: Optional[ProcessingConfig] = None,
) -> Trace:
    """
    Complete preprocessing pipeline for seismic waveforms.

    Follows Boore (2001) and Boore & Bommer (2005) standards:

    1. Pre-event baseline correction (mean of window before earthquake).
    2. Linear detrending.
    3. Cosine taper BEFORE integration (critical for long-period signals).
    4. Multi-step numerical integration (cumtrapz) with polynomial baseline correction:
       - ACC→VEL: Linear polynomial correction (degree=1).
       - VEL→DISP: Quadratic polynomial correction (degree=2) to preserve static offset.
    5. Bandpass filtering with zero-phase and padding.
    6. Trim to desired time window.

    Args:
        tr: ObsPy Trace.
        event_time: Event origin time (UTCDateTime).
        target_mode: Target output ('VEL' or 'DISP').
        config: ProcessingConfig instance. If None, uses defaults.

    Returns:
        Processed Trace (modifies in-place).

    Raises:
        ValueError: If target_mode not in ('VEL', 'DISP').
    """
    if config is None:
        config = ProcessingConfig()

    if target_mode not in ("VEL", "DISP"):
        raise ValueError(f"target_mode must be 'VEL' or 'DISP', got {target_mode}")

    # ========================================================================
    # 1. Pre-event baseline correction
    # ========================================================================
    pre_event_window = tr.slice(endtime=event_time - 0.1)
    if pre_event_window is not None and pre_event_window.stats.npts > 5:
        tr.data -= np.mean(pre_event_window.data)
    else:
        tr.detrend("demean")

    tr.detrend("linear")

    # ========================================================================
    # 2. Taper BEFORE integration (Boore 2001, §3)
    # ========================================================================
    # Critical for long-period content; prevents edge discontinuities from
    # being amplified through each integration step.
    tr.taper(max_percentage=config.taper_percentage, type=config.taper_type)

    # ========================================================================
    # 3. Determine integration steps
    # ========================================================================
    # Sensor type inferred from channel code:
    # - HN/HL channels: Accelerometer
    # - Others: Velocimeter
    sensor_type = "ACC" if tr.stats.channel[0:2] in ("HN", "HL") else "VEL"

    if sensor_type == "ACC":
        n_integration_steps = 1 if target_mode == "VEL" else 2
    else:
        n_integration_steps = 0 if target_mode == "VEL" else 1

    # ========================================================================
    # 4. Integration with polynomial baseline correction
    # ========================================================================
    for step_idx in range(n_integration_steps):
        tr.integrate(method="cumtrapz")

        if step_idx == 0:
            # First integration: ACC→VEL (or VEL→DISP if only 1 step)
            # Linear baseline correction sufficient here.
            tr.detrend("demean")
            polynomial_baseline_correction(tr, degree=config.polynomial_degree_vel)

        elif step_idx == 1:
            # Second integration: VEL→DISP
            # Quadratic baseline correction (Boore 2001) to preserve permanent offset.
            # DO NOT apply linear detrend; it could remove real static displacement.
            tr.detrend("demean")
            polynomial_baseline_correction(tr, degree=config.polynomial_degree_disp)

    return tr


def _infer_component(channel: str) -> str:
    """
    Extract component code (E/N/Z) from channel string.

    Maps channel conventions:
    - HN/HL/BH ending in E/N/Z → E/N/Z
    - HN/HL/BH ending in 1 → N
    - HN/HL/BH ending in 2 → E
    """
    ch = (channel or "").upper()
    if not ch:
        return ""
    last = ch[-1]
    if last in {"E", "N", "Z"}:
        return last
    if last == "1":
        return "N"
    if last == "2":
        return "E"
    return ""


def _find_inventory_files(raw_dir: Path) -> Optional:
    """
    Scan raw_dir for inventory files (XML, RESP, StationXML).

    Returns combined inventory or None if no files found.
    """
    if read_inventory is None:
        return None

    patterns = ("*.xml", "*.XML", "*.stationxml", "*.STATIONXML", "RESP*", "*.resp", "*.RESP")
    inv_files = []
    for pat in patterns:
        inv_files.extend(raw_dir.glob(pat))

    if not inv_files:
        return None

    inv = None
    for fpath in inv_files:
        try:
            cur = read_inventory(str(fpath))
            inv = cur if inv is None else (inv + cur)
        except Exception as e:
            warnings.warn(f"Could not read inventory file {fpath}: {e}")

    return inv


# ==============================================================================
# WAVEFORM LOADING & EXTRACTION
# ==============================================================================


class WaveformExtractor:
    """
    Extract and process waveforms from SAC/MiniSEED files in DATA/RAW/.

    Handles:
    - Reading multiple files (*.sac, *.mseed)
    - Grouping by station and component
    - Preprocessing with full pipeline
    - Validation and error handling
    """

    def __init__(
        self,
        raw_dir: Path,
        config: Optional[ProcessingConfig] = None,
    ):
        """
        Args:
            raw_dir: Path to DATA/RAW directory.
            config: ProcessingConfig for signal processing.
        """
        self.raw_dir = Path(raw_dir)
        self.config = config or ProcessingConfig()

        if not self.raw_dir.is_dir():
            raise FileNotFoundError(f"RAW directory not found: {self.raw_dir}")

        # Load all waveforms
        self.stream = self._load_waveforms()
        self.inventory = _find_inventory_files(self.raw_dir)

    def _load_waveforms(self) -> Stream:
        """Load all SAC/MiniSEED files from raw_dir."""
        patterns = ("*.sac", "*.SAC", "*.mseed", "*.MSEED")
        stream = Stream()

        for pattern in patterns:
            for fpath in self.raw_dir.glob(pattern):
                try:
                    stream += read(str(fpath))
                except Exception as e:
                    warnings.warn(f"Could not read {fpath}: {e}")

        if not stream:
            raise FileNotFoundError(f"No SAC/MiniSEED files found in {self.raw_dir}")

        return stream

    def extract_stations(self) -> List[str]:
        """Get sorted list of unique station names."""
        return sorted(list(set(tr.stats.station for tr in self.stream)))

    def get_station_waveforms(self, station: str) -> Dict[str, Trace]:
        """
        Get preprocessed waveforms for a station, keyed by component.

        Args:
            station: Station name (e.g., 'AC01').

        Returns:
            Dict mapping component ('E', 'N', 'Z') to preprocessed Trace.

        Raises:
            ValueError: If station not found or components incomplete.
        """
        st_select = self.stream.select(station=station)
        if not st_select:
            raise ValueError(f"Station {station} not found in stream")

        result = {}
        for comp_code in ("N", "E", "Z"):
            candidates = [
                tr for tr in st_select
                if _infer_component(tr.stats.channel) == comp_code
            ]

            if not candidates:
                raise ValueError(f"Station {station}: component {comp_code} not found")

            # Prefer highest sample rate
            tr = sorted(
                candidates,
                key=lambda x: float(getattr(x.stats, "sampling_rate", 0.0)),
                reverse=True,
            )[0]

            result[comp_code] = tr

        return result

    def preprocess_station(
        self,
        station: str,
        event_time: UTCDateTime,
        target_mode: str = "VEL",
        time_window: Optional[Tuple[float, float]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Preprocess all components for a station.

        Args:
            station: Station name.
            event_time: Event origin time.
            target_mode: 'VEL' or 'DISP'.
            time_window: (t_start, t_end) relative to event_time. If None, uses full trace.

        Returns:
            Dict mapping component to data array.
        """
        waveforms = self.get_station_waveforms(station)
        result = {}

        for comp, tr in waveforms.items():
            # Copy for processing
            tr_proc = tr.copy()

            # Apply full preprocessing pipeline
            process_trace(tr_proc, event_time, target_mode, self.config)

            # Apply bandpass filtering
            acausal_bandpass(
                tr_proc,
                freq_min=self.config.freq_min,
                freq_max=self.config.freq_max,
                corners=self.config.filter_corners,
                zerophase=self.config.zerophase,
                zeropad_duration=self.config.zeropad_duration,
            )

            # Remove instrument response if inventory available
            if self.inventory is not None:
                try:
                    out_type = "DISP" if target_mode == "DISP" else "VEL"
                    pre_filt = (
                        max(0.001, self.config.freq_min * 0.5),
                        max(0.002, self.config.freq_min * 0.8),
                        self.config.freq_max * 1.2,
                        self.config.freq_max * 1.8,
                    )
                    tr_proc.remove_response(
                        inventory=self.inventory,
                        output=out_type,
                        pre_filt=pre_filt,
                        water_level=60,
                    )
                except Exception as e:
                    warnings.warn(f"Could not remove response for {station}.{comp}: {e}")

            # Trim to time window if specified
            if time_window is not None:
                t_start, t_end = time_window
                tr_proc.trim(
                    starttime=event_time + t_start,
                    endtime=event_time + t_end,
                    pad=True,
                    fill_value=0.0,
                )

            result[comp] = np.asarray(tr_proc.data, dtype=float)

        return result
