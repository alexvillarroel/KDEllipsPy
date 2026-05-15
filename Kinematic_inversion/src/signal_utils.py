from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

try:
	from .config_parser import ConfigParser
except ImportError:
	from config_parser import ConfigParser

try:
	from obspy import Stream, Trace, read, read_inventory
except Exception:  # pragma: no cover - handled at runtime when RAW mode is requested
	Stream = None
	Trace = None
	read = None
	read_inventory = None

try:
	from obspy.taup import TauPyModel
	from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
except Exception:  # pragma: no cover - handled at runtime when azi_times is requested
	TauPyModel = None
	gps2dist_azimuth = None
	kilometers2degrees = None

try:
	from obspy.clients.iris import Client as IrisClient
except Exception:  # pragma: no cover - optional fallback to legacy web service
	IrisClient = None


def _component_suffix_for_units(units: int) -> str:
	if int(units) == 1:
		return "disp"
	if int(units) == 2:
		return "vel"
	raise ValueError(f"Unsupported units={units}. Expected 1 (disp) or 2 (vel).")


def _read_observed_flat_file(
	path: Path,
	t1: float,
	delta: float,
	npts: int,
	n_stations: int,
) -> Tuple[np.ndarray, np.ndarray]:
	"""
	Read flat observed waveform files in either format:
	- Two columns: time value
	- One column: value only (time is inferred from input.ctl)
	"""
	arr = np.loadtxt(path, dtype=float)
	arr = np.asarray(arr, dtype=float)

	# Case A: value-only file (1D or Nx1)
	if arr.ndim == 1:
		data = arr
		expected_rows = int(n_stations) * int(npts)
		if data.size != expected_rows:
			raise ValueError(
				f"{path}: value-only file has {data.size} rows, expected {expected_rows} "
				f"(npts={npts} x nstations={n_stations})"
			)
		time_block = float(t1) + np.arange(int(npts), dtype=float) * float(delta)
		time = np.tile(time_block, int(n_stations))
	elif arr.ndim == 2 and arr.shape[1] == 1:
		data = arr[:, 0]
		expected_rows = int(n_stations) * int(npts)
		if data.size != expected_rows:
			raise ValueError(
				f"{path}: value-only file has {data.size} rows, expected {expected_rows} "
				f"(npts={npts} x nstations={n_stations})"
			)
		time_block = float(t1) + np.arange(int(npts), dtype=float) * float(delta)
		time = np.tile(time_block, int(n_stations))
	# Case B: classic time/value file
	elif arr.ndim == 2 and arr.shape[1] >= 2:
		time = arr[:, 0]
		data = arr[:, 1]
	else:
		raise ValueError(
			f"{path} must contain either one column (value) or at least two columns (time value)"
		)

	if not np.all(np.isfinite(time)) or not np.all(np.isfinite(data)):
		raise ValueError(f"{path} contains NaN/Inf values")

	return time, data


def _reshape_flat_component(
	time: np.ndarray,
	data: np.ndarray,
	n_stations: int,
	npts: int,
	label: str,
) -> Tuple[np.ndarray, np.ndarray]:
	expected_rows = n_stations * npts
	if data.size != expected_rows:
		raise ValueError(
			f"{label}: expected {expected_rows} rows (npts={npts} x nstations={n_stations}), got {data.size}"
		)

	time_2d = time.reshape(n_stations, npts)
	data_2d = data.reshape(n_stations, npts)

	ref_t = time_2d[0]
	for i in range(1, n_stations):
		if not np.allclose(time_2d[i], ref_t, atol=1e-6, rtol=0.0):
			raise ValueError(
				f"{label}: time column mismatch between station blocks 0 and {i}"
			)

	return ref_t, data_2d


def _validate_time_axis(time: np.ndarray, t1: float, delta: float, npts: int) -> None:
	if time.size != npts:
		raise ValueError(f"Time axis has {time.size} samples, expected npts={npts}")

	if time.size > 1:
		dt_measured = float(np.median(np.diff(time)))
		if not np.isclose(dt_measured, float(delta), atol=1e-6, rtol=1e-3):
			raise ValueError(
				f"Time step mismatch: measured dt={dt_measured:.6g}, input.ctl delta={float(delta):.6g}"
			)

	if not np.isclose(float(time[0]), float(t1), atol=max(1e-6, abs(float(delta)) * 0.51), rtol=0.0):
		raise ValueError(
			f"Time start mismatch: file starts at {float(time[0]):.6g}, expected t1={float(t1):.6g}"
		)


def _load_from_flat_files(data_dir: Path, cfg: ConfigParser) -> Tuple[np.ndarray, np.ndarray]:
	n_stations = len(cfg.stations.stations)
	npts = int(cfg.observed_data.npts)
	units_suffix = _component_suffix_for_units(cfg.observed_data.units)

	fx = data_dir / f"real_{units_suffix}_x"
	fy = data_dir / f"real_{units_suffix}_y"
	fz = data_dir / f"real_{units_suffix}_z"

	missing = [str(p) for p in (fx, fy, fz) if not p.exists()]
	if missing:
		raise FileNotFoundError(
			"Missing flat observed files for selected units in input.ctl: "
			+ ", ".join(missing)
		)

	tx, x = _read_observed_flat_file(
		fx,
		t1=float(cfg.observed_data.t1),
		delta=float(cfg.observed_data.delta),
		npts=npts,
		n_stations=n_stations,
	)
	ty, y = _read_observed_flat_file(
		fy,
		t1=float(cfg.observed_data.t1),
		delta=float(cfg.observed_data.delta),
		npts=npts,
		n_stations=n_stations,
	)
	tz, z = _read_observed_flat_file(
		fz,
		t1=float(cfg.observed_data.t1),
		delta=float(cfg.observed_data.delta),
		npts=npts,
		n_stations=n_stations,
	)

	time, x2 = _reshape_flat_component(tx, x, n_stations, npts, fx.name)
	time_y, y2 = _reshape_flat_component(ty, y, n_stations, npts, fy.name)
	time_z, z2 = _reshape_flat_component(tz, z, n_stations, npts, fz.name)

	if not np.allclose(time_y, time, atol=1e-6, rtol=0.0):
		raise ValueError(f"Time mismatch between {fx.name} and {fy.name}")
	if not np.allclose(time_z, time, atol=1e-6, rtol=0.0):
		raise ValueError(f"Time mismatch between {fx.name} and {fz.name}")

	_validate_time_axis(
		time,
		float(cfg.observed_data.t1),
		float(cfg.observed_data.delta),
		npts,
	)

	# (n_stations, 3, npts)
	observed = np.stack([x2, y2, z2], axis=1)
	return observed, time


def _iter_raw_waveform_files(raw_dir: Path) -> Iterable[Path]:
	for ext in ("*.sac", "*.SAC", "*.mseed", "*.MSEED"):
		yield from raw_dir.glob(ext)


def _component_from_channel(channel: str) -> str:
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


def _find_inventory(raw_dir: Path):
	if read_inventory is None:
		return None
	inv_files = []
	for pat in ("*.xml", "*.XML", "*.stationxml", "*.STATIONXML", "RESP*", "*.resp", "*.RESP"):
		inv_files.extend(raw_dir.glob(pat))
	if not inv_files:
		return None
	inv = None
	for f in inv_files:
		cur = read_inventory(str(f))
		inv = cur if inv is None else (inv + cur)
	return inv


def _preprocess_trace(
	tr,
	freq1: float,
	freq2: float,
	delta: float,
	npts: int,
	units: int,
	inventory,
	starttime,
):
	tr = tr.copy()
	tr.detrend("demean")
	tr.detrend("linear")
	tr.taper(max_percentage=0.05, type="cosine")

	if inventory is not None:
		pre_filt = (
			max(0.001, float(freq1) * 0.5),
			max(0.002, float(freq1) * 0.8),
			float(freq2) * 1.2,
			float(freq2) * 1.8,
		)
		out = "DISP" if int(units) == 1 else "VEL"
		tr.remove_response(
			inventory=inventory,
			output=out,
			pre_filt=pre_filt,
			water_level=60,
		)

	tr.filter("bandpass", freqmin=float(freq1), freqmax=float(freq2), corners=4, zerophase=True)
	tr.interpolate(sampling_rate=1.0 / float(delta), starttime=starttime, method="linear")

	endtime = starttime + (int(npts) - 1) * float(delta)
	tr.trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0.0)

	data = np.asarray(tr.data, dtype=float)
	if data.size < int(npts):
		pad = np.zeros(int(npts), dtype=float)
		pad[: data.size] = data
		data = pad
	elif data.size > int(npts):
		data = data[: int(npts)]
	return data


def _load_from_raw(raw_dir: Path, cfg: ConfigParser, freq1: float, freq2: float) -> Tuple[np.ndarray, np.ndarray]:
	if read is None:
		raise ImportError(
			"ObsPy is required for RAW mode (.sac/.mseed). Install with: pip install obspy"
		)

	paths = sorted(_iter_raw_waveform_files(raw_dir))
	if not paths:
		raise FileNotFoundError(f"No .sac/.mseed files found in {raw_dir}")

	stream = Stream()
	for p in paths:
		stream += read(str(p))

	station_names = [s.name.strip().upper() for s in cfg.stations.stations]
	inventory = _find_inventory(raw_dir)

	if inventory is None:
		print("[signal_utils] WARNING: no instrument inventory found in DATA/RAW; response removal skipped.")

	# Reference start from earliest available trace. Then apply t1 offset from input.ctl.
	global_start = min(tr.stats.starttime for tr in stream)
	starttime = global_start + float(cfg.observed_data.t1)

	npts = int(cfg.observed_data.npts)
	delta = float(cfg.observed_data.delta)
	time = float(cfg.observed_data.t1) + np.arange(npts, dtype=float) * delta

	comp_idx = {"E": 0, "N": 1, "Z": 2}
	observed = np.zeros((len(station_names), 3, npts), dtype=float)

	missing = []
	for ista, sta in enumerate(station_names):
		st_sta = stream.select(station=sta)
		if len(st_sta) == 0:
			missing.append(f"{sta}: all components")
			continue

		for comp_name, icomp in comp_idx.items():
			candidates = [tr for tr in st_sta if _component_from_channel(getattr(tr.stats, "channel", "")) == comp_name]
			if not candidates:
				missing.append(f"{sta}: {comp_name}")
				continue

			# Prefer highest sample rate if duplicates exist.
			tr = sorted(candidates, key=lambda x: float(getattr(x.stats, "sampling_rate", 0.0)), reverse=True)[0]
			observed[ista, icomp, :] = _preprocess_trace(
				tr,
				freq1=float(freq1),
				freq2=float(freq2),
				delta=delta,
				npts=npts,
				units=int(cfg.observed_data.units),
				inventory=inventory,
				starttime=starttime,
			)

	if missing:
		raise ValueError(
			"Missing required station/components in RAW data: " + ", ".join(missing)
		)

	return observed, time


def _first_arrival_time(model, evdep_km: float, distdeg: float, phases: Tuple[str, str]) -> float:
	for phase in phases:
		arrivals = model.get_travel_times(float(evdep_km), float(distdeg), phase_list=[phase])
		if arrivals:
			return float(arrivals[0].time)
	raise ValueError(f"No arrivals found for phases {phases} at distance {distdeg:.3f} deg")


def _distance_azimuth_deg(
	stlat: float,
	stlon: float,
	evlat: float,
	evlon: float,
	use_iris: bool,
	iris_client,
) -> Tuple[float, float]:
	if use_iris:
		if iris_client is None:
			raise RuntimeError("IRIS client is unavailable. Install/enable obspy.clients.iris or set use_iris=False.")
		result = iris_client.distaz(stlat, stlon, evlat, evlon)
		return float(result["distance"]), float(result["azimuth"])

	if gps2dist_azimuth is None or kilometers2degrees is None:
		raise RuntimeError("ObsPy geodetics is unavailable. Install obspy to compute azimuth and distance.")

	dist_m, azimuth_deg, _ = gps2dist_azimuth(stlat, stlon, evlat, evlon)
	distdeg = float(kilometers2degrees(float(dist_m) / 1000.0))
	return distdeg, float(azimuth_deg)


def build_azi_times_array(
	input_ctl_path: str | Path = "input.ctl",
	model_name: str = "iasp91",
	p_shift_s: float = -1.0,
	s_shift_s: float = -1.0,
	use_iris: bool = False,
) -> np.ndarray:
	"""
	Build legacy-compatible azi_times values as an array with columns:
	[azimuth_rad, tP, tS].

	This reproduces the old script behavior while avoiding IRIS network
	requirements by default (use_iris=False).
	"""
	if TauPyModel is None:
		raise RuntimeError("ObsPy TauPyModel is unavailable. Install obspy to compute travel times.")

	input_ctl_path = Path(input_ctl_path).resolve()
	cfg = ConfigParser(str(input_ctl_path))

	evlat = float(cfg.source_position.latitude)
	evlon = float(cfg.source_position.longitude)
	evdep = float(cfg.source_position.depth)
	stations = cfg.stations.stations

	model = TauPyModel(model=model_name)
	iris_client = IrisClient() if use_iris and IrisClient is not None else None

	rows = []
	for st in stations:
		stlat = float(st.latitude)
		stlon = float(st.longitude)

		distdeg, azi_deg = _distance_azimuth_deg(
			stlat=stlat,
			stlon=stlon,
			evlat=evlat,
			evlon=evlon,
			use_iris=bool(use_iris),
			iris_client=iris_client,
		)

		pt = _first_arrival_time(model, evdep, distdeg, ("P", "p"))
		stt = _first_arrival_time(model, evdep, distdeg, ("S", "s"))

		rows.append([
			math.radians(float(azi_deg)),
			float(pt) + float(p_shift_s),
			float(stt) + float(s_shift_s),
		])

	return np.asarray(rows, dtype=float)


def write_azi_times_file(
	input_ctl_path: str | Path = "input.ctl",
	output_path: str | Path | None = None,
	model_name: str = "iasp91",
	p_shift_s: float = -1.0,
	s_shift_s: float = -1.0,
	use_iris: bool = False,
) -> Path:
	"""Write Event/azi_times.txt with legacy-compatible 3-column format."""
	input_ctl_path = Path(input_ctl_path).resolve()
	if output_path is None:
		output = input_ctl_path.parent / "Event" / "azi_times.txt"
	else:
		output = Path(output_path)
		if not output.is_absolute():
			output = (input_ctl_path.parent / output).resolve()

	arr = build_azi_times_array(
		input_ctl_path=input_ctl_path,
		model_name=model_name,
		p_shift_s=p_shift_s,
		s_shift_s=s_shift_s,
		use_iris=use_iris,
	)

	output.parent.mkdir(parents=True, exist_ok=True)
	np.savetxt(str(output), arr, fmt="%.4f %.1f %.1f")
	return output


def load_and_filter_observed_data(
	freq1: Optional[float] = None,
	freq2: Optional[float] = None,
	input_ctl_path: str | Path = "input.ctl",
	data_dir: str | Path = "DATA",
	prefer_raw: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
	"""
	Load observed data with two supported user workflows:

	1) Flat files in DATA:
	   - real_disp_x/y/z  (if input.ctl units=1)
	   - real_vel_x/y/z   (if input.ctl units=2)
	   Each file must contain 2 columns: time and data,
	   and total rows must be npts * nstations.

	2) RAW files in DATA/RAW:
	   - .sac or .mseed files
	   - station name extracted from waveform metadata
	   - bandpass filtering using input.ctl frequency band
	   - response removal attempted when inventory files are present in RAW

	Returns:
		observed_waveforms: ndarray with shape (n_stations, 3, npts)
		time_array: ndarray with shape (npts,)
	"""
	input_ctl_path = Path(input_ctl_path).resolve()
	base_dir = input_ctl_path.parent
	data_dir = Path(data_dir)
	if not data_dir.is_absolute():
		data_dir = (base_dir / data_dir).resolve()

	cfg = ConfigParser(str(input_ctl_path))
	freq1_eff = float(cfg.ellipse.freq1 if freq1 is None else freq1)
	freq2_eff = float(cfg.ellipse.freq2 if freq2 is None else freq2)

	if freq1_eff <= 0 or freq2_eff <= 0 or freq2_eff <= freq1_eff:
		raise ValueError(
			f"Invalid frequency band: freq1={freq1_eff}, freq2={freq2_eff}. Require 0 < freq1 < freq2."
		)

	raw_dir = data_dir / "RAW"

	if prefer_raw and raw_dir.exists():
		return _load_from_raw(raw_dir, cfg, freq1_eff, freq2_eff)

	# Try flat-file mode first (fast and deterministic).
	try:
		return _load_from_flat_files(data_dir, cfg)
	except Exception as flat_exc:
		if raw_dir.exists():
			print(f"[signal_utils] INFO: flat DATA mode unavailable ({flat_exc}); trying RAW mode...")
			return _load_from_raw(raw_dir, cfg, freq1_eff, freq2_eff)
		raise


def bandpass_filter_waveforms(
	synthetic: np.ndarray,
	time: np.ndarray,
	freq1: float,
	freq2: float,
	corners: int = 4,
	zerophase: bool = True,
) -> np.ndarray:
	"""
	Apply Butterworth bandpass filter to synthetic waveforms.
	
	This ensures synthetic waveforms match the frequency band of observed data,
	maintaining consistency in the misfit calculation.

	Args:
		synthetic: waveform array with shape (n_stations, 3, npts)
		time: time array with shape (npts,)
		freq1: low frequency cutoff (Hz)
		freq2: high frequency cutoff (Hz)
		corners: filter order (default 4 = 4th order)
		zerophase: if True, apply filter forward-backward for zero phase distortion

	Returns:
		filtered_synthetic: filtered waveforms with same shape as input

	Raises:
		ValueError: if time array has < 2 points or freq2 >= Nyquist frequency
	"""
	try:
		from scipy.signal import butter, sosfilt
	except ImportError:
		raise ImportError("scipy is required for bandpass filtering. Install with: pip install scipy")

	if len(time) < 2:
		raise ValueError(f"Need at least 2 time points, got {len(time)}")

	# Infer sampling rate from time
	dt = float(time[1] - time[0])
	if dt <= 0:
		raise ValueError(f"Invalid time step: {dt}")

	sr = 1.0 / dt
	nyquist = sr / 2.0

	# Validate frequencies
	if freq1 <= 0 or freq2 <= 0:
		raise ValueError(f"Frequencies must be positive: freq1={freq1}, freq2={freq2}")
	if freq1 >= freq2:
		raise ValueError(f"freq1 ({freq1}) must be < freq2 ({freq2})")
	if freq2 >= nyquist:
		raise ValueError(
			f"freq2 ({freq2}) >= Nyquist frequency ({nyquist:.2f}). "
			f"Increase sampling rate or reduce freq2."
		)

	# Design Butterworth filter in SOS (second-order sections) form for better numerics
	sos = butter(int(corners), [float(freq1), float(freq2)], btype='band', fs=sr, output='sos')

	# Apply to each station and component
	nstations, ncomps, npts = synthetic.shape
	filtered = np.zeros_like(synthetic, dtype=float)

	for ista in range(nstations):
		for icomp in range(ncomps):
			data = synthetic[ista, icomp, :]

			# Forward pass
			filtered_fwd = sosfilt(sos, data)

			if zerophase:
				# Reverse, filter again, reverse back (for zero phase)
				filtered_rev = sosfilt(sos, filtered_fwd[::-1])
				filtered[ista, icomp, :] = filtered_rev[::-1]
			else:
				filtered[ista, icomp, :] = filtered_fwd

	return filtered
