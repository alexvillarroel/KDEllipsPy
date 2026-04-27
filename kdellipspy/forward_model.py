from pathlib import Path
from typing import Optional, Tuple
import importlib
import sys
import io
import os
from contextlib import redirect_stderr, redirect_stdout

import numpy as np

# Handle both relative imports (when used as package) and direct imports (from notebooks)
try:
	from .config_parser import ConfigParser
	from .geometry import EllipticalSlipMapper, FaultGeometry, GeometryBuilder, UTMProjection
except ImportError:
	from config_parser import ConfigParser
	from geometry import EllipticalSlipMapper, FaultGeometry, GeometryBuilder, UTMProjection


class AxitraForwardModel:
	"""
	Forward model using AXITRA2024 python wrapper.

	Main workflow:
	1) Build invariant geometry from input.ctl
	2) Apply ellipse model (a1,a2,theta,np,tp,dmax,vr)
	3) Build Axitra instance with model/stations/sources
	4) Compute Green functions: moment.green(ap)
	5) Convolve: moment.conv(ap, hist, ...)
	"""

	def __init__(self, input_ctl_path: str, axitra_dir: Optional[str] = None):
		self.input_ctl_path = Path(input_ctl_path).resolve()
		self.base_dir = self.input_ctl_path.parent
		self.cfg = ConfigParser(str(self.input_ctl_path))
		self.geometry_builder = GeometryBuilder(self.cfg)
		self.ellipse_mapper = EllipticalSlipMapper(self.cfg)
		self._projection = UTMProjection(
			float(self.cfg.source_position.latitude),
			float(self.cfg.source_position.longitude),
		)
		src_east, src_north = self._projection.latlon_to_xy(
			float(self.cfg.source_position.latitude),
			float(self.cfg.source_position.longitude),
		)
		self._src_east_m = float(src_east)
		self._src_north_m = float(src_north)
		self.axitra_dir = (
			Path(axitra_dir).resolve()
			if axitra_dir
			else (self.base_dir / "AXITRA2024").resolve()
		)

	@classmethod
	def from_params(cls, params: dict, axitra_dir: Optional[str] = None) -> 'AxitraForwardModel':
		"""
		Create AxitraForwardModel from a dictionary of parameters instead of an input.ctl file.
		
		This allows building the forward model programmatically without file I/O.
		
		Parameters:
		-----------
		params : dict
			Dictionary containing all configuration sections:
			- 'observed_data': Dictionary with time window and sampling parameters
			- 'source_position': Dictionary with event location and focal mechanism
			- 'fault_plane': Dictionary with fault geometry discretization
			- 'ellipse': Dictionary with ellipse model and frequency band
			- 'inversion_params': Dictionary with inversion parameter ranges
			- 'inversion_process': Dictionary with algorithm settings
			- 'moment_tensor': Dictionary with moment tensor components
			- 'stations': List of station dictionaries (lat, lon, height, name)
			- 'velocity_model': List of velocity layer dictionaries
		
		axitra_dir : str, optional
			Path to AXITRA2024 directory. If not provided, assumes it's a sibling folder.
		
		Returns:
		--------
		AxitraForwardModel
			Initialized forward model instance with all configuration loaded from params.
		
		Example:
		--------
		>>> params = {
		...     'observed_data': {'Time window start (t1)': 0.0, ...},
		...     'source_position': {'Latitude': -20.0, 'Longitude': -70.0, ...},
		...     'fault_plane': {'Length along strike (Lx)': 100000.0, ...},
		...     'ellipse': {'Number of ellipses': 1, ...},
		...     'stations': [
		...         {'latitude': -19.5, 'longitude': -69.5, 'height': 0.0, 'name': 'SL01'},
		...         ...
		...     ],
		...     'velocity_model': [
		...         {'thickness': 10000.0, 'vp': 5500.0, 'vs': 3200.0, 'rho': 2700.0, 'qp': 100.0, 'qs': 50.0},
		...         ...
		...     ],
		...     'inversion_params': {...},
		...     'inversion_process': {...},
		...     'moment_tensor': {...},
		... }
		>>> fm = AxitraForwardModel.from_params(params, axitra_dir='/path/to/AXITRA2024')
		"""
		# Create ConfigParser from dict
		cfg = ConfigParser.from_dict(params)
		
		# Create instance without calling __init__ (avoid file parsing)
		instance = cls.__new__(cls)
		instance.input_ctl_path = Path('<from_params>')
		instance.base_dir = Path.cwd()
		instance.cfg = cfg
		instance.geometry_builder = GeometryBuilder(cfg)
		instance.ellipse_mapper = EllipticalSlipMapper(cfg)
		instance._projection = UTMProjection(
			float(cfg.source_position.latitude),
			float(cfg.source_position.longitude),
		)
		src_east, src_north = instance._projection.latlon_to_xy(
			float(cfg.source_position.latitude),
			float(cfg.source_position.longitude),
		)
		instance._src_east_m = float(src_east)
		instance._src_north_m = float(src_north)
		
		# Set AXITRA directory
		if axitra_dir:
			instance.axitra_dir = Path(axitra_dir).resolve()
		else:
			# Try to find AXITRA2024 relative to cwd
			axitra_candidate = Path.cwd() / "AXITRA2024"
			if axitra_candidate.exists():
				instance.axitra_dir = axitra_candidate.resolve()
			else:
				raise FileNotFoundError(
					"AXITRA2024 directory not found. Please provide axitra_dir parameter."
				)
		
		return instance

	def _import_axitra(self):
		if str(self.axitra_dir) not in sys.path:
			sys.path.append(str(self.axitra_dir))
		mod = importlib.import_module("axitra")
		return mod.Axitra, mod.moment

	def _call_with_optional_silence(self, fn, *args, quiet: bool = True, **kwargs):
		"""Call AXITRA wrapper function suppressing noisy stdout/stderr if requested."""
		if not quiet:
			return fn(*args, **kwargs)

		# First layer: Python-level redirection.
		py_out = io.StringIO()
		py_err = io.StringIO()

		# Second layer: low-level FD redirection, needed for Fortran/C prints.
		devnull_fd = os.open(os.devnull, os.O_WRONLY)
		stdout_fd = os.dup(1)
		stderr_fd = os.dup(2)
		try:
			os.dup2(devnull_fd, 1)
			os.dup2(devnull_fd, 2)
			with redirect_stdout(py_out), redirect_stderr(py_err):
				return fn(*args, **kwargs)
		finally:
			os.dup2(stdout_fd, 1)
			os.dup2(stderr_fd, 2)
			os.close(stdout_fd)
			os.close(stderr_fd)
			os.close(devnull_fd)

	def model_array(self) -> np.ndarray:
		layers = self.cfg.velocity_model.layers
		return np.array(
			[[l.thickness, l.vp, l.vs, l.rho, l.qp, l.qs] for l in layers],
			dtype="float64",
		)

	def _latlon_to_local_xy_m(self, lat: float, lon: float) -> Tuple[float, float]:
		"""Convert station lat/lon to local Cartesian x/y in meters.

		Local convention matches geometry module:
		x -> North, y -> East, relative to source reference point.
		"""
		sta_east_m, sta_north_m = self._projection.latlon_to_xy(float(lat), float(lon))
		x_north_m = float(sta_north_m - self._src_north_m)
		y_east_m = float(sta_east_m - self._src_east_m)
		return x_north_m, y_east_m

	def stations_array(self, latlon: bool = True) -> np.ndarray:
		rows = []
		for i, st in enumerate(self.cfg.stations.stations, start=1):
			if latlon:
				rows.append([i, st.latitude, st.longitude, st.height])
			else:
				x_m, y_m = self._latlon_to_local_xy_m(float(st.latitude), float(st.longitude))
				rows.append([i, x_m, y_m, st.height])
		return np.array(rows, dtype="float64")

	def suggested_fmax_duration(self) -> Tuple[float, float]:
		duration = max(
			float(self.cfg.observed_data.t2 - self.cfg.observed_data.t1),
			float(self.cfg.observed_data.delta),
		)
		fmax = 1.0 / float(self.cfg.observed_data.delta) / 2 # Nyquist frequency
		return fmax, duration

	def build_geometry(
		self,
		slip_geom: float = 1.0,
	) -> FaultGeometry:
		return self.geometry_builder.build(
			slip_geom=slip_geom,
		)

	def apply_ellipse_model_to_geometry(
		self,
		geometry: FaultGeometry,
		model: np.ndarray,
	) -> FaultGeometry:
		"""Apply 7-parameter ellipse model to an existing geometry instance."""
		if len(model) < 7:
			raise ValueError("Model must include 7 parameters: a1,a2,theta,np,tp,dmax,vr")

		return self.ellipse_mapper.apply_to_geometry(
			geom=geometry,
			model=model,
		)

	def build_geometry_with_ellipse_slip(self, model: np.ndarray) -> FaultGeometry:
		"""
		Build fault geometry with ellipse-to-subfault slip mapping.

		Model parameters (7):
		  model[0]: a1 - semi-axis 1 (km)
		  model[1]: a2 - semi-axis 2 (km)
		  model[2]: theta - rotation angle (fraction of pi)
		  model[3]: np - position fraction [0,1]
		  model[4]: tp - position angle (fraction of 2pi)
		  model[5]: dmax - maximum slip (m)
		  model[6]: vr - rupture velocity (km/s)

		Returns:
			FaultGeometry with slip distribution applied from ellipse parameters.
		"""
		geom = self.build_geometry()
		return self.apply_ellipse_model_to_geometry(geometry=geom, model=model)

	def estimate_total_moment_and_mw(
		self,
		model: np.ndarray,
		geometry: Optional[FaultGeometry] = None,
	) -> Tuple[float, float]:
		"""Estimate total scalar moment (N.m) and Mw for a 7-parameter ellipse model."""
		geom = geometry if geometry is not None else self.build_geometry_with_ellipse_slip(model)
		m0 = float(geom.total_moment_nm())
		mw = float(geom.moment_magnitude_mw())
		return m0, mw

	def build_axitra(
		self,
		geometry: FaultGeometry,
		fmax: Optional[float] = None,
		duration: Optional[float] = None,
		xl: float = 0.0,
		latlon: bool = True,
		freesurface: bool = True,
		ikmax: int = 10000,
		id: Optional[int] = None,
		aw: Optional[float] = None,
	):
		Axitra, _ = self._import_axitra()
		auto_fmax, auto_duration = self.suggested_fmax_duration()
		return Axitra(
			self.model_array(),
			self.stations_array(latlon=latlon),
			geometry.to_axitra_sources(latlon=latlon),
			fmax=auto_fmax if fmax is None else float(fmax),
			duration=auto_duration if duration is None else float(duration),
			xl=float(xl),
			latlon=latlon,
			freesurface=freesurface,
			ikmax=int(ikmax),
			axpath=str(self.axitra_dir),
			id=id,
			aw=aw,
		)

	def green(self, ap, quiet: bool = True):
		_, moment = self._import_axitra()
		return self._call_with_optional_silence(moment.green, ap, quiet=quiet)

	def conv(
		self,
		ap,
		geometry: FaultGeometry,
		source_type: int = 1,
		t0: Optional[float] = None,
		t1: float = 0.0,
		unit: Optional[int] = None,
		sfunc=None,
		quiet: bool = True,
	):
		_, moment = self._import_axitra()
		t0_val = float(self.cfg.ellipse.t0) if t0 is None else float(t0)
		# AXITRA Ricker source (source_type=1) becomes degenerate for t0<=0.
		# Keep a small positive fallback to avoid null synthetics.
		if int(source_type) == 1 and t0_val <= 0.0:
			t0_val = max(float(self.cfg.observed_data.delta), 0.5)
		unit_val = int(self.cfg.observed_data.units) if unit is None else int(unit)
		hist = geometry.to_axitra_hist()
		return self._call_with_optional_silence(
			moment.conv,
			ap,
			hist,
			source_type=source_type,
			t0=t0_val,
			t1=t1,
			unit=unit_val,
			sfunc=sfunc,
			quiet=quiet,
		)


def precompute_greens_functions(params, v_model):
	"""
	Compatibility shim for legacy calls.
	"""
	_ = params, v_model
	return {
		"mode": "wrapper",
		"status": "use AxitraForwardModel.build_geometry -> build_axitra -> green",
	}

