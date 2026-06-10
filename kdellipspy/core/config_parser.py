from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import re
import numpy as np
from pathlib import Path

@dataclass
class ObservedDataParams:
    """Section 1: Observed Data Parameters"""
    t1: float
    t2: float
    npts: int
    delta: float
    units: int  # 1: displacement, 2: velocity

    @classmethod
    def from_dict(cls, params: Dict) -> 'ObservedDataParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            t1=float(get('Time window start (t1)', 0.0)),
            t2=float(get('Time window end (t2)', 128.0)),
            npts=int(get('Number of points (Npts)', 512)),
            delta=float(get('Delta / Time step', 0.25)),
            units=int(get('Units', 1)) # Precise prefix
        )


@dataclass
class SourcePosition:
    """Section 2: Source Position & Focal Mechanism"""
    event_name: str
    origin_time: Optional[str]  # Origin Time (UTC)
    latitude: float
    longitude: float
    depth: float
    strike: float
    dip: float
    rake: float

    @classmethod
    def from_dict(cls, params: Dict) -> 'SourcePosition':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            event_name=str(get('Event Name', 'Unknown')),
            origin_time=get('Origin Time', None),
            latitude=float(get('Latitude', 0.0)),
            longitude=float(get('Longitude', 0.0)),
            depth=float(get('Depth', 0.0)),
            strike=float(get('Strike', 0.0)),
            dip=float(get('Dip', 0.0)),
            rake=float(get('Rake', 0.0))
        )


@dataclass
class FaultPlaneParams:
    """Section 3: Fault Plane Parameters"""
    lx: float  # Length along strike (m)
    ly: float  # Length along dip (m)
    hx: float  # Hypocenter position strike (m)
    hy: float  # Hypocenter position dip (m)
    nx: int    # Number of subfaults along strike
    ny: int    # Number of subfaults along dip

    @classmethod
    def from_dict(cls, params: Dict) -> 'FaultPlaneParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            lx=float(get('Length along strike (Lx)', 0.0)),
            ly=float(get('Length along dip (Ly)', 0.0)),
            hx=float(get('Hypocenter position strike (Hx)', 0.0)),
            hy=float(get('Hypocenter position dip (Hy)', 0.0)),
            nx=int(get('Number of subfaults along strike (Nx)', 1)),
            ny=int(get('Number of subfaults along dip (Ny)', 1))
        )


@dataclass
class EllipseParams:
    """Section 4: Ellipse Parameters & Frequency Band"""
    num_ellipses: int
    initial_slip: int  # 0: no, 1: yes
    slip_shape: int    # 0: constant, 1: gaussian, 2: ellipse
    freq1: float       # Hz
    freq2: float       # Hz
    t0: float          # Time shift (s)
    source_type: int   # Axitra source time function type
    # Fase del filtro pasabanda: True = acausal (filtfilt, fase cero),
    # False = causal (una pasada). Por defecto True para retrocompatibilidad.
    zerophase: bool = True

    @classmethod
    def from_dict(cls, params: Dict) -> 'EllipseParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            num_ellipses=int(get('Number of ellipses', 1)),
            initial_slip=int(get('Initial slip', 0)),
            slip_shape=int(get('Slip shape', 1)),
            freq1=float(get('Frequency 1 (Freq1)', 0.02)),
            freq2=float(get('Frequency 2 (Freq2)', 0.10)),
            t0=float(get('Time shift (T0)', 3.0)),
            source_type=int(get('Source type', params.get('source_type', 4))),
            # 1 = acausal (zero-phase), 0 = causal. Default 1 si no esta en input.ctl.
            zerophase=bool(int(get('Zerophase filter', 1))),
        )


@dataclass
class InversionParam:
    """Single parameter for inversion"""
    name: str
    min_val: float
    max_val: float
    flag: int  # 0: fixed, 1: invert

    def __repr__(self) -> str:
        status = "invert" if self.flag else "fixed"
        return f"{self.name}: [{self.min_val}, {self.max_val}] ({status})"


@dataclass
class InversionParams:
    """Section 5: Parameters to Invert"""
    parameters: List[InversionParam]

    @classmethod
    def from_dict(cls, params: Dict, param_lines: List[str]) -> 'InversionParams':
        params_list = []
        for line in param_lines:
            # Expected format: "Param N: <name text> : <min> <max> <flag>"
            # Use split by last colon for robust name capture
            label, sep, data = line.rpartition(':')
            if not sep: continue
            
            # Label is "Param N: Name text"
            # We want "Name text"
            name_part = label.split(':', 1)[-1].strip()
            
            parts = data.strip().split()
            if len(parts) >= 3:
                params_list.append(
                    InversionParam(
                        name=name_part,
                        min_val=float(parts[0]),
                        max_val=float(parts[1]),
                        flag=int(parts[2]),
                    )
                )
        
        return cls(parameters=params_list)


@dataclass
class InversionProcessParams:
    """Section 6: Inversion Process Parameters"""
    algorithm_type: int  # 0: NA, 1: MCMC (Metropolis-Hastings)
    num_iterations: int
    ss1: int             # Sample size for first iteration (NA) / unused default for MCMC
    ss_other: int        # Sample size for other iterations (NA)
    cells_resample: int
    misfit_time_window: float = 0.0 # 0=full signal, >0=window in seconds
    n_jobs: int = 1      # Parallel workers for NA/Forward evaluation
    # MCMC (only used when algorithm_type == 1); optional lines in input.ctl
    mcmc_total_steps: int = 500
    mcmc_burn_in: int = 0
    mcmc_proposal_scale: float = 0.08  # Gaussian RW std as fraction of (max-min) per parameter
    mcmc_thin: int = 1
    mcmc_chains: int = 1  # PyMC only; use 1 if axitra is not fork-safe

    @classmethod
    def from_dict(cls, params: Dict) -> 'InversionProcessParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            algorithm_type=int(get('Algorithm type', 0)),
            num_iterations=int(get('Number of iterations', 1)),
            ss1=int(get('Sample size for first iteration (SS1)', 30)),
            ss_other=int(get('Sample size for other iterations', 30)),
            cells_resample=int(get('Cells to resample', 7)),
            misfit_time_window=float(get('Misfit time window', 0.0)),
            n_jobs=int(get('Parallel workers', get('n_jobs', 1))),
            mcmc_total_steps=int(get('MCMC total steps', 500)),
            mcmc_burn_in=int(get('MCMC burn-in', 0)),
            mcmc_proposal_scale=float(get('MCMC proposal scale', 0.08)),
            mcmc_thin=int(get('MCMC thinning', 1)),
            mcmc_chains=int(get('MCMC chains', 1)),
        )


@dataclass
class MomentTensor:
    """Section 7: Moment Tensor (Full MT)"""
    flag: int
    mrr: float
    mtt: float
    mpp: float
    mrt: float
    mrp: float
    mtp: float
    exponent: float
    scaling_mode: str = "no_mt"

    @classmethod
    def from_dict(cls, params: Dict) -> 'MomentTensor':
        get = lambda key, default: _get_param_value(params, key, default)
        flag = int(get('Moment Tensor Flag', 0))
        raw_mode = str(get('MT Scaling Mode', '')).strip()
        scaling_mode = _normalize_mt_scaling_mode(flag=flag, raw_mode=raw_mode)
        return cls(
            flag=flag,
            mrr=float(get('Mrr', 0.0)),
            mtt=float(get('Mtt', 0.0)),
            mpp=float(get('Mpp', 0.0)),
            mrt=float(get('Mrt', 0.0)),
            mrp=float(get('Mrp', 0.0)),
            mtp=float(get('Mtp', 0.0)),
            exponent=float(get('Exponent (iexp)', 18.0)),
            scaling_mode=scaling_mode,
        )


@dataclass
class Station:
    """Individual station parameters"""
    latitude: float
    longitude: float
    height: float
    name: str
    use_n: bool = True
    use_e: bool = True 
    use_z: bool = True 

@dataclass
class StationParams:
    """Section 8: Station Parameters"""
    stations: List[Station] = field(default_factory=list)

    @classmethod
    def from_lines(cls, station_lines: List[str]) -> 'StationParams':
        stations = []
        for line in station_lines:
            parts = line.split()
            if len(parts) >= 4:
                use_n = use_e = use_z = True
                if len(parts) >= 7:
                    use_n = int(parts[4]) == 1
                    use_e = int(parts[5]) == 1
                    use_z = int(parts[6]) == 1
                stations.append(Station(
                    latitude=float(parts[0]),
                    longitude=float(parts[1]),
                    height=float(parts[2]),
                    name=parts[3],
                    use_n=use_n,
                    use_e=use_e,
                    use_z=use_z
                ))
        return cls(stations=stations)


@dataclass
class VelocityLayer:
    """Single velocity layer"""
    thickness: float
    vp: float
    vs: float
    rho: float
    qp: float
    qs: float


@dataclass
class VelocityModel:
    """Section 9: Velocity Model 1D"""
    layers: List[VelocityLayer] = field(default_factory=list)

    @classmethod
    def from_lines(cls, layer_lines: List[str]) -> 'VelocityModel':
        layers = []
        for line in layer_lines:
            parts = line.split()
            if len(parts) >= 6:
                layers.append(VelocityLayer(
                    thickness=float(parts[0]),
                    vp=float(parts[1]),
                    vs=float(parts[2]),
                    rho=float(parts[3]),
                    qp=float(parts[4]),
                    qs=float(parts[5])
                ))
        return cls(layers=layers)

    def to_numpy(self) -> np.ndarray:
        return np.array([
            [layer.thickness, layer.vp, layer.vs, layer.rho, layer.qp, layer.qs]
            for layer in self.layers
        ])


class ConfigParser:
    """Parses and stores the inversion configuration from the 'input.ctl' file."""

    @staticmethod
    def get_base_template_path() -> Path:
        return Path(__file__).parent.parent / "io" / "templates" / "input.ctl.base"

    def __init__(self, filepath: Optional[str | Path] = None):
        if filepath is None:
            filepath = "<manual>"
        self.filepath = str(filepath)
        self.observed_data = None
        self.source_position = None
        self.fault_plane = None
        self.ellipse = None
        self.inversion_params = None
        self.inversion_process = None
        self.moment_tensor = None
        self.stations = None
        self.velocity_model = None
        
        if self.filepath not in ("<manual>", "<from_dict>"):
            self.parse()

    def plot_stations(self, show: bool = True, save_path: Optional[str] = None, use_cartopy: bool = True) -> Tuple[Any, Any]:
        """
        Grafica la ubicación de las estaciones e hipocentro.
        """
        from .plotting import plot_stations_map
        return plot_stations_map(self, show=show, save_path=save_path, use_cartopy=use_cartopy)

    @classmethod
    def from_dict(cls, params: Dict[str, Any]) -> 'ConfigParser':
        """Create a ConfigParser instance from a dictionary of parameters."""
        instance = cls("<from_dict>")
        
        # Section 1
        instance.observed_data = ObservedDataParams.from_dict(params.get('observed_data', {}))
        # Section 2
        instance.source_position = SourcePosition.from_dict(params.get('source_position', {}))
        # Section 3
        instance.fault_plane = FaultPlaneParams.from_dict(params.get('fault_plane', {}))
        # Section 4
        instance.ellipse = EllipseParams.from_dict(params.get('ellipse', {}))
        # Section 5
        # InversionParams needs a list of lines for its from_dict, but we can bypass it
        # or implement a simpler from_dict for it.
        # For now, let's assume it's passed as a list of InversionParam objects or similar.
        inv_data = params.get('inversion_params', {})
        if isinstance(inv_data, list):
            instance.inversion_params = InversionParams(parameters=inv_data)
        else:
            # Try to handle it if it's a dict
            instance.inversion_params = InversionParams(parameters=[])
            
        # Section 6
        instance.inversion_process = InversionProcessParams.from_dict(params.get('inversion_process', {}))
        # Section 7
        instance.moment_tensor = MomentTensor.from_dict(params.get('moment_tensor', {}))
        
        # Section 8
        stations_data = params.get('stations', [])
        if isinstance(stations_data, list):
            stations_list = []
            for st in stations_data:
                if isinstance(st, dict):
                    stations_list.append(Station(
                        latitude=float(st.get('latitude', 0.0)),
                        longitude=float(st.get('longitude', 0.0)),
                        height=float(st.get('height', 0.0)),
                        name=str(st.get('name', 'ST')),
                        use_n=bool(st.get('use_n', True)),
                        use_e=bool(st.get('use_e', True)),
                        use_z=bool(st.get('use_z', True))
                    ))
                elif isinstance(st, Station):
                    stations_list.append(st)
            instance.stations = StationParams(stations=stations_list)
        
        # Section 9
        vel_data = params.get('velocity_model', [])
        if isinstance(vel_data, list):
            layers = []
            for l in vel_data:
                if isinstance(l, dict):
                    layers.append(VelocityLayer(
                        thickness=float(l.get('thickness', 0.0)),
                        vp=float(l.get('vp', 0.0)),
                        vs=float(l.get('vs', 0.0)),
                        rho=float(l.get('rho', 0.0)),
                        qp=float(l.get('qp', 0.0)),
                        qs=float(l.get('qs', 0.0))
                    ))
                elif isinstance(l, VelocityLayer):
                    layers.append(l)
            instance.velocity_model = VelocityModel(layers=layers)
            
        return instance

    def update_stations(self, keep_station_names: list[str]) -> None:
        if not self.stations or not self.stations.stations:
            return
        keep_upper = [s.strip().upper() for s in keep_station_names]
        self.stations.stations = [s for s in self.stations.stations if s.name.upper() in keep_upper]
        if self.filepath not in ("<manual>", "<from_dict>"):
            self.save()

    def update_stations_from_sac(self, raw_dir: str | Path, station_names: list[str]) -> None:
        try: from obspy import read
        except ImportError: return
        if not self.stations: self.stations = StationParams()
        raw_path = Path(raw_dir)
        sac_files = list(raw_path.glob("*.SAC")) + list(raw_path.glob("*.sac"))
        if not sac_files: return
        coords_map = {}
        for f in sac_files:
            try:
                st = read(str(f), headonly=True)
                tr = st[0]
                coords_map[tr.stats.station.strip().upper()] = (float(tr.stats.sac.stla), float(tr.stats.sac.stlo))
            except Exception: pass
        new_stations = []
        for name in station_names:
            nu = name.strip().upper()
            if nu in coords_map:
                lat, lon = coords_map[nu]
                new_stations.append(Station(latitude=lat, longitude=lon, height=0.0, name=nu))
        self.stations.stations = new_stations
        if self.filepath not in ("<manual>", "<from_dict>"):
            self.save()

    def save(self, output_path: Optional[str | Path] = None) -> None:
        path = Path(output_path if output_path else self.filepath)
        if str(path) in ("<manual>", "<from_dict>"):
            self._save_from_scratch(path)
            return

        with open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        simple_params = self._get_simple_params_map()
        new_lines = []
        current_section = 0
        skip_data = False
        COLON_POS = 48 # Standard alignment
        VALUE_WIDTH = 25

        for line in lines:
            stripped = line.strip()
            m_header = re.match(r"^#\s*(\d+)\.", stripped)
            if m_header:
                current_section = int(m_header.group(1))
                skip_data = False
                new_lines.append(line)
                continue
            
            if skip_data:
                if current_section in (5, 8, 9) and stripped and not stripped.startswith("#") and ":" not in stripped:
                    continue
                elif current_section == 5 and stripped.startswith("Param"):
                    continue
                else:
                    skip_data = False

            if ":" in line and not skip_data:
                label, value = _split_line_by_separator(line)
                norm_key = _normalize_key(label)
                
                matched_key = None
                if norm_key in simple_params:
                    matched_key = norm_key
                else:
                    for sk in simple_params:
                        if norm_key.startswith(sk) or sk.startswith(norm_key):
                            matched_key = sk
                            break
                
                if matched_key:
                    val = simple_params[matched_key]
                    # Fixed alignment logic
                    new_line = label.ljust(COLON_POS-1) + ":" + val.rjust(VALUE_WIDTH) + "\n"
                    new_lines.append(new_line)
                    
                    if matched_key in ("number of parameters", "number of stations", "number of layers"):
                        self._inject_list_data(new_lines, current_section)
                        skip_data = True
                    continue

            # Handle list data that doesn't have a "Number of" line
            if current_section in (5, 8, 9) and stripped and not stripped.startswith("#") and not skip_data:
                is_data = (current_section == 5 and stripped.startswith("Param")) or (":" not in stripped)
                if is_data:
                    self._inject_list_data(new_lines, current_section)
                    skip_data = True
                    continue

            new_lines.append(line)

        with open(path, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        print(f"✓ Configuration surgically updated: {path}")

    def _get_simple_params_map(self) -> Dict[str, str]:
        p = {}
        # Section 1
        p["time window start (t1)"] = f"{self.observed_data.t1:.6f}"
        p["time window end (t2)"] = f"{self.observed_data.t2:.6f}"
        p["number of points (npts)"] = str(self.observed_data.npts)
        p["delta / time step"] = f"{self.observed_data.delta:.6f}"
        p["units"] = str(self.observed_data.units)
        # Section 2
        p["event name"] = self.source_position.event_name
        p["origin time"] = self.source_position.origin_time if self.source_position.origin_time else ""
        p["latitude"] = f"{self.source_position.latitude:.6f}"
        p["longitude"] = f"{self.source_position.longitude:.6f}"
        p["depth"] = f"{self.source_position.depth:.6f}"
        p["strike"] = f"{self.source_position.strike:.6f}"
        p["dip"] = f"{self.source_position.dip:.6f}"
        p["rake"] = f"{self.source_position.rake:.6f}"
        # Section 3
        p["length along strike (lx)"] = f"{self.fault_plane.lx:.6f}"
        p["length along dip (ly)"] = f"{self.fault_plane.ly:.6f}"
        p["hypocenter position strike (hx)"] = f"{self.fault_plane.hx:.6f}"
        p["hypocenter position dip (hy)"] = f"{self.fault_plane.hy:.6f}"
        p["number of subfaults along strike (nx)"] = str(self.fault_plane.nx)
        p["number of subfaults along dip (ny)"] = str(self.fault_plane.ny)
        # Section 4
        p["number of ellipses"] = str(self.ellipse.num_ellipses)
        p["initial slip"] = str(self.ellipse.initial_slip)
        p["slip shape"] = str(self.ellipse.slip_shape)
        p["frequency 1 (freq1)"] = f"{self.ellipse.freq1:.6f}"
        p["frequency 2 (freq2)"] = f"{self.ellipse.freq2:.6f}"
        p["time shift (t0)"] = f"{self.ellipse.t0:.6f}"
        p["source type"] = str(self.ellipse.source_type)
        p["zerophase filter"] = str(int(self.ellipse.zerophase))
        # Section 6
        p["algorithm type"] = str(self.inversion_process.algorithm_type)
        p["number of iterations"] = str(self.inversion_process.num_iterations)
        p["sample size for first iteration (ss1)"] = str(self.inversion_process.ss1)
        p["sample size for other iterations"] = str(self.inversion_process.ss_other)
        p["cells to resample"] = str(self.inversion_process.cells_resample)
        p["misfit time window"] = f"{self.inversion_process.misfit_time_window:.6f}"
        p["parallel workers (n_jobs)"] = str(self.inversion_process.n_jobs)
        p["mcmc total steps"] = str(self.inversion_process.mcmc_total_steps)
        p["mcmc burn-in"] = str(self.inversion_process.mcmc_burn_in)
        p["mcmc proposal scale"] = f"{self.inversion_process.mcmc_proposal_scale:.6f}"
        p["mcmc thinning"] = str(self.inversion_process.mcmc_thin)
        p["mcmc chains"] = str(self.inversion_process.mcmc_chains)
        # Section 7
        p["moment tensor flag"] = str(self.moment_tensor.flag)
        p["mt scaling mode"] = self.moment_tensor.scaling_mode
        p["mrr"] = f"{self.moment_tensor.mrr:.6e}"; p["mtt"] = f"{self.moment_tensor.mtt:.6e}"; p["mpp"] = f"{self.moment_tensor.mpp:.6e}"
        p["mrt"] = f"{self.moment_tensor.mrt:.6e}"; p["mrp"] = f"{self.moment_tensor.mrp:.6e}"; p["mtp"] = f"{self.moment_tensor.mtp:.6e}"
        p["exponent (iexp)"] = str(self.moment_tensor.exponent)
        # Counts
        p["number of parameters"] = str(len(self.inversion_params.parameters))
        p["number of stations"] = str(len(self.stations.stations))
        p["number of layers"] = str(len(self.velocity_model.layers))
        return p

    def _inject_list_data(self, lines: List[str], section: int):
        if section == 5:
            for i, p in enumerate(self.inversion_params.parameters, 1):
                lines.append(f" Param {i} : {p.name.ljust(30)} : {p.min_val:10.6f} {p.max_val:10.6f} {p.flag:3}\n")
        elif section == 8:
            for st in self.stations.stations:
                un, ue, uz = (1 if st.use_n else 0, 1 if st.use_e else 0, 1 if st.use_z else 0)
                lines.append(f" {st.latitude:10.6f} {st.longitude:10.6f} {st.height:10.6f} {st.name:8} {un} {ue} {uz}\n")
        elif section == 9:
            for l in self.velocity_model.layers:
                lines.append(f" {l.thickness:10.6f} {l.vp:10.6f} {l.vs:10.6f} {l.rho:10.6f} {l.qp:10.6f} {l.qs:10.6f}\n")

    def parse(self) -> None:
        with open(self.filepath, 'r') as f: content = f.read()
        sections = self._split_sections(content)
        if 1 in sections: self.observed_data = ObservedDataParams.from_dict(self._extract_params(sections[1]))
        if 2 in sections: self.source_position = SourcePosition.from_dict(self._extract_params(sections[2]))
        if 3 in sections: self.fault_plane = FaultPlaneParams.from_dict(self._extract_params(sections[3]))
        if 4 in sections: self.ellipse = EllipseParams.from_dict(self._extract_params(sections[4]))
        if 5 in sections:
            p_dict = self._extract_params(sections[5])
            p_lines = [l.strip() for l in sections[5].split('\n') if l.strip().startswith('Param')]
            self.inversion_params = InversionParams.from_dict(p_dict, p_lines)
        if 6 in sections: self.inversion_process = InversionProcessParams.from_dict(self._extract_params(sections[6]))
        if 7 in sections: self.moment_tensor = MomentTensor.from_dict(self._extract_params(sections[7]))
        if 8 in sections: self.stations = StationParams.from_lines(self._extract_data_lines(sections[8]))
        if 9 in sections: self.velocity_model = VelocityModel.from_lines(self._extract_data_lines(sections[9]))

    def _split_sections(self, content: str) -> Dict[int, str]:
        sections = {}; cur = None; chunk = []
        for line in content.splitlines():
            m = re.match(r"^#\s*(\d+)\.", line)
            if m:
                if cur is not None: sections[cur] = "\n".join(chunk)
                cur = int(m.group(1)); chunk = []
                continue
            if cur is not None: chunk.append(line)
        if cur is not None: sections[cur] = "\n".join(chunk)
        return sections

    def _extract_params(self, section: str) -> Dict[str, str]:
        params = {}
        for line in section.split('\n'):
            line = line.strip()
            if ':' in line:
                label, val = _split_line_by_separator(line)
                if label: params[label.strip()] = val.strip()
        return params

    def _extract_data_lines(self, section: str) -> List[str]:
        return [l.strip() for l in section.split('\n') if l.strip() and not l.strip().startswith('#') and ":" not in l]

    def _save_from_scratch(self, path: Path):
        # Fallback if needed, but not really expected here
        pass

def _normalize_key(text: str) -> str:
    text = text.lower().strip().replace(":", "")
    return re.sub(r"\s+", " ", text)

def _get_param_value(params: Dict[str, str], key: str, default: Any) -> Any:
    wanted = _normalize_key(key)
    # 1. Try exact match or startswith (standard)
    for k, v in params.items():
        norm_k = _normalize_key(k)
        if norm_k == wanted or norm_k.startswith(wanted):
            return v
    
    # 2. Try 'contains' match for cases like "NA: Cells to resample"
    for k, v in params.items():
        if wanted in _normalize_key(k):
            return v
            
    return default

def _split_line_by_separator(line: str) -> Tuple[str, str]:
    """Choose the separator colon based on surrounding whitespace (prefer LAST surrounded)."""
    indices = [i for i, c in enumerate(line) if c == ':']
    if not indices: return line, ""
    
    # Standard choice: find colons surrounded by at least one space
    # (e.g. " ) : " or " (km) : ")
    spaced_indices = []
    for idx in indices:
        has_left = (idx > 0 and line[idx-1].isspace())
        has_right = (idx < len(line)-1 and line[idx+1].isspace())
        if has_left and has_right:
            spaced_indices.append(idx)
            
    if spaced_indices:
        # Use the LAST spaced colon as separator (handles labels with internal spaced colons)
        best_idx = spaced_indices[-1]
        return line[:best_idx], line[best_idx+1:]
        
    # Fallback 1: last colon NOT part of a time string \d\d:\d\d
    for idx in reversed(indices):
        is_time = False
        if idx > 1 and idx < len(line)-2:
            if line[idx-2:idx].isdigit() and line[idx+1:idx+3].isdigit():
                is_time = True
        if not is_time:
            return line[:idx], line[idx+1:]
            
    # Final fallback: just use the first colon
    return line[:indices[0]], line[indices[0]+1:]

def _normalize_mt_scaling_mode(flag: int, raw_mode: str) -> str:
    if int(flag) == 0: return "no_mt"
    m = str(raw_mode).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {"":"mt_factored", "mt":"mt_factored", "full_mt":"mt_factored", "factored":"mt_factored", "strict":"mt_strict", "fixed_m0":"mt_strict", "no_mt":"no_mt"}
    return aliases.get(m, "mt_factored")

def read_input_ctl(filepath: str) -> Dict[str, Any]:
    """Compatibility helper returning a flat dictionary used by the pipeline."""
    cfg = ConfigParser(filepath)

    out: Dict[str, Any] = {
        't1': cfg.observed_data.t1,
        't2': cfg.observed_data.t2,
        'npts': cfg.observed_data.npts,
        'delta': cfg.observed_data.delta,
        'units': cfg.observed_data.units,
        'event_name': cfg.source_position.event_name,
        'origin_time': cfg.source_position.origin_time,
        'lat': cfg.source_position.latitude,
        'lon': cfg.source_position.longitude,
        'depth': cfg.source_position.depth,
        'strike': cfg.source_position.strike,
        'dip': cfg.source_position.dip,
        'rake': cfg.source_position.rake,
        'lx': cfg.fault_plane.lx,
        'ly': cfg.fault_plane.ly,
        'hx': cfg.fault_plane.hx,
        'hy': cfg.fault_plane.hy,
        'nx': cfg.fault_plane.nx,
        'ny': cfg.fault_plane.ny,
        'num_ellipses': cfg.ellipse.num_ellipses,
        'initial_slip': cfg.ellipse.initial_slip,
        'slip_shape': cfg.ellipse.slip_shape,
        'freq1': cfg.ellipse.freq1,
        'freq2': cfg.ellipse.freq2,
        't0': cfg.ellipse.t0,
        'source_type': cfg.ellipse.source_type,
        'algorithm_type': cfg.inversion_process.algorithm_type,
        'num_iterations': cfg.inversion_process.num_iterations,
        'ss1': cfg.inversion_process.ss1,
        'ss_other': cfg.inversion_process.ss_other,
        'cells_resample': cfg.inversion_process.cells_resample,
        'misfit_time_window': cfg.inversion_process.misfit_time_window,
        'mcmc_total_steps': cfg.inversion_process.mcmc_total_steps,
        'mcmc_burn_in': cfg.inversion_process.mcmc_burn_in,
        'mcmc_proposal_scale': cfg.inversion_process.mcmc_proposal_scale,
        'mcmc_thin': cfg.inversion_process.mcmc_thin,
        'mcmc_chains': cfg.inversion_process.mcmc_chains,
        'moment_tensor_flag': cfg.moment_tensor.flag,
        'mrr': cfg.moment_tensor.mrr,
        'mtt': cfg.moment_tensor.mtt,
        'mpp': cfg.moment_tensor.mpp,
        'mrt': cfg.moment_tensor.mrt,
        'mrp': cfg.moment_tensor.mrp,
        'mtp': cfg.moment_tensor.mtp,
        'iexp': cfg.moment_tensor.exponent,
        'mt_scaling_mode': cfg.moment_tensor.scaling_mode,
        'inversion_params': cfg.inversion_params.parameters,
        'stations': cfg.stations.stations,
        'velocity_layers': cfg.velocity_model.layers,
    }
    return out


def parse_velocity_model(filepath: str) -> List[Dict[str, float]]:
    """Compatibility helper for existing code path."""
    cfg = ConfigParser(filepath)
    return [
        {
            'thickness': layer.thickness,
            'vp': layer.vp,
            'vs': layer.vs,
            'rho': layer.rho,
            'qp': layer.qp,
            'qs': layer.qs,
        }
        for layer in cfg.velocity_model.layers
    ]


def validate_input_ctl(filepath: str) -> Dict[str, Any]:
    """Return a concise parse summary for smoke testing."""
    cfg = ConfigParser(filepath)
    return {
        'sections_ok': all(
            x is not None
            for x in [
                cfg.observed_data,
                cfg.source_position,
                cfg.fault_plane,
                cfg.ellipse,
                cfg.inversion_params,
                cfg.inversion_process,
                cfg.moment_tensor,
                cfg.stations,
                cfg.velocity_model,
            ]
        ),
        'n_inversion_params': len(cfg.inversion_params.parameters),
        'n_stations': len(cfg.stations.stations),
        'n_layers': len(cfg.velocity_model.layers),
        'freq_band': (cfg.ellipse.freq1, cfg.ellipse.freq2),
    }


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:
        print('Usage: python -m src.config_parser <path_to_input.ctl>')
        raise SystemExit(2)

    summary = validate_input_ctl(sys.argv[1])
    print('Validation summary:', summary)
