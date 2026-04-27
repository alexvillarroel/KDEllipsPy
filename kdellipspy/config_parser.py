from dataclasses import dataclass
from typing import Any, Dict, List
import re
import numpy as np


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
            units=int(get('Units (1:disp, 2:vel)', 1))
        )


@dataclass
class SourcePosition:
    """Section 2: Source Position & Focal Mechanism"""
    event_name: str
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

    @classmethod
    def from_dict(cls, params: Dict) -> 'EllipseParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            num_ellipses=int(get('Number of ellipses', 1)),
            initial_slip=int(get('Initial slip', 0)),
            slip_shape=int(get('Slip shape', 1)),
            freq1=float(get('Frequency 1 (Freq1)', 0.02)),
            freq2=float(get('Frequency 2 (Freq2)', 0.10)),
            t0=float(get('Time shift (T0)', 3.0))
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
            m = re.match(r"^\s*Param\s+\d+\s*:\s*(.*?)\s*:\s*(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(\d+)\s*$", line)
            if not m:
                continue
            params_list.append(
                InversionParam(
                    name=m.group(1).strip(),
                    min_val=float(m.group(2)),
                    max_val=float(m.group(3)),
                    flag=int(m.group(4)),
                )
            )
        
        return cls(parameters=params_list)


@dataclass
class InversionProcessParams:
    """Section 6: Inversion Process Parameters"""
    algorithm_type: int  # 0: NA, 1: MC
    num_iterations: int
    ss1: int             # Sample size for first iteration
    ss_other: int        # Sample size for other iterations
    cells_resample: int

    @classmethod
    def from_dict(cls, params: Dict) -> 'InversionProcessParams':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            algorithm_type=int(get('Algorithm type', 0)),
            num_iterations=int(get('Number of iterations', 1)),
            ss1=int(get('Sample size for first iteration (SS1)', 30)),
            ss_other=int(get('Sample size for other iterations', 30)),
            cells_resample=int(get('Cells to resample', 7))
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

    @classmethod
    def from_dict(cls, params: Dict) -> 'MomentTensor':
        get = lambda key, default: _get_param_value(params, key, default)
        return cls(
            flag=int(get('Moment Tensor Flag', 0)),
            mrr=float(get('Mrr', 0.0)),
            mtt=float(get('Mtt', 0.0)),
            mpp=float(get('Mpp', 0.0)),
            mrt=float(get('Mrt', 0.0)),
            mrp=float(get('Mrp', 0.0)),
            mtp=float(get('Mtp', 0.0)),
            exponent=float(get('Exponent (iexp)', 18.0))
        )


@dataclass
class Station:
    """Individual station parameters"""
    latitude: float
    longitude: float
    height: float
    name: str


@dataclass
class StationParams:
    """Section 8: Station Parameters"""
    stations: List[Station]

    @classmethod
    def from_lines(cls, station_lines: List[str]) -> 'StationParams':
        stations = []
        for line in station_lines:
            parts = line.split()
            if len(parts) >= 4:
                stations.append(Station(
                    latitude=float(parts[0]),
                    longitude=float(parts[1]),
                    height=float(parts[2]),
                    name=parts[3]
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
    layers: List[VelocityLayer]

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
        """
        Exporta las capas directamente como una matriz de NumPy (Nx6).
        Este es exactamente el formato que requiere el wrapper de Axitra.
        """
        return np.array([
            [layer.thickness, layer.vp, layer.vs, layer.rho, layer.qp, layer.qs]
            for layer in self.layers
        ])


class ConfigParser:
    """
    Parses and stores the inversion configuration from the 'input.ctl' file.
    (Parsea y almacena la configuración de la inversión desde el archivo 'input.ctl'.)
    
    This class is responsible for reading seismic event parameters, station grids, 
    filtering frequencies, and Neighbourhood Algorithm (NA) parameters, structuring 
    them for use throughout the pipeline.
    (Esta clase es responsable de leer los parámetros del evento sísmico, la grilla 
    de estaciones, las frecuencias de filtrado y los parámetros del algoritmo NA, 
    estructurándolos para su uso en el pipeline.)

    Attributes:
        filepath (str or Path): Path to the 'input.ctl' control file. 
                                (Ruta al archivo de control 'input.ctl'.)
        source_position (SourceConfig): Seismic source parameters. 
                                        (Parámetros de la fuente sísmica.)
        stations (StationConfig): List and information of stations. 
                                  (Lista e información de las estaciones.)
        ellipse (EllipseConfig): Ellipse geometry parameters and frequencies. 
                                 (Parámetros geométricos de la elipse y frecuencias.)
        inversion_process (InversionConfig): NA algorithm parameters. 
                                             (Parámetros del algoritmo de inversión NA.)
    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.observed_data = None
        self.source_position = None
        self.fault_plane = None
        self.ellipse = None
        self.inversion_params = None
        self.inversion_process = None
        self.moment_tensor = None
        self.stations = None
        self.velocity_model = None
        
        self.parse()

    def parse(self) -> None:
        """Parse the entire control file"""
        with open(self.filepath, 'r') as f:
            content = f.read()

        sections = self._split_sections(content)
        self._parse_all_sections(sections)

    def _split_sections(self, content: str) -> Dict[int, str]:
        """Split content into 9 sections"""
        sections: Dict[int, str] = {}
        current_section = None
        chunk: List[str] = []

        for line in content.splitlines():
            m = re.match(r"^#\s*(\d+)\.", line)
            if m:
                if current_section is not None:
                    sections[current_section] = "\n".join(chunk)
                current_section = int(m.group(1))
                chunk = []
                continue
            if current_section is not None:
                chunk.append(line)

        if current_section is not None:
            sections[current_section] = "\n".join(chunk)

        return sections

    def _parse_all_sections(self, sections: Dict[int, str]) -> None:
        """Parse each section"""
        if 1 in sections:
            self.observed_data = ObservedDataParams.from_dict(
                self._extract_params(sections[1])
            )
        if 2 in sections:
            self.source_position = SourcePosition.from_dict(
                self._extract_params(sections[2])
            )
        if 3 in sections:
            self.fault_plane = FaultPlaneParams.from_dict(
                self._extract_params(sections[3])
            )
        if 4 in sections:
            self.ellipse = EllipseParams.from_dict(
                self._extract_params(sections[4])
            )
        if 5 in sections:
            params_dict = self._extract_params(sections[5])
            param_lines = self._extract_param_lines(sections[5])
            self.inversion_params = InversionParams.from_dict(params_dict, param_lines)
        if 6 in sections:
            self.inversion_process = InversionProcessParams.from_dict(
                self._extract_params(sections[6])
            )
        if 7 in sections:
            self.moment_tensor = MomentTensor.from_dict(
                self._extract_params(sections[7])
            )
        if 8 in sections:
            station_lines = self._extract_data_lines(sections[8])
            self.stations = StationParams.from_lines(station_lines)
        if 9 in sections:
            layer_lines = self._extract_data_lines(sections[9])
            self.velocity_model = VelocityModel.from_lines(layer_lines)

    def _extract_params(self, section: str) -> Dict[str, str]:
        """Extract key-value pairs from section"""
        params = {}
        for line in section.split('\n'):
            if ':' in line:
                # Use right split because labels can include ':' (e.g. 'Param 1: ... : val')
                parts = line.rsplit(':', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    if key and value:
                        params[key] = value
        return params

    def _extract_param_lines(self, section: str) -> List[str]:
        """Extract lines starting with 'Param'"""
        lines = []
        for line in section.split('\n'):
            if line.strip().startswith('Param'):
                lines.append(line.strip())
        return lines

    def _extract_data_lines(self, section: str) -> List[str]:
        """Extract data lines (non-header, non-empty)"""
        lines = []
        for line in section.split('\n'):
            line = line.strip()
            if line and not line.startswith('#') and ':' not in line:
                lines.append(line)
        return lines

    def __repr__(self) -> str:
        return f"ConfigParser({self.filepath})"


def _normalize_key(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"\s+", " ", text)
    return text


def _get_param_value(params: Dict[str, str], key: str, default: Any) -> Any:
    """Lookup by normalized prefix so parser is tolerant to spacing/units formatting."""
    wanted = _normalize_key(key)
    for k, v in params.items():
        nk = _normalize_key(k)
        if nk.startswith(wanted):
            return v
    return default


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
        'algorithm_type': cfg.inversion_process.algorithm_type,
        'num_iterations': cfg.inversion_process.num_iterations,
        'ss1': cfg.inversion_process.ss1,
        'ss_other': cfg.inversion_process.ss_other,
        'cells_resample': cfg.inversion_process.cells_resample,
        'moment_tensor_flag': cfg.moment_tensor.flag,
        'mrr': cfg.moment_tensor.mrr,
        'mtt': cfg.moment_tensor.mtt,
        'mpp': cfg.moment_tensor.mpp,
        'mrt': cfg.moment_tensor.mrt,
        'mrp': cfg.moment_tensor.mrp,
        'mtp': cfg.moment_tensor.mtp,
        'iexp': cfg.moment_tensor.exponent,
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