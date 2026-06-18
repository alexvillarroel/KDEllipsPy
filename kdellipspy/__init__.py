# kdellipspy/__init__.py

# --- Configuración y Parámetros ---
from .core.config_parser import (
    ConfigParser, ObservedDataParams, SourcePosition, FaultPlaneParams,
    EllipseParams, InversionParam, InversionParams, InversionProcessParams,
    MomentTensor, StationParams, VelocityLayer, VelocityModel,
    read_input_ctl, parse_velocity_model, validate_input_ctl
)
# Renombramos Station para evitar colisión con geometry.py
from .core.config_parser import Station as ConfigStation

# --- Procesamiento de señales sísmicas (SAC/MiniSEED) ---
from .core.sac_processor import (
    ProcessingConfig,
    WaveformExtractor,
    process_trace,
    acausal_bandpass,
    polynomial_baseline_correction,
)

# --- Carga de datos RAW para inversion ---
from .core.raw_loader import (
    load_observed_from_raw,
)

# --- Modelado Foward ---
from .core.forward_model import (
    AxitraForwardModel, precompute_greens_functions
)

# --- Geometría y Cinemática ---
from .core.geometry import (
    UTMProjection, StationGeometry, Subfault, SourcePoint, FaultGeometry,
    GeometryBuilder, EllipticalSlipMapper, build_geometry_from_input_ctl,
    build_station_geometry
)
# Renombramos Station de geometry
from .core.geometry import Station as GeometryStation

# --- Gráficos ---
from .core import plotting

# --- Inversión: módulo base (compartido) ---
from .inversion.base import (
    NAModel,
    MisfitCalculator,
    NAResult,
    BaseInversionModel,
)

# --- Inversión: Neighbourhood Algorithm ---
from .inversion.kinematic.model_na import (
    NAConfig,
    NAInversionModel,
)

# --- Inversión: MCMC (PyMC + ArviZ) ---
from .inversion.kinematic.model_mcmc import (
    MCMCConfig,
    MCMCInversionModel,
)

# --- Procesamiento de Señales ---
from .core.signal_utils import (
    build_azi_times_array, write_azi_times_file,
    load_and_filter_observed_data, bandpass_filter_waveforms,
    integrate_waveforms
)
from .core.data_preprocessor import DataPreprocessor

# Definimos exactamente qué se expone al exterior
__all__ = [
    # Config
    "ConfigParser", "ObservedDataParams", "SourcePosition", "FaultPlaneParams",
    "EllipseParams", "InversionParam", "InversionParams", "InversionProcessParams",
    "MomentTensor", "StationParams", "VelocityLayer", "VelocityModel",
    "ConfigStation", "read_input_ctl", "parse_velocity_model", "validate_input_ctl",
    
    # Forward Model
    "AxitraForwardModel", "precompute_greens_functions",
    
    # Geometría
    "UTMProjection", "GeometryStation", "StationGeometry", "Subfault", "SourcePoint", 
    "FaultGeometry", "GeometryBuilder", "EllipticalSlipMapper", "build_geometry_from_input_ctl",
    "build_station_geometry",
    
    # Gráficos
    "plotting",
    
    # Inversión base
    "NAModel",
    "MisfitCalculator",
    "NAResult",
    "BaseInversionModel",
    # NA
    "NAConfig",
    "NAInversionModel",
    # MCMC
    "MCMCConfig",
    "MCMCInversionModel",

    # Procesamiento de Señales
    "build_azi_times_array", "write_azi_times_file",
    "load_and_filter_observed_data", "bandpass_filter_waveforms", "integrate_waveforms",
    "DataPreprocessor"
]
