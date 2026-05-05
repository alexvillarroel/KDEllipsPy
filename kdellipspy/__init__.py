# kdellipspy/__init__.py

# --- Configuración y Parámetros ---
from .config_parser import (
    ConfigParser, ObservedDataParams, SourcePosition, FaultPlaneParams,
    EllipseParams, InversionParam, InversionParams, InversionProcessParams,
    MomentTensor, StationParams, VelocityLayer, VelocityModel,
    read_input_ctl, parse_velocity_model, validate_input_ctl
)
# Renombramos Station para evitar colisión con geometry.py
from .config_parser import Station as ConfigStation 

# --- Modelado Foward ---
from .forward_model import (
    AxitraForwardModel, precompute_greens_functions
)

# --- Geometría y Cinemática ---
from .geometry import (
    UTMProjection, StationGeometry, Subfault, SourcePoint, FaultGeometry,
    GeometryBuilder, EllipticalSlipMapper, build_geometry_from_input_ctl,
    build_station_geometry, EllipseDiagnosticsResult, EllipseDiagnostics
)
# Renombramos Station de geometry
from .geometry import Station as GeometryStation 

# --- Gráficos ---
from .graphics_suite import (
    GraphicsConfig, GraphicsSuite
)

# --- Inversión: módulo base (compartido) ---
from .inversion_base import (
    NAModel,
    MisfitCalculator,
    NAResult,
    BaseInversionModel,
)

# --- Inversión: Neighbourhood Algorithm ---
from .inversion_na import (
    NAConfig,
    NAInversionModel,
)

# --- Inversión: MCMC (PyMC + ArviZ) ---
from .inversion_mcmc import (
    MCMCConfig,
    MCMCInversionModel,
)

# --- Procesamiento de Señales ---
from .signal_utils import (
    build_azi_times_array, write_azi_times_file,
    load_and_filter_observed_data, bandpass_filter_waveforms
)

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
    "build_station_geometry", "EllipseDiagnosticsResult", "EllipseDiagnostics",
    
    # Gráficos
    "GraphicsConfig", "GraphicsSuite",
    
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
    "load_and_filter_observed_data", "bandpass_filter_waveforms"
]