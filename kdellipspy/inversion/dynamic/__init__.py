"""Dynamic inversion adapter modules."""

from .config import DynamicLegacyRunConfig
from .legacy_bridge import run_dynamic_inversion
from .model_dynamic import DynamicInversionConfig, DynamicInversionModel

# fd3d_TSN-based pipeline (Fases 1-4, 2026-07-13) — supersedes the
# subprocess/bash bridge above, kept only for backwards compatibility.
from .dynamic_convolution import (
    bin_slip_rate_to_subfaults,
    convolve_dynamic_sources,
    project_rake,
    resample_to_axitra_grid,
)
from .forward_model_dynamic import DynamicForwardModel
from .model_dynamic_na import DynamicNAInversionModel
from .tsn_bridge import TSNRunConfig, read_tsn_fault_field, run_tsn_forward

__all__ = [
    "DynamicLegacyRunConfig",
    "run_dynamic_inversion",
    "DynamicInversionConfig",
    "DynamicInversionModel",
    "bin_slip_rate_to_subfaults",
    "convolve_dynamic_sources",
    "project_rake",
    "resample_to_axitra_grid",
    "DynamicForwardModel",
    "DynamicNAInversionModel",
    "TSNRunConfig",
    "read_tsn_fault_field",
    "run_tsn_forward",
]
