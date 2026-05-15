"""Compatibility layer and public exports for inversion models."""

from .base import (
    NAModel,
    MisfitCalculator,
    NAResult,
    BaseInversionModel,
)
from .kinematic.model_na import (
    NAConfig,
    NAInversionModel,
)
from .kinematic.model_mcmc import (
    MCMCConfig,
    MCMCInversionModel,
)
from .dynamic.model_dynamic import (
    DynamicInversionConfig,
    DynamicInversionModel,
)

__all__ = [
    "NAModel",
    "MisfitCalculator",
    "NAResult",
    "BaseInversionModel",
    "NAConfig",
    "NAInversionModel",
    "MCMCConfig",
    "MCMCInversionModel",
    "DynamicInversionConfig",
    "DynamicInversionModel",
]
