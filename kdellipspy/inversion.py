"""Compatibility wrapper for the split inversion pipeline.

Legacy code that imported from ``kdellipspy.inversion`` can keep working.
The real implementations now live in:
  - inversion_base  : NAModel, MisfitCalculator, NAResult, BaseInversionModel
  - inversion_na    : NAConfig, NAInversionModel
  - inversion_mcmc  : MCMCConfig, MCMCInversionModel
"""

from .inversion_base import (
    NAModel,
    MisfitCalculator,
    NAResult,
    BaseInversionModel,
)
from .inversion_na import (
    NAConfig,
    NAInversionModel,
)
from .inversion_mcmc import (
    MCMCConfig,
    MCMCInversionModel,
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
]