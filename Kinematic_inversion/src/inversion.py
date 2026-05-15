"""Compatibility wrapper for the new neighpy-based inversion pipeline.

Legacy code that imported src.inversion can keep working, but the real
implementation now lives in src.inversion_na.
"""

from .inversion_na import NAConfig, NAModel, MisfitCalculator, NAResult, NAInversionModel

__all__ = [
    "NAConfig",
    "NAModel",
    "MisfitCalculator",
    "NAResult",
    "NAInversionModel",
]