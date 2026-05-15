"""Input parsing for kinematic and dynamic modes."""

from .kinematic_parser import read_kinematic_input
from .dynamic_parser import read_dynamic_input

__all__ = [
    "read_kinematic_input",
    "read_dynamic_input",
]
