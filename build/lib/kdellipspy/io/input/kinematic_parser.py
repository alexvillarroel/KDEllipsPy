from __future__ import annotations

from pathlib import Path

from ...core.config_parser import ConfigParser


def read_kinematic_input(input_ctl_path: str | Path) -> ConfigParser:
    return ConfigParser(str(Path(input_ctl_path).resolve()))

