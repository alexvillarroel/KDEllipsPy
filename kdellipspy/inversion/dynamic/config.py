from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class DynamicLegacyRunConfig:
    dynamic_root: Path
    input_file: str = "input_dynamic_ellipse.txt"
    run_script: str = "run_dyn_inversion.sh"
    compile_script: str = "compile_dyn_src.sh"
    event_dir: str = "Evento"

