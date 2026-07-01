from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import DynamicLegacyRunConfig
from .legacy_bridge import run_dynamic_inversion


@dataclass
class DynamicInversionConfig:
    dynamic_root: Optional[str | Path] = None
    compile_sources: bool = False
    input_file: str = "input_dynamic_ellipse.txt"
    run_script: str = "run_dyn_inversion.sh"
    compile_script: str = "compile_dyn_src.sh"
    event_dir: str = "Evento"


class DynamicInversionModel:
    def __init__(self, dynamic_root: Optional[str | Path] = None):
        self.dynamic_root = Path(dynamic_root).resolve() if dynamic_root else None

    def run_dynamic_search(self, config: Optional[DynamicInversionConfig] = None) -> dict:
        cfg = config or DynamicInversionConfig(dynamic_root=self.dynamic_root)
        selected_root = (
            Path(cfg.dynamic_root).resolve()
            if cfg.dynamic_root is not None
            else self.dynamic_root
        )

        legacy_cfg = None
        if selected_root is not None:
            legacy_cfg = DynamicLegacyRunConfig(
                dynamic_root=selected_root,
                input_file=cfg.input_file,
                run_script=cfg.run_script,
                compile_script=cfg.compile_script,
                event_dir=cfg.event_dir,
            )

        return run_dynamic_inversion(legacy_cfg, compile_sources=cfg.compile_sources)
