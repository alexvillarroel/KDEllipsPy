from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..inversion.dynamic.model_dynamic import DynamicInversionConfig, DynamicInversionModel


def run_dynamic_workflow(
    dynamic_root: Optional[str | Path] = None,
    *,
    compile_sources: bool = False,
) -> dict:
    model = DynamicInversionModel(dynamic_root=dynamic_root)
    config = DynamicInversionConfig(
        dynamic_root=dynamic_root,
        compile_sources=compile_sources,
    )
    return model.run_dynamic_search(config=config)

