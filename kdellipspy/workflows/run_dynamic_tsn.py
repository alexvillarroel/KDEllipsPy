"""Fase 5: build a DynamicNAInversionModel from an input.ctl (with an
optional Section 10 for the fd3d_TSN solver config) and run the NA search.

Mirrors ``kdellipspy/cli.py``'s kinematic workflow, minus plotting.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..core.config_parser import ConfigParser
from ..core.geometry import TSNFaultGridSpec
from ..core.signal_utils import build_azi_times_array, load_and_filter_observed_data
from ..inversion.base import NAResult
from ..inversion.dynamic.model_dynamic_na import DynamicNAInversionModel
from ..inversion.dynamic.tsn_bridge import TSNRunConfig


def run_dynamic_tsn_workflow(
    project_dir: str | Path,
    *,
    input_ctl: Optional[str | Path] = None,
    data_dir: Optional[str | Path] = None,
    model_name: str = "iasp91",
) -> NAResult:
    project_dir = Path(project_dir)
    input_ctl = Path(input_ctl) if input_ctl else project_dir / "input.ctl"
    data_dir = Path(data_dir) if data_dir else project_dir / "DATA"

    cfg = ConfigParser(filepath=str(input_ctl))
    if cfg.dynamic_solver is None:
        raise ValueError(
            f"{input_ctl} has no Section 10 (TSN dynamic solver config) — "
            "required for run_dynamic_tsn_workflow."
        )

    observed, time_array = load_and_filter_observed_data(
        input_ctl_path=str(input_ctl),
        data_dir=str(data_dir),
        freq1=cfg.ellipse.freq1,
        freq2=cfg.ellipse.freq2,
    )
    azi_times_array = build_azi_times_array(config=cfg, model_name=model_name)

    ds = cfg.dynamic_solver
    grid = TSNFaultGridSpec(dh=ds.dh, dip_deg=cfg.source_position.dip, nztT=ds.nztT, nabc=ds.nabc)
    run_cfg = TSNRunConfig(work_dir=ds.work_dir, nxtT=ds.nxtT, nztT=ds.nztT, dt_s=ds.dt_s, binary=ds.binary)

    model = DynamicNAInversionModel(
        config=cfg,
        observed_waveforms=observed,
        time_array=time_array,
        azi_times_array=azi_times_array,
        tsn_run_cfg=run_cfg,
        tsn_grid=grid,
    )
    try:
        return model.run_na_search()
    finally:
        model.clean()
