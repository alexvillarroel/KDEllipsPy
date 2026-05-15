from __future__ import annotations

from pathlib import Path
import subprocess

from .config import DynamicLegacyRunConfig


def _default_dynamic_root() -> Path:
    repo_root = Path(__file__).resolve().parents[3]
    return repo_root / "Dynamic_inversion"


def _assert_required_paths(cfg: DynamicLegacyRunConfig) -> None:
    required = [
        cfg.dynamic_root / cfg.input_file,
        cfg.dynamic_root / cfg.run_script,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Required dynamic file not found: {path}")


def run_dynamic_inversion(
    cfg: DynamicLegacyRunConfig | None = None,
    *,
    compile_sources: bool = False,
) -> dict:
    resolved = cfg or DynamicLegacyRunConfig(dynamic_root=_default_dynamic_root())
    resolved.dynamic_root = Path(resolved.dynamic_root).resolve()
    _assert_required_paths(resolved)

    compile_stdout = ""
    compile_stderr = ""
    if compile_sources:
        compile_path = resolved.dynamic_root / resolved.compile_script
        if not compile_path.exists():
            raise FileNotFoundError(f"Compile script not found: {compile_path}")
        compile_result = subprocess.run(
            ["bash", str(compile_path)],
            cwd=resolved.dynamic_root,
            check=True,
            capture_output=True,
            text=True,
        )
        compile_stdout = compile_result.stdout
        compile_stderr = compile_result.stderr

    run_result = subprocess.run(
        ["bash", str(resolved.dynamic_root / resolved.run_script)],
        cwd=resolved.dynamic_root,
        check=True,
        capture_output=True,
        text=True,
    )

    return {
        "dynamic_root": str(resolved.dynamic_root),
        "event_dir": str((resolved.dynamic_root / resolved.event_dir).resolve()),
        "stdout": run_result.stdout,
        "stderr": run_result.stderr,
        "compile_stdout": compile_stdout,
        "compile_stderr": compile_stderr,
    }

