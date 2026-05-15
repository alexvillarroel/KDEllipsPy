from __future__ import annotations

from pathlib import Path


def read_text_lines(path: str | Path) -> list[str]:
    file_path = Path(path).resolve()
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    return file_path.read_text(encoding="utf-8").splitlines()

