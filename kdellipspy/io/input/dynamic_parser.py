from __future__ import annotations

from pathlib import Path

from .shared_parser import read_text_lines


def read_dynamic_input(path: str | Path) -> dict:
    file_path = Path(path).resolve()
    lines = read_text_lines(file_path)
    content = [
        line.strip()
        for line in lines
        if line.strip() and not line.strip().startswith(("-", "="))
    ]
    return {
        "path": str(file_path),
        "line_count": len(lines),
        "content_lines": content,
    }

