from pathlib import Path
import sys

# Run from Scripts/ while importing from src/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.signal_utils import write_azi_times_file


if __name__ == "__main__":
    try:
        output = write_azi_times_file(input_ctl_path=ROOT / "input.ctl")
        print(f"azi_times generado: {output}")
    except Exception as e:
        raise RuntimeError(
            f"No se pudo generar azi_times: {e}. "
            "Instala obspy en este entorno para calcular tiempos P/S."
        ) from e
