import os
from pathlib import Path
import numpy as np
from obspy import read, UTCDateTime
from kdellipspy.signal_utils import integrate_waveforms

def process_sac_files(input_dir, output_subdir="VELOCITY"):
    input_path = Path(input_dir)
    output_path = input_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    sac_files = list(input_path.glob("*.SAC")) + list(input_path.glob("*.sac"))

    if not sac_files:
        print(f"No se encontraron archivos SAC en {input_dir}")
        return

    # Tiempo del evento para calcular el pre-evento (basado en kdeformat.py)
    event_time = UTCDateTime("2025-06-06T17:15:06.000")

    print(f"Procesando {len(sac_files)} archivos con integración consolidada (Matlab logic)...")

    for sac_file in sac_files:
        try:
            st = read(str(sac_file))
            tr = st[0]

            if tr.stats.starttime < event_time:
                n_baseline = int((event_time - 1.0 - tr.stats.starttime) / tr.stats.delta)
                n_baseline = max(100, min(n_baseline, len(tr.data) // 2))
            else:
                n_baseline = 500

            # Usamos la función principal ahora consolidada
            integrated_data = integrate_waveforms(
                tr.data, 
                delta=tr.stats.delta,
                baseline_samples=n_baseline,
                steps=1
            )

            tr.data = integrated_data

            new_filename = sac_file.stem + "_to_vel.SAC"
            save_path = output_path / new_filename
            tr.write(str(save_path), format="SAC")
            print(f"Guardado (baseline={n_baseline} pts): {save_path.name}")

        except Exception as e:
            print(f"Error procesando {sac_file.name}: {e}")

if __name__ == "__main__":
    raw_dir = "inversions/copiapo-06-06-2025/RAW"
    process_sac_files(raw_dir)

