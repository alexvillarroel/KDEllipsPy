import numpy as np
import shutil
from pathlib import Path


def resample_data(base_path, npts_orig, npts_new):
    base_path = Path(base_path)
    ratio = int(npts_orig / npts_new)
    files = [
        "real_vel_x",
        "real_vel_y",
        "real_vel_z",
        "real_disp_x",
        "real_disp_y",
        "real_disp_z",
    ]

    recovery_route = base_path / "original_data"
    recovery_route.mkdir(parents=True, exist_ok=True)

    # 1. VALIDACIÓN INICIAL: ¿Ya existen respaldos de una ejecución anterior?
    # Usamos any() sobre iterdir() que es el equivalente limpio en pathlib a os.listdir
    ya_respaldado = any(recovery_route.iterdir())

    if ya_respaldado:
        print("⚠️ Atención: La carpeta 'original_data' ya contiene archivos.")
        print("Se omitirá la copia de respaldo para no borrar tus archivos originales.")

    # 2. Congelamos la lista de archivos para procesar los 3 correctamente
    all_items = list(base_path.iterdir())

    for item in all_items:
        if item.is_file() and any(f in item.name for f in files):
            # 3. RESPALDO SEGURO: Solo copia si la carpeta de respaldo estaba vacía
            if not ya_respaldado:
                shutil.copy2(item, recovery_route / item.name)
                print(f"✓ Respaldo creado para: {item.name}")

            # 4. PROCESAMIENTO: Esto SIEMPRE se ejecuta (así se sobreescriben si cambias los npts)
            txt = np.loadtxt(item)
            data_resampled = txt[::ratio]

            print(
                f"Procesado: {item.name} - Guardando {len(data_resampled)} puntos nuevos."
            )

            # 5. Guardar los datos nuevos sobreescribiendo el archivo de la carpeta principal
            np.savetxt(item, data_resampled)


if __name__ == "__main__":
    ruta_datos = "/home/alex/KDEllipsPy/inversions/2020-06-03_Calama/banda_0.02-0.1/Herreraetal/DATA/"
    puntos_originales = 512
    puntos_nuevos = 128

    print("Iniciando el proceso de remuestreo...")
    resample_data(ruta_datos, puntos_originales, puntos_nuevos)
    print("¡Proceso finalizado con éxito!")
