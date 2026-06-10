# %% Este script sirve para convertir los .sac a desplazamiento y guardarlos en un formato
#  donde todos se encuentran en un solo archivo de una columna llamado real_disp_x,
#  real_disp_y y real_disp_z.
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import UTCDateTime, read

dir_files = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/DATA_ISOLA")
stream = read(str(dir_files / "*.sac"))

# settings input.ctl
origin_time = UTCDateTime("2025-06-06T17:15:06.000Z")
t1 = origin_time
t2 = origin_time + 204.8
dt = 0.20
f = 1.0 / dt
f1 = 0.06
f2 = 0.09

for tr in stream:
    # 1. Quitar tendencia y media antes de integrar para evitar drifts
    tr.detrend("linear")
    tr.detrend("demean")

    # 2. Integrar de velocidad a desplazamiento
    tr.integrate(method="cumtrapz")

    # 3. Aplicar filtro pasabanda en desplazamiento
    tr.filter("bandpass", freqmin=f1, freqmax=f2, corners=4, zerophase=True)

    # 4. Trim y resample
    tr.trim(starttime=t1, endtime=t2, fill_value=0.0, pad=True)
    tr.resample(f)

    if tr.stats.channel.endswith("Z"):
        print(
            f"{tr.stats.sac.stla:.4f} {tr.stats.sac.stlo:.4f} {0:.4f} {tr.stats.station}     1 1 1"
        )

# %% crear arreglos para guardar los datos concatenados
n_traces = len(stream.select(channel="*Z"))
n_pts = len(stream[0].data)
print(f"Número de trazas: {n_traces}")
print(f"Cada sismograma tiene {n_pts} muestras")

real_disp_x = np.zeros((n_pts * n_traces, 1))
real_disp_y = np.zeros((n_pts * n_traces, 1))
real_disp_z = np.zeros((n_pts * n_traces, 1))

for i, tr in enumerate(stream.select(channel="*N")):
    real_disp_x[i * n_pts : (i + 1) * n_pts] = tr.data[:, np.newaxis]
for i, tr in enumerate(stream.select(channel="*E")):
    real_disp_y[i * n_pts : (i + 1) * n_pts] = tr.data[:, np.newaxis]
for i, tr in enumerate(stream.select(channel="*Z")):
    real_disp_z[i * n_pts : (i + 1) * n_pts] = tr.data[:, np.newaxis]

# %% Guardar resultados
np.savetxt(dir_files / "real_disp_x", real_disp_x)
np.savetxt(dir_files / "real_disp_y", real_disp_y)
np.savetxt(dir_files / "real_disp_z", real_disp_z)

plt.figure(figsize=(12, 6))
plt.subplot(3, 1, 1)
plt.plot(real_disp_x, label="Desplazamiento X")
plt.legend()
plt.subplot(3, 1, 2)
plt.plot(real_disp_y, label="Desplazamiento Y")
plt.legend()
plt.subplot(3, 1, 3)
plt.plot(real_disp_z, label="Desplazamiento Z")
plt.legend()
plt.tight_layout()
plt.show()
print("Archivos real_disp_x, y, z generados exitosamente.")
