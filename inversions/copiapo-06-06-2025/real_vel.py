# %% Este script sirve para convertir los .sac a un formato en donde
#  todas las velocidades se encuentran en un solo archivo de una columna llamado real_vel_x,
#  real_vel_y y real_vel_z.
# para ello, se resamplea a lo ingresado en input.ctl
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read

dir_files = Path("/home/alex/KDEllipsPy/inversions/copiapo-06-06-2025/SELECTED_DATA")
stream = read(str(dir_files) + "/*.sac")
# settings input.ctl
origin_time = UTCDateTime("2025-06-06T17:15:06.000Z")
t1 = origin_time
t2 = origin_time + 128.0
dt = 0.25
f = 1.0 / dt
f1 = 0.1
f2 = 0.2

for tr in stream:
    # fill the trace with trim
    tr.trim(starttime=t1, endtime=t2, fill_value=0.0, pad=True)
    # resample to input.ctl frequency
    tr.filter("bandpass", freqmin=f1, freqmax=f2, corners=4, zerophase=True)
    tr.resample(f)
    if tr.stats.channel.endswith("Z"):
        print(
            f"{tr.stats.sac.stla:.4f} {tr.stats.sac.stlo:.4f} {0:.4f} {tr.stats.station}     1 1 1"
        )
print(f"Cada sismograma tiene {len(tr.data)} muestras")

# %% crear una columna de largo len(tr.data) x len(stream.select(channel='*Z'))
n_traces = len(stream.select(channel="*Z"))
print(f"Número de trazas: {n_traces}")
real_vel_x = np.zeros((len(stream[0].data) * n_traces, 1))
real_vel_y = np.zeros((len(stream[0].data) * n_traces, 1))
real_vel_z = np.zeros((len(stream[0].data) * n_traces, 1))

for i, tr in enumerate(stream.select(channel="*N")):
    real_vel_x[i * len(stream[0].data) : (i + 1) * len(stream[0].data)] = tr.data[
        :, np.newaxis
    ]
for i, tr in enumerate(stream.select(channel="*E")):
    real_vel_y[i * len(stream[0].data) : (i + 1) * len(stream[0].data)] = tr.data[
        :, np.newaxis
    ]
for i, tr in enumerate(stream.select(channel="*Z")):
    real_vel_z[i * len(stream[0].data) : (i + 1) * len(stream[0].data)] = tr.data[
        :, np.newaxis
    ]
    # %%

plt.figure()
plt.plot(real_vel_x)
plt.show()
# %% Save real_vel_x, real_vel_y, real_vel_z as text files
np.savetxt(dir_files / "real_vel_x", real_vel_x)
np.savetxt(dir_files / "real_vel_y", real_vel_y)
np.savetxt(dir_files / "real_vel_z", real_vel_z)
