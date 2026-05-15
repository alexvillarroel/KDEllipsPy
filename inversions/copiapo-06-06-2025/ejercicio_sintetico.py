import numpy as np
import os
import sys
from pathlib import Path

# Asegurar que podemos importar kdellipspy core
sys.path.append('/home/alex/KDEllipsPy/kdellipspy/core')

from forward_model import AxitraForwardModel
from config_parser import ConfigParser
from graphics_suite import GraphicsSuite

# 1. Inicializamos la configuración usando el nuevo método .build()
# Este método ya trae valores por defecto razonables
cfg = ConfigParser.build()

# 2. Personalizamos los parámetros básicos
cfg.observed_data.t1 = 0.0
cfg.observed_data.t2 = 128.0
cfg.observed_data.npts = 512
cfg.observed_data.delta = 0.25

cfg.source_position.event_name = 'Terremoto_Sintetico'
cfg.source_position.latitude = -33.45
cfg.source_position.longitude = -70.66
cfg.source_position.depth = 25.0
cfg.source_position.strike = 10.0
cfg.source_position.dip = 15.0
cfg.source_position.rake = 90.0

cfg.fault_plane.lx = 20000.0
cfg.fault_plane.ly = 20000.0
cfg.fault_plane.hx = 10000.0
cfg.fault_plane.hy = 10000.0
cfg.fault_plane.nx = 30
cfg.fault_plane.ny = 30

# 3. Agregamos estaciones de forma dinámica
# Estación con todas las componentes
cfg.stations.add(name='STA1', latitude=-33.0, longitude=-71.0)
# Estación sin componente Este
cfg.stations.add(name='STA2', latitude=-34.0, longitude=-70.0, use_e=False)
# Estación solo componente Vertical
cfg.stations.add(name='STA3', latitude=-33.5, longitude=-70.5, use_n=False, use_e=False)

# 4. Definimos el modelo de velocidad
cfg.velocity_model.add(thickness=0.0, vp=6000.0, vs=3500.0, rho=2700.0, qp=200.0, qs=400.0)

# 5. Inicializamos el modelo forward
# Usamos la ruta absoluta de AXITRA
axitra_path = '/home/alex/KDEllipsPy/kdellipspy/AXITRA2024'
fm = AxitraForwardModel.from_config(cfg, axitra_dir=axitra_path)

# 6. Definimos los 7 parámetros de nuestra elipse sintética
# [a1, a2, theta, np, tp, dmax, vr]
modelo_real = np.array([5.0, 10.0, 0.5, 0.2, 0.5, 3.0, 3.0])

print("Generando geometría de la elipse...")
geom = fm.build_geometry_with_ellipse_slip(modelo_real)

print("Calculando funciones de Green (Paso pesado)...")
ap = fm.build_axitra(geom)
fm.green(ap)

print("Calculando sismogramas (Paso rápido)...")
time, sismox, sismoy, sismoz = fm.conv(ap, geom)

# Apilamos para tener la forma (nsta, 3, npts) que esperabas
sismogramas = np.stack([sismox, sismoy, sismoz], axis=1)

print(f"¡Éxito! Sismogramas generados con forma: {sismogramas.shape}")

# 7. Graficamos la distribución de estaciones y los sismogramas
print("Generando gráficos...")
# Usamos el directorio actual para guardar las figuras
gs = GraphicsSuite(base_dir=Path.cwd(), show=True)

print("Mapa de estaciones...")
gs.plot_stations_map(cfg)

print("Sismogramas por componente...")
station_names = [s.name for s in cfg.stations.stations]
gs.plot_synthetic_components(time, sismogramas, station_names)

print("Gráfico de la fuente (Slip moderno)...")
gs.plot_source_ellipse(geom, title="Distribución de Slip Sintética (Estilo Moderno)")

print("Fin del ejercicio sintético.")
