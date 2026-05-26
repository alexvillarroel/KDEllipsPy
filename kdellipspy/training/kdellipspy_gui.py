"""
kdellipspy_gui.py
=================
Interfaz gráfica completa para el sistema de inversión cinemática KDEllipsPy.
Permite crear, cargar y editar archivos input.ctl, y visualizar datos sismológicos.

Dependencias:
  PyQt6, matplotlib, cartopy (opcional), numpy, obspy (opcional)

Uso:
  python kdellipspy_gui.py
"""

import sys
import os
import re
from pathlib import Path
from typing import Optional, List, Dict, Any

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QListWidget, QListWidgetItem, QStackedWidget, QLabel, QPushButton,
    QLineEdit, QFormLayout, QSpinBox, QDoubleSpinBox, QCheckBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea, QGridLayout,
    QFileDialog, QMessageBox, QGroupBox, QComboBox, QFrame,
    QSizePolicy, QSplitter, QTextEdit
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QFont, QIcon, QColor

# ─── Matplotlib integrado en PyQt6 ──────────────────────────────────────────
import matplotlib
matplotlib.use("QtAgg")  # Backend para PyQt6
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# ─── Valor centinela para detectar si cartopy está disponible ───────────────
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False


# ════════════════════════════════════════════════════════════════════════════
#  SECCIÓN 1: CTL IO MANAGER — Parseo y escritura del archivo input.ctl
# ════════════════════════════════════════════════════════════════════════════

class CTLIOManager:
    """
    Responsable de parsear y escribir el archivo input.ctl.

    El archivo está dividido en secciones #1. … #9.
    Cada sección contiene pares clave : valor o tablas de datos crudos.
    """

    # Plantilla de texto para un input.ctl vacío/por defecto
    DEFAULT_TEMPLATE = """\
#1. Observed Data Parameters
Time window start (t1)          :  -10.0
Time window end (t2)            :  120.0
Number of points (Npts)         :  512
Delta / Time step               :  0.25
Units (1:disp, 2:vel)           :  1

#2. Source Position & Focal Mechanism
Event Name                      :  MyEvent
Latitude                        :  -20.0
Longitude                       :  -70.0
Depth                           :  15.0
Strike                          :  180.0
Dip                             :  45.0
Rake                            :  -90.0

#3. Fault Plane Parameters
Length along strike (Lx)        :  60000.0
Length along dip (Ly)           :  40000.0
Hypocenter position strike (Hx) :  30000.0
Hypocenter position dip (Hy)    :  20000.0
Number of subfaults along strike (Nx) :  10
Number of subfaults along dip (Ny)    :  7

#4. Ellipse Parameters & Frequency Band
Number of ellipses              :  1
Initial slip                    :  0
Slip shape                      :  1
Frequency 1 (Freq1)             :  0.02
Frequency 2 (Freq2)             :  0.10
Time shift (T0)                 :  3.0

#5. Parameters to Invert
Param 1 : Semi-axis 1 (a1 km)          :  5.0   40.0  1
Param 2 : Semi-axis 2 (a2 km)          :  5.0   30.0  1
Param 3 : Rotation angle (theta)        :  0.0    1.0  1
Param 4 : Position fraction (np)        :  0.0    1.0  1
Param 5 : Position angle (tp)           :  0.0    1.0  1
Param 6 : Maximum slip (dmax m)         :  0.5   10.0  1
Param 7 : Rupture velocity (vr km/s)    :  1.5    4.0  1

#6. Inversion Process Parameters
Algorithm type                  :  0
Number of iterations            :  10
Sample size for first iteration (SS1) :  30
Sample size for other iterations      :  30
Cells to resample               :  7
Misfit time window              :  0.0
MCMC total steps                :  500
MCMC burn-in                    :  0
MCMC proposal scale             :  0.08
MCMC thinning                   :  1
MCMC chains                     :  1

#7. Moment Tensor (Full MT)
Moment Tensor Flag              :  0
MT Scaling Mode                 :  mt_factored
Mrr                             :  0.0
Mtt                             :  0.0
Mpp                             :  0.0
Mrt                             :  0.0
Mrp                             :  0.0
Mtp                             :  0.0
Exponent (iexp)                 :  18.0

#8. Station Parameters
# lat          lon         height   name   use_N  use_E  use_Z
-19.5000    -70.1000     0.000   SL01    1   1   1
-19.8000    -69.9000     0.000   SL02    1   1   1

#9. Velocity Model 1D
# thickness    vp        vs       rho      qp      qs
  10000.0    5500.0   3180.0   2700.0   200.0   100.0
  20000.0    6200.0   3580.0   2900.0   400.0   200.0
      0.0    8100.0   4680.0   3300.0   600.0   300.0
"""

    def parse(self, filepath: str) -> Dict[str, Any]:
        """
        Lee el archivo input.ctl y devuelve un dict con las 9 secciones.
        Retorna None si el archivo no existe o está vacío.
        """
        path = Path(filepath)
        if not path.exists():
            return self._empty_data()

        content = path.read_text(encoding="utf-8")
        sections = self._split_sections(content)
        return self._parse_all(sections)

    def write(self, filepath: str, data: Dict[str, Any]) -> None:
        """
        Escribe el dict `data` como un archivo input.ctl bien formateado.
        Las tablas (estaciones, velocidades) se escriben con columnas alineadas.
        """
        lines = []

        # ── Sección #1 ──────────────────────────────────────────────────────
        lines += [
            "#1. Observed Data Parameters",
            f"Time window start (t1)          :  {data['t1']}",
            f"Time window end (t2)            :  {data['t2']}",
            f"Number of points (Npts)         :  {data['npts']}",
            f"Delta / Time step               :  {data['delta']}",
            f"Units (1:disp, 2:vel)           :  {data['units']}",
            "",
        ]

        # ── Sección #2 ──────────────────────────────────────────────────────
        lines += [
            "#2. Source Position & Focal Mechanism",
            f"Event Name                      :  {data['event_name']}",
            f"Latitude                        :  {data['latitude']}",
            f"Longitude                       :  {data['longitude']}",
            f"Depth                           :  {data['depth']}",
            f"Strike                          :  {data['strike']}",
            f"Dip                             :  {data['dip']}",
            f"Rake                            :  {data['rake']}",
            "",
        ]

        # ── Sección #3 ──────────────────────────────────────────────────────
        lines += [
            "#3. Fault Plane Parameters",
            f"Length along strike (Lx)        :  {data['lx']}",
            f"Length along dip (Ly)           :  {data['ly']}",
            f"Hypocenter position strike (Hx) :  {data['hx']}",
            f"Hypocenter position dip (Hy)    :  {data['hy']}",
            f"Number of subfaults along strike (Nx) :  {data['nx']}",
            f"Number of subfaults along dip (Ny)    :  {data['ny']}",
            "",
        ]

        # ── Sección #4 ──────────────────────────────────────────────────────
        lines += [
            "#4. Ellipse Parameters & Frequency Band",
            f"Number of ellipses              :  {data['num_ellipses']}",
            f"Initial slip                    :  {data['initial_slip']}",
            f"Slip shape                      :  {data['slip_shape']}",
            f"Frequency 1 (Freq1)             :  {data['freq1']}",
            f"Frequency 2 (Freq2)             :  {data['freq2']}",
            f"Time shift (T0)                 :  {data['t0']}",
            "",
        ]

        # ── Sección #5 — parámetros de inversión ───────────────────────────
        lines.append("#5. Parameters to Invert")
        for i, p in enumerate(data.get('inversion_params', []), 1):
            # Formato: Param N : <nombre> : min max flag
            label = f"Param {i} : {p['name']}"
            values = f"{p['min_val']}   {p['max_val']}  {p['flag']}"
            lines.append(f"{label:<50} :  {values}")
        lines.append("")

        # ── Sección #6 ──────────────────────────────────────────────────────
        lines += [
            "#6. Inversion Process Parameters",
            f"Algorithm type                  :  {data['algorithm_type']}",
            f"Number of iterations            :  {data['num_iterations']}",
            f"Sample size for first iteration (SS1) :  {data['ss1']}",
            f"Sample size for other iterations      :  {data['ss_other']}",
            f"Cells to resample               :  {data['cells_resample']}",
            f"Misfit time window              :  {data['misfit_time_window']}",
            f"MCMC total steps                :  {data['mcmc_total_steps']}",
            f"MCMC burn-in                    :  {data['mcmc_burn_in']}",
            f"MCMC proposal scale             :  {data['mcmc_proposal_scale']}",
            f"MCMC thinning                   :  {data['mcmc_thin']}",
            f"MCMC chains                     :  {data['mcmc_chains']}",
            "",
        ]

        # ── Sección #7 ──────────────────────────────────────────────────────
        lines += [
            "#7. Moment Tensor (Full MT)",
            f"Moment Tensor Flag              :  {data['mt_flag']}",
            f"MT Scaling Mode                 :  {data['mt_scaling_mode']}",
            f"Mrr                             :  {data['mrr']}",
            f"Mtt                             :  {data['mtt']}",
            f"Mpp                             :  {data['mpp']}",
            f"Mrt                             :  {data['mrt']}",
            f"Mrp                             :  {data['mrp']}",
            f"Mtp                             :  {data['mtp']}",
            f"Exponent (iexp)                 :  {data['iexp']}",
            "",
        ]

        # ── Sección #8 — tabla de estaciones ───────────────────────────────
        # Escritura alineada por columnas para compatibilidad con lectores numéricos
        lines.append("#8. Station Parameters")
        lines.append("# lat          lon         height   name   use_N  use_E  use_Z")
        for st in data.get('stations', []):
            lat  = f"{float(st['lat']):12.4f}"
            lon  = f"{float(st['lon']):12.4f}"
            h    = f"{float(st['height']):9.3f}"
            name = f"{st['name']:<8}"
            n    = "1" if st['use_n'] else "0"
            e    = "1" if st['use_e'] else "0"
            z    = "1" if st['use_z'] else "0"
            lines.append(f"{lat}  {lon}  {h}  {name}  {n}  {e}  {z}")
        lines.append("")

        # ── Sección #9 — modelo de velocidades ─────────────────────────────
        lines.append("#9. Velocity Model 1D")
        lines.append("# thickness    vp        vs       rho      qp      qs")
        for lyr in data.get('velocity_layers', []):
            th  = f"{float(lyr['thickness']):10.1f}"
            vp  = f"{float(lyr['vp']):9.1f}"
            vs  = f"{float(lyr['vs']):9.1f}"
            rho = f"{float(lyr['rho']):9.1f}"
            qp  = f"{float(lyr['qp']):8.1f}"
            qs  = f"{float(lyr['qs']):8.1f}"
            lines.append(f"{th}  {vp}  {vs}  {rho}  {qp}  {qs}")
        lines.append("")

        Path(filepath).write_text("\n".join(lines), encoding="utf-8")

    # ── Helpers internos ────────────────────────────────────────────────────

    def _split_sections(self, content: str) -> Dict[int, str]:
        """Divide el contenido en dict {número_sección: texto_sección}."""
        sections: Dict[int, str] = {}
        current = None
        chunk: List[str] = []
        for line in content.splitlines():
            m = re.match(r"^#\s*(\d+)\.", line)
            if m:
                if current is not None:
                    sections[current] = "\n".join(chunk)
                current = int(m.group(1))
                chunk = []
            elif current is not None:
                chunk.append(line)
        if current is not None:
            sections[current] = "\n".join(chunk)
        return sections

    def _kv(self, section: str) -> Dict[str, str]:
        """Extrae pares clave:valor de una sección."""
        out = {}
        for line in section.splitlines():
            if ":" in line and not line.strip().startswith("#"):
                parts = line.rsplit(":", 1)
                key = parts[0].strip()
                val = parts[1].strip()
                if key:
                    out[key] = val
        return out

    def _data_lines(self, section: str) -> List[str]:
        """Extrae líneas de datos (no comentarios, no clave:valor)."""
        out = []
        for line in section.splitlines():
            s = line.strip()
            if s and not s.startswith("#") and ":" not in s:
                out.append(s)
        return out

    def _get(self, d: Dict, key: str, default):
        """Búsqueda tolerante por prefijo normalizado."""
        wanted = key.lower().strip()
        for k, v in d.items():
            if k.lower().strip().startswith(wanted):
                return v
        return default

    def _parse_all(self, sections: Dict[int, str]) -> Dict[str, Any]:
        data = self._empty_data()
        g = self._get

        # ── #1 ──────────────────────────────────────────────────────────────
        if 1 in sections:
            kv = self._kv(sections[1])
            data.update({
                't1':    float(g(kv, 'Time window start', data['t1'])),
                't2':    float(g(kv, 'Time window end', data['t2'])),
                'npts':  int(g(kv, 'Number of points', data['npts'])),
                'delta': float(g(kv, 'Delta', data['delta'])),
                'units': int(g(kv, 'Units', data['units'])),
            })

        # ── #2 ──────────────────────────────────────────────────────────────
        if 2 in sections:
            kv = self._kv(sections[2])
            data.update({
                'event_name': str(g(kv, 'Event Name', data['event_name'])),
                'latitude':   float(g(kv, 'Latitude', data['latitude'])),
                'longitude':  float(g(kv, 'Longitude', data['longitude'])),
                'depth':      float(g(kv, 'Depth', data['depth'])),
                'strike':     float(g(kv, 'Strike', data['strike'])),
                'dip':        float(g(kv, 'Dip', data['dip'])),
                'rake':       float(g(kv, 'Rake', data['rake'])),
            })

        # ── #3 ──────────────────────────────────────────────────────────────
        if 3 in sections:
            kv = self._kv(sections[3])
            data.update({
                'lx': float(g(kv, 'Length along strike', data['lx'])),
                'ly': float(g(kv, 'Length along dip', data['ly'])),
                'hx': float(g(kv, 'Hypocenter position strike', data['hx'])),
                'hy': float(g(kv, 'Hypocenter position dip', data['hy'])),
                'nx': int(g(kv, 'Number of subfaults along strike', data['nx'])),
                'ny': int(g(kv, 'Number of subfaults along dip', data['ny'])),
            })

        # ── #4 ──────────────────────────────────────────────────────────────
        if 4 in sections:
            kv = self._kv(sections[4])
            data.update({
                'num_ellipses': int(g(kv, 'Number of ellipses', data['num_ellipses'])),
                'initial_slip': int(g(kv, 'Initial slip', data['initial_slip'])),
                'slip_shape':   int(g(kv, 'Slip shape', data['slip_shape'])),
                'freq1': float(g(kv, 'Frequency 1', data['freq1'])),
                'freq2': float(g(kv, 'Frequency 2', data['freq2'])),
                't0':    float(g(kv, 'Time shift', data['t0'])),
            })

        # ── #5 — parámetros de inversión ────────────────────────────────────
        if 5 in sections:
            params = []
            for line in sections[5].splitlines():
                m = re.match(
                    r"^\s*Param\s+\d+\s*:\s*(.*?)\s*:\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(\d+)",
                    line,
                )
                if m:
                    params.append({
                        'name':    m.group(1).strip(),
                        'min_val': float(m.group(2)),
                        'max_val': float(m.group(3)),
                        'flag':    int(m.group(4)),
                    })
            data['inversion_params'] = params

        # ── #6 ──────────────────────────────────────────────────────────────
        if 6 in sections:
            kv = self._kv(sections[6])
            data.update({
                'algorithm_type':     int(g(kv, 'Algorithm type', data['algorithm_type'])),
                'num_iterations':     int(g(kv, 'Number of iterations', data['num_iterations'])),
                'ss1':                int(g(kv, 'Sample size for first', data['ss1'])),
                'ss_other':           int(g(kv, 'Sample size for other', data['ss_other'])),
                'cells_resample':     int(g(kv, 'Cells to resample', data['cells_resample'])),
                'misfit_time_window': float(g(kv, 'Misfit time window', data['misfit_time_window'])),
                'mcmc_total_steps':   int(g(kv, 'MCMC total steps', data['mcmc_total_steps'])),
                'mcmc_burn_in':       int(g(kv, 'MCMC burn-in', data['mcmc_burn_in'])),
                'mcmc_proposal_scale':float(g(kv, 'MCMC proposal scale', data['mcmc_proposal_scale'])),
                'mcmc_thin':          int(g(kv, 'MCMC thinning', data['mcmc_thin'])),
                'mcmc_chains':        int(g(kv, 'MCMC chains', data['mcmc_chains'])),
            })

        # ── #7 ──────────────────────────────────────────────────────────────
        if 7 in sections:
            kv = self._kv(sections[7])
            data.update({
                'mt_flag':         int(g(kv, 'Moment Tensor Flag', data['mt_flag'])),
                'mt_scaling_mode': str(g(kv, 'MT Scaling Mode', data['mt_scaling_mode'])),
                'mrr':  float(g(kv, 'Mrr', data['mrr'])),
                'mtt':  float(g(kv, 'Mtt', data['mtt'])),
                'mpp':  float(g(kv, 'Mpp', data['mpp'])),
                'mrt':  float(g(kv, 'Mrt', data['mrt'])),
                'mrp':  float(g(kv, 'Mrp', data['mrp'])),
                'mtp':  float(g(kv, 'Mtp', data['mtp'])),
                'iexp': float(g(kv, 'Exponent', data['iexp'])),
            })

        # ── #8 — estaciones ─────────────────────────────────────────────────
        if 8 in sections:
            stations = []
            for line in self._data_lines(sections[8]):
                parts = line.split()
                if len(parts) >= 4:
                    stations.append({
                        'lat':    float(parts[0]),
                        'lon':    float(parts[1]),
                        'height': float(parts[2]),
                        'name':   parts[3],
                        'use_n':  int(parts[4]) == 1 if len(parts) > 4 else True,
                        'use_e':  int(parts[5]) == 1 if len(parts) > 5 else True,
                        'use_z':  int(parts[6]) == 1 if len(parts) > 6 else True,
                    })
            data['stations'] = stations

        # ── #9 — modelo de velocidades ──────────────────────────────────────
        if 9 in sections:
            layers = []
            for line in self._data_lines(sections[9]):
                parts = line.split()
                if len(parts) >= 6:
                    layers.append({
                        'thickness': float(parts[0]),
                        'vp':        float(parts[1]),
                        'vs':        float(parts[2]),
                        'rho':       float(parts[3]),
                        'qp':        float(parts[4]),
                        'qs':        float(parts[5]),
                    })
            data['velocity_layers'] = layers

        return data

    def _empty_data(self) -> Dict[str, Any]:
        """Valores por defecto cuando no hay archivo."""
        return {
            # #1
            't1': -10.0, 't2': 120.0, 'npts': 512, 'delta': 0.25, 'units': 1,
            # #2
            'event_name': 'MyEvent', 'latitude': -20.0, 'longitude': -70.0,
            'depth': 15.0, 'strike': 180.0, 'dip': 45.0, 'rake': -90.0,
            # #3
            'lx': 60000.0, 'ly': 40000.0, 'hx': 30000.0, 'hy': 20000.0,
            'nx': 10, 'ny': 7,
            # #4
            'num_ellipses': 1, 'initial_slip': 0, 'slip_shape': 1,
            'freq1': 0.02, 'freq2': 0.10, 't0': 3.0,
            # #5
            'inversion_params': [
                {'name': 'Semi-axis 1 (a1 km)',       'min_val': 5.0,  'max_val': 40.0, 'flag': 1},
                {'name': 'Semi-axis 2 (a2 km)',       'min_val': 5.0,  'max_val': 30.0, 'flag': 1},
                {'name': 'Rotation angle (theta)',     'min_val': 0.0,  'max_val': 1.0,  'flag': 1},
                {'name': 'Position fraction (np)',     'min_val': 0.0,  'max_val': 1.0,  'flag': 1},
                {'name': 'Position angle (tp)',        'min_val': 0.0,  'max_val': 1.0,  'flag': 1},
                {'name': 'Maximum slip (dmax m)',      'min_val': 0.5,  'max_val': 10.0, 'flag': 1},
                {'name': 'Rupture velocity (vr km/s)', 'min_val': 1.5,  'max_val': 4.0,  'flag': 1},
            ],
            # #6
            'algorithm_type': 0, 'num_iterations': 10, 'ss1': 30, 'ss_other': 30,
            'cells_resample': 7, 'misfit_time_window': 0.0,
            'mcmc_total_steps': 500, 'mcmc_burn_in': 0, 'mcmc_proposal_scale': 0.08,
            'mcmc_thin': 1, 'mcmc_chains': 1,
            # #7
            'mt_flag': 0, 'mt_scaling_mode': 'mt_factored',
            'mrr': 0.0, 'mtt': 0.0, 'mpp': 0.0,
            'mrt': 0.0, 'mrp': 0.0, 'mtp': 0.0, 'iexp': 18.0,
            # #8
            'stations': [],
            # #9
            'velocity_layers': [],
        }


# ════════════════════════════════════════════════════════════════════════════
#  SECCIÓN 2: VISTAS (QWidget hijos de QStackedWidget)
# ════════════════════════════════════════════════════════════════════════════

def _make_dbl(value=0.0, minimum=-1e9, maximum=1e9, decimals=4, step=0.1) -> QDoubleSpinBox:
    """Factory para QDoubleSpinBox con configuración común."""
    w = QDoubleSpinBox()
    w.setRange(minimum, maximum)
    w.setDecimals(decimals)
    w.setSingleStep(step)
    w.setValue(value)
    return w

def _make_int(value=0, minimum=0, maximum=99999) -> QSpinBox:
    """Factory para QSpinBox."""
    w = QSpinBox()
    w.setRange(minimum, maximum)
    w.setValue(value)
    return w

def _section_title(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setStyleSheet("font-weight: bold; font-size: 13px; color: #2c7be5; "
                      "border-bottom: 1px solid #2c7be5; padding-bottom: 4px;")
    return lbl


# ── Vista 1: Configuración General (#1, #2, #3, #4) ─────────────────────────

class GeneralConfigView(QWidget):
    """Agrupa las secciones #1 (datos observados), #2 (fuente), #3 (falla) y #4 (elipse)."""

    def __init__(self):
        super().__init__()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        inner = QWidget()
        layout = QVBoxLayout(inner)
        layout.setSpacing(16)

        # ── #1 Observed Data ─────────────────────────────────────────────────
        g1 = QGroupBox("#1  Datos Observados")
        f1 = QFormLayout(g1)
        self.t1    = _make_dbl(-10.0, -9999, 9999, 3, 1.0)
        self.t2    = _make_dbl(120.0, -9999, 9999, 3, 1.0)
        self.npts  = _make_int(512, 1, 65536)
        self.delta = _make_dbl(0.25, 0.001, 100.0, 4, 0.01)
        self.units = QComboBox()
        self.units.addItems(["1 - Desplazamiento (m)", "2 - Velocidad (m/s)"])
        f1.addRow("Inicio ventana (t1) [s]:", self.t1)
        f1.addRow("Fin ventana (t2) [s]:",    self.t2)
        f1.addRow("Número de puntos (Npts):", self.npts)
        f1.addRow("Paso de tiempo (Δt) [s]:", self.delta)
        f1.addRow("Unidades:",                self.units)
        layout.addWidget(g1)

        # ── #2 Source Position ────────────────────────────────────────────────
        g2 = QGroupBox("#2  Posición de la Fuente y Mecanismo Focal")
        f2 = QFormLayout(g2)
        self.event_name = QLineEdit("MyEvent")
        self.lat    = _make_dbl(-20.0, -90.0, 90.0,    5, 0.01)
        self.lon    = _make_dbl(-70.0, -180.0, 180.0,  5, 0.01)
        self.depth  = _make_dbl(15.0,  0.0,   800.0,   3, 1.0)
        self.strike = _make_dbl(180.0, 0.0,   360.0,   2, 1.0)
        self.dip    = _make_dbl(45.0,  0.0,   90.0,    2, 1.0)
        self.rake   = _make_dbl(-90.0, -180.0, 180.0,  2, 1.0)
        f2.addRow("Nombre del Evento:",   self.event_name)
        f2.addRow("Latitud [°]:",          self.lat)
        f2.addRow("Longitud [°]:",         self.lon)
        f2.addRow("Profundidad [km]:",     self.depth)
        f2.addRow("Strike [°]:",           self.strike)
        f2.addRow("Dip [°]:",              self.dip)
        f2.addRow("Rake [°]:",             self.rake)
        layout.addWidget(g2)

        # ── #3 Fault Plane ───────────────────────────────────────────────────
        g3 = QGroupBox("#3  Parámetros del Plano de Falla")
        f3 = QFormLayout(g3)
        self.lx = _make_dbl(60000.0, 0, 1e7, 1, 1000.0)
        self.ly = _make_dbl(40000.0, 0, 1e7, 1, 1000.0)
        self.hx = _make_dbl(30000.0, 0, 1e7, 1, 1000.0)
        self.hy = _make_dbl(20000.0, 0, 1e7, 1, 1000.0)
        self.nx = _make_int(10, 1, 200)
        self.ny = _make_int(7,  1, 200)
        f3.addRow("Longitud en strike (Lx) [m]:",     self.lx)
        f3.addRow("Longitud en dip (Ly) [m]:",         self.ly)
        f3.addRow("Posición hipocentro strike (Hx):", self.hx)
        f3.addRow("Posición hipocentro dip (Hy):",     self.hy)
        f3.addRow("Nº subfallas en strike (Nx):",      self.nx)
        f3.addRow("Nº subfallas en dip (Ny):",         self.ny)
        layout.addWidget(g3)

        # ── #4 Ellipse ───────────────────────────────────────────────────────
        g4 = QGroupBox("#4  Parámetros de Elipse y Banda de Frecuencia")
        f4 = QFormLayout(g4)
        self.num_ellipses = _make_int(1, 1, 10)
        self.initial_slip = QComboBox()
        self.initial_slip.addItems(["0 - No", "1 - Sí"])
        self.slip_shape = QComboBox()
        self.slip_shape.addItems(["0 - Constante", "1 - Gaussiano", "2 - Elipse"])
        self.freq1 = _make_dbl(0.02, 0.001, 100.0, 4, 0.01)
        self.freq2 = _make_dbl(0.10, 0.001, 100.0, 4, 0.01)
        self.t0_src = _make_dbl(3.0, 0.0, 1000.0, 3, 0.5)
        f4.addRow("Nº de Elipses:",        self.num_ellipses)
        f4.addRow("Slip inicial:",          self.initial_slip)
        f4.addRow("Forma del Slip:",        self.slip_shape)
        f4.addRow("Frecuencia 1 [Hz]:",     self.freq1)
        f4.addRow("Frecuencia 2 [Hz]:",     self.freq2)
        f4.addRow("Tiempo origen (T0) [s]:", self.t0_src)
        layout.addWidget(g4)

        layout.addStretch()
        scroll.setWidget(inner)

        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)

    def load(self, d: Dict):
        self.t1.setValue(d['t1'])
        self.t2.setValue(d['t2'])
        self.npts.setValue(d['npts'])
        self.delta.setValue(d['delta'])
        self.units.setCurrentIndex(d['units'] - 1)

        self.event_name.setText(d['event_name'])
        self.lat.setValue(d['latitude'])
        self.lon.setValue(d['longitude'])
        self.depth.setValue(d['depth'])
        self.strike.setValue(d['strike'])
        self.dip.setValue(d['dip'])
        self.rake.setValue(d['rake'])

        self.lx.setValue(d['lx'])
        self.ly.setValue(d['ly'])
        self.hx.setValue(d['hx'])
        self.hy.setValue(d['hy'])
        self.nx.setValue(d['nx'])
        self.ny.setValue(d['ny'])

        self.num_ellipses.setValue(d['num_ellipses'])
        self.initial_slip.setCurrentIndex(d['initial_slip'])
        self.slip_shape.setCurrentIndex(d['slip_shape'])
        self.freq1.setValue(d['freq1'])
        self.freq2.setValue(d['freq2'])
        self.t0_src.setValue(d['t0'])

    def dump(self, d: Dict):
        d['t1']    = self.t1.value()
        d['t2']    = self.t2.value()
        d['npts']  = self.npts.value()
        d['delta'] = self.delta.value()
        d['units'] = self.units.currentIndex() + 1

        d['event_name'] = self.event_name.text()
        d['latitude']   = self.lat.value()
        d['longitude']  = self.lon.value()
        d['depth']      = self.depth.value()
        d['strike']     = self.strike.value()
        d['dip']        = self.dip.value()
        d['rake']       = self.rake.value()

        d['lx'] = self.lx.value()
        d['ly'] = self.ly.value()
        d['hx'] = self.hx.value()
        d['hy'] = self.hy.value()
        d['nx'] = self.nx.value()
        d['ny'] = self.ny.value()

        d['num_ellipses'] = self.num_ellipses.value()
        d['initial_slip'] = self.initial_slip.currentIndex()
        d['slip_shape']   = self.slip_shape.currentIndex()
        d['freq1']        = self.freq1.value()
        d['freq2']        = self.freq2.value()
        d['t0']           = self.t0_src.value()


# ── Vista 2: Inversión (#5 y #6) ─────────────────────────────────────────────

class InversionView(QWidget):
    """Parámetros de inversión: rango por parámetro (#5) y proceso (#6)."""

    def __init__(self):
        super().__init__()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        layout = QVBoxLayout(inner)

        # ── #5 Tabla de parámetros a invertir ──────────────────────────────
        g5 = QGroupBox("#5  Parámetros a Invertir")
        v5 = QVBoxLayout(g5)
        self.param_table = QTableWidget(7, 4)
        self.param_table.setHorizontalHeaderLabels(["Nombre", "Mínimo", "Máximo", "Flag (0=fijo, 1=invertir)"])
        self.param_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for col in (1, 2, 3):
            self.param_table.horizontalHeader().setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)

        add_row_btn = QPushButton("➕ Agregar Parámetro")
        del_row_btn = QPushButton("🗑 Eliminar Fila")
        add_row_btn.clicked.connect(lambda: self.param_table.insertRow(self.param_table.rowCount()))
        del_row_btn.clicked.connect(self._del_param_row)
        btn_row = QHBoxLayout()
        btn_row.addWidget(add_row_btn)
        btn_row.addWidget(del_row_btn)
        btn_row.addStretch()
        v5.addWidget(self.param_table)
        v5.addLayout(btn_row)
        layout.addWidget(g5)

        # ── #6 Proceso de inversión ─────────────────────────────────────────
        g6 = QGroupBox("#6  Proceso de Inversión")
        f6 = QFormLayout(g6)
        self.algorithm_type = QComboBox()
        self.algorithm_type.addItems(["0 - Neighbourhood Algorithm (NA)", "1 - MCMC (Metropolis-Hastings)"])
        self.algorithm_type.currentIndexChanged.connect(self._toggle_mcmc)

        self.num_iterations  = _make_int(10, 1, 9999)
        self.ss1             = _make_int(30, 1, 9999)
        self.ss_other        = _make_int(30, 1, 9999)
        self.cells_resample  = _make_int(7,  1, 9999)
        self.misfit_tw       = _make_dbl(0.0, 0, 9999, 2, 1.0)

        f6.addRow("Tipo de algoritmo:",           self.algorithm_type)
        f6.addRow("Nº de iteraciones:",            self.num_iterations)
        f6.addRow("SS1 (muestras iniciales):",     self.ss1)
        f6.addRow("SS_other (por iteración):",     self.ss_other)
        f6.addRow("Celdas resampleadas (nr):",     self.cells_resample)
        f6.addRow("Ventana misfit [s] (0=todo):", self.misfit_tw)

        # ── MCMC ─────────────────────────────────────────────────────────────
        self.mcmc_group = QGroupBox("Parámetros MCMC (solo si algoritmo = 1)")
        f_mc = QFormLayout(self.mcmc_group)
        self.mcmc_total_steps    = _make_int(500, 1, 999999)
        self.mcmc_burn_in        = _make_int(0, 0, 999999)
        self.mcmc_proposal_scale = _make_dbl(0.08, 1e-6, 10.0, 4, 0.01)
        self.mcmc_thin           = _make_int(1, 1, 9999)
        self.mcmc_chains         = _make_int(1, 1, 64)
        f_mc.addRow("Total de pasos:",     self.mcmc_total_steps)
        f_mc.addRow("Burn-in:",            self.mcmc_burn_in)
        f_mc.addRow("Escala de propuesta:", self.mcmc_proposal_scale)
        f_mc.addRow("Thinning:",           self.mcmc_thin)
        f_mc.addRow("Cadenas:",            self.mcmc_chains)
        self.mcmc_group.setEnabled(False)

        layout.addWidget(g6)
        layout.addWidget(self.mcmc_group)
        layout.addStretch()
        scroll.setWidget(inner)
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)

    def _toggle_mcmc(self, idx):
        self.mcmc_group.setEnabled(idx == 1)

    def _del_param_row(self):
        row = self.param_table.currentRow()
        if row >= 0:
            self.param_table.removeRow(row)

    def load(self, d: Dict):
        params = d.get('inversion_params', [])
        self.param_table.setRowCount(len(params))
        for r, p in enumerate(params):
            self.param_table.setItem(r, 0, QTableWidgetItem(str(p['name'])))
            self.param_table.setItem(r, 1, QTableWidgetItem(str(p['min_val'])))
            self.param_table.setItem(r, 2, QTableWidgetItem(str(p['max_val'])))
            self.param_table.setItem(r, 3, QTableWidgetItem(str(p['flag'])))

        self.algorithm_type.setCurrentIndex(d['algorithm_type'])
        self.num_iterations.setValue(d['num_iterations'])
        self.ss1.setValue(d['ss1'])
        self.ss_other.setValue(d['ss_other'])
        self.cells_resample.setValue(d['cells_resample'])
        self.misfit_tw.setValue(d['misfit_time_window'])
        self.mcmc_total_steps.setValue(d['mcmc_total_steps'])
        self.mcmc_burn_in.setValue(d['mcmc_burn_in'])
        self.mcmc_proposal_scale.setValue(d['mcmc_proposal_scale'])
        self.mcmc_thin.setValue(d['mcmc_thin'])
        self.mcmc_chains.setValue(d['mcmc_chains'])

    def dump(self, d: Dict):
        params = []
        for r in range(self.param_table.rowCount()):
            def cell(c): return self.param_table.item(r, c)
            if cell(0):
                params.append({
                    'name':    cell(0).text() if cell(0) else '',
                    'min_val': float(cell(1).text()) if cell(1) else 0.0,
                    'max_val': float(cell(2).text()) if cell(2) else 1.0,
                    'flag':    int(cell(3).text()) if cell(3) else 1,
                })
        d['inversion_params']      = params
        d['algorithm_type']        = self.algorithm_type.currentIndex()
        d['num_iterations']        = self.num_iterations.value()
        d['ss1']                   = self.ss1.value()
        d['ss_other']              = self.ss_other.value()
        d['cells_resample']        = self.cells_resample.value()
        d['misfit_time_window']    = self.misfit_tw.value()
        d['mcmc_total_steps']      = self.mcmc_total_steps.value()
        d['mcmc_burn_in']          = self.mcmc_burn_in.value()
        d['mcmc_proposal_scale']   = self.mcmc_proposal_scale.value()
        d['mcmc_thin']             = self.mcmc_thin.value()
        d['mcmc_chains']           = self.mcmc_chains.value()


# ── Vista 3: Fuente / Source (#7) ────────────────────────────────────────────

class SourceView(QWidget):
    """Sección #7: Tensor de Momento. Los campos MT se habilitan con un checkbox."""

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setSpacing(12)

        layout.addWidget(_section_title("#7  Momento Tensor (Full MT)"))

        g7 = QGroupBox("Parámetros de la Fuente")
        f7 = QFormLayout(g7)

        self.mt_flag = QCheckBox("Activar Tensor de Momento (MT)")
        self.mt_flag.stateChanged.connect(self._toggle_mt)

        self.mt_scaling_mode = QComboBox()
        self.mt_scaling_mode.addItems(["no_mt", "mt_factored", "mt_strict"])

        self.iexp = _make_dbl(18.0, 0, 30, 1, 1.0)

        f7.addRow("", self.mt_flag)
        f7.addRow("Modo de escalado MT:", self.mt_scaling_mode)
        f7.addRow("Exponente (iexp):",    self.iexp)
        layout.addWidget(g7)

        # ── Grupo de componentes MT (habilitado / deshabilitado) ──────────
        self.mt_group = QGroupBox("Componentes del Tensor de Momento")
        f_mt = QFormLayout(self.mt_group)
        self.mrr = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        self.mtt = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        self.mpp = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        self.mrt = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        self.mrp = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        self.mtp = _make_dbl(0, -1e20, 1e20, 4, 1.0)
        for name, widget in [("Mrr:", self.mrr), ("Mtt:", self.mtt),
                              ("Mpp:", self.mpp), ("Mrt:", self.mrt),
                              ("Mrp:", self.mrp), ("Mtp:", self.mtp)]:
            f_mt.addRow(name, widget)
        # Todos los campos MT deshabilitados por defecto
        self.mt_group.setEnabled(False)
        layout.addWidget(self.mt_group)
        layout.addStretch()

    def _toggle_mt(self, state):
        """Habilita/deshabilita los campos del tensor de momento."""
        enabled = state == Qt.CheckState.Checked.value
        self.mt_group.setEnabled(enabled)

    def load(self, d: Dict):
        flag = d.get('mt_flag', 0)
        self.mt_flag.setChecked(bool(flag))
        mode = d.get('mt_scaling_mode', 'mt_factored')
        idx = self.mt_scaling_mode.findText(mode)
        if idx >= 0:
            self.mt_scaling_mode.setCurrentIndex(idx)
        self.iexp.setValue(d.get('iexp', 18.0))
        self.mrr.setValue(d.get('mrr', 0.0))
        self.mtt.setValue(d.get('mtt', 0.0))
        self.mpp.setValue(d.get('mpp', 0.0))
        self.mrt.setValue(d.get('mrt', 0.0))
        self.mrp.setValue(d.get('mrp', 0.0))
        self.mtp.setValue(d.get('mtp', 0.0))

    def dump(self, d: Dict):
        d['mt_flag']         = 1 if self.mt_flag.isChecked() else 0
        d['mt_scaling_mode'] = self.mt_scaling_mode.currentText()
        d['iexp']            = self.iexp.value()
        d['mrr']             = self.mrr.value()
        d['mtt']             = self.mtt.value()
        d['mpp']             = self.mpp.value()
        d['mrt']             = self.mrt.value()
        d['mrp']             = self.mrp.value()
        d['mtp']             = self.mtp.value()


# ── Vista 4: Modelo de Velocidades (#9) ──────────────────────────────────────

class VelocityModelView(QWidget):
    """Tabla editable para el modelo de velocidades 1D. Máximo 20 capas."""

    MAX_ROWS = 20
    HEADERS = ["Espesor (m)", "Vp (m/s)", "Vs (m/s)", "Rho (kg/m³)", "Qp", "Qs"]

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.addWidget(_section_title("#9  Modelo de Velocidades 1D"))
        lbl = QLabel(f"Máximo {self.MAX_ROWS} capas. La última capa debe tener espesor = 0 (semiespacio).")
        lbl.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(lbl)

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        btn_row = QHBoxLayout()
        btn_add = QPushButton("➕ Agregar Capa")
        btn_del = QPushButton("🗑 Eliminar Capa")
        btn_add.clicked.connect(self._add_row)
        btn_del.clicked.connect(self._del_row)
        btn_row.addWidget(btn_add)
        btn_row.addWidget(btn_del)
        btn_row.addStretch()
        layout.addLayout(btn_row)

    def _add_row(self):
        if self.table.rowCount() >= self.MAX_ROWS:
            QMessageBox.warning(self, "Límite", f"Máximo {self.MAX_ROWS} capas permitidas.")
            return
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c, val in enumerate([10000.0, 5500.0, 3180.0, 2700.0, 200.0, 100.0]):
            self.table.setItem(r, c, QTableWidgetItem(str(val)))

    def _del_row(self):
        row = self.table.currentRow()
        if row >= 0:
            self.table.removeRow(row)

    def load(self, d: Dict):
        layers = d.get('velocity_layers', [])
        self.table.setRowCount(0)
        for lyr in layers[:self.MAX_ROWS]:
            r = self.table.rowCount()
            self.table.insertRow(r)
            for c, key in enumerate(['thickness', 'vp', 'vs', 'rho', 'qp', 'qs']):
                self.table.setItem(r, c, QTableWidgetItem(str(lyr.get(key, 0.0))))

    def dump(self, d: Dict):
        layers = []
        for r in range(self.table.rowCount()):
            def cell(c):
                it = self.table.item(r, c)
                return float(it.text()) if it else 0.0
            layers.append({
                'thickness': cell(0), 'vp': cell(1), 'vs': cell(2),
                'rho': cell(3), 'qp': cell(4), 'qs': cell(5),
            })
        d['velocity_layers'] = layers


# ── Vista 5: Estaciones (#8) ─────────────────────────────────────────────────

class _CheckBoxWidget(QWidget):
    """QCheckBox centrado dentro de una celda de QTableWidget."""
    def __init__(self, checked=True):
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.setContentsMargins(0, 0, 0, 0)
        self.cb = QCheckBox()
        self.cb.setChecked(checked)
        layout.addWidget(self.cb)

    def isChecked(self) -> bool:
        return self.cb.isChecked()


class StationsView(QWidget):
    """
    Tabla de estaciones (#8) + mapa con Cartopy/Matplotlib embebido.

    Columnas 0-3: Nombre, Lat, Lon, Z (texto/numérico)
    Columnas 4-6: use_N, use_E, use_Z (QCheckBox centrados)
    """

    HEADERS = ["Nombre", "Latitud", "Longitud", "Altura (m)", "Usar N", "Usar E", "Usar Z"]

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.addWidget(_section_title("#8  Estaciones Sismológicas"))

        # ── Tabla ───────────────────────────────────────────────────────────
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for c in range(1, 7):
            self.table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.table)

        btn_row = QHBoxLayout()
        btn_add  = QPushButton("➕ Agregar Estación")
        btn_del  = QPushButton("🗑 Eliminar Seleccionada")
        btn_plot = QPushButton("🗺 Plotear Estaciones")
        btn_add.clicked.connect(self._add_row)
        btn_del.clicked.connect(self._del_row)
        btn_plot.clicked.connect(self.plot_stations)
        for b in (btn_add, btn_del, btn_plot):
            btn_row.addWidget(b)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        # ── Mapa embebido (Matplotlib + Cartopy) ────────────────────────────
        # El canvas se crea vacío y se actualiza al hacer clic en "Plotear".
        self.map_canvas_container = QWidget()
        self.map_canvas_container.setMinimumHeight(320)
        self.map_layout = QVBoxLayout(self.map_canvas_container)
        self.map_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.map_canvas_container)

        self._map_canvas = None  # Se crea al primer uso

    def _add_row(self):
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(f"STA{r+1:02d}"))
        self.table.setItem(r, 1, QTableWidgetItem("-20.0"))
        self.table.setItem(r, 2, QTableWidgetItem("-70.0"))
        self.table.setItem(r, 3, QTableWidgetItem("0.0"))
        for c in (4, 5, 6):
            cb = _CheckBoxWidget(True)
            self.table.setCellWidget(r, c, cb)

    def _del_row(self):
        row = self.table.currentRow()
        if row >= 0:
            self.table.removeRow(row)

    def _get_stations(self):
        stations = []
        for r in range(self.table.rowCount()):
            try:
                name = self.table.item(r, 0).text() if self.table.item(r, 0) else f"STA{r}"
                lat  = float(self.table.item(r, 1).text()) if self.table.item(r, 1) else 0.0
                lon  = float(self.table.item(r, 2).text()) if self.table.item(r, 2) else 0.0
                stations.append((name, lat, lon))
            except ValueError:
                pass
        return stations

    def plot_stations(self):
        """
        Genera un mapa de las estaciones incrustado en la GUI.
        Usa Cartopy si está disponible, si no usa un scatter simple de Matplotlib.
        """
        stations = self._get_stations()
        if not stations:
            QMessageBox.information(self, "Sin datos", "No hay estaciones para graficar.")
            return

        # Limpia el canvas anterior si existe
        if self._map_canvas:
            self.map_layout.removeWidget(self._map_canvas)
            self._map_canvas.deleteLater()
            self._map_canvas = None

        lats = [s[1] for s in stations]
        lons = [s[2] for s in stations]
        names = [s[0] for s in stations]

        pad = max(2.0, (max(lats)-min(lats))*0.5 + 1, (max(lons)-min(lons))*0.5 + 1)
        lon0, lon1 = min(lons)-pad, max(lons)+pad
        lat0, lat1 = min(lats)-pad, max(lats)+pad

        # ── Figura de Matplotlib ─────────────────────────────────────────────
        if HAS_CARTOPY:
            # Con Cartopy: proyección PlateCarree + features geográficos
            fig = Figure(figsize=(7, 4), tight_layout=True)
            ax = fig.add_subplot(1, 1, 1,
                                 projection=ccrs.PlateCarree())
            ax.set_extent([lon0, lon1, lat0, lat1], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND,   facecolor="#f0ede0")
            ax.add_feature(cfeature.OCEAN,  facecolor="#cce5ff")
            ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=":")
            ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5)
            ax.scatter(lons, lats, s=60, c="#e05c2a", zorder=5,
                       transform=ccrs.PlateCarree(), label="Estaciones")
            for name, la, lo in stations:
                ax.text(lo + 0.08, la + 0.08, name, fontsize=7,
                        transform=ccrs.PlateCarree(), zorder=6)
            ax.set_title("Mapa de Estaciones (Cartopy)")
        else:
            # Sin Cartopy: scatter simple
            fig = Figure(figsize=(7, 4), tight_layout=True)
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(lons, lats, s=60, c="#e05c2a", zorder=3)
            for name, la, lo in stations:
                ax.annotate(name, (lo, la), textcoords="offset points",
                            xytext=(4, 4), fontsize=7)
            ax.set_xlim(lon0, lon1)
            ax.set_ylim(lat0, lat1)
            ax.set_xlabel("Longitud [°]")
            ax.set_ylabel("Latitud [°]")
            ax.set_title("Mapa de Estaciones (sin Cartopy)")
            ax.grid(True, alpha=0.3)

        # ── Integra el canvas de Matplotlib en el layout de PyQt6 ────────────
        self._map_canvas = FigureCanvas(fig)
        self._map_canvas.setMinimumHeight(300)
        self.map_layout.addWidget(self._map_canvas)
        self._map_canvas.draw()

    def load(self, d: Dict):
        self.table.setRowCount(0)
        for st in d.get('stations', []):
            r = self.table.rowCount()
            self.table.insertRow(r)
            self.table.setItem(r, 0, QTableWidgetItem(str(st.get('name', ''))))
            self.table.setItem(r, 1, QTableWidgetItem(str(st.get('lat', 0.0))))
            self.table.setItem(r, 2, QTableWidgetItem(str(st.get('lon', 0.0))))
            self.table.setItem(r, 3, QTableWidgetItem(str(st.get('height', 0.0))))
            for col_idx, key in [(4, 'use_n'), (5, 'use_e'), (6, 'use_z')]:
                cb = _CheckBoxWidget(bool(st.get(key, True)))
                self.table.setCellWidget(r, col_idx, cb)

    def dump(self, d: Dict):
        stations = []
        for r in range(self.table.rowCount()):
            def txt(c):
                it = self.table.item(r, c)
                return it.text() if it else ""
            def chk(c):
                w = self.table.cellWidget(r, c)
                return w.isChecked() if w else True
            stations.append({
                'name':   txt(0),
                'lat':    float(txt(1)) if txt(1) else 0.0,
                'lon':    float(txt(2)) if txt(2) else 0.0,
                'height': float(txt(3)) if txt(3) else 0.0,
                'use_n':  chk(4),
                'use_e':  chk(5),
                'use_z':  chk(6),
            })
        d['stations'] = stations


# ── Vista 6: Visualización de Datos ──────────────────────────────────────────

class DataVisualizationView(QWidget):
    """
    Lee sismogramas desde la carpeta DATA y los muestra como una grilla scrollable
    con 3 columnas fijas por componente (N, E, Z) y una fila por estación.

    Formato esperado: archivos .sac, .mseed, o binarios numpy (.npy)
    en subcarpetas o directamente en DATA/.
    """

    def __init__(self):
        super().__init__()
        self._work_dir: Optional[Path] = None
        self._base_waveforms: Optional[np.ndarray] = None
        self._view_waveforms: Optional[np.ndarray] = None
        self._station_names: List[str] = []
        self._station_dist_km: Dict[str, float] = {}
        self._event_lat: Optional[float] = None
        self._event_lon: Optional[float] = None
        self._dt: float = 0.25

        layout = QVBoxLayout(self)
        layout.addWidget(_section_title("📊  Visualización de Sismogramas"))

        ctrl = QHBoxLayout()
        self.load_btn = QPushButton("📂 Cargar Datos desde DATA/")
        self.load_btn.clicked.connect(self.load_data)
        self.info_lbl = QLabel("Sin datos cargados")
        self.info_lbl.setStyleSheet("color: #888;")
        ctrl.addWidget(self.load_btn)
        ctrl.addWidget(self.info_lbl)
        ctrl.addStretch()
        layout.addLayout(ctrl)

        proc = QHBoxLayout()
        proc.addWidget(QLabel("fmin [Hz]:"))
        self.freq_min = _make_dbl(0.04, 0.001, 50.0, 4, 0.01)
        self.freq_min.setFixedWidth(100)
        proc.addWidget(self.freq_min)

        proc.addWidget(QLabel("fmax [Hz]:"))
        self.freq_max = _make_dbl(0.15, 0.001, 50.0, 4, 0.01)
        self.freq_max.setFixedWidth(100)
        proc.addWidget(self.freq_max)

        self.btn_bandpass = QPushButton("Aplicar Pasabanda")
        self.btn_bandpass.clicked.connect(self.apply_bandpass)
        proc.addWidget(self.btn_bandpass)

        self.btn_integrate = QPushButton("Integrar Señal")
        self.btn_integrate.clicked.connect(self.integrate_signal)
        proc.addWidget(self.btn_integrate)

        self.btn_reset = QPushButton("Reiniciar")
        self.btn_reset.clicked.connect(self.reset_processing)
        proc.addWidget(self.btn_reset)

        proc.addStretch()
        layout.addLayout(proc)

        # ── ScrollArea que contiene una grilla por estación ────────────────
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.canvas_container = QWidget()
        self.canvas_layout = QVBoxLayout(self.canvas_container)
        self.canvas_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas_layout.setSpacing(12)
        self.scroll.setWidget(self.canvas_container)
        layout.addWidget(self.scroll)

        self._station_cards: List[QWidget] = []

    def set_work_dir(self, work_dir: Path):
        self._work_dir = work_dir

    def set_event_location(self, evlat: Optional[float], evlon: Optional[float]):
        """Define coordenadas del evento desde input.ctl para ordenar por distancia."""
        self._event_lat = float(evlat) if evlat is not None else None
        self._event_lon = float(evlon) if evlon is not None else None

    def load_data(self):
        """
        Carga sismogramas desde DATA/.
        Soporta: archivos .npy (numpy) y SAC/MiniSEED via obspy (si está instalado).
        Muestra una advertencia si no hay datos válidos.
        """
        if not self._work_dir:
            QMessageBox.warning(self, "Sin directorio", "Selecciona un directorio de trabajo primero.")
            return

        data_dir = self._work_dir / "DATA"
        raw_dir = data_dir / "RAW"
        if raw_dir.exists():
            data_dir = raw_dir
        elif not data_dir.exists():
            QMessageBox.warning(self, "Sin DATA/", "La carpeta DATA/ no existe en el directorio de trabajo.")
            return

        waveforms, station_names, dt = self._try_load(data_dir)

        if waveforms is None or len(waveforms) == 0:
            QMessageBox.information(self, "Sin datos",
                "No se encontraron datos compatibles en DATA/.\n"
                "Formatos soportados: .npy (shape: [nsta, 3, npts]), .sac, .mseed")
            return

        self._base_waveforms = np.asarray(waveforms, dtype=float).copy()
        self._view_waveforms = self._base_waveforms.copy()
        self._station_names = list(station_names)
        self._station_dist_km = dict(getattr(self, "_last_station_dist_km", {}))
        self._dt = float(dt)
        self._plot_waveforms(self._view_waveforms, self._station_names, self._dt)

    def _has_loaded_waveforms(self) -> bool:
        return self._view_waveforms is not None and len(self._station_names) > 0

    def integrate_signal(self):
        """Integra una vez en el tiempo las señales actuales (por estación/componente)."""
        if not self._has_loaded_waveforms():
            QMessageBox.information(self, "Sin datos", "Carga datos antes de integrar.")
            return

        self._view_waveforms = np.cumsum(self._view_waveforms, axis=2) * self._dt
        self._plot_waveforms(self._view_waveforms, self._station_names, self._dt)
        self.info_lbl.setText(self.info_lbl.text() + " · integración aplicada")

    def apply_bandpass(self):
        """Aplica filtro pasabanda con fmin/fmax seleccionadas en la UI."""
        if not self._has_loaded_waveforms():
            QMessageBox.information(self, "Sin datos", "Carga datos antes de filtrar.")
            return

        fmin = float(self.freq_min.value())
        fmax = float(self.freq_max.value())
        if fmin <= 0 or fmax <= 0 or fmin >= fmax:
            QMessageBox.warning(self, "Frecuencias inválidas", "Debe cumplirse: 0 < fmin < fmax")
            return

        nyquist = 0.5 / max(self._dt, 1e-12)
        if fmax >= nyquist:
            QMessageBox.warning(
                self,
                "Frecuencia inválida",
                f"fmax ({fmax:.3f} Hz) debe ser menor que Nyquist ({nyquist:.3f} Hz).",
            )
            return

        try:
            from obspy.signal.filter import bandpass
        except Exception:
            QMessageBox.warning(
                self,
                "Dependencia faltante",
                "No se pudo importar obspy.signal.filter.bandpass.",
            )
            return

        filtered = self._view_waveforms.copy()
        n_sta, n_comp, _ = filtered.shape
        df = 1.0 / self._dt
        for i in range(n_sta):
            for j in range(n_comp):
                tr = filtered[i, j]
                tr = tr - np.mean(tr)
                filtered[i, j] = bandpass(tr, fmin, fmax, df=df, corners=4, zerophase=True)

        self._view_waveforms = filtered
        self._plot_waveforms(self._view_waveforms, self._station_names, self._dt)
        self.info_lbl.setText(
            self.info_lbl.text() + f" · pasabanda {fmin:.3f}-{fmax:.3f} Hz"
        )

    def reset_processing(self):
        """Restaura las trazas originales cargadas desde disco."""
        if self._base_waveforms is None:
            QMessageBox.information(self, "Sin datos", "No hay datos para reiniciar.")
            return
        self._view_waveforms = self._base_waveforms.copy()
        self._plot_waveforms(self._view_waveforms, self._station_names, self._dt)
        self.info_lbl.setText("Procesamiento reiniciado")

    def _clear_station_cards(self):
        while self.canvas_layout.count():
            item = self.canvas_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._station_cards.clear()

    def _component_figure(self, t: np.ndarray, signal: np.ndarray, title: str, color: str):
        fig = Figure(figsize=(4.2, 1.7), tight_layout=True)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(t, signal, color=color, linewidth=0.8)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_xlabel("Tiempo [s]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.25)
        return fig

    def _station_card(self, station_name: str, row_data: np.ndarray, t: np.ndarray):
        card = QGroupBox(station_name)
        card.setObjectName("stationCard")
        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(10, 12, 10, 10)
        card_layout.setSpacing(8)

        grid = QGridLayout()
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(6)

        headers = [("N", 0, "#2c7be5"), ("E", 1, "#e8623a"), ("Z", 2, "#3cb179")]
        for label, col, color in headers:
            title = QLabel(label)
            title.setAlignment(Qt.AlignmentFlag.AlignCenter)
            title.setStyleSheet("font-weight: bold; color: #8eb4e3;")
            grid.addWidget(title, 0, col)

            fig = self._component_figure(t, row_data[col], f"{station_name} - {label}", color)
            canvas = FigureCanvas(fig)
            canvas.setMinimumHeight(170)
            canvas.setMaximumHeight(200)
            grid.addWidget(canvas, 1, col)

        card_layout.addLayout(grid)
        return card

    def _inject_overview(self, waveforms: np.ndarray, station_names: List[str], dt: float):
        summary = QGroupBox("Resumen de Datos")
        layout = QVBoxLayout(summary)
        n_sta, _, npts = waveforms.shape
        layout.addWidget(QLabel(f"Estaciones cargadas: {n_sta}"))
        layout.addWidget(QLabel(f"Componentes: N / E / Z"))
        layout.addWidget(QLabel(f"Puntos por traza: {npts}"))
        layout.addWidget(QLabel(f"Delta: {dt:.3f} s"))
        if self._station_dist_km:
            ordered = []
            for s in station_names:
                dkm = self._station_dist_km.get(s, np.nan)
                if np.isfinite(dkm):
                    ordered.append(f"{s} ({dkm:.1f} km)")
                else:
                    ordered.append(s)
            layout.addWidget(QLabel("Orden por distancia al evento:"))
            layout.addWidget(QLabel(", ".join(ordered)))
        else:
            layout.addWidget(QLabel(f"Orden: {', '.join(station_names)}"))
        self.canvas_layout.addWidget(summary)

    def _try_load(self, data_dir: Path):
        """
        Intenta cargar datos en varios formatos:
        1. data.npy / waveforms.npy: array (nsta, 3, npts)
        2. Archivos .sac / .mseed via obspy
        """
        # ── Opción 1: numpy ─────────────────────────────────────────────────
        for fname in ("data.npy", "waveforms.npy", "observed.npy"):
            fpath = data_dir / fname
            if fpath.exists():
                arr = np.load(str(fpath))
                if arr.ndim == 3 and arr.shape[1] == 3:
                    n_sta = arr.shape[0]
                    names = [f"STA{i+1:02d}" for i in range(n_sta)]
                    return arr, names, 0.25
                    
        # ── Opción 2: obspy ─────────────────────────────────────────────────
        try:
            import obspy
            from obspy.geodetics import gps2dist_azimuth
            stream = obspy.Stream()
            for ext in ("*.sac", "*.SAC", "*.mseed", "*.MiniSEED"):
                for fpath in data_dir.glob(ext):
                    try:
                        stream += obspy.read(str(fpath))
                    except Exception:
                        pass
            if len(stream) > 0:
                stations_dict: Dict[str, list] = {}
                for tr in stream:
                    sta = str(tr.stats.station).strip().upper()
                    if sta not in stations_dict:
                        stations_dict[sta] = []
                    stations_dict[sta].append(tr)
                
                # Compute station distances (km), prioritizing event coordinates from input.ctl.
                dist_km: Dict[str, float] = {}
                for sname, traces in stations_dict.items():
                    dvals = []
                    for tr in traces:
                        d_m = getattr(tr.stats, "distance", None)
                        if d_m is not None:
                            try:
                                dvals.append(float(d_m) / 1000.0)
                                continue
                            except Exception:
                                pass
                        sac = getattr(tr.stats, "sac", None)
                        if sac is None:
                            continue
                        stla = getattr(sac, "stla", None)
                        stlo = getattr(sac, "stlo", None)

                        # Prefer event coordinates provided by input.ctl.
                        evla = self._event_lat
                        evlo = self._event_lon

                        # Fallback to SAC event coordinates if input.ctl values are unavailable.
                        if evla is None or evlo is None:
                            evla = getattr(sac, "evla", None)
                            evlo = getattr(sac, "evlo", None)

                        if None in (stla, stlo, evla, evlo):
                            continue
                        try:
                            d_m_calc = gps2dist_azimuth(float(evla), float(evlo), float(stla), float(stlo))[0]
                            dvals.append(float(d_m_calc) / 1000.0)
                        except Exception:
                            pass
                    if dvals:
                        dist_km[sname] = float(np.min(dvals))

                station_names = sorted(
                    stations_dict.keys(),
                    key=lambda s: (dist_km.get(s, float("inf")), s),
                )
                npts_max = max(tr.stats.npts for sta in stations_dict for tr in stations_dict[sta])
                waveforms = np.zeros((len(station_names), 3, npts_max))
                dt = stream[0].stats.delta
                for i, sname in enumerate(station_names):
                    for tr in stations_dict[sname]:
                        comp = str(getattr(tr.stats, "component", "")).upper().strip()
                        if not comp:
                            channel = str(getattr(tr.stats, "channel", "")).upper().strip()
                            comp = channel[-1:] if channel else "Z"
                        ch = {"N": 0, "E": 1, "Z": 2}.get(comp, 2)
                        n = min(tr.stats.npts, npts_max)
                        waveforms[i, ch, :n] = tr.data[:n]
                self._last_station_dist_km = dist_km
                return waveforms, station_names, dt
        except ImportError:
            pass

        # ── Opción 3: datos sintéticos de ejemplo ──────────────────────────
        n_sta, npts, dt_ex = 3, 512, 0.25
        t = np.linspace(0, npts * dt_ex, npts)
        waveforms = np.zeros((n_sta, 3, npts))
        for i in range(n_sta):
            for c in range(3):
                f0 = 0.05 + i * 0.01 + c * 0.005
                amp = np.exp(-0.02 * t) * np.sin(2 * np.pi * f0 * t)
                waveforms[i, c] = amp
        names = [f"EJ{i+1:02d}" for i in range(n_sta)]
        return waveforms, names, dt_ex

    def _plot_waveforms(self, waveforms: np.ndarray, station_names: List[str], dt: float):
        """
        Genera una grilla scrollable con una fila por estación y tres columnas por componente.
        """
        self._clear_station_cards()

        n_sta = waveforms.shape[0]
        npts  = waveforms.shape[2]
        t     = np.arange(npts) * dt
        self._inject_overview(waveforms, station_names, dt)

        for row, sname in enumerate(station_names):
            station_card = self._station_card(sname, waveforms[row], t)
            self.canvas_layout.addWidget(station_card)
            self._station_cards.append(station_card)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.canvas_layout.addWidget(spacer)

        visible_rows = min(5, n_sta)
        self.info_lbl.setText(
            f"{n_sta} estaciones cargadas · mostrando {visible_rows} visibles antes del scroll · Δt={dt:.3f}s"
        )


# ════════════════════════════════════════════════════════════════════════════
#  SECCIÓN 3: VENTANA PRINCIPAL
# ════════════════════════════════════════════════════════════════════════════

NAV_ITEMS = [
    ("⚙️  Configuración General", "Secciones #1, #2, #3, #4"),
    ("🔄  Inversión",             "Secciones #5, #6"),
    ("💥  Fuente / Source",       "Sección #7"),
    ("🌊  Modelo de Velocidades", "Sección #9"),
    ("📡  Estaciones",            "Sección #8"),
    ("📊  Visualización Data",    "Sismogramas"),
]


class MainWindow(QMainWindow):
    """
    Ventana principal de KDEllipsPy GUI.

    Layout:
      ┌─────────────┬────────────────────────────────────┐
      │  sidebar    │         content area               │
      │  QListWidget│    QStackedWidget con 6 vistas     │
      └─────────────┴────────────────────────────────────┘
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("KDEllipsPy — Editor de Configuración Sísmica")
        self.resize(1200, 750)

        self._work_dir: Optional[Path] = None
        self._ctl_path: Optional[Path] = None
        self._io = CTLIOManager()
        self._data: Dict[str, Any] = self._io._empty_data()

        self._build_ui()
        self._apply_stylesheet()
        self._ask_work_dir()

    # ── Construcción de la UI ────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Barra lateral izquierda ──────────────────────────────────────────
        sidebar = QWidget()
        sidebar.setFixedWidth(210)
        sidebar.setObjectName("sidebar")
        side_layout = QVBoxLayout(sidebar)
        side_layout.setContentsMargins(0, 0, 0, 0)
        side_layout.setSpacing(0)

        logo = QLabel("KDEllipsPy")
        logo.setObjectName("logo")
        logo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        logo.setFixedHeight(54)
        side_layout.addWidget(logo)

        self.nav = QListWidget()
        self.nav.setObjectName("navList")
        self.nav.setSpacing(2)
        for title, subtitle in NAV_ITEMS:
            item = QListWidgetItem(title)
            item.setToolTip(subtitle)
            item.setSizeHint(QSize(200, 44))
            self.nav.addItem(item)
        self.nav.setCurrentRow(0)
        self.nav.currentRowChanged.connect(self._on_nav)
        side_layout.addWidget(self.nav)
        side_layout.addStretch()

        # ── Botones inferiores del sidebar ──────────────────────────────────
        for label, slot in [("💾  Guardar CTL", self.save_ctl),
                             ("📂  Cargar CTL",  self.load_ctl)]:
            btn = QPushButton(label)
            btn.setObjectName("sideBtn")
            btn.clicked.connect(slot)
            side_layout.addWidget(btn)

        self.dir_label = QLabel("Sin directorio")
        self.dir_label.setObjectName("dirLabel")
        self.dir_label.setWordWrap(True)
        self.dir_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        side_layout.addWidget(self.dir_label)

        # ── Área de contenido: QStackedWidget ───────────────────────────────
        self.stack = QStackedWidget()

        self.view_general  = GeneralConfigView()
        self.view_inv      = InversionView()
        self.view_source   = SourceView()
        self.view_velmodel = VelocityModelView()
        self.view_stations = StationsView()
        self.view_data     = DataVisualizationView()

        for v in (self.view_general, self.view_inv, self.view_source,
                  self.view_velmodel, self.view_stations, self.view_data):
            self.stack.addWidget(v)

        root.addWidget(sidebar)
        # Separador visual
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setObjectName("separator")
        root.addWidget(sep)
        root.addWidget(self.stack)

        # ── Barra de estado ──────────────────────────────────────────────────
        self.statusBar().showMessage("Listo")

    def _on_nav(self, row: int):
        self.stack.setCurrentIndex(row)

    # ── Directorio de trabajo ────────────────────────────────────────────────

    def _ask_work_dir(self):
        """Al inicio, pide al usuario que seleccione la carpeta de trabajo."""
        QMessageBox.information(
            self, "Bienvenido a KDEllipsPy",
            "Por favor selecciona una carpeta de trabajo.\n"
            "El programa verificará o creará input.ctl y la carpeta DATA/.",
        )
        self._select_work_dir()

    def _select_work_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Seleccionar Directorio de Trabajo")
        if not path:
            return
        self._work_dir = Path(path)
        self._setup_work_dir()

    def _setup_work_dir(self):
        """Verifica o crea DATA/ e input.ctl dentro del directorio de trabajo."""
        data_dir = self._work_dir / "DATA"
        data_dir.mkdir(parents=True, exist_ok=True)

        self._ctl_path = self._work_dir / "input.ctl"
        if not self._ctl_path.exists():
            # Crea un input.ctl por defecto desde la plantilla
            self._ctl_path.write_text(CTLIOManager.DEFAULT_TEMPLATE, encoding="utf-8")
            self.statusBar().showMessage(f"Creado nuevo input.ctl en {self._work_dir}")
        else:
            self.statusBar().showMessage(f"Directorio cargado: {self._work_dir}")

        self._data = self._io.parse(str(self._ctl_path))
        self._populate_all_views()
        self.view_data.set_work_dir(self._work_dir)
        self.view_data.set_event_location(self._data.get("latitude"), self._data.get("longitude"))

        short = str(self._work_dir)
        if len(short) > 28:
            short = "…" + short[-26:]
        self.dir_label.setText(short)

    # ── Carga / Guardado ────────────────────────────────────────────────────

    def load_ctl(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Abrir input.ctl", str(self._work_dir or ""),
            "Control Files (*.ctl);;All Files (*)"
        )
        if path:
            self._ctl_path = Path(path)
            self._work_dir = self._ctl_path.parent
            self._data = self._io.parse(path)
            self._populate_all_views()
            self.view_data.set_work_dir(self._work_dir)
            self.view_data.set_event_location(self._data.get("latitude"), self._data.get("longitude"))
            self.statusBar().showMessage(f"Cargado: {path}")

    def save_ctl(self):
        if not self._ctl_path:
            self._ctl_path = Path(
                QFileDialog.getSaveFileName(
                    self, "Guardar input.ctl", "input.ctl",
                    "Control Files (*.ctl)"
                )[0]
            )
        if not self._ctl_path or str(self._ctl_path) == ".":
            return
        # Recoge datos de todas las vistas
        self._collect_all_views()
        try:
            self._io.write(str(self._ctl_path), self._data)
            self.statusBar().showMessage(f"✔ Guardado: {self._ctl_path}")
            QMessageBox.information(self, "Guardado", f"Archivo guardado exitosamente:\n{self._ctl_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error al guardar", str(e))

    def _populate_all_views(self):
        """Distribuye `self._data` hacia todas las vistas."""
        self.view_general.load(self._data)
        self.view_inv.load(self._data)
        self.view_source.load(self._data)
        self.view_velmodel.load(self._data)
        self.view_stations.load(self._data)
        self.view_data.set_event_location(self._data.get("latitude"), self._data.get("longitude"))

    def _collect_all_views(self):
        """Recoge los valores de todas las vistas hacia `self._data`."""
        self.view_general.dump(self._data)
        self.view_inv.dump(self._data)
        self.view_source.dump(self._data)
        self.view_velmodel.dump(self._data)
        self.view_stations.dump(self._data)

    # ── Stylesheet ──────────────────────────────────────────────────────────

    def _apply_stylesheet(self):
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #1e2433;
                color: #dce3f0;
                font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
                font-size: 12px;
            }
            QWidget#sidebar {
                background-color: #151a27;
            }
            QLabel#logo {
                font-size: 17px;
                font-weight: bold;
                color: #4fa3e8;
                background-color: #0e1220;
                letter-spacing: 1px;
            }
            QListWidget#navList {
                background-color: #151a27;
                border: none;
                color: #c0cade;
                font-size: 12px;
            }
            QListWidget#navList::item {
                padding: 8px 14px;
                border-radius: 6px;
                margin: 2px 6px;
            }
            QListWidget#navList::item:selected {
                background-color: #2c7be5;
                color: #ffffff;
            }
            QListWidget#navList::item:hover:!selected {
                background-color: #222d45;
            }
            QPushButton#sideBtn {
                background-color: #2c7be5;
                color: white;
                border: none;
                padding: 9px 12px;
                font-size: 12px;
                margin: 4px 8px;
                border-radius: 5px;
            }
            QPushButton#sideBtn:hover {
                background-color: #1a60c8;
            }
            QLabel#dirLabel {
                color: #6b7a99;
                font-size: 10px;
                padding: 6px;
            }
            QFrame#separator {
                color: #2a3045;
                max-width: 1px;
            }
            QGroupBox {
                border: 1px solid #2a3450;
                border-radius: 6px;
                margin-top: 10px;
                padding-top: 6px;
                font-weight: bold;
                color: #8eb4e3;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
            }
            QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox {
                background-color: #252e45;
                border: 1px solid #2e3c5a;
                border-radius: 4px;
                padding: 4px 7px;
                color: #dce3f0;
                min-width: 100px;
            }
            QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {
                border-color: #2c7be5;
            }
            QTableWidget {
                background-color: #1e2740;
                border: 1px solid #2a3450;
                gridline-color: #2a3450;
                color: #dce3f0;
            }
            QHeaderView::section {
                background-color: #151a27;
                color: #8eb4e3;
                border: none;
                padding: 5px;
                font-weight: bold;
            }
            QPushButton {
                background-color: #252e45;
                color: #dce3f0;
                border: 1px solid #2e3c5a;
                padding: 5px 12px;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #2c7be5;
                color: white;
                border-color: #2c7be5;
            }
            QScrollBar:vertical {
                background: #1a2036;
                width: 8px;
                border-radius: 4px;
            }
            QScrollBar::handle:vertical {
                background: #2e3c5a;
                border-radius: 4px;
                min-height: 20px;
            }
            QStatusBar {
                background: #0e1220;
                color: #6b7a99;
                font-size: 11px;
            }
        """)


# ════════════════════════════════════════════════════════════════════════════
#  ENTRYPOINT
# ════════════════════════════════════════════════════════════════════════════

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("KDEllipsPy")
    app.setOrganizationName("KDEllipsPy")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
