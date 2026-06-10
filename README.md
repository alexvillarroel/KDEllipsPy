# KDEllipsPy

**Kinematic Ellipse Inversion Model for Earthquakes**

KDEllipsPy is a robust Python framework for seismic source inversion using elliptical slip patches. It combines kinematic inversion (rupture kinematics and slip distribution) with the Neighbourhood Algorithm (NA) or Bayesian MCMC methods.

---

## 🎯 Features

- **Ellipse Patch Model**: Flexible representation of slip distribution using elliptical patches.
- **Multi-Parameterization**: Invert for rupture velocity, maximum slip, ellipse geometry, and position.
- **Integration with Axitra**: Rapid Green's functions via the auto-compiled Axitra engine.
- **Signal Processing**: Native SAC/MiniSEED support with robust filtering and integration pipelines.
- **Zero-Configuration**: Automatic resolution of Fortran binaries; works out-of-the-box in local environments and Google Colab.

---

## 🚀 Installation

The recommended way to install KDEllipsPy is using a Conda environment to manage scientific dependencies and compilers.

### 1. Create Conda Environment
```bash
conda env create -f environment.yml
conda activate kdellipspy
```

### 2. Install the Package
This will install the Python dependencies and automatically compile the Axitra Fortran engine.
```bash
pip install -e .
```

### Google Colab (Quick Start)
You can run a synthetic simulation directly in your browser without any local installation:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/alexvillarroel/KDEllipsPy/blob/main/Colab_Example_Synthetic.ipynb)

---

## 📖 Usage

### CLI Interface
You can run a full inversion pipeline using the command-line interface:
```bash
kdellipspy path/to/run_folder
```
*Note: The run folder should contain an `input.ctl` file and a `DATA/` subdirectory.*

### Python API
```python
import numpy as np
from kdellipspy import ConfigParser, AxitraForwardModel

# Load configuration and initialize model
cfg = ConfigParser("input.ctl")
fm = AxitraForwardModel("input.ctl")

# Setup synthetic model parameters (a1, a2, theta, np, tp, dmax, vr)
model_params = np.array([10.0, 7.0, 0.0, 0.5, 0.5, 3.0, 2.8])

# Build geometry and run simulation
geom = fm.build_geometry_with_ellipse_slip(model_params)
ap = fm.build_axitra(geom)
fm.green(ap)
time, sx, sy, sz = fm.conv(ap, geom)

# Visualize results directly from objects
fm.plot(np.array([sx, sy, sz]).transpose(1, 0, 2), time=time)
geom.plot(title="Synthetic Slip Distribution")
```

---

## 🎓 References

1. **Ellipse Inversion**: Hernandez, B., et al. (1999). Rupture history of 1997 Umbria-Marche earthquakes. *J. Geophys. Res.*
2. **Forward Modeling**: Sánchez-Sesma, F. J., & Campillo, M. (2006). *J. Acoust. Soc. Am.*
3. **NA Algorithm**: Sambridge, M. (1999). *Geophys. J. Int.*

---

## 📄 License
This project is licensed under the MIT License.
