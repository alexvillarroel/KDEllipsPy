# KDEllipsPy

**Kinematic and Dynamic Ellipse Inversion Model for Earthquakes**

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A robust Python framework for seismic source inversion using elliptical slip patches. Combines kinematic inversion (rupture kinematics and slip distribution) with dynamic inversion (moment tensor) using Neighbourhood Algorithm (NA) or MCMC methods.

---

## 🎯 Features

### Kinematic Inversion
- **Ellipse Patch Model**: Flexible representation of slip distribution using elliptical patches
- **Multi-Parameterization**: Invert for rupture velocity, maximum slip, ellipse geometry, and position
- **Neighbourhood Algorithm (NA)**: Fast global optimization for efficient parameter space exploration
- **Event-Based Workflows**: Organize inversions by earthquake with dedicated run folders

### Data Processing
- **SAC/MiniSEED Support**: Native reading and preprocessing of waveform data
- **Robust Signal Processing**:
  - Zero-phase bandpass filtering with manual padding (Boore 2001)
  - Polynomial baseline correction (linear for VEL, quadratic for DISP)
  - Multi-step integration (ACC→VEL) with detrending and tapering
  - Instrument response removal (if inventory available)
- **Velocity Output**: All waveforms processed to velocity (m/s)

### Forward Modeling
- **AXITRA Integration**: Rapid Green's functions via AXITRA (Sánchez-Sesma & Campillo 2006)
- **1D Velocity Models**: Flexible stratified Earth model support

### Inversion Algorithms
- **Neighbourhood Algorithm (NA)**: Resampling algorithm for efficient global search
- **MCMC**: PyMC + ArviZ for Bayesian parameter estimation
- **Custom Misfits**: Compute waveform misfit with flexible time windows and component selection

### Visualization
- **Result Plots**: Misfit evolution, parameter distributions, slip distributions
- **Waveform Comparison**: Observed vs. synthetic waveform overlays
- **Focal Mechanism Plots**: Strike/dip/rake visualization

---

## 📋 Project Structure

```
KDEllipsPy/
├── README.md                              # This file
├── kdellipspy/                            # Main Python package
│   ├── __init__.py                        # Public API
│   ├── config_parser.py                   # input.ctl parsing
│   ├── sac_processor.py                   # SAC/MiniSEED processing pipeline
│   ├── raw_loader.py                      # DATA/RAW integration
│   ├── forward_model.py                   # AXITRA forward modeling
│   ├── geometry.py                        # Fault geometry and ellipse mapping
│   ├── inversion_base.py                  # Base inversion classes
│   ├── inversion_na.py                    # Neighbourhood Algorithm
│   ├── inversion_mcmc.py                  # MCMC inversion
│   ├── graphics_suite.py                  # Plotting and visualization
│   ├── cli.py                             # Command-line interface
│   └── AXITRA2024/                        # AXITRA binaries and wrappers
│
├── Kinematic_inversion/                   # Event-run workflow
│   ├── main.py                            # Pipeline entry point
│   ├── input.ctl                          # Example configuration
│   ├── src/                               # Source modules (duplicated from kdellipspy)
│   ├── Scripts/                           # Utility scripts
│   └── examples/                          # Example notebooks and scripts
│
├── Dynamic_inversion/                     # Legacy dynamic inversion code
├── examples/                              # Jupyter notebooks and examples
└── docs/                                  # API reference and guides
```

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/alexvillarroel/KDEllipsPy.git
cd KDEllipsPy

# Install package and dependencies
pip install -e .

# Install optional dependencies for full functionality
pip install -e ".[dev,mcmc]"
```

### Required Dependencies
- `obspy>=1.4.0` — Seismic data I/O and processing
- `numpy>=1.20` — Numerical computing
- `scipy>=1.7` — Scientific algorithms
- `pandas>=1.3` — Data handling
- `matplotlib>=3.5` — Visualization

### Optional Dependencies
- `neighpy` — Neighbourhood Algorithm (NA)
- `pymc` — MCMC sampling
- `arviz` — MCMC diagnostics visualization

### Run Your First Inversion

```bash
# Navigate to kinematic inversion
cd Kinematic_inversion

# Create event run folder
mkdir my_event_2025
cd my_event_2025

# Set up structure
mkdir DATA
cp ../input.ctl .

# Copy your SAC files to DATA/RAW/
mkdir DATA/RAW
cp /path/to/data/*.sac DATA/RAW/

# Edit input.ctl with your event parameters
# - Event name, origin date/time
# - Latitude, longitude, depth
# - Strike, dip, rake
# - Frequency band (e.g., 0.04–0.2 Hz)
# - Inversion parameters

# Run inversion
cd ..
python main.py --run my_event_2025 --data-mode raw
```

Results will be saved in `my_event_2025/output/`:
- `na_results.json` — Full inversion results
- `na_results.csv` — Parameter estimates
- `Figures/NA_results_summary.png` — Summary plot

---

## 📖 Workflows

### 1. Kinematic Inversion (Ellipse Patches)

**Input**: Observed waveforms (SAC/MiniSEED) + event parameters

**Output**: Optimal slip distribution, rupture velocity, and ellipse parameters

**Process**:
1. Load and preprocess waveforms (bandpass filtering, integration ACC→VEL)
2. Build fault geometry with elliptical patches
3. Run forward modeling to compute synthetic waveforms
4. Execute NA to optimize parameters
5. Export results and visualization

**Configuration**: `input.ctl` (Section 4–5 control ellipse parameters and inversion settings)

### 2. Data Processing from SAC Files

**Example**:
```python
from kdellipspy import load_observed_from_raw

# Load SAC files from DATA/RAW/
observed, time = load_observed_from_raw(
    input_ctl_path="input.ctl",
    data_dir="DATA",  # Looks for DATA/RAW/
)
# Returns: observed (n_stations, 3, npts), time (npts,)
```

**Processing Pipeline**:
- Detrending and baseline correction (pre-event window)
- Cosine tapering (5% maximum)
- Bandpass filtering: zero-phase with 100 s padding
- Instrument response removal (if StationXML/RESP available)
- Integration: ACC → VEL (polynomial correction, degree=1)
- Resampling and time windowing

---

## ⚙️ Configuration (`input.ctl`)

### Section 1: Observed Data Parameters
- Time window (t1, t2)
- Number of points and sampling interval
- Units: 1=displacement, 2=velocity

### Section 2: Source Position & Focal Mechanism
- **Event Name**, **Origin Date** (UTC, YYYY-MM-DD), **Origin Time** (HH:MM:SS)
- Latitude, longitude, depth
- Strike, dip, rake

### Section 3: Fault Plane Parameters
- Fault length and width
- Hypocenter position
- Subfault discretization (Nx, Ny)

### Section 4: Ellipse Parameters & Frequency Band
- Number of ellipses
- Initial slip and slip shape
- **Frequency band** (Freq1, Freq2) — critical for signal processing
- Time shift (T0)

### Section 5: Parameters to Invert
- Which parameters are free vs. fixed
- Min/max bounds for each parameter
- Inversion flags (0=fixed, 1=invert)

### Sections 6–9
- Inversion process settings (algorithm, iterations, sampling)
- Optional moment tensor
- Station coordinates
- 1D velocity model

---

## 📊 Data Organization

### Flat File Mode (Pre-processed)
```
RUN_DIR/
├── input.ctl
└── DATA/
    ├── real_vel_x
    ├── real_vel_y
    └── real_vel_z
```
Files: 2 columns (time, value) or value-only, npts × nstations rows total.

### RAW Mode (Unprocessed SAC/MiniSEED)
```
RUN_DIR/
├── input.ctl
└── DATA/
    └── RAW/
        ├── STATION1.HNE.sac
        ├── STATION1.HNN.sac
        ├── STATION1.HNZ.sac
        ├── ...
        └── [inventory.xml]  (optional, for response removal)
```

**Preprocessing is applied automatically** when using `--data-mode raw`.

---

## 🔍 Examples

### Example 1: Simple NA Inversion
```python
from pathlib import Path
from kdellipspy import ConfigParser, load_observed_flat, NAInversionModel, NAConfig

run_dir = Path("my_event_2025")
input_ctl = run_dir / "input.ctl"

# Load configuration
cfg = ConfigParser(str(input_ctl))

# Load observed waveforms
observed, time = load_observed_flat(
    input_ctl_path=input_ctl,
    data_dir=run_dir / "DATA",
)

# Initialize inversion
inversion = NAInversionModel(
    str(input_ctl),
    observed_waveforms=observed,
    time_array=time,
)

# Run NA search
config = NAConfig(
    n_samples_initial=100,
    n_samples_iteration=30,
    n_iterations=10,
)
result = inversion.run_na_search(config)

# Export results
result.export_results(run_dir / "output" / "results.json")
```

### Example 2: Load and Preprocess SAC Files
```python
from kdellipspy import load_observed_from_raw

observed, time = load_observed_from_raw(
    input_ctl_path="input.ctl",
    data_dir="DATA",
)

print(f"Waveform shape: {observed.shape}")  # (n_stations, 3, npts)
print(f"Time shape: {time.shape}")            # (npts,)
```

See `examples/` directory for Jupyter notebooks with full workflows.

---

## 📚 Documentation

- [Kinematic Inversion Guide](Kinematic_inversion/RUN_FOLDER_GUIDE.md) — Detailed workflow
- [API Reference](docs/api_reference.rst) — Full API documentation
- [Quick Start](Kinematic_inversion/QUICK_START.md) — 5-minute introduction
- [Configuration Guide](Kinematic_inversion/input.ctl) — Annotated example

---

## 🔧 Advanced Usage

### Custom Signal Processing Configuration
```python
from kdellipspy import ProcessingConfig, load_observed_from_raw

config = ProcessingConfig(
    freq_min=0.05,
    freq_max=0.25,
    filter_corners=4,
    zerophase=True,
    zeropad_duration=150.0,
    polynomial_degree_vel=1,
    polynomial_degree_disp=2,
)

observed, time = load_observed_from_raw(
    input_ctl_path="input.ctl",
    data_dir="DATA",
    config=config,
)
```

### Running Multiple Events in Parallel
```bash
# Create event directories
for date in 2025-06-06 2025-06-07; do
    mkdir event_$date/DATA/RAW
    # Copy SAC files...
done

# Run inversions in parallel
for event in event_*/; do
    python Kinematic_inversion/main.py --run "$event" --no-plot &
done
wait
```

---

## 🎓 References

1. **Ellipse Inversion**:
   - Hernandez, B., Cocco, M., Cotton, F., & Campillo, M. (1999). Rupture history of 1997 Umbria-Marche (Central Italy) earthquakes from teleseismic and regional seismograms. *J. Geophys. Res.*

2. **Signal Processing**:
   - Boore, D. M. (2001). Effect of baseline corrections on displacements and response spectra. *Bull. Seism. Soc. Am.*
   - Boore, D. M., & Bommer, J. J. (2005). Processing of peak and pseudo-absolute acceleration recordings. *Soil Dyn. Earthq. Eng.*

3. **Neighbourhood Algorithm**:
   - Sambridge, M. (1999). Geophysical inversion with a Neighbourhood Algorithm. *Geophys. J. Int.*

4. **Forward Modeling**:
   - Sánchez-Sesma, F. J., & Campillo, M. (2006). Retrieval of the Green's function from cross correlation. *J. Acoust. Soc. Am.*

---

## 📝 Citation

If you use KDEllipsPy in your research, please cite:

```bibtex
@software{kdellipspy2025,
  author = {Villarroel, A.},
  title = {KDEllipsPy: Kinematic and Dynamic Ellipse Inversion for Earthquakes},
  year = {2025},
  url = {https://github.com/alexvillarroel/KDEllipsPy}
}
```

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit changes (`git commit -am 'Add feature'`)
4. Push to branch (`git push origin feature/your-feature`)
5. Open a Pull Request

---

## 📧 Contact & Support

For issues, questions, or suggestions:
- Open an [Issue](https://github.com/alexvillarroel/KDEllipsPy/issues)
- Contact: [Alejandro Villarroel] — `alexvillarroel@udec.cl`

---

## 📄 License

This project is licensed under the MIT License — see [LICENSE](LICENSE) file for details.