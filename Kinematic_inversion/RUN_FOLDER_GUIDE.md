# Event-Run Folder Workflow Guide

## Overview

The kinematic inversion pipeline now supports an **event-run folder** workflow, allowing you to run multiple earthquake inversions while keeping inputs and outputs organized in dedicated directories.

## Quick Start

### Basic Usage

```bash
# Run inversion in a specific folder
python main.py --run /path/to/my_event

# Run inversion in current directory (default)
python main.py

# Skip plot generation
python main.py --run /path/to/my_event --no-plot
```

## Required Run Folder Structure

Each run folder (`<RUN_DIR>`) must contain:

```
<RUN_DIR>/
├── input.ctl                    # ✓ REQUIRED: Configuration file
├── DATA/                        # ✓ REQUIRED: Observed data directory
│   ├── real_disp_x              #   Option 1: Flat files (if units=1)
│   ├── real_disp_y
│   ├── real_disp_z
│   └── RAW/                     #   Option 2: SAC/MiniSEED files
│       ├── *.sac or *.mseed
│       └── [inventory files]    #   (optional, for response removal)
└── [optional] AXITRA2024/       # Optional: Local copy of AXITRA2024
```

### Files

| File/Folder | Required? | Purpose |
|---|---|---|
| `input.ctl` | YES | Control file with configuration parameters |
| `DATA/` | YES | Observed data (flat files or RAW/) |
| `AXITRA2024/` | NO | AXITRA executables; defaults to `src/AXITRA2024` if not present |

### Outputs

The pipeline creates an `output/` subdirectory:

```
<RUN_DIR>/
└── output/
    ├── na_results.json          # Results in JSON format
    ├── na_results.csv           # Results in CSV format
    └── Figures/
        └── NA_results_summary.png   # Summary plot
```

## Configuration File (input.ctl)

The `input.ctl` file contains all inversion parameters. Each run folder should have its own, tuned for that specific earthquake event.

### Basic Sections:

1. **Observed Data Parameters**: time window, sampling, units
2. **Source Position & Focal Mechanism**: latitude, longitude, depth, strike, dip, rake
3. **Fault Plane Parameters**: dimensions, hypocenter location, subfault grid
4. **Ellipse Parameters**: number of ellipses, frequency band, time shift
5. **Parameters to Invert**: which fault parameters to vary and their bounds
6. **Inversion Process**: algorithm, iterations, sampling strategy
7. **Moment Tensor** (optional): if inverting for full moment tensor
8. **Station Parameters**: seismic station coordinates

See [input.ctl](input.ctl) for a complete example.

## Data Organization

### Option 1: Flat Files (Recommended for Reproducibility)

Place flat ASCII files in `<RUN_DIR>/DATA/`:

```
real_disp_x    # For units=1 (displacement)
real_disp_y
real_disp_z
```

Or:

```
real_vel_x     # For units=2 (velocity)
real_vel_y
real_vel_z
```

**Format**: 2 columns (time, value) with `npts * nstations` rows.

### Option 2: SAC/MiniSEED Files (Raw Format)

Place waveform files in `<RUN_DIR>/DATA/RAW/`:

```
DATA/
└── RAW/
    ├── 20200603T073534-A05C-HNE.sac
    ├── 20200603T073534-A05C-HNN.sac
    ├── 20200603T073534-A05C-HNZ.sac
    ├── ...
    └── [optional_inventory.xml]
```

**Notes:**
- Files are automatically bandpass-filtered to the frequency range in `input.ctl`
- Station name is extracted from file metadata
- If inventory files present, instrument response is removed
- Requires: `obspy` for SAC/MiniSEED reading

## Example Workflow

### 1. Create a Run Folder

```bash
mkdir my_event_2025
cd my_event_2025
```

### 2. Copy Configuration and Data

```bash
# Copy and customize config for your event
cp ../input.ctl .
vim input.ctl        # Edit for your event

# Copy or link observed data
mkdir DATA
cp /path/to/observed_data/real_disp_* DATA/
```

### 3. Run Inversion

```bash
cd ..
python main.py --run my_event_2025
```

### 4. Check Results

```bash
ls my_event_2025/output/
cat my_event_2025/output/na_results.json | jq .best_model
```

## Advanced Options

### Using AXITRA from Run Directory

If you have a local copy of AXITRA2024 in your run folder:

```
my_event_2025/
├── input.ctl
├── DATA/
│   └── ...
└── AXITRA2024/          # Preferred over src/AXITRA2024
    ├── axitra
    └── ...
```

The pipeline will use it automatically.

### Disabling Plots

To skip plot generation (useful for batch processing):

```bash
python main.py --run my_event_2025 --no-plot
```

## Troubleshooting

### Error: "Invalid run directory structure"

**Check:**
- ✓ `<RUN_DIR>/input.ctl` exists and is readable
- ✓ `<RUN_DIR>/DATA/` exists and is a directory
- ✓ Data files are present in `DATA/` (either flat files or `RAW/` subdirectory)

### Error: "AXITRA2024 directory not found"

**Solutions:**
- Option 1: Copy `src/AXITRA2024` to your run folder
- Option 2: Ensure `src/AXITRA2024` exists in the project root
- Option 3: Use absolute paths if running from different working directory

### Error: "Missing flat observed files"

**Check:**
- `input.ctl` specifies correct units (1=displacement, 2=velocity)
- File names match expected format: `real_disp_x/y/z` or `real_vel_x/y/z`
- Files contain correct number of rows: `npts * nstations`

### Error: "No se pudo generar azi_times"

**Solution:** Install obspy for P/S time calculations:

```bash
pip install obspy
```

## Performance Tips

1. **Pre-process data**: Generate flat files from SAC/MiniSEED once, reuse for multiple runs
2. **Parallel runs**: Run inversions on different events simultaneously (each in own folder)
3. **Resource limits**: Adjust `ss1` and `ss_other` in `input.ctl` if runs are memory-intensive

## Migration from Old Workflow

If you have projects using the old single-directory structure:

```bash
# Old structure
Kinematic_inversion/
├── main.py
├── input.ctl          ← Move to run folder
├── DATA/              ← Move to run folder
└── src/

# New structure
Kinematic_inversion/
├── main.py
├── src/
│   └── AXITRA2024/
└── my_event/          ← New run folder
    ├── input.ctl
    ├── DATA/
    └── output/        ← Auto-created
```

Then run:

```bash
python main.py --run my_event
```

## See Also

- [input.ctl](input.ctl) - Complete configuration file example
- `src/config_parser.py` - Configuration parsing logic
- `src/signal_utils.py` - Data loading and filtering
