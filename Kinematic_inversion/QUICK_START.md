# Event-Run Folder Workflow - Quick Reference

## 🚀 Quick Start

```bash
# Create a run folder
mkdir my_event_run
cd my_event_run

# Set up structure
mkdir DATA

# Copy configuration for your event
cp ../input.ctl .

# Copy observed data
cp /path/to/data/real_disp_* DATA/

# Go back to main directory and run
cd ..
python main.py --run my_event_run
```

## 📋 Requirements

Each run folder (`--run` argument) **must** contain:

1. **input.ctl** - Configuration file with all inversion parameters
2. **DATA/** - Directory with observed waveforms
   - Either: `real_disp_x`, `real_disp_y`, `real_disp_z` (flat files)
   - Or: `DATA/RAW/` with SAC/MiniSEED files

## 📊 Output

Results automatically created in `<RUN_DIR>/output/`:
- `na_results.json` - Full results in JSON format
- `na_results.csv` - Results in CSV format
- `Figures/NA_results_summary.png` - Summary plot

## 🎯 Common Workflows

### Single Event Analysis
```bash
python main.py --run event_2025_01
```

### Batch Processing
```bash
for event in events/event_*; do
    python main.py --run "$event" --no-plot
done
```

### From Specific Directory
```bash
cd /my/data/folder
python /path/to/kinematic_inversion/main.py --run .
```

## ⚙️ Options

| Option | Purpose |
|---|---|
| `--run <DIR>` | Event run folder (default: current dir) |
| `--no-plot` | Skip plot generation (faster for batch) |

## 📖 For More Details

See [RUN_FOLDER_GUIDE.md](RUN_FOLDER_GUIDE.md) for comprehensive documentation including:
- Data format specifications
- Configuration file parameters
- Troubleshooting
- Advanced options
