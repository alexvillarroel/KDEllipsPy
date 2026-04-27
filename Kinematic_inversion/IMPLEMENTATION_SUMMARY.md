# Event-Run Folder Workflow - Implementation Summary

## ✅ Implementation Complete

The kinematic inversion pipeline has been successfully refactored to support an **event-run folder** workflow, enabling users to:

- 🎯 Run multiple earthquake inversions from different directories
- 📁 Keep inputs and outputs organized per event
- 🔄 Batch process multiple events
- 🔁 Reuse shared AXITRA installation while maintaining separate runs
- ⚙️ Maintain backward compatibility (defaults to current directory)

## 🎯 Key Changes

### main.py Refactored (Core Implementation)

**New Features:**
- ✅ CLI argument `--run <RUN_DIR>` - specify event run folder
- ✅ CLI argument `--no-plot` - skip plot generation for batch processing
- ✅ `validate_run_dir()` - validates required files (input.ctl, DATA/)
- ✅ `get_axitra_dir()` - intelligent AXITRA directory resolution
- ✅ All paths now resolve relative to run directory

**Pipeline Flow (6 Steps):**
1. Validate run directory structure
2. Load configuration from input.ctl
3. Prepare observed data
4. Initialize inversion engine
5. Execute inversion search
6. Export results and generate graphics

### Documentation Created

1. **RUN_FOLDER_GUIDE.md** (Comprehensive)
   - Overview and quick start
   - Required run folder structure
   - Data organization (flat files vs. RAW)
   - Configuration file sections
   - Example workflows
   - Troubleshooting guide
   - Performance tips

2. **QUICK_START.md** (Quick Reference)
   - Fast setup instructions
   - Common workflows
   - CLI options table

3. **examples/run_example_event.sh** (Executable Example)
   - Step-by-step shell script
   - Shows complete workflow

## 📋 Required Run Folder Structure

```
<RUN_DIR>/                          # User specifies with --run
├── input.ctl                        # ✓ REQUIRED: Config file
├── DATA/                            # ✓ REQUIRED: Observed data
│   ├── real_disp_x (or real_vel_x)  #   Option 1: Flat files
│   ├── real_disp_y (or real_vel_y)
│   ├── real_disp_z (or real_vel_z)
│   └── RAW/                         #   Option 2: SAC/MiniSEED
│       └── *.sac or *.mseed
└── [Optional] AXITRA2024/           #   Local AXITRA (else src/AXITRA2024)

# Outputs (auto-created)
output/
├── na_results.json
├── na_results.csv
└── Figures/NA_results_summary.png
```

## 🚀 Usage

### Basic Usage

```bash
# Run inversion in specified event folder
python main.py --run /path/to/my_event

# Use current directory (backward compatible)
python main.py

# Skip plot generation
python main.py --run /path/to/my_event --no-plot
```

### Batch Processing

```bash
# Process multiple events
for event_dir in events/event_*; do
    echo "Processing $event_dir..."
    python main.py --run "$event_dir" --no-plot
done

# Process from within event directory
cd my_event
../../kinematic_inversion/main.py --run .
```

## 🔍 Validation & Error Handling

**Automatic Validation:**
- ✓ Checks `<RUN_DIR>/input.ctl` exists
- ✓ Checks `<RUN_DIR>/DATA/` exists
- ✓ Validates AXITRA2024 availability

**Clear Error Messages:**
```
ERROR: Invalid run directory structure
Run directory: /path/to/my_event
Missing required files/directories:
  ✗ /path/to/my_event/input.ctl (not found or not a file)
  ✗ /path/to/my_event/DATA (not found or not a directory)

Usage: python main.py --run <RUN_DIR>
       <RUN_DIR> must contain: input.ctl and DATA/
```

## 🔧 Technical Details

### Path Resolution Strategy

| Path | Resolution |
|---|---|
| Config | `<RUN_DIR>/input.ctl` |
| Observed data | `<RUN_DIR>/DATA/` |
| AXITRA | Try `<RUN_DIR>/AXITRA2024/` → fallback `src/AXITRA2024/` |
| Output | `<RUN_DIR>/output/` (created if needed) |

### Unchanged Components

These modules already supported flexible paths - no changes needed:

- **config_parser.py** - Accepts `input_ctl_path` parameter
- **signal_utils.py** - Handles `input_ctl_path` and `data_dir` parameters
- **inversion_na.py** - Accepts `axitra_dir` parameter
- **graphics_suite.py** - Uses `base_dir` for output directories

### Backward Compatibility

✅ **Fully backward compatible:**
- If `--run` not specified, defaults to current directory
- Works with existing single-directory projects
- No breaking changes to configuration format
- All existing scripts continue to work

## 📊 Example Workflows

### Setup for New Event

```bash
# Create run folder
mkdir events/copiapo_2025
cd events/copiapo_2025
mkdir DATA

# Configure
cp ../../input.ctl .
# Edit input.ctl for your event...

# Add data
cp /raw/data/real_disp_* DATA/

# Run from project root
cd ../..
python Kinematic_inversion/main.py --run events/copiapo_2025
```

### Reuse Results Directory Structure

```bash
# Copy previous run as template
cp -r events/event_2024 events/event_2025

# Update only what changed
vim events/event_2025/input.ctl
cp new_data/real_disp_* events/event_2025/DATA/

# Run
python Kinematic_inversion/main.py --run events/event_2025
```

## ✨ Benefits

- 🎯 **Organization**: Each event self-contained with inputs/outputs
- 🔄 **Scalability**: Batch processing multiple events
- 🔁 **Reusability**: Share AXITRA, custom configs per event
- 📊 **Traceability**: Full history per event in separate folders
- ⚡ **Performance**: Process events in parallel (different --run folders)
- 🔧 **Maintenance**: Easy to archive, compare, or reproduce results

## 📝 Files Modified

### Core Implementation
- ✅ `main.py` - Complete refactor with CLI and validation

### New Documentation
- ✅ `RUN_FOLDER_GUIDE.md` - 300+ line comprehensive guide
- ✅ `QUICK_START.md` - Quick reference
- ✅ `examples/run_example_event.sh` - Executable example
- ✅ `IMPLEMENTATION_SUMMARY.md` - This file

## 🧪 Testing Checklist

- ✅ Python syntax validation
- ✅ Path resolution logic tested
- ✅ Validation error messages verified
- ✅ Backward compatibility confirmed
- ✅ AXITRA directory fallback tested
- ✅ Output directory creation verified

## 📚 Documentation Files

### For Users Starting Out
→ Read [QUICK_START.md](QUICK_START.md) first

### For Comprehensive Reference
→ See [RUN_FOLDER_GUIDE.md](RUN_FOLDER_GUIDE.md)

### For Implementation Details
→ Check [main.py](main.py) source code

## 🎓 Migration Guide

**Old workflow:**
```bash
Kinematic_inversion/
├── main.py
├── input.ctl        ← Hard-coded location
└── DATA/            ← Hard-coded location
```

**New workflow:**
```bash
# Project structure unchanged
Kinematic_inversion/main.py

# But run from event folders
events/
├── event_1/         ← python ../Kinematic_inversion/main.py --run .
│   ├── input.ctl
│   └── DATA/
└── event_2/
    ├── input.ctl
    └── DATA/
```

**Migration:** Just move your `input.ctl` and `DATA/` to event folders and run with `--run`!

---

**Status**: ✅ Ready for production use
**Branch**: feature/moment-tensor
**Last Updated**: 2025-04-24
