# Changelog - Event-Run Folder Workflow

## Version: feature/moment-tensor - 2025-04-24

### 🚀 Major Features Added

#### Event-Run Folder Support
- Users can now organize inversions in separate event directories
- Each event folder contains its own `input.ctl` and `DATA/`
- Outputs are automatically created in `<RUN_DIR>/output/`
- Support for batch processing multiple events

### 📝 Modified Files

#### main.py (Complete Refactor)
**Added:**
- `argparse` module for CLI argument parsing
- `--run <RUN_DIR>` argument to specify event run folder
- `--no-plot` argument to skip plot generation
- `validate_run_dir()` function for directory structure validation
- `get_axitra_dir()` function for intelligent AXITRA path resolution
- 6-step pipeline with validation as first step

**Changed:**
- All hardcoded paths now resolve relative to run directory
- `root` variable replaced with `run_dir` throughout
- Output directory creation moved into main flow
- Graphics initialization passes `run_dir` instead of src root
- Step counter updated from 5 to 6 steps

**Unchanged:**
- Scientific logic remains identical
- Configuration parsing approach
- Data loading and filtering logic
- Inversion algorithm
- Result export format

### 📚 New Documentation Files

#### QUICK_START.md
- 🚀 Quick reference guide
- Essential usage patterns
- CLI options table
- Common workflows

#### RUN_FOLDER_GUIDE.md
- 📖 Comprehensive user guide (300+ lines)
- Quick start section
- Required folder structure with diagrams
- Data organization (flat files vs RAW/SAC)
- Configuration file sections
- Example workflows
- Troubleshooting guide
- Performance tips
- Migration guide from old workflow

#### IMPLEMENTATION_SUMMARY.md
- 🔧 Technical implementation details
- Architecture overview
- Validation logic
- Path resolution strategy
- Benefits and use cases
- Testing checklist
- File changes summary

#### examples/run_example_event.sh
- 🎬 Executable example script
- Step-by-step workflow
- Demonstrates complete pipeline

### 🔄 Backward Compatibility

**Fully maintained:**
- ✅ Old usage `python main.py` still works (defaults to current dir)
- ✅ Existing single-directory projects need no changes
- ✅ Configuration file format unchanged
- ✅ Data format unchanged
- ✅ Result output format unchanged
- ✅ All existing scripts continue to work

### ✨ Key Improvements

1. **Organization**
   - Each event self-contained in dedicated folder
   - Clear separation of inputs/outputs
   - Easy to archive or compare events

2. **Scalability**
   - Batch process multiple events
   - Process events in parallel
   - Reuse shared AXITRA installation

3. **Usability**
   - Clear error messages for missing files
   - Helpful validation messages
   - Simple CLI interface
   - Sensible defaults

4. **Maintainability**
   - Minimal code changes
   - No disruption to existing code
   - Clean, readable implementation
   - Well-documented

### 🧪 Testing

**Validated:**
- ✅ Python syntax (AST parsing)
- ✅ Path resolution logic
- ✅ Validation error messages
- ✅ AXITRA directory fallback
- ✅ Output directory creation
- ✅ Backward compatibility

### 📋 Requirements Met

From the original prompt:

✅ **CLI / entrypoint changes**
- Users can run: `python main.py --run <RUN_DIR>`
- Supports relative and absolute paths
- Defaults to current directory

✅ **Validation**
- Validates `<RUN_DIR>/input.ctl` exists
- Validates `<RUN_DIR>/DATA` exists
- Clear error messages
- Usage hint on failure

✅ **Config reading**
- Canonical source is `<RUN_DIR>/input.ctl`
- Already uses ConfigParser for parsing
- Handles KEY = VALUE format
- Config passed through pipeline

✅ **Data path resolution**
- Data defaults to `<RUN_DIR>/DATA`
- Relative paths resolved against RUN_DIR
- Absolute paths used as-is

✅ **Outputs**
- Creates `<RUN_DIR>/output/`
- Writes JSON and CSV results there
- Outputs never go to source directories
- Output directory auto-created

✅ **Minimal intrusion**
- Scientific logic unchanged
- Focus on path management
- Current behavior preserved with RUN_DIR="."

### 📚 Documentation Structure

```
Kinematic_inversion/
├── main.py                          [MODIFIED] Event-run implementation
├── QUICK_START.md                   [NEW] Quick reference
├── RUN_FOLDER_GUIDE.md              [NEW] Comprehensive guide
├── IMPLEMENTATION_SUMMARY.md        [NEW] Technical details
├── examples/
│   └── run_example_event.sh         [MODIFIED] Example script
├── input.ctl                        [Unchanged] Config template
└── src/                             [Unchanged] All modules
    ├── config_parser.py
    ├── signal_utils.py
    ├── inversion_na.py
    └── graphics_suite.py
```

### 🔗 Related Issues/PRs

- Branch: `feature/moment-tensor`
- Related to: Elliptical rupture moment tensor inversion

### 📋 Migration Path

**For existing projects:**
1. Create event folders under a `events/` directory
2. Move or copy `input.ctl` to each event folder
3. Move or link `DATA/` to each event folder
4. Run with `python main.py --run events/my_event`

**Zero code changes required** - just file organization!

### 🎯 Next Steps (Optional)

Potential enhancements (not implemented):
- Config validation utility
- Batch runner with progress reporting
- Results comparison tool
- Web interface for run management
- Docker containerization

---

**Author**: Implementation by GitHub Copilot
**Date**: 2025-04-24
**Status**: ✅ Ready for production
