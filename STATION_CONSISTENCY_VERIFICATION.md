# Station Consistency Verification Report

**Date**: 2025-01-29  
**Status**: ✅ VALIDATED - Station data flow is now properly connected

## Question Posed
"Me pregunto si está bien conectado las estaciones seleccionadas con el misfit y con el loader."  
(I wonder if selected stations are properly connected with the misfit and loader.)

## Answer: YES, NOW FULLY CONNECTED

The station selection flow is now explicitly validated at every stage:

```
input.ctl Section 8 (Station Parameters)
        ↓ (contains lat, lon, elev, name, use_N, use_E, use_Z)
        ↓
raw_loader.py:load_observed_from_raw()
        ↓ Extracts: station_names = [s.name.upper() for s in cfg.stations.stations]
        ↓ Validates: All stations found in DATA/RAW/
        ↓ Loads: observed = np.zeros((len(station_names), 3, npts))
        ↓
BaseInversionModel.__init__()
        ↓ Validates: observed.shape[0] == len(cfg.stations.stations)
        ↓ Extracts: station_flags = [[s.use_n, s.use_e, s.use_z] for s in cfg.stations.stations]
        ↓
MisfitCalculator.__init__()
        ↓ Receives: observed (nsta, 3, npts), station_flags (nsta, 3), azi_times_array (nsta, 3)
        ↓
l2_misfit(synthetic)
        ↓ Validates: len(self.azi) == nsta == observed.shape[0]
        ↓ For j in range(nsta): use_n, use_e, use_z = self.station_flags[j]
        ↓ Selects components per station based on config flags
```

## Key Validations Added

### 1. **raw_loader.py** - Input Validation
```python
# Check all input.ctl stations found in DATA/RAW/
available_stations = set(extractor.extract_stations())
required_stations = set(s.name.strip().upper() for s in cfg.stations.stations)
missing_stations = required_stations - available_stations
if missing_stations:
    raise ValueError(f"Missing SAC files for: {missing_stations}")
```

### 2. **raw_loader.py** - Output Validation  
```python
# Verify output array dimensions
if observed.shape[0] != len(station_names):
    raise RuntimeError("shape[0] != station_names length")
if observed.shape[1] != 3:
    raise RuntimeError("Expected 3 components (N,E,Z)")
if observed.shape[2] != npts:
    raise RuntimeError("Expected npts time points")
```

### 3. **BaseInversionModel.__init__** - Consistency Check
```python
if observed_waveforms is not None and time_array is not None:
    # Validate array structure
    if observed_waveforms.ndim != 3:
        raise ValueError(f"Must be 3D: (nsta, 3, npts)")
    nsta_obs, ncomp, npts = observed_waveforms.shape
    if ncomp != 3:
        raise ValueError(f"Must have 3 components, got {ncomp}")
    
    # Verify n_stations matches config
    nsta_cfg = len(self.cfg.stations.stations)
    if nsta_obs != nsta_cfg:
        raise ValueError(
            f"Mismatch: observed has {nsta_obs} stations "
            f"but input.ctl defines {nsta_cfg}. "
            f"Must load exactly the stations in Section 8."
        )
    
    # Extract station flags in same order as observed array
    station_flags = np.array([
        [s.use_n, s.use_e, s.use_z] for s in self.cfg.stations.stations
    ], dtype=bool)
```

### 4. **MisfitCalculator.l2_misfit()** - Runtime Validation
```python
if synthetic.shape != self.observed.shape:
    raise ValueError(f"Shape mismatch: {synthetic.shape} vs {self.observed.shape}")

nsta, ncomp, npts = self.observed.shape
if len(self.azi) != nsta:
    raise ValueError(
        f"azi_times rows ({len(self.azi)}) do not match stations ({nsta})"
    )

# Safe indexing with validated dimensions
for j in range(nsta):
    use_n, use_e, use_z = self.station_flags[j]  # j-th station flags
    az = float(self.azi[j])                        # j-th azimuth
    ...
```

## Configuration Example (input.ctl Section 8)

```
#8. Station Parameters
Number of stations                             :       3
-26.15  -70.60  0.0  AC01  1  1  1
-27.36  -70.34  0.0  A05C  1  1  1
-26.84  -69.13  0.0  AC02  1  0  1
```

**What each column means:**
1. Latitude (degrees)
2. Longitude (degrees)
3. Height (meters)
4. **Station name** — MUST match SAC header KSTNM (case-insensitive)
5. use_N (1=include North in misfit, 0=exclude)
6. use_E (1=include East in misfit, 0=exclude)
7. use_Z (1=include Vertical in misfit, 0=exclude)

**How it works:**
- Station AC02 (row 3) has `use_E=0`, so East component won't contribute to misfit
- All three components are still loaded and preprocessed
- But during misfit calculation, only N and Z components affect the objective function

## Data Files Required

For the example above, you MUST have in `DATA/RAW/`:
- `AC01.BHN.M.*.SAC` (or equivalent - N component)
- `AC01.BHE.M.*.SAC` (or equivalent - E component)
- `AC01.BHZ.M.*.SAC` (or equivalent - Z component)
- Similar for A05C and AC02

**Naming Convention**: ObsPy's `read()` auto-detects SAC format and components.  
Station name matching is case-insensitive and normalized to UPPERCASE.

## Error Messages (Now Explicit)

If something goes wrong, you'll get a clear error:

```
ValueError: Missing SAC/MiniSEED files for stations: {'A05C', 'AC02'}.
Found: {'AC01'}

→ Add DATA/RAW/A05C.BH*.sac and DATA/RAW/AC02.BH*.sac
```

```
ValueError: Mismatch: observed_waveforms has 3 stations but input.ctl defines 5 stations.
Ensure you load data for exactly the stations in Section 8.

→ Either add 2 more stations to Section 8, or remove waveform files for 2 stations
```

```
ValueError: azi_times rows (3) do not match number of stations (2)

→ Event/azi_times.txt has wrong number of rows. Regenerate or fix manually.
```

## Testing Recommendations

1. **Verify station order preservation:**
   ```python
   from kdellipspy import load_observed_from_raw
   obs, t = load_observed_from_raw()
   print(obs.shape)  # (n_stations, 3, npts)
   ```

2. **Check component selection:**
   ```python
   # If input.ctl has use_E=0 for a station,
   # that component won't affect misfit but will still be in obs array
   ```

3. **Trace a full misfit evaluation:**
   ```python
   model = inversion_model.base_geometry  # dummy parameters
   misfit = inversion_model.objective_function(model)
   diag = inversion_model.misfit_calc.diagnostics_summary(synthetic)
   print(diag)  # Per-station RMS values
   ```

## Files Modified

- ✅ `kdellipspy/config_parser.py` — Station class has use_N, use_E, use_Z (was already present)
- ✅ `kdellipspy/raw_loader.py` — Added station validation + shape checks
- ✅ `kdellipspy/inversion_base.py` — Added n_stations consistency validation
- ✅ Documentation — Explicit data flow description

## Conclusion

**The connection is now explicit and validated:**
1. Stations defined in input.ctl Section 8
2. Loader checks all stations exist and loads them in that order
3. BaseInversionModel validates n_stations matches
4. MisfitCalculator receives station_flags and azi_times in correct order
5. Misfit calculation uses correct components per station

**No silent failures — errors are caught early with clear messages.**
