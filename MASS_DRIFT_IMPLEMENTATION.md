# Mass Drift Diagnostic Implementation

## Overview

Implemented a comprehensive mass drift diagnostic to validate that mass intake rate Ṁ is negligible in the Solar System regime. This diagnostic tracks mass conservation over the trajectory and provides pass/fail validation against strict thresholds.

## Implementation Details

### 1. New Function: `compute_mass_drift()` in `slab/diagnostics.py`

**Location:** Lines 884-1019 in `slab/diagnostics.py`

**Purpose:** Computes mass drift diagnostics for each body over the entire trajectory.

**Inputs:**
- `trajectory`: Dictionary from `integrate_orbit()` with keys:
  - `'t'`: time array, shape (n_steps,)
  - `'M'`: mass array, shape (n_steps, n_bodies)
- `body_names`: List of body names for labeling

**Outputs:**
Dictionary with:
- `'bodies'`: list of body names
- `'initial_mass'`: M(0) for each body
- `'final_mass'`: M(final) for each body
- `'abs_drift'`: |M(final) - M(0)| for each body
- `'rel_drift'`: |M(final) - M(0)| / M(0) for each body
- `'max_rel_drift'`: max(|M(t) - M(0)| / M(0)) over all timesteps

**Key Design Decisions:**

1. **Relative drift is the primary metric**: `|ΔM|/M(0)` is mass-independent and physically meaningful
2. **Max drift over entire trajectory**: Catches transient mass fluctuations, not just final drift
3. **Handles edge cases**:
   - Zero initial mass → returns `np.inf` for relative drift
   - Single timestep → drift = 0 by definition
   - Empty trajectory → returns empty arrays

### 2. Integration into `slab/run.py`

**Changes made:**

1. **Import statement** (line 56):
   ```python
   from slab.diagnostics import compute_precession, find_periapsis_passages, compute_mass_drift
   ```

2. **Call in `compute_summary()`** (lines 314-316):
   ```python
   # Mass drift diagnostics
   body_names = [body.name for body in bodies]
   mass_drift_data = compute_mass_drift(trajectory, body_names)
   ```

3. **Added to summary dict** (line 348):
   ```python
   'mass_drift': mass_drift_data,
   ```

4. **Summary printout in `print_summary()`** (lines 425-458):
   - Prints mass drift for each body with format: `Body_name: ΔM/M = X.XXe-XX`
   - Shows max drift over trajectory
   - Displays pass/fail status for each body:
     - **PASS**: drift < 1e-12 (negligible)
     - **WARN**: 1e-12 ≤ drift < 1e-9 (acceptable)
     - **FAIL**: drift ≥ 1e-9 (significant)
   - Overall status assessment

5. **JSON output**: Mass drift data automatically included in `diagnostics.json` via existing serialization

### 3. Pass/Fail Thresholds

**Per-body thresholds:**
- `drift < 1e-12`: **PASS** (negligible drift, expected for Solar System regime)
- `1e-12 ≤ drift < 1e-9`: **WARN** (acceptable but non-negligible)
- `drift ≥ 1e-9`: **FAIL** (significant drift, check `use_flux_mass` setting)

**Overall status thresholds:**
- `max_drift < 1e-12`: **EXCELLENT** (negligible drift)
- `1e-12 ≤ max_drift < 1e-9`: **GOOD** (acceptable drift)
- `1e-9 ≤ max_drift < 1e-6`: **WARNING** (non-negligible drift)
- `max_drift ≥ 1e-6`: **FAIL** (significant drift)

These thresholds ensure:
- Constant-mass systems (Solar System regime with `use_flux_mass=False`) should show **PASS**
- Variable-mass systems (with `use_flux_mass=True`) are distinguished from numerical errors
- Numerical roundoff errors (≈1e-15) are well below the **PASS** threshold

## Validation and Testing

### Test Suite: `test_mass_drift.py`

Created comprehensive test suite with 4 test categories:

1. **Constant-mass trajectory** (Solar System regime)
   - Validates drift = 0 for constant masses
   - Expected: All bodies **PASS**

2. **Variable-mass trajectory** (with intake)
   - Validates non-zero drift is correctly computed
   - Tests 10% mass increase scenario

3. **Edge cases**
   - Zero initial mass (returns `np.inf`)
   - Single timestep (drift = 0)
   - Empty trajectory (returns empty arrays)

4. **Pass/fail criteria**
   - Validates threshold logic at each level
   - Tests: 1e-15, 1e-13, 1e-10, 1e-7, 1e-3

**All tests PASS** ✓

### Realistic Simulation Test

Tested with `examples/mercury_quick_test.yaml` and `examples/audit_test_minimal.yaml`:

**Results:**
```
Mass conservation:
  Sun         : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Mercury     : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Status:             EXCELLENT (negligible drift)
```

This confirms:
- Mass is perfectly conserved in constant-mass regime (`use_flux_mass=False`)
- Diagnostic correctly reports negligible drift
- Output formatting is clear and informative

## Output Examples

### Console Output

```
Mass conservation:
  Sun         : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Mercury     : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Status:             EXCELLENT (negligible drift)
```

### JSON Output (`diagnostics.json`)

```json
{
  "summary": {
    "mass_drift": {
      "bodies": ["Sun", "Mercury"],
      "initial_mass": [1.0, 3.3e-07],
      "final_mass": [1.0, 3.3e-07],
      "abs_drift": [0.0, 0.0],
      "rel_drift": [0.0, 0.0],
      "max_rel_drift": [0.0, 0.0]
    }
  }
}
```

## Usage

### Command Line

```bash
# Standard run (shows mass drift in summary)
python -m slab.run examples/mercury_orbit.yaml --verbose

# Quick test (minimal steps)
python -m slab.run examples/mercury_quick_test.yaml --verbose

# Mass drift data saved to output/diagnostics.json
```

### Python API

```python
from slab.diagnostics import compute_mass_drift

# trajectory from integrate_orbit()
mass_drift = compute_mass_drift(trajectory, body_names)

# Access results
print(f"Relative drift: {mass_drift['rel_drift']}")
print(f"Max drift: {mass_drift['max_rel_drift']}")

# Check if drift is negligible
if mass_drift['max_rel_drift'].max() < 1e-12:
    print("Mass conserved (Solar System regime validated)")
```

## Design Rationale

### Why relative drift `|ΔM|/M(0)` instead of absolute drift?

- **Mass-independent**: Works for bodies with vastly different masses (Sun vs. planets)
- **Physically meaningful**: Fractional change is more interpretable than absolute change
- **Threshold-friendly**: Single threshold (1e-12) works across all mass scales

### Why track max drift over entire trajectory?

- **Catches transient errors**: Final drift might be small even if intermediate drift was large
- **More robust diagnostic**: Reveals numerical instabilities that cancel out
- **Conservative metric**: Ensures mass conservation throughout, not just at endpoints

### Why threshold at 1e-12?

- **Well above machine precision**: float64 has ~1e-16 precision, so 1e-12 has margin
- **Solar System regime expectation**: With `use_flux_mass=False`, masses never change → drift should be exactly 0 (or roundoff ~1e-15)
- **Distinguishes numerical errors from physics**: Real mass intake would have drift >> 1e-12

## Files Modified

1. **`slab/diagnostics.py`**: Added `compute_mass_drift()` function (136 lines)
2. **`slab/run.py`**: Integrated mass drift into pipeline (4 locations):
   - Import statement
   - Compute in `compute_summary()`
   - Add to summary dict
   - Print in `print_summary()`
   - CSV output bug fix (2 locations)
3. **`test_mass_drift.py`**: Test suite (206 lines)
4. **`MASS_DRIFT_IMPLEMENTATION.md`**: This documentation

## Verification

Run the test suite:
```bash
python test_mass_drift.py
```

Expected output:
```
ALL TESTS PASSED

The compute_mass_drift() diagnostic is working correctly!
```

Run a realistic simulation:
```bash
python -m slab.run examples/mercury_quick_test.yaml --verbose
```

Expected mass drift output:
```
Mass conservation:
  Sun         : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Mercury     : ΔM/M = 0.00e+00  (max: 0.00e+00)  [PASS]
  Status:             EXCELLENT (negligible drift)
```

## Summary

The mass drift diagnostic successfully validates that:
1. **Mass is conserved** in Solar System regime (constant-mass systems)
2. **Drift is negligible** (< 1e-12) as expected when `use_flux_mass=False`
3. **Implementation is robust** (handles edge cases, tested comprehensively)
4. **Output is clear** (pass/fail status, formatted values, included in JSON)

This diagnostic provides high confidence that the simulation correctly implements constant-mass dynamics in the Solar System regime, where mass intake rate Ṁ should be negligible.
