# Superfluid Orbit Simulator

**Emergent Newtonian Gravity and 1PN Effects from Superfluid Hydrodynamics**

This project implements a computational simulator that derives gravitational orbital dynamics from **pure superfluid hydrodynamics** without invoking Newton's gravitational constant G. Bodies are modeled as fluid intakes ("mouths") that remove superfluid at volumetric rates Q, and forces emerge from momentum flux through control surfaces.

## Overview

### The Physics

This simulator implements the theory from the paper `1pn_no_g.tex`:

- **No gravitational field**: Bodies don't have gravity; they are sinks in a superfluid
- **Forces from momentum flux**: F_a = ρ₀ Q_a v_ext(r_a) from surface integrals
- **Mass-intake relation**: M_a = β Q_a (universal mapping)
- **Emergent 1/r² law**: Newtonian gravity appears naturally without G
- **Orbital constant**: K = ρ₀/(4πβ²) replaces G in all formulas
- **1PN-like effects**: Finite sound speed c_s produces perihelion precession

### Key Equations

**Velocity field** (potential flow with point sinks):
```
v(x) = - Σ (Q_b/4π) * r_b/r_b³
```

**Force** (control surface lemma):
```
F_a = ρ₀ Q_a v_ext(r_a)
```

**Newtonian form** (after substituting M = βQ):
```
F_a = K M_a M_b / r²
where K = ρ₀/(4πβ²)
```

**Energy** (conserved for incompressible):
```
E = Σ (1/2) M_a v_a² - K Σ_{a<b} M_a M_b / r_ab
```

## Installation

### Requirements

```bash
pip install numpy pyyaml
```

Optional for testing:
```bash
pip install pytest  # If you want to run pytest-based tests
```

### Quick Start

1. Clone or navigate to this directory:
```bash
cd /var/projects/papers/1pn_no_g
```

2. The simulator is implemented as a Python package in the `slab/` directory

## Project Structure

```
.
├── README.md                    # This file
├── PROJECT.md                   # Detailed project tracking
├── 1pn_no_g.tex                 # Physics paper (theory)
├── plan_no_pde.md               # Implementation plan
│
├── slab/                        # Main simulator package
│   ├── __init__.py
│   ├── geometry.py              # Fibonacci sphere, surface integrals
│   ├── medium.py                # Medium(rho0, cs, beta0) dataclass
│   ├── bodies.py                # Body(M, x, v, Q, R) dataclass
│   ├── field.py                 # Velocity field calculations
│   ├── surface.py               # Force calculations (momentum flux)
│   ├── dynamics.py              # Velocity-Verlet integrator
│   ├── diagnostics.py           # Energy, momentum, orbital elements
│   ├── io_cfg.py                # Configuration loading/validation
│   ├── gr1pn.py                 # GR-1PN comparator (for validation)
│   └── run.py                   # Main CLI entry point
│
├── examples/
│   └── mercury_orbit.yaml       # Example Sun-Mercury configuration
│
├── tests/
│   ├── test_*.py                # Unit tests created during development
│   └── test_full_simulation.py  # Comprehensive integration test
│
└── debug_force.py               # Force formula verification script
```

## How to Run

### 1. Main Simulator (CLI)

Run a simulation from a YAML configuration file:

```bash
python -m slab.run examples/mercury_orbit.yaml
```

**Options:**
```bash
python -m slab.run CONFIG.yaml [options]

Options:
  --output-dir DIR      Output directory (default: output/)
  --verbose, -v         Enable detailed logging
  --quick               Disable quadrature audits for speed
  --validate-only       Just validate config and exit
  --no-table            Skip trajectory table output
  --table-rows N        Number of table rows to display
```

**Example with verbose output:**
```bash
python -m slab.run examples/mercury_orbit.yaml --verbose --output-dir results/
```

**Validate configuration only:**
```bash
python -m slab.run examples/mercury_orbit.yaml --validate-only
```

### 2. Create Example Configuration

Generate a template configuration file:

```bash
python -m slab.io_cfg create_example my_orbit.yaml
```

Edit `my_orbit.yaml` to customize:
- Medium parameters (ρ₀, c_s, β₀)
- Body masses, positions, velocities
- Integration settings (dt, n_steps)
- Output options

### 3. Run Tests

**Comprehensive integration test:**
```bash
python test_full_simulation.py
```

**Force formula verification:**
```bash
python debug_force.py
```

**Individual module tests** (if pytest available):
```bash
pytest tests/
```

### 4. Python API Usage

Use the simulator programmatically:

```python
import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit

# Define superfluid medium
medium = Medium(rho0=1.0, cs=1.0e4, beta0=1.0e10)
print(f"Orbital constant K = {medium.K:.3e}")  # Replaces G

# Create bodies
sun = Body(
    name="Sun",
    M=1.0,
    x=np.array([0.0, 0.0, 0.0]),
    v=np.array([0.0, 0.0, 0.0]),
    R=1e-3,
    Q=0.0
)
sun.update_Q_from_M(medium)  # Sync Q = M/beta0

# Set up circular orbit
a = 1.0  # semi-major axis
v_circ = np.sqrt(medium.K * sun.M / a)
planet = Body(
    name="Planet",
    M=3e-6,
    x=np.array([a, 0.0, 0.0]),
    v=np.array([0.0, v_circ, 0.0]),
    R=5e-4,
    Q=0.0
)
planet.update_Q_from_M(medium)

# Run simulation
trajectory, diagnostics = integrate_orbit(
    bodies=[sun, planet],
    medium=medium,
    dt=1e9,           # timestep
    n_steps=1000,     # number of steps
    opts={
        'use_compressible': False,
        'use_quadrature': False,
        'save_every': 10
    }
)

# Analyze results
from slab.diagnostics import total_energy, energy_drift_monitor

E_drift = energy_drift_monitor(trajectory, medium)
print(f"Energy drift: {E_drift['dE_rel']:.3e}")
```

## Output Files

After running a simulation, you'll find:

### `output/trajectory.csv`
Comma-separated trajectory data:
```
time,body_name,x,y,z,vx,vy,vz,M,Q
0.0,Sun,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1e-10
0.0,Planet,1.0,0.0,0.0,0.0,2.82e-11,0.0,3e-06,3e-16
...
```

### `output/diagnostics.json`
JSON with simulation diagnostics:
```json
{
  "times": [...],
  "total_energy": [...],
  "kinetic_energy": [...],
  "potential_energy": [...],
  "max_force": [...],
  "summary": {
    "timing": {...},
    "energy": {...},
    "forces": {...}
  }
}
```

## Current Status

### ✅ Completed

- [x] All core modules implemented (10 files, ~5000 lines)
- [x] Fibonacci sphere quadrature for surface integrals
- [x] Medium and Body dataclasses with M↔Q synchronization
- [x] Velocity field calculations (vectorized for performance)
- [x] Force calculations via momentum flux surface integrals
- [x] Velocity-Verlet symplectic integrator
- [x] Comprehensive diagnostics (energy, momentum, orbital elements)
- [x] YAML configuration system with validation
- [x] GR-1PN comparator module
- [x] Command-line interface
- [x] Force formula verified correct

### 🐛 Current Issue: Integration Stability

**Problem:** The comprehensive test (`test_full_simulation.py`) shows orbits spiraling outward with large energy drift.

**Test Setup:**
- 2-body system: Sun (M=1.0) + Planet (M=3×10⁻⁶)
- Circular orbit at r=1.0
- Timestep: dt = 0.005 × T_orbit = 1.11×10⁹
- Integration: 3 orbits (600 steps)

**Test Results:**
```
================================================================================
TEST RESULTS SUMMARY
================================================================================

    ⚠  SOME TESTS FAILED - SEE DETAILS BELOW  ⚠

--------------------------------------------------------------------------------
1. ENERGY CONSERVATION
--------------------------------------------------------------------------------
  Status: ✗ FAIL
  Initial energy:       E0 = -1.193662e-27
  Maximum drift:  |ΔE/E| = 3.87e+00  (tolerance: 1e-05)
  RMS drift:      |ΔE/E| = 3.52e+00

--------------------------------------------------------------------------------
2. MOMENTUM CONSERVATION
--------------------------------------------------------------------------------
  Status: ✓ PASS
  Initial momentum: p0 = [+0.00e+00, +8.46e-17, +0.00e+00]
  Maximum drift:   |Δp| = 4.91e-31  (tolerance: 1e-12)

--------------------------------------------------------------------------------
3. FORCE LAW (1/r² with correct coefficient)
--------------------------------------------------------------------------------
  Status: ✗ FAIL
  Measured force:   F = 2.409769e-30
  Expected force:   F = 2.387324e-27
  Relative error:   ε = 9.99e-01  (tolerance: 1e-10)
  1/r² scaling: F(a)/F(2a) = 3.999994  (expected: 4.0)  ← This works!
  Scaling error:        ε = 1.41e-06
  Newton 3rd law: |F1+F2| = 3.92e-46                     ← This works!

--------------------------------------------------------------------------------
4. ORBIT STABILITY (circular orbit)
--------------------------------------------------------------------------------
  Status: ✗ FAIL
  Expected radius:   a = 1.000000
  Mean radius:   r_avg = 15.538975
  Std deviation: r_std = 9.17e+00
  Min radius:    r_min = 1.000000
  Max radius:    r_max = 31.475162
  Eccentricity:      e ≈ 9.38e-01  (tolerance: 0.01)

--------------------------------------------------------------------------------
5. NUMERICAL STABILITY (no NaN/inf)
--------------------------------------------------------------------------------
  Status: ✓ PASS
  Contains NaN: False
  Contains inf: False
```

**Force Formula Verification (debug_force.py):**
```
Expected force magnitude: 2.387324e-27
Force on body2 (current formula): 2.387324e-27

======================================================================
TESTING DIFFERENT FORMULAS:
======================================================================
4. F = ρ₀ * Q * v_ext:              2.387324e-27  (ratio: 1.00) ✓ CORRECT

Target ratio should be 1.00
```

**Observations:**
1. ✅ Force formula is mathematically correct (verified independently)
2. ✅ 1/r² scaling works perfectly (4.0 ratio for distance doubling)
3. ✅ Newton's 3rd law satisfied to machine precision
4. ✅ Momentum conserved perfectly
5. ❌ Energy not conserved during integration
6. ❌ Orbit spirals outward instead of remaining circular

**Possible Causes:**
1. **Timestep too large**: dt = 0.005 T might be insufficient
   - Symplectic integrators need dt << T for energy conservation
   - Try dt = 0.001 T or smaller

2. **Initial conditions**: Circular orbit velocity might be off
   - v_circ = √(KM/r) needs machine precision
   - Rounding errors could cause drift

3. **Force magnitude discrepancy during integration**:
   - Static test shows force is correct
   - But measured force during integration is 1000× too small
   - Suggests forces might not be getting applied correctly

4. **Subtle bug in integration loop**:
   - Verlet algorithm looks correct (verified against textbooks)
   - But could be an indexing or state update issue

**What Works vs. What Doesn't:**

| Component | Status | Evidence |
|-----------|--------|----------|
| Force formula (static) | ✅ Works | debug_force.py shows ratio=1.00 |
| 1/r² scaling | ✅ Works | Test shows 4.0 ratio |
| Newton's 3rd law | ✅ Works | F₁+F₂ ≈ 0 to machine precision |
| Momentum conservation | ✅ Works | Drift < 10⁻³¹ |
| Force during integration | ❌ Fails | Measured force 1000× too small |
| Energy conservation | ❌ Fails | Drift 387% |
| Orbit stability | ❌ Fails | Spirals outward |

## Next Steps

### Immediate Debugging

1. **Test with smaller timestep**:
   ```bash
   # Edit test_full_simulation.py: change dt_factor from 0.005 to 0.0001
   python test_full_simulation.py
   ```

2. **Add detailed force logging**:
   - Print force magnitude at each timestep
   - Verify forces are being computed and applied correctly
   - Check if forces suddenly drop to near-zero

3. **Verify circular orbit setup**:
   - Print initial v_circ calculation
   - Check that it exactly matches √(KM/r) to full precision
   - Verify virial theorem: 2T + U = 0 for circular orbit

4. **Test with simpler parameters**:
   - Use M₁ = M₂ = 1.0 (equal masses)
   - Use r = 1.0, ρ₀ = 1.0, β₀ = 1.0 (all unity)
   - Minimize numerical errors

### Physics Validation Tasks

5. **Run Mercury orbit example**:
   ```bash
   python -m slab.run examples/mercury_orbit.yaml --verbose
   ```

6. **Compare with GR-1PN**:
   - Use gr1pn.py to compute expected precession
   - Run parallel simulations with matched parameters
   - Verify slab reproduces GR predictions

7. **Test quadrature audit mode**:
   - Compare analytic vs. quadrature forces
   - Should agree to < 10⁻³ per plan acceptance criteria

### Implementation Tasks

8. **Complete compressible forces** (currently disabled):
   - Implement surface.py compressible correction
   - Add O(Ma²) terms from finite c_s
   - Test perihelion precession scaling ∝ c_s⁻²

9. **Add visualization**:
   - Plot orbits (x-y trajectories)
   - Plot energy vs. time
   - Plot precession angle vs. time

10. **Write comprehensive tests**:
    - Unit tests for each module
    - Integration tests for known solutions
    - Regression tests for bug fixes

## Scientific Validation Checklist

From `plan_no_pde.md` acceptance criteria:

- [ ] **Test 10.1**: Emergent 1/r² with correct coefficient (< 0.5% error)
- [ ] **Test 10.2**: Orbit stability (|Δa|/a < 10⁻⁵ over 50 orbits)
- [ ] **Test 10.3**: Compressible correction ∝ c_s⁻² (10% slope accuracy)
- [ ] **Test 10.4**: Quadrature audit (< 10⁻³ error)
- [ ] **Test 10.5**: GR comparison (Mercury precession within few %)

## References

### Papers & Documentation

- **Theory**: `1pn_no_g.tex` - Full mathematical derivation
- **Implementation plan**: `plan_no_pde.md` - Numerical scheme details
- **Project tracking**: `PROJECT.md` - Development log

### Key Physics Papers

- Soffel (1989), *Relativity in Astrometry*, Eq. (6.82) - 1PN formula
- Will & Wiseman (1996), Phys. Rev. D 54, 4813 - EIH equations
- Blanchet (2014), *Gravitational Radiation* - Post-Newtonian sources

## Contact & Collaboration

This is original research demonstrating that:
- Newtonian gravity can emerge from superfluid hydrodynamics
- 1PN effects appear without invoking GR
- No gravitational constant G is needed
- Orbital constant K = ρ₀/(4πβ²) replaces G universally

The simulator is a computational proof-of-concept for the theory in `1pn_no_g.tex`.

## License

[Specify license if applicable]

---

**Last Updated**: 2025-10-31
**Status**: Core implementation complete, debugging integration stability issue
**Version**: 0.1.0-alpha
