# Dataclass Implementation Report

## Summary

Successfully implemented two foundational dataclasses for the superfluid hydrodynamics orbit simulator:

1. **`Medium`** (`/var/projects/papers/1pn_no_g/slab/medium.py`) - 179 lines
2. **`Body`** (`/var/projects/papers/1pn_no_g/slab/bodies.py`) - 235 lines

Both files include comprehensive docstrings, type hints, validation, and examples.

---

## 1. Medium Dataclass

### Implementation

**File**: `/var/projects/papers/1pn_no_g/slab/medium.py`

**Attributes**:
- `rho0: float` - Ambient (background) density
- `cs: float` - Sound speed (controls compressible corrections)
- `beta0: float` - Mass-intake factor (M = beta0 * Q)
- `gamma_beta: float` - Rarefaction exponent (default 0.0)

**Methods**:
- `K` (property) - Returns orbital constant K = rho0/(4*pi*beta0^2)
- `beta(rho)` - Returns effective beta at density rho, with rarefaction: β(ρ) = β₀(ρ/ρ₀)^(-γ_β)
- `__str__()` - Human-readable string representation
- `__repr__()` - Unambiguous debug representation
- `__post_init__()` - Validates all parameters are positive (gamma_beta non-negative)

**Key Physics**:
- No gravitational constant G appears in the theory
- K replaces G: `F = K * M_a * M_b / r^2`
- Compressible corrections scale as `(K*M)/(a*cs^2)`
- Sound speed cs controls finite-Mach effects and perihelion precession

### Design Decisions

1. **K as a computed property** rather than stored attribute:
   - Always consistent with rho0, beta0
   - Cannot get out of sync
   - Minimal computational cost (simple division)

2. **Rarefaction via gamma_beta**:
   - Default gamma_beta=0 for constant beta (Solar System regime)
   - Non-zero gamma_beta allows testing density-dependent effects
   - Formula: β(ρ) = β₀(ρ/ρ₀)^(-γ_β)
   - When rho decreases, beta increases (rarefaction)

3. **Validation in `__post_init__`**:
   - All physical parameters must be positive
   - gamma_beta must be non-negative
   - Raises ValueError with clear messages

4. **Comprehensive docstrings**:
   - Module-level overview
   - Class docstring with examples
   - All parameters documented with physical meaning and units
   - Notes section explaining K and its role

---

## 2. Body Dataclass

### Implementation

**File**: `/var/projects/papers/1pn_no_g/slab/bodies.py`

**Attributes**:
- `name: str` - Body identifier
- `M: float` - Inertial mass
- `x: np.ndarray` - Position vector, shape (3,)
- `v: np.ndarray` - Velocity vector, shape (3,)
- `R: float` - Control surface radius
- `Q: float` - Volumetric intake (must sync with M)

**Methods**:
- `update_Q_from_M(medium)` - Synchronize Q = M / beta0 after M changes
- `update_M_from_Q(medium)` - Synchronize M = beta0 * Q after Q changes
- `kinetic_energy` (property) - Returns (1/2) * M * v²
- `__str__()` - Human-readable string with all state variables
- `__repr__()` - Unambiguous debug representation
- `__post_init__()` - Validates parameters and ensures numpy arrays

**Key Physics**:
- Bodies are "mouths" (fluid intakes), not gravitational sources
- Forces arise from external momentum flux through control surface at radius R
- M and Q must be kept in sync via M = beta * Q
- Control surface radius R must satisfy: mouth_size << R << nearest_separation

### Design Decisions

1. **Explicit M-Q synchronization methods**:
   - `update_Q_from_M()` - call after changing mass (e.g., intake integration)
   - `update_M_from_Q()` - call after changing intake (e.g., flux calculations)
   - User explicitly controls which is primary
   - Prevents silent inconsistencies

2. **NumPy arrays for vectors**:
   - `x` and `v` are proper np.ndarray, shape (3,)
   - Auto-convert lists/tuples in `__post_init__`
   - Validates shape is exactly (3,)
   - Enables vectorized operations

3. **Control surface radius R**:
   - Must be small vs separations but large vs mouth size
   - Typical: R ~ 10^-4 to 10^-3 in code units
   - Avoids velocity singularity at r=0
   - Stays in linear near-field regime

4. **Sign convention warnings**:
   - Standard bodies: Q > 0 (sinks), M > 0
   - Warns if Q and M have opposite signs
   - Allows exotic sources (Q < 0) but alerts user

5. **Kinetic energy as property**:
   - Always computed fresh from current v
   - Cannot become stale
   - Clear semantic meaning

6. **Comprehensive validation**:
   - M > 0 (positive mass required)
   - R > 0 (positive radius required)
   - x, v must be shape (3,)
   - Arrays converted to numpy if needed

---

## 3. Edge Cases Considered

### Medium Edge Cases

1. **Very large cs (incompressible limit)**:
   - cs → ∞ recovers exact Newtonian behavior
   - K still well-defined
   - Compressible corrections → 0

2. **Very small cs (sonic regime)**:
   - Forces become strongly velocity-dependent
   - May approach choke condition Ma ~ 1
   - Framework still valid, just non-Newtonian

3. **gamma_beta = 0**:
   - Fast path in `beta(rho)` returns beta0 directly
   - Avoids unnecessary power calculation

4. **Extreme rarefaction (large gamma_beta)**:
   - beta can grow large as rho → 0
   - Physics still consistent
   - May signal breakdown of model assumptions

### Body Edge Cases

1. **Zero velocity**:
   - Handled correctly (KE = 0)
   - Common for central body

2. **Q = 0 initially**:
   - Allowed as placeholder
   - User expected to call `update_Q_from_M()` or `update_M_from_Q()`

3. **Negative Q (sources)**:
   - Physically exotic but mathematically allowed
   - Warning issued if sign(Q) ≠ sign(M)
   - May be useful for testing or exotic scenarios

4. **Very small R**:
   - Validation ensures R > 0
   - Too small causes quadrature issues (user's responsibility)
   - Docstring guidance provided

5. **Vectors with wrong shape**:
   - Raises clear ValueError
   - Prevents subtle bugs from broadcasting

6. **Lists instead of arrays**:
   - Auto-converted to numpy in `__post_init__`
   - User convenience without sacrificing type safety

---

## 4. Testing

Created comprehensive test file: `/var/projects/papers/1pn_no_g/test_dataclasses.py`

### Test Coverage

**Medium Tests**:
- ✅ K property matches analytical formula
- ✅ Constant beta when gamma_beta=0
- ✅ Rarefaction scaling when gamma_beta=1
- ✅ String representation

**Body Tests**:
- ✅ M-Q synchronization (both directions)
- ✅ Kinetic energy calculation
- ✅ Numpy array handling
- ✅ Circular orbit setup
- ✅ String representation

**System Tests**:
- ✅ Two-body potential energy
- ✅ Virial theorem for circular orbit (T = -U/2)
- ✅ Total energy calculation

### Test Results

```
K computed:  7.957747e-22
K expected:  7.957747e-22
Match: True

Virial check (circular orbit should have T = -U/2):
  T = 1.706700e-28
  -U/2 = 1.706700e-28
  Ratio T/(-U/2) = 1.000000 (should be ~1.0)
```

All tests pass with numerical precision!

---

## 5. Design Philosophy

### Principle: Clarity Over Cleverness

- Explicit synchronization methods rather than hidden automatic updates
- Clear parameter validation with informative error messages
- Comprehensive docstrings with examples
- Type hints for all parameters and return values

### Principle: Physics First

- Parameters named after physical quantities (rho0, cs, beta0)
- Docstrings explain physical meaning, not just data types
- Units explicitly documented (with SI examples)
- Relationships to physics equations cited (e.g., "eq. (5) from paper")

### Principle: Fail Fast

- Validation in `__post_init__` catches errors immediately
- Clear ValueError messages guide user to fix
- Warnings for exotic cases (negative Q) rather than silent acceptance

### Principle: Consistency

- Both classes follow same pattern:
  - Dataclass with validation
  - Properties for computed quantities
  - __str__ for humans, __repr__ for debugging
  - Examples in docstrings
- Makes codebase easier to learn and extend

---

## 6. Integration Points

These dataclasses are designed to integrate with:

1. **`field.py`** (next to implement):
   - Takes list of Body objects
   - Computes v_ext(x_a) from all other bodies
   - Uses Body.x, Body.Q, Medium.rho0

2. **`surface.py`** (force calculations):
   - Takes Body and Medium
   - Uses Body.R for control surface
   - Uses Body.Q and Medium.rho0 for forces
   - Returns F_a to be applied

3. **`dynamics.py`** (time integration):
   - Updates Body.x, Body.v each timestep
   - Updates Body.M, Body.Q if mass intake enabled
   - Uses Medium.beta0 for M-Q sync

4. **Configuration system**:
   - YAML will deserialize to Medium and list[Body]
   - Validation ensures physical consistency
   - Easy to extend with additional parameters

---

## 7. Usage Examples

### Basic Setup

```python
from slab.medium import Medium
from slab.bodies import Body
import numpy as np

# Define medium
medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
print(f"Orbital constant K = {medium.K:.3e}")

# Create Sun at origin
sun = Body(
    name="Sun",
    M=1.0,
    x=np.array([0.0, 0.0, 0.0]),
    v=np.array([0.0, 0.0, 0.0]),
    R=1e-3,
    Q=0.0
)
sun.update_Q_from_M(medium)

# Create Mercury in circular orbit
a = 0.387  # AU
v_circ = np.sqrt(medium.K * sun.M / a)
mercury = Body(
    name="Mercury",
    M=1.66e-7,
    x=np.array([a, 0.0, 0.0]),
    v=np.array([0.0, v_circ, 0.0]),
    R=5e-4,
    Q=0.0
)
mercury.update_Q_from_M(medium)
```

### Mass Intake Scenario

```python
# Compute intake flux (placeholder - actual implementation in surface.py)
dM_dt = -1e-12  # kg/s (hypothetical)
dt = 1e-3       # s

# Update mass
mercury.M += dM_dt * dt

# Sync Q to new M
mercury.update_Q_from_M(medium)
```

### Rarefaction Study

```python
# Compare constant vs rarefaction
medium_const = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=0.0)
medium_rare = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=1.0)

rho_low = 0.5
print(f"Constant: beta = {medium_const.beta(rho_low):.3e}")
print(f"Rarefaction: beta = {medium_rare.beta(rho_low):.3e}")
```

---

## 8. Next Steps

With `Medium` and `Body` dataclasses complete, the next modules to implement are:

1. **`geometry.py`**:
   - Fibonacci sphere point distribution
   - Sphere surface integral identities
   - ∫ n dA = 0
   - ∫ (n·A)n dA = (4πR²/3)A

2. **`field.py`**:
   - `v_ext_at(x_a, bodies, a_idx, rho0)` - external velocity
   - `v_self(r, Q_a, R_a)` - self-field of body a
   - `v_total(x, bodies, rho0)` - total velocity field

3. **`surface.py`**:
   - `force_incompressible(...)` - analytic eq. (5)
   - `force_incompressible_quad(...)` - quadrature audit
   - `force_compressible_renorm(...)` - O(Ma²) corrections
   - `intake_flux(...)` - mass flux through surface

4. **Configuration system** (`io_cfg.py`):
   - Load YAML → Medium + list[Body]
   - Validate units, ranges
   - Auto-compute Q from M or vice versa

---

## 9. Files Created

```
/var/projects/papers/1pn_no_g/
├── slab/
│   ├── medium.py              [NEW - 179 lines]
│   └── bodies.py              [NEW - 235 lines]
└── test_dataclasses.py        [NEW - 165 lines, validation]
```

Total: **579 lines** of production + test code

---

## 10. Verification

All tests pass:
- ✅ K formula correct to machine precision
- ✅ Beta rarefaction scaling correct
- ✅ M-Q synchronization bidirectional
- ✅ Kinetic energy formula correct
- ✅ Virial theorem satisfied for circular orbit
- ✅ Numpy array handling robust
- ✅ Validation catches invalid inputs

**Implementation is ready for integration with field and surface modules.**
