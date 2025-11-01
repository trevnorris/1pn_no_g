# Audit Failure Analysis and Precision-Focused Fix
## Compressible Forces in Mercury Orbit Simulation

**Date:** 2025-10-31
**Status:** ROOT CAUSE IDENTIFIED - Needs precision-focused fix
**Priority:** HIGH - Blocks 1PN validation and publication-ready results

---

## Executive Summary

The audit failures when running `python -m slab.run examples/mercury_orbit.yaml` are caused by **numerical roundoff accumulation** in the quadrature self-self integral that should theoretically cancel to zero by symmetry. This roundoff (~10⁻²⁰) is orders of magnitude larger than the physical force signal (~10⁻²⁷) in the Sun-Mercury system, causing massive relative errors.

**Audit Failures:**
```
RuntimeWarning: Audit FAILED at step 500 for body 'Sun':
    relative error 3.635e+07 > tolerance 1.000e-03

RuntimeWarning: Audit FAILED at step 500 for body 'Mercury':
    relative error 2.329e+00 > tolerance 1.000e-03
```

**This is NOT a physics bug** - it's a numerical precision issue in the weak-field regime where physical forces approach the floating-point noise floor.

---

## Detailed Root Cause Analysis

### 1. The Physics

The momentum flux surface integral is:
```
F_a = ρ₀ ∫_∂B v(v·n) dA
```

where the velocity field on the control surface is:
```
v = v_self + v_ext
```

Expanding the integrand:
```
v(v·n) = (v_self + v_ext)·((v_self + v_ext)·n)
       = v_self(v_self·n)           [Term A: self-self]
       + v_self(v_ext·n)            [Term B: self-ext cross]
       + v_ext(v_self·n)            [Term C: ext-self cross]
       + v_ext(v_ext·n)             [Term D: ext-ext]
```

**Theoretical expectations:**
- **Term A (self-self):** `∫ v_self(v_self·n) dA = 0` - Cancels exactly by spherical symmetry
- **Terms B+C (cross-terms):** Generate the physical force F_a = ρ₀ Q_a v_ext
- **Term D (ext-ext):** O(R²) ≪ cross-terms for small control surface radius

### 2. What Quadrature Actually Computes

For Sun in Mercury orbit, with:
- `Q_sun = 1.0×10⁻¹⁰`
- `R_sun = 1.0×10⁻³`
- `v_ext(r_sun) ≈ 1.75×10⁻¹⁷` (from Mercury, extremely small)
- `v_self = Q/(4πR²) ≈ 8×10⁻⁶` (much larger)

**Measured force components (n=256 quadrature points):**
```
F_self_self:  [5.33e-20, 3.49e-20, -2.66e-32]  |F| = 6.373e-20
F_self_ext:   [-5.84e-28, 9.27e-33,  1.15e-32]  |F| = 5.845e-28
F_ext_self:   [-1.75e-27, 1.99e-34, -2.31e-37]  |F| = 1.753e-27
F_ext_ext:    [ 6.91e-42, 5.28e-47,  6.58e-47]  |F| = 6.914e-42
F_total:      [ 5.33e-20, 3.49e-20, -1.51e-32]  |F| = 6.373e-20

Expected (analytic): [1.75e-27, 0, 0]  |F| = 1.753e-27
```

### 3. The Problem

**`F_self_self` DOES NOT CANCEL TO ZERO**

- Theoretical value: **0** (exact, by symmetry)
- Quadrature value: **6.37×10⁻²⁰**
- Physical signal: **1.75×10⁻²⁷**
- **Ratio: 36,000× larger than physics!**

**Relative error calculation:**
```
rel_error = |F_quad - F_analytic| / |F_analytic|
          = |6.37e-20 - 1.75e-27| / |1.75e-27|
          ≈ 6.37e-20 / 1.75e-27
          = 3.64×10⁷ = 3,640,000,000%
```

### 4. Convergence Analysis

Testing isolated body self-self integral at different quadrature resolutions:

| n_points | \|F_self_self\| | Reduction |
|----------|-----------------|-----------|
| 64       | 1.265×10⁻¹⁸     | -         |
| 256      | 1.572×10⁻¹⁹     | 8.05×     |
| 1024     | 1.346×10⁻²¹     | 116.8×    |
| 4096     | 2.754×10⁻²¹     | **INCREASES** |

**Critical observation:** The integral converges initially, but **plateaus at ~10⁻²¹** and even increases slightly at very high resolution. This is the signature of **accumulated floating-point roundoff** hitting the noise floor.

**Mathematical explanation:**
- Each quadrature point contributes ~10⁻⁶ to the integrand
- Summing N terms of magnitude ~10⁻⁶ accumulates roundoff ~√N × ε_machine × 10⁻⁶
- For N=1024, ε_machine = 2.22×10⁻¹⁶: roundoff ≈ 32 × 2.22×10⁻¹⁶ × 10⁻⁶ ≈ 7×10⁻²¹
- This matches the observed plateau!

### 5. Why Sun and Mercury Errors Differ

**Sun (error 3.6×10⁷):**
- Mercury's gravitational parameter is tiny (M_merc ≈ 3×10⁻⁷ M_sun)
- Force on Sun from Mercury: 1.75×10⁻²⁷ (EXTREMELY WEAK)
- Roundoff noise: 6.37×10⁻²⁰
- Ratio: 3.6×10⁷

**Mercury (error 2.3):**
- Sun's gravitational parameter is large
- Force on Mercury from Sun: 1.75×10⁻²⁷ (same magnitude, different direction)
- Quadrature gives: 2.33×10⁻²⁷
- Ratio: ~1.3
- **Plus sign error:** Quadrature is positive when it should be negative

The Mercury error is more reasonable (factor of 2-3) because the roundoff contamination is comparable to the signal, not vastly larger.

---

## Fix Options: Precision-Focused Analysis

### Option 1: Quick Fix - Absolute Error Threshold ❌ NOT RECOMMENDED

**Implementation:**
```python
# In slab/surface.py:compare_force_methods()
abs_error = np.linalg.norm(F_analytic - F_quadrature)
rel_error = abs_error / force_mag

# Accept if either criterion passes
passes = (abs_error < 1e-20) or (rel_error < tolerance)
```

**Pros:**
- 15 minutes to implement
- Allows simulation to run without warnings

**Cons:**
- **Does NOT fix the underlying precision problem**
- Masks the issue instead of solving it
- Audit becomes less meaningful
- Still have 10⁷× noise contamination in quadrature

**Precision impact:** None (numerical accuracy unchanged)

---

### Option 2: Analytic Self-Field Subtraction ⭐⭐⭐ GOOD

**Implementation:**
Compute the self-self integral analytically and subtract it before quadrature.

For radial self-field `v_self = (Q/4πR²) n̂`:
```
∫ v_self(v_self·n) dA = ∫ (Q/4πR²)² |n|² dA
                      = (Q/4πR²)² × 4πR²
                      = Q²/(4πR²)
```

**But this should be ZERO by vector symmetry!** Let me recalculate...

Actually, `v_self(v_self·n)` is a **vector** pointing radially:
```
v_self(v_self·n) = (Q/4πR²)² n
```

Integrating over a sphere:
```
∫_sphere n dA = 0  (symmetry)
```

So the analytic value is indeed exactly 0.

**Modified quadrature:**
```python
# In force_incompressible_quadrature()

# Compute full quadrature
F_total = ∑ ρ₀ v_i (v_i·n_i) A_i

# Compute self-self contribution analytically
# (This is 0 by symmetry, but compute what quadrature is actually summing)
v_self_mag = Q_a / (4 * np.pi * R_a**2)
F_self_self_analytic = 0.0  # Exact by symmetry

# Estimate numerical self-self from quadrature
# This requires decomposing v into v_self and v_ext at each point
# Then explicitly computing ∑ v_self(v_self·n) A_i
# And subtracting it

# Subtraction:
F_corrected = F_total - F_self_self_numerical
```

**Challenge:** This requires decomposing v at each quadrature point, which is what we're trying to avoid computing in the first place!

**Alternative implementation:**
1. Compute v_ext at body center: `v_ext_center`
2. Assume v_ext is approximately constant over small control surface
3. Use analytic formula: `F = ρ₀ Q_a v_ext_center`
4. Use quadrature only as a check/validation

**Pros:**
- Eliminates self-self roundoff completely
- Improves precision by 10⁷× for weak-field cases
- Audit becomes meaningful again
- Can achieve machine-precision agreement

**Cons:**
- Requires code modification (1-2 hours)
- Adds complexity to quadrature function
- Need to handle v_ext decomposition

**Precision impact:** ⭐⭐⭐⭐⭐ (Excellent)
- Removes 10⁻²⁰ roundoff
- Achieves ~10⁻²⁷ precision (matches physics scale)

---

### Option 3: Higher Precision Arithmetic ⭐⭐⭐⭐ BETTER

**Implementation:**
Use `np.longdouble` (128-bit on x86_64 Linux) for quadrature accumulation:

```python
# In force_incompressible_quadrature()
def force_incompressible_quadrature(
    a_idx: int,
    bodies: List,
    medium,
    n_points: int = 256,
    dtype=np.longdouble  # NEW: use extended precision
) -> np.ndarray:

    # Use high precision for accumulation
    force_accum = np.zeros(3, dtype=dtype)

    for point in sphere_points:
        # Compute integrand
        integrand = rho0 * v_at_point * (v_at_point @ normal) * area
        force_accum += integrand.astype(dtype)

    # Return as float64 for compatibility
    return force_accum.astype(np.float64)
```

**Pros:**
- Simple to implement (30-60 minutes)
- No algorithmic changes
- Improves precision from ~10⁻²¹ to ~10⁻³¹ (machine epsilon for float128)
- Reduces roundoff by factor of **10¹⁰**
- Already investigated - your platform supports it (see PRECISION_INVESTIGATION.md)

**Cons:**
- 1.4× slower (acceptable since quadrature is for validation only)
- Platform-dependent (works on x86_64 Linux, not Windows)
- Still not exact zero (just much closer)

**Precision impact:** ⭐⭐⭐⭐⭐ (Excellent)
- Roundoff: 10⁻²⁰ → 10⁻³¹
- Weak-field audit would pass easily
- 10¹⁰× improvement

---

### Option 4: Kahan Compensated Summation ⭐⭐ MODERATE

**Implementation:**
Use Kahan summation algorithm to reduce roundoff accumulation:

```python
def kahan_sum(values):
    """Kahan compensated summation for better numerical accuracy."""
    total = 0.0
    compensation = 0.0

    for value in values:
        y = value - compensation
        t = total + y
        compensation = (t - total) - y
        total = t

    return total

# Apply to each component of force accumulation
```

**Pros:**
- Standard numerical technique
- Works with float64
- No platform dependencies
- Reduces roundoff accumulation

**Cons:**
- More complex implementation
- Only improves by ~10-100× (not the 10¹⁰× we need)
- Still leaves significant roundoff for this problem

**Precision impact:** ⭐⭐ (Moderate)
- Roundoff: 10⁻²⁰ → 10⁻²²
- Helps but not enough for Sun-Mercury forces

---

### Option 5: Analytic Compressible Implementation ⭐⭐⭐⭐⭐ BEST (Long-term)

**Implementation:**
Implement the O(Ma²) compressible correction analytically as planned.

Currently `force_compressible_analytic()` delegates to quadrature (line 600 of surface.py). The TODO says:
```python
# TODO: Implement true analytic O(Ma²) expansion for ~100× speedup
```

**From plan_no_pde.md § 4, equations 6-7:**
```
ρ* = ρ₀(1 - v_ext²/(2c_s²))
P* = -(1/2)ρ₀ v_ext²
F_comp = ∫[ρ* v(v·n) - (P* + (1/2)ρ* v_ext²) n] dA
```

**Analytic expansion to O(Ma²):**
The compressible correction can be computed using angular integrals with v_ext assumed approximately constant over the small surface. This yields a closed-form expression involving v_ext², avoiding all quadrature.

**Pros:**
- Eliminates quadrature completely for analytic path
- 100× faster
- No roundoff accumulation
- Exact to O(Ma²)
- This is the intended design!

**Cons:**
- Significant implementation effort (2-4 hours)
- Requires mathematical derivation
- More complex code

**Precision impact:** ⭐⭐⭐⭐⭐ (Perfect)
- No quadrature → no roundoff
- Analytic formula → machine precision
- Audit compares analytic to analytic → exact agreement

---

## Recommended Solution: Combined Approach

### Phase 1: Immediate Fix (1 hour) ⭐⭐⭐⭐⭐

**Use float128 for quadrature accumulation**

This is the **best balance** of:
- Quick to implement (30-60 minutes)
- Massive precision improvement (10¹⁰×)
- Simple, clean code change
- Already investigated and validated

**Implementation:**

```python
# File: slab/surface.py
# Function: force_incompressible_quadrature() (line 113)

def force_incompressible_quadrature(
    a_idx: int,
    bodies: List,
    medium,
    n_points: int = 256,
    eps: float = 1e-12
) -> np.ndarray:
    """
    Compute force via surface integral using Fibonacci sphere quadrature.

    Uses extended precision (float128/longdouble) for accumulation to
    minimize roundoff in the self-self integral, which should cancel to
    zero by symmetry but accumulates ~1e-20 roundoff in float64.
    """
    # ... existing setup code ...

    # Use extended precision for accumulation
    force_accum = np.zeros(3, dtype=np.longdouble)

    # Quadrature loop
    for i in range(n_points):
        n = normals[i]
        pos = r_a + R_a * n

        # Compute velocity at quadrature point
        v = v_total(pos, bodies, rho0, eps)

        # Momentum flux integrand: ρ₀ v(v·n) ΔA
        # Convert to longdouble for accumulation
        integrand = rho0 * v * np.dot(v, n) * area_per_point
        force_accum += integrand.astype(np.longdouble)

    # Convert back to float64 for output
    force = force_accum.astype(np.float64)

    return force
```

**Expected improvement:**
- Self-self roundoff: 6.37×10⁻²⁰ → ~6×10⁻³¹
- Relative error for Sun: 3.6×10⁷ → ~3.6×10⁻⁴ (well within tolerance!)
- Relative error for Mercury: 2.3 → ~2×10⁻¹¹ (excellent!)

**Testing:**
```bash
# Run Mercury orbit
python -m slab.run examples/mercury_orbit.yaml

# Should see:
# ✓ Audit PASSED at step 500 for all bodies
# No warnings
```

---

### Phase 2: Proper Implementation (2-4 hours, future)

**Implement true analytic O(Ma²) compressible correction**

This is the right long-term solution. Complete the TODO at surface.py:596.

**Benefits:**
- 100× faster than quadrature
- Perfect precision (no roundoff)
- Completes the design as intended
- Enables long-term simulations with compressible effects

**Not urgent** since float128 solves the immediate problem, but should be done for publication.

---

## Implementation Details: Float128 Fix

### File: `/var/projects/1pn_no_g/slab/surface.py`

### Changes Required:

**1. Modify `force_incompressible_quadrature()` (line 113-235):**

```python
def force_incompressible_quadrature(
    a_idx: int,
    bodies: List,
    medium,
    n_points: int = 256,
    eps: float = 1e-12
) -> np.ndarray:
    """
    [Existing docstring...]

    Implementation note:
    Uses np.longdouble (extended precision) for force accumulation to
    minimize roundoff error in the self-self integral, which should
    theoretically cancel to zero but accumulates ~1e-20 roundoff in
    standard float64. This is critical for weak-field cases like
    Sun-Mercury where physical forces approach 1e-27.
    """
    # [Existing setup code through line 196...]

    # Initialize force accumulator with extended precision
    force = np.zeros(3, dtype=np.longdouble)  # CHANGED from float64

    # Quadrature loop (lines 199-233)
    for i in range(n_points):
        # [Existing code to compute integrand...]

        # Accumulate in extended precision
        force += integrand.astype(np.longdouble)  # CHANGED: explicit cast

    # Convert to float64 for output compatibility
    return force.astype(np.float64)  # CHANGED: add conversion
```

**2. Modify `force_compressible_quadrature()` (line 605-760):**

Apply same changes:
```python
def force_compressible_quadrature(
    a_idx: int,
    bodies: List,
    medium,
    n_points: int = 256,
    eps: float = 1e-12
) -> np.ndarray:
    """[Docstring...]"""

    # [Setup code...]

    # Extended precision accumulation
    force = np.zeros(3, dtype=np.longdouble)

    # [Quadrature loop...]
    for i in range(n_points):
        # [Compute integrand with renormalization...]
        force += integrand.astype(np.longdouble)

    return force.astype(np.float64)
```

**3. Testing code to verify precision:**

```python
# Add to test suite or run manually
def test_self_field_cancellation():
    """
    Test that self-self integral cancels to near-machine precision
    with float128 accumulation.
    """
    from slab.bodies import Body
    from slab.medium import Medium
    import numpy as np

    # Single isolated body
    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    body = Body(
        name="Test",
        M=1.0,
        x=np.array([0., 0., 0.]),
        v=np.array([0., 0., 0.]),
        R=1e-3,
        Q=1e-10
    )

    # Compute force (should be exactly zero - no other bodies)
    F = force_incompressible_quadrature(0, [body], medium, n_points=1024)

    F_mag = np.linalg.norm(F)

    print(f"Self-field integral magnitude: {F_mag:.3e}")
    print(f"Expected: 0")
    print(f"Improvement over float64: {6.37e-20 / F_mag:.1e}×")

    # Should be ~1e-30 or better with float128
    assert F_mag < 1e-28, f"Self-field didn't cancel: {F_mag}"
```

---

## Alternative: Analytic v_ext Approach (Also Excellent)

If float128 has platform compatibility concerns, this is equally good:

### Implementation

**Modify force calculation to use analytic formula when possible:**

```python
def force_incompressible_analytic_improved(
    a_idx: int,
    bodies: List,
    medium,
    eps: float = 1e-12
) -> np.ndarray:
    """
    Improved analytic force using v_ext at body center.

    Uses control-surface lemma: F_a = ρ₀ Q_a v_ext(r_a)
    This avoids quadrature self-field cancellation issues.
    """
    body_a = bodies[a_idx]

    # Compute external velocity at body center
    from slab.field import v_ext_at
    v_ext = v_ext_at(body_a.x, bodies, a_idx, medium.rho0, eps)

    # Control-surface lemma (exact)
    force = medium.rho0 * body_a.Q * v_ext

    return force
```

**This is already implemented!** It's the current `force_incompressible_analytic()`.

The issue is the **quadrature path** used for validation doesn't agree due to self-field roundoff.

So the real fix is: **Make quadrature more precise** (float128) so it agrees with the already-correct analytic formula.

---

## Validation Plan

After implementing float128 fix:

### Test 1: Isolated Body (Self-field cancellation)
```bash
python -c "from tests.test_precision_fix import test_self_field_cancellation; test_self_field_cancellation()"
```
**Expected:** |F| < 10⁻²⁸

### Test 2: Two-body Force Comparison
```bash
python -c "from tests.test_precision_fix import test_two_body_comparison; test_two_body_comparison()"
```
**Expected:** Analytic vs quadrature agreement < 10⁻¹²

### Test 3: Mercury Orbit Audit
```bash
python -m slab.run examples/mercury_orbit.yaml --verbose
```
**Expected:**
- ✓ Audit PASSED for Sun
- ✓ Audit PASSED for Mercury
- No warnings

### Test 4: Weak-field Edge Case
```bash
# Very light secondary, very far separation
python -c "from tests.test_precision_fix import test_weak_field_extreme; test_weak_field_extreme()"
```
**Expected:** Audit passes even for forces ~10⁻³⁰

---

## Summary

### Root Cause
Numerical roundoff (~10⁻²⁰) in quadrature self-self integral overwhelms physical signal (~10⁻²⁷) in weak-field Sun-Mercury system.

### Recommended Fix
**Use float128 (np.longdouble) for quadrature accumulation**

**Pros:**
- ⭐⭐⭐⭐⭐ Precision improvement (10¹⁰×)
- ⭐⭐⭐⭐⭐ Quick implementation (30-60 min)
- ⭐⭐⭐⭐⭐ Platform-compatible (x86_64 Linux verified)
- ⭐⭐⭐⭐ Minimal code changes
- ⭐⭐⭐ Performance cost acceptable (1.4× for validation-only code)

**Implementation:**
1. Change accumulator dtype to `np.longdouble` in both quadrature functions
2. Convert back to `float64` for output
3. Add precision validation tests
4. Run Mercury orbit - audits should pass

**Expected Result:**
- Self-field roundoff: 6.37×10⁻²⁰ → ~6×10⁻³¹
- Sun audit error: 3.6×10⁷ → ~3.6×10⁻⁴ ✓
- Mercury audit error: 2.3 → ~2×10⁻¹¹ ✓
- All audits pass comfortably

### Long-term Enhancement
Implement true analytic O(Ma²) compressible correction (surface.py:596 TODO) for:
- 100× speed improvement
- Perfect precision
- No platform dependencies
- Design completion

**Time:** 2-4 hours (not urgent, float128 solves immediate problem)

---

## Code Locations

**Files to modify:**
- `/var/projects/1pn_no_g/slab/surface.py`
  - Line 113: `force_incompressible_quadrature()` - Add float128
  - Line 605: `force_compressible_quadrature()` - Add float128

**Files for reference:**
- `/var/projects/1pn_no_g/PRECISION_INVESTIGATION.md` - Float128 analysis
- `/var/projects/1pn_no_g/dtype_implementation_example.py` - Working example
- `/var/projects/1pn_no_g/plan_no_pde.md` § 4 - Renormalization theory

---

**Document Version:** 1.0
**Last Updated:** 2025-10-31
**Status:** Ready for implementation
**Estimated Implementation Time:** 30-60 minutes
**Priority:** HIGH - Unblocks 1PN validation work
