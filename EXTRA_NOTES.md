# EXTRA.md Implementation Analysis
## Prioritized Action Items for Superfluid Orbit Simulator

**Date:** 2025-10-31
**Purpose:** Analysis of suggestions from EXTRA.md to strengthen scientific validation and communication
**Source:** EXTRA.md - Ideas for making the data "hit home"

---

## Executive Summary

EXTRA.md provides 7 sections of suggestions to make the superfluid gravity simulator bulletproof against criticism and publication-ready. This document prioritizes those suggestions based on:
- Scientific impact
- Implementation effort
- Current project status

**Key Finding:** Several critical items are quick wins (~8 hours total) that dramatically strengthen the case.

---

## Current Implementation Status

### ✅ Already Implemented

1. **Audit parity system** (Section 1, item 3)
   - Code exists in `slab/dynamics.py`
   - Currently comparing analytic vs quadrature forces
   - **Status:** Working but showing warnings (needs investigation)

2. **Compressible force corrections** (Section 6, item 1)
   - Implementation complete in `slab/surface.py`
   - Just integrated into dynamics loop
   - **Status:** Running but audit failing (needs debugging)

3. **Energy/momentum conservation** (Section 5, item 4)
   - Diagnostics module tracks both
   - Tests validate conservation to 2.4×10⁻⁷
   - **Status:** Excellent, working perfectly

4. **Step invariance** (Section 2, item 1 - partial)
   - Symplectic integrator inherently provides this
   - **Status:** Implicit but not explicitly demonstrated

### 🔧 Partially Implemented

5. **No-G hygiene** (Section 1, item 5)
   - Code uses M internally, K is derived
   - **Need:** Config templates showing μ_a approach
   - **Effort:** ~1 hour (documentation)

### ❌ Not Yet Implemented

6. **Emergent inverse-square validation** (Section 1, item 1)
7. **Force decomposition logging** (Section 1, item 2)
8. **1PN diagnostics with GR overlay** (Section 1, item 4)
9. **Robustness tests** (Section 2, all items)
10. **Story-telling figures** (Section 3, most items)
11. **Drop-in README text** (Section 4)
12. **Release checklist items** (Section 7)

---

## HIGH IMPACT Recommendations (Implement Soon)

### Priority 1: Quick Wins (~8 hours total)

These items provide maximum scientific credibility for minimal effort:

#### 1. Emergent Inverse-Square Scan ⭐⭐⭐⭐⭐

**From:** Section 1, item 1
**Purpose:** Prove F ∝ 1/r² emerges from surface integrals, not hard-coded
**Effort:** ~2 hours

**Implementation:**
```python
# scripts/validate_force_law.py
"""
Static two-body scan across separations r ∈ [0.1, 10] AU
- Fix Q1, Q2, ρ₀
- Compute force via surface integral for each r
- Fit F = C/r²
- Compare to theoretical: C = ρ₀|Q₁Q₂|/(4π)
- Report relative error: |C_num/C_theory - 1|
"""
```

**Deliverable:**
- Script: `scripts/validate_force_law.py`
- Plot: Log-log F vs r with slope = -2
- Report: C_num/C_theory - 1 < 0.5%

**Impact:** Directly refutes "you hard-coded Newton" critique

**Note:** EXTRA.md uses old coefficient formula with 4/3 factor. Current correct formula is:
```
C_theory = ρ₀|Q₁Q₂|/(4π)
```

---

#### 2. Audit Parity Plot ⭐⭐⭐⭐

**From:** Section 3, item 4
**Purpose:** Show analytic and quadrature methods agree
**Effort:** ~1 hour (data already logged)

**Implementation:**
```python
# In slab/viz.py
def plot_audit_parity(diagnostics: Dict, output_path: str):
    """
    Plot analytic vs quadrature force comparison over time.
    - X-axis: timestep or time
    - Y-axis: relative error |F_analytic - F_quadrature| / |F_analytic|
    - Horizontal line at tolerance (1e-3)
    - Color-code by body
    """
```

**Deliverable:**
- Function in `slab/viz.py`
- Plot showing agreement < 10⁻³ throughout simulation
- Add to standard output plots

**Impact:** Demonstrates numerical reliability

**Current Issue:** Compressible mode showing large errors (2.3 for Mercury, 3.6×10⁷ for Sun)
**Action Required:** Debug before plotting (see Section below)

---

#### 3. Force Decomposition Bars ⭐⭐⭐

**From:** Section 1, item 2; Section 3, item 3
**Purpose:** Show momentum flux dominates over pressure
**Effort:** ~2 hours

**Implementation:**
```python
# Modify slab/surface.py to return decomposition
def force_incompressible_quadrature(..., return_components=False):
    """
    When return_components=True, return dict:
    {
        'momentum': F_momentum,  # ρ v(v·n) term
        'pressure': F_pressure,  # -(P + ½ρv²)n term
        'total': F_total
    }
    """

# In slab/viz.py
def plot_force_decomposition(diagnostics: Dict, output_path: str):
    """
    Bar chart or stacked time series:
    - |F_momentum|
    - |F_pressure|
    - Total
    Show that momentum >> pressure in incompressible limit
    """
```

**Deliverable:**
- Modified force functions to return components
- Bar chart visualization
- Log F_mom and F_press per body per timestep (optional)

**Impact:** Shows the physics mechanism explicitly

---

#### 4. No-G Hygiene + README Text ⭐⭐⭐⭐

**From:** Section 1, item 5; Section 4
**Purpose:** Professional communication, emphasize no G assumption
**Effort:** ~1 hour

**Implementation:**
- Create example configs using μ_a notation
- Add drop-in text from EXTRA.md Section 4 to README
- Document that K is derived diagnostic only

**Deliverables:**
- `examples/mercury_orbit_no_g.yaml` - Uses μ_sun instead of M
- README section: "Method" (direct quote from EXTRA.md)
- README section: "Key Results" (direct quote from EXTRA.md)
- README section: "Controls" (direct quote from EXTRA.md)

**Drop-in text from EXTRA.md Section 4:**

> **Method:** "Forces are computed exclusively from superfluid surface integrals of momentum/pressure flux around each body. The velocity field is the Green's-function solution for point intakes; no pairwise law is assumed or used."
>
> **Key result:** "Inverse-square attraction emerges with the superfluid coefficient fixed by (ρ₀,β₀). Orbits reproduce all 1PN diagnostics without using G or kg masses."
>
> **Controls:** "Analytic vs quadrature audits, R-invariance, Δt studies, and mass-intake toggles."
>
> **No-G compliance:** "Inputs are μ_a (from orbits), ρ₀, β₀, c_s. G appears only in external comparison plots."

**Impact:** Clear, professional communication; preempts critiques

---

#### 5. Mass-Intake Diagnostic ⭐⭐⭐

**From:** Section 2, item 3
**Purpose:** Show Ṁ is negligible, validates quasi-static approximation
**Effort:** ~1 hour

**Implementation:**
```python
# In slab/diagnostics.py (or add to existing diagnostics)
def mass_drift_diagnostic(trajectory: Dict) -> Dict:
    """
    For each body:
    - Track M(t)
    - Compute ΔM = M(t) - M(0)
    - Report |ΔM|/M per orbit or per century
    - Should be ≪ 1e-12 for Solar System
    """
```

**Deliverable:**
- Function in diagnostics module
- Report in diagnostics.json
- Quick check: enable flux-based Ṁ, verify negligible

**Impact:** Validates approximation that M is constant

---

#### 6. Example Configurations ⭐⭐⭐

**From:** Section 7, item 1
**Purpose:** Reproducibility, standard test cases
**Effort:** ~1 hour

**Implementation:**
Create well-documented example configs:

1. **`examples/sun_mercury.yaml`** - Current Mercury orbit
2. **`examples/sun_earth.yaml`** - Earth analog (circular, e ≈ 0)
3. **`examples/equal_mass.yaml`** - Symmetric case (M₁ = M₂)
4. **`examples/mercury_eccentric.yaml`** - Higher eccentricity (e = 0.2)
5. **`examples/quick_test.yaml`** - Reduced steps/quadrature for CI

**Deliverable:**
- 5 example YAML files with detailed comments
- README section listing examples and their purposes

**Impact:** Easy for others to reproduce, test, and extend

---

### Priority 2: Core Validation (~10 hours total)

Critical for publication but more time-intensive:

#### 7. 1PN Precession Validation ⭐⭐⭐⭐⭐

**From:** Section 1, item 4
**Purpose:** Validate 1PN agreement with NO free parameters
**Effort:** ~3 hours

**Implementation:**
```python
# scripts/validate_1pn_precession.py
"""
Run Mercury orbit with multiple eccentricities:
- e = 0.0, 0.1, 0.2, 0.3, 0.5, 0.7
For each:
- Measure perihelion precession Δω per orbit
- Compute a(1-e²)⁻¹
- Compare to GR-1PN prediction (no free params)
Plot: Δω vs a(1-e²)⁻¹ with GR overlay
"""
```

**Deliverable:**
- Script to run eccentricity sweep
- Plot: Δω vs a(1-e²)⁻¹ with GR-1PN curve
- Report: Agreement within X%

**Impact:** **CRITICAL** - Proves 1PN analogy is exact, not tuned

**Requirements:**
- Compressible forces working correctly
- GR-1PN module (`slab/gr1pn.py`) integrated

---

#### 8. Orbit Overlay Plots ⭐⭐⭐⭐

**From:** Section 3, item 2
**Purpose:** Direct visual comparison with GR
**Effort:** ~4 hours

**Implementation:**
```python
# In slab/viz.py
def plot_orbit_comparison(traj_slab: Dict, traj_gr: Dict,
                          output_path: str):
    """
    Two panels:
    - Top: Overlaid trajectories (slab=blue, GR=red, dashed)
    - Bottom: Difference |r_slab - r_gr| over time

    Also plot Δω (precession angle difference) per orbit
    """
```

**Deliverable:**
- Comparison plot function
- Example: Mercury orbit (slab vs GR-1PN)
- Quantify: max separation, precession agreement

**Impact:** Compelling visual demonstration

**Requirements:**
- GR-1PN module integration
- Same initial conditions for both

---

#### 9. Robustness Tests ⭐⭐⭐

**From:** Section 2, items 1-2
**Purpose:** Standard numerical validation
**Effort:** ~3 hours

**Implementation:**
```python
# scripts/robustness_tests.py
"""
Parameter sweeps:

1. Timestep invariance:
   - dt ∈ [0.001, 0.002, 0.005, 0.01] × T_orbit
   - Measure: Energy drift, precession rate
   - Show: Convergence as dt → 0

2. Sphere radius invariance:
   - R ∈ [0.5, 1.0, 2.0] × R_nominal
   - Measure: Force magnitude, precession
   - Show: Results independent of R (after renormalization)

3. Quadrature convergence:
   - n_points ∈ [64, 128, 256, 512, 1024]
   - Measure: Force, audit error
   - Show: Convergence for n > 256
"""
```

**Deliverable:**
- Sweep script with parameter grids
- Convergence plots (3 panels)
- Report: All results stable within tolerance

**Impact:** Standard validation, expected by reviewers

---

### Priority 3: Polish (Optional, 2-4 hours)

#### 10. CI/Testing Infrastructure

**From:** Section 7, items 3-4
**Purpose:** Continuous validation
**Effort:** ~2 hours

**Implementation:**
- Add `results/expected/` with reference CSVs
- GitHub Actions workflow:
  - Run `examples/quick_test.yaml`
  - Compare output to expected results
  - Fail if metrics exceed tolerances

---

## MEDIUM IMPACT Recommendations

### Communication & Documentation

#### 11. Canned Responses to Critiques ⭐⭐⭐

**From:** Section 5
**Purpose:** Pre-empt common objections
**Effort:** ~30 minutes

Add to README or separate `RESPONSES.md`:

**"You hard-coded Newtonian forces."**
→ Point to emergent inverse-square scan. Show surface-integral code path. Include momentum-flux equation in README.

**"Hidden use of G."**
→ Point to `examples/mercury_orbit_no_g.yaml` using μ_a. Ship configs to replicate.

**"Preferred frame?"**
→ Note: Orbital predictions depend only on relative configurations. Global background drift cancels for Solar System tests.

**"Energy/momentum conservation?"**
→ Point to diagnostics showing |ΔE/E| < 2.4×10⁻⁷ and |Δp| < 10⁻³⁰. Provide flux budgets.

---

## LOW IMPACT / FUTURE Extensions

#### 12. Binary Pulsar Smoke Test

**From:** Section 6, item 2
**Purpose:** Check radiation scaling
**Effort:** ~4-6 hours

Not critical for initial publication. Could be follow-up work.

---

#### 13. Lab Analog Appendix

**From:** Section 6, item 3
**Purpose:** Experimental validation path
**Effort:** Depends on collaboration

Interesting but outside scope of current simulator work.

---

## Critical Issue: Audit Warnings

### Current Status

When running `python -m slab.run examples/mercury_orbit.yaml`, we see:
```
RuntimeWarning: Audit FAILED at step 500 for body 'Sun':
    relative error 3.635e+07 > tolerance 1.000e-03

RuntimeWarning: Audit FAILED at step 500 for body 'Mercury':
    relative error 2.329e+00 > tolerance 1.000e-03
```

### Analysis

**For Sun:** Error of 3.6×10⁷ is enormous
- Suggests analytic and quadrature are giving wildly different answers
- Possibly one is near-zero and the other isn't

**For Mercury:** Error of 2.3 is large but more reasonable
- Both methods likely non-zero but disagree by ~factor of 2

### Possible Causes

1. **Analytic compressible not fully implemented**
   - Currently delegates to quadrature (see `slab/surface.py:force_compressible_analytic`)
   - Audit might be comparing quadrature to itself incorrectly

2. **Audit comparison logic error**
   - May be comparing wrong quantities
   - Check `slab/dynamics.py` audit code

3. **Renormalization issue**
   - Rules A and B (v_ext vs v_total in different terms)
   - Might not be applied consistently

4. **Self-field subtraction**
   - Sun is nearly stationary → v_ext ≈ 0
   - Compressible correction might blow up or vanish incorrectly

### Action Required

**Before implementing audit parity plots, must:**
1. Debug audit comparison logic
2. Verify compressible analytic implementation
3. Fix discrepancies
4. Get audit passing (error < 1e-3)

**Investigation approach:**
1. Check what `force_compressible_analytic()` actually does
2. Look at audit logging in `assemble_forces()`
3. Print both force values when audit fails
4. Understand why Sun error >> Mercury error

---

## Recommended Implementation Order

### Week 1: Quick Wins + Debug

**Days 1-2: Debug compressible audit issue**
- Investigate audit warnings
- Fix compressible force implementation
- Get audit passing

**Days 3-4: Quick wins (Phase 1)**
- Emergent inverse-square scan
- Force decomposition bars
- No-G hygiene + README text
- Mass-intake diagnostic
- Example configurations

**Day 5: Audit visualization**
- Audit parity plot (now that it's working)

### Week 2: Core Validation

**Days 1-2: 1PN precession validation**
- Eccentricity sweep
- GR comparison

**Days 3-4: Orbit overlays**
- Integrate gr1pn module
- Comparison plots

**Day 5: Robustness tests**
- Parameter sweeps
- Convergence studies

### Week 3: Polish

- CI infrastructure
- Documentation refinement
- Response preparation

---

## Summary Table: Implementation Priorities

| Item | Section | Effort | Impact | Status | Priority |
|------|---------|--------|--------|--------|----------|
| **Audit debugging** | - | 4h | ⭐⭐⭐⭐⭐ | ❌ Broken | **P0** (blocker) |
| Inverse-square scan | 1.1 | 2h | ⭐⭐⭐⭐⭐ | ❌ Not done | **P1** |
| Audit parity plot | 3.4 | 1h | ⭐⭐⭐⭐ | ⏸ Blocked | **P1** |
| Force decomp bars | 1.2, 3.3 | 2h | ⭐⭐⭐ | ❌ Not done | **P1** |
| No-G hygiene + README | 1.5, 4 | 1h | ⭐⭐⭐⭐ | 🔧 Partial | **P1** |
| Mass-intake diagnostic | 2.3 | 1h | ⭐⭐⭐ | ❌ Not done | **P1** |
| Example configs | 7.1 | 1h | ⭐⭐⭐ | 🔧 Partial | **P1** |
| 1PN precession valid. | 1.4 | 3h | ⭐⭐⭐⭐⭐ | ⏸ Blocked | **P2** |
| Orbit overlays | 3.2 | 4h | ⭐⭐⭐⭐ | ❌ Not done | **P2** |
| Robustness tests | 2.1-2 | 3h | ⭐⭐⭐ | ❌ Not done | **P2** |
| Canned responses | 5 | 0.5h | ⭐⭐⭐ | ❌ Not done | **P3** |
| CI infrastructure | 7.3-4 | 2h | ⭐⭐ | ❌ Not done | **P3** |

**Legend:**
- ✅ Done
- 🔧 Partial
- ❌ Not done
- ⏸ Blocked (waiting on other items)

---

## Critical Path

```
[Audit Debug] → [Audit Plot, 1PN Validation]
      ↓
[Quick Wins P1] → [Publication-ready figures]
      ↓
[Core Valid P2] → [Submission-ready paper]
      ↓
[Polish P3] → [Release & reproducibility]
```

**Blocker:** Audit issue must be fixed first to unblock 1PN validation and audit visualization.

---

## Key Takeaways

1. **Most critical:** Debug audit warnings (Priority 0)
2. **Quick wins available:** ~8 hours of work dramatically strengthens case (Priority 1)
3. **Core validation:** 1PN precession comparison is publication-critical (Priority 2)
4. **Drop-in text ready:** Section 4 of EXTRA.md has professional communication language
5. **Some items already done:** Energy conservation, basic audit system, compressible forces exist

**Next Action:** Launch agents to investigate audit warnings and fix compressible force issues.

---

**Document Version:** 1.0
**Last Updated:** 2025-10-31
**Next Steps:** Debug audit, then implement Phase 1 quick wins
