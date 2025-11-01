# Newtonian Gravity Validation Framework

This directory contains a comprehensive validation framework for establishing baseline precision requirements before measuring 1PN (post-Newtonian) corrections.

## Overview

Before we can reliably measure tiny GR effects (~0.1 arcsec/orbit for Mercury), we must first understand what timestep is required to eliminate spurious numerical precession from the Newtonian baseline. This framework:

1. Implements a pure Newtonian reference integrator (F = -GMm/r²)
2. Compares it with the superfluid implementation (F = ρ₀Qv, incompressible)
3. Measures spurious precession vs timestep across 9 resolutions
4. Establishes minimum precision requirements for 1PN studies

## Key Files

### Scripts

- **`validate_newtonian.py`** - Main validation script
  - Runs convergence study comparing Newtonian vs Superfluid
  - Tests multiple timesteps (default: 8 values from 0.002 to 5×10⁻⁶ yr)
  - Outputs formatted tables with precession, energy drift, and timing
  - Usage: `python validate_newtonian.py --quick` for fast test

- **`plot_newtonian_convergence.py`** - Generates publication figures
  - Creates log-log plots showing dt² scaling
  - Visualizes GR signal vs spurious precession
  - Highlights recommended timestep region
  - Outputs PNG and PDF formats

- **`precession_helpers.py`** - Shared utilities
  - Barycentric Mercury configuration
  - Precession analysis from trajectories
  - Orbital element computation

### Documentation

- **`../NEWTONIAN_BASELINE.md`** - Main results document
  - Executive summary of findings
  - Detailed convergence tables
  - Recommendations for 1PN studies
  - Computational cost estimates

## Quick Start

### Run the validation study:
```bash
python scripts/validate_newtonian.py
```

### Quick test (3 timesteps):
```bash
python scripts/validate_newtonian.py --quick
```

### Custom timesteps:
```bash
python scripts/validate_newtonian.py --dt 0.001 0.0001 1e-5 5e-6
```

### Generate plots:
```bash
python scripts/plot_newtonian_convergence.py
```

## Key Results

### Spurious Precession vs Timestep

| dt [yr]   | Steps/orbit | Spurious Prec [arcsec/orbit] | vs GR Signal |
|-----------|-------------|------------------------------|--------------|
| 0.002     | 121         | ~899                        | 9000×        |
| 0.001     | 241         | ~227                        | 2300×        |
| 0.0001    | 2413        | ~2.3                        | 23×          |
| 0.00001   | 24,130      | ~0.023                      | 0.23×        |
| **5×10⁻⁶** | **48,261** | **~0.006**                 | **0.06×**    |

**GR Signal**: 0.0997 arcsec/orbit (Mercury, cs = 63,240 AU/yr)

### Main Findings

1. **Spurious precession scales as dt²** (velocity-Verlet is 2nd order)
2. **Recommended minimum: dt ≤ 5×10⁻⁶ yr** (48,261 steps/orbit)
   - Achieves <0.01 arcsec/orbit (10% of GR signal)
   - Excellent energy conservation (ΔE/E ~ 10⁻¹³)
3. **Superfluid = Newtonian** (agreement to machine precision)
   - Validates incompressible superfluid force law
   - No systematic offsets detected

### Computational Cost

For 1 century of Mercury orbit (414 orbits):
- **Recommended (dt = 5×10⁻⁶ yr)**: 20 million steps
  - Newtonian: ~10 minutes
  - Superfluid: ~50 minutes

## Validation Methodology

### Test Configuration
- **Orbit**: Mercury barycentric (a = 0.387 AU, e = 0.1)
- **Integrator**: Velocity-Verlet (symplectic, 2nd order)
- **Integration time**: 8 orbits (~1.93 yr)
- **Save frequency**: ~2000 snapshots per run

### Methods Compared

1. **Pure Newtonian**
   ```python
   F_i = -Σ_j (G M_i M_j / r_ij²) r̂_ij
   ```
   Reference implementation with textbook force law.

2. **Incompressible Superfluid**
   ```python
   F_i = ρ₀ Q_i v_ext(x_i)   where v_ext = Σ_j Q_j/(4πr_ij²) r̂_ij
   ```
   Momentum flux through control surface, use_compressible=False.

3. **GR 1PN** (optional)
   ```python
   F_i = F_Newton + (1/c²) × (1PN corrections)
   ```
   Einstein-Infeld-Hoffmann equations for comparison.

### Metrics Measured

1. **Perihelion precession** (primary metric)
   - Measured from periapsis-to-periapsis Δω
   - Should be zero for closed Newtonian orbits
   - Non-zero values indicate numerical errors

2. **Energy conservation**
   - ΔE/E should be O(dt²) for symplectic integrator
   - No secular drift expected

3. **Orbital element drift**
   - Semi-major axis: Δa/a
   - Eccentricity: Δe

## Physical Interpretation

### Why does spurious precession occur?

The velocity-Verlet integrator:
- Conserves energy to O(dt²) ✓
- Preserves phase-space volume (symplectic) ✓
- Does **not** exactly conserve angular momentum for eccentric orbits ✗

Small violations in L̂ direction cause the orbital plane to precess, which manifests as perihelion precession in the osculating elements.

### Why dt² scaling?

- Velocity-Verlet is 2nd-order accurate: local error is O(dt³) per step
- Global error over N steps: O(N × dt³) = O(T_total × dt²)
- Precession is a cumulative effect → scales as dt²

Confirmed empirically: reducing dt by factor of 10 reduces spurious precession by factor of ~100.

## Recommendations for 1PN Studies

### Minimum Requirements

✅ **Timestep**: dt ≤ 5×10⁻⁶ yr
- Achieves <0.01 arcsec/orbit spurious precession
- This is <10% of the GR signal (0.1 arcsec/orbit)

✅ **Integration time**: ≥10 orbits
- Needed for statistical significance
- Improves signal-to-noise in precession measurement

✅ **Always run Newtonian baseline first**
- Measure spurious precession in incompressible case
- Verify it's <10% of expected 1PN signal
- If not, refine timestep further

### Best Practices

1. **Use barycentric initial conditions**
   - Prevents CoM drift from corrupting measurements
   - See `create_barycentric_mercury_config()` in `precession_helpers.py`

2. **Check energy conservation**
   - Should have ΔE/E < 10⁻¹⁰ for acceptable precision
   - Larger drift indicates numerical problems

3. **Measure precession multiple ways**
   - Periapsis-to-periapsis: robust for short runs
   - Linear fit to ω(t): better for long integrations
   - Compare both to check consistency

4. **Run convergence tests**
   - Test 2-3 timesteps to verify dt² scaling
   - Extrapolate to dt→0 if needed (Richardson extrapolation)

## Future Improvements

### Adaptive Timestepping
- Use smaller dt near periapsis (where forces are largest)
- Could reduce cost by factor of 2-3
- Requires careful implementation to maintain symplectic property

### Higher-Order Integrators
- 4th-order symplectic methods (Yoshida, Forest-Ruth)
- Would scale as dt⁴ → potentially 100× fewer steps
- Trade-off: more force evaluations per step

### Richardson Extrapolation
- Run at 2-3 timesteps, extrapolate to dt→0
- Can improve accuracy without full refinement
- Requires dt² scaling to hold (verified here)

## References

1. **Hairer, Lubich & Wanner** (2006). *Geometric Numerical Integration*.
   - Standard reference on symplectic integrators

2. **Quinn, Tremaine & Duncan** (1991). AJ 101, 2287.
   - Long-term stability of leapfrog integrator for planetary orbits

3. **Wisdom & Holman** (1991). AJ 102, 1528.
   - Symplectic maps for N-body problem

## Troubleshooting

### "Precession is still too large at finest timestep"
- Try dt = 2×10⁻⁶ yr or 1×10⁻⁶ yr
- Check for bugs in initial conditions (should be barycentric)
- Verify eccentricity (higher e needs smaller dt)

### "Superfluid disagrees with Newtonian"
- Check that use_compressible=False in opts
- Verify K = ρ₀/(4πβ₀²) matches the G value used
- Check force calculation (should use analytic v_ext)

### "Energy drift is too large"
- Reduce timestep
- Check for collision (r_ij → 0)
- Verify masses are positive

### "Integration is too slow"
- Use --quick flag for testing
- Consider coarser timestep for exploratory runs
- Remember: Superfluid is ~5× slower than pure Newtonian

## Contact

For questions about this validation framework, see:
- Main documentation: `../NEWTONIAN_BASELINE.md`
- Implementation: `validate_newtonian.py`
- Original diagnosis: `final_diagnosis.py`

---
*Generated: 2025-11-01*
*Framework: Newtonian Gravity Validation for 1PN Studies*
