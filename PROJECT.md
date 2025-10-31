# Superfluid Orbit Simulator - Project Tracking

## Project Goal
Implement a software simulator that uses superfluid hydrodynamics to model 1PN general relativity effects, based on the "Superfluid Slab Gravity Without G" paper.

## Core Principles
- **No gravitational constant G**: All forces emerge from superfluid momentum flux
- **Bodies as fluid intakes**: Mass objects are "mouths" that remove fluid at rate Q_a
- **Forces from surface integrals**: Momentum and pressure flux around control surfaces
- **Emergent 1/r² law**: Newtonian gravity appears without being assumed
- **1PN-like corrections**: From finite sound speed, healing length, mouth size

## Physics Summary

### Key Equations
- **Mass-intake relation**: M_a = β Q_a
- **Orbital constant**: K = ρ₀/(4πβ²) (replaces G)
- **Velocity field**: v(x) = Σ_b (Q_b/4πρ₀) r_b/r_b³
- **Incompressible force**: F_a = (4/3)(Q_a/4π) v_ext(x_a)
- **Control surface lemma**: F_a = ρ₀ Q_a v_ext(r_a)

### Parameters
- ρ₀: ambient density
- c_s: sound speed (controls retardation corrections)
- β₀: mass-intake factor (Q/M = 1/β₀)
- ξ: healing length (dispersion)
- a₀: mouth size

## Implementation Status

### Phase 1: Core Infrastructure ⏳
- [ ] Directory structure (slab/ modules)
- [ ] geometry.py - Fibonacci sphere, surface integrals
- [ ] medium.py - Medium dataclass
- [ ] bodies.py - Body dataclass with M↔Q sync
- [ ] Configuration system (YAML)

### Phase 2: Physics Engine ⏳
- [ ] field.py - v_total(), v_ext(), v_self()
- [ ] surface.py - Incompressible forces (analytic + quadrature)
- [ ] surface.py - Compressible forces with renormalization
- [ ] surface.py - Intake flux calculation
- [ ] dynamics.py - Velocity-Verlet integrator

### Phase 3: Validation & Comparison ⏳
- [ ] diagnostics.py - Orbital elements, energy checks
- [ ] gr1pn.py - EIH 1PN comparator (separate)
- [ ] Unit tests for physics identities

### Phase 4: Runner & Output ⏳
- [ ] run.py - CLI entry point
- [ ] io_cfg.py - Config validation
- [ ] Table output for debugging
- [ ] CSV export
- [ ] (Later) Visualization plots

## Test Acceptance Criteria

### Test 1: Emergent 1/r² law
- Place two equal sinks, measure F via quadrature
- **Pass**: |C/C_theory - 1| < 0.005
- C_theory = (4/3)|Q₁Q₂|/(4πρ₀)

### Test 2: Orbit Stability
- Two-body, low eccentricity
- **Pass**: |Δa|/a < 10⁻⁵ over 50 orbits

### Test 3: Compressible Scaling
- Extract perihelion precession
- **Pass**: (a) finite, (b) scales ∝ c_s⁻²

### Test 4: Quadrature Audit
- Compare analytic vs quadrature forces
- **Pass**: relative diff < 10⁻³ (incomp), < 5×10⁻³ (comp)

### Test 5: GR-1PN Comparator
- Independent module reproduces Mercury precession
- Within few % of standard result

## Design Decisions

### Language & Tools
- **Python** for implementation
- **SymPy** for symbolic verification
- **NumPy** for numerics
- **PyYAML** for configuration
- **pytest** for testing

### Units
- Convenient code units (AU, yr, code masses)
- K reported as diagnostic only

### Numerical Strategy
- **Fast path**: Analytic force equation (5) every step
- **Audit path**: Quadrature every N_audit steps
- **Renormalization**: Use v_ext only for ρ* and P* in compressible mode
- **Integration**: Velocity-Verlet (leapfrog)

## Module Structure

```
slab/
  __init__.py
  geometry.py          # Fibonacci sphere, sphere identities
  medium.py            # Medium(rho0, cs, beta0, gamma_beta)
  bodies.py            # Body dataclass; Q=M/beta; sync M↔Q
  field.py             # v_total(), v_ext(), v_self()
  surface.py           # force_incompressible(), force_compressible_renorm(), intake_flux()
  dynamics.py          # step_verlet(), assemble_forces()
  diagnostics.py       # osculating_elements(), peri_finder(), energy checks
  gr1pn.py             # eih_1pn_accel(); separate comparator
  io_cfg.py            # load/validate config; unit checks
  run.py               # CLI entry
tests/
  test_geometry.py
  test_surface.py
  test_forces.py
  test_integration.py
```

## Key Implementation Notes

### Force Calculation
1. **Incompressible** (baseline):
   - Use analytic F_a = (4/3)(Q_a/4π) v_ext
   - Quadrature audit with eq. (2) periodically

2. **Compressible** (small-Mach):
   - **Rule A**: Full v in momentum term ρ v(v·n)
   - **Rule B**: Use v_ext only in ρ*, P* (near-field renormalization)
   - Prevents singular self-field blowup

### Surface Integrals
- Fibonacci sphere sampling (N=256-1024 points)
- Cache normals per body
- Key identity: ∫(n·A)n dA = (4πR²/3)A for constant A

### Performance
- Vectorize pairwise v_ext sums
- Skip heavy quadrature most steps (analytic only)
- Audit quadrature every ~100-1000 steps

## Session Progress

### Session 1 (Current)
- [x] Project planning and tracking file created
- [ ] Directory structure setup
- [ ] Core modules implementation (via agents)
- [ ] Basic 2-body test

## Questions & Decisions Log

### Q1: Sign convention for forces
- **Decision**: Use theoretical exact values, prioritize accuracy
- Factor 4/3 is validated to 4×10⁻⁴, use exact theory

### Q2: Output format
- **Decision**: Start with table output for debugging
- Add CSV export
- Defer plots to later phase

## References
- Main paper: `1pn_no_g.tex`
- Implementation plan: `plan_no_pde.md`
- Working directory: `/var/projects/papers/1pn_no_g/`

## Next Steps
1. Create slab/ directory structure
2. Implement geometry module (fibonacci sphere)
3. Implement medium and bodies modules
4. Implement field calculations
5. Implement surface force calculations
6. Build up to working 2-body orbit

---
*Last updated: 2025-10-31*
