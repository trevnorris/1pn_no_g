# Visualization Implementation Plan
## Superfluid Orbit Simulator

**Date:** 2025-10-31
**Status:** Ready for implementation
**Dependencies:** matplotlib >= 3.5.0, pandas >= 1.3.0 (already in pyproject.toml)

---

## Executive Summary

The superfluid orbit simulator has **all the infrastructure needed** for visualization but **no plotting code exists yet**. The data output (CSV trajectories, JSON diagnostics) and configuration hooks are in place. This document outlines a prioritized plan to implement essential physics validation plots.

**Quick Start:** Implement Priority 1 plots (~9 hours) for publication-ready validation.

---

## Current State

### ✅ What Exists

1. **Data Infrastructure:**
   - CSV output: `output/trajectory.csv` (time, position, velocity, mass, intake)
   - JSON diagnostics: `output/diagnostics.json` (energy, forces, orbital elements, precession)
   - Comprehensive diagnostics module: `slab/diagnostics.py`

2. **Dependencies Already Defined:**
   ```toml
   [project.optional-dependencies]
   viz = [
       "matplotlib>=3.5.0",
       "pandas>=1.3.0",
   ]
   ```

3. **Configuration Hooks Ready:**
   The YAML config already expects plot types:
   ```yaml
   outputs:
     plots:
       - orbit
       - precession_vs_time
       - force_decomp
       - energy
       - elements
   ```

### ❌ What's Missing

- No `slab/viz.py` module
- No matplotlib/plotting code anywhere
- CLI mentions plots in config but doesn't generate them
- No velocity field visualization
- No 3D trajectory plots
- No animation capabilities

---

## Physics Context: What Needs Visualization

This simulator derives **Newtonian gravity from superfluid hydrodynamics**:
- Bodies are **sinks** (fluid intakes) with Q > 0 meaning inward flow
- Forces emerge from **momentum flux** through control surfaces
- **No gravitational constant G** - replaced by K = ρ₀/(4πβ²)
- 1PN-like effects from finite sound speed, dispersion, etc.

**Key claims to visualize:**
1. Stable orbits emerge without assuming gravity
2. 1/r² force law emerges naturally from fluid dynamics
3. Energy conservation validates symplectic integrator
4. Perihelion precession from finite sound speed (not GR curvature)
5. Velocity field shows radial inflow patterns around sinks

---

## Priority 1: Essential Physics Validation Plots ⭐

**Time Estimate:** ~9 hours
**Impact:** Critical for publication and scientific validation

### 1.1 Orbital Trajectory Plot (`plot_orbit`)

**Purpose:** Verify stable orbits, demonstrate emergent Newtonian dynamics

**Format:** 2D plot (x-y plane projection)
- Multi-body trajectories with different colors/styles
- Show control surface radii as circles at initial positions
- Mark periapsis/apoapsis passages
- Equal aspect ratio to avoid distortion

**Physics to demonstrate:**
- Stable elliptical orbits without assuming gravity
- Closed vs. precessing orbits (compressible vs. incompressible)

**Implementation:** Straightforward matplotlib scatter/line plots
**Estimated time:** 2 hours

**Example code structure:**
```python
def plot_orbit(trajectory: Dict, config: Dict, output_path: str) -> None:
    """Plot 2D orbital trajectories in x-y plane."""
    fig, ax = plt.subplots(figsize=(10, 10))

    for i, body in enumerate(bodies):
        x = positions[:, i, 0]
        y = positions[:, i, 1]
        ax.plot(x, y, label=body.name, linewidth=2)
        ax.scatter(x[0], y[0], s=100, marker='o')  # Start position

        # Draw control surface
        circle = Circle((x[0], y[0]), body.R, fill=False,
                       linestyle='--', alpha=0.3)
        ax.add_patch(circle)

    ax.set_xlabel('x [AU]', fontsize=14)
    ax.set_ylabel('y [AU]', fontsize=14)
    ax.set_aspect('equal')
    ax.legend()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
```

---

### 1.2 Energy Conservation Plot (`plot_energy`)

**Purpose:** Validate symplectic integrator quality

**Format:** Time series with 2 subplots
- Top: Total energy E(t), kinetic T(t), potential U(t)
- Bottom: Relative energy drift |ΔE/E₀|

**Acceptance criteria:**
- |ΔE/E| < 10⁻⁵ for incompressible (from plan_no_pde.md § 10.2)
- Should show bounded oscillation, not secular drift

**Implementation:** Direct plot from diagnostics JSON
**Estimated time:** 1 hour

**Example code:**
```python
def plot_energy(diagnostics: Dict, output_path: str) -> None:
    """Plot energy conservation over time."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    times = diagnostics['times']
    E_total = diagnostics['total_energy']
    E_kin = diagnostics['kinetic_energy']
    E_pot = diagnostics['potential_energy']

    # Top: Energy components
    ax1.plot(times, E_total, label='Total E', linewidth=2)
    ax1.plot(times, E_kin, label='Kinetic T', alpha=0.7)
    ax1.plot(times, E_pot, label='Potential U', alpha=0.7)
    ax1.set_ylabel('Energy')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Bottom: Relative drift
    E0 = E_total[0]
    drift = np.abs((E_total - E0) / E0)
    ax2.semilogy(times, drift, linewidth=2)
    ax2.axhline(1e-5, color='r', linestyle='--',
                label='Target: 1e-5')
    ax2.set_xlabel('Time [yr]')
    ax2.set_ylabel('|ΔE/E₀|')
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
```

---

### 1.3 Force Law Validation (`plot_force_law`)

**Purpose:** Demonstrate emergent 1/r² scaling

**Format:** Log-log plot
- X-axis: log(r) - separation distance
- Y-axis: log(F) - force magnitude
- Overlay: Theoretical 1/r² line with slope = -2

**Physics to demonstrate:**
- Slope ≈ -2.0 (validate 1/r² scaling)
- Coefficient matches K·M₁·M₂ (no G needed!)

**Implementation:** Extract force vs. distance from trajectory, fit power law
**Estimated time:** 3 hours

**Example code:**
```python
def plot_force_law(trajectory: Dict, medium, output_path: str) -> None:
    """Validate 1/r² force scaling."""
    # Extract separations and compute forces
    r_values = []
    F_values = []

    for snapshot in trajectory:
        r = np.linalg.norm(snapshot['x'][1] - snapshot['x'][0])
        # Compute theoretical force: K M1 M2 / r²
        F = medium.K * M1 * M2 / (r * r)
        r_values.append(r)
        F_values.append(F)

    fig, ax = plt.subplots(figsize=(10, 8))

    # Log-log plot
    ax.loglog(r_values, F_values, 'o', alpha=0.5,
              label='Measured')

    # Theoretical 1/r² line
    r_theory = np.logspace(np.log10(min(r_values)),
                          np.log10(max(r_values)), 100)
    F_theory = medium.K * M1 * M2 / (r_theory ** 2)
    ax.loglog(r_theory, F_theory, 'r--', linewidth=2,
              label='Theory: K M₁M₂/r² (slope=-2)')

    # Fit and display slope
    from scipy.stats import linregress
    slope, intercept, r_val, p_val, std_err = linregress(
        np.log(r_values), np.log(F_values))
    ax.text(0.05, 0.95, f'Fitted slope: {slope:.3f}',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat'))

    ax.set_xlabel('Separation r [AU]', fontsize=14)
    ax.set_ylabel('Force magnitude F', fontsize=14)
    ax.set_title('Emergent 1/r² Force Law (No G!)', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, which='both')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
```

---

### 1.4 Perihelion Precession Plot (`plot_precession`)

**Purpose:** Show 1PN-like effects from finite sound speed

**Format:** Dual plot
- Top: Precession angle ω(t) vs. time (linear trend for secular precession)
- Bottom: Precession rate Δω vs. 1/c_s² (demonstrate scaling)

**Physics to demonstrate:**
- Secular precession for compressible mode
- Scaling ∝ c_s⁻² validates 1PN analogy

**Implementation:** Use `diagnostics.find_periapsis_passages()` and `compute_precession()`
**Estimated time:** 3 hours

**Example code:**
```python
def plot_precession(diagnostics: Dict, output_path: str) -> None:
    """Plot perihelion precession vs time."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Extract periapsis data from diagnostics
    peri_times = diagnostics.get('periapsis_times', [])
    peri_angles = diagnostics.get('periapsis_angles', [])

    # Top: Precession angle vs time
    ax1.plot(peri_times, np.rad2deg(peri_angles),
             'o-', markersize=6, linewidth=2)

    # Fit linear trend
    if len(peri_times) > 2:
        slope, intercept = np.polyfit(peri_times, peri_angles, 1)
        ax1.plot(peri_times,
                np.rad2deg(slope * np.array(peri_times) + intercept),
                'r--', linewidth=2,
                label=f'Δω/orbit = {np.rad2deg(slope):.2e}°')

    ax1.set_xlabel('Time [yr]', fontsize=14)
    ax1.set_ylabel('Argument of Periapsis ω [deg]', fontsize=14)
    ax1.set_title('Perihelion Precession (from Superfluid, not GR!)',
                  fontsize=16)
    ax1.legend(fontsize=12)
    ax1.grid(alpha=0.3)

    # Bottom: Scaling with sound speed (if multiple runs available)
    # This would require multiple simulation runs with different cs
    # For now, show the single measurement
    ax2.text(0.5, 0.5,
             'Run multiple simulations with varying cs\n'
             'to demonstrate Δω ∝ cs⁻² scaling',
             ha='center', va='center', transform=ax2.transAxes,
             fontsize=12, bbox=dict(boxstyle='round',
                                   facecolor='lightblue'))
    ax2.set_xlabel('1/cs² [yr²/AU²]', fontsize=14)
    ax2.set_ylabel('Precession rate Δω [deg/orbit]', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
```

---

## Priority 2: Field Visualization 🌊

**Time Estimate:** ~6 hours
**Impact:** High value for demonstrating superfluid physics

### 2.1 Velocity Field Quiver Plot (`plot_velocity_field`)

**Purpose:** Show radial inflow patterns from superfluid sinks

**Format:** 2D quiver plot on x-y plane slice
- Background color: velocity magnitude |v(x)|
- Arrows: velocity direction
- Overlays: Body positions as circles with radius R

**Physics to demonstrate:**
- Potential flow: radial inflow toward sinks
- Linear superposition of sink fields
- 1/r² velocity decay

**Implementation:**
```python
def plot_velocity_field(bodies: List, medium, t: float,
                        output_path: str) -> None:
    """2D quiver plot of velocity field."""
    # Sample velocity field on grid
    x_range = np.linspace(-2, 2, 50)
    y_range = np.linspace(-2, 2, 50)
    X, Y = np.meshgrid(x_range, y_range)

    # Compute v(x, y, z=0) using slab.field.v_total()
    from slab.field import v_total

    Vx = np.zeros_like(X)
    Vy = np.zeros_like(Y)
    Vmag = np.zeros_like(X)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            pos = np.array([X[i,j], Y[i,j], 0.0])
            v = v_total(pos, bodies, medium.rho0)
            Vx[i,j] = v[0]
            Vy[i,j] = v[1]
            Vmag[i,j] = np.linalg.norm(v)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Background: velocity magnitude (log scale)
    im = ax.contourf(X, Y, np.log10(Vmag + 1e-20),
                     levels=20, cmap='viridis')
    plt.colorbar(im, ax=ax, label='log₁₀|v|')

    # Arrows: velocity direction
    skip = 3  # Show every 3rd arrow
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
              Vx[::skip, ::skip], Vy[::skip, ::skip],
              color='white', alpha=0.6)

    # Overlay body positions
    for body in bodies:
        circle = Circle((body.x[0], body.x[1]), body.R,
                       fill=True, color='red', alpha=0.7)
        ax.add_patch(circle)
        ax.text(body.x[0], body.x[1] + body.R + 0.1,
                body.name, ha='center', color='white')

    ax.set_xlabel('x [AU]', fontsize=14)
    ax.set_ylabel('y [AU]', fontsize=14)
    ax.set_title(f'Superfluid Velocity Field at t={t:.2f} yr',
                 fontsize=16)
    ax.set_aspect('equal')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
```

**Estimated time:** 4 hours

---

### 2.2 Velocity Field Magnitude Heatmap (`plot_velocity_magnitude`)

**Purpose:** Visualize field structure around bodies

**Format:** 2D heatmap with logarithmic color scale
- log₁₀|v(x)| to show 1/r² structure
- Multiple distance scales visible

**Implementation:** Sample `v_total()` on 2D grid, use `plt.imshow()` with log norm
**Estimated time:** 2 hours

---

## Priority 3: Diagnostics & Comparison 📊

**Time Estimate:** ~9 hours
**Impact:** Medium (useful for debugging and validation)

### 3.1 Force Decomposition (`plot_force_decomp`)

**Purpose:** Break down force contributions

**Format:** Stacked time series
- Incompressible baseline F_inc(t)
- Compressible correction F_comp(t)
- Total F_total(t)
- Show relative magnitudes

**Physics to demonstrate:**
- Compressible correction is small (few %) but measurable
- Correction grows as c_s decreases

**Estimated time:** 2 hours

---

### 3.2 Orbital Elements Evolution (`plot_elements`)

**Purpose:** Track osculating elements over time

**Format:** Multi-panel time series
- Semi-major axis a(t)
- Eccentricity e(t)
- Argument of periapsis ω(t) (shows precession)
- Inclination i(t)

**Implementation:** Call `diagnostics.osculating_elements()` at each saved step
**Estimated time:** 2 hours

---

### 3.3 GR-1PN Comparison (`plot_gr_comparison`)

**Purpose:** Validate against independent GR calculation

**Format:** Side-by-side or overlay plots
- Left: Slab trajectory
- Right: GR-1PN trajectory
- Bottom: Difference |Δr(t)|

**Physics to demonstrate:**
- Slab reproduces GR predictions
- Precession rates match within error bars

**Estimated time:** 5 hours (requires gr1pn module integration)

---

## Priority 4: Advanced & Interactive 🚀

**Time Estimate:** ~18-30 hours
**Impact:** Low (nice to have, not essential for publication)

### 4.1 3D Trajectory Visualization
- Use `mpl_toolkits.mplot3d` or Plotly
- Interactive rotation for multi-body systems
- Time-slider animation
- **Estimated time:** 4 hours

### 4.2 Phase Space Plots
- Position vs. velocity phase diagrams
- Poincaré sections for periodic orbits
- **Estimated time:** 4 hours

### 4.3 Animation
- Bodies moving along trajectories
- Velocity field updating in real-time
- Export to MP4 via matplotlib.animation
- **Estimated time:** 6 hours

### 4.4 Interactive Dashboards (Optional)
- Plotly Dash or Bokeh for web interface
- Sliders for cs, β₀ parameters
- Real-time simulation + visualization
- **Estimated time:** 8-12 hours

---

## Implementation Roadmap

### Phase 1: Core Module Setup (Week 1 - Day 1)

**Create `/var/projects/1pn_no_g/slab/viz.py`**

```python
"""Visualization module for superfluid orbit simulator.

This module provides plotting functions to visualize orbital dynamics,
energy conservation, force laws, and velocity fields emerging from
superfluid hydrodynamics.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle
from typing import Dict, List, Optional
import pandas as pd


def plot_orbit(trajectory: Dict, config: Dict, output_path: str) -> None:
    """Plot 2D orbital trajectories in x-y plane."""
    pass  # Implement


def plot_energy(diagnostics: Dict, output_path: str) -> None:
    """Plot energy conservation over time."""
    pass  # Implement


def plot_force_law(trajectory: Dict, medium, output_path: str) -> None:
    """Validate 1/r² force scaling."""
    pass  # Implement


def plot_precession(diagnostics: Dict, output_path: str) -> None:
    """Plot perihelion precession vs time."""
    pass  # Implement


def plot_velocity_field(bodies: List, medium, t: float,
                        output_path: str) -> None:
    """2D quiver plot of velocity field."""
    pass  # Implement


def generate_plots(config: Dict, results: Dict,
                   output_dir: str) -> None:
    """Generate all requested plots from config.

    Parameters
    ----------
    config : dict
        Configuration dictionary with 'outputs']['plots'] list
    results : dict
        Dictionary containing 'trajectory' and 'diagnostics' data
    output_dir : str
        Directory to save plot files
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    plot_list = config.get('outputs', {}).get('plots', [])

    for plot_name in plot_list:
        if plot_name == 'orbit':
            plot_orbit(results['trajectory'], config,
                      f'{output_dir}/orbit.png')
        elif plot_name == 'energy':
            plot_energy(results['diagnostics'],
                       f'{output_dir}/energy.png')
        elif plot_name == 'force_law':
            plot_force_law(results['trajectory'], config['medium'],
                          f'{output_dir}/force_law.png')
        elif plot_name == 'precession_vs_time':
            plot_precession(results['diagnostics'],
                           f'{output_dir}/precession.png')
        elif plot_name == 'velocity_field':
            plot_velocity_field(results['bodies'], config['medium'],
                               results['t_final'],
                               f'{output_dir}/velocity_field.png')
        else:
            print(f"Warning: Unknown plot type '{plot_name}'")
```

---

### Phase 2: CLI Integration (Week 1 - Day 1)

**Modify `/var/projects/1pn_no_g/slab/run.py`:**

```python
# Add import at top
from slab import viz

# In main() or save_outputs() function, after saving CSV/JSON:
def save_outputs(...):
    # ... existing CSV/JSON save code ...

    # Generate plots if requested
    if config.get('outputs', {}).get('plots'):
        print("\n" + "="*70)
        print("GENERATING PLOTS")
        print("="*70)

        results = {
            'trajectory': trajectory_data,
            'diagnostics': diagnostics_data,
            'bodies': bodies,
            't_final': t_final,
        }

        viz.generate_plots(config, results, output_dir)
        print(f"\n✓ Plots saved to {output_dir}/")
```

---

### Phase 3: Implement Priority 1 Plots (Week 1 - Days 2-3)

**Timeline:**
- Day 2 Morning: `plot_orbit()` (2 hours)
- Day 2 Afternoon: `plot_energy()` (1 hour)
- Day 3 Morning: `plot_force_law()` (3 hours)
- Day 3 Afternoon: `plot_precession()` (3 hours)

**Testing after each plot:**
```bash
python -m slab.run examples/mercury_orbit.yaml
# Check output/orbit.png, output/energy.png, etc.
```

---

### Phase 4: Implement Priority 2 Plots (Week 2)

**Timeline:**
- Day 1: `plot_velocity_field()` (4 hours)
- Day 2: `plot_velocity_magnitude()` (2 hours)

---

### Phase 5: Polish & Documentation (Week 2 - Days 3-4)

**Features:**
- Consistent style (colors, fonts, LaTeX labels)
- Figure saving (PNG, PDF, SVG)
- Multi-format support via `savefig()` kwargs
- Publication-quality defaults (300 DPI, tight layout)

**Documentation:**
- Add docstrings with Examples sections
- Document expected data formats
- Add visualization section to README.md

**Example README addition:**
```markdown
## Visualization

Generate plots of orbital trajectories, energy conservation, and
velocity fields:

```yaml
# config.yaml
outputs:
  plots:
    - orbit              # 2D trajectories
    - energy             # Energy conservation
    - force_law          # 1/r² validation
    - precession_vs_time # Perihelion precession
    - velocity_field     # Superfluid flow patterns
```

```bash
python -m slab.run config.yaml
ls output/figures/
# orbit.png, energy.png, force_law.png, ...
```
```

---

## Module Structure

```
slab/
├── viz.py              # Main plotting module (NEW)
│   ├── plot_orbit()
│   ├── plot_energy()
│   ├── plot_force_law()
│   ├── plot_precession()
│   ├── plot_velocity_field()
│   ├── plot_velocity_magnitude()
│   ├── plot_force_decomp()
│   ├── plot_elements()
│   ├── plot_gr_comparison()
│   └── generate_plots()  # Main dispatcher
│
└── run.py              # CLI integration (MODIFY)
    └── Call viz.generate_plots() in save_outputs()
```

---

## Output Structure

```
output/
├── trajectory.csv           # Existing
├── diagnostics.json         # Existing
└── figures/                 # NEW
    ├── orbit.png
    ├── orbit.pdf
    ├── energy_conservation.png
    ├── force_law_log.png
    ├── precession_time.png
    ├── velocity_field_t0.png
    └── summary_multi.png    # Multi-panel overview
```

---

## Configuration Options

**Extend YAML config:**

```yaml
outputs:
  save_every: 1000
  write_csv: true

  # NEW visualization options
  plots:
    - orbit
    - energy
    - force_law
    - precession_vs_time
    - velocity_field

  plot_formats: [png, pdf]  # NEW
  plot_dpi: 300              # NEW
  plot_style: 'seaborn'      # NEW (matplotlib style)
  figures_dir: 'output/figures'  # NEW (optional, defaults to output/)
```

---

## Dependencies

### Current (Already Installed)
```toml
[project.optional-dependencies]
viz = [
    "matplotlib>=3.5.0",
    "pandas>=1.3.0",
]
```

### Suggested Additions
```toml
viz = [
    "matplotlib>=3.5.0",
    "pandas>=1.3.0",
    "seaborn>=0.11.0",        # Better default styles (optional)
]

viz-extended = [
    "plotly>=5.0.0",          # Interactive plots (optional)
    "imageio>=2.9.0",         # Animation export (optional)
]
```

**Installation:**
```bash
pip install -e ".[viz]"              # Basic visualization
pip install -e ".[viz-extended]"     # With interactive/animation
```

---

## Testing Strategy

### Unit Tests

**Create `/var/projects/1pn_no_g/tests/test_viz.py`:**

```python
import pytest
import numpy as np
from slab import viz

def test_plot_orbit_creates_file(tmp_path):
    """Test that plot_orbit generates output file."""
    # Create mock trajectory
    trajectory = {
        't': np.linspace(0, 1, 100),
        'x': np.random.randn(100, 2, 3),  # 2 bodies
    }
    config = {
        'bodies': [
            {'name': 'Sun', 'R': 0.001},
            {'name': 'Planet', 'R': 0.0005}
        ]
    }
    output = tmp_path / "orbit.png"

    viz.plot_orbit(trajectory, config, str(output))

    assert output.exists()
    assert output.stat().st_size > 0


def test_plot_energy_handles_empty_data():
    """Test graceful handling of empty diagnostics."""
    diagnostics = {'times': [], 'total_energy': []}

    with pytest.raises(ValueError, match="No data"):
        viz.plot_energy(diagnostics, "output.png")


def test_generate_plots_skips_unknown():
    """Test that unknown plot types are skipped gracefully."""
    config = {
        'outputs': {
            'plots': ['orbit', 'unknown_plot', 'energy']
        }
    }
    # Should warn but not crash
    # Implementation depends on error handling strategy
```

### Integration Tests

```bash
# Run full simulation with plotting
python -m slab.run examples/mercury_orbit.yaml --verbose

# Check outputs exist
ls output/figures/orbit.png
ls output/figures/energy.png
ls output/figures/force_law.png
ls output/figures/precession.png
```

---

## Performance Considerations

1. **Field sampling:**
   - Cache velocity field grid calculations
   - Use reasonable grid resolution (50×50 for quiver, 200×200 for heatmaps)
   - Avoid recomputing for multiple plot types

2. **Large datasets:**
   - For 200,000 step simulations, subsample trajectory for plotting
   - Use `save_every` to control output frequency
   - Consider separate high-res output for critical sections

3. **Multi-format output:**
   - Generate PNG by default (fast, good quality)
   - Add PDF/SVG as options for publications

---

## Summary Table

| Plot Type | Priority | Time | Physics Value | Implementation Difficulty |
|-----------|----------|------|---------------|--------------------------|
| Orbital trajectory | P1 | 2h | ⭐⭐⭐⭐⭐ | Low |
| Energy conservation | P1 | 1h | ⭐⭐⭐⭐⭐ | Low |
| Force law (1/r²) | P1 | 3h | ⭐⭐⭐⭐⭐ | Medium |
| Perihelion precession | P1 | 3h | ⭐⭐⭐⭐⭐ | Medium |
| Velocity field quiver | P2 | 4h | ⭐⭐⭐⭐ | Medium |
| Velocity heatmap | P2 | 2h | ⭐⭐⭐ | Low |
| Force decomposition | P3 | 2h | ⭐⭐⭐ | Low |
| Orbital elements | P3 | 2h | ⭐⭐⭐ | Low |
| GR comparison | P3 | 5h | ⭐⭐⭐⭐ | High |
| 3D trajectories | P4 | 4h | ⭐⭐ | Medium |
| Animation | P4 | 6h | ⭐⭐⭐ | High |
| Interactive dashboard | P4 | 8-12h | ⭐⭐ | High |

**Total for Priority 1 (Essential):** ~9 hours
**Total for Priority 1-2 (Recommended):** ~15 hours
**Total for Priority 1-3 (Complete):** ~24 hours

---

## Next Actions

### Immediate (This Week)

1. ✅ Review this document
2. ⬜ Create `slab/viz.py` with function stubs
3. ⬜ Integrate into CLI (`slab/run.py`)
4. ⬜ Implement `plot_orbit()` and test
5. ⬜ Implement `plot_energy()` and test
6. ⬜ Implement `plot_force_law()` and test
7. ⬜ Implement `plot_precession()` and test
8. ⬜ Update README.md with visualization examples

### Medium Term (Next Week)

9. ⬜ Add velocity field visualization
10. ⬜ Implement diagnostic plots
11. ⬜ Add multi-format export
12. ⬜ Write unit tests

### Long Term (Future)

13. ⬜ Interactive plots with Plotly
14. ⬜ Animation support
15. ⬜ Comparison with observations

---

## Key Insights

1. **Infrastructure is ready** - All data is computed and saved, just needs plotting
2. **Priority 1 plots are essential** - These directly validate the core physics claims
3. **Implementation is straightforward** - No complex algorithms, mostly matplotlib basics
4. **Visual demonstration is critical** - The simulator's value is showing gravity emerges from fluids
5. **9 hours to publication-ready** - Priority 1 suite provides everything needed for papers

**Recommendation:** Implement Priority 1 plots immediately to unlock the full scientific impact of this work.

---

**Document Version:** 1.0
**Last Updated:** 2025-10-31
**Status:** Ready for implementation
