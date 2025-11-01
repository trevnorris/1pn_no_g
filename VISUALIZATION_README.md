# Visualization Tools for Slab-GR Comparison

This document describes the visualization tools implemented to compare slab trajectory simulations with GR 1PN predictions.

## Overview

The visualization system consists of two main components:

1. **`slab/viz.py`**: Plotting functions for trajectory comparison
2. **`scripts/compare_with_gr.py`**: Standalone script to run comparisons

## Files Created

### 1. `slab/viz.py`

Visualization module with three main plotting functions:

#### `plot_orbit_comparison(slab_trajectory, gr_trajectory, output_path, body_idx=1, dpi=150)`

Creates a two-panel comparison plot:
- **Top panel**: Overlaid x-y trajectories (slab=blue solid, GR=red dashed)
- **Bottom panel**: Separation |r_slab - r_gr| over time

**Features**:
- Publication-ready formatting
- Clear labels and legends
- Grid lines for readability
- Statistics box showing max/mean separation
- Initial position marker (green circle)
- Central body marker (black star)

#### `plot_precession_comparison(slab_trajectory, gr_trajectory, body_idx, M_central, K_slab, G_gr, output_path, dpi=150)`

Plots perihelion precession comparison:
- Computes ω(t) (argument of periapsis) for both trajectories
- Unwraps 2π discontinuities
- Fits linear trends to extract precession rates
- Shows both data points and fitted lines
- Reports precession per orbit in rad/orbit and arcsec/orbit
- Displays agreement ratio

**Features**:
- Color-coded data (slab=blue, GR=red)
- Linear fits with slope displayed
- Text box with detailed statistics
- Monospace font for alignment

#### `plot_trajectory_3d(trajectory, body_names, output_path, dpi=150)`

Creates 3D visualization of orbital trajectories:
- Color-coded paths for each body
- Initial positions marked with circles
- Final positions marked with squares
- Central body marked with gold star
- Equal aspect ratio on all axes

**Features**:
- Interactive 3D perspective
- Legend with body names
- Proper axis scaling
- Grid lines for depth perception

### 2. `scripts/compare_with_gr.py`

Standalone script that orchestrates the full comparison workflow:

**Command-line arguments**:
- `--config`: Configuration file path (default: `examples/quick_validation.yaml`)
- `--steps`: Override number of steps from config
- `--output-dir`: Output directory (default: `output/comparison`)
- `--dpi`: Plot resolution (default: 150)
- `--body-idx`: Index of orbiting body (default: 1)

**Workflow**:
1. Load and validate configuration
2. Run slab simulation with configured parameters
3. Run GR 1PN simulation with same initial conditions
4. Compute comparison statistics
5. Generate four plots:
   - `orbit_comparison.png`
   - `precession_comparison.png`
   - `trajectory_3d_slab.png`
   - `trajectory_3d_gr.png`
6. Save summary statistics to `summary.txt`

## Usage Examples

### Basic Usage

```bash
# Run comparison with default config
python scripts/compare_with_gr.py

# Specify a different config
python scripts/compare_with_gr.py --config examples/mercury_orbit.yaml

# Quick test with fewer steps
python scripts/compare_with_gr.py --steps 100 --output-dir output/quick_test

# High-resolution plots for publication
python scripts/compare_with_gr.py --config examples/mercury_orbit.yaml --dpi 300
```

### Using Visualization Functions Directly

```python
from slab.dynamics import integrate_orbit
from slab.gr1pn import integrate_gr1pn_orbit
from slab.viz import plot_orbit_comparison, plot_precession_comparison

# After running both simulations...
slab_traj = integrate_orbit(bodies, medium, dt, n_steps, opts)
gr_traj = integrate_gr1pn_orbit(bodies, c_light, G, dt, n_steps)

# Generate plots
plot_orbit_comparison(
    slab_traj, gr_traj,
    output_path="output/orbit.png",
    body_idx=1
)

plot_precession_comparison(
    slab_traj, gr_traj,
    body_idx=1,
    M_central=1.0,
    K_slab=medium.K,
    G_gr=6.67e-11,
    output_path="output/precession.png"
)
```

## Design Decisions

### Color Scheme
- **Slab**: Blue (consistent with "cold" superfluid concept)
- **GR**: Red (standard for comparison/reference)
- **Initial positions**: Green circles
- **Central body**: Gold/yellow star
- **Markers**: Black edges for visibility

### Units
- Positions: AU (Astronomical Units)
- Velocities: AU/yr
- Time: years
- Precession: radians and arcseconds
- All labels clearly specify units

### Plot Quality
- Default DPI: 150 (good for screen viewing)
- Optional DPI: 300 (publication quality)
- Grid lines: alpha=0.3 for visibility without distraction
- Font sizes: 10-14pt for readability
- Tight layout to avoid label clipping

### Data Handling
- Automatic interpolation if time grids differ
- Handles missing or NaN values gracefully
- Unwraps angular data to avoid 2π jumps
- Uses numpy for efficient array operations

## Output Files

### Generated Plots

1. **orbit_comparison.png**: Two-panel trajectory comparison
   - Size: ~100KB at 150 DPI
   - Shows spatial agreement between slab and GR

2. **precession_comparison.png**: Precession rate analysis
   - Size: ~100KB at 150 DPI
   - Shows temporal evolution of ω(t)

3. **trajectory_3d_slab.png**: 3D view of slab trajectory
   - Size: ~280KB at 150 DPI
   - Useful for visualizing out-of-plane motion

4. **trajectory_3d_gr.png**: 3D view of GR trajectory
   - Size: ~280KB at 150 DPI
   - Reference for comparison

### Summary Statistics

The `summary.txt` file contains:
- Orbital parameters (a, e, period, number of orbits)
- Position differences (max, mean, final separation)
- Precession rates (both trajectories, in rad/orbit and arcsec/orbit)
- Agreement ratio and interpretation

## Testing

### Test with Quick Validation Config

```bash
python scripts/compare_with_gr.py \
    --config examples/quick_validation.yaml \
    --output-dir output/test
```

**Expected results**:
- All plots generated successfully
- Energy conservation: ΔE/E < 10⁻⁵
- Audit tests pass: relative error < 10⁻³
- Separation remains small (< 10⁻¹⁰ AU for identical ICs)

### Interpreting Results

The `quick_validation.yaml` config uses arbitrary code units where:
- K = 7.96e-22 (orbital constant)
- Orbital period ≈ 5.4e10 years (in code units)

**Important**: These are not physical Mercury orbits! The config uses:
- Extremely small K value
- Very long time scales
- Short simulation time (1000 steps = 2 years in code units)

This means precession measurements from short runs may not be meaningful. For meaningful precession comparison, use:
1. Longer simulations (>100 orbits)
2. More realistic parameters (see `examples/mercury_orbit.yaml`)
3. Or accept qualitative comparison only

### Validation Checklist

✓ Visualization module created: `slab/viz.py`
✓ Three plotting functions implemented:
  - `plot_orbit_comparison`
  - `plot_precession_comparison`
  - `plot_trajectory_3d`
✓ Comparison script created: `scripts/compare_with_gr.py`
✓ Command-line interface with arguments
✓ Integration with existing modules:
  - Uses `slab.dynamics.integrate_orbit`
  - Uses `slab.gr1pn.integrate_gr1pn_orbit`
  - Uses `slab.gr1pn.compute_orbit_elements`
✓ Error handling for missing matplotlib
✓ Publication-ready plots (grid, labels, legends)
✓ High-DPI output support
✓ Summary statistics generation
✓ Tested with `quick_validation.yaml`

## Dependencies

Required:
- `numpy` (already in dependencies)
- `matplotlib>=3.5.0` (in optional viz dependencies)

Install visualization dependencies:
```bash
pip install matplotlib
# or
pip install -e ".[viz]"
```

## Integration with Existing Code

The visualization tools integrate seamlessly with:

1. **Configuration system** (`slab.io_cfg`):
   - Loads config from YAML
   - Validates parameters
   - Extracts GR comparison settings

2. **Dynamics module** (`slab.dynamics`):
   - Uses `integrate_orbit()` for slab simulation
   - Returns trajectory in compatible format

3. **GR 1PN module** (`slab.gr1pn`):
   - Uses `integrate_gr1pn_orbit()` for GR reference
   - Uses `compute_orbit_elements()` for precession analysis
   - Returns trajectory in same format as slab

4. **Body and Medium classes**:
   - Works with existing `Body` and `Medium` dataclasses
   - No modifications needed

## Future Enhancements

Potential improvements:
1. Animated trajectories (gif/mp4 output)
2. Interactive plots (plotly/bokeh)
3. Comparison of multiple cs values
4. Energy conservation plots
5. Phase space plots (x-v diagrams)
6. Automated report generation (LaTeX/PDF)
7. Statistical analysis (chi-squared tests)
8. Batch comparison across parameter sweeps

## Troubleshooting

### "matplotlib not available" error
```bash
pip install matplotlib
```

### Plots are blurry
Increase DPI:
```bash
python scripts/compare_with_gr.py --dpi 300
```

### Precession agreement is poor
- Check simulation time: need >100 orbits for good statistics
- Verify initial conditions match between slab and GR
- Check that cs is appropriately chosen
- Review orbital parameters in summary.txt

### Memory issues with large trajectories
Reduce `save_every` in config:
```yaml
outputs:
  save_every: 1000  # Save less frequently
```

## References

Implementation based on:
- GR 1PN formalism: Soffel (1989), Will & Wiseman (1996)
- Orbital elements: Murray & Dermott (1999)
- Precession measurement: Standard astrometry techniques
- Matplotlib best practices: https://matplotlib.org/

## Authors

Implemented as part of the superfluid orbit simulator project.
