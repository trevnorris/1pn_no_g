#!/usr/bin/env python3
"""Compare slab simulation with GR 1PN predictions.

This script runs both slab and GR-1PN simulations with identical initial
conditions and generates comparison plots to validate that the slab model
correctly reproduces General Relativity effects.

Usage:
------
    python scripts/compare_with_gr.py --config examples/quick_validation.yaml
    python scripts/compare_with_gr.py --steps 1000 --output-dir output/comparison
    python scripts/compare_with_gr.py --config examples/mercury_orbit.yaml --dpi 300

The script performs:
1. Load configuration from YAML file
2. Run slab simulation with configured parameters
3. Run GR 1PN simulation with same initial conditions
4. Generate comparison plots (orbit, precession, 3D)
5. Print summary statistics (max separation, precession agreement)

Output files:
-------------
- orbit_comparison.png: 2-panel plot (trajectories + separation)
- precession_comparison.png: omega vs time with linear fits
- trajectory_3d_slab.png: 3D visualization of slab trajectory
- trajectory_3d_gr.png: 3D visualization of GR trajectory
- summary.txt: Text file with statistics

Requirements:
-------------
- matplotlib (for plotting): pip install matplotlib
- Configuration file with GR comparison enabled
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import argparse
from typing import Dict
import copy

from slab.io_cfg import load_config, validate_config
from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit
from slab.gr1pn import integrate_gr1pn_orbit, compute_precession_rate

# Import visualization functions
try:
    from slab.viz import (
        plot_orbit_comparison,
        plot_precession_comparison,
        plot_trajectory_3d
    )
    HAS_VIZ = True
except ImportError:
    HAS_VIZ = False
    print("Warning: Visualization module not available. Install matplotlib.")


def run_slab_simulation(config: Dict) -> Dict:
    """Run slab simulation from configuration.

    Parameters
    ----------
    config : Dict
        Configuration dictionary from load_config()

    Returns
    -------
    trajectory : Dict
        Trajectory data with 't', 'x', 'v', 'M', 'Q' arrays
    """
    medium = config['medium']
    bodies = config['bodies']
    numerics = config['numerics']

    # Prepare integration options
    opts = {
        'use_compressible': numerics['use_compressible'],
        'use_quadrature': False,  # Use analytic forces for speed
        'n_points': numerics['npts_audit'],
        'use_flux_mass': numerics['use_flux_mass'],
        'flux_every': numerics['intake_every'],
        'audit_every': numerics.get('audit_every', 0),
        'audit_tolerance': 1e-3,
        'save_every': config['outputs']['save_every'],
        'verbose': True,
        'progress_every': max(1, numerics['steps'] // 10),
    }

    print("=" * 70)
    print("RUNNING SLAB SIMULATION")
    print("=" * 70)
    print(f"Medium: rho0={medium.rho0}, cs={medium.cs}, beta0={medium.beta0}")
    print(f"  K = {medium.K:.6e} (orbital constant)")
    print(f"Bodies: {len(bodies)}")
    for body in bodies:
        print(f"  {body.name}: M={body.M:.3e}, x={body.x}, v={body.v}")
    print(f"Integration: dt={numerics['dt']}, steps={numerics['steps']}")
    print(f"  Total time: {numerics['dt'] * numerics['steps']:.3f} yr")
    print(f"Compressible corrections: {opts['use_compressible']}")
    print()

    # Run integration
    trajectory, diagnostics = integrate_orbit(
        bodies, medium, numerics['dt'], numerics['steps'], opts
    )

    print()
    print("Slab simulation complete!")
    print()

    return trajectory


def run_gr_simulation(config: Dict) -> Dict:
    """Run GR 1PN simulation with same initial conditions.

    Parameters
    ----------
    config : Dict
        Configuration dictionary from load_config()

    Returns
    -------
    trajectory : Dict
        GR trajectory data with 't', 'x', 'v', 'a' arrays
    """
    # Get GR parameters
    gr_config = config['compare_gr_1pn']
    c_light = gr_config['c_light']

    # Standard gravitational constant (SI units would be 6.6743e-11)
    # But we use code units where G matches K from slab
    # To compare like-for-like, we use K from slab as G for GR
    medium = config['medium']
    G_newton = medium.K  # Use same orbital constant

    # Copy bodies to avoid modifying original
    bodies = [copy.deepcopy(body) for body in config['bodies']]

    # Integration parameters
    numerics = config['numerics']
    dt = numerics['dt']
    n_steps = numerics['steps']
    save_every = config['outputs']['save_every']

    print("=" * 70)
    print("RUNNING GR 1PN SIMULATION")
    print("=" * 70)
    print(f"GR parameters: G={G_newton:.6e}, c={c_light:.6e}")
    print(f"Bodies: {len(bodies)}")
    for body in bodies:
        print(f"  {body.name}: M={body.M:.3e}, x={body.x}, v={body.v}")
    print(f"Integration: dt={dt}, steps={n_steps}")
    print(f"  Total time: {dt * n_steps:.3f} yr")
    print()

    # Run integration
    trajectory = integrate_gr1pn_orbit(
        bodies, c_light, G_newton, dt, n_steps, save_every
    )

    print()
    print("GR simulation complete!")
    print()

    return trajectory


def compute_statistics(
    slab_trajectory: Dict,
    gr_trajectory: Dict,
    body_idx: int,
    M_central: float,
    K_slab: float,
    G_gr: float
) -> Dict:
    """Compute comparison statistics.

    Parameters
    ----------
    slab_trajectory : Dict
        Slab trajectory
    gr_trajectory : Dict
        GR trajectory
    body_idx : int
        Index of orbiting body
    M_central : float
        Mass of central body
    K_slab : float
        Slab orbital constant
    G_gr : float
        GR gravitational constant

    Returns
    -------
    stats : Dict
        Dictionary with statistics
    """
    from slab.gr1pn import compute_orbit_elements

    # Extract positions
    t_slab = slab_trajectory['t']
    x_slab = slab_trajectory['x'][:, body_idx, :]
    v_slab = slab_trajectory['v'][:, body_idx, :]

    t_gr = gr_trajectory['t']
    x_gr = gr_trajectory['x'][:, body_idx, :]
    v_gr = gr_trajectory['v'][:, body_idx, :]

    # Compute separation (interpolate if needed)
    if len(t_slab) != len(t_gr) or not np.allclose(t_slab, t_gr):
        x_gr_interp = np.zeros_like(x_slab)
        for dim in range(3):
            x_gr_interp[:, dim] = np.interp(t_slab, t_gr, x_gr[:, dim])
        separation = np.linalg.norm(x_slab - x_gr_interp, axis=1)
    else:
        separation = np.linalg.norm(x_slab - x_gr, axis=1)

    max_separation = np.max(separation)
    mean_separation = np.mean(separation)
    final_separation = separation[-1]

    # Compute precession rates
    # Slab precession
    omega_slab = []
    a_slab = []
    for i in range(len(t_slab)):
        elems = compute_orbit_elements(x_slab[i], v_slab[i], M_central, K_slab)
        omega_slab.append(elems['omega'])
        if np.isfinite(elems['a']):
            a_slab.append(elems['a'])

    omega_slab = np.array(omega_slab)
    omega_slab_unwrapped = np.unwrap(omega_slab)
    coeffs_slab = np.polyfit(t_slab, omega_slab_unwrapped, deg=1)
    domega_dt_slab = coeffs_slab[0]

    # GR precession
    omega_gr = []
    a_gr = []
    for i in range(len(t_gr)):
        elems = compute_orbit_elements(x_gr[i], v_gr[i], M_central, G_gr)
        omega_gr.append(elems['omega'])
        if np.isfinite(elems['a']):
            a_gr.append(elems['a'])

    omega_gr = np.array(omega_gr)
    omega_gr_unwrapped = np.unwrap(omega_gr)
    coeffs_gr = np.polyfit(t_gr, omega_gr_unwrapped, deg=1)
    domega_dt_gr = coeffs_gr[0]

    # Orbital period (use slab)
    a_mean = np.mean(a_slab) if a_slab else 1.0
    P_orb = 2 * np.pi * np.sqrt(a_mean**3 / (K_slab * M_central))
    n_orbits = (t_slab[-1] - t_slab[0]) / P_orb

    # Precession per orbit
    domega_per_orbit_slab = domega_dt_slab * P_orb
    domega_per_orbit_gr = domega_dt_gr * P_orb

    # Agreement ratio
    if domega_per_orbit_gr != 0:
        precession_ratio = domega_per_orbit_slab / domega_per_orbit_gr
    else:
        precession_ratio = np.nan

    return {
        'max_separation': max_separation,
        'mean_separation': mean_separation,
        'final_separation': final_separation,
        'slab_precession_per_orbit': domega_per_orbit_slab,
        'gr_precession_per_orbit': domega_per_orbit_gr,
        'precession_ratio': precession_ratio,
        'orbital_period': P_orb,
        'n_orbits': n_orbits,
        'semi_major_axis': a_mean,
    }


def save_summary(stats: Dict, output_path: str) -> None:
    """Save summary statistics to text file.

    Parameters
    ----------
    stats : Dict
        Statistics dictionary
    output_path : str
        Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("SLAB vs GR 1PN COMPARISON SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write("ORBITAL PARAMETERS\n")
        f.write("-" * 70 + "\n")
        f.write(f"Semi-major axis: {stats['semi_major_axis']:.6e} AU\n")
        f.write(f"Orbital period:  {stats['orbital_period']:.6e} yr\n")
        f.write(f"                 ({stats['orbital_period'] * 365.25:.2f} days)\n")
        f.write(f"Number of orbits: {stats['n_orbits']:.2f}\n")
        f.write("\n")

        f.write("POSITION DIFFERENCE\n")
        f.write("-" * 70 + "\n")
        f.write(f"Maximum separation:  {stats['max_separation']:.6e} AU\n")
        f.write(f"Mean separation:     {stats['mean_separation']:.6e} AU\n")
        f.write(f"Final separation:    {stats['final_separation']:.6e} AU\n")
        f.write("\n")

        f.write("PERIHELION PRECESSION\n")
        f.write("-" * 70 + "\n")
        f.write(f"Slab precession:  {stats['slab_precession_per_orbit']:.6e} rad/orbit\n")
        f.write(f"                  {stats['slab_precession_per_orbit'] * 206265:.3f} arcsec/orbit\n")
        f.write(f"GR precession:    {stats['gr_precession_per_orbit']:.6e} rad/orbit\n")
        f.write(f"                  {stats['gr_precession_per_orbit'] * 206265:.3f} arcsec/orbit\n")
        f.write(f"Agreement ratio:  {stats['precession_ratio']:.6f}\n")
        f.write(f"Percent difference: {abs(stats['precession_ratio'] - 1.0) * 100:.2f}%\n")
        f.write("\n")

        f.write("INTERPRETATION\n")
        f.write("-" * 70 + "\n")
        if abs(stats['precession_ratio'] - 1.0) < 0.05:
            f.write("EXCELLENT: Precession agreement within 5%\n")
        elif abs(stats['precession_ratio'] - 1.0) < 0.10:
            f.write("GOOD: Precession agreement within 10%\n")
        elif abs(stats['precession_ratio'] - 1.0) < 0.20:
            f.write("FAIR: Precession agreement within 20%\n")
        else:
            f.write("POOR: Precession differs by >20%\n")
        f.write("\n")

    print(f"Saved summary statistics to {output_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Compare slab simulation with GR 1PN predictions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--config',
        type=str,
        default='examples/quick_validation.yaml',
        help='Configuration file path (default: examples/quick_validation.yaml)'
    )
    parser.add_argument(
        '--steps',
        type=int,
        default=None,
        help='Override number of steps from config'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output/comparison',
        help='Output directory for plots and summary (default: output/comparison)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='Plot resolution (default: 150)'
    )
    parser.add_argument(
        '--body-idx',
        type=int,
        default=1,
        help='Index of orbiting body to analyze (default: 1)'
    )

    args = parser.parse_args()

    # Check if visualization is available
    if not HAS_VIZ:
        print("ERROR: matplotlib not installed. Install with: pip install matplotlib")
        sys.exit(1)

    # Load configuration
    print("=" * 70)
    print("LOADING CONFIGURATION")
    print("=" * 70)
    print(f"Config file: {args.config}")

    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"ERROR loading config: {e}")
        sys.exit(1)

    # Validate configuration
    is_valid, warnings = validate_config(config)
    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  - {w}")
    if not is_valid:
        print("\nERROR: Configuration is invalid!")
        sys.exit(1)

    print("\nConfiguration is valid.")

    # Check that GR comparison is enabled
    if not config['compare_gr_1pn']['enable']:
        print("\nWARNING: GR comparison is not enabled in config!")
        print("Enabling it for this run...")
        config['compare_gr_1pn']['enable'] = True

    # Override steps if requested
    if args.steps is not None:
        print(f"\nOverriding steps: {config['numerics']['steps']} -> {args.steps}")
        config['numerics']['steps'] = args.steps

    print()

    # Run slab simulation
    slab_trajectory = run_slab_simulation(config)

    # Run GR simulation
    gr_trajectory = run_gr_simulation(config)

    # Compute statistics
    print("=" * 70)
    print("COMPUTING COMPARISON STATISTICS")
    print("=" * 70)

    medium = config['medium']
    body_idx = args.body_idx
    M_central = config['bodies'][0].M  # Assume first body is central
    K_slab = medium.K
    G_gr = medium.K  # Use same value for fair comparison

    stats = compute_statistics(
        slab_trajectory, gr_trajectory,
        body_idx, M_central, K_slab, G_gr
    )

    print(f"\nMax separation: {stats['max_separation']:.3e} AU")
    print(f"Mean separation: {stats['mean_separation']:.3e} AU")
    print(f"Slab precession: {stats['slab_precession_per_orbit'] * 206265:.3f} arcsec/orbit")
    print(f"GR precession:   {stats['gr_precession_per_orbit'] * 206265:.3f} arcsec/orbit")
    print(f"Agreement ratio: {stats['precession_ratio']:.6f}")
    print()

    # Generate plots
    print("=" * 70)
    print("GENERATING COMPARISON PLOTS")
    print("=" * 70)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Orbit comparison plot
    print("\n1. Orbit comparison (2-panel)...")
    plot_orbit_comparison(
        slab_trajectory, gr_trajectory,
        output_path=output_dir / "orbit_comparison.png",
        body_idx=body_idx,
        dpi=args.dpi
    )

    # Precession comparison plot
    print("\n2. Precession comparison...")
    plot_precession_comparison(
        slab_trajectory, gr_trajectory,
        body_idx=body_idx,
        M_central=M_central,
        K_slab=K_slab,
        G_gr=G_gr,
        output_path=output_dir / "precession_comparison.png",
        dpi=args.dpi
    )

    # 3D trajectory plots
    print("\n3. 3D trajectory (slab)...")
    body_names = [body.name for body in config['bodies']]
    plot_trajectory_3d(
        slab_trajectory,
        body_names=body_names,
        output_path=output_dir / "trajectory_3d_slab.png",
        dpi=args.dpi
    )

    print("\n4. 3D trajectory (GR)...")
    plot_trajectory_3d(
        gr_trajectory,
        body_names=body_names,
        output_path=output_dir / "trajectory_3d_gr.png",
        dpi=args.dpi
    )

    # Save summary
    print("\n5. Summary statistics...")
    save_summary(stats, output_dir / "summary.txt")

    print()
    print("=" * 70)
    print("COMPARISON COMPLETE!")
    print("=" * 70)
    print(f"\nOutput files saved to: {output_dir}")
    print(f"  - orbit_comparison.png")
    print(f"  - precession_comparison.png")
    print(f"  - trajectory_3d_slab.png")
    print(f"  - trajectory_3d_gr.png")
    print(f"  - summary.txt")
    print()


if __name__ == "__main__":
    main()
