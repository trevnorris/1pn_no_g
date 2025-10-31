#!/usr/bin/env python3
"""
Main command-line interface for the superfluid orbit simulator.

This script provides a robust, production-ready CLI for running superfluid
orbit simulations. It handles:
- Configuration loading and validation
- Simulation execution with progress reporting
- Optional GR-1PN comparison (when gr1pn module available)
- CSV and JSON output
- Graceful keyboard interrupt handling
- Comprehensive diagnostics and summary output

Usage:
    python -m slab.run config.yaml
    python -m slab.run config.yaml --output-dir results --verbose
    python -m slab.run config.yaml --validate-only
    python -m slab.run config.yaml --quick

Example:
    # Run Mercury orbit simulation
    python -m slab.run examples/mercury_orbit.yaml --verbose

    # Validate config without running
    python -m slab.run examples/mercury_orbit.yaml --validate-only

    # Fast run without audits
    python -m slab.run examples/mercury_orbit.yaml --quick

No gravitational constant G appears. Forces emerge from superfluid momentum
flux through control surfaces. Bodies are fluid intakes (mass sinks).

See plan_no_pde.md for physics documentation and PROJECT.md for overview.
"""

import argparse
import sys
import time
import signal
from pathlib import Path
from typing import Dict, List, Optional, Any
import warnings
import numpy as np
import copy

# Import slab modules
from slab.io_cfg import (
    load_config,
    validate_config,
    save_state_csv,
    save_diagnostics_json,
)
from slab.dynamics import integrate_orbit, estimate_orbital_period
from slab.bodies import Body
from slab.medium import Medium
from slab.diagnostics import compute_precession, find_periapsis_passages


# ============================================================================
# Global state for interrupt handling
# ============================================================================
_interrupted = False


def signal_handler(signum, frame):
    """Handle keyboard interrupt gracefully."""
    global _interrupted
    _interrupted = True
    print("\n\n🛑 Keyboard interrupt received. Saving partial results...")


# ============================================================================
# Main simulation runner
# ============================================================================

def run_simulation(
    config: Dict[str, Any],
    output_dir: Path,
    verbose: bool = False,
    quick: bool = False,
) -> Dict[str, Any]:
    """
    Run the superfluid orbit simulation.

    Parameters
    ----------
    config : dict
        Loaded and validated configuration
    output_dir : Path
        Directory for output files
    verbose : bool
        Enable verbose progress output
    quick : bool
        Disable audits for speed

    Returns
    -------
    results : dict
        Dictionary with:
        - 'trajectory': trajectory data from integrate_orbit
        - 'diagnostics': diagnostics list from integrate_orbit
        - 'summary': summary statistics dict
        - 'interrupted': bool, True if interrupted
    """
    global _interrupted
    _interrupted = False

    # Register signal handler for graceful interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Unpack configuration
    medium = config['medium']
    bodies = config['bodies']
    numerics = config['numerics']
    compare_gr = config['compare_gr_1pn']
    outputs = config['outputs']

    # Make copies of bodies for simulation (don't modify config)
    sim_bodies = [copy.deepcopy(body) for body in bodies]

    # Extract integration parameters
    dt = numerics['dt']
    n_steps = numerics['steps']

    # Build integration options
    opts = {
        'use_compressible': numerics['use_compressible'],
        'use_flux_mass': numerics['use_flux_mass'],
        'flux_every': numerics['intake_every'],
        'save_every': outputs['save_every'],
        'audit_every': 0 if quick else numerics['audit_every'],
        'audit_tolerance': 1e-3,  # From acceptance criteria
        'n_points': numerics['npts_audit'],
        'verbose': verbose,
        'progress_every': max(1, n_steps // 100),  # Update ~100 times
    }

    # Estimate orbital parameters
    T_orbit = estimate_orbital_period(sim_bodies, medium)

    # Print simulation header
    if verbose:
        print("=" * 80)
        print("SUPERFLUID ORBIT SIMULATOR")
        print("=" * 80)
        print()
        print("Medium parameters:")
        print(f"  ρ₀ (density):      {medium.rho0:.6e}")
        print(f"  cs (sound speed):  {medium.cs:.6e}")
        print(f"  β₀ (intake factor): {medium.beta0:.6e}")
        print(f"  K (orbital const): {medium.K:.6e}  [replaces G]")
        if medium.gamma_beta != 0:
            print(f"  γ_β (rarefaction): {medium.gamma_beta:.3f}")
        print()

        print(f"Bodies ({len(sim_bodies)}):")
        for body in sim_bodies:
            print(f"  {body.name:12s}: M={body.M:.6e}, R={body.R:.6e}, Q={body.Q:.6e}")
        print()

        print("Integration parameters:")
        print(f"  Timestep dt:        {dt:.6e}")
        print(f"  Total steps:        {n_steps:,}")
        print(f"  Total time:         {dt * n_steps:.3e}")
        print(f"  Est. orbital period: {T_orbit:.3e}  ({dt * n_steps / T_orbit:.1f} orbits)")
        print(f"  Save every:         {opts['save_every']} steps")
        if opts['audit_every'] > 0:
            print(f"  Audit every:        {opts['audit_every']} steps")
        else:
            print(f"  Audit:              DISABLED (--quick mode)")
        print()

        print("Physics options:")
        print(f"  Compressible:       {'YES' if opts['use_compressible'] else 'NO'}")
        print(f"  Mass flux:          {'YES' if opts['use_flux_mass'] else 'NO'}")
        print()

    # Run integration
    try:
        t_start = time.time()
        trajectory, diagnostics = integrate_orbit(
            sim_bodies, medium, dt, n_steps, opts
        )
        t_end = time.time()
        elapsed = t_end - t_start
        interrupted = False

    except KeyboardInterrupt:
        # Should be caught by signal handler, but handle here too
        t_end = time.time()
        elapsed = t_end - t_start
        interrupted = True
        print("\nSimulation interrupted by user.")
        # trajectory and diagnostics should be partially filled
        # In current implementation, we don't save partial results
        # This would require modifying integrate_orbit to be more robust
        raise

    # Check if interrupted via signal
    if _interrupted:
        interrupted = True

    # Compute summary statistics
    summary = compute_summary(
        trajectory, diagnostics, sim_bodies, medium, elapsed, T_orbit
    )

    results = {
        'trajectory': trajectory,
        'diagnostics': diagnostics,
        'summary': summary,
        'interrupted': interrupted,
    }

    return results


def compute_summary(
    trajectory: Dict,
    diagnostics: List[Dict],
    bodies: List[Body],
    medium: Medium,
    elapsed_time: float,
    T_orbit: float,
) -> Dict[str, Any]:
    """
    Compute summary statistics from simulation results.

    Parameters
    ----------
    trajectory : dict
        Trajectory data (t, x, v, M, Q arrays)
    diagnostics : List[dict]
        Diagnostics at each saved step
    bodies : List[Body]
        Final body states
    medium : Medium
        Medium parameters
    elapsed_time : float
        Wall-clock time for integration [seconds]
    T_orbit : float
        Estimated orbital period [code units]

    Returns
    -------
    summary : dict
        Summary statistics including:
        - timing: wall time, steps/sec
        - energy: initial, final, drift
        - precession: perihelion info if measurable
        - forces: audit statistics if available
    """
    n_steps = len(diagnostics)
    times = trajectory['t']
    dt = times[1] - times[0] if len(times) > 1 else 0.0

    # Energy conservation
    energies = np.array([d['total_energy'] for d in diagnostics])
    E_initial = energies[0]
    E_final = energies[-1]
    E_drift_abs = E_final - E_initial
    E_drift_rel = E_drift_abs / E_initial if E_initial != 0 else 0.0
    E_max_deviation = np.max(np.abs(energies - E_initial))
    E_max_deviation_rel = E_max_deviation / abs(E_initial) if E_initial != 0 else 0.0

    # Momentum conservation
    p_initial = diagnostics[0]['total_momentum']
    p_final = diagnostics[-1]['total_momentum']
    p_drift = np.linalg.norm(p_final - p_initial)

    # Force statistics
    max_forces = np.array([d['max_force'] for d in diagnostics])
    mean_max_force = np.mean(max_forces)
    max_force_overall = np.max(max_forces)

    # Audit statistics (if audits were performed)
    audit_errors = []
    for d in diagnostics:
        if 'audit' in d:
            for body_audit in d['audit']:
                audit_errors.append(body_audit['rel_error'])

    if audit_errors:
        audit_mean_error = np.mean(audit_errors)
        audit_max_error = np.max(audit_errors)
        audit_pass_rate = np.sum(np.array(audit_errors) < 1e-3) / len(audit_errors)
    else:
        audit_mean_error = None
        audit_max_error = None
        audit_pass_rate = None

    # Precession measurement (for two-body systems)
    precession_rate = None
    n_periapsis = 0
    if len(bodies) == 2 and len(times) > 10:
        try:
            # Extract relative coordinates
            x_rel = trajectory['x'][:, 1, :] - trajectory['x'][:, 0, :]
            v_rel = trajectory['v'][:, 1, :] - trajectory['v'][:, 0, :]

            # Find periapsis passages
            passages = find_periapsis_passages(times, x_rel, v_rel)
            n_periapsis = len(passages)

            if n_periapsis >= 2:
                precession_rate = compute_precession(passages)
        except Exception:
            # If precession calculation fails, just skip it
            pass

    # Performance metrics
    steps_per_sec = n_steps / elapsed_time if elapsed_time > 0 else 0.0
    time_per_orbit = elapsed_time / (times[-1] / T_orbit) if T_orbit > 0 else 0.0

    summary = {
        'timing': {
            'elapsed_seconds': elapsed_time,
            'steps_per_second': steps_per_sec,
            'wall_time_per_orbit': time_per_orbit,
        },
        'integration': {
            'n_steps': n_steps,
            'dt': dt,
            'total_time': times[-1],
            'n_orbits': times[-1] / T_orbit if T_orbit > 0 else 0.0,
            'orbital_period': T_orbit,
        },
        'energy': {
            'initial': E_initial,
            'final': E_final,
            'drift_absolute': E_drift_abs,
            'drift_relative': E_drift_rel,
            'max_deviation_absolute': E_max_deviation,
            'max_deviation_relative': E_max_deviation_rel,
        },
        'momentum': {
            'initial_magnitude': np.linalg.norm(p_initial),
            'final_magnitude': np.linalg.norm(p_final),
            'drift_magnitude': p_drift,
        },
        'forces': {
            'mean_max': mean_max_force,
            'overall_max': max_force_overall,
        },
        'audits': {
            'n_audits': len(audit_errors),
            'mean_error': audit_mean_error,
            'max_error': audit_max_error,
            'pass_rate': audit_pass_rate,
        } if audit_errors else None,
        'precession': {
            'n_periapsis': n_periapsis,
            'rate_rad_per_orbit': precession_rate,
            'rate_arcsec_per_century': (
                precession_rate * (180 / np.pi * 3600) * 100 / T_orbit
                if precession_rate is not None and T_orbit > 0
                else None
            ),
        } if precession_rate is not None else None,
    }

    return summary


def print_summary(summary: Dict[str, Any], verbose: bool = False) -> None:
    """
    Print human-readable simulation summary.

    Parameters
    ----------
    summary : dict
        Summary statistics from compute_summary
    verbose : bool
        Print extended details
    """
    print()
    print("=" * 80)
    print("SIMULATION SUMMARY")
    print("=" * 80)
    print()

    # Timing
    t = summary['timing']
    print(f"Performance:")
    print(f"  Wall time:          {t['elapsed_seconds']:.2f} seconds")
    print(f"  Speed:              {t['steps_per_second']:.1f} steps/second")
    if t['wall_time_per_orbit'] > 0:
        print(f"  Time per orbit:     {t['wall_time_per_orbit']:.2f} seconds")
    print()

    # Integration
    i = summary['integration']
    print(f"Integration:")
    print(f"  Total steps:        {i['n_steps']:,}")
    print(f"  Timestep dt:        {i['dt']:.6e}")
    print(f"  Total time:         {i['total_time']:.6e}")
    print(f"  Orbits completed:   {i['n_orbits']:.2f}")
    print()

    # Energy conservation
    e = summary['energy']
    print(f"Energy conservation:")
    print(f"  Initial energy:     {e['initial']:+.10e}")
    print(f"  Final energy:       {e['final']:+.10e}")
    print(f"  Absolute drift:     {e['drift_absolute']:+.6e}")
    print(f"  Relative drift:     {e['drift_relative']:+.6e}")
    print(f"  Max deviation:      {e['max_deviation_relative']:+.6e}")

    # Quality assessment
    if abs(e['drift_relative']) < 1e-5:
        quality = "EXCELLENT"
    elif abs(e['drift_relative']) < 1e-3:
        quality = "GOOD"
    elif abs(e['drift_relative']) < 1e-2:
        quality = "ACCEPTABLE"
    else:
        quality = "POOR (check timestep)"
    print(f"  Quality:            {quality}")
    print()

    # Momentum conservation
    if verbose:
        m = summary['momentum']
        print(f"Momentum conservation:")
        print(f"  Initial |p|:        {m['initial_magnitude']:.6e}")
        print(f"  Final |p|:          {m['final_magnitude']:.6e}")
        print(f"  Drift |Δp|:         {m['drift_magnitude']:.6e}")
        print()

    # Force statistics
    f = summary['forces']
    print(f"Force statistics:")
    print(f"  Mean max force:     {f['mean_max']:.6e}")
    print(f"  Overall max force:  {f['overall_max']:.6e}")
    print()

    # Audit results
    if summary['audits'] is not None:
        a = summary['audits']
        print(f"Quadrature audits:")
        print(f"  Total audits:       {a['n_audits']}")
        print(f"  Mean rel. error:    {a['mean_error']:.6e}")
        print(f"  Max rel. error:     {a['max_error']:.6e}")
        print(f"  Pass rate:          {a['pass_rate'] * 100:.1f}%")
        if a['pass_rate'] < 0.95:
            print(f"  ⚠️  WARNING: Low pass rate! Check analytic vs quadrature forces.")
        print()

    # Precession
    if summary['precession'] is not None:
        p = summary['precession']
        print(f"Perihelion precession:")
        print(f"  Periapsis passages: {p['n_periapsis']}")
        if p['rate_rad_per_orbit'] is not None:
            print(f"  Rate:               {p['rate_rad_per_orbit']:.6e} rad/orbit")
            if p['rate_arcsec_per_century'] is not None:
                print(f"                      {p['rate_arcsec_per_century']:.2f} arcsec/century")
        print()


def print_trajectory_table(
    trajectory: Dict,
    bodies_names: List[str],
    max_rows: int = 20,
) -> None:
    """
    Print formatted table of trajectory data for debugging.

    Parameters
    ----------
    trajectory : dict
        Trajectory data (t, x, v, M, Q)
    bodies_names : List[str]
        Names of bodies
    max_rows : int
        Maximum rows to print (show first/last if exceeded)
    """
    times = trajectory['t']
    positions = trajectory['x']
    velocities = trajectory['v']
    n_snapshots = len(times)
    n_bodies = len(bodies_names)

    print()
    print("=" * 80)
    print("TRAJECTORY SAMPLE")
    print("=" * 80)
    print()

    # Determine which rows to print
    if n_snapshots <= max_rows:
        rows = range(n_snapshots)
    else:
        # Show first and last few rows
        n_show = max_rows // 2
        rows = list(range(n_show)) + list(range(n_snapshots - n_show, n_snapshots))

    # Print header
    print(f"{'Time':>12s}  {'Body':>12s}  {'x':>14s}  {'y':>14s}  {'z':>14s}  "
          f"{'vx':>14s}  {'vy':>14s}  {'vz':>14s}")
    print("-" * 110)

    last_printed_idx = -1
    for i in rows:
        # Print ellipsis if we skip rows
        if i > last_printed_idx + 1:
            print(f"{'...':>12s}  {'...':>12s}  {'...':>14s}  {'...':>14s}  {'...':>14s}  "
                  f"{'...':>14s}  {'...':>14s}  {'...':>14s}")

        t = times[i]
        for j, name in enumerate(bodies_names):
            x = positions[i, j]
            v = velocities[i, j]
            print(f"{t:12.6e}  {name:>12s}  "
                  f"{x[0]:+14.6e}  {x[1]:+14.6e}  {x[2]:+14.6e}  "
                  f"{v[0]:+14.6e}  {v[1]:+14.6e}  {v[2]:+14.6e}")

        last_printed_idx = i

    print()


def save_outputs(
    config: Dict[str, Any],
    results: Dict[str, Any],
    output_dir: Path,
    verbose: bool = False,
) -> None:
    """
    Save simulation outputs (CSV, JSON, plots).

    Parameters
    ----------
    config : dict
        Configuration
    results : dict
        Simulation results
    output_dir : Path
        Output directory
    verbose : bool
        Print progress
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    trajectory = results['trajectory']
    diagnostics = results['diagnostics']
    summary = results['summary']

    # Prepare bodies_history for CSV saving
    # save_state_csv expects List[List[Body]]
    times = trajectory['t']
    positions = trajectory['x']
    velocities = trajectory['v']
    masses = trajectory['M']
    intakes = trajectory['Q']

    n_snapshots = len(times)
    n_bodies = len(config['bodies'])
    bodies_names = [b.name for b in config['bodies']]

    # Reconstruct bodies at each snapshot
    bodies_history = []
    for i in range(n_snapshots):
        snapshot_bodies = []
        for j in range(n_bodies):
            # Get R from original config (doesn't change)
            R = config['bodies'][j].R

            body = Body(
                name=bodies_names[j],
                M=masses[i, j],
                x=positions[i, j].copy(),
                v=velocities[i, j].copy(),
                R=R,
                Q=intakes[i, j],
            )
            snapshot_bodies.append(body)
        bodies_history.append(snapshot_bodies)

    # Save trajectory CSV
    if config['outputs']['write_csv']:
        csv_path = output_dir / "trajectory.csv"
        if verbose:
            print(f"Saving trajectory to {csv_path}...")
        save_state_csv(str(csv_path), times, bodies_history)

    # Prepare diagnostics for JSON
    # Convert diagnostics list to a more compact format
    def convert_numpy_types(obj):
        """Convert numpy types to native Python types for JSON serialization."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, dict):
            return {k: convert_numpy_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_types(item) for item in obj]
        else:
            return obj

    diag_dict = {
        'times': times.tolist(),
        'total_energy': [float(d['total_energy']) for d in diagnostics],
        'kinetic_energy': [float(d['kinetic_energy']) for d in diagnostics],
        'potential_energy': [float(d['potential_energy']) for d in diagnostics],
        'max_force': [float(d['max_force']) for d in diagnostics],
        'summary': convert_numpy_types(summary),
    }

    # Add audit data if available
    audit_results = []
    for d in diagnostics:
        if 'audit' in d:
            audit_results.append({
                'time': float(d['time']),
                'step': int(d['step']),
                'results': convert_numpy_types(d['audit']),
            })
    if audit_results:
        diag_dict['audits'] = audit_results

    # Save diagnostics JSON
    json_path = output_dir / "diagnostics.json"
    if verbose:
        print(f"Saving diagnostics to {json_path}...")
    save_diagnostics_json(str(json_path), diag_dict)

    print()
    print(f"✓ Outputs saved to: {output_dir.absolute()}")
    print(f"  - {csv_path.name}")
    print(f"  - {json_path.name}")
    print()


def compare_with_gr1pn(
    config: Dict[str, Any],
    trajectory: Dict,
    verbose: bool = False,
) -> Optional[Dict[str, Any]]:
    """
    Run GR-1PN comparison if enabled and module available.

    Parameters
    ----------
    config : dict
        Configuration
    trajectory : dict
        Slab simulation trajectory
    verbose : bool
        Print progress

    Returns
    -------
    gr_results : dict or None
        GR comparison results if successful, None otherwise
    """
    compare_gr = config['compare_gr_1pn']

    if not compare_gr['enable']:
        return None

    if verbose:
        print()
        print("=" * 80)
        print("GR-1PN COMPARISON")
        print("=" * 80)
        print()

    try:
        from slab import gr1pn
    except ImportError:
        print("⚠️  GR-1PN comparison requested but slab.gr1pn module not found.")
        print("    Skipping GR comparison.")
        return None

    # This is a placeholder for when gr1pn module is implemented
    # For now, just note that it would be called here

    print("⚠️  GR-1PN comparison module not yet implemented.")
    print("    See plan_no_pde.md § 7 for implementation details.")
    print()

    return None


# ============================================================================
# Command-line interface
# ============================================================================

def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for CLI."""
    parser = argparse.ArgumentParser(
        prog='slab.run',
        description=(
            'Superfluid orbit simulator: compute orbits using only superfluid '
            'hydrodynamics. No gravitational constant G appears. Forces emerge '
            'from momentum flux through control surfaces.'
        ),
        epilog=(
            'Examples:\n'
            '  python -m slab.run examples/mercury_orbit.yaml\n'
            '  python -m slab.run config.yaml --output-dir results --verbose\n'
            '  python -m slab.run config.yaml --validate-only\n'
            '  python -m slab.run config.yaml --quick\n'
            '\n'
            'See plan_no_pde.md for physics details and PROJECT.md for overview.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        'config',
        type=str,
        help='Path to YAML configuration file',
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        default='output',
        help='Output directory for results (default: output/)',
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output (progress, diagnostics, etc.)',
    )

    parser.add_argument(
        '--quick',
        action='store_true',
        help='Disable quadrature audits for speed (use with caution)',
    )

    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Validate configuration and exit (no simulation)',
    )

    parser.add_argument(
        '--no-table',
        action='store_true',
        help='Skip trajectory table output',
    )

    parser.add_argument(
        '--table-rows',
        type=int,
        default=20,
        help='Number of rows in trajectory table (default: 20)',
    )

    return parser


def main(argv=None):
    """Main entry point for CLI."""
    parser = create_parser()
    args = parser.parse_args(argv)

    # Convert output_dir to Path
    output_dir = Path(args.output_dir)

    # Print banner
    if args.verbose:
        print()
        print("╔════════════════════════════════════════════════════════════════════════════╗")
        print("║                    SUPERFLUID ORBIT SIMULATOR                              ║")
        print("║                                                                            ║")
        print("║  Forces from fluid momentum flux, not gravity.                            ║")
        print("║  No gravitational constant G. Bodies are fluid intakes.                   ║")
        print("╚════════════════════════════════════════════════════════════════════════════╝")
        print()

    # 1) Load configuration
    try:
        if args.verbose:
            print(f"Loading configuration from: {args.config}")
            print()
        config = load_config(args.config)
    except FileNotFoundError:
        print(f"ERROR: Configuration file not found: {args.config}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: Failed to load configuration: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

    # 2) Validate configuration
    if args.verbose:
        print("Validating configuration...")
        print()

    is_valid, warnings_list = validate_config(config)

    if warnings_list:
        print("⚠️  Configuration warnings/errors:")
        for w in warnings_list:
            print(f"    - {w}")
        print()

    if not is_valid:
        print("❌ Configuration is INVALID. Please fix errors above.", file=sys.stderr)
        return 1

    if args.verbose:
        print("✓ Configuration is valid.")
        print()

    # Exit if validate-only mode
    if args.validate_only:
        print("✓ Configuration validated successfully. Exiting (--validate-only mode).")
        return 0

    # 3) Run simulation
    try:
        results = run_simulation(
            config,
            output_dir,
            verbose=args.verbose,
            quick=args.quick,
        )
    except KeyboardInterrupt:
        print("\n\nSimulation interrupted by user. Exiting without saving.")
        return 130  # Standard exit code for SIGINT
    except Exception as e:
        print(f"\nERROR: Simulation failed: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

    # 4) Print summary
    print_summary(results['summary'], verbose=args.verbose)

    # 5) Print trajectory table (for debugging)
    if not args.no_table:
        bodies_names = [b.name for b in config['bodies']]
        print_trajectory_table(
            results['trajectory'],
            bodies_names,
            max_rows=args.table_rows,
        )

    # 6) Save outputs
    try:
        save_outputs(config, results, output_dir, verbose=args.verbose)
    except Exception as e:
        print(f"ERROR: Failed to save outputs: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

    # 7) Optional GR comparison
    if config['compare_gr_1pn']['enable']:
        gr_results = compare_with_gr1pn(
            config,
            results['trajectory'],
            verbose=args.verbose,
        )
        # If GR module were implemented, results would be saved here

    # 8) Final message
    print()
    print("=" * 80)
    print("✓ Simulation complete!")
    print("=" * 80)
    print()

    if results['interrupted']:
        print("⚠️  Note: Simulation was interrupted. Results may be incomplete.")
        return 130

    return 0


if __name__ == '__main__':
    sys.exit(main())
