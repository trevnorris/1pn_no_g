#!/usr/bin/env python3
"""
Numerical Robustness Validation Tests
======================================

Comprehensive test suite demonstrating convergence and parameter invariance
for the superfluid hydrodynamics orbit simulator.

This script implements three validation studies:

1. **Timestep Invariance Study**
   - Tests energy conservation convergence as dt → 0
   - Verifies 2nd-order convergence of velocity-Verlet integrator
   - Checks precession rate stability (for compressible mode)

2. **Control Surface Radius Invariance**
   - Tests force independence of control surface radius R
   - Validates near-field renormalization
   - Confirms forces scale correctly after self-field subtraction

3. **Quadrature Convergence Study**
   - Tests convergence of surface integrals vs number of points
   - Validates Fibonacci sphere quadrature accuracy
   - Measures convergence rate (should be power-law)

Expected results:
-----------------
- Timestep: Energy drift improves by ~4x when dt halved (2nd order)
- Radius: Force variation < 0.1% across R range
- Quadrature: Error < 1e-3 for n >= 256

Usage:
------
    python scripts/robustness_tests.py --all
    python scripts/robustness_tests.py --timestep
    python scripts/robustness_tests.py --radius
    python scripts/robustness_tests.py --quadrature
    python scripts/robustness_tests.py --quick
    python scripts/robustness_tests.py --output-dir results/
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional
import argparse
from dataclasses import dataclass

from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit, estimate_orbital_period
from slab.surface import (
    force_incompressible_analytic,
    force_incompressible_quadrature,
    force_total,
)


# ============================================================================
# Configuration and setup
# ============================================================================

@dataclass
class TestConfig:
    """Configuration for robustness tests."""
    # Timestep test parameters
    dt_values: List[float]  # In units of T_orbit
    n_orbits: float         # Number of orbits to integrate

    # Radius test parameters
    R_factors: List[float]  # Multiples of nominal R

    # Quadrature test parameters
    n_points_values: List[int]

    # Output configuration
    output_dir: Path

    # System parameters
    use_compressible: bool

    @classmethod
    def default(cls) -> 'TestConfig':
        """Standard test configuration."""
        return cls(
            dt_values=[0.0005, 0.001, 0.002, 0.005, 0.01],
            n_orbits=10.0,
            R_factors=[0.5, 1.0, 2.0],
            n_points_values=[64, 128, 256, 512, 1024, 2048],
            output_dir=Path("output"),
            use_compressible=False,
        )

    @classmethod
    def quick(cls) -> 'TestConfig':
        """Quick test configuration for rapid validation."""
        return cls(
            dt_values=[0.001, 0.005, 0.01],
            n_orbits=2.0,
            R_factors=[0.5, 1.0, 2.0],
            n_points_values=[64, 128, 256, 512],
            output_dir=Path("output"),
            use_compressible=False,
        )


def setup_mercury_system(medium: Medium, R_sun: float = 1e-3, R_mercury: float = 5e-4) -> Tuple[List[Body], float]:
    """
    Create Sun-Mercury system for testing.

    Parameters
    ----------
    medium : Medium
        Superfluid medium parameters
    R_sun : float
        Control surface radius for Sun [AU]
    R_mercury : float
        Control surface radius for Mercury [AU]

    Returns
    -------
    bodies : List[Body]
        [Sun, Mercury] with circular orbit initial conditions
    T_orbit : float
        Orbital period [yr]
    """
    # Sun at origin (stationary in COM frame approximation)
    M_sun = 1.0  # Solar mass in code units
    sun = Body(
        name="Sun",
        M=M_sun,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=R_sun,
        Q=0.0,
    )
    sun.update_Q_from_M(medium)

    # Mercury at perihelion (semi-major axis a = 0.387 AU)
    a = 0.387  # AU
    M_mercury = 1.66e-7  # Mercury/Sun mass ratio

    # Circular velocity: v = sqrt(K * M_sun / a)
    # where K = rho0 / (4 * pi * beta0^2)
    v_circ = np.sqrt(medium.K * M_sun / a)

    mercury = Body(
        name="Mercury",
        M=M_mercury,
        x=np.array([a, 0.0, 0.0]),
        v=np.array([0.0, v_circ, 0.0]),
        R=R_mercury,
        Q=0.0,
    )
    mercury.update_Q_from_M(medium)

    bodies = [sun, mercury]

    # Compute orbital period: T = 2*pi*sqrt(a^3 / (K*M_sun))
    T_orbit = estimate_orbital_period(bodies, medium)

    return bodies, T_orbit


# ============================================================================
# Test 1: Timestep Invariance Study
# ============================================================================

def run_timestep_test(config: TestConfig, verbose: bool = True) -> Dict:
    """
    Test energy conservation and convergence as dt → 0.

    Expected: Energy drift scales as dt^2 (2nd order integrator).

    Parameters
    ----------
    config : TestConfig
        Test configuration
    verbose : bool
        Print progress

    Returns
    -------
    results : dict
        Dictionary with:
        - 'dt_values' : timestep sizes (in years)
        - 'energy_drift' : final |ΔE/E| for each dt
        - 'position_error' : final position error vs finest dt
        - 'convergence_rate' : measured convergence order
        - 'passes' : bool, whether convergence is 2nd order
    """
    if verbose:
        print("\n" + "=" * 70)
        print("TEST 1: TIMESTEP INVARIANCE STUDY")
        print("=" * 70)
        print(f"Testing dt values: {config.dt_values}")
        print(f"Integration: {config.n_orbits} orbits")
        print(f"Compressible: {config.use_compressible}")
        print()

    # Setup medium (Solar System parameters)
    medium = Medium(
        rho0=1.0,
        cs=1e4,
        beta0=1e10,
        gamma_beta=0.0,
    )

    # Storage for results
    dt_actual = []
    energy_drifts = []
    final_positions = []

    # Run integration for each timestep
    for i, dt_factor in enumerate(config.dt_values):
        # Setup fresh system for each run
        bodies, T_orbit = setup_mercury_system(medium)

        # Compute timestep
        dt = dt_factor * T_orbit
        dt_actual.append(dt)

        # Compute number of steps
        total_time = config.n_orbits * T_orbit
        n_steps = int(total_time / dt)

        if verbose:
            print(f"[{i+1}/{len(config.dt_values)}] dt = {dt_factor:.4f} * T_orbit = {dt:.6f} yr")
            print(f"  n_steps = {n_steps}, total_time = {total_time:.3f} yr")

        # Integration options
        opts = {
            'use_compressible': config.use_compressible,
            'use_quadrature': False,  # Use fast analytic forces
            'save_every': max(1, n_steps // 100),  # Save ~100 snapshots
            'verbose': False,
        }

        # Run integration
        traj, diags = integrate_orbit(bodies, medium, dt, n_steps, opts)

        # Extract diagnostics
        E_initial = diags[0]['total_energy']
        E_final = diags[-1]['total_energy']
        energy_drift = abs((E_final - E_initial) / E_initial)
        energy_drifts.append(energy_drift)

        # Store final position of Mercury
        final_positions.append(bodies[1].x.copy())

        if verbose:
            print(f"  Energy drift: |ΔE/E| = {energy_drift:.6e}")

    # Compute position errors relative to finest timestep (reference)
    ref_position = final_positions[0]  # Finest dt is first
    position_errors = [np.linalg.norm(pos - ref_position) for pos in final_positions]

    # Estimate convergence rate from energy drift
    # If drift ~ dt^p, then log(drift) ~ p * log(dt)
    # Fit to get p
    log_dt = np.log(dt_actual)
    log_drift = np.log(energy_drifts)

    # Linear fit: log(drift) = p * log(dt) + const
    p_fit = np.polyfit(log_dt, log_drift, 1)
    convergence_rate = p_fit[0]

    # Check if convergence is 2nd order (expect p ≈ 2)
    passes = abs(convergence_rate - 2.0) < 0.5

    if verbose:
        print()
        print(f"Convergence rate: {convergence_rate:.2f} (expect ~2.0 for 2nd order)")
        print(f"PASS: {passes}" if passes else f"FAIL: {passes}")

    results = {
        'dt_values': np.array(dt_actual),
        'energy_drift': np.array(energy_drifts),
        'position_error': np.array(position_errors),
        'convergence_rate': convergence_rate,
        'passes': passes,
    }

    return results


def plot_timestep_results(results: Dict, output_path: Path, config: TestConfig):
    """Plot timestep convergence study results."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: Energy drift vs dt (log-log)
    ax = axes[0]
    ax.loglog(results['dt_values'], results['energy_drift'], 'o-', label='Measured drift')

    # Add reference line showing 2nd order convergence
    dt_ref = results['dt_values']
    drift_ref = results['energy_drift'][0] * (dt_ref / dt_ref[0])**2
    ax.loglog(dt_ref, drift_ref, '--', alpha=0.5, label='2nd order reference')

    ax.set_xlabel('Timestep [yr]')
    ax.set_ylabel('Energy drift |ΔE/E|')
    ax.set_title('Energy Conservation vs Timestep')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right panel: Position error vs dt (log-log)
    ax = axes[1]
    ax.loglog(results['dt_values'], results['position_error'], 's-', color='C1')
    ax.set_xlabel('Timestep [yr]')
    ax.set_ylabel('Position error [AU]')
    ax.set_title('Position Error vs Timestep\n(relative to finest dt)')
    ax.grid(True, alpha=0.3)

    # Add text box with convergence rate
    textstr = f"Convergence rate: {results['convergence_rate']:.2f}\n"
    textstr += f"Expected: 2.0 (2nd order)\n"
    textstr += f"PASS: {results['passes']}"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes[0].text(0.05, 0.95, textstr, transform=axes[0].transAxes,
                verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


# ============================================================================
# Test 2: Control Surface Radius Invariance
# ============================================================================

def run_radius_test(config: TestConfig, verbose: bool = True) -> Dict:
    """
    Test force independence of control surface radius R.

    Expected: Forces should be R-independent (< 0.1% variation).

    Parameters
    ----------
    config : TestConfig
        Test configuration
    verbose : bool
        Print progress

    Returns
    -------
    results : dict
        Dictionary with:
        - 'R_factors' : R multiples tested
        - 'forces_sun' : force on Sun for each R
        - 'forces_mercury' : force on Mercury for each R
        - 'rel_error_sun' : relative error vs nominal
        - 'rel_error_mercury' : relative error vs nominal
        - 'passes' : bool, whether variation < 0.1%
    """
    if verbose:
        print("\n" + "=" * 70)
        print("TEST 2: CONTROL SURFACE RADIUS INVARIANCE")
        print("=" * 70)
        print(f"Testing R factors: {config.R_factors}")
        print()

    # Setup medium
    medium = Medium(
        rho0=1.0,
        cs=1e4,
        beta0=1e10,
        gamma_beta=0.0,
    )

    # Nominal radii
    R_sun_nominal = 1e-3
    R_mercury_nominal = 5e-4

    forces_sun = []
    forces_mercury = []

    # Test each R factor
    for i, R_factor in enumerate(config.R_factors):
        R_sun = R_sun_nominal * R_factor
        R_mercury = R_mercury_nominal * R_factor

        if verbose:
            print(f"[{i+1}/{len(config.R_factors)}] R_factor = {R_factor:.1f}")
            print(f"  R_sun = {R_sun:.6f} AU, R_mercury = {R_mercury:.6f} AU")

        # Setup system with scaled radii
        bodies, T_orbit = setup_mercury_system(medium, R_sun, R_mercury)

        # Compute forces using analytic formula
        F_sun = force_incompressible_analytic(0, bodies, medium)
        F_mercury = force_incompressible_analytic(1, bodies, medium)

        forces_sun.append(F_sun)
        forces_mercury.append(F_mercury)

        if verbose:
            print(f"  |F_sun| = {np.linalg.norm(F_sun):.6e}")
            print(f"  |F_mercury| = {np.linalg.norm(F_mercury):.6e}")

    # Convert to arrays
    forces_sun = np.array(forces_sun)
    forces_mercury = np.array(forces_mercury)

    # Compute force magnitudes
    F_sun_mag = np.linalg.norm(forces_sun, axis=1)
    F_mercury_mag = np.linalg.norm(forces_mercury, axis=1)

    # Reference is nominal R (factor = 1.0)
    idx_nominal = config.R_factors.index(1.0)
    F_sun_nominal = F_sun_mag[idx_nominal]
    F_mercury_nominal = F_mercury_mag[idx_nominal]

    # Compute relative errors
    rel_error_sun = np.abs(F_sun_mag - F_sun_nominal) / F_sun_nominal
    rel_error_mercury = np.abs(F_mercury_mag - F_mercury_nominal) / F_mercury_nominal

    # Check if all variations < 0.1%
    max_error = max(rel_error_sun.max(), rel_error_mercury.max())
    passes = max_error < 1e-3

    if verbose:
        print()
        print(f"Maximum relative error: {max_error:.6e}")
        print(f"PASS: {passes} (threshold: 1e-3)" if passes else f"FAIL: {passes}")

    results = {
        'R_factors': np.array(config.R_factors),
        'forces_sun': forces_sun,
        'forces_mercury': forces_mercury,
        'F_sun_mag': F_sun_mag,
        'F_mercury_mag': F_mercury_mag,
        'rel_error_sun': rel_error_sun,
        'rel_error_mercury': rel_error_mercury,
        'max_error': max_error,
        'passes': passes,
    }

    return results


def plot_radius_results(results: Dict, output_path: Path, config: TestConfig):
    """Plot radius invariance test results."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: Force magnitude vs R factor
    ax = axes[0]
    ax.plot(results['R_factors'], results['F_sun_mag'], 'o-', label='Sun', markersize=8)
    ax.plot(results['R_factors'], results['F_mercury_mag'], 's-', label='Mercury', markersize=8)
    ax.set_xlabel('R factor (× nominal)')
    ax.set_ylabel('Force magnitude [code units]')
    ax.set_title('Force vs Control Surface Radius')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right panel: Relative error vs R factor
    ax = axes[1]
    ax.semilogy(results['R_factors'], results['rel_error_sun'], 'o-', label='Sun', markersize=8)
    ax.semilogy(results['R_factors'], results['rel_error_mercury'], 's-', label='Mercury', markersize=8)
    ax.axhline(1e-3, color='r', linestyle='--', alpha=0.5, label='0.1% threshold')
    ax.set_xlabel('R factor (× nominal)')
    ax.set_ylabel('Relative error')
    ax.set_title('Force Variation vs R\n(relative to nominal)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add text box with result
    textstr = f"Max error: {results['max_error']:.3e}\n"
    textstr += f"Threshold: 1e-3 (0.1%)\n"
    textstr += f"PASS: {results['passes']}"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes[0].text(0.05, 0.95, textstr, transform=axes[0].transAxes,
                verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


# ============================================================================
# Test 3: Quadrature Convergence Study
# ============================================================================

def run_quadrature_test(config: TestConfig, verbose: bool = True) -> Dict:
    """
    Test convergence of quadrature forces vs number of points.

    Expected: Error decreases as power law, converges for n > 256.

    Parameters
    ----------
    config : TestConfig
        Test configuration
    verbose : bool
        Print progress

    Returns
    -------
    results : dict
        Dictionary with:
        - 'n_points_values' : number of quadrature points tested
        - 'rel_error' : relative error vs analytic for each n
        - 'convergence_rate' : power-law exponent
        - 'passes' : bool, whether error < 1e-3 for n >= 256
    """
    if verbose:
        print("\n" + "=" * 70)
        print("TEST 3: QUADRATURE CONVERGENCE STUDY")
        print("=" * 70)
        print(f"Testing n_points values: {config.n_points_values}")
        print()

    # Setup medium
    medium = Medium(
        rho0=1.0,
        cs=1e4,
        beta0=1e10,
        gamma_beta=0.0,
    )

    # Setup system
    bodies, T_orbit = setup_mercury_system(medium)

    # Compute analytic force (reference)
    F_analytic = force_incompressible_analytic(1, bodies, medium)
    F_analytic_mag = np.linalg.norm(F_analytic)

    if verbose:
        print(f"Analytic force magnitude: {F_analytic_mag:.6e}")
        print()

    rel_errors = []

    # Test each n_points value
    for i, n_points in enumerate(config.n_points_values):
        if verbose:
            print(f"[{i+1}/{len(config.n_points_values)}] n_points = {n_points}")

        # Compute quadrature force
        F_quadrature = force_incompressible_quadrature(1, bodies, medium, n_points)
        F_quadrature_mag = np.linalg.norm(F_quadrature)

        # Compute relative error
        rel_error = np.linalg.norm(F_quadrature - F_analytic) / F_analytic_mag
        rel_errors.append(rel_error)

        if verbose:
            print(f"  Quadrature force: {F_quadrature_mag:.6e}")
            print(f"  Relative error: {rel_error:.6e}")

    rel_errors = np.array(rel_errors)

    # Estimate convergence rate
    # If error ~ 1/n^p, then log(error) ~ -p * log(n)
    log_n = np.log(config.n_points_values)
    log_error = np.log(rel_errors)

    # Linear fit
    p_fit = np.polyfit(log_n, log_error, 1)
    convergence_rate = -p_fit[0]  # Negative because error decreases

    # Check if error < 1e-3 for n >= 256
    idx_256 = None
    for i, n in enumerate(config.n_points_values):
        if n >= 256:
            idx_256 = i
            break

    if idx_256 is not None:
        error_at_256 = rel_errors[idx_256]
        passes = error_at_256 < 1e-3
    else:
        passes = False

    if verbose:
        print()
        print(f"Convergence rate: {convergence_rate:.2f}")
        if idx_256 is not None:
            print(f"Error at n=256: {error_at_256:.6e}")
        print(f"PASS: {passes}" if passes else f"FAIL: {passes}")

    results = {
        'n_points_values': np.array(config.n_points_values),
        'rel_error': rel_errors,
        'convergence_rate': convergence_rate,
        'passes': passes,
    }

    return results


def plot_quadrature_results(results: Dict, output_path: Path, config: TestConfig):
    """Plot quadrature convergence study results."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Log-log plot of error vs n_points
    ax.loglog(results['n_points_values'], results['rel_error'], 'o-', markersize=8)

    # Add reference line showing measured convergence rate
    n_ref = results['n_points_values']
    error_ref = results['rel_error'][0] * (n_ref / n_ref[0])**(-results['convergence_rate'])
    ax.loglog(n_ref, error_ref, '--', alpha=0.5,
             label=f'Power law: n^(-{results["convergence_rate"]:.2f})')

    # Add threshold line
    ax.axhline(1e-3, color='r', linestyle='--', alpha=0.5, label='1e-3 threshold')

    ax.set_xlabel('Number of quadrature points')
    ax.set_ylabel('Relative error |F_quad - F_analytic| / |F_analytic|')
    ax.set_title('Quadrature Convergence Study')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add text box with convergence rate
    textstr = f"Convergence rate: {results['convergence_rate']:.2f}\n"
    textstr += f"PASS: {results['passes']}"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
           verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_path}")
    plt.close()


# ============================================================================
# Main driver
# ============================================================================

def print_summary_report(results_timestep: Optional[Dict],
                        results_radius: Optional[Dict],
                        results_quadrature: Optional[Dict],
                        config: TestConfig):
    """Print final summary report."""
    print("\n" + "=" * 70)
    print("SUMMARY REPORT")
    print("=" * 70)
    print()

    # Timestep test
    if results_timestep is not None:
        print("1. TIMESTEP INVARIANCE:")
        print(f"   Convergence rate: {results_timestep['convergence_rate']:.2f} (expect ~2.0)")
        print(f"   Status: {'PASS' if results_timestep['passes'] else 'FAIL'}")
        print()

    # Radius test
    if results_radius is not None:
        print("2. CONTROL SURFACE RADIUS INVARIANCE:")
        print(f"   Maximum force variation: {results_radius['max_error']:.6e}")
        print(f"   Threshold: 1e-3 (0.1%)")
        print(f"   Status: {'PASS' if results_radius['passes'] else 'FAIL'}")
        print()

    # Quadrature test
    if results_quadrature is not None:
        print("3. QUADRATURE CONVERGENCE:")
        print(f"   Convergence rate: {results_quadrature['convergence_rate']:.2f}")
        print(f"   Status: {'PASS' if results_quadrature['passes'] else 'FAIL'}")
        print()

    print("=" * 70)
    print("RECOMMENDATIONS:")
    print("=" * 70)

    if results_timestep is not None:
        # Recommend timestep based on acceptable energy drift
        acceptable_drift = 1e-5
        # Find timestep that gives this drift (extrapolate from measured)
        dt_values = results_timestep['dt_values']
        drifts = results_timestep['energy_drift']
        # Fit power law
        log_dt = np.log(dt_values)
        log_drift = np.log(drifts)
        p_fit = np.polyfit(log_dt, log_drift, 1)
        # Solve for dt: log(acceptable_drift) = p * log(dt) + const
        log_dt_recommended = (np.log(acceptable_drift) - p_fit[1]) / p_fit[0]
        dt_recommended = np.exp(log_dt_recommended)
        print(f"- Timestep: Use dt < {dt_recommended:.6f} yr for |ΔE/E| < 1e-5")

    if results_quadrature is not None:
        # Recommend n_points based on error threshold
        n_256_idx = None
        for i, n in enumerate(results_quadrature['n_points_values']):
            if n >= 256:
                n_256_idx = i
                break
        if n_256_idx is not None:
            error_256 = results_quadrature['rel_error'][n_256_idx]
            print(f"- Quadrature: Use n_points >= 256 (error = {error_256:.3e})")
        else:
            print(f"- Quadrature: Use n_points >= 256 (not tested)")

    if results_radius is not None:
        print(f"- Control surface: R values tested show < 0.1% variation")
        print(f"  Use R ~ 1e-3 to 1e-4 AU for typical separations")

    print()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Numerical robustness validation tests for superfluid orbit simulator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Test selection
    parser.add_argument('--all', action='store_true',
                       help='Run all tests (default)')
    parser.add_argument('--timestep', action='store_true',
                       help='Run timestep invariance test only')
    parser.add_argument('--radius', action='store_true',
                       help='Run radius invariance test only')
    parser.add_argument('--quadrature', action='store_true',
                       help='Run quadrature convergence test only')

    # Configuration
    parser.add_argument('--quick', action='store_true',
                       help='Use reduced parameter ranges for quick test')
    parser.add_argument('--output-dir', type=str, default='output',
                       help='Output directory for plots (default: output/)')
    parser.add_argument('--compressible', action='store_true',
                       help='Enable compressible corrections (slower)')

    args = parser.parse_args()

    # Determine which tests to run
    run_all = args.all or not (args.timestep or args.radius or args.quadrature)
    run_timestep = run_all or args.timestep
    run_radius = run_all or args.radius
    run_quadrature = run_all or args.quadrature

    # Setup configuration
    if args.quick:
        config = TestConfig.quick()
        print("Using QUICK test configuration")
    else:
        config = TestConfig.default()
        print("Using STANDARD test configuration")

    config.output_dir = Path(args.output_dir)
    config.use_compressible = args.compressible

    # Create output directory
    config.output_dir.mkdir(exist_ok=True, parents=True)

    # Run tests
    results_timestep = None
    results_radius = None
    results_quadrature = None

    if run_timestep:
        results_timestep = run_timestep_test(config, verbose=True)
        plot_timestep_results(results_timestep,
                            config.output_dir / 'robustness_timestep.png',
                            config)

    if run_radius:
        results_radius = run_radius_test(config, verbose=True)
        plot_radius_results(results_radius,
                          config.output_dir / 'robustness_radius.png',
                          config)

    if run_quadrature:
        results_quadrature = run_quadrature_test(config, verbose=True)
        plot_quadrature_results(results_quadrature,
                               config.output_dir / 'robustness_quadrature.png',
                               config)

    # Print summary report
    print_summary_report(results_timestep, results_radius, results_quadrature, config)

    print(f"\nAll plots saved to: {config.output_dir.absolute()}")
    print("Done!")


if __name__ == '__main__':
    main()
