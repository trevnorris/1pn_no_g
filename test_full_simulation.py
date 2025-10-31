#!/usr/bin/env python3
"""
Comprehensive integration test for the superfluid orbit simulator.

This test verifies that the complete simulation pipeline works correctly
after the bug fix. It creates a simple 2-body system and checks:
1. Energy conservation (< 10^-5 drift for symplectic integrator)
2. Momentum conservation (center of mass stationary)
3. Force law accuracy (1/r² scaling with correct coefficient)
4. Orbit stability (remains approximately circular)
5. No numerical errors (NaN, inf)
6. All modules integrate correctly (medium, bodies, dynamics, diagnostics)

Test philosophy:
- Simple, deterministic setup (no randomness)
- Clean parameters for easy verification
- Multiple orbits to check long-term stability
- Clear pass/fail criteria
- Comprehensive diagnostics table
"""

import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit, assemble_forces
from slab.diagnostics import (
    total_energy,
    total_kinetic_energy,
    fluid_potential_energy,
    angular_momentum,
)


def create_simple_circular_orbit():
    """Create a simple 2-body circular orbit for testing.

    Returns
    -------
    bodies : List[Body]
        Two-body system (Sun + test body)
    medium : Medium
        Ambient superfluid medium
    orbital_params : dict
        Expected orbital parameters for verification
    """
    # Medium with simple parameters
    medium = Medium(
        rho0=1.0,
        cs=1e4,
        beta0=1e10,
        gamma_beta=0.0  # No rarefaction
    )

    # Two bodies: heavy primary + light secondary (like Sun-Earth)
    M_sun = 1.0
    M_planet = 3e-6  # ~Earth mass ratio
    a = 1.0  # Semi-major axis (separation)

    # Circular orbit velocity: v = sqrt(K*M_sun/a)
    # where K = rho0/(4*pi*beta0^2)
    K = medium.K
    v_circ = np.sqrt(K * M_sun / a)

    # Create bodies
    sun = Body(
        name="Sun",
        M=M_sun,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=1e-3,  # Control surface radius
        Q=0.0
    )
    sun.update_Q_from_M(medium)

    planet = Body(
        name="Planet",
        M=M_planet,
        x=np.array([a, 0.0, 0.0]),  # Start on x-axis
        v=np.array([0.0, v_circ, 0.0]),  # Velocity in y-direction
        R=5e-4,
        Q=0.0
    )
    planet.update_Q_from_M(medium)

    bodies = [sun, planet]

    # Calculate expected orbital parameters
    # Period: T = 2*pi*sqrt(a^3 / (K*M_sun))
    T_orbit = 2.0 * np.pi * np.sqrt(a**3 / (K * M_sun))

    # Expected energy (circular orbit): E = -K*M1*M2/(2*a)
    E_expected = -K * M_sun * M_planet / (2.0 * a)

    # Expected angular momentum magnitude: L = M_planet * v_circ * a
    L_expected = M_planet * v_circ * a

    # Expected force magnitude on planet: F = K*M1*M2/a^2
    F_expected = K * M_sun * M_planet / (a**2)

    orbital_params = {
        'a': a,
        'K': K,
        'v_circ': v_circ,
        'T_orbit': T_orbit,
        'E_expected': E_expected,
        'L_expected': L_expected,
        'F_expected': F_expected,
    }

    return bodies, medium, orbital_params


def run_integration(bodies, medium, orbital_params, n_orbits=3):
    """Run the integration for several orbits.

    Parameters
    ----------
    bodies : List[Body]
        Bodies to integrate
    medium : Medium
        Superfluid medium
    orbital_params : dict
        Expected orbital parameters
    n_orbits : int
        Number of orbits to simulate

    Returns
    -------
    trajectory : dict
        Trajectory data
    diagnostics : List[dict]
        Diagnostics at each saved step
    """
    T_orbit = orbital_params['T_orbit']

    # Choose timestep: dt = 0.005 * T (1/200 of orbital period)
    # This should give excellent energy conservation
    dt = 0.005 * T_orbit

    # Total integration time
    t_total = n_orbits * T_orbit
    n_steps = int(t_total / dt)

    print(f"Integration setup:")
    print(f"  Orbital period: T = {T_orbit:.6e}")
    print(f"  Timestep: dt = {dt:.6e} ({dt/T_orbit:.4f} T)")
    print(f"  Number of orbits: {n_orbits}")
    print(f"  Total steps: {n_steps}")
    print(f"  Total time: {t_total:.6e}")
    print()

    # Integration options
    opts = {
        'use_compressible': False,  # Incompressible (simplest)
        'use_quadrature': False,    # Analytic forces (fast)
        'save_every': max(1, n_steps // 100),  # ~100 snapshots
        'verbose': False,
    }

    # Run integration
    print("Running integration...")
    trajectory, diagnostics = integrate_orbit(bodies, medium, dt, n_steps, opts)
    print("Integration complete!")
    print()

    return trajectory, diagnostics


def check_energy_conservation(diagnostics, orbital_params, tolerance=1e-5):
    """Check energy conservation over the trajectory.

    Parameters
    ----------
    diagnostics : List[dict]
        Diagnostics from integration
    orbital_params : dict
        Expected orbital parameters
    tolerance : float
        Relative energy drift tolerance

    Returns
    -------
    result : dict
        Test results
    """
    energies = np.array([d['total_energy'] for d in diagnostics])
    E0 = energies[0]
    E_final = energies[-1]

    # Energy drift
    dE = energies - E0
    rel_dE = dE / np.abs(E0)
    max_drift = np.max(np.abs(rel_dE))
    rms_drift = np.sqrt(np.mean(rel_dE**2))

    # Check against expected energy
    E_expected = orbital_params['E_expected']
    E_error = np.abs(E0 - E_expected) / np.abs(E_expected)

    # Pass/fail
    passed = max_drift < tolerance

    return {
        'E0': E0,
        'E_final': E_final,
        'E_expected': E_expected,
        'E_error': E_error,
        'max_drift': max_drift,
        'rms_drift': rms_drift,
        'tolerance': tolerance,
        'passed': passed,
    }


def check_momentum_conservation(diagnostics, tolerance=1e-12):
    """Check total momentum conservation.

    For an isolated system, total momentum should be exactly conserved.

    Parameters
    ----------
    diagnostics : List[dict]
        Diagnostics from integration
    tolerance : float
        Absolute momentum drift tolerance

    Returns
    -------
    result : dict
        Test results
    """
    p_total = np.array([d['total_momentum'] for d in diagnostics])
    p0 = p_total[0]

    # Momentum drift
    dp = np.linalg.norm(p_total - p0[None, :], axis=1)
    max_drift = np.max(dp)

    # For circular orbit starting from rest in CM frame, p0 should be ~0
    p0_mag = np.linalg.norm(p0)

    passed = max_drift < tolerance

    return {
        'p0': p0,
        'p0_magnitude': p0_mag,
        'max_drift': max_drift,
        'tolerance': tolerance,
        'passed': passed,
    }


def check_force_law(bodies, medium, orbital_params, tolerance=1e-10):
    """Check that forces obey the 1/r² law with correct coefficient.

    For two bodies at separation r:
        F = K * M1 * M2 / r²

    where K = rho0/(4*pi*beta0^2).

    Parameters
    ----------
    bodies : List[Body]
        Current body states
    medium : Medium
        Superfluid medium
    orbital_params : dict
        Expected orbital parameters
    tolerance : float
        Relative error tolerance

    Returns
    -------
    result : dict
        Test results
    """
    # Compute forces
    forces = assemble_forces(bodies, medium)

    # Force on planet (body 1)
    F_planet = forces[1]
    F_mag = np.linalg.norm(F_planet)

    # Expected force
    F_expected = orbital_params['F_expected']

    # Relative error
    rel_error = np.abs(F_mag - F_expected) / F_expected

    # Check 1/r² scaling: compute force at different separation
    # Move planet to 2*a and check force drops by 1/4
    original_x = bodies[1].x.copy()
    bodies[1].x = 2.0 * original_x
    forces_2a = assemble_forces(bodies, medium)
    F_2a_mag = np.linalg.norm(forces_2a[1])
    bodies[1].x = original_x  # Restore

    # Should have F(2a) = F(a) / 4
    scaling_ratio = F_mag / F_2a_mag
    scaling_expected = 4.0
    scaling_error = np.abs(scaling_ratio - scaling_expected) / scaling_expected

    # Check Newton's third law: F1 = -F2
    F_sum = forces[0] + forces[1]
    newton3_error = np.linalg.norm(F_sum)

    passed = (rel_error < tolerance and
              scaling_error < tolerance and
              newton3_error < 1e-14)

    return {
        'F_measured': F_mag,
        'F_expected': F_expected,
        'rel_error': rel_error,
        'scaling_ratio': scaling_ratio,
        'scaling_expected': scaling_expected,
        'scaling_error': scaling_error,
        'newton3_error': newton3_error,
        'tolerance': tolerance,
        'passed': passed,
    }


def check_orbit_circularity(trajectory, orbital_params, tolerance=0.01):
    """Check that orbit remains approximately circular.

    For a circular orbit, the distance from primary should be constant.

    Parameters
    ----------
    trajectory : dict
        Trajectory data
    orbital_params : dict
        Expected orbital parameters
    tolerance : float
        Relative variation tolerance (e.g., 0.01 = 1%)

    Returns
    -------
    result : dict
        Test results
    """
    # Extract positions
    positions = trajectory['x']  # Shape: (n_saved, 2, 3)

    # Compute separation between bodies
    separations = np.linalg.norm(positions[:, 1, :] - positions[:, 0, :], axis=1)

    # Expected separation (semi-major axis)
    a = orbital_params['a']

    # Statistics
    r_mean = np.mean(separations)
    r_std = np.std(separations)
    r_min = np.min(separations)
    r_max = np.max(separations)

    # Eccentricity estimate: e ≈ (r_max - r_min) / (r_max + r_min)
    e_estimate = (r_max - r_min) / (r_max + r_min)

    # Relative variation
    rel_variation = r_std / r_mean

    passed = rel_variation < tolerance and e_estimate < tolerance

    return {
        'a_expected': a,
        'r_mean': r_mean,
        'r_std': r_std,
        'r_min': r_min,
        'r_max': r_max,
        'eccentricity': e_estimate,
        'rel_variation': rel_variation,
        'tolerance': tolerance,
        'passed': passed,
    }


def check_no_numerical_errors(trajectory, diagnostics):
    """Check for NaN or inf values in trajectory and diagnostics.

    Parameters
    ----------
    trajectory : dict
        Trajectory data
    diagnostics : List[dict]
        Diagnostics from integration

    Returns
    -------
    result : dict
        Test results
    """
    # Check trajectory arrays
    positions = trajectory['x']
    velocities = trajectory['v']

    has_nan_pos = np.any(np.isnan(positions))
    has_inf_pos = np.any(np.isinf(positions))
    has_nan_vel = np.any(np.isnan(velocities))
    has_inf_vel = np.any(np.isinf(velocities))

    # Check diagnostics
    energies = np.array([d['total_energy'] for d in diagnostics])
    has_nan_energy = np.any(np.isnan(energies))
    has_inf_energy = np.any(np.isinf(energies))

    passed = not (has_nan_pos or has_inf_pos or
                  has_nan_vel or has_inf_vel or
                  has_nan_energy or has_inf_energy)

    return {
        'has_nan': has_nan_pos or has_nan_vel or has_nan_energy,
        'has_inf': has_inf_pos or has_inf_vel or has_inf_energy,
        'passed': passed,
    }


def print_results_table(results):
    """Print a formatted summary table of test results.

    Parameters
    ----------
    results : dict
        Dictionary of test results
    """
    print("=" * 80)
    print("TEST RESULTS SUMMARY")
    print("=" * 80)
    print()

    # Overall status
    all_passed = all(r['passed'] for r in results.values())

    if all_passed:
        print("    ███████╗██╗   ██╗ ██████╗ ██████╗███████╗███████╗███████╗")
        print("    ██╔════╝██║   ██║██╔════╝██╔════╝██╔════╝██╔════╝██╔════╝")
        print("    ███████╗██║   ██║██║     ██║     █████╗  ███████╗███████╗")
        print("    ╚════██║██║   ██║██║     ██║     ██╔══╝  ╚════██║╚════██║")
        print("    ███████║╚██████╔╝╚██████╗╚██████╗███████╗███████║███████║")
        print("    ╚══════╝ ╚═════╝  ╚═════╝ ╚═════╝╚══════╝╚══════╝╚══════╝")
        print()
        print("                 ALL TESTS PASSED!")
        print()
    else:
        print("    ⚠  SOME TESTS FAILED - SEE DETAILS BELOW  ⚠")
        print()

    # Energy conservation
    print("-" * 80)
    print("1. ENERGY CONSERVATION")
    print("-" * 80)
    r = results['energy']
    status = "✓ PASS" if r['passed'] else "✗ FAIL"
    print(f"  Status: {status}")
    print(f"  Initial energy:       E0 = {r['E0']:+.6e}")
    print(f"  Expected energy:  E_exp = {r['E_expected']:+.6e}")
    print(f"  Initial error:        ε = {r['E_error']:.2e}")
    print(f"  Maximum drift:  |ΔE/E| = {r['max_drift']:.2e}  (tolerance: {r['tolerance']:.0e})")
    print(f"  RMS drift:      |ΔE/E| = {r['rms_drift']:.2e}")
    print()

    # Momentum conservation
    print("-" * 80)
    print("2. MOMENTUM CONSERVATION")
    print("-" * 80)
    r = results['momentum']
    status = "✓ PASS" if r['passed'] else "✗ FAIL"
    print(f"  Status: {status}")
    print(f"  Initial momentum: p0 = [{r['p0'][0]:+.2e}, {r['p0'][1]:+.2e}, {r['p0'][2]:+.2e}]")
    print(f"  |p0|:                  {r['p0_magnitude']:.2e}")
    print(f"  Maximum drift:   |Δp| = {r['max_drift']:.2e}  (tolerance: {r['tolerance']:.0e})")
    print()

    # Force law
    print("-" * 80)
    print("3. FORCE LAW (1/r² with correct coefficient)")
    print("-" * 80)
    r = results['force']
    status = "✓ PASS" if r['passed'] else "✗ FAIL"
    print(f"  Status: {status}")
    print(f"  Measured force:   F = {r['F_measured']:.6e}")
    print(f"  Expected force:   F = {r['F_expected']:.6e}")
    print(f"  Relative error:   ε = {r['rel_error']:.2e}  (tolerance: {r['tolerance']:.0e})")
    print(f"  1/r² scaling: F(a)/F(2a) = {r['scaling_ratio']:.6f}  (expected: {r['scaling_expected']:.1f})")
    print(f"  Scaling error:        ε = {r['scaling_error']:.2e}")
    print(f"  Newton 3rd law: |F1+F2| = {r['newton3_error']:.2e}")
    print()

    # Orbit circularity
    print("-" * 80)
    print("4. ORBIT STABILITY (circular orbit)")
    print("-" * 80)
    r = results['orbit']
    status = "✓ PASS" if r['passed'] else "✗ FAIL"
    print(f"  Status: {status}")
    print(f"  Expected radius:   a = {r['a_expected']:.6f}")
    print(f"  Mean radius:   r_avg = {r['r_mean']:.6f}")
    print(f"  Std deviation: r_std = {r['r_std']:.2e}")
    print(f"  Min radius:    r_min = {r['r_min']:.6f}")
    print(f"  Max radius:    r_max = {r['r_max']:.6f}")
    print(f"  Eccentricity:      e ≈ {r['eccentricity']:.2e}  (tolerance: {r['tolerance']:.2f})")
    print(f"  Rel. variation:    σ = {r['rel_variation']:.2e}  (tolerance: {r['tolerance']:.2f})")
    print()

    # Numerical errors
    print("-" * 80)
    print("5. NUMERICAL STABILITY (no NaN/inf)")
    print("-" * 80)
    r = results['numerical']
    status = "✓ PASS" if r['passed'] else "✗ FAIL"
    print(f"  Status: {status}")
    print(f"  Contains NaN: {r['has_nan']}")
    print(f"  Contains inf: {r['has_inf']}")
    print()

    print("=" * 80)
    if all_passed:
        print("OVERALL: ✓ ALL CHECKS PASSED")
    else:
        failed_tests = [name for name, r in results.items() if not r['passed']]
        print(f"OVERALL: ✗ FAILED TESTS: {', '.join(failed_tests)}")
    print("=" * 80)
    print()


def main():
    """Run the comprehensive integration test."""
    print()
    print("=" * 80)
    print("COMPREHENSIVE INTEGRATION TEST FOR SUPERFLUID ORBIT SIMULATOR")
    print("=" * 80)
    print()
    print("This test verifies the complete simulation pipeline:")
    print("  - Medium setup")
    print("  - Body creation and initialization")
    print("  - Force calculation (analytic)")
    print("  - Symplectic integration (Verlet)")
    print("  - Energy/momentum conservation")
    print("  - Long-term orbit stability")
    print()

    # Create test system
    print("-" * 80)
    print("STEP 1: Create simple 2-body circular orbit")
    print("-" * 80)
    bodies, medium, orbital_params = create_simple_circular_orbit()

    print("Medium parameters:")
    print(f"  rho0  = {medium.rho0:.3e}")
    print(f"  cs    = {medium.cs:.3e}")
    print(f"  beta0 = {medium.beta0:.3e}")
    print(f"  K     = {medium.K:.6e}  (orbital constant, replaces G)")
    print()

    print("Body setup:")
    for body in bodies:
        print(f"  {body.name:8s}: M = {body.M:.3e}, Q = {body.Q:.3e}, R = {body.R:.3e}")
    print()

    print("Orbital parameters:")
    print(f"  Semi-major axis: a = {orbital_params['a']:.6f}")
    print(f"  Circular velocity: v = {orbital_params['v_circ']:.6e}")
    print(f"  Orbital period: T = {orbital_params['T_orbit']:.6e}")
    print()

    # Run integration
    print("-" * 80)
    print("STEP 2: Run integration for 3 orbits")
    print("-" * 80)
    trajectory, diagnostics = run_integration(bodies, medium, orbital_params, n_orbits=3)

    # Perform checks
    print("-" * 80)
    print("STEP 3: Perform verification checks")
    print("-" * 80)
    print()

    results = {}

    print("Checking energy conservation...")
    results['energy'] = check_energy_conservation(diagnostics, orbital_params)

    print("Checking momentum conservation...")
    results['momentum'] = check_momentum_conservation(diagnostics)

    print("Checking force law...")
    results['force'] = check_force_law(bodies, medium, orbital_params)

    print("Checking orbit circularity...")
    results['orbit'] = check_orbit_circularity(trajectory, orbital_params)

    print("Checking for numerical errors...")
    results['numerical'] = check_no_numerical_errors(trajectory, diagnostics)

    print("All checks complete!")
    print()

    # Print results
    print_results_table(results)

    # Return status code
    all_passed = all(r['passed'] for r in results.values())
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
