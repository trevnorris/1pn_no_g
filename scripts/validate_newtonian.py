#!/usr/bin/env python3
"""
Newtonian Gravity Validation Framework
========================================

This script establishes baseline precision requirements for Newtonian gravity
before measuring tiny 1PN corrections. It compares three integration methods:

1. **Pure Newtonian gravity**: F = -G M₁ M₂ r̂ / r²
2. **Incompressible superfluid**: F = ρ₀ Q v_ext (use_compressible=False)
3. **GR 1PN with cs→∞**: Einstein-Infeld-Hoffmann equations (optional)

Key Question
------------
At what timestep does spurious precession become smaller than the GR signal?
- GR prediction: ~0.1 arcsec/orbit for Mercury
- Target: <0.01 arcsec/orbit (10% of signal)
- Stretch goal: <0.001 arcsec/orbit (1% of signal)

Context
-------
The project currently shows large spurious precession artifacts:
- dt = 0.002 yr: ~890 arcsec/orbit spurious precession
- dt = 1e-5 yr: ~0.02 arcsec/orbit spurious precession

This validation determines what minimum timestep is needed for reliable
measurements of 1PN effects.

Usage
-----
    python scripts/validate_newtonian.py
    python scripts/validate_newtonian.py --dt 0.001 0.0005 0.0001
    python scripts/validate_newtonian.py --n-orbits 10
    python scripts/validate_newtonian.py --quick  # Fast test with 3 timesteps
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import argparse
import numpy as np
import copy
from typing import Dict, List, Tuple
from dataclasses import dataclass

from slab.bodies import Body
from slab.medium import Medium
from slab.dynamics import integrate_orbit
from slab.gr1pn import integrate_gr1pn_orbit, compute_orbit_elements
from scripts.precession_helpers import (
    create_barycentric_mercury_config,
    analyse_precession,
    TwoBodyParams,
)


# ============================================================================
# Pure Newtonian Reference Implementation
# ============================================================================

def newtonian_accel(bodies: List[Body], G: float) -> np.ndarray:
    """
    Compute pure Newtonian accelerations: a_i = -Σ_j G M_j r_ij / |r_ij|³

    This is the textbook F = -G M₁ M₂ r̂ / r² force law, for comparison
    with the superfluid implementation.

    Parameters
    ----------
    bodies : List[Body]
        List of bodies with M, x, v attributes
    G : float
        Gravitational constant (in code units, replaces K)

    Returns
    -------
    accels : np.ndarray
        Accelerations for each body, shape (N, 3)
    """
    N = len(bodies)
    accels = np.zeros((N, 3), dtype=np.float64)

    for i in range(N):
        for j in range(N):
            if i == j:
                continue

            # Separation vector from j to i
            r_ij = bodies[i].x - bodies[j].x
            r_mag = np.linalg.norm(r_ij)

            if r_mag < 1e-15:
                raise ValueError(f"Bodies {i} and {j} have collided!")

            # Newtonian force: F = -G M_i M_j / r² * r_hat
            # Acceleration: a_i = F/M_i = -G M_j / r² * r_hat
            accels[i] -= G * bodies[j].M * r_ij / (r_mag ** 3)

    return accels


def integrate_newtonian_orbit(
    bodies: List[Body],
    G: float,
    dt: float,
    n_steps: int,
    save_every: int = 1,
) -> Dict[str, np.ndarray]:
    """
    Integrate N-body system using pure Newtonian gravity.

    Uses velocity-Verlet integrator (same as slab and GR integrators)
    for fair comparison.

    Parameters
    ----------
    bodies : List[Body]
        Initial conditions. Bodies are modified in-place.
    G : float
        Gravitational constant (replaces K from superfluid)
    dt : float
        Timestep size [yr]
    n_steps : int
        Number of integration steps
    save_every : int
        Save trajectory every N steps

    Returns
    -------
    trajectory : Dict[str, np.ndarray]
        't': times, shape (n_saved,)
        'x': positions, shape (n_saved, N, 3)
        'v': velocities, shape (n_saved, N, 3)
        'a': accelerations, shape (n_saved, N, 3)
    """
    N = len(bodies)
    n_saved = (n_steps // save_every) + 1

    trajectory = {
        't': np.zeros(n_saved),
        'x': np.zeros((n_saved, N, 3)),
        'v': np.zeros((n_saved, N, 3)),
        'a': np.zeros((n_saved, N, 3)),
    }

    # Initial state
    t = 0.0
    save_idx = 0

    # Compute initial accelerations
    a_current = newtonian_accel(bodies, G)

    # Save initial state
    trajectory['t'][save_idx] = t
    for i, body in enumerate(bodies):
        trajectory['x'][save_idx, i] = body.x.copy()
        trajectory['v'][save_idx, i] = body.v.copy()
        trajectory['a'][save_idx, i] = a_current[i].copy()
    save_idx += 1

    # Velocity-Verlet integration
    for step in range(1, n_steps + 1):
        # Half-step velocity
        for i, body in enumerate(bodies):
            body.v += 0.5 * a_current[i] * dt

        # Full-step position
        for body in bodies:
            body.x += body.v * dt

        # New accelerations
        a_new = newtonian_accel(bodies, G)

        # Complete velocity update
        for i, body in enumerate(bodies):
            body.v += 0.5 * a_new[i] * dt

        t += dt
        a_current = a_new

        # Save if requested
        if step % save_every == 0:
            trajectory['t'][save_idx] = t
            for i, body in enumerate(bodies):
                trajectory['x'][save_idx, i] = body.x.copy()
                trajectory['v'][save_idx, i] = body.v.copy()
                trajectory['a'][save_idx, i] = a_current[i].copy()
            save_idx += 1

    # Trim if needed
    if save_idx < n_saved:
        for key in trajectory:
            trajectory[key] = trajectory[key][:save_idx]

    return trajectory


# ============================================================================
# Convergence Analysis
# ============================================================================

@dataclass
class ConvergenceResult:
    """Results for a single timestep in convergence study."""
    dt: float
    steps_per_orbit: float
    method: str
    precession_arcsec: float
    precession_rad: float
    precession_std_arcsec: float
    n_periapsis: int
    energy_drift: float
    a_drift: float
    e_drift: float
    integration_time: float


def run_convergence_study(
    dt_values: List[float],
    n_orbits: float,
    eccentricity: float = 0.1,
    include_gr1pn: bool = False,
    verbose: bool = True,
) -> Tuple[List[ConvergenceResult], TwoBodyParams]:
    """
    Run timestep convergence study comparing integration methods.

    Parameters
    ----------
    dt_values : List[float]
        Timestep values to test [yr]
    n_orbits : float
        Number of orbits to integrate
    eccentricity : float
        Orbital eccentricity
    include_gr1pn : bool
        Include GR 1PN comparison (slower)
    verbose : bool
        Print progress updates

    Returns
    -------
    results : List[ConvergenceResult]
        Convergence results for all methods and timesteps
    params : TwoBodyParams
        Orbital parameters
    """
    # Set up Mercury-Sun barycentric system
    rho0 = 1.0
    beta0 = 0.045
    cs = 63239.7263  # Large cs for "Newtonian" limit
    a = 0.387  # AU

    medium, bodies_template, params = create_barycentric_mercury_config(
        eccentricity=eccentricity,
        rho0=rho0,
        beta0=beta0,
        cs=cs,
        semi_major_axis=a,
    )

    K = medium.K
    T_orbit = params.T_orbit

    if verbose:
        print("=" * 80)
        print("NEWTONIAN GRAVITY VALIDATION FRAMEWORK")
        print("=" * 80)
        print()
        print("Orbital Parameters:")
        print(f"  Semi-major axis:  a = {a:.4f} AU")
        print(f"  Eccentricity:     e = {eccentricity:.3f}")
        print(f"  Orbital period:   T = {T_orbit:.4f} yr")
        print(f"  Total mass:       M = {params.M_total:.3e}")
        print(f"  Orbital constant: K = {K:.3e}")
        print()
        print(f"Integration: {n_orbits:.1f} orbits")
        print()
        print("Methods:")
        print("  1. Pure Newtonian (F = -GMm/r²)")
        print("  2. Incompressible superfluid (F = ρ₀Qv, use_compressible=False)")
        if include_gr1pn:
            print("  3. GR 1PN with large cs (Einstein-Infeld-Hoffmann)")
        print()
        print("=" * 80)
        print()

    results = []

    for dt in dt_values:
        n_steps = int(n_orbits * T_orbit / dt)
        steps_per_orbit = T_orbit / dt
        save_every = max(1, n_steps // 2000)  # Save ~2000 points max

        if verbose:
            print(f"\nTimestep: dt = {dt:.6f} yr")
            print(f"  Steps: {n_steps} ({steps_per_orbit:.1f} per orbit)")
            print()

        # ====================================================================
        # Method 1: Pure Newtonian
        # ====================================================================
        if verbose:
            print("  [1/3] Pure Newtonian integration...")

        import time
        t_start = time.time()

        bodies_newton = [copy.deepcopy(b) for b in bodies_template]
        traj_newton = integrate_newtonian_orbit(
            bodies_newton, K, dt, n_steps, save_every
        )

        t_elapsed = time.time() - t_start

        # Analyze precession
        analysis_newton = analyse_precession(traj_newton, params)
        prec_rad = analysis_newton['peri_per_orbit']
        prec_std = analysis_newton['peri_std']
        n_peri = analysis_newton['n_peri_cycles']

        # Compute drifts
        t_newton, x_rel, v_rel = traj_newton['t'], None, None
        x_rel = traj_newton['x'][:, 1, :] - traj_newton['x'][:, 0, :]
        v_rel = traj_newton['v'][:, 1, :] - traj_newton['v'][:, 0, :]

        elems_list = []
        for xr, vr in zip(x_rel, v_rel):
            elems = compute_orbit_elements(xr, vr, params.M_total, K)
            elems_list.append(elems)

        a_values = np.array([e['a'] for e in elems_list if np.isfinite(e['a'])])
        e_values = np.array([e['e'] for e in elems_list if np.isfinite(e['e'])])
        E_values = np.array([e['E'] for e in elems_list if np.isfinite(e['E'])])

        a_drift = (a_values[-1] - a_values[0]) / a_values[0] if len(a_values) > 0 else np.nan
        e_drift = (e_values[-1] - e_values[0]) if len(e_values) > 0 else np.nan
        E_drift = (E_values[-1] - E_values[0]) / abs(E_values[0]) if len(E_values) > 0 else np.nan

        results.append(ConvergenceResult(
            dt=dt,
            steps_per_orbit=steps_per_orbit,
            method='Newtonian',
            precession_arcsec=prec_rad * 206265.0,
            precession_rad=prec_rad,
            precession_std_arcsec=prec_std * 206265.0,
            n_periapsis=n_peri,
            energy_drift=E_drift,
            a_drift=a_drift,
            e_drift=e_drift,
            integration_time=t_elapsed,
        ))

        if verbose:
            print(f"        Precession: {prec_rad * 206265:.6f} ± {prec_std * 206265:.6f} arcsec/orbit ({n_peri} cycles)")
            print(f"        Energy drift: {E_drift:.3e}")
            print(f"        Time: {t_elapsed:.2f} s")

        # ====================================================================
        # Method 2: Incompressible Superfluid
        # ====================================================================
        if verbose:
            print("  [2/3] Incompressible superfluid integration...")

        t_start = time.time()

        bodies_slab = [copy.deepcopy(b) for b in bodies_template]
        opts = {
            'use_compressible': False,
            'save_every': save_every,
            'verbose': False,
        }
        traj_slab, _ = integrate_orbit(bodies_slab, medium, dt, n_steps, opts)

        t_elapsed = time.time() - t_start

        # Analyze precession
        analysis_slab = analyse_precession(traj_slab, params)
        prec_rad = analysis_slab['peri_per_orbit']
        prec_std = analysis_slab['peri_std']
        n_peri = analysis_slab['n_peri_cycles']

        # Compute drifts
        x_rel = traj_slab['x'][:, 1, :] - traj_slab['x'][:, 0, :]
        v_rel = traj_slab['v'][:, 1, :] - traj_slab['v'][:, 0, :]

        elems_list = []
        for xr, vr in zip(x_rel, v_rel):
            elems = compute_orbit_elements(xr, vr, params.M_total, K)
            elems_list.append(elems)

        a_values = np.array([e['a'] for e in elems_list if np.isfinite(e['a'])])
        e_values = np.array([e['e'] for e in elems_list if np.isfinite(e['e'])])
        E_values = np.array([e['E'] for e in elems_list if np.isfinite(e['E'])])

        a_drift = (a_values[-1] - a_values[0]) / a_values[0] if len(a_values) > 0 else np.nan
        e_drift = (e_values[-1] - e_values[0]) if len(e_values) > 0 else np.nan
        E_drift = (E_values[-1] - E_values[0]) / abs(E_values[0]) if len(E_values) > 0 else np.nan

        results.append(ConvergenceResult(
            dt=dt,
            steps_per_orbit=steps_per_orbit,
            method='Superfluid',
            precession_arcsec=prec_rad * 206265.0,
            precession_rad=prec_rad,
            precession_std_arcsec=prec_std * 206265.0,
            n_periapsis=n_peri,
            energy_drift=E_drift,
            a_drift=a_drift,
            e_drift=e_drift,
            integration_time=t_elapsed,
        ))

        if verbose:
            print(f"        Precession: {prec_rad * 206265:.6f} ± {prec_std * 206265:.6f} arcsec/orbit ({n_peri} cycles)")
            print(f"        Energy drift: {E_drift:.3e}")
            print(f"        Time: {t_elapsed:.2f} s")

        # ====================================================================
        # Method 3: GR 1PN (optional, slower)
        # ====================================================================
        if include_gr1pn:
            if verbose:
                print("  [3/3] GR 1PN integration...")

            t_start = time.time()

            bodies_gr = [copy.deepcopy(b) for b in bodies_template]
            traj_gr = integrate_gr1pn_orbit(
                bodies_gr, cs, K, dt, n_steps, save_every
            )

            t_elapsed = time.time() - t_start

            # Analyze precession
            analysis_gr = analyse_precession(traj_gr, params)
            prec_rad = analysis_gr['peri_per_orbit']
            prec_std = analysis_gr['peri_std']
            n_peri = analysis_gr['n_peri_cycles']

            # Compute drifts (same as above)
            x_rel = traj_gr['x'][:, 1, :] - traj_gr['x'][:, 0, :]
            v_rel = traj_gr['v'][:, 1, :] - traj_gr['v'][:, 0, :]

            elems_list = []
            for xr, vr in zip(x_rel, v_rel):
                elems = compute_orbit_elements(xr, vr, params.M_total, K)
                elems_list.append(elems)

            a_values = np.array([e['a'] for e in elems_list if np.isfinite(e['a'])])
            e_values = np.array([e['e'] for e in elems_list if np.isfinite(e['e'])])
            E_values = np.array([e['E'] for e in elems_list if np.isfinite(e['e'])])

            a_drift = (a_values[-1] - a_values[0]) / a_values[0] if len(a_values) > 0 else np.nan
            e_drift = (e_values[-1] - e_values[0]) if len(e_values) > 0 else np.nan
            E_drift = (E_values[-1] - E_values[0]) / abs(E_values[0]) if len(E_values) > 0 else np.nan

            results.append(ConvergenceResult(
                dt=dt,
                steps_per_orbit=steps_per_orbit,
                method='GR1PN',
                precession_arcsec=prec_rad * 206265.0,
                precession_rad=prec_rad,
                precession_std_arcsec=prec_std * 206265.0,
                n_periapsis=n_peri,
                energy_drift=E_drift,
                a_drift=a_drift,
                e_drift=e_drift,
                integration_time=t_elapsed,
            ))

            if verbose:
                print(f"        Precession: {prec_rad * 206265:.6f} ± {prec_std * 206265:.6f} arcsec/orbit ({n_peri} cycles)")
                print(f"        Energy drift: {E_drift:.3e}")
                print(f"        Time: {t_elapsed:.2f} s")

    return results, params


def print_results_table(results: List[ConvergenceResult], params: TwoBodyParams):
    """Print formatted results table."""

    # Compute GR prediction
    cs = 63239.7263
    K = params.K
    M_total = params.M_total
    a = params.a
    e = params.e
    GR_precession_arcsec = (6 * np.pi * K * M_total) / (a * cs**2 * (1 - e**2)) * 206265.0

    print()
    print("=" * 120)
    print("CONVERGENCE STUDY RESULTS")
    print("=" * 120)
    print()
    print(f"GR Prediction: {GR_precession_arcsec:.6f} arcsec/orbit")
    print(f"Target (10% of GR): {0.1 * GR_precession_arcsec:.6f} arcsec/orbit")
    print(f"Stretch goal (1% of GR): {0.01 * GR_precession_arcsec:.6f} arcsec/orbit")
    print()

    # Group by timestep
    dt_values = sorted(set(r.dt for r in results))

    print(f"{'dt [yr]':<12} {'Steps/orbit':<12} {'Method':<15} {'Prec [arcsec/orb]':<20} "
          f"{'ΔE/E':<12} {'Δa/a':<12} {'Time [s]':<10}")
    print("-" * 120)

    for dt in dt_values:
        dt_results = [r for r in results if r.dt == dt]
        for i, res in enumerate(dt_results):
            if i == 0:
                dt_str = f"{res.dt:.6f}"
                steps_str = f"{res.steps_per_orbit:.1f}"
            else:
                dt_str = ""
                steps_str = ""

            prec_str = f"{res.precession_arcsec:.6f} ± {res.precession_std_arcsec:.6f}"

            print(f"{dt_str:<12} {steps_str:<12} {res.method:<15} {prec_str:<20} "
                  f"{res.energy_drift:< 12.2e} {res.a_drift:< 12.2e} {res.integration_time:<10.2f}")

        # Add blank line between timesteps
        if dt != dt_values[-1]:
            print()

    print("-" * 120)
    print()

    # Analysis section
    print("=" * 120)
    print("ANALYSIS")
    print("=" * 120)
    print()

    # Find timesteps that meet criteria
    newton_results = [r for r in results if r.method == 'Newtonian']
    slab_results = [r for r in results if r.method == 'Superfluid']

    print("Precision Requirements:")
    print()

    for target_frac, target_name in [(0.1, "10% of GR"), (0.01, "1% of GR")]:
        target_prec = target_frac * GR_precession_arcsec
        print(f"  Target: <{target_prec:.6f} arcsec/orbit ({target_name})")

        # Check Newtonian
        newton_meets = [r for r in newton_results if abs(r.precession_arcsec) < target_prec]
        if newton_meets:
            best = min(newton_meets, key=lambda r: r.steps_per_orbit)
            print(f"    Newtonian:   dt ≤ {best.dt:.6f} yr ({best.steps_per_orbit:.0f} steps/orbit)")
        else:
            print(f"    Newtonian:   NONE of the tested timesteps meet this requirement")

        # Check Superfluid
        slab_meets = [r for r in slab_results if abs(r.precession_arcsec) < target_prec]
        if slab_meets:
            best = min(slab_meets, key=lambda r: r.steps_per_orbit)
            print(f"    Superfluid:  dt ≤ {best.dt:.6f} yr ({best.steps_per_orbit:.0f} steps/orbit)")
        else:
            print(f"    Superfluid:  NONE of the tested timesteps meet this requirement")
        print()

    # Method agreement
    print("Method Agreement:")
    print()
    for dt in dt_values:
        newton_res = [r for r in newton_results if r.dt == dt][0]
        slab_res = [r for r in slab_results if r.dt == dt][0]

        diff = abs(newton_res.precession_arcsec - slab_res.precession_arcsec)
        rel_diff = diff / max(abs(newton_res.precession_arcsec), 1e-10)

        print(f"  dt = {dt:.6f} yr: |Newtonian - Superfluid| = {diff:.6f} arcsec/orbit "
              f"({rel_diff*100:.2f}% relative)")

    print()
    print("=" * 120)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Validate Newtonian gravity precision for 1PN studies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--dt', nargs='+', type=float,
        default=[0.002, 0.001, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5],
        help='Timestep values to test [yr] (default: 8 values from 0.002 to 1e-5)'
    )

    parser.add_argument(
        '--n-orbits', type=float, default=10.0,
        help='Number of orbits to integrate (default: 10)'
    )

    parser.add_argument(
        '--eccentricity', type=float, default=0.1,
        help='Orbital eccentricity (default: 0.1, Mercury-like)'
    )

    parser.add_argument(
        '--include-gr', action='store_true',
        help='Include GR 1PN comparison (slower)'
    )

    parser.add_argument(
        '--quick', action='store_true',
        help='Quick test with only 3 timesteps: [0.001, 0.0001, 1e-5]'
    )

    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress output'
    )

    args = parser.parse_args()

    # Handle quick mode
    if args.quick:
        dt_values = [0.001, 0.0001, 1e-5]
    else:
        dt_values = sorted(args.dt, reverse=True)  # Coarse to fine

    # Run convergence study
    results, params = run_convergence_study(
        dt_values=dt_values,
        n_orbits=args.n_orbits,
        eccentricity=args.eccentricity,
        include_gr1pn=args.include_gr,
        verbose=not args.quiet,
    )

    # Print results table
    print_results_table(results, params)

    return 0


if __name__ == '__main__':
    sys.exit(main())
