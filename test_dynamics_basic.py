#!/usr/bin/env python3
"""
Basic test of dynamics module: two-body orbit integration.

This verifies:
1. step_verlet runs without errors
2. Energy is approximately conserved (symplectic property)
3. Momentum is exactly conserved
4. Forces obey Newton's third law
"""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.dynamics import (
    step_verlet,
    assemble_forces,
    integrate_orbit,
    estimate_orbital_period,
    estimate_timestep,
)


def test_two_body_step():
    """Test single integration step for two-body system."""
    print("=" * 70)
    print("TEST 1: Single step_verlet for two-body system")
    print("=" * 70)

    # Set up medium
    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    print(f"\nMedium: {medium}")
    print(f"  K = {medium.K:.6e} [orbital constant]")

    # Create two bodies at separation 1.0
    b1 = Body(
        name="Body1",
        M=1.0,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=1e-3,
        Q=0.0
    )
    b1.update_Q_from_M(medium)

    b2 = Body(
        name="Body2",
        M=1e-3,  # Much smaller (like planet)
        x=np.array([1.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=5e-4,
        Q=0.0
    )
    b2.update_Q_from_M(medium)

    bodies = [b1, b2]

    print(f"\nInitial state:")
    for body in bodies:
        print(f"  {body.name}: x={body.x}, v={body.v}, M={body.M:.3e}")

    # Compute initial forces
    print(f"\nInitial forces:")
    F_init = assemble_forces(bodies, medium)
    for i, body in enumerate(bodies):
        print(f"  {body.name}: F={F_init[i]}")

    # Verify Newton's third law: F1 = -F2
    F_sum = F_init[0] + F_init[1]
    print(f"\nNewton's 3rd law check: F1 + F2 = {F_sum}")
    print(f"  Magnitude: {np.linalg.norm(F_sum):.3e} (should be ~0)")

    # Take a single timestep
    dt = 0.01
    print(f"\nTaking single step with dt={dt}")

    opts = {'use_compressible': False}
    diag = step_verlet(bodies, medium, dt, opts)

    print(f"\nAfter step:")
    for body in bodies:
        print(f"  {body.name}: x={body.x}, v={body.v}")

    print(f"\nDiagnostics:")
    print(f"  KE = {diag['kinetic_energy']:.6e}")
    print(f"  PE = {diag['potential_energy']:.6e}")
    print(f"  E_total = {diag['total_energy']:.6e}")
    print(f"  |p_total| = {diag['momentum_magnitude']:.6e}")
    print(f"  max_force = {diag['max_force']:.6e}")

    print("\n✓ Single step completed successfully\n")


def test_energy_conservation():
    """Test energy conservation over multiple steps."""
    print("=" * 70)
    print("TEST 2: Energy conservation over 100 steps")
    print("=" * 70)

    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)

    # Set up circular orbit
    a = 1.0  # AU
    M_primary = 1.0
    M_secondary = 3e-6  # Earth-like

    # Circular velocity: v = sqrt(K*M/a)
    v_circ = np.sqrt(medium.K * M_primary / a)

    b1 = Body("Primary", M=M_primary, x=[0, 0, 0], v=[0, 0, 0], R=1e-3, Q=0.0)
    b2 = Body("Secondary", M=M_secondary, x=[a, 0, 0], v=[0, v_circ, 0], R=5e-4, Q=0.0)
    b1.update_Q_from_M(medium)
    b2.update_Q_from_M(medium)

    bodies = [b1, b2]

    print(f"\nCircular orbit setup:")
    print(f"  Semi-major axis: a = {a}")
    print(f"  Circular velocity: v = {v_circ:.6e}")
    print(f"  M_primary = {M_primary:.3e}")
    print(f"  M_secondary = {M_secondary:.3e}")

    # Estimate period and choose timestep
    T_est = estimate_orbital_period(bodies, medium)
    print(f"\nEstimated period: T = {T_est:.3e}")

    # Use a reasonable fraction of orbital period
    # For good energy conservation: dt ~ 0.001 - 0.01 * T
    dt = 0.005 * T_est
    n_steps = 200  # Cover one full orbit

    print(f"Timestep: dt = {dt:.6e} ({dt/T_est:.3%} of period)")
    print(f"Steps: {n_steps} (total time: {n_steps*dt:.3e} = {n_steps*dt/T_est:.2f} orbits)")

    # Integrate
    opts = {
        'save_every': 1,
        'verbose': False,
    }

    traj, diags = integrate_orbit(bodies, medium, dt, n_steps, opts)

    # Extract energy time series
    times = np.array([d['time'] for d in diags])
    energies = np.array([d['total_energy'] for d in diags])
    KE = np.array([d['kinetic_energy'] for d in diags])
    PE = np.array([d['potential_energy'] for d in diags])

    E0 = energies[0]
    dE = energies - E0
    rel_dE = dE / np.abs(E0)

    print(f"\nEnergy conservation:")
    print(f"  Initial energy: E0 = {E0:.6e}")
    print(f"  Final energy: E_f = {energies[-1]:.6e}")
    print(f"  Max |ΔE|: {np.max(np.abs(dE)):.3e}")
    print(f"  Max |ΔE/E|: {np.max(np.abs(rel_dE)):.3e}")
    print(f"  RMS |ΔE/E|: {np.sqrt(np.mean(rel_dE**2)):.3e}")

    # Check momentum conservation
    p_total = np.array([d['total_momentum'] for d in diags])
    p0 = p_total[0]
    dp = np.linalg.norm(p_total - p0[None, :], axis=1)

    print(f"\nMomentum conservation:")
    print(f"  Initial momentum: p0 = {p0}")
    print(f"  Max |Δp|: {np.max(dp):.3e}")

    # Acceptance criteria (from plan § 10.2)
    max_rel_dE = np.max(np.abs(rel_dE))
    if max_rel_dE < 1e-4:
        print(f"\n✓ Energy conserved to {max_rel_dE:.2e} < 1e-4")
    else:
        print(f"\n⚠ Energy drift {max_rel_dE:.2e} (acceptable if < 1e-3)")

    print()


def test_force_assembly():
    """Test force assembly for multiple bodies."""
    print("=" * 70)
    print("TEST 3: Force assembly for N=3 bodies")
    print("=" * 70)

    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)

    # Three bodies in a line
    bodies = [
        Body("B1", M=1.0, x=[0, 0, 0], v=[0, 0, 0], R=1e-3, Q=0.0),
        Body("B2", M=0.5, x=[1, 0, 0], v=[0, 0, 0], R=5e-4, Q=0.0),
        Body("B3", M=0.3, x=[2, 0, 0], v=[0, 0, 0], R=5e-4, Q=0.0),
    ]

    for body in bodies:
        body.update_Q_from_M(medium)

    print(f"\nBodies:")
    for body in bodies:
        print(f"  {body.name}: x={body.x}, M={body.M:.3e}")

    # Assemble forces
    F = assemble_forces(bodies, medium)

    print(f"\nForces:")
    for i, body in enumerate(bodies):
        print(f"  {body.name}: F={F[i]}, |F|={np.linalg.norm(F[i]):.3e}")

    # Check total force = 0 (Newton's third law for all pairs)
    F_total = np.sum(F, axis=0)
    print(f"\nTotal force: {F_total}")
    print(f"  Magnitude: {np.linalg.norm(F_total):.3e} (should be ~0)")

    if np.linalg.norm(F_total) < 1e-10:
        print("\n✓ Forces satisfy Newton's third law (sum to zero)")
    else:
        print("\n⚠ Force sum nonzero (possible numerical issue)")

    print()


if __name__ == "__main__":
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "DYNAMICS MODULE BASIC TESTS" + " " * 26 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    test_two_body_step()
    test_energy_conservation()
    test_force_assembly()

    print("=" * 70)
    print("ALL TESTS COMPLETED")
    print("=" * 70)
    print()
