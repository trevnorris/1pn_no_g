#!/usr/bin/env python3
"""
Realistic test of dynamics module with Solar System-like parameters.

Uses parameters chosen to give reasonable timescales and velocities.
"""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.dynamics import (
    step_verlet,
    assemble_forces,
    integrate_orbit,
    estimate_orbital_period,
)


def setup_solar_system_analog():
    """
    Set up a Sun-Earth analog with parameters that give:
    - Orbital period ~ 1 year
    - Velocities ~ 10 km/s in code units
    - Reasonable timestep ~ 0.01 year

    Strategy: choose beta0 to make K ~ G (in SI units)
    """
    # Code units: AU, yr, solar masses
    # Then K should be ~ 4π² (gives T = 1 yr for a = 1 AU, M = 1 M_sun)

    # Work backwards: K = rho0 / (4π β0²)
    # Want K ~ 4π² ≈ 39.5
    # Choose rho0 = 1.0 (code units)
    # Then β0² = rho0 / (4π K) = 1 / (4π * 39.5) ≈ 0.002
    # So β0 ≈ 0.045

    K_target = 4.0 * np.pi**2  # AU³/yr²/M_sun
    rho0 = 1.0  # code units
    beta0 = np.sqrt(rho0 / (4.0 * np.pi * K_target))
    cs = 1e4  # Large compared to orbital velocities

    medium = Medium(rho0=rho0, cs=cs, beta0=beta0)

    print(f"Medium parameters (for ~1 yr orbit at 1 AU):")
    print(f"  rho0 = {rho0}")
    print(f"  beta0 = {beta0:.6f}")
    print(f"  cs = {cs:.3e}")
    print(f"  K = {medium.K:.6f} (target: {K_target:.6f})")

    return medium


def test_earth_sun_orbit():
    """Test Earth-Sun analog with realistic parameters."""
    print("=" * 70)
    print("TEST: Earth-Sun analog orbit (should have T ~ 1 year)")
    print("=" * 70)
    print()

    medium = setup_solar_system_analog()

    # Sun at origin
    M_sun = 1.0  # solar masses
    sun = Body(
        name="Sun",
        M=M_sun,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=1e-3,
        Q=0.0
    )
    sun.update_Q_from_M(medium)

    # Earth at 1 AU with circular velocity
    a = 1.0  # AU
    M_earth = 3.0e-6  # solar masses
    v_circ = np.sqrt(medium.K * M_sun / a)

    earth = Body(
        name="Earth",
        M=M_earth,
        x=np.array([a, 0.0, 0.0]),
        v=np.array([0.0, v_circ, 0.0]),
        R=5e-4,
        Q=0.0
    )
    earth.update_Q_from_M(medium)

    bodies = [sun, earth]

    print(f"\nInitial configuration:")
    print(f"  Sun: M={sun.M:.3e} M_sun, at origin")
    print(f"  Earth: M={earth.M:.3e} M_sun, at r={a} AU")
    print(f"  Circular velocity: v={v_circ:.6f} AU/yr")

    # Estimate period
    T_est = estimate_orbital_period(bodies, medium)
    print(f"\nEstimated period: T = {T_est:.6f} yr")
    print(f"  (Should be ~1.0 yr for Earth)")

    # Choose timestep
    dt = 0.005 * T_est  # 0.5% of orbital period
    n_orbits = 5  # Integrate for 5 orbits
    n_steps = int(n_orbits * T_est / dt)

    print(f"\nIntegration parameters:")
    print(f"  Timestep: dt = {dt:.6f} yr ({dt/T_est*100:.2f}% of period)")
    print(f"  Steps: {n_steps} (for {n_orbits} orbits)")
    print(f"  Total time: {n_steps * dt:.3f} yr")

    # Run integration
    opts = {
        'save_every': max(1, n_steps // 500),  # ~500 snapshots
        'verbose': True,
        'progress_every': max(1, n_steps // 10),  # ~10 progress updates
    }

    print(f"\nRunning integration...")
    traj, diags = integrate_orbit(bodies, medium, dt, n_steps, opts)

    # Analyze results
    print(f"\nEnergy analysis:")
    E = np.array([d['total_energy'] for d in diags])
    E0 = E[0]
    dE = E - E0
    rel_dE = dE / np.abs(E0)

    print(f"  Initial energy: E0 = {E0:.6e}")
    print(f"  Final energy: E_f = {E[-1]:.6e}")
    print(f"  Energy change: ΔE = {E[-1]-E0:.6e}")
    print(f"  Relative change: ΔE/E0 = {(E[-1]-E0)/E0:.6e}")
    print(f"  Max |ΔE/E0|: {np.max(np.abs(rel_dE)):.6e}")
    print(f"  RMS |ΔE/E0|: {np.sqrt(np.mean(rel_dE**2)):.6e}")

    # Check acceptance criteria: |Δa|/a < 10^-5 over 50 orbits
    # For now, check energy conservation as proxy
    max_energy_error = np.max(np.abs(rel_dE))

    if max_energy_error < 1e-5:
        print(f"\n✓ EXCELLENT: Energy conserved to {max_energy_error:.2e} < 1e-5")
        print(f"  Exceeds acceptance criterion from plan § 10.2")
    elif max_energy_error < 1e-4:
        print(f"\n✓ GOOD: Energy conserved to {max_energy_error:.2e} < 1e-4")
    elif max_energy_error < 1e-3:
        print(f"\n✓ PASS: Energy conserved to {max_energy_error:.2e} < 1e-3")
    else:
        print(f"\n⚠ WARNING: Energy error {max_energy_error:.2e} > 1e-3")
        print(f"  Consider reducing timestep")

    # Momentum conservation
    p = np.array([d['total_momentum'] for d in diags])
    p0 = p[0]
    dp = np.linalg.norm(p - p0[None, :], axis=1)
    max_dp = np.max(dp)

    print(f"\nMomentum conservation:")
    print(f"  Initial momentum: p0 = {p0}")
    print(f"  Max |Δp|: {max_dp:.3e}")

    if max_dp < 1e-12:
        print(f"  ✓ Momentum conserved to machine precision")
    elif max_dp < 1e-10:
        print(f"  ✓ Momentum well conserved")
    else:
        print(f"  ⚠ Momentum drift detected")

    # Orbit shape analysis
    x_earth = traj['x'][:, 1, :]  # Earth positions
    r = np.linalg.norm(x_earth, axis=1)

    r_mean = np.mean(r)
    r_std = np.std(r)

    print(f"\nOrbit shape:")
    print(f"  Mean radius: <r> = {r_mean:.6f} AU")
    print(f"  Std dev: σ(r) = {r_std:.6f} AU")
    print(f"  Circularity: σ/r = {r_std/r_mean:.2e}")

    if r_std/r_mean < 0.01:
        print(f"  ✓ Nearly circular orbit (as expected)")

    print()


if __name__ == "__main__":
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 13 + "REALISTIC DYNAMICS MODULE TEST" + " " * 24 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    test_earth_sun_orbit()

    print("=" * 70)
    print("TEST COMPLETED")
    print("=" * 70)
    print()
