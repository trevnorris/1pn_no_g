#!/usr/bin/env python3
"""Quick test/demo of the Medium and Body dataclasses."""

import numpy as np
from slab.medium import Medium
from slab.bodies import Body


def test_medium():
    """Test Medium dataclass."""
    print("=" * 70)
    print("TESTING MEDIUM DATACLASS")
    print("=" * 70)

    # Create a medium with typical parameters
    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=0.0)
    print(medium)
    print()

    # Test K property
    K_expected = 1.0 / (4.0 * np.pi * 1e20)
    print(f"K computed:  {medium.K:.6e}")
    print(f"K expected:  {K_expected:.6e}")
    print(f"Match: {np.isclose(medium.K, K_expected)}")
    print()

    # Test constant beta
    print("Testing beta(rho) with gamma_beta=0:")
    for rho_frac in [0.5, 1.0, 2.0]:
        rho = rho_frac * medium.rho0
        beta_val = medium.beta(rho)
        print(f"  rho={rho:.1f} -> beta={beta_val:.3e} (should be constant)")
    print()

    # Test rarefaction
    print("Testing beta(rho) with gamma_beta=1:")
    medium_rare = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=1.0)
    for rho_frac in [0.5, 1.0, 2.0]:
        rho = rho_frac * medium_rare.rho0
        beta_val = medium_rare.beta(rho)
        beta_expected = 1e10 * rho_frac**(-1.0)
        print(f"  rho={rho:.1f} -> beta={beta_val:.3e} (expected {beta_expected:.3e})")
    print()


def test_bodies():
    """Test Body dataclass."""
    print("=" * 70)
    print("TESTING BODY DATACLASS")
    print("=" * 70)

    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    print(f"Medium: K = {medium.K:.3e}")
    print()

    # Create Sun at origin
    sun = Body(
        name="Sun",
        M=1.0,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=1e-3,
        Q=0.0
    )
    sun.update_Q_from_M(medium)
    print("Sun (after update_Q_from_M):")
    print(sun)
    print()

    # Create Mercury in circular orbit
    a = 0.387  # semi-major axis in AU (code units)
    v_circ = np.sqrt(medium.K * sun.M / a)
    print(f"Circular orbit at a={a} AU requires v={v_circ:.6f}")

    mercury = Body(
        name="Mercury",
        M=1.66e-7,
        x=np.array([a, 0.0, 0.0]),
        v=np.array([0.0, v_circ, 0.0]),
        R=5e-4,
        Q=0.0
    )
    mercury.update_Q_from_M(medium)
    print("Mercury (after update_Q_from_M):")
    print(mercury)
    print()

    # Test kinetic energy
    print("Kinetic energies:")
    print(f"  Sun: KE = {sun.kinetic_energy:.6e}")
    print(f"  Mercury: KE = {mercury.kinetic_energy:.6e}")
    print()

    # Test M-Q synchronization
    print("Testing M-Q synchronization:")
    old_M = mercury.M
    old_Q = mercury.Q
    mercury.Q *= 1.1  # Increase Q by 10%
    print(f"  Increased Q from {old_Q:.3e} to {mercury.Q:.3e}")
    mercury.update_M_from_Q(medium)
    print(f"  M updated from {old_M:.3e} to {mercury.M:.3e}")
    print(f"  Relative change: {(mercury.M/old_M - 1)*100:.1f}%")
    print()


def test_system_energy():
    """Test two-body system energy."""
    print("=" * 70)
    print("TESTING TWO-BODY SYSTEM ENERGY")
    print("=" * 70)

    medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)

    # Sun at origin
    sun = Body(
        name="Sun",
        M=1.0,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=1e-3,
        Q=1e-10
    )

    # Mercury in circular orbit
    a = 0.387
    v_circ = np.sqrt(medium.K * sun.M / a)
    mercury = Body(
        name="Mercury",
        M=1.66e-7,
        x=np.array([a, 0.0, 0.0]),
        v=np.array([0.0, v_circ, 0.0]),
        R=5e-4,
        Q=1.66e-17
    )

    # Compute energies
    T_total = sun.kinetic_energy + mercury.kinetic_energy

    r_vec = mercury.x - sun.x
    r = np.linalg.norm(r_vec)
    U_fluid = -medium.K * sun.M * mercury.M / r

    E_total = T_total + U_fluid

    print(f"Separation: r = {r:.6f}")
    print(f"Kinetic energy: T = {T_total:.6e}")
    print(f"Fluid potential energy: U = {U_fluid:.6e}")
    print(f"Total energy: E = {E_total:.6e}")
    print()

    # For circular orbit: T = -U/2 (virial theorem)
    print("Virial check (circular orbit should have T = -U/2):")
    print(f"  T = {T_total:.6e}")
    print(f"  -U/2 = {-U_fluid/2:.6e}")
    print(f"  Ratio T/(-U/2) = {T_total / (-U_fluid/2):.6f} (should be ~1.0)")
    print()


if __name__ == "__main__":
    test_medium()
    test_bodies()
    test_system_energy()

    print("=" * 70)
    print("ALL TESTS COMPLETED")
    print("=" * 70)
