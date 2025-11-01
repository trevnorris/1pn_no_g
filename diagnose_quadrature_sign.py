#!/usr/bin/env python3
"""Diagnose why quadrature has wrong sign."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.field import v_self, v_ext_at
from slab.geometry import fibonacci_sphere

# Create simple two-body system (Sun-Mercury)
sun = Body(
    name="Sun",
    x=np.array([0.0, 0.0, 0.0]),
    v=np.array([0.0, 0.0, 0.0]),
    M=1.0,
    Q=1.0e-10,
    R=0.01,
)

mercury = Body(
    name="Mercury",
    x=np.array([0.387, 0.0, 0.0]),
    v=np.array([0.0, 4.53e-11, 0.0]),
    M=3.3e-7,
    Q=3.3e-17,
    R=0.001,
)

bodies = [sun, mercury]
medium = Medium(rho0=1.0, cs=10000.0, beta0=1.0e10)

print("="*70)
print("DIAGNOSING QUADRATURE SIGN ERROR")
print("="*70)

a_idx = 0
body_a = sun
x_a = body_a.x
R_a = body_a.R
Q_a = body_a.Q
rho0 = medium.rho0

# Test at a single point: normal pointing in +x direction
normal = np.array([1.0, 0.0, 0.0])
x_surface = x_a + R_a * normal

print(f"\nTest point:")
print(f"  Normal (outward): {normal}")
print(f"  Surface position: {x_surface}")

# Compute velocities
v_self_vec = v_self(x_surface, x_a, Q_a, rho0)
v_ext_vec = v_ext_at(x_surface, bodies, a_idx, rho0)

print(f"\nVelocity components:")
print(f"  v_self = {v_self_vec}")
print(f"  v_ext = {v_ext_vec}")

# Compute normal components
v_self_dot_n = np.dot(v_self_vec, normal)
v_ext_dot_n = np.dot(v_ext_vec, normal)

print(f"\nNormal components:")
print(f"  v_self·n = {v_self_dot_n:.6e}")
print(f"  v_ext·n = {v_ext_dot_n:.6e}")

# The cross-term integrand is:
# ρ₀ [v_self(v_ext·n) + v_ext(v_self·n)]

cross_term = rho0 * (v_self_vec * v_ext_dot_n + v_ext_vec * v_self_dot_n)

print(f"\nCross-term integrand:")
print(f"  ρ₀ [v_self(v_ext·n) + v_ext(v_self·n)] = {cross_term}")

print("\n" + "="*70)
print("UNDERSTANDING THE PHYSICS")
print("="*70)

print(f"\nFor a sink (Q > 0), v_self points INWARD:")
print(f"  At r = R*x̂, v_self should point in -x̂ direction")
print(f"  v_self·n̂ (with n̂ = +x̂) should be NEGATIVE")
print(f"  Actual: v_self·n̂ = {v_self_dot_n:.6e} ✓" if v_self_dot_n < 0 else f"  Actual: v_self·n̂ = {v_self_dot_n:.6e} ✗")

print(f"\nFor attractive force, v_ext should point TOWARD Mercury:")
print(f"  Mercury is at +x, so v_ext should point in +x̂ direction")
print(f"  v_ext·n̂ (with n̂ = +x̂) should be POSITIVE")
print(f"  Actual: v_ext·n̂ = {v_ext_dot_n:.6e} ✓" if v_ext_dot_n > 0 else f"  Actual: v_ext·n̂ = {v_ext_dot_n:.6e} ✗")

print("\n" + "="*70)
print("ANALYZING THE CROSS-TERM")
print("="*70)

term1 = rho0 * v_self_vec * v_ext_dot_n
term2 = rho0 * v_ext_vec * v_self_dot_n

print(f"\nTerm 1: ρ₀ v_self (v_ext·n)")
print(f"  v_self = {v_self_vec}")
print(f"  v_ext·n = {v_ext_dot_n:.6e}")
print(f"  Term 1 = {term1}")

print(f"\nTerm 2: ρ₀ v_ext (v_self·n)")
print(f"  v_ext = {v_ext_vec}")
print(f"  v_self·n = {v_self_dot_n:.6e}")
print(f"  Term 2 = {term2}")

print(f"\nSum: {term1 + term2}")

print("\n" + "="*70)
print("EXPECTED BEHAVIOR")
print("="*70)

print("\nThe cross-term has TWO parts:")
print("  1. v_self (v_ext·n): negative vector × positive = NEGATIVE contribution")
print("  2. v_ext (v_self·n): positive vector × negative = NEGATIVE contribution")
print("\nBoth terms contribute NEGATIVELY!")
print("But we want a POSITIVE force (toward Mercury)!")

print("\n" + "="*70)
print("THE BUG")
print("="*70)

print("\nThe formula in decompose_momentum_integrand uses:")
print("  F = ∫ ρ₀ v(v·n) dA")
print("\nBut the normal n points OUTWARD.")
print("For momentum flux, we want the integral over INWARD flow!")
print("\nThe correct formula should have a MINUS SIGN:")
print("  F = -∫ ρ₀ v(v·n) dA  (with outward n)")
print("OR use inward normal:")
print("  F = ∫ ρ₀ v(v·(-n)) dA")

# Let's check what sign we need
F_expected = rho0 * Q_a * v_ext_at(x_a, bodies, a_idx, rho0)
print(f"\nExpected force: {F_expected}")
print(f"Cross-term (current): {cross_term} (wrong sign!)")
print(f"Cross-term (negated): {-cross_term} (correct sign!)")
