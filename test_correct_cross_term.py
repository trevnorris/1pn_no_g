#!/usr/bin/env python3
"""Test the correctly-signed cross-term formula."""

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
print("TESTING CORRECTED CROSS-TERM WITH MINUS SIGN")
print("="*70)

a_idx = 0
body_a = sun
x_a = body_a.x
R_a = body_a.R
Q_a = body_a.Q
rho0 = medium.rho0

n_points = 256
normals = fibonacci_sphere(n_points)
dA = (4.0 * np.pi * R_a**2) / n_points

# Compute cross-term with CORRECT sign (negative)
F_cross_correct = np.zeros(3, dtype=np.longdouble)
F_cross_wrong = np.zeros(3, dtype=np.longdouble)

rho0_ld = np.longdouble(rho0)
dA_ld = np.longdouble(dA)

for i in range(n_points):
    n_i = normals[i]
    x_i = x_a + R_a * n_i

    v_self_i = v_self(x_i, x_a, Q_a, rho0)
    v_ext_i = v_ext_at(x_i, bodies, a_idx, rho0)

    v_self_dot_n = np.dot(v_self_i, n_i)
    v_ext_dot_n = np.dot(v_ext_i, n_i)

    # Wrong formula (from decompose_momentum_integrand)
    integrand_wrong = rho0 * (v_self_i * v_ext_dot_n + v_ext_i * v_self_dot_n)

    # Correct formula (with minus sign for outward normal)
    integrand_correct = -rho0 * (v_self_i * v_ext_dot_n + v_ext_i * v_self_dot_n)

    F_cross_wrong += integrand_wrong.astype(np.longdouble) * dA_ld
    F_cross_correct += integrand_correct.astype(np.longdouble) * dA_ld

F_cross_wrong = F_cross_wrong.astype(np.float64)
F_cross_correct = F_cross_correct.astype(np.float64)

# Get analytic force
from slab.surface import force_incompressible_analytic
F_analytic = force_incompressible_analytic(a_idx, bodies, medium)

print(f"\nAnalytic force:              {F_analytic}")
print(f"Cross-term (wrong sign):     {F_cross_wrong}")
print(f"Cross-term (correct sign):   {F_cross_correct}")

error_wrong = np.linalg.norm(F_cross_wrong - F_analytic) / np.linalg.norm(F_analytic)
error_correct = np.linalg.norm(F_cross_correct - F_analytic) / np.linalg.norm(F_analytic)

print(f"\nRelative error (wrong):      {error_wrong:.6e}")
print(f"Relative error (correct):    {error_correct:.6e}")

print(f"\n{'✓ PASSES AUDIT!' if error_correct < 1e-3 else '✗ Still fails'}")

if error_correct < 1e-3:
    print("\n" + "="*70)
    print("SOLUTION FOUND!")
    print("="*70)
    print("\nThe fix is to:")
    print("1. Use ONLY the cross-term (v_self × v_ext + v_ext × v_self)")
    print("2. Add a MINUS SIGN: F = -∫ ρ (cross-term) dA")
    print("\nThis avoids the catastrophic self×self cancellation")
    print("and uses the correct sign convention for momentum flux.")
