#!/usr/bin/env python3
"""Compute all terms of the momentum flux integral."""

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
print("COMPUTING ALL TERMS OF MOMENTUM FLUX INTEGRAL")
print("="*70)

a_idx = 0
body_a = sun
x_a = body_a.x
R_a = body_a.R
Q_a = body_a.Q
rho0 = medium.rho0

n_points = 512
normals = fibonacci_sphere(n_points)
dA = (4.0 * np.pi * R_a**2) / n_points

# Compute all four terms separately WITH MINUS SIGN
F_self_self = np.zeros(3, dtype=np.longdouble)
F_self_ext = np.zeros(3, dtype=np.longdouble)
F_ext_self = np.zeros(3, dtype=np.longdouble)
F_ext_ext = np.zeros(3, dtype=np.longdouble)

rho0_ld = np.longdouble(rho0)
dA_ld = np.longdouble(dA)

for i in range(n_points):
    n_i = normals[i]
    x_i = x_a + R_a * n_i

    v_self_i = v_self(x_i, x_a, Q_a, rho0)
    v_ext_i = v_ext_at(x_i, bodies, a_idx, rho0)

    v_self_dot_n = np.dot(v_self_i, n_i)
    v_ext_dot_n = np.dot(v_ext_i, n_i)

    # WITH MINUS SIGN (correct for outward normal)
    F_self_self += (-rho0 * v_self_i * v_self_dot_n * dA).astype(np.longdouble)
    F_self_ext += (-rho0 * v_self_i * v_ext_dot_n * dA).astype(np.longdouble)
    F_ext_self += (-rho0 * v_ext_i * v_self_dot_n * dA).astype(np.longdouble)
    F_ext_ext += (-rho0 * v_ext_i * v_ext_dot_n * dA).astype(np.longdouble)

F_self_self = F_self_self.astype(np.float64)
F_self_ext = F_self_ext.astype(np.float64)
F_ext_self = F_ext_self.astype(np.float64)
F_ext_ext = F_ext_ext.astype(np.float64)

F_total_decomposed = F_self_self + F_self_ext + F_ext_self + F_ext_ext
F_cross_only = F_self_ext + F_ext_self

# Get analytic force
from slab.surface import force_incompressible_analytic
F_analytic = force_incompressible_analytic(a_idx, bodies, medium)

print("\nIndividual terms (with correct sign):")
print(f"  F_self×self = {F_self_self}")
print(f"  F_self×ext  = {F_self_ext}")
print(f"  F_ext×self  = {F_ext_self}")
print(f"  F_ext×ext   = {F_ext_ext}")

print(f"\nSums:")
print(f"  F_cross (self×ext + ext×self) = {F_cross_only}")
print(f"  F_total (all four terms)      = {F_total_decomposed}")

print(f"\nExpected:")
print(f"  F_analytic = {F_analytic}")

print(f"\nMagnitude comparison:")
print(f"  |F_self×self| = {np.linalg.norm(F_self_self):.6e}")
print(f"  |F_cross|     = {np.linalg.norm(F_cross_only):.6e}")
print(f"  |F_ext×ext|   = {np.linalg.norm(F_ext_ext):.6e}")
print(f"  |F_total|     = {np.linalg.norm(F_total_decomposed):.6e}")
print(f"  |F_analytic|  = {np.linalg.norm(F_analytic):.6e}")

error_cross = np.linalg.norm(F_cross_only - F_analytic) / np.linalg.norm(F_analytic)
error_total = np.linalg.norm(F_total_decomposed - F_analytic) / np.linalg.norm(F_analytic)

print(f"\nRelative errors:")
print(f"  Cross-term only:  {error_cross:.6e}")
print(f"  All terms (total): {error_total:.6e}")

if np.linalg.norm(F_self_self) > 1e-20:
    print(f"\n⚠ WARNING: self×self term is NOT negligible!")
    print(f"  This explains the huge error in the original quadrature.")
