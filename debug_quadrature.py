#!/usr/bin/env python3
"""Debug the quadrature calculation in detail."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.field import v_total, v_self, v_ext_at
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
print("DEBUGGING QUADRATURE IN DETAIL")
print("="*70)

# We're computing the force on body 0 (Sun)
a_idx = 0
body_a = bodies[a_idx]
x_a = body_a.x
R_a = body_a.R
rho0 = medium.rho0

print(f"\nComputing force on body {a_idx} ({body_a.name})")
print(f"Position: {x_a}")
print(f"Control surface radius: R = {R_a}")

# Generate a few sample points on the sphere
n_points = 8  # Just a few for debugging
normals = fibonacci_sphere(n_points)
sphere_area = 4.0 * np.pi * R_a**2
dA = sphere_area / n_points

print(f"\nQuadrature setup:")
print(f"  n_points = {n_points}")
print(f"  Sphere area = 4π R² = {sphere_area:.6e}")
print(f"  dA = {dA:.6e}")

print("\n" + "="*70)
print("SAMPLING POINTS ON CONTROL SURFACE")
print("="*70)

force_accumulator = np.zeros(3)

for i in range(min(5, n_points)):  # Show first 5 points
    n_i = normals[i]
    x_i = x_a + R_a * n_i

    print(f"\nPoint {i}:")
    print(f"  Normal: {n_i}")
    print(f"  Position on surface: {x_i}")

    # Compute velocity components
    v_self_i = v_self(x_i, x_a, body_a.Q, rho0)
    v_ext_i = v_ext_at(x_i, bodies, a_idx, rho0)
    v_total_i = v_total(x_i, bodies, rho0)

    print(f"  v_self:  {v_self_i}")
    print(f"  v_ext:   {v_ext_i}")
    print(f"  v_total: {v_total_i}")
    print(f"  v_self + v_ext: {v_self_i + v_ext_i}")
    print(f"  Agreement: {np.allclose(v_total_i, v_self_i + v_ext_i)}")

    # Compute normal component
    v_dot_n = np.dot(v_total_i, n_i)
    print(f"  v·n: {v_dot_n:.6e}")

    # Momentum flux contribution
    flux = rho0 * v_total_i * v_dot_n * dA
    print(f"  ρ₀ v(v·n) dA: {flux}")

    force_accumulator += flux

print(f"\n" + "="*70)
print(f"Force from {min(5, n_points)} sample points: {force_accumulator}")

# Now compute full quadrature with 256 points
print("\n" + "="*70)
print("FULL QUADRATURE WITH n_points=256")
print("="*70)

n_points = 256
normals = fibonacci_sphere(n_points)
dA = (4.0 * np.pi * R_a**2) / n_points

force_full = np.zeros(3, dtype=np.longdouble)
rho0_ld = np.longdouble(rho0)
dA_ld = np.longdouble(dA)

for i in range(n_points):
    n_i = normals[i]
    x_i = x_a + R_a * n_i

    v_i = v_total(x_i, bodies, rho0)
    v_dot_n = np.dot(v_i, n_i)

    v_i_ld = v_i.astype(np.longdouble)
    v_dot_n_ld = np.longdouble(v_dot_n)
    force_full += rho0_ld * v_i_ld * v_dot_n_ld * dA_ld

force_full = force_full.astype(np.float64)

print(f"\nQuadrature force: {force_full}")

# Compare to analytic
from slab.surface import force_incompressible_analytic
F_analytic = force_incompressible_analytic(a_idx, bodies, medium)

print(f"Analytic force:   {F_analytic}")
print(f"Ratio (quad/analytic): {np.linalg.norm(force_full) / np.linalg.norm(F_analytic):.2e}")

# Check what v_ext at the center should be
v_ext_center = v_ext_at(x_a, bodies, a_idx, rho0)
print(f"\nv_ext at Sun center: {v_ext_center}")
print(f"Expected force = rho0 * Q_a * v_ext = {rho0 * body_a.Q * v_ext_center}")
