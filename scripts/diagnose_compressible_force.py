#!/usr/bin/env python3
"""
Diagnostic script to trace compressible force calculation.

This will help identify where the factor of ~10^5 error is coming from
in the 1PN precession validation.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.field import v_ext_at, v_total, v_self
from slab.surface import force_incompressible_analytic, force_compressible_quadrature
from slab.geometry import fibonacci_sphere

# Create Mercury configuration at perihelion with e=0.1
rho0 = 1.0
beta0 = 0.045
cs = 63239.7263  # Speed of light [AU/yr]
a = 0.387  # Semi-major axis [AU]
e = 0.1

medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
K = medium.K

# Masses
M_sun = 1.0
M_mercury = 3.3e-7

# Perihelion setup
r_peri = a * (1 - e)
v_peri = np.sqrt(K * M_sun * (1 + e) / (a * (1 - e)))

sun = Body(
    name="Sun",
    M=M_sun,
    x=np.array([0.0, 0.0, 0.0]),
    v=np.array([0.0, 0.0, 0.0]),
    R=0.001,
    Q=M_sun / beta0,
)

mercury = Body(
    name="Mercury",
    M=M_mercury,
    x=np.array([r_peri, 0.0, 0.0]),
    v=np.array([0.0, v_peri, 0.0]),
    R=0.0005,
    Q=M_mercury / beta0,
)

bodies = [sun, mercury]

print("=" * 70)
print("COMPRESSIBLE FORCE DIAGNOSTIC")
print("=" * 70)
print()
print("Configuration:")
print(f"  cs = {cs:.4f} AU/yr (speed of light)")
print(f"  a = {a:.4f} AU")
print(f"  e = {e:.3f}")
print(f"  r_peri = {r_peri:.4f} AU")
print(f"  v_peri = {v_peri:.4f} AU/yr")
print(f"  Mercury v_body = {mercury.v}")
print(f"  K = {K:.6e}")
print()

# Compute external velocity at Mercury
v_ext = v_ext_at(mercury.x, bodies, 1, rho0)
v_ext_mag = np.linalg.norm(v_ext)

print("Velocity field analysis:")
print(f"  v_ext at Mercury = {v_ext}")
print(f"  |v_ext| = {v_ext_mag:.6e} AU/yr")
print(f"  v_ext / cs = {v_ext_mag / cs:.6e} (should be ~1.76e-4)")
print()

# Compute Mercury's velocity relative to the field
v_rel = v_ext - mercury.v
v_rel_mag = np.linalg.norm(v_rel)

print("Relative velocity (for thermodynamics):")
print(f"  v_rel = v_ext - v_body = {v_rel}")
print(f"  |v_rel| = {v_rel_mag:.6e} AU/yr")
print(f"  v_rel / cs = {v_rel_mag / cs:.6e}")
print(f"  (v_rel / cs)^2 = {(v_rel_mag / cs)**2:.6e}")
print()

# Expected density and pressure perturbations
Delta_rho = -rho0 * v_rel_mag**2 / (2 * cs**2)
P_star = -0.5 * rho0 * v_rel_mag**2
rho_star = rho0 + Delta_rho

print("Thermodynamic perturbations:")
print(f"  Δρ/ρ₀ = {Delta_rho / rho0:.6e} (should be ~Ma^2 ~ 3e-8)")
print(f"  P* = {P_star:.6e}")
print(f"  ρ* = {rho_star:.6e}")
print()

# Compute incompressible force
F_inc = force_incompressible_analytic(1, bodies, medium)
F_inc_mag = np.linalg.norm(F_inc)

print("Incompressible force:")
print(f"  F_inc = {F_inc}")
print(f"  |F_inc| = {F_inc_mag:.6e}")
print()

# Compute compressible correction
F_comp = force_compressible_quadrature(1, bodies, medium, n_points=512)
F_comp_mag = np.linalg.norm(F_comp)

print("Compressible correction:")
print(f"  F_comp = {F_comp}")
print(f"  |F_comp| = {F_comp_mag:.6e}")
print(f"  |F_comp| / |F_inc| = {F_comp_mag / F_inc_mag:.6e}")
print(f"  Expected ratio ~ Ma^2 = {(v_ext_mag / cs)**2:.6e}")
print()

# Manual integrand inspection
print("=" * 70)
print("SURFACE INTEGRAND INSPECTION")
print("=" * 70)
print()

# Sample one point on the control surface
normals = fibonacci_sphere(8)  # Just a few points for inspection
R_a = mercury.R

print(f"Sampling {len(normals)} points on Mercury's control surface (R={R_a} AU)")
print()

for i in range(min(4, len(normals))):
    n_i = normals[i]
    x_i = mercury.x + R_a * n_i

    print(f"Point {i}: normal = {n_i}")

    # Compute velocities at this surface point
    v_ext_i = v_ext_at(x_i, bodies, 1, rho0)
    v_total_i = v_total(x_i, bodies, rho0)
    v_self_i = v_self(x_i, mercury.x, mercury.Q, rho0)

    print(f"  v_ext(x_i) = {v_ext_i}")
    print(f"  |v_ext(x_i)| = {np.linalg.norm(v_ext_i):.6e} AU/yr")

    # Boost to body frame
    v_rel_i = v_ext_i - mercury.v
    v_rel_mag_sq_i = np.dot(v_rel_i, v_rel_i)

    print(f"  v_rel(x_i) = v_ext - v_body = {v_rel_i}")
    print(f"  |v_rel|^2 = {v_rel_mag_sq_i:.6e} (AU/yr)^2")

    # Thermodynamic quantities
    Delta_rho_i = -rho0 * v_rel_mag_sq_i / (2 * cs**2)
    P_star_i = -0.5 * rho0 * v_rel_mag_sq_i

    print(f"  Δρ = {Delta_rho_i:.6e}")
    print(f"  P* = {P_star_i:.6e}")

    # Momentum flux components
    v_total_i_boosted = v_total_i - mercury.v
    v_dot_n = np.dot(v_total_i_boosted, n_i)

    momentum_correction = Delta_rho_i * v_total_i_boosted * v_dot_n

    print(f"  v_total(x_i) boosted = {v_total_i_boosted}")
    print(f"  v·n = {v_dot_n:.6e}")
    print(f"  Δρ * v * (v·n) = {momentum_correction}")
    print(f"  |momentum term| = {np.linalg.norm(momentum_correction):.6e}")

    print()

print("=" * 70)
print("DIMENSIONAL ANALYSIS")
print("=" * 70)
print()

# Check if cs has the wrong units or scale factor
print("Checking for dimensional issues:")
print(f"  cs (provided) = {cs:.4f} AU/yr")
print(f"  c (physical) = 299792 km/s = 63239.7263 AU/yr")
print(f"  v_peri (computed) = {v_peri:.4f} AU/yr")
print(f"  v_peri / cs = {v_peri / cs:.6e} (Mach number)")
print()

# Check force dimensions
print("Force dimensions:")
print(f"  rho0 = {rho0} (dimensionless)")
print(f"  Q_mercury = {mercury.Q:.6e} (mass units / beta0)")
print(f"  v_ext ~ {v_ext_mag:.6e} AU/yr")
print(f"  F_inc ~ rho0 * Q * v_ext = {rho0 * mercury.Q * v_ext_mag:.6e}")
print(f"  |F_inc| (actual) = {F_inc_mag:.6e}")
print()

# Estimate what compressible force SHOULD be
expected_F_comp = F_inc_mag * (v_peri / cs)**2
print(f"Expected |F_comp| ~ |F_inc| * (v/cs)^2:")
print(f"  = {F_inc_mag:.6e} * {(v_peri/cs)**2:.6e}")
print(f"  = {expected_F_comp:.6e}")
print(f"Actual |F_comp| = {F_comp_mag:.6e}")
print(f"Ratio (actual/expected) = {F_comp_mag / expected_F_comp:.6e}")
print()

if F_comp_mag / expected_F_comp > 100:
    print("ERROR: Compressible force is ~{:.0f}× too large!".format(
        F_comp_mag / expected_F_comp))
    print()
    print("Likely causes:")
    print("  1. v_rel = v_ext - v_body is being computed incorrectly")
    print("  2. cs is in wrong units or has wrong value")
    print("  3. Boost to body frame is missing or incorrect")
    print("  4. Integrand has sign error or extra factor")
