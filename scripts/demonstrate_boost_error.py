#!/usr/bin/env python3
"""
Demonstrate the error in the velocity boost application.

The compressible force should be the DIFFERENCE between:
  F_comp_full = ∫ ρ* v(v·n) - (P* + 0.5 ρ* v_rel²) n dA
  F_incomp = ∫ ρ₀ v(v·n) dA

Where v is the velocity in the LAB frame (not boosted).

But the current implementation appears to boost v in the momentum term,
which creates an artificial large force.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.field import v_ext_at, v_total, v_self
from slab.geometry import fibonacci_sphere

# Create Mercury configuration
rho0 = 1.0
beta0 = 0.045
cs = 63239.7263  # AU/yr
a = 0.387
e = 0.1

medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
K = medium.K

M_sun = 1.0
M_mercury = 3.3e-7

r_peri = a * (1 - e)
v_peri = np.sqrt(K * M_sun * (1 + e) / (a * (1 - e)))

sun = Body(
    name="Sun", M=M_sun, x=np.array([0.0, 0.0, 0.0]),
    v=np.array([0.0, 0.0, 0.0]), R=0.001, Q=M_sun / beta0,
)

mercury = Body(
    name="Mercury", M=M_mercury,
    x=np.array([r_peri, 0.0, 0.0]),
    v=np.array([0.0, v_peri, 0.0]),
    R=0.0005, Q=M_mercury / beta0,
)

bodies = [sun, mercury]

print("=" * 70)
print("DEMONSTRATING VELOCITY BOOST ERROR")
print("=" * 70)
print()

# Sample the control surface
n_points = 512
normals = fibonacci_sphere(n_points)
R_a = mercury.R
dA = 4.0 * np.pi * R_a**2 / n_points

v_body = mercury.v
print(f"Mercury velocity: v_body = {v_body}")
print(f"|v_body| = {np.linalg.norm(v_body):.4f} AU/yr")
print()

# Compute forces using different prescriptions
F_correct = np.zeros(3, dtype=np.longdouble)
F_wrong = np.zeros(3, dtype=np.longdouble)
F_incomp_check = np.zeros(3, dtype=np.longdouble)

rho0_ld = np.longdouble(rho0)
cs_ld = np.longdouble(cs)
dA_ld = np.longdouble(dA)

for i in range(n_points):
    n_i = normals[i]
    x_i = mercury.x + R_a * n_i

    # Velocities at surface point
    v_ext_i = v_ext_at(x_i, bodies, 1, rho0)
    v_total_i = v_total(x_i, bodies, rho0)

    # Boost to body frame FOR THERMODYNAMICS ONLY
    v_rel_i = v_ext_i - v_body
    v_rel_mag_sq_i = np.dot(v_rel_i, v_rel_i)

    # Thermodynamic quantities (using v_rel)
    Delta_rho_i = -rho0 * v_rel_mag_sq_i / (2 * cs**2)
    P_star_i = -0.5 * rho0 * v_rel_mag_sq_i
    rho_star_i = rho0 + Delta_rho_i

    # CORRECT: Use v_total (LAB frame) for momentum flux
    v_dot_n_correct = np.dot(v_total_i, n_i)

    momentum_correction_correct = Delta_rho_i * v_total_i * v_dot_n_correct
    bracket_correct = P_star_i + 0.5 * rho_star_i * v_rel_mag_sq_i
    pressure_term_correct = -bracket_correct * n_i

    F_correct += (momentum_correction_correct + pressure_term_correct) * dA

    # WRONG: Boost v_total before using in momentum flux
    v_total_boosted = v_total_i - v_body
    v_dot_n_wrong = np.dot(v_total_boosted, n_i)

    momentum_correction_wrong = Delta_rho_i * v_total_boosted * v_dot_n_wrong
    bracket_wrong = P_star_i + 0.5 * rho_star_i * v_rel_mag_sq_i
    pressure_term_wrong = -bracket_wrong * n_i

    F_wrong += (momentum_correction_wrong + pressure_term_wrong) * dA

    # Also check incompressible for reference
    F_incomp_check += rho0 * v_total_i * v_dot_n_correct * dA

print("RESULTS:")
print("=" * 70)
print()

F_correct_64 = F_correct.astype(np.float64)
F_wrong_64 = F_wrong.astype(np.float64)
F_incomp_check_64 = F_incomp_check.astype(np.float64)

print("Incompressible force (for reference):")
print(f"  F_inc = {F_incomp_check_64}")
print(f"  |F_inc| = {np.linalg.norm(F_incomp_check_64):.6e}")
print()

print("Compressible correction (CORRECT - no boost in momentum):")
print(f"  F_comp = {F_correct_64}")
print(f"  |F_comp| = {np.linalg.norm(F_correct_64):.6e}")
print(f"  |F_comp| / |F_inc| = {np.linalg.norm(F_correct_64) / np.linalg.norm(F_incomp_check_64):.6e}")
print()

print("Compressible correction (WRONG - boost in momentum):")
print(f"  F_comp = {F_wrong_64}")
print(f"  |F_comp| = {np.linalg.norm(F_wrong_64):.6e}")
print(f"  |F_comp| / |F_inc| = {np.linalg.norm(F_wrong_64) / np.linalg.norm(F_incomp_check_64):.6e}")
print()

ratio = np.linalg.norm(F_wrong_64) / np.linalg.norm(F_correct_64)
print(f"Ratio (wrong / correct) = {ratio:.2e}")
print()

print("=" * 70)
print("ANALYSIS:")
print("=" * 70)
print()

if ratio > 1000:
    print(f"The WRONG method gives a force {ratio:.0e}× larger!")
    print()
    print("This explains the catastrophic precession error.")
    print()
    print("The issue: boosting v in the momentum flux creates")
    print("a term ~ Δρ * v_body * (v·n) which scales as:")
    print(f"  ~ (v_ext²/cs²) * v_body * v_ext")
    print(f"  ~ v_body * v_ext³ / cs²")
    print()
    print(f"But v_body ~ v_ext ~ v_orbital ~ {v_peri:.1f} AU/yr")
    print(f"So this gives a huge spurious force!")
    print()
    print("FIX: Remove the '- v_body' from line 791 in surface.py")
    print("     v_total_i = v_total(x_i, bodies, rho0)  # NO boost here!")
elif 0.5 < ratio < 2:
    print(f"The two methods give similar results (ratio {ratio:.2f})")
    print("The boost in momentum flux is not the main issue.")
else:
    print(f"Unexpected ratio: {ratio:.2e}")
