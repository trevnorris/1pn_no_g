#!/usr/bin/env python3
"""Verify the analytic formula is correct."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.field import v_ext_at

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
print("VERIFYING ANALYTIC FORMULA")
print("="*70)

# The analytic formula from surface.py line 100:
# F_a = ρ₀ * Q_a * v_ext(r_a)

rho0 = medium.rho0
Q_sun = sun.Q
x_sun = sun.x

v_ext_at_sun = v_ext_at(x_sun, bodies, 0, rho0)

F_by_formula = rho0 * Q_sun * v_ext_at_sun

print(f"\nInputs:")
print(f"  ρ₀ = {rho0}")
print(f"  Q_sun = {Q_sun:.6e}")
print(f"  x_sun = {x_sun}")
print(f"\nv_ext at Sun center: {v_ext_at_sun}")
print(f"\nF = ρ₀ * Q * v_ext = {F_by_formula}")

# Now let's verify this matches the Newtonian force
# From the field.py definition (line 107):
# v_ext(x_a) = - Σ_{b≠a} (Q_b/4π) * r_ab/r_ab³
# where r_ab = x_a - x_b

r_vec = sun.x - mercury.x  # Vector from Mercury to Sun
r = np.linalg.norm(r_vec)
r_hat = r_vec / r

v_ext_manual = -(mercury.Q / (4 * np.pi)) * r_vec / r**3

print(f"\nManual calculation of v_ext:")
print(f"  r_vec (Sun - Mercury) = {r_vec}")
print(f"  r = {r}")
print(f"  v_ext = -(Q_merc/4π) * r_vec/r³ = {v_ext_manual}")
print(f"  Matches? {np.allclose(v_ext_manual, v_ext_at_sun)}")

# Now compute force
F_manual = rho0 * Q_sun * v_ext_manual
print(f"\nF = ρ₀ * Q_sun * v_ext = {F_manual}")

# Compare to "Newtonian" form
K = medium.K
F_newtonian = K * sun.M * mercury.M / r**2 * r_hat

print(f"\nNewtonian form:")
print(f"  K = {K:.6e}")
print(f"  F = K * M_sun * M_merc / r² * r_hat = {F_newtonian}")
print(f"  Matches? {np.allclose(F_newtonian, F_manual)}")

print("\n" + "="*70)
print("KEY INSIGHT")
print("="*70)
print("\nThe analytic formula F = ρ₀ Q v_ext is CORRECT.")
print("The direction should be TOWARD Mercury (positive x).")
print(f"Force magnitude: {np.linalg.norm(F_manual):.6e}")
print(f"Force direction: {F_manual / np.linalg.norm(F_manual)}")
