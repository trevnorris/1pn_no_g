#!/usr/bin/env python3
"""
Estimate what the 1PN force magnitude should be to produce GR precession.

The GR perihelion precession is:
    Δω = (6π G M) / (a c² (1 - e²))

This should arise from velocity-dependent corrections to the force of order v²/c².
Let's check what force magnitude is needed and compare to what we're computing.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

# Mercury orbit parameters
a = 0.387  # AU
e = 0.1
M_sun = 1.0
beta0 = 0.045
rho0 = 1.0
cs = 63239.7263  # AU/yr

# Compute K
K = rho0 / (4 * np.pi * beta0**2)

# Orbital parameters
r_peri = a * (1 - e)
v_peri = np.sqrt(K * M_sun * (1 + e) / (a * (1 - e)))
T_orbit = 2 * np.pi * np.sqrt(a**3 / (K * M_sun))

print("=" * 70)
print("ESTIMATING REQUIRED 1PN FORCE MAGNITUDE")
print("=" * 70)
print()
print("Orbital parameters:")
print(f"  a = {a:.4f} AU")
print(f"  e = {e:.3f}")
print(f"  r_peri = {r_peri:.4f} AU")
print(f"  v_peri = {v_peri:.4f} AU/yr")
print(f"  T_orbit = {T_orbit:.4f} yr")
print(f"  K = {K:.6e}")
print(f"  cs = {cs:.4f} AU/yr")
print()

# GR precession formula
GR_precession_per_orbit = (6 * np.pi * K * M_sun) / (a * cs**2 * (1 - e**2))
GR_precession_rate = GR_precession_per_orbit / T_orbit  # rad/yr

print("GR predictions:")
print(f"  Δω_GR = {GR_precession_per_orbit:.6e} rad/orbit")
print(f"  Δω_GR = {GR_precession_per_orbit * 206265:.6e} arcsec/orbit")
print(f"  dω/dt = {GR_precession_rate:.6e} rad/yr")
print()

# Newtonian force at perihelion
F_newtonian = K * M_sun * 3.3e-7 / r_peri**2
print(f"Newtonian force at perihelion: F_N ~ {F_newtonian:.6e}")
print()

# The 1PN force correction should be ~ (v²/c²) F_N
F_1pn_expected = (v_peri / cs)**2 * F_newtonian
print(f"Expected 1PN force: F_1PN ~ (v/cs)² F_N")
print(f"  = ({v_peri:.4f} / {cs:.4f})² × {F_newtonian:.6e}")
print(f"  = {(v_peri/cs)**2:.6e} × {F_newtonian:.6e}")
print(f"  = {F_1pn_expected:.6e}")
print()

# Typical orbit-averaged force needed to produce precession
# The precession comes from the azimuthal component of the 1PN force
# For a rotating orbit, Δω ~ (1/L) ∫ τ dt where τ is the torque
# This is roughly: Δω ~ (r × F_1PN) / L ~ (r F_1PN) / (m v r) ~ F_1PN / (m v)

# Angular momentum
M_mercury = 3.3e-7
L = M_mercury * np.sqrt(K * M_sun * a * (1 - e**2))
print(f"Angular momentum: L = {L:.6e}")
print()

# Characteristic torque needed
# Δω = ∫(dω/dt) dt = ∫(τ/L) dt = (1/L) ∫ τ dt
# Over one orbit: Δω ~ (τ_avg × T) / L
tau_needed = GR_precession_per_orbit * L / T_orbit
print(f"Torque needed: τ = (Δω × L) / T")
print(f"  = {GR_precession_per_orbit:.6e} × {L:.6e} / {T_orbit:.4f}")
print(f"  = {tau_needed:.6e}")
print()

# Torque ~ r × F, so F ~ τ / r
F_from_torque = tau_needed / r_peri
print(f"Force from torque estimate: F ~ τ / r")
print(f"  = {tau_needed:.6e} / {r_peri:.4f}")
print(f"  = {F_from_torque:.6e}")
print()

print("=" * 70)
print("COMPARISON WITH ACTUAL COMPRESSIBLE FORCE")
print("=" * 70)
print()

# From the diagnostic script, we found:
F_comp_actual = 7.5e-12
print(f"Actual compressible force: {F_comp_actual:.6e}")
print(f"Expected from (v/cs)² scaling: {F_1pn_expected:.6e}")
print(f"Expected from torque balance: {F_from_torque:.6e}")
print()

print("Ratio (actual / expected):")
print(f"  vs (v/cs)² estimate: {F_comp_actual / F_1pn_expected:.2e}")
print(f"  vs torque estimate: {F_comp_actual / F_from_torque:.2e}")
print()

# The actual precession observed
precession_observed = 0.15  # rad/orbit (magnitude)
print(f"Observed precession: {precession_observed:.6e} rad/orbit")
print(f"GR prediction: {GR_precession_per_orbit:.6e} rad/orbit")
print(f"Ratio: {precession_observed / GR_precession_per_orbit:.2e}×")
print()

# What force would produce the observed precession?
F_for_observed = F_from_torque * (precession_observed / GR_precession_per_orbit)
print(f"Force needed for observed precession:")
print(f"  = {F_from_torque:.6e} × {precession_observed / GR_precession_per_orbit:.2e}")
print(f"  = {F_for_observed:.6e}")
print()

print("=" * 70)
print("DIAGNOSIS:")
print("=" * 70)
print()

ratio_force = F_comp_actual / F_from_torque
ratio_precession = precession_observed / GR_precession_per_orbit

print(f"Force ratio (actual/expected): {ratio_force:.2e}")
print(f"Precession ratio (observed/GR): {ratio_precession:.2e}")
print()

if abs(ratio_force - 1) < 10 and ratio_precession > 1000:
    print("CONCLUSION: Force magnitude is CORRECT!")
    print("The problem is NOT in the force calculation.")
    print()
    print("Possible issues:")
    print("  1. Precession measurement error (omega unwrapping, fitting)")
    print("  2. Initial conditions not exactly at perihelion")
    print("  3. Integration errors accumulating")
    print("  4. cs value in wrong units somewhere")
    print("  5. Force direction issue (causing large, non-precessing torques)")
else:
    print(f"The compressible force is ~{ratio_force:.0e}× the expected magnitude.")
    print(f"But the precession is ~{ratio_precession:.0e}× too large.")
    print()
    print("This suggests a mismatch between force and precession,")
    print("possibly in how the force is being integrated or applied.")
