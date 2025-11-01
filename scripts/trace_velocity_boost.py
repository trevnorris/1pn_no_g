#!/usr/bin/env python3
"""
Trace how the velocity boost (v_ext - v_body) changes around the orbit.

The huge precession suggests the compressible force may be incorrectly
applying the body velocity boost, causing forces that scale with v_body^2
instead of v_rel^2.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.field import v_ext_at

# Create Mercury configuration at various orbital phases
rho0 = 1.0
beta0 = 0.045
cs = 63239.7263  # Speed of light [AU/yr]
a = 0.387
e = 0.1

medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
K = medium.K

M_sun = 1.0
M_mercury = 3.3e-7

# Test at perihelion and aphelion
r_peri = a * (1 - e)
r_apo = a * (1 + e)

v_peri = np.sqrt(K * M_sun * (1 + e) / (a * (1 - e)))
v_apo = np.sqrt(K * M_sun * (1 - e) / (a * (1 + e)))

print("=" * 70)
print("VELOCITY BOOST ANALYSIS AROUND ORBIT")
print("=" * 70)
print()
print(f"cs = {cs:.4f} AU/yr")
print(f"e = {e:.3f}")
print()

# Perihelion (x = r_peri, v = [0, v_peri, 0])
print("PERIHELION:")
print(f"  r = {r_peri:.4f} AU")
print(f"  v_body = [0, {v_peri:.4f}, 0] AU/yr")
print(f"  |v_body| = {v_peri:.4f} AU/yr")

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

v_ext = v_ext_at(mercury.x, bodies, 1, rho0)
v_ext_mag = np.linalg.norm(v_ext)
v_rel = v_ext - mercury.v
v_rel_mag = np.linalg.norm(v_rel)

print(f"  v_ext = {v_ext}")
print(f"  |v_ext| = {v_ext_mag:.4f} AU/yr")
print(f"  v_rel = v_ext - v_body = {v_rel}")
print(f"  |v_rel| = {v_rel_mag:.4f} AU/yr")
print(f"  |v_rel|^2 = {v_rel_mag**2:.4f} (AU/yr)^2")
print(f"  |v_rel|^2 / cs^2 = {(v_rel_mag / cs)**2:.6e}")
print()

# Aphelion (x = -r_apo, v = [0, -v_apo, 0])
print("APHELION:")
print(f"  r = {r_apo:.4f} AU")
print(f"  v_body = [0, {-v_apo:.4f}, 0] AU/yr")
print(f"  |v_body| = {v_apo:.4f} AU/yr")

mercury.x = np.array([-r_apo, 0.0, 0.0])
mercury.v = np.array([0.0, -v_apo, 0.0])

v_ext = v_ext_at(mercury.x, bodies, 1, rho0)
v_ext_mag = np.linalg.norm(v_ext)
v_rel = v_ext - mercury.v
v_rel_mag = np.linalg.norm(v_rel)

print(f"  v_ext = {v_ext}")
print(f"  |v_ext| = {v_ext_mag:.4f} AU/yr")
print(f"  v_rel = v_ext - v_body = {v_rel}")
print(f"  |v_rel| = {v_rel_mag:.4f} AU/yr")
print(f"  |v_rel|^2 = {v_rel_mag**2:.4f} (AU/yr)^2")
print(f"  |v_rel|^2 / cs^2 = {(v_rel_mag / cs)**2:.6e}")
print()

print("=" * 70)
print("KEY OBSERVATION:")
print("=" * 70)
print()
print("If the compressible force is missing the velocity boost:")
print(f"  - It would use |v_ext|^2 / cs^2 ~ {(v_ext_mag/cs)**2:.6e}")
print(f"  - Instead of |v_rel|^2 / cs^2 ~ {(v_rel_mag/cs)**2:.6e}")
print(f"  - Ratio: {(v_rel_mag/v_ext_mag)**2:.2f}×")
print()
print("But we see |v_rel| ~ |v_ext + v_body| >> |v_ext|")
print("because v_ext and v_body are nearly perpendicular.")
print()
print("This means:")
print("  |v_rel|^2 ~ v_ext^2 + v_body^2 ~ v_body^2 (since v_body > v_ext)")
print(f"  v_body^2 = {v_peri**2:.4f} (AU/yr)^2")
print(f"  v_ext^2 = {v_ext_mag**2:.4f} (AU/yr)^2")
print(f"  v_rel^2 = {v_rel_mag**2:.4f} (AU/yr)^2")
print()
print("Ratio v_body^2 / v_ext^2 = {:.2f}".format(v_peri**2 / v_ext_mag**2))
print()
print("HYPOTHESIS: The boost is CORRECT, but it's making the force LARGER")
print("because |v_rel|^2 > |v_ext|^2 when v_body is significant.")
print()

# Now check what happens with NO boost (old behavior before fix)
print("=" * 70)
print("COMPARISON: WITH vs WITHOUT BOOST")
print("=" * 70)
print()

mercury.x = np.array([r_peri, 0.0, 0.0])
mercury.v = np.array([0.0, v_peri, 0.0])
v_ext = v_ext_at(mercury.x, bodies, 1, rho0)

print("At perihelion:")
print(f"  WITHOUT boost: |v_thermodynamic|^2 = |v_ext|^2 = {np.dot(v_ext, v_ext):.4f}")
print(f"  WITH boost:    |v_thermodynamic|^2 = |v_rel|^2 = {np.dot(v_rel, v_rel):.4f}")
print(f"  Ratio: {np.dot(v_rel, v_rel) / np.dot(v_ext, v_ext):.2f}×")
print()

factor = np.dot(v_rel, v_rel) / np.dot(v_ext, v_ext)
print(f"If the boost was ADDED as part of the fix, the compressible force")
print(f"would be ~{factor:.1f}× larger than before!")
print()
print(f"Expected precession increase: ~{factor:.1f}×")
print(f"Observed: from ~0.1 rad/orbit to ~0.15 rad/orbit")
print(f"That's a factor of ~{0.15 / 0.1:.1f}×")
print()

# Check what cs value would be needed to get the right precession
print("=" * 70)
print("WHAT cs VALUE WOULD FIX THE PRECESSION?")
print("=" * 70)
print()

# GR prediction at e=0.1
GR_precession = (6 * np.pi * K * M_sun) / (a * cs**2 * (1 - e**2))
observed_precession = -0.151  # rad/orbit from log

print(f"GR prediction: {GR_precession:.6e} rad/orbit")
print(f"Observed: {observed_precession:.6e} rad/orbit")
print(f"Ratio: {abs(observed_precession / GR_precession):.0f}×")
print()

# If precession ~ 1/cs^2, what cs would give the observed value?
cs_needed = cs * np.sqrt(abs(observed_precession / GR_precession))
print(f"cs needed for correct precession: {cs_needed:.4f} AU/yr")
print(f"Current cs: {cs:.4f} AU/yr")
print(f"Ratio: {cs_needed / cs:.2f}×")
print()

# Check if this corresponds to using v_ext instead of v_rel
cs_eff_with_boost = cs * np.sqrt(factor)
print(f"If we're using |v_rel|^2 instead of |v_ext|^2:")
print(f"  Effective cs^2 is reduced by factor of {factor:.2f}")
print(f"  Effective cs = {cs / np.sqrt(factor):.4f} AU/yr")
print(f"  This would increase precession by factor of {factor:.2f}")
