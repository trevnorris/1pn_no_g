#!/usr/bin/env python3
"""Verify the magnitude ratio between v_self and v_ext."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.field import v_self, v_ext_at

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

# Pick a point on Sun's control surface
x_a = sun.x
R_a = sun.R
x_surface = x_a + np.array([R_a, 0.0, 0.0])  # Point on surface in +x direction

v_self_mag = np.linalg.norm(v_self(x_surface, x_a, sun.Q, medium.rho0))
v_ext_mag = np.linalg.norm(v_ext_at(x_surface, bodies, 0, medium.rho0))

print("="*70)
print("MAGNITUDE COMPARISON: v_self vs v_ext on Sun's control surface")
print("="*70)
print(f"\nSun control surface radius: R = {R_a}")
print(f"Point on surface: {x_surface}")
print(f"\n|v_self| = {v_self_mag:.6e}")
print(f"|v_ext| = {v_ext_mag:.6e}")
print(f"\nRatio: |v_self| / |v_ext| = {v_self_mag / v_ext_mag:.2e}")
print(f"\nThis means v_self is {v_self_mag / v_ext_mag:.0f}x larger than v_ext!")

# Estimate the self×self contribution
v_self_vec = v_self(x_surface, x_a, sun.Q, medium.rho0)
n = np.array([1.0, 0.0, 0.0])  # Normal at this point
v_self_dot_n = np.dot(v_self_vec, n)

self_self_integrand_mag = medium.rho0 * v_self_mag * abs(v_self_dot_n)
print(f"\nMagnitude of self×self integrand: ρ₀ |v_self| |v_self·n| = {self_self_integrand_mag:.6e}")

# The self×self integral should give zero, but needs ~10^10 precision to cancel
print(f"\nTo get 1% accuracy on the force ({1e-27:.0e}), we need cancellation")
print(f"of the self×self term ({self_self_integrand_mag:.0e}) to precision:")
print(f"  {self_self_integrand_mag / 1e-27:.0e} (about 10^{np.log10(self_self_integrand_mag / 1e-27):.0f})")
print(f"\nThis is IMPOSSIBLE with float64 precision (~10^-16 relative)!")
