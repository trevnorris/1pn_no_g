#!/usr/bin/env python3
"""Test that the self×ext cross-term gives the correct force."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.surface import decompose_momentum_integrand, force_incompressible_analytic

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
print("TESTING CROSS-TERM FIX")
print("="*70)

# Get decomposition
decomp = decompose_momentum_integrand(0, bodies, medium, n_points=256)

print("\nDecomposed momentum flux integrals:")
print(f"  F_total (using v_total):  {decomp['total']}")
print(f"  F_self_ext (cross-term):  {decomp['self_ext']}")
print(f"  F_self_self:              {decomp['self_self']}")
print(f"  F_ext_ext:                {decomp['ext_ext']}")

# Get analytic force
F_analytic = force_incompressible_analytic(0, bodies, medium)
print(f"\nAnalytic force:             {F_analytic}")

# Compare
print("\n" + "="*70)
print("COMPARISON TO ANALYTIC")
print("="*70)

error_total = np.linalg.norm(decomp['total'] - F_analytic) / np.linalg.norm(F_analytic)
error_cross = np.linalg.norm(decomp['self_ext'] - F_analytic) / np.linalg.norm(F_analytic)

print(f"\nF_total vs analytic:")
print(f"  Relative error: {error_total:.6e}")
print(f"  Passes audit (< 1e-3)? {error_total < 1e-3}")

print(f"\nF_self_ext vs analytic:")
print(f"  Relative error: {error_cross:.6e}")
print(f"  Passes audit (< 1e-3)? {error_cross < 1e-3}")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if error_cross < 1e-3:
    print("\n✓ The self×ext CROSS-TERM gives the correct force!")
    print("✓ The fix is to use ONLY the cross-term in the quadrature.")
else:
    print("\n✗ Something else is wrong...")

if decomp['total'][0] / F_analytic[0] > 100:
    print(f"\n✗ F_total is {decomp['total'][0] / F_analytic[0]:.0f}x too large!")
    print("✗ This is because of catastrophic cancellation in self×self term.")
