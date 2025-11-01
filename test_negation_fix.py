#!/usr/bin/env python3
"""Test if negating the quadrature result fixes the sign."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.surface import force_incompressible_analytic, force_incompressible_quadrature

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
print("TESTING NEGATION FIX")
print("="*70)

F_analytic = force_incompressible_analytic(0, bodies, medium)
F_quadrature = force_incompressible_quadrature(0, bodies, medium, n_points=256)
F_quadrature_negated = -F_quadrature

print(f"\nAnalytic force:           {F_analytic}")
print(f"Quadrature (as-is):       {F_quadrature}")
print(f"Quadrature (negated):     {F_quadrature_negated}")

error_as_is = np.linalg.norm(F_quadrature - F_analytic) / np.linalg.norm(F_analytic)
error_negated = np.linalg.norm(F_quadrature_negated - F_analytic) / np.linalg.norm(F_analytic)

print(f"\nRelative error (as-is):   {error_as_is:.6e}")
print(f"Relative error (negated): {error_negated:.6e}")

print(f"\nDoes negation fix it? {error_negated < 1e-3}")

if error_negated < 1e-3:
    print("\n✓ YES! Negating the quadrature result gives the correct answer!")
    print("✓ The bug is a SIGN ERROR in force_incompressible_quadrature")
else:
    print("\n✗ NO. The problem is more complex than just a sign error.")
