#!/usr/bin/env python3
"""Debug script to understand why audit is failing."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.surface import (
    force_incompressible_analytic,
    force_incompressible_quadrature,
    force_compressible_analytic,
    force_compressible_quadrature,
    force_total,
)

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
print("DEBUGGING AUDIT FAILURE")
print("="*70)
print(f"\nSun: M={sun.M}, Q={sun.Q}, R={sun.R}, x={sun.x}")
print(f"Mercury: M={mercury.M}, Q={mercury.Q}, R={mercury.R}, x={mercury.x}")
print(f"\nMedium: rho0={medium.rho0}, cs={medium.cs}, beta0={medium.beta0}")

print("\n" + "="*70)
print("INCOMPRESSIBLE FORCES (body 0 = Sun)")
print("="*70)

# Direct calls to force functions
F_inc_analytic = force_incompressible_analytic(0, bodies, medium)
F_inc_quadrature = force_incompressible_quadrature(0, bodies, medium, n_points=256)

print(f"\nIncompressible Analytic:   {F_inc_analytic}")
print(f"Incompressible Quadrature: {F_inc_quadrature}")
print(f"Difference:                {F_inc_quadrature - F_inc_analytic}")
print(f"Ratio (quad/analytic):     {np.linalg.norm(F_inc_quadrature) / np.linalg.norm(F_inc_analytic):.2e}")

print("\n" + "="*70)
print("COMPRESSIBLE FORCES (body 0 = Sun)")
print("="*70)

F_comp_analytic = force_compressible_analytic(0, bodies, medium, n_points=256)
F_comp_quadrature = force_compressible_quadrature(0, bodies, medium, n_points=256)

print(f"\nCompressible Analytic:   {F_comp_analytic}")
print(f"Compressible Quadrature: {F_comp_quadrature}")
print(f"Difference:              {F_comp_quadrature - F_comp_analytic}")

print("\n" + "="*70)
print("FORCE_TOTAL WITH use_compressible=True")
print("="*70)

F_total_analytic = force_total(
    a_idx=0,
    bodies=bodies,
    medium=medium,
    use_compressible=True,
    use_quadrature=False,
    n_points=256,
)

F_total_quadrature = force_total(
    a_idx=0,
    bodies=bodies,
    medium=medium,
    use_compressible=True,
    use_quadrature=True,
    n_points=256,
)

print(f"\nforce_total (analytic, use_compressible=True):   {F_total_analytic}")
print(f"force_total (quadrature, use_compressible=True): {F_total_quadrature}")
print(f"Difference:                                       {F_total_quadrature - F_total_analytic}")
print(f"Ratio (quad/analytic):                            {np.linalg.norm(F_total_quadrature) / np.linalg.norm(F_total_analytic):.2e}")

print("\n" + "="*70)
print("FORCE_TOTAL WITH use_compressible=False")
print("="*70)

F_total_analytic_inc = force_total(
    a_idx=0,
    bodies=bodies,
    medium=medium,
    use_compressible=False,
    use_quadrature=False,
    n_points=256,
)

F_total_quadrature_inc = force_total(
    a_idx=0,
    bodies=bodies,
    medium=medium,
    use_compressible=False,
    use_quadrature=True,
    n_points=256,
)

print(f"\nforce_total (analytic, use_compressible=False):   {F_total_analytic_inc}")
print(f"force_total (quadrature, use_compressible=False): {F_total_quadrature_inc}")
print(f"Difference:                                        {F_total_quadrature_inc - F_total_analytic_inc}")
print(f"Ratio (quad/analytic):                             {np.linalg.norm(F_total_quadrature_inc) / np.linalg.norm(F_total_analytic_inc):.2e}")

print("\n" + "="*70)
print("EXPECTED FORCE MAGNITUDE")
print("="*70)

# Calculate expected force using F = K * M1 * M2 / r^2
r = np.linalg.norm(mercury.x - sun.x)
K = medium.K
F_expected = K * sun.M * mercury.M / r**2

print(f"\nSeparation r = {r:.3f} AU")
print(f"K = rho0/(4*pi*beta0^2) = {K:.6e}")
print(f"Expected force magnitude: F = K * M1 * M2 / r^2 = {F_expected:.6e}")
print(f"\nCompare to:")
print(f"  Analytic incompressible:   {np.linalg.norm(F_inc_analytic):.6e}")
print(f"  Quadrature incompressible: {np.linalg.norm(F_inc_quadrature):.6e}")
