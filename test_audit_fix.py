#!/usr/bin/env python3
"""Quick test to verify the audit fix works correctly."""

import numpy as np
from slab.bodies import Body
from slab.medium import Medium
from slab.surface import compare_force_methods

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

print("Testing audit with incompressible forces...")
result_inc = compare_force_methods(
    0, bodies, medium, n_points=256,
    use_compressible=False,
    verbose=True
)

print("\nTesting audit with compressible forces...")
result_comp = compare_force_methods(
    0, bodies, medium, n_points=256,
    use_compressible=True,
    verbose=True
)

print(f"\nIncompressible audit: {'PASS' if result_inc['passes_audit'] else 'FAIL'}")
print(f"Compressible audit: {'PASS' if result_comp['passes_audit'] else 'FAIL'}")
