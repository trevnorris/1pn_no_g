"""Minimal script to debug force calculation formula."""
import numpy as np
from slab.medium import Medium
from slab.bodies import Body
from slab.field import v_ext_at
from slab.surface import force_incompressible_analytic

# Simple parameters
rho0 = 1.0
beta0 = 1.0e10
cs = 1.0e4

# Create medium
medium = Medium(rho0=rho0, cs=cs, beta0=beta0)
print(f"K = {medium.K:.6e}")

# Two bodies separated by r=1.0
M1 = 1.0
M2 = 3.0e-6
r = 1.0

Q1 = M1 / beta0
Q2 = M2 / beta0

print(f"\nBody 1: M={M1}, Q={Q1:.6e}")
print(f"Body 2: M={M2}, Q={Q2:.6e}")
print(f"Separation: r={r}")

# Create bodies
body1 = Body("Body1", M=M1, x=np.array([0.0, 0.0, 0.0]), v=np.array([0.0, 0.0, 0.0]), R=1e-3, Q=Q1)
body2 = Body("Body2", M=M2, x=np.array([r, 0.0, 0.0]), v=np.array([0.0, 0.0, 0.0]), R=5e-4, Q=Q2)

bodies = [body1, body2]

# Compute v_ext at body2 due to body1
v_ext = v_ext_at(body2.x, bodies, 1, rho0)
print(f"\nv_ext at body2: {v_ext}")
print(f"|v_ext| = {np.linalg.norm(v_ext):.6e}")

# Compute force on body2 using current formula
F2 = force_incompressible_analytic(1, bodies, medium)
print(f"\nForce on body2 (current formula): {F2}")
print(f"|F2| = {np.linalg.norm(F2):.6e}")

# Expected force from Newton-like law
F_expected = medium.K * M1 * M2 / r**2
print(f"\nExpected force magnitude: {F_expected:.6e}")

# Test various formulas
print("\n" + "="*70)
print("TESTING DIFFERENT FORMULAS:")
print("="*70)

# Formula 1: Current (4/3) * (Q/4π) * v_ext
F1 = (4.0/3.0) * (Q2/(4.0*np.pi)) * v_ext
print(f"1. F = (4/3) * (Q/4π) * v_ext:     {np.linalg.norm(F1):.6e}  (ratio: {F_expected/np.linalg.norm(F1):.2f})")

# Formula 2: (4/3) * Q * v_ext
F2 = (4.0/3.0) * Q2 * v_ext
print(f"2. F = (4/3) * Q * v_ext:          {np.linalg.norm(F2):.6e}  (ratio: {F_expected/np.linalg.norm(F2):.2f})")

# Formula 3: ρ₀ * (4/3) * (Q/4π) * v_ext
F3 = rho0 * (4.0/3.0) * (Q2/(4.0*np.pi)) * v_ext
print(f"3. F = ρ₀ * (4/3) * (Q/4π) * v_ext: {np.linalg.norm(F3):.6e}  (ratio: {F_expected/np.linalg.norm(F3):.2f})")

# Formula 4: ρ₀ * Q * v_ext
F4 = rho0 * Q2 * v_ext
print(f"4. F = ρ₀ * Q * v_ext:              {np.linalg.norm(F4):.6e}  (ratio: {F_expected/np.linalg.norm(F4):.2f})")

# Formula 5: ρ₀² * (4/3) * (Q/4π) * v_ext
F5 = rho0**2 * (4.0/3.0) * (Q2/(4.0*np.pi)) * v_ext
print(f"5. F = ρ₀² * (4/3) * (Q/4π) * v_ext: {np.linalg.norm(F5):.6e}  (ratio: {F_expected/np.linalg.norm(F5):.2f})")

# Formula 6: Direct from paper's formula with Q/(4πρ₀) definition
# If v_ext = Σ Q/(4πρ₀) r/r³, and F = ρ₀ Q v_ext from paper
F6 = rho0 * Q2 * v_ext
print(f"6. F = ρ₀ * Q * v_ext (paper):      {np.linalg.norm(F6):.6e}  (ratio: {F_expected/np.linalg.norm(F6):.2f})")

print(f"\nTarget ratio should be 1.00")
