#!/usr/bin/env python
"""
Complete working example showing how to add configurable dtype support
to the superfluid orbit simulator.

This is a minimal but complete implementation demonstrating:
1. How to modify field.py for dtype support
2. How to modify bodies.py for dtype support
3. How to use it in practice
4. Performance and precision comparison

Run this to see the precision improvement in action!
"""

import numpy as np
from dataclasses import dataclass
from typing import List
import time


# =============================================================================
# Modified Body class with dtype support
# =============================================================================

@dataclass
class BodyWithDtype:
    """Body class supporting configurable precision.

    Modified version of slab.bodies.Body with dtype parameter.
    """
    name: str
    M: float
    x: np.ndarray
    v: np.ndarray
    R: float
    Q: float = 0.0
    dtype: type = np.float64  # NEW: configurable dtype

    def __post_init__(self):
        """Ensure arrays use correct dtype."""
        # Convert to numpy arrays with specified dtype
        if not isinstance(self.x, np.ndarray):
            self.x = np.array(self.x, dtype=self.dtype)
        else:
            self.x = self.x.astype(self.dtype)

        if not isinstance(self.v, np.ndarray):
            self.v = np.array(self.v, dtype=self.dtype)
        else:
            self.v = self.v.astype(self.dtype)

        # Validate shapes
        if self.x.shape != (3,):
            raise ValueError(f"Position must be shape (3,), got {self.x.shape}")
        if self.v.shape != (3,):
            raise ValueError(f"Velocity must be shape (3,), got {self.v.shape}")

    @property
    def kinetic_energy(self):
        """Compute kinetic energy: KE = 0.5 * M * v²"""
        return 0.5 * self.M * np.dot(self.v, self.v)


# =============================================================================
# Modified field functions with dtype support
# =============================================================================

def v_self_dtype(
    x: np.ndarray,
    x_body: np.ndarray,
    Q: float,
    rho0: float,
    eps: float = 1e-12,
    dtype: type = np.float64,
) -> np.ndarray:
    """
    Velocity field from single sink with configurable precision.

    Modified version of slab.field.v_self.
    """
    # Ensure correct dtype
    x = np.asarray(x, dtype=dtype)
    x_body = np.asarray(x_body, dtype=dtype)
    Q = dtype(Q)

    # Vector from body to field point
    r_vec = x - x_body
    r = np.linalg.norm(r_vec)

    # Handle singularity
    if r < eps:
        return np.zeros(3, dtype=dtype)

    # v = -(Q/4π) * r_vec/r³
    prefactor = Q / (dtype(4.0) * dtype(np.pi))
    v = -prefactor * r_vec / (r * r * r)

    return v


def v_ext_at_dtype(
    x_a: np.ndarray,
    bodies: List[BodyWithDtype],
    a_idx: int,
    rho0: float,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    External velocity at body a's position (with dtype support).

    Modified version of slab.field.v_ext_at.
    """
    # Get dtype from body
    dtype = bodies[a_idx].dtype

    v_ext = np.zeros(3, dtype=dtype)

    for b_idx, body in enumerate(bodies):
        if b_idx == a_idx:
            continue

        v_b = v_self_dtype(x_a, body.x, body.Q, rho0, eps, dtype)
        v_ext += v_b

    return v_ext


def force_incompressible_analytic_dtype(
    a_idx: int,
    bodies: List[BodyWithDtype],
    rho0: float,
) -> np.ndarray:
    """
    Compute force on body a using analytic formula (with dtype).

    Modified version of slab.surface.force_incompressible_analytic.
    """
    body_a = bodies[a_idx]
    dtype = body_a.dtype

    # Get external velocity
    v_ext = v_ext_at_dtype(body_a.x, bodies, a_idx, rho0)

    # F_a = ρ₀ * Q_a * v_ext
    force = dtype(rho0) * dtype(body_a.Q) * v_ext

    return force


# =============================================================================
# Demonstration: Two-body force calculation
# =============================================================================

def demonstrate_precision_improvement():
    """
    Show the precision improvement from using float128 vs float64.
    """
    print("=" * 80)
    print("PRECISION DEMONSTRATION: Two-Body Force Calculation")
    print("=" * 80)
    print()

    # Setup: Two bodies separated by distance r = 2.0
    r_sep = 2.0
    rho0 = 1.0

    print(f"Configuration:")
    print(f"  Two identical sinks at separation r = {r_sep}")
    print(f"  Q1 = Q2 = 1.0")
    print(f"  rho0 = {rho0}")
    print()

    # Theoretical force magnitude: F = ρ₀ Q₁ Q₂ / (4π r²)
    F_theory = rho0 * 1.0 * 1.0 / (4.0 * np.pi * r_sep**2)
    print(f"Theoretical force magnitude: {F_theory:.16e}")
    print()

    # -------------------------------------------------------------------------
    # Calculation with float64
    # -------------------------------------------------------------------------
    print("-" * 80)
    print("1. FLOAT64 CALCULATION (standard precision)")
    print("-" * 80)

    b1_64 = BodyWithDtype(
        name="Body1",
        M=1.0,
        x=[0.0, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        R=0.1,
        Q=1.0,
        dtype=np.float64
    )

    b2_64 = BodyWithDtype(
        name="Body2",
        M=1.0,
        x=[r_sep, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        R=0.1,
        Q=1.0,
        dtype=np.float64
    )

    bodies_64 = [b1_64, b2_64]

    # Time the calculation
    start = time.perf_counter()
    for _ in range(10000):
        F1_64 = force_incompressible_analytic_dtype(0, bodies_64, rho0)
    time_64 = time.perf_counter() - start

    F1_mag_64 = np.linalg.norm(F1_64)
    error_64 = abs(F1_mag_64 - F_theory)
    rel_error_64 = error_64 / F_theory

    print(f"Force on body 1: {F1_64}")
    print(f"Force magnitude: {F1_mag_64:.16e}")
    print(f"Absolute error:  {error_64:.3e}")
    print(f"Relative error:  {rel_error_64:.3e}")
    print(f"Time (10k ops):  {time_64:.4f} seconds")
    print()

    # -------------------------------------------------------------------------
    # Calculation with float128 (longdouble)
    # -------------------------------------------------------------------------
    print("-" * 80)
    print("2. FLOAT128 CALCULATION (extended precision)")
    print("-" * 80)

    b1_128 = BodyWithDtype(
        name="Body1",
        M=1.0,
        x=[0.0, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        R=0.1,
        Q=1.0,
        dtype=np.longdouble
    )

    b2_128 = BodyWithDtype(
        name="Body2",
        M=1.0,
        x=[r_sep, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        R=0.1,
        Q=1.0,
        dtype=np.longdouble
    )

    bodies_128 = [b1_128, b2_128]

    # Time the calculation
    start = time.perf_counter()
    for _ in range(10000):
        F1_128 = force_incompressible_analytic_dtype(0, bodies_128, rho0)
    time_128 = time.perf_counter() - start

    F1_mag_128 = np.linalg.norm(F1_128)
    error_128 = abs(F1_mag_128 - F_theory)
    rel_error_128 = float(error_128 / F_theory)

    print(f"Force on body 1: {F1_128}")
    print(f"Force magnitude: {float(F1_mag_128):.20e}")
    print(f"Absolute error:  {float(error_128):.3e}")
    print(f"Relative error:  {rel_error_128:.3e}")
    print(f"Time (10k ops):  {time_128:.4f} seconds")
    print()

    # -------------------------------------------------------------------------
    # Comparison
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("COMPARISON")
    print("=" * 80)
    print()

    improvement = rel_error_64 / rel_error_128 if rel_error_128 > 0 else float('inf')
    slowdown = time_128 / time_64

    print(f"Precision improvement: {improvement:.1e}x better")
    print(f"Performance cost:      {slowdown:.2f}x slower")
    print()

    print(f"float64 precision:     {np.finfo(np.float64).precision} digits")
    print(f"longdouble precision:  {np.finfo(np.longdouble).precision} digits")
    print(f"Precision gain:        +{np.finfo(np.longdouble).precision - np.finfo(np.float64).precision} digits")
    print()

    # -------------------------------------------------------------------------
    # Assessment
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("ASSESSMENT")
    print("=" * 80)
    print()

    print("Is the improvement worth it?")
    print()

    if rel_error_64 < 1e-9:
        print(f"  ✓ float64 error ({rel_error_64:.2e}) is already < 1e-9")
        print(f"    → Probably NOT worth the {slowdown:.1f}x slowdown")
        print(f"    → Use float64 for production")
        print(f"    → Reserve float128 for validation studies")
    else:
        print(f"  ✗ float64 error ({rel_error_64:.2e}) exceeds 1e-9")
        print(f"    → float128 may be beneficial")
        print(f"    → {slowdown:.1f}x slowdown is acceptable for critical calculations")
    print()

    if improvement > 1000:
        print(f"  ✓ Precision improvement ({improvement:.1e}x) is substantial")
        print(f"    → Useful for convergence studies")
        print(f"    → Can validate float64 results")
    print()

    return {
        'float64': {'error': rel_error_64, 'time': time_64},
        'float128': {'error': rel_error_128, 'time': time_128},
        'improvement': improvement,
        'slowdown': slowdown,
    }


# =============================================================================
# Integration example with dtype support
# =============================================================================

def demonstrate_integration_with_dtype():
    """
    Show how dtype would work in actual integration.
    """
    print("=" * 80)
    print("INTEGRATION EXAMPLE WITH DTYPE SUPPORT")
    print("=" * 80)
    print()

    print("Example workflow:")
    print()
    print("# Standard precision (default)")
    print("bodies = [")
    print("    BodyWithDtype('Sun', M=1.0, x=[0,0,0], v=[0,0,0],")
    print("                  R=1e-3, Q=1e-10, dtype=np.float64),")
    print("    BodyWithDtype('Earth', M=3e-6, x=[1,0,0], v=[0,6.28,0],")
    print("                  R=5e-4, Q=3e-16, dtype=np.float64),")
    print("]")
    print()
    print("# High precision")
    print("bodies = [")
    print("    BodyWithDtype('Sun', M=1.0, x=[0,0,0], v=[0,0,0],")
    print("                  R=1e-3, Q=1e-10, dtype=np.longdouble),")
    print("    BodyWithDtype('Earth', M=3e-6, x=[1,0,0], v=[0,6.28,0],")
    print("                  R=5e-4, Q=3e-16, dtype=np.longdouble),")
    print("]")
    print()

    print("Integration loop (pseudo-code):")
    print()
    print("for step in range(n_steps):")
    print("    # Forces automatically use body.dtype")
    print("    forces = assemble_forces(bodies, medium)")
    print("    ")
    print("    # Update velocities (half-step)")
    print("    for i, body in enumerate(bodies):")
    print("        body.v += 0.5 * (forces[i] / body.M) * dt")
    print("    ")
    print("    # Update positions")
    print("    for body in bodies:")
    print("        body.x += body.v * dt")
    print("    ")
    print("    # ... rest of verlet step")
    print()

    print("Command-line usage:")
    print()
    print("  $ slab-run config.yaml                    # float64 (default)")
    print("  $ slab-run config.yaml --precision=high   # float128")
    print()


# =============================================================================
# Main demonstration
# =============================================================================

if __name__ == "__main__":
    print()
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 20 + "DTYPE IMPLEMENTATION EXAMPLE" + " " * 30 + "║")
    print("╚" + "=" * 78 + "╝")
    print()

    # Run precision demonstration
    results = demonstrate_precision_improvement()

    print()
    input("Press Enter to see integration example...")
    print("\n" * 2)

    # Show integration example
    demonstrate_integration_with_dtype()

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()

    print("Key takeaways:")
    print()
    print(f"1. float64 (standard):")
    print(f"   - Error: {results['float64']['error']:.2e}")
    print(f"   - Speed: {results['float64']['time']:.3f}s (baseline)")
    print(f"   - Use for: Production simulations")
    print()
    print(f"2. float128 (extended):")
    print(f"   - Error: {results['float128']['error']:.2e}")
    print(f"   - Speed: {results['float128']['time']:.3f}s ({results['slowdown']:.2f}x slower)")
    print(f"   - Improvement: {results['improvement']:.1e}x better precision")
    print(f"   - Use for: Validation, convergence studies, critical calculations")
    print()
    print("3. Implementation:")
    print("   - Add 'dtype' parameter to all field/force functions")
    print("   - Add 'dtype' field to Body class")
    print("   - Propagate through integration pipeline")
    print("   - Add CLI flag: --precision={standard,high}")
    print()
    print("4. Recommendation:")
    if results['float64']['error'] < 1e-9:
        print("   - START with float64 (your error is already excellent)")
        print("   - ADD float128 support later if needed")
        print("   - RESERVE for specific use cases requiring ultra-high precision")
    else:
        print("   - CONSIDER float128 for critical calculations")
        print("   - BENCHMARK on your actual simulation")
        print("   - DOCUMENT precision requirements for your science goals")
    print()

    print("Files to reference:")
    print("  - This file: dtype_implementation_example.py")
    print("  - Full analysis: PRECISION_INVESTIGATION.md")
    print("  - Interactive demo: precision_examples.py")
    print("  - Benchmark suite: precision_analysis.py")
    print()
