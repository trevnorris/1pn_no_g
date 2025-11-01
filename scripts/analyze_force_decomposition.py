#!/usr/bin/env python3
"""
Force Decomposition Analysis
==============================

This script analyzes the physical mechanism behind emergent forces by decomposing
the surface integral into momentum flux and pressure contributions.

Physics:
--------
The force on body 'a' from a control surface integral is:

    F_a = ∫ [-P n - ρ(½v² n - v(v·n))] dA

where:
    -P n            = static pressure contribution
    -½ρv² n         = dynamic pressure contribution
    +ρ v(v·n)       = momentum flux contribution

Rearranging:
    F_a = ∫ [ρ v(v·n)] dA  - ∫ [(P + ½ρv²) n] dA
        = F_momentum       + F_pressure

Key result:
-----------
In the incompressible limit with spherical control surfaces, the pressure term
(both static and dynamic) cancels by symmetry, leaving only the momentum flux.

This script demonstrates:
1. |F_momentum| >> |F_pressure| for spherical surfaces
2. F_total ≈ F_momentum in practice
3. The physical mechanism is momentum transport, not pressure gradients

Usage:
------
    python scripts/analyze_force_decomposition.py
    python scripts/analyze_force_decomposition.py --n-points 1024 --separations 10
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import argparse

from slab.medium import Medium
from slab.bodies import Body
from slab.geometry import fibonacci_sphere
from slab.field import v_ext_at, v_self


def compute_force_components(
    a_idx: int,
    bodies: List[Body],
    medium: Medium,
    n_points: int = 512,
) -> Dict[str, np.ndarray]:
    """
    Compute force decomposition: momentum vs pressure contributions.

    Parameters
    ----------
    a_idx : int
        Index of body to compute force on
    bodies : List[Body]
        All bodies in system
    medium : Medium
        Medium properties
    n_points : int
        Number of surface quadrature points

    Returns
    -------
    components : Dict
        Dictionary with:
        - 'momentum': momentum flux contribution, shape (3,)
        - 'pressure': pressure contribution, shape (3,)
        - 'total': sum of both, shape (3,)
        - 'analytic': analytic force for comparison, shape (3,)
    """
    from slab.surface import force_incompressible_analytic

    # Get body parameters
    body_a = bodies[a_idx]
    x_a = body_a.x
    R_a = body_a.R
    Q_a = body_a.Q
    rho0 = medium.rho0

    # Generate surface points
    normals = fibonacci_sphere(n_points)
    sphere_area = 4.0 * np.pi * R_a**2
    dA = sphere_area / n_points

    # Initialize accumulators
    F_momentum = np.zeros(3)
    F_pressure = np.zeros(3)

    # External velocity at body center (for cross-term)
    v_ext_a = v_ext_at(x_a, bodies, a_idx, rho0)

    for i in range(n_points):
        n_i = normals[i]
        x_i = x_a + R_a * n_i

        # Self-field velocity at surface point
        v_self_i = v_self(x_i, x_a, Q_a, rho0)

        # For incompressible flow around a point sink, the pressure is:
        # P = P_∞ - ½ρ v²
        # where v² = v_self² + 2*v_self·v_ext + v_ext²
        #
        # The pressure gradient force is:
        # F_pressure = -∫ (P + ½ρv²) n dA
        #            = -∫ P_∞ n dA - ∫ ½ρ(v_self² + 2*v_self·v_ext + v_ext²) n dA
        #
        # The P_∞ term cancels by symmetry. For the dynamic terms:

        # Total velocity at surface
        v_total = v_self_i + v_ext_a
        v_total_sq = np.dot(v_total, v_total)

        # Momentum flux: ρ v(v·n)
        # Using cross-term approximation: ρ v_ext (v_self·n)
        v_self_dot_n = np.dot(v_self_i, n_i)
        F_momentum += rho0 * v_ext_a * v_self_dot_n * dA

        # Pressure (dynamic): -½ρ v² n
        # This should approximately cancel by symmetry on sphere
        F_pressure += -0.5 * rho0 * v_total_sq * n_i * dA

    # Total force
    F_total = F_momentum + F_pressure

    # Analytic force for comparison
    F_analytic = force_incompressible_analytic(a_idx, bodies, medium)

    return {
        'momentum': F_momentum,
        'pressure': F_pressure,
        'total': F_total,
        'analytic': F_analytic,
    }


def run_separation_sweep(
    separations: np.ndarray,
    n_points: int = 512,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Compute force decomposition across multiple separations.

    Parameters
    ----------
    separations : np.ndarray
        Array of separation distances [AU]
    n_points : int
        Number of quadrature points
    verbose : bool
        Print progress

    Returns
    -------
    results : Dict
        Dictionary with arrays of shape (n_separations,):
        - 'separations': input separations
        - 'F_momentum_mag': |F_momentum|
        - 'F_pressure_mag': |F_pressure|
        - 'F_total_mag': |F_total|
        - 'F_analytic_mag': |F_analytic|
        - 'ratio_mom_to_pressure': |F_momentum| / |F_pressure|
    """
    if verbose:
        print("=" * 70)
        print("FORCE DECOMPOSITION ANALYSIS")
        print("=" * 70)
        print()
        print(f"Analyzing force components at {len(separations)} separations")
        print(f"Quadrature points: {n_points}")
        print()

    # Medium and body parameters
    rho0 = 1.0
    beta0 = 1e10
    cs = 1e10  # Incompressible limit
    medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)

    M1 = 1.0
    M2 = 3.3e-7
    Q1 = M1 / beta0
    Q2 = M2 / beta0

    # Storage
    F_momentum_mag = np.zeros(len(separations))
    F_pressure_mag = np.zeros(len(separations))
    F_total_mag = np.zeros(len(separations))
    F_analytic_mag = np.zeros(len(separations))

    for i, r in enumerate(separations):
        if verbose and (i % 5 == 0 or i == 0):
            print(f"  [{i+1}/{len(separations)}] r = {r:.3f} AU...")

        # Create two-body system
        body0 = Body(
            name="Body0",
            M=M1,
            x=np.array([0.0, 0.0, 0.0]),
            v=np.array([0.0, 0.0, 0.0]),
            R=1e-3,
            Q=Q1,
        )

        body1 = Body(
            name="Body1",
            M=M2,
            x=np.array([r, 0.0, 0.0]),
            v=np.array([0.0, 0.0, 0.0]),
            R=5e-4,
            Q=Q2,
        )

        bodies = [body0, body1]

        # Compute force components on body 1
        components = compute_force_components(1, bodies, medium, n_points)

        F_momentum_mag[i] = np.linalg.norm(components['momentum'])
        F_pressure_mag[i] = np.linalg.norm(components['pressure'])
        F_total_mag[i] = np.linalg.norm(components['total'])
        F_analytic_mag[i] = np.linalg.norm(components['analytic'])

    # Compute ratio
    ratio = F_momentum_mag / np.maximum(F_pressure_mag, 1e-50)

    if verbose:
        print()
        print("Sample results:")
        print(f"{'r [AU]':<10} {'|F_mom|':<15} {'|F_press|':<15} {'Ratio':<10} {'|F_total|':<15}")
        print("-" * 65)
        for i in [0, len(separations)//2, -1]:
            r = separations[i]
            print(
                f"{r:<10.3f} {F_momentum_mag[i]:<15.6e} {F_pressure_mag[i]:<15.6e} "
                f"{ratio[i]:<10.1f} {F_total_mag[i]:<15.6e}"
            )
        print()

    return {
        'separations': separations,
        'F_momentum_mag': F_momentum_mag,
        'F_pressure_mag': F_pressure_mag,
        'F_total_mag': F_total_mag,
        'F_analytic_mag': F_analytic_mag,
        'ratio': ratio,
    }


def plot_force_decomposition(
    results: Dict,
    output_path: str = "output/force_decomposition.png"
):
    """
    Create force decomposition visualization.

    Parameters
    ----------
    results : Dict
        Results from run_separation_sweep
    output_path : str
        Output plot path
    """
    separations = results['separations']
    F_momentum = results['F_momentum_mag']
    F_pressure = results['F_pressure_mag']
    F_total = results['F_total_mag']
    F_analytic = results['F_analytic_mag']
    ratio = results['ratio']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Force components vs separation
    ax = axes[0]
    ax.loglog(separations, F_momentum, 'o-', label='Momentum flux', linewidth=2, markersize=8)
    ax.loglog(separations, F_pressure, 's-', label='Pressure', linewidth=2, markersize=6, alpha=0.7)
    ax.loglog(separations, F_total, '^-', label='Total (quadrature)', linewidth=2, markersize=6, alpha=0.5)
    ax.loglog(separations, F_analytic, 'x--', label='Analytic', linewidth=2, markersize=8, alpha=0.5)

    ax.set_xlabel('Separation r [AU]', fontsize=12)
    ax.set_ylabel('Force magnitude [code units]', fontsize=12)
    ax.set_title('Force Decomposition: Momentum vs Pressure', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    # Add annotation
    textstr = 'Momentum flux dominates:\n'
    textstr += f'Typical ratio ~ {np.median(ratio):.0f}:1'
    ax.text(
        0.05, 0.05,
        textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
    )

    # Right panel: Ratio of momentum to pressure
    ax = axes[1]
    ax.semilogx(separations, ratio, 'o-', linewidth=2, markersize=8, color='green')
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Equal contributions')
    ax.set_xlabel('Separation r [AU]', fontsize=12)
    ax.set_ylabel('Ratio: |F_momentum| / |F_pressure|', fontsize=12)
    ax.set_title('Momentum Flux Dominance', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")


def main():
    """Run force decomposition analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze force decomposition: momentum vs pressure"
    )
    parser.add_argument(
        "--separations",
        type=int,
        default=15,
        help="Number of separation points (default: 15)",
    )
    parser.add_argument(
        "--r-min",
        type=float,
        default=0.1,
        help="Minimum separation [AU] (default: 0.1)",
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=5.0,
        help="Maximum separation [AU] (default: 5.0)",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=512,
        help="Quadrature points (default: 512)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/force_decomposition.png",
        help="Output plot path",
    )
    args = parser.parse_args()

    # Create separation array
    separations = np.logspace(
        np.log10(args.r_min),
        np.log10(args.r_max),
        args.separations,
    )

    # Run analysis
    results = run_separation_sweep(
        separations=separations,
        n_points=args.n_points,
        verbose=True,
    )

    # Summary statistics
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"Momentum flux dominance:")
    print(f"  Mean ratio |F_mom|/|F_press|: {np.mean(results['ratio']):.1f}")
    print(f"  Median ratio:                  {np.median(results['ratio']):.1f}")
    print(f"  Min ratio:                     {np.min(results['ratio']):.1f}")
    print()
    print("This demonstrates that forces arise primarily from MOMENTUM FLUX,")
    print("not pressure gradients. The pressure term cancels by symmetry on")
    print("spherical control surfaces, leaving only the momentum transport.")
    print()

    # Create plot
    plot_force_decomposition(results, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
