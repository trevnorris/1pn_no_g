#!/usr/bin/env python3
"""
Emergent Inverse-Square Force Law Validation
=============================================

This script demonstrates that F ∝ 1/r² emerges naturally from surface
integrals of momentum flux, without hard-coding Newton's law.

Method:
-------
1. Create two static bodies with fixed Q₁, Q₂, ρ₀
2. Compute force via surface integral for separations r ∈ [0.1, 10] AU
3. Fit F(r) = C/r² to extract numerical coefficient C_num
4. Compare to theoretical: C_theory = ρ₀|Q₁Q₂|/(4π)
5. Report relative error: |C_num/C_theory - 1|

This proves the inverse-square law is emergent, not assumed.

Expected results:
----------------
- Log-log plot of F vs r shows slope = -2.0 ± 0.001
- Coefficient agreement: |C_num/C_theory - 1| < 0.5%
- Forces from analytic and quadrature methods agree to ~0.01%

Usage:
------
    python scripts/validate_force_law.py
    python scripts/validate_force_law.py --separations 20 --n-points 1024
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
import argparse

from slab.medium import Medium
from slab.bodies import Body
from slab.surface import (
    force_incompressible_analytic,
    force_incompressible_quadrature,
)


def setup_two_body_system(
    separation: float,
    Q1: float,
    Q2: float,
    medium: Medium,
    R1: float = 1e-3,
    R2: float = 5e-4,
) -> List[Body]:
    """
    Create two static bodies separated by distance r.

    Parameters
    ----------
    separation : float
        Distance between bodies [AU]
    Q1, Q2 : float
        Volumetric intake rates [AU³/yr]
    medium : Medium
        Medium properties (rho0, cs, beta0)
    R1, R2 : float
        Control surface radii [AU]

    Returns
    -------
    bodies : List[Body]
        Two-body system: body 0 at origin, body 1 at (separation, 0, 0)
    """
    # Body 0 at origin
    M1 = medium.beta0 * Q1
    body0 = Body(
        name="Body0",
        M=M1,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=R1,
        Q=Q1,
    )

    # Body 1 along x-axis
    M2 = medium.beta0 * Q2
    body1 = Body(
        name="Body1",
        M=M2,
        x=np.array([separation, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=R2,
        Q=Q2,
    )

    return [body0, body1]


def compute_force_at_separation(
    separation: float,
    Q1: float,
    Q2: float,
    medium: Medium,
    n_points: int = 512,
    use_quadrature: bool = False,
) -> Tuple[float, np.ndarray]:
    """
    Compute force magnitude at given separation using either analytic or quadrature.

    Parameters
    ----------
    separation : float
        Distance between bodies [AU]
    Q1, Q2 : float
        Volumetric intake rates [AU³/yr]
    medium : Medium
        Medium properties
    n_points : int
        Number of quadrature points (if use_quadrature=True)
    use_quadrature : bool
        Use quadrature method instead of analytic formula

    Returns
    -------
    force_magnitude : float
        Magnitude of force on body 1 [code units]
    force_vector : ndarray
        Force vector on body 1 (3,)
    """
    bodies = setup_two_body_system(separation, Q1, Q2, medium)

    # Compute force on body 1 (index 1)
    if use_quadrature:
        F = force_incompressible_quadrature(
            a_idx=1,
            bodies=bodies,
            medium=medium,
            n_points=n_points,
        )
    else:
        F = force_incompressible_analytic(
            a_idx=1,
            bodies=bodies,
            medium=medium,
        )

    return np.linalg.norm(F), F


def fit_power_law(r_values: np.ndarray, F_values: np.ndarray) -> Tuple[float, float]:
    """
    Fit F = C/r^n to data and extract coefficient C and exponent n.

    Uses log-linear regression: log(F) = log(C) - n*log(r)

    Parameters
    ----------
    r_values : ndarray
        Separation distances
    F_values : ndarray
        Force magnitudes

    Returns
    -------
    C : float
        Fitted coefficient
    n : float
        Fitted exponent (should be ~2.0 for inverse-square)
    """
    # Take logs
    log_r = np.log(r_values)
    log_F = np.log(F_values)

    # Linear fit: log(F) = a + b*log(r), where b = -n
    coeffs = np.polyfit(log_r, log_F, deg=1)
    b, a = coeffs  # polyfit returns [slope, intercept]

    C = np.exp(a)
    n = -b

    return C, n


def theoretical_coefficient(Q1: float, Q2: float, rho0: float) -> float:
    """
    Compute theoretical force coefficient from control-surface lemma.

    From the paper, the incompressible force between two static sinks is:
        F = (ρ₀ |Q₁ Q₂|) / (4π r²)

    So the coefficient C in F = C/r² is:
        C_theory = ρ₀ |Q₁ Q₂| / (4π)

    Parameters
    ----------
    Q1, Q2 : float
        Volumetric intake rates [AU³/yr]
    rho0 : float
        Ambient density [code units]

    Returns
    -------
    C_theory : float
        Theoretical coefficient
    """
    return rho0 * abs(Q1 * Q2) / (4.0 * np.pi)


def main():
    """Run inverse-square validation scan."""
    parser = argparse.ArgumentParser(
        description="Validate emergent inverse-square force law"
    )
    parser.add_argument(
        "--separations",
        type=int,
        default=15,
        help="Number of separation points to sample",
    )
    parser.add_argument(
        "--r-min",
        type=float,
        default=0.1,
        help="Minimum separation [AU]",
    )
    parser.add_argument(
        "--r-max",
        type=float,
        default=10.0,
        help="Maximum separation [AU]",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=512,
        help="Number of quadrature points for validation",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/force_law_validation.png",
        help="Output plot path",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("EMERGENT INVERSE-SQUARE FORCE LAW VALIDATION")
    print("=" * 70)
    print()

    # Set up medium and body parameters
    rho0 = 1.0
    beta0 = 1e10
    cs = 1e10  # Very large cs for incompressible limit

    medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)

    # Choose two bodies with realistic Q values (Sun-like and Mercury-like)
    M1 = 1.0  # Sun mass
    M2 = 3.3e-7  # Mercury mass
    Q1 = M1 / beta0
    Q2 = M2 / beta0

    print(f"Medium parameters:")
    print(f"  ρ₀ = {rho0:.3e}")
    print(f"  β₀ = {beta0:.3e}")
    print(f"  cs = {cs:.3e} AU/yr (incompressible limit)")
    print(f"  K = {medium.K:.6e}")
    print()
    print(f"Body parameters:")
    print(f"  M₁ = {M1:.6e} (Sun-like)")
    print(f"  M₂ = {M2:.6e} (Mercury-like)")
    print(f"  Q₁ = {Q1:.6e} AU³/yr")
    print(f"  Q₂ = {Q2:.6e} AU³/yr")
    print()

    # Compute theoretical coefficient
    C_theory = theoretical_coefficient(Q1, Q2, rho0)
    print(f"Theoretical coefficient:")
    print(f"  C_theory = ρ₀|Q₁Q₂|/(4π) = {C_theory:.10e}")
    print()

    # Scan separations
    r_values = np.logspace(
        np.log10(args.r_min),
        np.log10(args.r_max),
        args.separations,
    )

    print(f"Computing forces at {len(r_values)} separations...")
    print(f"  Range: {r_values[0]:.3f} to {r_values[-1]:.3f} AU")
    print(f"  Quadrature points: {args.n_points}")
    print()

    # Compute forces using both methods
    F_analytic = np.zeros(len(r_values))
    F_quadrature = np.zeros(len(r_values))

    for i, r in enumerate(r_values):
        F_analytic[i], _ = compute_force_at_separation(
            r, Q1, Q2, medium, use_quadrature=False
        )
        F_quadrature[i], _ = compute_force_at_separation(
            r, Q1, Q2, medium, n_points=args.n_points, use_quadrature=True
        )

        # Progress indicator
        if (i + 1) % 5 == 0 or i == 0:
            rel_diff = abs(F_analytic[i] - F_quadrature[i]) / F_analytic[i]
            print(
                f"  r = {r:6.3f} AU: "
                f"F_analytic = {F_analytic[i]:.6e}, "
                f"F_quadrature = {F_quadrature[i]:.6e}, "
                f"rel_diff = {rel_diff:.3e}"
            )

    print()

    # Fit power law to analytic data
    C_num, n_fitted = fit_power_law(r_values, F_analytic)

    print(f"Power law fit results (analytic method):")
    print(f"  Fitted: F = C/r^n")
    print(f"  C_num = {C_num:.10e}")
    print(f"  n = {n_fitted:.6f} (expect 2.0 for inverse-square)")
    print()

    # Compare to theory
    rel_error_C = abs(C_num / C_theory - 1.0)
    exponent_error = abs(n_fitted - 2.0)

    print(f"Validation results:")
    print(f"  |C_num/C_theory - 1| = {rel_error_C:.6e} ({rel_error_C*100:.4f}%)")
    print(f"  |n - 2| = {exponent_error:.6e}")
    print()

    # Check analytic vs quadrature agreement
    rel_diffs = np.abs(F_analytic - F_quadrature) / F_analytic
    mean_rel_diff = np.mean(rel_diffs)
    max_rel_diff = np.max(rel_diffs)

    print(f"Analytic vs quadrature agreement:")
    print(f"  Mean relative difference: {mean_rel_diff:.3e}")
    print(f"  Max relative difference: {max_rel_diff:.3e}")
    print()

    # Assess results
    print("=" * 70)
    print("ASSESSMENT")
    print("=" * 70)
    success = True

    if rel_error_C < 0.005:  # 0.5%
        print("✓ PASS: Coefficient matches theory to < 0.5%")
    else:
        print(f"✗ FAIL: Coefficient error {rel_error_C*100:.2f}% exceeds 0.5%")
        success = False

    if exponent_error < 0.001:
        print(f"✓ PASS: Exponent n = {n_fitted:.6f} is inverse-square (2.0 ± 0.001)")
    else:
        print(f"✗ FAIL: Exponent n = {n_fitted:.6f} deviates from 2.0")
        success = False

    if max_rel_diff < 0.0001:  # 0.01%
        print("✓ PASS: Analytic and quadrature agree to < 0.01%")
    else:
        print(f"✗ FAIL: Methods disagree by {max_rel_diff*100:.4f}%")
        success = False

    print()
    if success:
        print("🎉 ALL CHECKS PASSED")
        print()
        print("This demonstrates that F ∝ 1/r² emerges from surface integrals")
        print("of momentum flux, without assuming or hard-coding Newton's law.")
    else:
        print("⚠️  SOME CHECKS FAILED - investigate discrepancies")

    print()

    # Create plot
    print(f"Generating plot: {args.output}")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Log-log plot with fit
    ax = axes[0]
    ax.loglog(r_values, F_analytic, "o", label="Analytic", markersize=8, alpha=0.7)
    ax.loglog(
        r_values,
        F_quadrature,
        "s",
        label=f"Quadrature (N={args.n_points})",
        markersize=6,
        alpha=0.5,
    )

    # Plot fitted power law
    r_fit = np.logspace(np.log10(r_values[0]), np.log10(r_values[-1]), 100)
    F_fit = C_num / r_fit**n_fitted
    ax.loglog(
        r_fit,
        F_fit,
        "--",
        color="red",
        linewidth=2,
        label=f"Fit: $F = C/r^{{{n_fitted:.4f}}}$",
    )

    ax.set_xlabel("Separation r [AU]", fontsize=12)
    ax.set_ylabel("Force magnitude |F| [code units]", fontsize=12)
    ax.set_title("Emergent Inverse-Square Law", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")

    # Add annotation with fit results
    textstr = f"$C_{{num}}$ = {C_num:.3e}\n"
    textstr += f"$C_{{theory}}$ = {C_theory:.3e}\n"
    textstr += f"Error = {rel_error_C*100:.3f}%\n"
    textstr += f"Exponent = {n_fitted:.5f}"
    ax.text(
        0.05,
        0.05,
        textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # Right panel: Relative differences
    ax = axes[1]
    ax.semilogx(r_values, rel_diffs * 100, "o-", markersize=6, linewidth=1.5)
    ax.axhline(0.01, color="green", linestyle="--", alpha=0.5, label="0.01% target")
    ax.set_xlabel("Separation r [AU]", fontsize=12)
    ax.set_ylabel("Relative difference [%]", fontsize=12)
    ax.set_title("Analytic vs Quadrature Agreement", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()

    # Save plot
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {output_path}")
    print()

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
