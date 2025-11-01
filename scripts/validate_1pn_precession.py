#!/usr/bin/env python3
"""
1PN Precession Validation
==========================

This script validates that the superfluid model reproduces General Relativity's
perihelion precession at the 1st post-Newtonian (1PN) order, WITHOUT free parameters.

Method:
-------
1. Run Mercury orbit simulations at multiple eccentricities: e = 0.1, 0.2, 0.3, 0.5, 0.7
2. For each eccentricity, measure perihelion precession Δω per orbit
3. Compute GR-1PN theoretical prediction:
       Δω_GR = (6π G M_sun) / (a c² (1 - e²))
4. Plot Δω vs a(1-e²)⁻¹ with GR overlay
5. Report agreement: |Δω_slab / Δω_GR - 1|

This is a CRITICAL test: GR precession depends on (G, M, a, e, c) with NO tunable
parameters. If the slab matches GR, it proves the 1PN effects emerge from the
superfluid dynamics with the sound speed cs playing the role of c.

Expected results:
----------------
- Precession scales as 1/(1-e²) for fixed semi-major axis
- Agreement with GR: |Δω_slab / Δω_GR - 1| < 5-10%
- cs parameter maps to speed of light c

Usage:
------
    python scripts/validate_1pn_precession.py
    python scripts/validate_1pn_precession.py --eccentricities 0.1 0.2 0.3 --steps 50000
    python scripts/validate_1pn_precession.py --output output/1pn_precession.png
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
import argparse
import copy

from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit, estimate_orbital_period
from slab.gr1pn import compute_orbit_elements


def create_mercury_config(
    eccentricity: float,
    rho0: float = 1.0,
    beta0: float = 1e10,
    cs: float = 1e4,
    semi_major_axis: float = 0.387,
) -> Tuple[Medium, List[Body], Dict]:
    """
    Create a Mercury-like orbit configuration with specified eccentricity.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity (0 = circular, <1 = elliptical)
    rho0 : float
        Ambient density [code units]
    beta0 : float
        Mass-intake factor [code units]
    cs : float
        Sound speed [AU/yr], controls 1PN-like corrections
    semi_major_axis : float
        Semi-major axis [AU]

    Returns
    -------
    medium : Medium
        Medium object
    bodies : List[Body]
        Two bodies: Sun and Mercury
    params : Dict
        Orbital parameters (a, e, v_perihelion, etc.)
    """
    # Masses (code units: M_sun = 1, M_mercury = 3.3e-7)
    M_sun = 1.0
    M_mercury = 3.3e-7

    # Create medium
    medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)

    # For a given semi-major axis a and eccentricity e,
    # the perihelion distance is: r_peri = a(1-e)
    # and the orbital speed at perihelion is: v_peri = sqrt(K*M*(1+e)/(a(1-e)))
    # where K is the "gravitational constant" analog

    K = medium.K
    a = semi_major_axis
    e = eccentricity

    # Perihelion distance and velocity
    r_peri = a * (1 - e)
    v_peri = np.sqrt(K * M_sun * (1 + e) / (a * (1 - e)))

    # Place Sun at origin
    sun = Body(
        name="Sun",
        M=M_sun,
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        R=0.001,  # Control surface radius
        Q=M_sun / beta0,
    )

    # Place Mercury at perihelion (along +x axis) with velocity in +y direction
    mercury = Body(
        name="Mercury",
        M=M_mercury,
        x=np.array([r_peri, 0.0, 0.0]),
        v=np.array([0.0, v_peri, 0.0]),
        R=0.0005,
        Q=M_mercury / beta0,
    )

    bodies = [sun, mercury]

    # Compute orbital period using Kepler's 3rd law: T = 2π sqrt(a³/(K*M))
    T_orbit = 2 * np.pi * np.sqrt(a**3 / (K * M_sun))

    params = {
        'a': a,
        'e': e,
        'r_peri': r_peri,
        'r_apo': a * (1 + e),
        'v_peri': v_peri,
        'T_orbit': T_orbit,
        'K': K,
    }

    return medium, bodies, params


def gr_1pn_precession(
    M_sun: float,
    a: float,
    e: float,
    G: float,
    c: float
) -> float:
    """
    Compute GR-1PN perihelion precession per orbit.

    Formula from Soffel (1989), Eq. (8.42):
        Δω = (6π G M) / (a c² (1 - e²))   [radians per orbit]

    Parameters
    ----------
    M_sun : float
        Mass of central body
    a : float
        Semi-major axis
    e : float
        Eccentricity
    G : float
        Gravitational constant (or K for slab)
    c : float
        Speed of light (or cs for slab)

    Returns
    -------
    domega : float
        Precession angle per orbit [radians]
    """
    return (6 * np.pi * G * M_sun) / (a * c**2 * (1 - e**2))


def run_eccentricity_sweep(
    eccentricities: List[float],
    n_steps: int = 100000,
    dt: float = 0.002,
    n_points: int = 256,
    use_compressible: bool = True,
    verbose: bool = True,
) -> Dict[str, List]:
    """
    Run simulations across multiple eccentricities and extract precession rates.

    Parameters
    ----------
    eccentricities : List[float]
        List of eccentricities to test
    n_steps : int
        Number of integration steps per simulation
    dt : float
        Time step [yr]
    n_points : int
        Number of quadrature points for audits
    use_compressible : bool
        Enable compressible forces (1PN-like corrections)
    verbose : bool
        Print progress

    Returns
    -------
    results : Dict
        Dictionary with lists:
        - 'eccentricities': input eccentricities
        - 'precession_slab': measured slab precession [rad/orbit]
        - 'precession_gr': theoretical GR precession [rad/orbit]
        - 'a': semi-major axes
        - 'params': orbital parameters for each case
    """
    if verbose:
        print("=" * 70)
        print("1PN PRECESSION VALIDATION: ECCENTRICITY SWEEP")
        print("=" * 70)
        print()
        print(f"Testing {len(eccentricities)} eccentricities:")
        print(f"  e = {eccentricities}")
        print(f"  Integration: {n_steps} steps × dt={dt} yr = {n_steps*dt:.1f} yr")
        print(f"  Compressible forces: {use_compressible}")
        print()

    # Medium parameters (same for all runs)
    rho0 = 1.0
    # beta0 = 0.045 gives K ≈ 39.4 AU³/(M_sun yr²), matching realistic solar system
    # dynamics where Mercury's orbital period is ~0.24 years (not billions of years).
    # This enables precession measurement over reasonable integration times.
    beta0 = 0.045

    # cs = speed of light in AU/yr for proper 1PN comparison with GR
    # c = 299792 km/s = 63239.7263 AU/yr
    # In the slab model, cs (sound speed in the superfluid) plays the role of c (speed of light)
    # For the model to reproduce GR 1PN effects, we set cs = c.
    cs = 63239.7263  # Speed of light [AU/yr], for GR 1PN validation

    # For GR comparison, we need to map K → G and cs → c
    # We'll use the fact that K = rho0/(4*pi*beta0^2)
    medium_ref = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
    K = medium_ref.K

    # For GR, use K as "effective G" and cs as "effective c"
    G_eff = K
    c_eff = cs

    # Semi-major axis (same for all runs, for fair comparison)
    a = 0.387  # AU, Mercury's value

    # Storage
    results = {
        'eccentricities': [],
        'precession_slab': [],
        'precession_gr': [],
        'a': [],
        'params': [],
        'n_orbits': [],
    }

    for i, ecc in enumerate(eccentricities):
        if verbose:
            print(f"[{i+1}/{len(eccentricities)}] Running e = {ecc:.3f}...")

        # Create configuration
        medium, bodies, params = create_mercury_config(
            eccentricity=ecc,
            rho0=rho0,
            beta0=beta0,
            cs=cs,
            semi_major_axis=a,
        )

        if verbose:
            print(f"    a = {params['a']:.4f} AU, e = {params['e']:.3f}")
            print(f"    r_peri = {params['r_peri']:.4f} AU, r_apo = {params['r_apo']:.4f} AU")
            print(f"    v_peri = {params['v_peri']:.6e} AU/yr")
            print(f"    T_orbit = {params['T_orbit']:.4f} yr")
            print(f"    Expected n_orbits = {n_steps * dt / params['T_orbit']:.1f}")

        # Integration options
        opts = {
            'use_compressible': use_compressible,
            'use_flux_mass': False,
            'flux_every': 1000,
            'save_every': 100,
            'audit_every': 0,  # Disable audits for speed
            'audit_tolerance': 1e-3,
            'n_points': n_points,
            'verbose': False,
            'progress_every': n_steps // 10,
        }

        # Run simulation
        bodies_copy = [copy.deepcopy(b) for b in bodies]
        trajectory, diagnostics = integrate_orbit(
            bodies_copy, medium, dt, n_steps, opts
        )

        # Compute precession from trajectory by tracking argument of periapsis
        # Extract Mercury's trajectory (index 1)
        t_traj = trajectory['t']
        x_merc = trajectory['x'][:, 1, :]  # Mercury positions
        v_merc = trajectory['v'][:, 1, :]  # Mercury velocities

        # Compute osculating omega (argument of periapsis) at each time
        omega_list = []
        for i in range(len(t_traj)):
            elems = compute_orbit_elements(x_merc[i], v_merc[i], bodies[0].M, K)
            omega_list.append(elems['omega'])

        omega = np.array(omega_list)

        # Unwrap 2π discontinuities
        omega_unwrapped = np.unwrap(omega)

        # Fit linear trend: omega(t) = omega_0 + (dω/dt) * t
        if len(t_traj) > 10:
            coeffs = np.polyfit(t_traj, omega_unwrapped, deg=1)
            domega_dt = coeffs[0]  # radians per unit time

            # Compute precession per orbit
            T_orbit = params['T_orbit']
            domega_per_orbit_slab = domega_dt * T_orbit
            n_orbits = (t_traj[-1] - t_traj[0]) / T_orbit
        else:
            domega_per_orbit_slab = np.nan
            n_orbits = 0

        # Compute GR theoretical prediction
        M_sun = bodies[0].M
        domega_per_orbit_gr = gr_1pn_precession(M_sun, a, ecc, G_eff, c_eff)

        # Store results
        results['eccentricities'].append(ecc)
        results['precession_slab'].append(domega_per_orbit_slab)
        results['precession_gr'].append(domega_per_orbit_gr)
        results['a'].append(a)
        results['params'].append(params)
        results['n_orbits'].append(n_orbits)

        if verbose:
            print(f"    Precession (slab): {domega_per_orbit_slab:.6e} rad/orbit")
            print(f"    Precession (GR):   {domega_per_orbit_gr:.6e} rad/orbit")
            if np.isfinite(domega_per_orbit_slab) and domega_per_orbit_gr != 0:
                ratio = domega_per_orbit_slab / domega_per_orbit_gr
                print(f"    Ratio (slab/GR):   {ratio:.4f}")
                print(f"    Error:             {abs(ratio - 1.0) * 100:.2f}%")
            print()

    return results


def plot_precession_validation(
    results: Dict,
    output_path: str = "output/1pn_precession_validation.png"
):
    """
    Create validation plot: Δω vs a(1-e²)⁻¹ with GR overlay.

    Parameters
    ----------
    results : Dict
        Results from run_eccentricity_sweep
    output_path : str
        Path to save plot
    """
    eccentricities = np.array(results['eccentricities'])
    precession_slab = np.array(results['precession_slab'])
    precession_gr = np.array(results['precession_gr'])
    a_values = np.array(results['a'])

    # Compute x-axis: a/(1-e²)
    x_values = a_values / (1 - eccentricities**2)

    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Δω vs a(1-e²)
    ax = axes[0]
    ax.plot(
        x_values,
        precession_slab * 206265,  # Convert to arcsec
        'o',
        markersize=10,
        label='Slab simulation',
        color='blue',
        alpha=0.7,
    )
    ax.plot(
        x_values,
        precession_gr * 206265,
        's',
        markersize=8,
        label='GR 1PN (theory)',
        color='red',
        alpha=0.5,
    )

    # Add GR prediction line
    x_line = np.linspace(x_values.min() * 0.9, x_values.max() * 1.1, 100)
    # Δω = (6π G M) / (c² a(1-e²))
    # So Δω ∝ 1/(a(1-e²)) = 1/x
    # Get slope from first GR point
    if len(precession_gr) > 0 and x_values[0] != 0:
        slope = precession_gr[0] * x_values[0]  # Δω * x = constant
        y_line = slope / x_line
        ax.plot(
            x_line,
            y_line * 206265,
            '--',
            color='red',
            linewidth=2,
            alpha=0.3,
            label='GR 1PN (1/x trend)',
        )

    ax.set_xlabel(r'$a/(1-e^2)$ [AU]', fontsize=12)
    ax.set_ylabel(r'Perihelion precession $\Delta\omega$ [arcsec/orbit]', fontsize=12)
    ax.set_title('1PN Precession vs Eccentricity', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add eccentricity labels
    for i, (x, y, e) in enumerate(zip(x_values, precession_slab * 206265, eccentricities)):
        ax.annotate(
            f'e={e:.2f}',
            xy=(x, y),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8,
            alpha=0.6,
        )

    # Right panel: Agreement ratio
    ax = axes[1]
    ratios = precession_slab / precession_gr
    rel_errors = np.abs(ratios - 1.0) * 100  # Percent

    ax.plot(
        eccentricities,
        ratios,
        'o-',
        markersize=8,
        linewidth=2,
        color='green',
        label='Slab / GR',
    )
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Perfect agreement')
    ax.fill_between(
        [eccentricities.min(), eccentricities.max()],
        0.95, 1.05,
        color='green',
        alpha=0.1,
        label='±5% band',
    )

    ax.set_xlabel('Eccentricity e', fontsize=12)
    ax.set_ylabel('Precession ratio (Slab / GR)', fontsize=12)
    ax.set_title('Agreement with GR 1PN', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0.8, 1.2])

    # Add error annotations
    for i, (e, ratio, err) in enumerate(zip(eccentricities, ratios, rel_errors)):
        ax.annotate(
            f'{err:.1f}%',
            xy=(e, ratio),
            xytext=(5, -10 if ratio > 1.0 else 10),
            textcoords='offset points',
            fontsize=8,
            alpha=0.6,
        )

    plt.tight_layout()

    # Save plot
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")


def main():
    """Run 1PN precession validation."""
    parser = argparse.ArgumentParser(
        description="Validate 1PN precession emergence from superfluid dynamics"
    )
    parser.add_argument(
        "--eccentricities",
        type=float,
        nargs='+',
        default=[0.1, 0.2, 0.3, 0.5, 0.7],
        help="Eccentricities to test (default: 0.1 0.2 0.3 0.5 0.7)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=100000,
        help="Integration steps per simulation (default: 100000)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.002,
        help="Time step in years (default: 0.002)",
    )
    parser.add_argument(
        "--no-compressible",
        action="store_true",
        help="Disable compressible forces (should give zero precession)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/1pn_precession_validation.png",
        help="Output plot path",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("1PN PERIHELION PRECESSION VALIDATION")
    print("=" * 70)
    print()
    print("This test validates that perihelion precession emerges from")
    print("superfluid compressible corrections, matching GR 1PN predictions")
    print("WITHOUT free parameters.")
    print()

    # Run sweep
    results = run_eccentricity_sweep(
        eccentricities=args.eccentricities,
        n_steps=args.steps,
        dt=args.dt,
        use_compressible=not args.no_compressible,
        verbose=True,
    )

    # Print summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"{'Eccentricity':<15} {'Δω_slab [arcsec]':<20} {'Δω_GR [arcsec]':<20} {'Ratio':<10} {'Error':<10}")
    print("-" * 75)

    eccentricities = results['eccentricities']
    precession_slab = results['precession_slab']
    precession_gr = results['precession_gr']

    for e, p_slab, p_gr in zip(eccentricities, precession_slab, precession_gr):
        p_slab_arcsec = p_slab * 206265
        p_gr_arcsec = p_gr * 206265
        ratio = p_slab / p_gr if p_gr != 0 else np.nan
        error = abs(ratio - 1.0) * 100 if np.isfinite(ratio) else np.nan

        print(
            f"{e:<15.3f} {p_slab_arcsec:<20.4e} {p_gr_arcsec:<20.4e} "
            f"{ratio:<10.4f} {error:<10.2f}%"
        )

    print()

    # Assess results
    ratios = np.array(precession_slab) / np.array(precession_gr)
    valid_ratios = ratios[np.isfinite(ratios)]

    if len(valid_ratios) > 0:
        mean_ratio = np.mean(valid_ratios)
        std_ratio = np.std(valid_ratios)
        max_error = np.max(np.abs(valid_ratios - 1.0)) * 100

        print(f"Mean ratio (slab/GR):  {mean_ratio:.4f} ± {std_ratio:.4f}")
        print(f"Max error:             {max_error:.2f}%")
        print()

        if max_error < 5:
            print("✓ EXCELLENT: Agreement within 5%")
        elif max_error < 10:
            print("✓ GOOD: Agreement within 10%")
        elif max_error < 20:
            print("⚠ FAIR: Agreement within 20% - investigate discrepancies")
        else:
            print("✗ POOR: Agreement exceeds 20% - major issues")

    print()

    # Create plot
    plot_precession_validation(results, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
