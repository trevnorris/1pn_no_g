#!/usr/bin/env python3
"""
Final comprehensive diagnosis of the precession issue.

We know:
1. The compressible force MAGNITUDE is correct (~7.5e-12, matches torque estimate)
2. The precession is ~300,000× too large (~0.15 rad/orbit vs ~5e-7 rad/orbit)
3. Only ~8 orbits are being simulated

Hypothesis: The huge precession is NOT from the compressible force at all!
It might be from:
- Initial conditions not being exactly at perihelion
- Apsidal oscillations in an eccentric orbit
- Integration errors
- Something else entirely

Let's check by comparing:
A. Precession WITH compressible forces
B. Precession WITHOUT compressible forces

If both show huge precession, the problem is NOT the compressible force!
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import copy
from slab.dynamics import integrate_orbit

from scripts.precession_helpers import (
    analyse_precession,
    create_barycentric_mercury_config,
)

# Create Mercury configuration (barycentric initial conditions)
rho0 = 1.0
beta0 = 0.045
cs = 63239.7263
a = 0.387
e = 0.1

medium, bodies, params = create_barycentric_mercury_config(
    eccentricity=e,
    rho0=rho0,
    beta0=beta0,
    cs=cs,
    semi_major_axis=a,
)

sun, mercury = bodies
T_orbit = params.T_orbit
K = params.K
M_total = params.M_total

print("=" * 70)
print("FINAL DIAGNOSIS: COMPARING WITH vs WITHOUT COMPRESSIBLE FORCES")
print("=" * 70)
print()
print(f"Orbit parameters:")
print(f"  a = {a:.4f} AU, e = {e:.3f}")
print(f"  T_orbit = {T_orbit:.4f} yr")
print(f"  cs = {cs:.4f} AU/yr")
print()

# Integration parameters matching the production runs.  These are intentionally
# coarse so we can reproduce the large apparent precession that motivated this
# investigation, then quantify the discretisation error by re-running with a
# finer timestep below.
n_steps = 1000
dt = 0.002
total_time = n_steps * dt
n_orbits = total_time / T_orbit
steps_per_orbit = T_orbit / dt

print(f"Integration: {n_steps} steps × {dt} yr = {total_time} yr ~ {n_orbits:.1f} orbits")
print(f"Steps per orbit ≈ {steps_per_orbit:.1f}")
print()

# Run WITHOUT compressible forces
print("=" * 70)
print("RUN 1: WITHOUT COMPRESSIBLE FORCES (pure incompressible)")
print("=" * 70)
print()

bodies_incomp = [copy.deepcopy(sun), copy.deepcopy(mercury)]
opts_incomp = {
    'use_compressible': False,
    'use_quadrature': False,
    'save_every': max(1, n_steps // 2000),
    'verbose': False,
}

traj_incomp, diag_incomp = integrate_orbit(bodies_incomp, medium, dt, n_steps, opts_incomp)

analysis_incomp = analyse_precession(traj_incomp, params)
domega_per_orbit_incomp = analysis_incomp['slope_per_orbit']
peri_per_orbit_incomp = analysis_incomp['peri_per_orbit']
peri_std_incomp = analysis_incomp['peri_std']
n_peri_incomp = analysis_incomp['n_peri_cycles']
t_incomp = analysis_incomp['t']
elems_incomp = analysis_incomp['elements']

print(f"Precession (incompressible, slope fit): {domega_per_orbit_incomp:.6e} rad/orbit")
print(
    f"                                      = {domega_per_orbit_incomp * 206265:.6e} arcsec/orbit"
)
if np.isfinite(peri_per_orbit_incomp):
    err_arcsec = peri_std_incomp * 206265.0
    print(
        f"Periapsis-to-periapsis: {peri_per_orbit_incomp:.6e} rad/orbit"
        f" ± {peri_std_incomp:.2e} (n = {n_peri_incomp})"
    )
    print(
        f"                       = {peri_per_orbit_incomp * 206265:.6e} arcsec/orbit"
        f" ± {err_arcsec:.2e}"
    )
print()

# Run WITH compressible forces
print("=" * 70)
print("RUN 2: WITH COMPRESSIBLE FORCES")
print("=" * 70)
print()

bodies_comp = [copy.deepcopy(sun), copy.deepcopy(mercury)]
opts_comp = {
    'use_compressible': True,
    'use_quadrature': False,
    'save_every': max(1, n_steps // 2000),
    'verbose': False,
    'n_points': 128,
}

traj_comp, diag_comp = integrate_orbit(bodies_comp, medium, dt, n_steps, opts_comp)

analysis_comp = analyse_precession(traj_comp, params)
domega_per_orbit_comp = analysis_comp['slope_per_orbit']
peri_per_orbit_comp = analysis_comp['peri_per_orbit']
peri_std_comp = analysis_comp['peri_std']
n_peri_comp = analysis_comp['n_peri_cycles']
t_comp = analysis_comp['t']
elems_comp = analysis_comp['elements']

print(f"Precession (compressible, slope fit): {domega_per_orbit_comp:.6e} rad/orbit")
print(
    f"                                   = {domega_per_orbit_comp * 206265:.6e} arcsec/orbit"
)
if np.isfinite(peri_per_orbit_comp):
    err_arcsec = peri_std_comp * 206265.0
    print(
        f"Periapsis-to-periapsis: {peri_per_orbit_comp:.6e} rad/orbit"
        f" ± {peri_std_comp:.2e} (n = {n_peri_comp})"
    )
    print(
        f"                       = {peri_per_orbit_comp * 206265:.6e} arcsec/orbit"
        f" ± {err_arcsec:.2e}"
    )
print()

# GR prediction
GR_precession = (6 * np.pi * K * M_total) / (a * cs**2 * (1 - e**2))
print(f"GR prediction: {GR_precession:.6e} rad/orbit")
print(f"             = {GR_precession * 206265:.6e} arcsec/orbit")
print()

# Difference between compressible and incompressible
delta_precession_slope = domega_per_orbit_comp - domega_per_orbit_incomp
if np.isfinite(peri_per_orbit_comp) and np.isfinite(peri_per_orbit_incomp):
    delta_precession_peri = peri_per_orbit_comp - peri_per_orbit_incomp
else:
    delta_precession_peri = np.nan
print("=" * 70)
print("ANALYSIS")
print("=" * 70)
print()

print(
    f"Precession difference (slope fit): {delta_precession_slope:.6e} rad/orbit"
)
print(
    f"                                   = {delta_precession_slope * 206265:.6e} arcsec/orbit"
)
if GR_precession != 0:
    print(
        f"                                   ≈ {delta_precession_slope / GR_precession:.2f} × GR"
    )
if np.isfinite(delta_precession_peri):
    print(
        f"Periapsis-to-periapsis difference: {delta_precession_peri:.6e} rad/orbit"
    )
    print(
        f"                                   = {delta_precession_peri * 206265:.6e} arcsec/orbit"
    )
    if GR_precession != 0:
        print(
            f"                                   ≈ {delta_precession_peri / GR_precession:.2f} × GR"
        )
print()

if np.isfinite(domega_per_orbit_incomp) and abs(domega_per_orbit_incomp) > 1e-3:
    print(f"WARNING: Incompressible case shows HUGE precession!")
    print(f"  |Δω_incomp| = {abs(domega_per_orbit_incomp):.6e} rad/orbit")
    print()
    print("This means the problem is NOT the compressible force correction.")
    print("Possible causes:")
    print("  1. Initial conditions are not exactly at perihelion")
    print("  2. Orbit is not exactly closed (energy drift)")
    print("  3. Measurement artifact (too few orbits, noise)")
    print("  4. Apsidal precession from non-Keplerian forces")
    print()
    print("The compressible correction is:")
    print(f"  Δ(Δω)_slope = {delta_precession_slope:.6e} rad/orbit")
    print(f"             = {delta_precession_slope * 206265:.6e} arcsec/orbit")
    if np.isfinite(delta_precession_peri):
        print(f"  Δ(Δω)_peri  = {delta_precession_peri:.6e} rad/orbit")
        print(f"             = {delta_precession_peri * 206265:.6e} arcsec/orbit")
    if abs(delta_precession_slope) < abs(domega_per_orbit_incomp) * 0.1:
        print("  This is SMALL compared to the incompressible precession.")
        print("  So the compressible force is NOT the main issue!")
else:
    print(f"Incompressible precession is small: {domega_per_orbit_incomp:.6e} rad/orbit")
    print(f"Compressible precession is: {domega_per_orbit_comp:.6e} rad/orbit")
    print()
    if abs(domega_per_orbit_comp) > abs(GR_precession) * 100:
        print("Compressible precession is ~{:.0f}× larger than GR!".format(
            abs(domega_per_orbit_comp / GR_precession)))
        print("This suggests an error in the compressible force implementation.")

print()
print("=" * 70)
print("TIMESTEP CONVERGENCE CHECK (incompressible)")
print("=" * 70)
print()

refinement_factor = 40
dt_refined = dt / refinement_factor
n_steps_refined = n_steps * refinement_factor
opts_incomp_refined = opts_incomp.copy()
opts_incomp_refined['save_every'] = max(1, n_steps_refined // 200)

print(
    "Refining timestep by ×{:.0f}: dt = {:.3e} yr, steps per orbit ≈ {:.0f}".format(
        refinement_factor, dt_refined, steps_per_orbit * refinement_factor
    )
)

bodies_refined = [copy.deepcopy(sun), copy.deepcopy(mercury)]
traj_refined, _ = integrate_orbit(bodies_refined, medium, dt_refined, n_steps_refined, opts_incomp_refined)

analysis_refined = analyse_precession(traj_refined, params)
domega_per_orbit_refined = analysis_refined['slope_per_orbit']
peri_per_orbit_refined = analysis_refined['peri_per_orbit']
peri_std_refined = analysis_refined['peri_std']
n_peri_refined = analysis_refined['n_peri_cycles']

print(
    f"Refined precession (incompressible, slope fit): {domega_per_orbit_refined:.6e} rad/orbit"
)
print(
    f"                                           = {domega_per_orbit_refined * 206265:.6e} arcsec/orbit"
)
if np.isfinite(peri_per_orbit_refined):
    err_arcsec = peri_std_refined * 206265.0
    print(
        f"Periapsis-to-periapsis: {peri_per_orbit_refined:.6e} rad/orbit"
        f" ± {peri_std_refined:.2e} (n = {n_peri_refined})"
    )
    print(
        f"                       = {peri_per_orbit_refined * 206265:.6e} arcsec/orbit"
        f" ± {err_arcsec:.2e}"
    )
if domega_per_orbit_refined != 0:
    ratio = abs(domega_per_orbit_incomp / domega_per_orbit_refined)
else:
    ratio = np.inf
print(f"Ratio coarse/refined ≈ {ratio:.1f}")
if (
    np.isfinite(peri_per_orbit_incomp)
    and np.isfinite(peri_per_orbit_refined)
    and peri_per_orbit_refined != 0
):
    ratio_peri = abs(peri_per_orbit_incomp / peri_per_orbit_refined)
    print(f"Ratio coarse/refined (peri) ≈ {ratio_peri:.1f}")
print()

print("=" * 70)
print("ORBIT ELEMENT EVOLUTION")
print("=" * 70)
print()

# Check if the orbit elements are stable
if len(t_comp) > 5:
    print("Orbit elements over time (compressible case):")
    print(f"{'Time [yr]':<12} {'a [AU]':<12} {'e':<12} {'omega [deg]':<12}")
    print("-" * 48)
    for i in [0, len(t_comp)//4, len(t_comp)//2, 3*len(t_comp)//4, -1]:
        elems = elems_comp[i]
        print(f"{t_comp[i]:<12.4f} {elems['a']:<12.6f} {elems['e']:<12.6f} {np.degrees(elems['omega']):<12.2f}")
