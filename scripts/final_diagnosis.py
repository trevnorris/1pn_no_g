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

print("DIAGNOSTIC: Starting final_diagnosis.py", flush=True)
print("DIAGNOSTIC: Importing modules...", flush=True)

import numpy as np
import copy
from slab.dynamics import integrate_orbit

print("DIAGNOSTIC: Importing precession helpers...", flush=True)
from scripts.precession_helpers import (
    analyse_precession,
    create_barycentric_mercury_config,
)

print("DIAGNOSTIC: All imports complete!", flush=True)

# Create Mercury configuration (barycentric initial conditions)
print("DIAGNOSTIC: Creating Mercury configuration...", flush=True)
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

print("DIAGNOSTIC: Configuration created successfully!", flush=True)
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

# Integration parameters: Test if quadrature resolution (n_points) is more
# important than timestep resolution for capturing 1PN physics.
#
# Hypothesis: Surface integral accuracy >> timestep accuracy for 1PN signal
#
# Using dt=1×10⁻⁴ yr (2413 steps/orbit, still 20× better than original)
# with n_points=256 (2× higher quadrature than previous 85% result)
# Expected runtime: ~15-20 minutes
dt = 1.0e-4  # yr - moderate timestep
n_orbits_target = 8  # Number of orbits to integrate
total_time = n_orbits_target * T_orbit
n_steps = int(total_time / dt)
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
    'save_every': max(1, n_steps // 500),  # Save ~500 points (enough for analysis)
    'verbose': False,
}

print(f"DIAGNOSTIC: Starting incompressible integration ({n_steps} steps)...", flush=True)
print(f"DIAGNOSTIC: This may take ~30 seconds...", flush=True)
traj_incomp, diag_incomp = integrate_orbit(bodies_incomp, medium, dt, n_steps, opts_incomp)
print(f"DIAGNOSTIC: Incompressible integration complete!", flush=True)
print(f"DIAGNOSTIC: Trajectory has {len(traj_incomp['t'])} saved points", flush=True)
print(f"DIAGNOSTIC: Starting precession analysis (computing orbital elements)...", flush=True)

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
print("RUN 2: WITH COMPRESSIBLE FORCES (same dt as incompressible)")
print("=" * 70)
print()

bodies_comp = [copy.deepcopy(sun), copy.deepcopy(mercury)]
opts_comp = {
    'use_compressible': True,
    'use_quadrature': False,
    'save_every': max(1, n_steps // 500),
    'verbose': False,
    'n_points': 256,  # Higher quadrature (testing if surface integral accuracy matters more)
}

print(f"DIAGNOSTIC: Starting compressible integration ({n_steps} steps, n_points=256)...", flush=True)
print(f"DIAGNOSTIC: Expected time: ~15-20 minutes...", flush=True)
print(f"DIAGNOSTIC: Testing hypothesis: quadrature resolution > timestep resolution for 1PN", flush=True)
traj_comp, diag_comp = integrate_orbit(bodies_comp, medium, dt, n_steps, opts_comp)
print(f"DIAGNOSTIC: Compressible integration complete!", flush=True)
print(f"DIAGNOSTIC: Trajectory has {len(traj_comp['t'])} saved points", flush=True)
print(f"DIAGNOSTIC: Starting precession analysis...", flush=True)

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
print(f"NOTE: Both runs use dt = {dt:.3e} yr (~{steps_per_orbit:.0f} steps/orbit)")
print(f"      Difference isolates the 1PN compressible correction from matched baselines.")
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
print("COARSE TIMESTEP COMPARISON (for reference)")
print("=" * 70)
print()
print("The main runs above used dt = {:.3e} yr (~{:.0f} steps/orbit).".format(dt, steps_per_orbit))
print("For comparison, here's the same run at coarse resolution (dt = 0.002 yr, ~121 steps/orbit):")
print()

# Run at the old coarse timestep for comparison
dt_coarse = 0.002
n_steps_coarse = int(total_time / dt_coarse)
steps_per_orbit_coarse = T_orbit / dt_coarse
opts_incomp_coarse = opts_incomp.copy()
opts_incomp_coarse['save_every'] = max(1, n_steps_coarse // 2000)

print(
    "Coarse timestep: dt = {:.3e} yr, steps per orbit ≈ {:.0f}".format(
        dt_coarse, steps_per_orbit_coarse
    )
)

bodies_coarse = [copy.deepcopy(sun), copy.deepcopy(mercury)]
traj_coarse, _ = integrate_orbit(bodies_coarse, medium, dt_coarse, n_steps_coarse, opts_incomp_coarse)

analysis_coarse = analyse_precession(traj_coarse, params)
domega_per_orbit_coarse = analysis_coarse['slope_per_orbit']
peri_per_orbit_coarse = analysis_coarse['peri_per_orbit']
peri_std_coarse = analysis_coarse['peri_std']
n_peri_coarse = analysis_coarse['n_peri_cycles']

print(
    f"Coarse precession (incompressible, slope fit): {domega_per_orbit_coarse:.6e} rad/orbit"
)
print(
    f"                                          = {domega_per_orbit_coarse * 206265:.6e} arcsec/orbit"
)
if np.isfinite(peri_per_orbit_coarse):
    err_arcsec = peri_std_coarse * 206265.0
    print(
        f"Periapsis-to-periapsis: {peri_per_orbit_coarse:.6e} rad/orbit"
        f" ± {peri_std_coarse:.2e} (n = {n_peri_coarse})"
    )
    print(
        f"                       = {peri_per_orbit_coarse * 206265:.6e} arcsec/orbit"
        f" ± {err_arcsec:.2e}"
    )
if domega_per_orbit_coarse != 0:
    ratio = abs(domega_per_orbit_coarse / domega_per_orbit_incomp)
else:
    ratio = np.inf
print(f"\nRatio fine/coarse ≈ {1.0/ratio:.1f} (fine timestep reduces artifact by {ratio:.1f}×)")
if (
    np.isfinite(peri_per_orbit_incomp)
    and np.isfinite(peri_per_orbit_coarse)
    and peri_per_orbit_coarse != 0
):
    ratio_peri = abs(peri_per_orbit_coarse / peri_per_orbit_incomp)
    print(f"Ratio fine/coarse (peri) ≈ {1.0/ratio_peri:.1f}")
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
