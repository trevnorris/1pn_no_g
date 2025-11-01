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
from slab.medium import Medium
from slab.bodies import Body
from slab.dynamics import integrate_orbit
from slab.gr1pn import compute_orbit_elements


def relative_orbit_series(traj, M_total, K):
    """Return times, periapsis angles, and orbital elements for the relative orbit."""

    t = traj['t']
    x_sun = traj['x'][:, 0, :]
    x_mercury = traj['x'][:, 1, :]
    v_sun = traj['v'][:, 0, :]
    v_mercury = traj['v'][:, 1, :]

    omega = []
    elems_list = []
    for xs, xm, vs, vm in zip(x_sun, x_mercury, v_sun, v_mercury):
        x_rel = xm - xs
        v_rel = vm - vs
        elems = compute_orbit_elements(x_rel, v_rel, M_total, K)
        omega.append(elems['omega'])
        elems_list.append(elems)

    return t, np.array(omega), elems_list

# Create Mercury configuration (barycentric initial conditions)
rho0 = 1.0
beta0 = 0.045
cs = 63239.7263
a = 0.387
e = 0.1

medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
K = medium.K

M_sun = 1.0
M_mercury = 3.3e-7
M_total = M_sun + M_mercury

r_peri = a * (1 - e)
mu = K * M_total
v_rel_peri = np.sqrt(mu * (1 + e) / (a * (1 - e)))

# Place the Sun and Mercury on opposite sides of the barycentre so that the
# centre-of-mass remains at rest.  This avoids an artificial drift that can be
# misinterpreted as apsidal motion when orbital elements are extracted from the
# lab-frame trajectory.
sun = Body(
    name="Sun",
    M=M_sun,
    x=np.array([-M_mercury / M_total * r_peri, 0.0, 0.0]),
    v=np.array([0.0, -M_mercury / M_total * v_rel_peri, 0.0]),
    R=0.001,
    Q=M_sun / beta0,
)

mercury = Body(
    name="Mercury",
    M=M_mercury,
    x=np.array([M_sun / M_total * r_peri, 0.0, 0.0]),
    v=np.array([0.0, M_sun / M_total * v_rel_peri, 0.0]),
    R=0.0005,
    Q=M_mercury / beta0,
)

T_orbit = 2 * np.pi * np.sqrt(a**3 / (K * M_total))

bodies = [sun, mercury]

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
    'save_every': 100,
    'verbose': False,
}

traj_incomp, diag_incomp = integrate_orbit(bodies_incomp, medium, dt, n_steps, opts_incomp)

# Extract relative orbital elements
t_incomp, omega_incomp, elems_incomp = relative_orbit_series(traj_incomp, M_total, K)
omega_incomp_unwrapped = np.unwrap(omega_incomp)

if len(t_incomp) > 10:
    coeffs = np.polyfit(t_incomp, omega_incomp_unwrapped, deg=1)
    domega_dt_incomp = coeffs[0]
    domega_per_orbit_incomp = domega_dt_incomp * T_orbit
else:
    domega_per_orbit_incomp = np.nan

print(f"Precession (incompressible): {domega_per_orbit_incomp:.6e} rad/orbit")
print(f"                           = {domega_per_orbit_incomp * 206265:.6e} arcsec/orbit")
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
    'save_every': 100,
    'verbose': False,
    'n_points': 128,
}

traj_comp, diag_comp = integrate_orbit(bodies_comp, medium, dt, n_steps, opts_comp)

# Extract relative orbital elements
t_comp, omega_comp, elems_comp = relative_orbit_series(traj_comp, M_total, K)
omega_comp_unwrapped = np.unwrap(omega_comp)

if len(t_comp) > 10:
    coeffs = np.polyfit(t_comp, omega_comp_unwrapped, deg=1)
    domega_dt_comp = coeffs[0]
    domega_per_orbit_comp = domega_dt_comp * T_orbit
else:
    domega_per_orbit_comp = np.nan

print(f"Precession (compressible): {domega_per_orbit_comp:.6e} rad/orbit")
print(f"                         = {domega_per_orbit_comp * 206265:.6e} arcsec/orbit")
print()

# GR prediction
GR_precession = (6 * np.pi * K * M_total) / (a * cs**2 * (1 - e**2))
print(f"GR prediction: {GR_precession:.6e} rad/orbit")
print(f"             = {GR_precession * 206265:.6e} arcsec/orbit")
print()

# Difference between compressible and incompressible
delta_precession = domega_per_orbit_comp - domega_per_orbit_incomp
print("=" * 70)
print("ANALYSIS")
print("=" * 70)
print()

print(f"Precession difference (comp - incomp): {delta_precession:.6e} rad/orbit")
print(f"                                     = {delta_precession * 206265:.6e} arcsec/orbit")
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
    print(f"  Δ(Δω) = {delta_precession:.6e} rad/orbit")
    print(f"        = {delta_precession * 206265:.6e} arcsec/orbit")
    if abs(delta_precession) < abs(domega_per_orbit_incomp) * 0.1:
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

t_refined, omega_refined, _ = relative_orbit_series(traj_refined, M_total, K)
omega_refined_unwrapped = np.unwrap(omega_refined)
coeffs_refined = np.polyfit(t_refined, omega_refined_unwrapped, deg=1)
domega_dt_refined = coeffs_refined[0]
domega_per_orbit_refined = domega_dt_refined * T_orbit

print(f"Refined precession (incompressible): {domega_per_orbit_refined:.6e} rad/orbit")
print(f"                                     = {domega_per_orbit_refined * 206265:.6e} arcsec/orbit")
if domega_per_orbit_refined != 0:
    ratio = abs(domega_per_orbit_incomp / domega_per_orbit_refined)
else:
    ratio = np.inf
print(f"Ratio coarse/refined ≈ {ratio:.1f}")
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
