#!/usr/bin/env python3
"""Fast precession sanity-check runner.

This tool integrates a two-body Sun–Mercury system with configurable timestep
and orbit count, reporting incompressible versus compressible perihelion
precession alongside the General Relativity 1PN prediction.  It is intended for
rapid iteration when tuning integrator settings: the defaults resolve ~5 orbits
with ~4800 steps/orbit, which is sufficient to see whether the numerical drift
has collapsed before running the full validation sweep.
"""

from __future__ import annotations

import argparse
import copy
import math
from typing import Dict

import numpy as np

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.precession_helpers import (
    analyse_precession,
    create_barycentric_mercury_config,
)


def run_single_case(
    *,
    dt: float,
    n_steps: int,
    use_compressible: bool,
    save_every: int,
    n_points: int,
    medium,
    bodies,
    params,
) -> Dict[str, float]:
    from slab.dynamics import integrate_orbit

    opts = {
        "use_compressible": use_compressible,
        "save_every": save_every,
        "use_quadrature": False,
        "verbose": False,
        "n_points": n_points,
        "audit_every": 0,
    }
    traj, _ = integrate_orbit([copy.deepcopy(b) for b in bodies], medium, dt, n_steps, opts)
    return analyse_precession(traj, params)


def format_precession(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    return f"{value:.6e} rad/orbit = {value * 206265:.6e} arcsec/orbit"


def main() -> None:
    parser = argparse.ArgumentParser(description="Quick incompressible vs compressible precession check")
    parser.add_argument("--eccentricity", type=float, default=0.1, help="Orbital eccentricity (default: 0.1)")
    parser.add_argument(
        "--dt",
        type=float,
        default=5e-5,
        help="Timestep in years (default: 5e-5, ≈4800 steps/orbit)",
    )
    parser.add_argument(
        "--orbits",
        type=float,
        default=5.0,
        help="Number of orbital periods to integrate (default: 5)",
    )
    parser.add_argument("--save-every", type=int, default=1, help="Trajectory save interval (default: 1)")
    parser.add_argument("--n-points", type=int, default=64, help="Quadrature points when audits trigger (default: 64)")
    parser.add_argument("--rho0", type=float, default=1.0, help="Background density (default: 1.0)")
    parser.add_argument("--beta0", type=float, default=0.045, help="Mass-intake factor (default: 0.045)")
    parser.add_argument(
        "--cs",
        type=float,
        default=63239.7263,
        help="Sound speed in AU/yr (default: 63239.7263 ≈ c)",
    )
    parser.add_argument(
        "--semi-major-axis",
        type=float,
        default=0.387,
        help="Semi-major axis in AU (default: Mercury's 0.387)",
    )
    parser.add_argument(
        "--refine-factor",
        type=int,
        default=1,
        help="Optional extra run with dt/refinement (>=1). 1 disables refinement.",
    )
    args = parser.parse_args()

    medium, bodies, params = create_barycentric_mercury_config(
        eccentricity=args.eccentricity,
        rho0=args.rho0,
        beta0=args.beta0,
        cs=args.cs,
        semi_major_axis=args.semi_major_axis,
    )

    steps_per_orbit = params.T_orbit / args.dt
    total_time = args.orbits * params.T_orbit
    n_steps = math.ceil(total_time / args.dt)

    print("=" * 70)
    print("QUICK PRECESSION PROBE")
    print("=" * 70)
    print()
    print(f"a = {params.a:.6f} AU, e = {params.e:.3f}")
    print(f"T_orbit = {params.T_orbit:.6f} yr, cs = {args.cs:.4f} AU/yr")
    print(f"dt = {args.dt:.3e} yr → {steps_per_orbit:.1f} steps/orbit")
    print(f"Integrating {args.orbits:.2f} orbits = {total_time:.6f} yr → {n_steps} steps")
    print()

    save_every = max(1, args.save_every)

    results = {}
    for label, flag in [
        ("Incompressible", False),
        ("Compressible", True),
    ]:
        analysis = run_single_case(
            dt=args.dt,
            n_steps=n_steps,
            use_compressible=flag,
            save_every=save_every,
            n_points=args.n_points,
            medium=medium,
            bodies=bodies,
            params=params,
        )
        results[label] = analysis
        print(f"{label} run:")
        print(f"  slope-fit precession:   {format_precession(analysis['slope_per_orbit'])}")
        if np.isfinite(analysis["peri_per_orbit"]):
            print(
                "  peri-to-peri precession: "
                f"{format_precession(analysis['peri_per_orbit'])}"
                f" (n={analysis['n_peri_cycles']})"
            )
        print()

    delta = results["Compressible"]["slope_per_orbit"] - results["Incompressible"]["slope_per_orbit"]
    delta_arcsec = delta * 206265
    gr = (6 * np.pi * params.K * params.M_total) / (
        params.a * args.cs**2 * (1 - params.e**2)
    )
    print("Δω_compressible - Δω_incompressible:")
    print(f"  {delta:.6e} rad/orbit = {delta_arcsec:.6e} arcsec/orbit")
    if gr != 0:
        print(f"  GR 1PN expectation: {gr * 206265:.6e} arcsec/orbit")
        print(f"  Ratio to GR: {delta / gr:.4f}")
    print()

    if args.refine_factor > 1:
        refined_dt = args.dt / args.refine_factor
        refined_steps_per_orbit = params.T_orbit / refined_dt
        refined_steps = n_steps * args.refine_factor
        print("Refined incompressible rerun:")
        print(
            f"  dt = {refined_dt:.3e} yr → {refined_steps_per_orbit:.1f} steps/orbit,"
            f" n_steps = {refined_steps}"
        )
        refined = run_single_case(
            dt=refined_dt,
            n_steps=refined_steps,
            use_compressible=False,
            save_every=max(1, refined_steps // 5000),
            n_points=args.n_points,
            medium=medium,
            bodies=bodies,
            params=params,
        )
        print(f"  slope-fit precession: {format_precession(refined['slope_per_orbit'])}")
        if np.isfinite(refined["peri_per_orbit"]):
            print(
                "  peri-to-peri precession: "
                f"{format_precession(refined['peri_per_orbit'])}"
                f" (n={refined['n_peri_cycles']})"
            )
        coarse = results["Incompressible"]["slope_per_orbit"]
        if refined["slope_per_orbit"] != 0:
            print(f"  coarse/refined ratio: {coarse / refined['slope_per_orbit']:.1f}")
        print()


if __name__ == "__main__":
    main()
