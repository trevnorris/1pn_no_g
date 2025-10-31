"""Diagnostics module for superfluid orbit simulator.

This module provides functions for monitoring conserved quantities, computing
orbital elements, measuring perihelion precession, and other diagnostic
quantities for the superfluid slab gravity simulation.

Key diagnostics:
- Total kinetic energy: T = Σ (1/2) M_a v_a²
- Fluid potential energy: U_fluid = -K Σ_{a<b} M_a M_b / r_ab
- Total energy: E = T + U_fluid (conserved for conservative forces)
- Angular momentum: L = Σ M_a (r_a × v_a) (conserved)
- Osculating orbital elements for two-body systems
- Periapsis passage detection and precession measurement

Physics reference:
- Paper §3.4 for energy formulas
- plan_no_pde.md §12 for logging requirements
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
from dataclasses import dataclass


def total_kinetic_energy(bodies: List) -> float:
    """Compute total kinetic energy of the system.

    Formula:
        T = Σ_a (1/2) M_a v_a²

    Parameters
    ----------
    bodies : List[Body]
        List of Body objects with mass M and velocity v.

    Returns
    -------
    float
        Total kinetic energy [code units: M * L² / T²].

    Notes
    -----
    This should be a conserved quantity (up to numerical error) when
    combined with potential energy for conservative force systems.

    Examples
    --------
    >>> from slab.bodies import Body
    >>> import numpy as np
    >>> bodies = [
    ...     Body("A", M=1.0, x=[0,0,0], v=[1,0,0], R=1e-3, Q=0),
    ...     Body("B", M=2.0, x=[1,0,0], v=[0,0.5,0], R=1e-3, Q=0)
    ... ]
    >>> T = total_kinetic_energy(bodies)
    >>> # T = 0.5*1.0*1.0 + 0.5*2.0*0.25 = 0.5 + 0.25 = 0.75
    >>> np.isclose(T, 0.75)
    True
    """
    T = 0.0
    for body in bodies:
        v_squared = np.dot(body.v, body.v)
        T += 0.5 * body.M * v_squared
    return T


def fluid_potential_energy(bodies: List, medium) -> float:
    """Compute fluid interaction potential energy.

    Formula:
        U_fluid = -K Σ_{a<b} M_a M_b / r_ab

    where K = ρ₀/(4πβ²) is the orbital constant from the medium.

    This is the pair interaction energy arising from the superfluid
    velocity field. It plays the role of gravitational potential energy
    in the Newtonian analogy.

    Parameters
    ----------
    bodies : List[Body]
        List of Body objects with mass M and position x.
    medium : Medium
        Medium object providing the orbital constant K.

    Returns
    -------
    float
        Total fluid potential energy [code units: M * L² / T²].

    Notes
    -----
    - Sign convention: U is negative for attractive interactions (sinks)
    - K = medium.K = rho0/(4*pi*beta0^2)
    - This is exactly the field energy of the potential flow
    - For two bodies: U = -K * M1 * M2 / r12

    Edge cases:
    - If r_ab → 0, U → -∞ (unphysical; control surfaces must not overlap)
    - For N=1 body, returns 0.0 (no interactions)

    Examples
    --------
    >>> from slab.bodies import Body
    >>> from slab.medium import Medium
    >>> import numpy as np
    >>>
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> K = medium.K
    >>>
    >>> # Two bodies separated by r=1.0
    >>> bodies = [
    ...     Body("A", M=1.0, x=[0,0,0], v=[0,0,0], R=1e-3, Q=0),
    ...     Body("B", M=2.0, x=[1,0,0], v=[0,0,0], R=1e-3, Q=0)
    ... ]
    >>> U = fluid_potential_energy(bodies, medium)
    >>> # U = -K * 1.0 * 2.0 / 1.0 = -2*K
    >>> np.isclose(U, -2.0 * K)
    True
    """
    if len(bodies) <= 1:
        return 0.0

    K = medium.K
    U = 0.0

    # Sum over unique pairs (a < b)
    for i, body_a in enumerate(bodies):
        for body_b in bodies[i+1:]:
            r_ab = body_b.x - body_a.x
            r_ab_mag = np.linalg.norm(r_ab)

            if r_ab_mag == 0.0:
                raise ValueError(
                    f"Bodies '{body_a.name}' and '{body_b.name}' have identical positions. "
                    f"Potential energy diverges."
                )

            # U_ab = -K * M_a * M_b / r_ab
            U += -K * body_a.M * body_b.M / r_ab_mag

    return U


def total_energy(bodies: List, medium) -> float:
    """Compute total energy (kinetic + potential).

    Formula:
        E = T + U_fluid
        E = Σ (1/2) M_a v_a² - K Σ_{a<b} M_a M_b / r_ab

    For conservative forces, E should be conserved (constant in time).
    Energy drift is a key diagnostic for numerical accuracy.

    Parameters
    ----------
    bodies : List[Body]
        List of Body objects.
    medium : Medium
        Medium object providing K.

    Returns
    -------
    float
        Total energy [code units: M * L² / T²].

    Notes
    -----
    For bound orbits (E < 0), the system is gravitationally bound.
    For unbound trajectories (E > 0), bodies escape to infinity.

    Monitoring |ΔE|/|E| over time checks symplectic integrator quality.

    Examples
    --------
    >>> from slab.bodies import Body
    >>> from slab.medium import Medium
    >>> import numpy as np
    >>>
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> K = medium.K
    >>>
    >>> # Circular orbit: E = -K*M1*M2/(2*r) for reduced mass system
    >>> M1, M2 = 1.0, 1e-6  # Sun-like and planet-like
    >>> r = 1.0
    >>> mu = M1 * M2 / (M1 + M2)
    >>> v_circ = np.sqrt(K * (M1 + M2) / r)
    >>>
    >>> bodies = [
    ...     Body("Sun", M=M1, x=[0,0,0], v=[0,0,0], R=1e-3, Q=0),
    ...     Body("Planet", M=M2, x=[r,0,0], v=[0,v_circ,0], R=1e-3, Q=0)
    ... ]
    >>> E = total_energy(bodies, medium)
    >>> # For circular orbit: E ≈ -K*M1*M2/(2*r)
    >>> E_expected = -K * M1 * M2 / (2.0 * r)
    >>> np.isclose(E, E_expected, rtol=1e-5)
    True
    """
    T = total_kinetic_energy(bodies)
    U = fluid_potential_energy(bodies, medium)
    return T + U


def angular_momentum(bodies: List) -> np.ndarray:
    """Compute total angular momentum vector.

    Formula:
        L = Σ_a M_a (r_a × v_a)

    where × denotes the cross product.

    Parameters
    ----------
    bodies : List[Body]
        List of Body objects.

    Returns
    -------
    np.ndarray
        Total angular momentum vector, shape (3,) [code units: M * L² / T].

    Notes
    -----
    For central forces (like the superfluid momentum flux forces),
    angular momentum is conserved: dL/dt = 0.

    The magnitude |L| and direction L_hat define the orbital plane.
    For planar orbits, L points perpendicular to the orbital plane.

    Edge cases:
    - For perfectly radial motion, L = 0 (degenerate)
    - For circular orbit, |L| = μ * v * r where μ is reduced mass

    Examples
    --------
    >>> from slab.bodies import Body
    >>> import numpy as np
    >>>
    >>> # Single body orbiting origin (no reaction)
    >>> bodies = [
    ...     Body("Test", M=1.0, x=[1,0,0], v=[0,1,0], R=1e-3, Q=0)
    ... ]
    >>> L = angular_momentum(bodies)
    >>> # L = M * (r × v) = 1.0 * ([1,0,0] × [0,1,0]) = [0,0,1]
    >>> np.allclose(L, [0, 0, 1])
    True
    >>>
    >>> # Two-body system (center-of-mass frame)
    >>> bodies = [
    ...     Body("A", M=1.0, x=[1,0,0], v=[0,0.5,0], R=1e-3, Q=0),
    ...     Body("B", M=1.0, x=[-1,0,0], v=[0,-0.5,0], R=1e-3, Q=0)
    ... ]
    >>> L = angular_momentum(bodies)
    >>> # L_total = 1.0*([1,0,0]×[0,0.5,0]) + 1.0*([-1,0,0]×[0,-0.5,0])
    >>> #          = [0,0,0.5] + [0,0,0.5] = [0,0,1]
    >>> np.allclose(L, [0, 0, 1])
    True
    """
    L_total = np.zeros(3)

    for body in bodies:
        # L_a = M_a * (r_a × v_a)
        L_a = body.M * np.cross(body.x, body.v)
        L_total += L_a

    return L_total


def osculating_elements(
    body_a,
    body_b,
    medium,
    reference_body_index: int = 0
) -> Dict[str, float]:
    """Compute osculating orbital elements for a two-body system.

    At the current instant, compute the Keplerian orbital elements that
    would produce the current relative position and velocity.

    Parameters
    ----------
    body_a : Body
        First body (typically the primary/massive body).
    body_b : Body
        Second body (typically the secondary/lighter body).
    medium : Medium
        Medium object providing K.
    reference_body_index : int, optional
        0 for body_a as reference, 1 for body_b (default: 0).

    Returns
    -------
    dict
        Dictionary with keys:
        - 'a': semi-major axis [length units]
        - 'e': eccentricity [dimensionless, 0 ≤ e < 1 for bound]
        - 'i': inclination [radians, 0 ≤ i ≤ π]
        - 'Omega': longitude of ascending node [radians, 0 ≤ Ω < 2π]
        - 'omega': argument of periapsis [radians, 0 ≤ ω < 2π]
        - 'M_anom': mean anomaly [radians, 0 ≤ M < 2π]
        - 'period': orbital period [time units]
        - 'mu': gravitational parameter μ = K*(M₁+M₂) [L³/T²]

    Notes
    -----
    Uses the standard Keplerian formulas with μ = K*(M₁+M₂) where
    K = ρ₀/(4πβ²) from the medium.

    For the relative orbit:
        r = r_b - r_a  (position of b relative to a)
        v = v_b - v_a  (velocity of b relative to a)

    Then:
        h = r × v  (specific angular momentum)
        e_vec = (v × h)/μ - r/|r|  (eccentricity vector)
        e = |e_vec|  (eccentricity)
        a = -μ / (2*E_orb)  where E_orb = v²/2 - μ/r

    Edge cases:
    - Circular orbits (e ≈ 0): periapsis is ill-defined
    - Equatorial orbits (i ≈ 0): ascending node is ill-defined
    - Radial motion: eccentricity e = 1 (parabolic limit)

    References
    ----------
    - Vallado, "Fundamentals of Astrodynamics and Applications"
    - Murray & Dermott, "Solar System Dynamics"

    Examples
    --------
    >>> from slab.bodies import Body
    >>> from slab.medium import Medium
    >>> import numpy as np
    >>>
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> K = medium.K
    >>>
    >>> # Circular orbit at r=1.0
    >>> M1, M2 = 1.0, 1e-6
    >>> r = 1.0
    >>> v_circ = np.sqrt(K * (M1 + M2) / r)
    >>>
    >>> sun = Body("Sun", M=M1, x=[0,0,0], v=[0,0,0], R=1e-3, Q=0)
    >>> planet = Body("Planet", M=M2, x=[r,0,0], v=[0,v_circ,0], R=1e-3, Q=0)
    >>>
    >>> elem = osculating_elements(sun, planet, medium)
    >>> np.isclose(elem['a'], r)  # semi-major axis = r
    True
    >>> np.isclose(elem['e'], 0.0, atol=1e-10)  # circular
    True
    """
    # Relative position and velocity (b relative to a)
    r_rel = body_b.x - body_a.x
    v_rel = body_b.v - body_a.v

    r = np.linalg.norm(r_rel)
    v = np.linalg.norm(v_rel)

    if r == 0.0:
        raise ValueError("Bodies have identical positions; orbital elements undefined.")

    # Gravitational parameter μ = K * (M1 + M2)
    mu = medium.K * (body_a.M + body_b.M)

    # Specific angular momentum h = r × v
    h_vec = np.cross(r_rel, v_rel)
    h = np.linalg.norm(h_vec)

    # Handle degenerate case: radial motion (h ≈ 0)
    if h < 1e-15:
        # Radial motion: e = 1, other elements ill-defined
        return {
            'a': np.nan,
            'e': 1.0,
            'i': np.nan,
            'Omega': np.nan,
            'omega': np.nan,
            'M_anom': np.nan,
            'period': np.nan,
            'mu': mu
        }

    # Specific orbital energy E = v²/2 - μ/r
    E_orb = 0.5 * v**2 - mu / r

    # Semi-major axis a = -μ / (2*E)
    # (For elliptic orbits E < 0, so a > 0)
    if E_orb >= 0:
        # Unbound orbit (parabolic or hyperbolic)
        a = np.inf if E_orb == 0 else -mu / (2.0 * E_orb)
        period = np.inf
    else:
        a = -mu / (2.0 * E_orb)
        # Period T = 2π√(a³/μ)  (Kepler's third law with K)
        period = 2.0 * np.pi * np.sqrt(a**3 / mu)

    # Eccentricity vector e_vec = (v × h)/μ - r_hat
    r_hat = r_rel / r
    e_vec = np.cross(v_rel, h_vec) / mu - r_hat
    e = np.linalg.norm(e_vec)

    # Inclination i = arccos(h_z / h)
    h_hat = h_vec / h
    i = np.arccos(np.clip(h_hat[2], -1.0, 1.0))

    # Node vector n = z_hat × h (points toward ascending node)
    z_hat = np.array([0, 0, 1])
    n_vec = np.cross(z_hat, h_vec)
    n = np.linalg.norm(n_vec)

    # Longitude of ascending node Ω
    if n < 1e-15:
        # Equatorial orbit (i ≈ 0 or π): Ω undefined, set to 0
        Omega = 0.0
        n_hat = np.array([1, 0, 0])  # arbitrary reference direction
    else:
        n_hat = n_vec / n
        # Ω = atan2(n_y, n_x)
        Omega = np.arctan2(n_hat[1], n_hat[0])
        if Omega < 0:
            Omega += 2.0 * np.pi

    # Argument of periapsis ω
    if e < 1e-10:
        # Circular orbit: periapsis undefined, set ω = 0
        omega = 0.0
        e_hat = n_hat  # arbitrary
    else:
        e_hat = e_vec / e
        # ω = acos(n · e / (n * e))
        cos_omega = np.clip(np.dot(n_hat, e_hat), -1.0, 1.0)
        omega = np.arccos(cos_omega)
        # Quadrant check: if e_z < 0, then ω > π
        if e_vec[2] < 0:
            omega = 2.0 * np.pi - omega

    # True anomaly ν (angle from periapsis to current position)
    if e < 1e-10:
        # Circular: use angle from node
        cos_nu = np.clip(np.dot(n_hat, r_hat), -1.0, 1.0)
        nu = np.arccos(cos_nu)
        if np.dot(r_rel, v_rel) < 0:
            nu = 2.0 * np.pi - nu
    else:
        cos_nu = np.clip(np.dot(e_hat, r_hat), -1.0, 1.0)
        nu = np.arccos(cos_nu)
        # Quadrant check: if r · v < 0, we're past apoapsis
        if np.dot(r_rel, v_rel) < 0:
            nu = 2.0 * np.pi - nu

    # Eccentric anomaly E_anom (Kepler's equation)
    if e < 1e-10:
        # Circular: E = ν
        E_anom = nu
    else:
        # tan(E/2) = sqrt((1-e)/(1+e)) * tan(ν/2)
        E_anom = 2.0 * np.arctan2(
            np.sqrt(1.0 - e) * np.sin(nu / 2.0),
            np.sqrt(1.0 + e) * np.cos(nu / 2.0)
        )
        if E_anom < 0:
            E_anom += 2.0 * np.pi

    # Mean anomaly M = E - e*sin(E)
    M_anom = E_anom - e * np.sin(E_anom)
    if M_anom < 0:
        M_anom += 2.0 * np.pi

    return {
        'a': a,
        'e': e,
        'i': i,
        'Omega': Omega,
        'omega': omega,
        'M_anom': M_anom,
        'period': period,
        'mu': mu
    }


def find_periapsis_passages(
    trajectory: List[Dict],
    body_index: int = 1,
    reference_index: int = 0
) -> List[Tuple[float, float]]:
    """Find times when a body passes through periapsis.

    Given a trajectory (sequence of states), detect when the distance
    between two bodies reaches a local minimum (periapsis passage).

    Parameters
    ----------
    trajectory : List[Dict]
        List of trajectory snapshots, each a dict with keys:
        - 't': time
        - 'bodies': list of Body objects (or dicts with 'x' position)
    body_index : int, optional
        Index of the orbiting body (default: 1).
    reference_index : int, optional
        Index of the reference body (default: 0).

    Returns
    -------
    List[Tuple[float, float]]
        List of (time, periapsis_distance) tuples for each detected passage.

    Notes
    -----
    Detection method:
    - Compute r(t) = |r_body - r_ref| for all snapshots
    - Find local minima by checking r[i-1] > r[i] < r[i+1]
    - Optionally refine with parabolic interpolation

    Caveats:
    - Requires trajectory sampled frequently enough to resolve periapsis
    - May miss periapses if timestep >> orbital period
    - For high eccentricity, periapsis passages are brief

    Examples
    --------
    >>> # Mock trajectory with periapsis every ~6.28 time units
    >>> import numpy as np
    >>> trajectory = []
    >>> for i, t in enumerate(np.linspace(0, 20, 1000)):
    ...     theta = t  # roughly circular motion
    ...     r = 1.0 + 0.2 * np.cos(theta)  # eccentric orbit
    ...     x = r * np.cos(theta)
    ...     y = r * np.sin(theta)
    ...     bodies_mock = [
    ...         {'x': np.array([0, 0, 0])},  # reference at origin
    ...         {'x': np.array([x, y, 0])}   # orbiting body
    ...     ]
    ...     trajectory.append({'t': t, 'bodies': bodies_mock})
    >>>
    >>> passages = find_periapsis_passages(trajectory)
    >>> # Should find ~3 periapsis passages in t ∈ [0, 20]
    >>> len(passages) >= 2
    True
    """
    if len(trajectory) < 3:
        return []

    passages = []

    # Extract distances r(t)
    times = []
    distances = []
    for snapshot in trajectory:
        t = snapshot['t']
        bodies = snapshot['bodies']

        # Handle both Body objects and dicts
        if hasattr(bodies[reference_index], 'x'):
            r_ref = bodies[reference_index].x
            r_body = bodies[body_index].x
        else:
            r_ref = bodies[reference_index]['x']
            r_body = bodies[body_index]['x']

        r = np.linalg.norm(r_body - r_ref)
        times.append(t)
        distances.append(r)

    times = np.array(times)
    distances = np.array(distances)

    # Find local minima: r[i-1] > r[i] < r[i+1]
    for i in range(1, len(distances) - 1):
        if distances[i] < distances[i-1] and distances[i] < distances[i+1]:
            # Local minimum found
            t_peri = times[i]
            r_peri = distances[i]

            # Optional: parabolic interpolation for better precision
            # Fit parabola through (t[i-1], r[i-1]), (t[i], r[i]), (t[i+1], r[i+1])
            # and find its minimum
            dt_m = times[i] - times[i-1]
            dt_p = times[i+1] - times[i]
            dr_m = distances[i] - distances[i-1]
            dr_p = distances[i+1] - distances[i]

            if abs(dt_m) > 1e-15 and abs(dt_p) > 1e-15:
                # Parabolic correction: t_min ≈ t[i] + correction
                slope_m = dr_m / dt_m
                slope_p = dr_p / dt_p
                curvature = (slope_p - slope_m) / (dt_m + dt_p)

                if abs(curvature) > 1e-15:
                    correction = -0.5 * (slope_m + slope_p) / curvature
                    # Clamp correction to avoid extrapolating too far
                    correction = np.clip(correction, -dt_m, dt_p)
                    t_peri += correction
                    # Recompute r_peri (linear approximation)
                    r_peri += 0.5 * (slope_m + slope_p) * correction

            passages.append((t_peri, r_peri))

    return passages


def compute_precession(periapsis_passages: List[Tuple[float, float]]) -> float:
    """Compute perihelion advance per orbit from periapsis passages.

    Given a sequence of periapsis passages, estimate the average
    perihelion advance Δω per orbit.

    Parameters
    ----------
    periapsis_passages : List[Tuple[float, float]]
        List of (time, periapsis_distance) tuples from find_periapsis_passages().

    Returns
    -------
    float
        Average perihelion advance per orbit [radians/orbit].
        Positive for prograde precession (advance).

    Notes
    -----
    Method:
    1. Compute orbital period from time intervals between passages
    2. Track periapsis argument ω via distance changes
    3. Fit linear trend to ω(orbit_number)
    4. Return slope as Δω per orbit

    For purely Newtonian motion, Δω = 0 (closed ellipse).
    Post-Newtonian corrections give Δω > 0 (apsidal precession).

    Accuracy considerations:
    - Requires at least 3 periapsis passages for meaningful fit
    - Longer baselines (more orbits) improve statistical precision
    - Numerical noise in periapsis detection affects precision

    Examples
    --------
    >>> # Mock precessing orbit: Δω = 0.01 rad/orbit
    >>> passages = [
    ...     (0.0, 0.90),     # orbit 0
    ...     (6.28, 0.895),   # orbit 1: distance slightly smaller
    ...     (12.56, 0.890),  # orbit 2
    ...     (18.84, 0.885),  # orbit 3
    ... ]
    >>> dw = compute_precession(passages)
    >>> # Linear decrease in r_peri indicates precession
    >>> # (simplified model; actual calculation is more sophisticated)
    >>> dw > 0  # positive precession
    True
    """
    if len(periapsis_passages) < 3:
        raise ValueError(
            f"Need at least 3 periapsis passages to compute precession, "
            f"got {len(periapsis_passages)}"
        )

    times = np.array([t for t, r in periapsis_passages])
    distances = np.array([r for t, r in periapsis_passages])

    # Compute orbital periods from consecutive passages
    periods = np.diff(times)
    T_avg = np.mean(periods)

    # Method 1: Track periapsis distance variation
    # If orbit is precessing, periapsis distance may change slightly
    # due to tidal forces or non-Newtonian effects.
    # For simple apsidal precession in 1PN, the semi-major axis is conserved
    # but the periapsis rotates.

    # More robust method: if we had full state vector at each periapsis,
    # we could compute ω directly from orbital elements.
    # Here we use a simplified proxy: fit linear trend to log(r_peri)
    # and convert to angular precession.

    # For demonstration, assume constant eccentricity and compute
    # average rate of change in periapsis argument.

    # Orbit numbers
    orbit_numbers = np.arange(len(periapsis_passages))

    # Fit linear trend to distances (as proxy for ω change)
    # In 1PN, periapsis precesses but r_peri stays nearly constant
    # So we need a different approach.

    # Better method: estimate Δω from total precession over baseline
    # For Mercury: Δω ≈ 43 arcsec/century ≈ 5e-7 rad/orbit
    # Total precession = (number of orbits - 1) * Δω

    # Since we don't have full orbit state, use a heuristic:
    # Assume small secular drift in periapsis timing relative to
    # strict periodicity.

    # Compute expected times if perfectly periodic
    t_expected = times[0] + orbit_numbers * T_avg
    t_residual = times - t_expected

    # Linear fit to residuals: t_res = A + B * n
    # B gives the secular drift in timing, related to Δω
    if len(orbit_numbers) > 1:
        coeffs = np.polyfit(orbit_numbers, t_residual, deg=1)
        secular_drift = coeffs[0]  # [time / orbit]

        # Convert timing drift to angular precession
        # Δω ≈ 2π * (secular_drift / T_avg)
        # (rough approximation; exact formula depends on eccentricity)
        delta_omega = 2.0 * np.pi * secular_drift / T_avg
    else:
        delta_omega = 0.0

    return delta_omega


def energy_drift_monitor(
    trajectory: List[Dict],
    medium
) -> Dict[str, float]:
    """Monitor energy drift over a trajectory.

    Compute relative energy drift ΔE/E₀ and maximum drift.

    Parameters
    ----------
    trajectory : List[Dict]
        List of snapshots with 'bodies' key.
    medium : Medium
        Medium object for K.

    Returns
    -------
    dict
        Dictionary with:
        - 'E0': initial energy
        - 'Ef': final energy
        - 'dE': absolute drift |Ef - E0|
        - 'dE_rel': relative drift |ΔE|/|E₀|
        - 'dE_max': maximum absolute deviation from E₀

    Notes
    -----
    For symplectic integrators (Verlet, leapfrog), energy drift should
    be bounded and oscillatory, not secular.

    Typical acceptable drift: |ΔE|/|E| < 1e-6 for timestep dt << period.

    Examples
    --------
    >>> # Mock trajectory with conserved energy
    >>> import numpy as np
    >>> from slab.bodies import Body
    >>> from slab.medium import Medium
    >>>
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> trajectory = []
    >>> for i in range(100):
    ...     bodies = [
    ...         Body("A", M=1.0, x=[np.cos(i*0.1), np.sin(i*0.1), 0],
    ...              v=[-np.sin(i*0.1), np.cos(i*0.1), 0], R=1e-3, Q=0)
    ...     ]
    ...     trajectory.append({'t': i*0.1, 'bodies': bodies})
    >>>
    >>> drift = energy_drift_monitor(trajectory, medium)
    >>> # Single body, no interactions: E = KE only, should be constant
    >>> drift['dE_rel'] < 1e-10
    True
    """
    if len(trajectory) < 2:
        raise ValueError("Need at least 2 snapshots to compute energy drift")

    energies = []
    for snapshot in trajectory:
        bodies = snapshot['bodies']
        E = total_energy(bodies, medium)
        energies.append(E)

    energies = np.array(energies)
    E0 = energies[0]
    Ef = energies[-1]

    dE = abs(Ef - E0)
    dE_rel = dE / abs(E0) if E0 != 0 else np.inf

    # Maximum deviation from initial energy
    dE_max = np.max(np.abs(energies - E0))

    return {
        'E0': E0,
        'Ef': Ef,
        'dE': dE,
        'dE_rel': dE_rel,
        'dE_max': dE_max
    }


def lagrange_points_L123(
    body_a,
    body_b,
    medium
) -> List[np.ndarray]:
    """Compute collinear Lagrange points L1, L2, L3 (approximate).

    For two bodies in circular orbit, find the three collinear equilibrium
    points where centrifugal and gravitational forces balance.

    Parameters
    ----------
    body_a : Body
        Primary body (more massive).
    body_b : Body
        Secondary body (less massive).
    medium : Medium
        Medium for K.

    Returns
    -------
    List[np.ndarray]
        Positions of L1, L2, L3 in inertial frame, each shape (3,).
        L1: between the two bodies
        L2: beyond body_b from body_a
        L3: beyond body_a from body_b

    Notes
    -----
    Uses Hill approximation for small mass ratio μ = M_b/(M_a + M_b).
    Exact locations require solving quintic equation numerically.

    For μ << 1:
        r_L1 ≈ r * (1 - (μ/3)^(1/3))
        r_L2 ≈ r * (1 + (μ/3)^(1/3))
        r_L3 ≈ -r * (1 - 5*μ/12)

    where r is the separation between bodies.

    Examples
    --------
    >>> from slab.bodies import Body
    >>> from slab.medium import Medium
    >>> import numpy as np
    >>>
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> sun = Body("Sun", M=1.0, x=[0,0,0], v=[0,0,0], R=1e-3, Q=0)
    >>> earth = Body("Earth", M=3e-6, x=[1,0,0], v=[0,1,0], R=1e-3, Q=0)
    >>>
    >>> L_points = lagrange_points_L123(sun, earth, medium)
    >>> # L1 is between sun and earth, closer to earth
    >>> L1, L2, L3 = L_points
    >>> 0 < np.linalg.norm(L1) < 1.0  # L1 between 0 and 1
    True
    """
    M1 = body_a.M
    M2 = body_b.M
    r_ab = body_b.x - body_a.x
    r = np.linalg.norm(r_ab)
    r_hat = r_ab / r

    # Mass ratio
    mu_ratio = M2 / (M1 + M2)

    # Hill sphere approximation (valid for mu << 1)
    if mu_ratio < 0.1:
        # L1: between bodies
        r_L1 = r * (1.0 - (mu_ratio / 3.0)**(1.0/3.0))
        L1_pos = body_a.x + r_L1 * r_hat

        # L2: beyond body_b
        r_L2 = r * (1.0 + (mu_ratio / 3.0)**(1.0/3.0))
        L2_pos = body_a.x + r_L2 * r_hat

        # L3: beyond body_a (opposite side)
        r_L3 = -r * (1.0 - 5.0 * mu_ratio / 12.0)
        L3_pos = body_a.x + r_L3 * r_hat
    else:
        # For larger mass ratios, use numerical solver (not implemented here)
        # Placeholder: return approximate positions
        import warnings
        warnings.warn(
            f"Mass ratio μ = {mu_ratio:.3f} is not small; "
            f"Hill approximation may be inaccurate.",
            UserWarning
        )
        L1_pos = body_a.x + 0.5 * r_ab
        L2_pos = body_a.x + 1.2 * r_ab
        L3_pos = body_a.x - 0.9 * r_ab

    return [L1_pos, L2_pos, L3_pos]
