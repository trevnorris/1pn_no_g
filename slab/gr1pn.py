"""GR 1PN (Einstein-Infeld-Hoffmann) comparison module.

This module implements standard General Relativity at the 1st post-Newtonian (1PN)
order for validation purposes ONLY. It is completely independent from the slab
simulation and serves to:

1. Validate that the slab can reproduce GR-like effects with appropriate parameters
2. Compare perihelion precession between slab and GR
3. Provide a reference implementation of EIH equations in harmonic gauge

**CRITICAL**: This module uses the gravitational constant G (not K from the slab).
It never feeds back into slab force calculations. The slab simulation is entirely
based on superfluid hydrodynamics and never imports this module for physics.

Theory
------
The 1PN acceleration in harmonic gauge (Einstein-Infeld-Hoffmann) for body i is:

    a_i = a_N + a_1PN/c²

where a_N is the Newtonian acceleration and a_1PN contains post-Newtonian corrections:

    a_N = Σ_{j≠i} G M_j n_ij / r_ij²

    a_1PN/c² = Σ_{j≠i} {
        + (4 G M_i/r_ij + 4 G M_j/r_ij - v_i² - 2 v_j²
           + 4 v_i·v_j + 3/2 (n_ij·v_j)²) (G M_j n_ij / r_ij²)
        + (4 n_ij·v_i - 3 n_ij·v_j) (G M_j v_ij / r_ij²)
    }

where:
    n_ij = (x_i - x_j) / r_ij    (unit vector from j to i)
    v_ij = v_i - v_j              (relative velocity)
    r_ij = |x_i - x_j|            (separation)

This is the standard 1PN formula from:
- Soffel (1989), "Relativity in Astrometry, Celestial Mechanics and Geodesy"
- Will & Wiseman (1996), "Gravitational radiation from compact binary systems"
- Blanchet (2014), "Gravitational Radiation from Post-Newtonian Sources"

Parameter Mapping
-----------------
To compare slab with GR:

**Slab parameters**:
    K = ρ₀/(4π β₀²)    [orbital constant, replaces G]
    cs                  [sound speed, controls 1PN-like corrections]

**GR parameters**:
    G = 6.6743e-11 m³/(kg·s²)    [gravitational constant]
    c = 299792458 m/s             [speed of light]

**Mapping for validation**:
    K_slab ↔ G_GR     [both control Newtonian orbital dynamics]
    cs_slab ↔ c_GR    [both control 1PN correction magnitude]

For a given set of initial conditions (masses, positions, velocities):
1. Run slab with (K, cs) chosen to match desired orbital properties
2. Run GR with (G, c) using standard Solar System values
3. Compare precession rates, orbital elements, etc.

The slab should match GR precession (up to model differences) when cs is chosen
appropriately, demonstrating that 1PN effects emerge from superfluid dynamics.

Examples
--------
>>> import numpy as np
>>> from slab.bodies import Body
>>>
>>> # Create Sun and Mercury for GR simulation
>>> sun = Body(
...     name="Sun",
...     M=1.989e30,  # kg
...     x=np.array([0.0, 0.0, 0.0]),
...     v=np.array([0.0, 0.0, 0.0]),
...     R=0.0,  # Not used in GR
...     Q=0.0   # Not used in GR
... )
>>>
>>> # Mercury: a = 5.79e10 m, v = 47.87 km/s
>>> mercury = Body(
...     name="Mercury",
...     M=3.301e23,  # kg
...     x=np.array([5.79e10, 0.0, 0.0]),
...     v=np.array([0.0, 4.787e4, 0.0]),
...     R=0.0,
...     Q=0.0
... )
>>>
>>> bodies = [sun, mercury]
>>> c_light = 299792458.0  # m/s
>>> G_newton = 6.6743e-11  # m³/(kg·s²)
>>>
>>> # Compute 1PN accelerations
>>> accels = eih_1pn_accel(bodies, c_light, G_newton)
>>> print(f"Mercury 1PN accel: {accels[1]}")
>>>
>>> # Integrate orbit
>>> dt = 100.0  # seconds
>>> n_steps = 1000000  # ~11.6 days
>>> trajectory = integrate_gr1pn_orbit(bodies, c_light, G_newton, dt, n_steps)
>>> print(f"Integrated {n_steps} steps, final time = {trajectory['t'][-1]/86400:.2f} days")
"""

from typing import List, Dict
import numpy as np
from dataclasses import dataclass

# Import Body from slab, but only use it as a data container
from slab.bodies import Body


def eih_1pn_accel(
    bodies: List[Body],
    c_light: float,
    G_newton: float = 6.6743e-11
) -> np.ndarray:
    """Compute Einstein-Infeld-Hoffmann 1PN accelerations in harmonic gauge.

    This implements the standard post-Newtonian expansion of General Relativity
    to O(1/c²), also known as the EIH equations. These are the reference GR
    predictions for comparison with the slab simulation.

    Parameters
    ----------
    bodies : List[Body]
        List of N bodies with masses M, positions x, velocities v.
        Only M, x, v fields are used; Q, R are ignored.
    c_light : float
        Speed of light [m/s in SI, or code units].
        Standard value: c = 299792458 m/s.
    G_newton : float, optional
        Gravitational constant [m³/(kg·s²) in SI, or code units].
        Standard value: G = 6.6743e-11 m³/(kg·s²).
        Default: 6.6743e-11 (SI value).

    Returns
    -------
    accels : np.ndarray
        Array of shape (N, 3) containing the 1PN acceleration vectors
        for each body [m/s² in SI, or code units].

    Notes
    -----
    **1PN formula (harmonic gauge)**:

    For body i, the acceleration is:

        a_i = Σ_{j≠i} [ (G M_j / r_ij²) n_ij ] * [ 1 + (1/c²) × (1PN terms) ]

    where the 1PN correction factor includes:

        (1/c²) × {
            + 4 G M_i / r_ij         [self-potential]
            + 4 G M_j / r_ij         [external potential, doubled]
            - v_i²                    [kinetic energy of i]
            - 2 v_j²                  [kinetic energy of j, doubled]
            + 4 v_i · v_j             [velocity coupling]
            + 3/2 (n_ij · v_j)²       [tangential velocity squared]
        }

    plus a velocity-dependent vector term:

        + (1/c²) × (G M_j / r_ij²) × [4 (n_ij·v_i) - 3 (n_ij·v_j)] v_ij

    This is the standard EIH formula. For derivation, see:
    - Soffel (1989), Eq. (6.82)
    - Will & Wiseman (1996), Eq. (2.1)
    - Blanchet (2014), Eq. (6.25)

    **Perihelion precession**:

    For a two-body orbit with M_1 >> M_2, semi-major axis a, eccentricity e:

        Δφ ≈ (6π G M_1) / (a c² (1 - e²))    [radians per orbit]

    Mercury: Δφ ≈ 5.0×10⁻⁷ rad/orbit ≈ 43 arcsec/century.

    **Harmonic gauge**:

    The EIH equations are valid in harmonic (de Donder) gauge, where the
    coordinate condition is ∂_μ(√(-g) g^{μν}) = 0. This is the most common
    gauge for post-Newtonian calculations and matches standard textbooks.

    **N-body case**:

    The formula generalizes naturally to N bodies via pairwise summation.
    Three-body and higher interactions enter at 2PN order (O(1/c⁴)).

    Examples
    --------
    >>> # Mercury orbit (SI units)
    >>> sun = Body("Sun", M=1.989e30, x=[0,0,0], v=[0,0,0], R=0, Q=0)
    >>> mercury = Body("Mercury", M=3.301e23,
    ...                x=[5.79e10, 0, 0], v=[0, 4.787e4, 0], R=0, Q=0)
    >>> bodies = [sun, mercury]
    >>>
    >>> c_si = 299792458.0  # m/s
    >>> G_si = 6.6743e-11   # m³/(kg·s²)
    >>> accels = eih_1pn_accel(bodies, c_si, G_si)
    >>>
    >>> # Mercury's acceleration magnitude (mostly Newtonian)
    >>> a_merc = np.linalg.norm(accels[1])
    >>> print(f"Mercury accel: {a_merc:.3e} m/s²")  # ~0.04 m/s²
    >>>
    >>> # 1PN correction is small: ~(v/c)² ~ (50 km/s / 300000 km/s)² ~ 3e-8
    >>> # so a_1PN / a_N ~ 3e-8

    References
    ----------
    .. [1] M. Soffel, "Relativity in Astrometry, Celestial Mechanics and
           Geodesy" (1989), Springer.
    .. [2] C. Will & A. Wiseman, "Gravitational radiation from compact binary
           systems: Computational methods", Phys. Rev. D 54, 4813 (1996).
    .. [3] L. Blanchet, "Gravitational Radiation from Post-Newtonian Sources
           and Inspiralling Compact Binaries", Living Rev. Relativity 17, 2 (2014).
    """
    N = len(bodies)
    c2 = c_light ** 2
    accels = np.zeros((N, 3))

    # Loop over all bodies i
    for i in range(N):
        M_i = bodies[i].M
        x_i = bodies[i].x
        v_i = bodies[i].v
        v_i_sq = np.dot(v_i, v_i)

        a_i = np.zeros(3)

        # Sum over all other bodies j ≠ i
        for j in range(N):
            if i == j:
                continue

            M_j = bodies[j].M
            x_j = bodies[j].x
            v_j = bodies[j].v

            # Separation vector: r_ij = x_i - x_j  (from j to i)
            r_ij_vec = x_i - x_j
            r_ij = np.linalg.norm(r_ij_vec)

            if r_ij == 0:
                raise ValueError(
                    f"Bodies {i} ('{bodies[i].name}') and {j} ('{bodies[j].name}') "
                    f"have identical positions: collision detected."
                )

            # Unit vector from j to i
            n_ij = r_ij_vec / r_ij

            # Relative velocity: v_ij = v_i - v_j
            v_ij = v_i - v_j

            # Velocity magnitudes and dot products
            v_j_sq = np.dot(v_j, v_j)
            v_i_dot_v_j = np.dot(v_i, v_j)
            n_ij_dot_v_i = np.dot(n_ij, v_i)
            n_ij_dot_v_j = np.dot(n_ij, v_j)

            # Newtonian acceleration from body j
            # a_N = (G M_j / r_ij²) n_ij
            a_N_prefactor = G_newton * M_j / (r_ij ** 2)
            a_N = a_N_prefactor * n_ij

            # --- 1PN correction terms (all O(1/c²)) ---

            # Scalar 1PN correction factor (multiplies a_N direction)
            # This is the coefficient in front of n_ij
            gamma_1pn = (
                + 4.0 * G_newton * M_i / r_ij       # self-potential
                + 4.0 * G_newton * M_j / r_ij       # external potential (×2)
                - v_i_sq                             # -v_i²
                - 2.0 * v_j_sq                       # -2 v_j²
                + 4.0 * v_i_dot_v_j                  # +4 v_i·v_j
                + 1.5 * (n_ij_dot_v_j ** 2)          # +3/2 (n_ij·v_j)²
            )

            # Scalar 1PN correction (radial)
            a_1pn_radial = (a_N_prefactor / c2) * gamma_1pn * n_ij

            # Vector 1PN correction (velocity-dependent term)
            # Coefficient: [4 (n_ij·v_i) - 3 (n_ij·v_j)]
            # Direction: v_ij = v_i - v_j
            velocity_coeff = 4.0 * n_ij_dot_v_i - 3.0 * n_ij_dot_v_j
            a_1pn_velocity = (a_N_prefactor / c2) * velocity_coeff * v_ij

            # Total acceleration from body j
            a_j = a_N + a_1pn_radial + a_1pn_velocity

            # Accumulate
            a_i += a_j

        accels[i] = a_i

    return accels


def integrate_gr1pn_orbit(
    bodies: List[Body],
    c_light: float,
    G_newton: float,
    dt: float,
    n_steps: int,
    save_every: int = 1
) -> Dict[str, np.ndarray]:
    """Integrate N-body orbits using GR 1PN forces (EIH equations).

    This is a standalone GR integrator for comparison with the slab simulation.
    It uses the same velocity-Verlet algorithm but with EIH 1PN accelerations
    instead of slab forces.

    Parameters
    ----------
    bodies : List[Body]
        List of N bodies. Only M, x, v are used and updated; Q, R are ignored.
        Bodies are modified in-place during integration.
    c_light : float
        Speed of light [same units as positions and velocities].
    G_newton : float
        Gravitational constant [compatible with M, x, v units].
    dt : float
        Time step [same units as time].
    n_steps : int
        Number of integration steps.
    save_every : int, optional
        Save trajectory data every N steps (default: 1 = every step).

    Returns
    -------
    trajectory : Dict[str, np.ndarray]
        Dictionary containing:
        - 't': time array, shape (n_saved,)
        - 'x': positions, shape (n_saved, N, 3)
        - 'v': velocities, shape (n_saved, N, 3)
        - 'a': accelerations, shape (n_saved, N, 3)

    Notes
    -----
    **Integration scheme**: Velocity-Verlet (symplectic leapfrog)

        1. v(t+dt/2) = v(t) + (1/2) a(t) dt
        2. x(t+dt) = x(t) + v(t+dt/2) dt
        3. Compute a(t+dt) from new positions x(t+dt)
        4. v(t+dt) = v(t+dt/2) + (1/2) a(t+dt) dt

    This is the same algorithm used in the slab integrator, ensuring that
    any differences are due to forces, not numerics.

    **Energy conservation**:

    At 1PN order, the conserved energy is:

        E = Σ_i (1/2) M_i v_i² - Σ_{i<j} G M_i M_j / r_ij
            + (1/c²) × (1PN corrections to kinetic and potential energy)

    For small v/c, energy should be conserved to ~(v/c)² accuracy by Verlet.

    **Stability**:

    Choose dt such that:
        dt << P_orb / (2π × N_steps_per_orbit)

    where P_orb is the orbital period. Typically N_steps_per_orbit ~ 100-1000.

    Examples
    --------
    >>> # Mercury orbit: ~88 days, sample 1000 steps per orbit
    >>> sun = Body("Sun", M=1.989e30, x=[0,0,0], v=[0,0,0], R=0, Q=0)
    >>> mercury = Body("Mercury", M=3.301e23,
    ...                x=[5.79e10, 0, 0], v=[0, 4.787e4, 0], R=0, Q=0)
    >>>
    >>> P_orb = 88 * 86400  # seconds
    >>> dt = P_orb / 1000   # ~7600 s per step
    >>> n_steps = 10000     # 10 orbits
    >>>
    >>> c_si = 299792458.0
    >>> G_si = 6.6743e-11
    >>> traj = integrate_gr1pn_orbit([sun, mercury], c_si, G_si, dt, n_steps)
    >>>
    >>> # Extract Mercury's trajectory
    >>> x_merc = traj['x'][:, 1, :]  # shape (n_steps, 3)
    >>> # Compute precession, etc.
    """
    N = len(bodies)
    n_saved = (n_steps // save_every) + 1

    # Allocate trajectory storage
    trajectory = {
        't': np.zeros(n_saved),
        'x': np.zeros((n_saved, N, 3)),
        'v': np.zeros((n_saved, N, 3)),
        'a': np.zeros((n_saved, N, 3)),
    }

    # Initial state
    t = 0.0
    save_idx = 0

    # Compute initial accelerations
    a_current = eih_1pn_accel(bodies, c_light, G_newton)

    # Save initial state
    trajectory['t'][save_idx] = t
    for i, body in enumerate(bodies):
        trajectory['x'][save_idx, i] = body.x.copy()
        trajectory['v'][save_idx, i] = body.v.copy()
        trajectory['a'][save_idx, i] = a_current[i].copy()
    save_idx += 1

    # Velocity-Verlet integration loop
    for step in range(1, n_steps + 1):
        # 1. Half-step velocity update: v(t+dt/2) = v(t) + (1/2) a(t) dt
        for i, body in enumerate(bodies):
            body.v += 0.5 * a_current[i] * dt

        # 2. Full-step position update: x(t+dt) = x(t) + v(t+dt/2) dt
        for body in bodies:
            body.x += body.v * dt

        # 3. Compute new accelerations a(t+dt) at new positions
        a_new = eih_1pn_accel(bodies, c_light, G_newton)

        # 4. Complete velocity update: v(t+dt) = v(t+dt/2) + (1/2) a(t+dt) dt
        for i, body in enumerate(bodies):
            body.v += 0.5 * a_new[i] * dt

        # Update time and accelerations
        t += dt
        a_current = a_new

        # Save trajectory data
        if step % save_every == 0:
            trajectory['t'][save_idx] = t
            for i, body in enumerate(bodies):
                trajectory['x'][save_idx, i] = body.x.copy()
                trajectory['v'][save_idx, i] = body.v.copy()
                trajectory['a'][save_idx, i] = a_current[i].copy()
            save_idx += 1

    # Trim trajectory arrays if final step wasn't saved
    if save_idx < n_saved:
        for key in trajectory:
            trajectory[key] = trajectory[key][:save_idx]

    return trajectory


def compute_orbit_elements(
    x: np.ndarray,
    v: np.ndarray,
    M_central: float,
    G_newton: float
) -> Dict[str, float]:
    """Compute osculating orbital elements for a test particle around a central mass.

    Given position and velocity of a test particle orbiting a dominant central
    mass M_central, compute the instantaneous (osculating) Keplerian elements.

    Parameters
    ----------
    x : np.ndarray
        Position vector of test particle, shape (3,).
    v : np.ndarray
        Velocity vector of test particle, shape (3,).
    M_central : float
        Mass of central body (assumed >> test particle mass).
    G_newton : float
        Gravitational constant.

    Returns
    -------
    elements : Dict[str, float]
        Dictionary containing:
        - 'a': semi-major axis
        - 'e': eccentricity
        - 'i': inclination [radians]
        - 'Omega': longitude of ascending node [radians]
        - 'omega': argument of periapsis [radians]
        - 'f': true anomaly [radians]
        - 'E': specific orbital energy
        - 'h': specific angular momentum magnitude

    Notes
    -----
    **Assumptions**:
    - Central body is much more massive than test particle
    - Two-body approximation is valid
    - Osculating elements: instantaneous Keplerian ellipse tangent to true orbit

    **Precession**:
    In GR 1PN, the orbital plane and perihelion precess. Osculating elements
    capture these slow variations: omega increases by ~43 arcsec/century for Mercury.

    **Standard Keplerian formulas**:

        E = v²/2 - G M / r          (specific energy)
        h = x × v                    (specific angular momentum vector)
        e_vec = (v × h) / (G M) - x/r   (eccentricity vector, points to periapsis)
        e = |e_vec|                  (eccentricity)
        a = -G M / (2 E)             (semi-major axis)
        i = arccos(h_z / |h|)        (inclination)

    Examples
    --------
    >>> # Mercury at perihelion (simplified)
    >>> x = np.array([4.60e10, 0.0, 0.0])  # m
    >>> v = np.array([0.0, 5.87e4, 0.0])   # m/s
    >>> M_sun = 1.989e30  # kg
    >>> G_si = 6.6743e-11
    >>> elems = compute_orbit_elements(x, v, M_sun, G_si)
    >>> print(f"a = {elems['a']/1e9:.2f} Gm")
    >>> print(f"e = {elems['e']:.4f}")
    """
    # Position and velocity magnitudes
    r = np.linalg.norm(x)
    v_mag = np.linalg.norm(v)

    # Specific orbital energy: E = v²/2 - GM/r
    E = 0.5 * v_mag**2 - G_newton * M_central / r

    # Specific angular momentum: h = r × v
    h_vec = np.cross(x, v)
    h_mag = np.linalg.norm(h_vec)

    # Semi-major axis: a = -GM / (2E)
    if E >= 0:
        a = np.inf  # Unbound orbit
    else:
        a = -G_newton * M_central / (2.0 * E)

    # Eccentricity vector: e_vec = (v × h)/(GM) - r_hat
    # where r_hat = x / |x|
    e_vec = np.cross(v, h_vec) / (G_newton * M_central) - x / r
    e = np.linalg.norm(e_vec)

    # Inclination: i = arccos(h_z / |h|)
    if h_mag > 0:
        i = np.arccos(np.clip(h_vec[2] / h_mag, -1.0, 1.0))
    else:
        i = 0.0  # Undefined for zero angular momentum

    # Node vector: n = z_hat × h  (points to ascending node)
    z_hat = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(z_hat, h_vec)
    n_mag = np.linalg.norm(n_vec)

    # Longitude of ascending node: Omega
    if n_mag > 0:
        Omega = np.arctan2(n_vec[1], n_vec[0])
        if Omega < 0:
            Omega += 2 * np.pi
    else:
        Omega = 0.0  # Undefined for equatorial orbits

    # Argument of periapsis: omega
    if n_mag > 0 and e > 0:
        omega_cos = np.dot(n_vec, e_vec) / (n_mag * e)
        omega = np.arccos(np.clip(omega_cos, -1.0, 1.0))
        if e_vec[2] < 0:  # Check z-component of e_vec
            omega = 2 * np.pi - omega
    else:
        omega = 0.0  # Undefined for circular or equatorial orbits

    # True anomaly: f
    if e > 0:
        f_cos = np.dot(e_vec, x) / (e * r)
        f = np.arccos(np.clip(f_cos, -1.0, 1.0))
        if np.dot(x, v) < 0:  # Check if moving towards periapsis
            f = 2 * np.pi - f
    else:
        f = 0.0  # Undefined for circular orbits

    return {
        'a': a,
        'e': e,
        'i': i,
        'Omega': Omega,
        'omega': omega,
        'f': f,
        'E': E,
        'h': h_mag,
    }


def find_periastron_passages(
    t: np.ndarray,
    x: np.ndarray
) -> np.ndarray:
    """Find times of periastron passage (closest approach to origin).

    Parameters
    ----------
    t : np.ndarray
        Time array, shape (n_steps,).
    x : np.ndarray
        Position trajectory, shape (n_steps, 3).

    Returns
    -------
    t_peri : np.ndarray
        Array of times at periastron passages.

    Notes
    -----
    Periastron is detected by finding local minima in r(t) = |x(t)|.
    Uses simple finite differencing: r'(t) changes sign from - to +.

    For eccentric orbits, periastron is well-defined. For circular orbits,
    this will find spurious "passages" due to noise.
    """
    r = np.linalg.norm(x, axis=1)

    # Find local minima: r[i-1] > r[i] < r[i+1]
    is_minimum = (r[1:-1] < r[:-2]) & (r[1:-1] < r[2:])
    peri_indices = np.where(is_minimum)[0] + 1  # Offset by 1 due to slicing

    return t[peri_indices]


def compute_precession_rate(
    trajectory: Dict[str, np.ndarray],
    body_idx: int,
    M_central: float,
    G_newton: float
) -> Dict[str, float]:
    """Compute perihelion precession rate from a GR trajectory.

    Parameters
    ----------
    trajectory : Dict[str, np.ndarray]
        Trajectory data from integrate_gr1pn_orbit.
    body_idx : int
        Index of the orbiting body (e.g., 1 for Mercury if Sun is index 0).
    M_central : float
        Mass of the central body.
    G_newton : float
        Gravitational constant.

    Returns
    -------
    precession : Dict[str, float]
        Dictionary containing:
        - 'omega_mean': mean argument of periapsis [radians]
        - 'domega_dt': precession rate [radians per unit time]
        - 'domega_per_orbit': precession per orbit [radians]
        - 'n_orbits': number of complete orbits analyzed

    Notes
    -----
    **Method**:
    1. Compute osculating elements at each saved time step
    2. Track argument of periapsis omega(t)
    3. Unwrap 2π discontinuities
    4. Fit linear trend: omega(t) ≈ omega_0 + (dω/dt) t
    5. Compute precession per orbit from orbital period

    **GR prediction (two-body, M_1 >> M_2)**:

        Δω_GR = (6π G M_1) / (a c² (1 - e²))    [radians per orbit]

    For Mercury: Δω_GR ≈ 5.0×10⁻⁷ rad/orbit ≈ 43 arcsec/century.

    Examples
    --------
    >>> traj = integrate_gr1pn_orbit(bodies, c_light, G, dt, n_steps)
    >>> prec = compute_precession_rate(traj, body_idx=1, M_central=M_sun, G_newton=G)
    >>> print(f"Precession: {prec['domega_per_orbit']*206265:.2f} arcsec/orbit")
    """
    t = trajectory['t']
    x = trajectory['x'][:, body_idx, :]
    v = trajectory['v'][:, body_idx, :]

    # Compute osculating omega at each time step
    omega_list = []
    for i in range(len(t)):
        elems = compute_orbit_elements(x[i], v[i], M_central, G_newton)
        omega_list.append(elems['omega'])

    omega = np.array(omega_list)

    # Unwrap 2π discontinuities
    omega_unwrapped = np.unwrap(omega)

    # Fit linear trend: omega(t) = omega_0 + (dω/dt) * t
    coeffs = np.polyfit(t, omega_unwrapped, deg=1)
    domega_dt = coeffs[0]  # radians per unit time
    omega_mean = coeffs[1]

    # Estimate orbital period from mean semi-major axis
    a_list = []
    for i in range(len(t)):
        elems = compute_orbit_elements(x[i], v[i], M_central, G_newton)
        if np.isfinite(elems['a']):
            a_list.append(elems['a'])

    a_mean = np.mean(a_list) if a_list else 0.0

    if a_mean > 0:
        # Kepler's third law: P = 2π sqrt(a³ / (G M))
        P_orb = 2 * np.pi * np.sqrt(a_mean**3 / (G_newton * M_central))
        n_orbits = (t[-1] - t[0]) / P_orb
        domega_per_orbit = domega_dt * P_orb
    else:
        P_orb = np.nan
        n_orbits = 0.0
        domega_per_orbit = np.nan

    return {
        'omega_mean': omega_mean,
        'domega_dt': domega_dt,
        'domega_per_orbit': domega_per_orbit,
        'n_orbits': n_orbits,
        'P_orb': P_orb,
    }


def compare_precession(
    slab_trajectory: Dict[str, np.ndarray],
    gr_trajectory: Dict[str, np.ndarray],
    slab_medium,
    body_idx: int,
    M_central_slab: float,
    M_central_gr: float,
    G_newton: float
) -> Dict[str, float]:
    """Compare perihelion precession between slab and GR simulations.

    This is the key validation function: it measures precession rates from
    both simulations and computes the ratio, assessing whether the slab can
    reproduce GR-like 1PN effects.

    Parameters
    ----------
    slab_trajectory : Dict[str, np.ndarray]
        Trajectory from slab simulation (from dynamics.integrate_orbit).
    gr_trajectory : Dict[str, np.ndarray]
        Trajectory from GR simulation (from integrate_gr1pn_orbit).
    slab_medium : Medium
        Slab medium object (to extract K for orbital period calculation).
    body_idx : int
        Index of orbiting body (e.g., 1 for planet if central body is 0).
    M_central_slab : float
        Mass of central body in slab units.
    M_central_gr : float
        Mass of central body in GR/SI units.
    G_newton : float
        Gravitational constant for GR.

    Returns
    -------
    comparison : Dict[str, float]
        Dictionary containing:
        - 'slab_domega_per_orbit': slab precession [rad/orbit]
        - 'gr_domega_per_orbit': GR precession [rad/orbit]
        - 'ratio': slab / GR precession ratio (should be ~1 if matching)
        - 'slab_n_orbits': number of orbits in slab simulation
        - 'gr_n_orbits': number of orbits in GR simulation

    Notes
    -----
    **Interpretation**:
    - ratio ≈ 1.0: Slab matches GR precession (validates 1PN emergence)
    - ratio < 1.0: Slab under-precesses (cs too large, or compressible effects weak)
    - ratio > 1.0: Slab over-precesses (cs too small, or other corrections dominate)

    **Parameter mapping**:
    For a given slab medium (ρ₀, cs, β₀), the effective "gravitational constant" is:
        K_slab = ρ₀ / (4π β₀²)

    To match GR with (G, c), choose:
        K_slab ↔ G
        cs_slab ↔ c

    Then compare precession rates. The slab should reproduce GR precession
    within model uncertainties (~few percent) if cs is chosen correctly.

    Examples
    --------
    >>> # After running both simulations with same initial conditions
    >>> comparison = compare_precession(
    ...     slab_traj, gr_traj, medium,
    ...     body_idx=1, M_central_slab=1.0, M_central_gr=1.989e30, G_newton=6.67e-11
    ... )
    >>> print(f"Slab precession: {comparison['slab_domega_per_orbit']*206265:.2f} arcsec/orbit")
    >>> print(f"GR precession: {comparison['gr_domega_per_orbit']*206265:.2f} arcsec/orbit")
    >>> print(f"Ratio (slab/GR): {comparison['ratio']:.4f}")
    """
    # For slab: need to compute precession using K instead of G
    # The slab trajectory has same structure as GR trajectory
    t_slab = slab_trajectory['t']
    x_slab = slab_trajectory['x'][:, body_idx, :]
    v_slab = slab_trajectory['v'][:, body_idx, :]

    # Use K from slab medium as "effective G"
    K_slab = slab_medium.K

    # Compute osculating omega for slab
    omega_slab_list = []
    for i in range(len(t_slab)):
        elems = compute_orbit_elements(x_slab[i], v_slab[i], M_central_slab, K_slab)
        omega_slab_list.append(elems['omega'])

    omega_slab = np.array(omega_slab_list)
    omega_slab_unwrapped = np.unwrap(omega_slab)

    # Fit linear trend for slab
    coeffs_slab = np.polyfit(t_slab, omega_slab_unwrapped, deg=1)
    domega_dt_slab = coeffs_slab[0]

    # Estimate orbital period for slab
    a_list_slab = []
    for i in range(len(t_slab)):
        elems = compute_orbit_elements(x_slab[i], v_slab[i], M_central_slab, K_slab)
        if np.isfinite(elems['a']):
            a_list_slab.append(elems['a'])

    a_mean_slab = np.mean(a_list_slab) if a_list_slab else 0.0

    if a_mean_slab > 0:
        P_orb_slab = 2 * np.pi * np.sqrt(a_mean_slab**3 / (K_slab * M_central_slab))
        n_orbits_slab = (t_slab[-1] - t_slab[0]) / P_orb_slab
        domega_per_orbit_slab = domega_dt_slab * P_orb_slab
    else:
        n_orbits_slab = 0.0
        domega_per_orbit_slab = np.nan

    # Compute GR precession
    prec_gr = compute_precession_rate(gr_trajectory, body_idx, M_central_gr, G_newton)
    domega_per_orbit_gr = prec_gr['domega_per_orbit']
    n_orbits_gr = prec_gr['n_orbits']

    # Compute ratio
    if np.isfinite(domega_per_orbit_slab) and np.isfinite(domega_per_orbit_gr) and domega_per_orbit_gr != 0:
        ratio = domega_per_orbit_slab / domega_per_orbit_gr
    else:
        ratio = np.nan

    return {
        'slab_domega_per_orbit': domega_per_orbit_slab,
        'gr_domega_per_orbit': domega_per_orbit_gr,
        'ratio': ratio,
        'slab_n_orbits': n_orbits_slab,
        'gr_n_orbits': n_orbits_gr,
    }
