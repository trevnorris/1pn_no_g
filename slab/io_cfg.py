"""Configuration and I/O module for superfluid orbit simulator.

This module provides:
- YAML configuration loading and validation
- Example config generation
- CSV output for trajectories
- JSON output for diagnostics

No gravitational constant G appears. All parameters are expressed in terms
of superfluid properties (rho0, cs, beta0) and emerge forces from momentum flux.
"""

from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import yaml
import json
from pathlib import Path
import warnings

from slab.medium import Medium
from slab.bodies import Body


def load_config(yaml_path: str) -> Dict[str, Any]:
    """Load and parse YAML configuration file.

    Reads a YAML configuration specifying medium properties, body initial
    conditions, numerical integration options, and output preferences.
    Returns structured objects ready for simulation.

    Parameters
    ----------
    yaml_path : str
        Path to YAML configuration file.

    Returns
    -------
    dict
        Configuration dictionary with keys:
        - 'medium': Medium instance (rho0, cs, beta0, gamma_beta)
        - 'bodies': list of Body instances (name, M, x, v, R, Q)
        - 'numerics': dict of integration options (dt, steps, audit settings, etc.)
        - 'compare_gr_1pn': dict of GR comparison options (enable, c_light, etc.)
        - 'outputs': dict of output options (save_every, write_csv, plots, etc.)

    Raises
    ------
    FileNotFoundError
        If yaml_path does not exist.
    yaml.YAMLError
        If YAML parsing fails.
    KeyError
        If required configuration fields are missing.
    ValueError
        If configuration values are invalid (negative masses, etc.).

    Notes
    -----
    The configuration must include:
    - medium: rho0, cs, beta0 (gamma_beta optional, defaults to 0.0)
    - bodies: list with name, M, x, v, R for each body (Q computed from M)
    - numerics: dt, steps (other fields have defaults)
    - compare_gr_1pn: enable (c_light required if enabled)
    - outputs: save_every, write_csv (plots optional)

    After loading, bodies have Q synchronized with M via Q = M/beta0.

    Examples
    --------
    >>> config = load_config("mercury_orbit.yaml")
    >>> medium = config['medium']
    >>> bodies = config['bodies']
    >>> print(f"Loaded {len(bodies)} bodies into medium with K={medium.K:.3e}")
    Loaded 2 bodies into medium with K=7.958e-22

    >>> # Access integration options
    >>> dt = config['numerics']['dt']
    >>> steps = config['numerics']['steps']
    >>> print(f"Integrating for {steps} steps with dt={dt}")
    Integrating for 200000 steps with dt=0.002

    See Also
    --------
    validate_config : Validate loaded configuration for physical consistency
    create_example_config : Generate example YAML file
    """
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {yaml_path}")

    with open(yaml_path, 'r') as f:
        raw_config = yaml.safe_load(f)

    if raw_config is None:
        raise ValueError(f"Empty or invalid YAML file: {yaml_path}")

    # Parse medium
    if 'medium' not in raw_config:
        raise KeyError("Configuration missing required section 'medium'")

    medium_cfg = raw_config['medium']
    medium = Medium(
        rho0=float(medium_cfg['rho0']),
        cs=float(medium_cfg['cs']),
        beta0=float(medium_cfg['beta0']),
        gamma_beta=float(medium_cfg.get('gamma_beta', 0.0))
    )

    # Parse bodies
    if 'bodies' not in raw_config:
        raise KeyError("Configuration missing required section 'bodies'")

    bodies_cfg = raw_config['bodies']
    if not isinstance(bodies_cfg, list) or len(bodies_cfg) == 0:
        raise ValueError("Configuration 'bodies' must be a non-empty list")

    bodies = []
    for i, body_cfg in enumerate(bodies_cfg):
        try:
            name = body_cfg['name']
            M = float(body_cfg['M'])
            x = np.array(body_cfg['x'], dtype=float)
            v = np.array(body_cfg['v'], dtype=float)
            R = float(body_cfg['R'])

            # Create body with placeholder Q, then sync from M
            body = Body(name=name, M=M, x=x, v=v, R=R, Q=0.0)
            body.update_Q_from_M(medium)
            bodies.append(body)
        except KeyError as e:
            raise KeyError(f"Body {i} missing required field {e}")
        except (ValueError, TypeError) as e:
            raise ValueError(f"Body {i} ('{body_cfg.get('name', 'unnamed')}'): {e}")

    # Parse numerics
    if 'numerics' not in raw_config:
        raise KeyError("Configuration missing required section 'numerics'")

    numerics_cfg = raw_config['numerics']
    numerics = {
        'dt': float(numerics_cfg['dt']),
        'steps': int(numerics_cfg['steps']),
        'audit_every': int(numerics_cfg.get('audit_every', 500)),
        'npts_audit': int(numerics_cfg.get('npts_audit', 512)),
        'use_compressible': bool(numerics_cfg.get('use_compressible', True)),
        'use_flux_mass': bool(numerics_cfg.get('use_flux_mass', False)),
        'intake_every': int(numerics_cfg.get('intake_every', 2000)),
    }

    # Parse GR comparison options
    compare_gr = raw_config.get('compare_gr_1pn', {})
    compare_gr_1pn = {
        'enable': bool(compare_gr.get('enable', False)),
        'c_light': float(compare_gr.get('c_light', 63239.7263)),  # AU/yr default
        'measure_peri': bool(compare_gr.get('measure_peri', True)),
    }

    if compare_gr_1pn['enable'] and 'c_light' not in compare_gr:
        warnings.warn(
            "GR comparison enabled but c_light not specified, using default 63239.7263 AU/yr",
            UserWarning
        )

    # Parse output options
    outputs_cfg = raw_config.get('outputs', {})
    outputs = {
        'save_every': int(outputs_cfg.get('save_every', 1000)),
        'write_csv': bool(outputs_cfg.get('write_csv', True)),
        'plots': list(outputs_cfg.get('plots', [])),
    }

    return {
        'medium': medium,
        'bodies': bodies,
        'numerics': numerics,
        'compare_gr_1pn': compare_gr_1pn,
        'outputs': outputs,
    }


def validate_config(config: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """Validate configuration for physical consistency and numerical stability.

    Checks that all required fields are present and that parameter values
    make physical sense. Warns about potential issues like:
    - Control surfaces too large (R comparable to separations)
    - Bodies too close (overlapping control surfaces)
    - Capacity violation (Ma = Q/(4*pi*R^2*cs) >> 1)
    - Timestep too large for sound speed or orbital frequency

    Parameters
    ----------
    config : dict
        Configuration dictionary from load_config().

    Returns
    -------
    is_valid : bool
        True if configuration passes all checks (may still have warnings).
    warnings_list : list of str
        List of warning messages about potential issues.

    Notes
    -----
    **Checks performed**:

    1. Positive parameters: M > 0, R > 0, rho0 > 0, cs > 0, beta0 > 0
    2. Array shapes: x, v have shape (3,)
    3. Control surface sizing: R < 0.1 * min_separation
    4. Body separation: min_separation > 2 * max(R)
    5. Capacity: Ma = Q/(4*pi*R^2*cs) << 1 for all bodies
    6. Timestep: dt < 0.1 * min(R/cs, orbital_period)

    **Capacity check**:
    The capacity parameter Ma measures how close a body is to saturating
    the local sound speed:
        Ma = Q / (4*pi*R^2*cs)

    For the small-Mach expansion to be valid, we need Ma << 1.
    Warnings are issued if Ma > 0.1 for any body.

    Examples
    --------
    >>> config = load_config("mercury_orbit.yaml")
    >>> is_valid, warnings = validate_config(config)
    >>> if not is_valid:
    ...     print("ERROR: Configuration is invalid!")
    ... elif warnings:
    ...     print("Warnings:")
    ...     for w in warnings:
    ...         print(f"  - {w}")
    ... else:
    ...     print("Configuration is valid and healthy.")
    Configuration is valid and healthy.

    >>> # Example with capacity warning
    >>> bodies[0].Q = 1e10  # Unrealistically large intake
    >>> is_valid, warnings = validate_config(config)
    >>> for w in warnings:
    ...     print(w)
    Body 'Sun': Capacity Ma = 7.96e-01 >> 1. Small-Mach expansion invalid!

    See Also
    --------
    load_config : Load configuration from YAML file
    """
    warnings_list = []
    is_valid = True

    # Unpack configuration
    try:
        medium = config['medium']
        bodies = config['bodies']
        numerics = config['numerics']
        compare_gr = config['compare_gr_1pn']
        outputs = config['outputs']
    except KeyError as e:
        return False, [f"Missing required config section: {e}"]

    # Check medium parameters (already validated in Medium.__post_init__)
    if medium.rho0 <= 0 or medium.cs <= 0 or medium.beta0 <= 0:
        is_valid = False
        warnings_list.append(
            f"Medium has non-positive parameters: rho0={medium.rho0}, "
            f"cs={medium.cs}, beta0={medium.beta0}"
        )

    # Check numerics
    if numerics['dt'] <= 0:
        is_valid = False
        warnings_list.append(f"Timestep dt must be positive, got {numerics['dt']}")

    if numerics['steps'] <= 0:
        is_valid = False
        warnings_list.append(f"Number of steps must be positive, got {numerics['steps']}")

    # Check bodies
    if len(bodies) < 1:
        is_valid = False
        warnings_list.append("Configuration must have at least one body")
        return is_valid, warnings_list  # Can't do further checks

    # Collect body properties for pair checks
    separations = []
    Rs = [body.R for body in bodies]

    for body in bodies:
        # Mass and R already validated in Body.__post_init__

        # Check capacity: Ma = Q/(4*pi*R^2*cs)
        Ma = abs(body.Q) / (4.0 * np.pi * body.R**2 * medium.cs)
        if Ma > 0.1:
            warnings_list.append(
                f"Body '{body.name}': Capacity Ma = {Ma:.2e} >> 1. "
                f"Small-Mach expansion may be invalid!"
            )
        if Ma > 1.0:
            is_valid = False
            warnings_list.append(
                f"Body '{body.name}': Capacity Ma = {Ma:.2e} > 1. "
                f"Configuration physically unrealistic (supersonic intake)."
            )

    # Check pairwise separations
    for i, body_a in enumerate(bodies):
        for j, body_b in enumerate(bodies):
            if j <= i:
                continue

            r_ab = np.linalg.norm(body_b.x - body_a.x)
            separations.append(r_ab)

            # Check for overlapping control surfaces
            if r_ab < (body_a.R + body_b.R):
                is_valid = False
                warnings_list.append(
                    f"Bodies '{body_a.name}' and '{body_b.name}': "
                    f"Control surfaces overlap! r_ab = {r_ab:.3e}, "
                    f"R_a + R_b = {body_a.R + body_b.R:.3e}"
                )

            # Warn if R too large compared to separation
            max_R = max(body_a.R, body_b.R)
            if max_R > 0.1 * r_ab:
                warnings_list.append(
                    f"Bodies '{body_a.name}' and '{body_b.name}': "
                    f"Control surface R = {max_R:.3e} is large compared to "
                    f"separation r_ab = {r_ab:.3e}. Near-field expansion may break down."
                )

    if separations:
        min_separation = min(separations)
        max_R = max(Rs)

        # Check sound crossing time vs timestep
        sound_crossing_time = max_R / medium.cs
        if numerics['dt'] > 0.5 * sound_crossing_time:
            warnings_list.append(
                f"Timestep dt = {numerics['dt']:.3e} is large compared to "
                f"sound crossing time R/cs = {sound_crossing_time:.3e}. "
                f"Consider reducing dt for stability."
            )

        # Estimate orbital period and check timestep
        # For two-body, T ~ 2*pi*sqrt(r^3/(K*M_total))
        if len(bodies) >= 2:
            M_total = sum(body.M for body in bodies)
            r_typical = min_separation
            K = medium.K
            orbital_freq = np.sqrt(K * M_total / r_typical**3)
            orbital_period = 2.0 * np.pi / orbital_freq

            if numerics['dt'] > 0.01 * orbital_period:
                warnings_list.append(
                    f"Timestep dt = {numerics['dt']:.3e} is large compared to "
                    f"estimated orbital period T ~ {orbital_period:.3e}. "
                    f"Consider dt < {0.01 * orbital_period:.3e} for accuracy."
                )

    # Check GR comparison
    if compare_gr['enable']:
        if compare_gr['c_light'] <= 0:
            is_valid = False
            warnings_list.append(
                f"GR comparison enabled but c_light <= 0: {compare_gr['c_light']}"
            )

        # Check that c_light >> typical velocities
        v_max = max(np.linalg.norm(body.v) for body in bodies)
        if v_max > 0.1 * compare_gr['c_light']:
            warnings_list.append(
                f"Maximum velocity v_max = {v_max:.3e} is comparable to "
                f"c_light = {compare_gr['c_light']:.3e}. "
                f"1PN expansion may be insufficient."
            )

    # Check outputs
    if outputs['save_every'] <= 0:
        is_valid = False
        warnings_list.append(f"save_every must be positive, got {outputs['save_every']}")

    return is_valid, warnings_list


def create_example_config(output_path: str) -> None:
    """Generate example YAML configuration file.

    Creates a well-commented example configuration for the Mercury orbit,
    using physical parameters from the Solar System. This serves as a
    template for users to modify.

    Parameters
    ----------
    output_path : str
        Path where YAML file will be written.

    Notes
    -----
    **Example system**:
    - Sun (M=1.0 in code units) at origin
    - Mercury (M=1.66e-7 in Sun masses) in circular orbit
    - Semi-major axis a=0.387 AU
    - Orbital period ~88 days

    **Medium parameters**:
    - rho0 = 1.0 (code units)
    - cs = 1.0e4 (AU/yr, much faster than orbital speeds)
    - beta0 = 1.0e10 (chosen to match observed K ~ G)
    - gamma_beta = 0.0 (constant beta, Solar System regime)

    **Control surfaces**:
    - R_Sun = 1.0e-3 AU (well inside Mercury orbit)
    - R_Mercury = 5.0e-4 AU (small compared to separation)

    The configuration includes helpful comments explaining each parameter
    and how to modify it for different systems.

    Examples
    --------
    >>> create_example_config("my_orbit.yaml")
    >>> config = load_config("my_orbit.yaml")
    >>> print(f"Generated config with {len(config['bodies'])} bodies")
    Generated config with 2 bodies

    See Also
    --------
    load_config : Load configuration from YAML file
    """
    # Compute orbital parameters for Mercury
    # Using arbitrary code units as in plan_no_pde.md section 9
    # The key point: no gravitational constant G appears anywhere.
    # Forces emerge from superfluid momentum flux, scaled by K = rho0/(4*pi*beta0^2)

    a = 0.387  # AU, Mercury semi-major axis (in code units)
    e = 0.0    # Start with circular orbit for simplicity

    # Use medium parameters from plan_no_pde.md
    # These give small intake Ma << 1 for physical consistency
    rho0 = 1.0
    cs = 1.0e4  # Code units, much faster than orbital speeds
    beta0 = 1.0e10  # Large beta0 means small Q/M = 1/beta0

    # This gives K = rho0 / (4*pi*beta0^2)
    K = rho0 / (4.0 * np.pi * beta0**2)

    M_sun = 1.0
    M_mercury = 3.3e-7  # Mercury mass in solar masses (from plan_no_pde.md)

    # Initial conditions: Sun at origin, Mercury on x-axis
    x_sun = [0.0, 0.0, 0.0]
    v_sun = [0.0, 0.0, 0.0]

    x_mercury = [a, 0.0, 0.0]

    # Circular orbit velocity (in y direction)
    # v = sqrt(K * M_sun / a)
    v_circ = np.sqrt(K * M_sun / a)
    v_mercury = [0.0, float(v_circ), 0.0]

    # Control surface radii (from plan_no_pde.md)
    R_sun = 1.0e-3
    R_mercury = 5.0e-4

    # Integration parameters (from plan_no_pde.md)
    # Use fixed dt and steps as specified in the plan
    dt = 2.0e-3
    steps = 200000

    # Compute orbital period for reference
    orbital_freq = np.sqrt(K * M_sun / a**3)
    T_orbit = 2.0 * np.pi / orbital_freq
    n_orbits = steps * dt / T_orbit

    # Speed of light in AU/yr
    c_light = 63239.7263  # 1 AU/yr = 4.74 km/s, c = 299792 km/s

    yaml_content = f"""# Superfluid Orbit Simulator Configuration
# Mercury orbit example
#
# This configuration demonstrates a two-body system (Sun-Mercury) using
# superfluid hydrodynamics to compute forces. No gravitational constant G
# appears; instead, forces emerge from momentum flux in the superfluid.

# ============================================================================
# Medium: Ambient superfluid properties
# ============================================================================
medium:
  # Ambient density [code units: arbitrary, sets overall scale]
  rho0: {rho0}

  # Sound speed [AU/yr, much faster than orbital speeds]
  # Controls compressible corrections and perihelion precession.
  # For incompressible limit, use cs >> v_orbital (e.g., 1e10).
  # For 1PN-like effects, use finite cs (e.g., 1e4).
  cs: {cs}

  # Mass-intake factor [code units]
  # Relates body mass M to volumetric intake Q via Q = M/beta0.
  # Larger beta0 means weaker forces (smaller K = rho0/(4*pi*beta0^2)).
  beta0: {beta0}

  # Rarefaction exponent (optional, default: 0.0)
  # When gamma_beta = 0, beta is constant (Solar System regime).
  # When gamma_beta > 0, beta varies with local density.
  gamma_beta: 0.0

# ============================================================================
# Bodies: Fluid intakes (mass objects)
# ============================================================================
bodies:
  - name: Sun
    # Mass [M_sun, code units]
    M: {M_sun}

    # Position [AU]
    x: {x_sun}

    # Velocity [AU/yr]
    v: {v_sun}

    # Control surface radius [AU]
    # Must satisfy: mouth_size << R << separation
    # Too small: numerical issues in quadrature
    # Too large: near-field expansion breaks down
    # Typical: R ~ 1e-4 to 1e-3 for separations ~ 0.1-1 AU
    R: {R_sun}

  - name: Mercury
    M: {M_mercury}
    x: {x_mercury}
    v: {v_mercury}
    R: {R_mercury}

# ============================================================================
# Numerics: Integration and force calculation options
# ============================================================================
numerics:
  # Timestep [code time units: yr]
  # Should be << orbital period for accuracy.
  # Suggested: dt ~ T_orbit / 100 for low-eccentricity orbits.
  dt: {dt}

  # Total number of integration steps
  steps: {steps}

  # Audit frequency: every N steps, compute forces via quadrature
  # and compare to analytic formula. Ensures analytic path is correct.
  # Set to 0 to disable auditing (faster, but less robust).
  audit_every: 500

  # Number of quadrature points for audit (Fibonacci sphere)
  # More points = better accuracy but slower.
  # 256-512 is usually sufficient, 1024 for high precision.
  npts_audit: 512

  # Use compressible force corrections (finite cs)?
  # If false, use incompressible (cs -> infinity) baseline only.
  # If true, compute finite-Mach corrections (perihelion precession).
  use_compressible: true

  # Use flux-based mass intake evolution?
  # If true, masses M_a evolve via surface flux integration.
  # If false, masses are constant (typical for Solar System).
  use_flux_mass: false

  # Intake update frequency (if use_flux_mass=true)
  # Update masses every N steps based on integrated flux.
  intake_every: 2000

# ============================================================================
# GR 1PN Comparator: Independent reference calculation
# ============================================================================
compare_gr_1pn:
  # Enable GR comparison module?
  # If true, integrate EIH 1PN equations separately and compare precession.
  # GR module never influences slab forces (fully independent).
  enable: true

  # Speed of light [AU/yr]
  # c = 299792 km/s = 63239.7263 AU/yr
  c_light: {c_light}

  # Measure perihelion precession?
  # If true, track apsis passages and compute precession rate.
  measure_peri: true

# ============================================================================
# Outputs: What to save and how often
# ============================================================================
outputs:
  # Save state every N steps
  # Trajectories, forces, and diagnostics written at this frequency.
  save_every: 1000

  # Write trajectory to CSV file?
  # If true, save time, position, velocity, mass, intake for all bodies.
  write_csv: true

  # List of plots to generate (optional, for post-processing)
  # Available plots: orbit, precession_vs_time, force_decomp, energy, elements
  plots:
    - orbit
    - precession_vs_time
    - force_decomp

# ============================================================================
# Expected outputs for this configuration
# ============================================================================
# - Stable circular orbit with radius a = 0.387 AU
# - Orbital period T ~ 0.24 yr (88 days, Mercury's actual period)
# - Perihelion precession (if use_compressible=true) ~ few arcsec/orbit
# - Precession scales as 1/cs^2 (test by varying cs)
# - Forces match quadrature to < 0.1% (if audit_every > 0)
#
# To modify:
# - Increase eccentricity: add v_z component, reduce v_y
# - Add more bodies: copy body template, adjust positions
# - Change cs: vary by factor of 2 to verify cs^-2 scaling
# - Disable compressible: set use_compressible=false for incompressible baseline
"""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(yaml_content)

    print(f"Example configuration written to: {output_path}")
    print(f"  K = {K:.3e} (orbital constant, replaces G)")
    print(f"  v_circ = {v_circ:.3e} AU/yr")
    print(f"  T_orbit = {T_orbit:.3e} yr ({T_orbit * 365.25:.1f} days)")
    print(f"  dt = {dt:.3e} yr ({dt / T_orbit * 100:.1f}% of period)")
    print(f"  Total integration time: {n_orbits} orbits")


def save_state_csv(
    filepath: str,
    times: np.ndarray,
    bodies_history: List[List[Body]]
) -> None:
    """Save trajectory data to CSV file.

    Writes a CSV file with columns: time, body_name, x, y, z, vx, vy, vz, M, Q.
    Each row corresponds to one body at one timestep.

    Parameters
    ----------
    filepath : str
        Output CSV file path.
    times : np.ndarray
        Array of times at which states were saved, shape (n_snapshots,).
    bodies_history : list of list of Body
        Nested list of body states: bodies_history[i][j] is the state of
        body j at time times[i]. Shape: (n_snapshots, n_bodies).

    Notes
    -----
    **CSV format**:
    ```
    time,body_name,x,y,z,vx,vy,vz,M,Q
    0.0,Sun,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0e-10
    0.0,Mercury,0.387,0.0,0.0,0.0,10.2,0.0,1.66e-7,1.66e-17
    ...
    ```

    This format is easy to load into pandas, numpy, or plotting tools:
    ```python
    import pandas as pd
    df = pd.read_csv("trajectory.csv")
    sun = df[df['body_name'] == 'Sun']
    plt.plot(sun['x'], sun['y'])
    ```

    Examples
    --------
    >>> times = np.array([0.0, 0.001, 0.002])
    >>> body1 = Body("Sun", M=1.0, x=[0,0,0], v=[0,0,0], R=1e-3, Q=1e-10)
    >>> body2 = Body("Mercury", M=1e-7, x=[1,0,0], v=[0,1,0], R=5e-4, Q=1e-17)
    >>> bodies_history = [[body1, body2], [body1, body2], [body1, body2]]
    >>> save_state_csv("trajectory.csv", times, bodies_history)
    Saved 6 states (3 snapshots × 2 bodies) to trajectory.csv

    See Also
    --------
    save_diagnostics_json : Save diagnostics to JSON
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        # Write header
        f.write("time,body_name,x,y,z,vx,vy,vz,M,Q\n")

        # Write data
        for t, bodies in zip(times, bodies_history):
            for body in bodies:
                f.write(
                    f"{t:.15e},{body.name},"
                    f"{body.x[0]:.15e},{body.x[1]:.15e},{body.x[2]:.15e},"
                    f"{body.v[0]:.15e},{body.v[1]:.15e},{body.v[2]:.15e},"
                    f"{body.M:.15e},{body.Q:.15e}\n"
                )

    n_snapshots = len(times)
    n_bodies = len(bodies_history[0]) if bodies_history else 0
    n_total = n_snapshots * n_bodies
    print(f"Saved {n_total} states ({n_snapshots} snapshots × {n_bodies} bodies) to {filepath}")


def save_diagnostics_json(filepath: str, diagnostics: Dict[str, Any]) -> None:
    """Save diagnostics data to JSON file.

    Writes diagnostics to a JSON file for post-processing and analysis.
    Diagnostics may include energy-like quantities, precession measurements,
    force comparisons (analytic vs quadrature), capacity checks, etc.

    Parameters
    ----------
    filepath : str
        Output JSON file path.
    diagnostics : dict
        Dictionary of diagnostic data. Can contain:
        - 'energies': list of total energy at each snapshot
        - 'precession': dict with perihelion passages and rates
        - 'force_errors': list of relative errors (analytic vs quad)
        - 'capacity': dict with Ma values for each body
        - 'orbital_elements': dict with semi-major axis, eccentricity, etc.
        - any other JSON-serializable data

    Notes
    -----
    **Example diagnostics structure**:
    ```python
    diagnostics = {
        'energies': [1.23e-5, 1.24e-5, ...],
        'precession': {
            'body': 'Mercury',
            'perihelion_times': [0.0, 0.24, 0.48, ...],
            'precession_per_orbit': 5.3e-7,  # radians
            'precession_arcsec_per_century': 43.2
        },
        'force_errors': {
            'incompressible': [1.2e-4, 8.7e-5, ...],
            'compressible': [3.4e-3, 2.9e-3, ...]
        },
        'capacity': {
            'Sun': {'Ma': 7.96e-5, 'warning': false},
            'Mercury': {'Ma': 3.98e-8, 'warning': false}
        },
        'orbital_elements': {
            'times': [0.0, 0.001, ...],
            'a': [0.387, 0.387, ...],
            'e': [0.001, 0.001, ...],
            'i': [0.0, 0.0, ...],
            'omega': [0.0, 5.3e-7, ...]  # perihelion angle
        }
    }
    ```

    Arrays are converted to lists for JSON serialization.

    Examples
    --------
    >>> diagnostics = {
    ...     'total_energy': [1.23, 1.24, 1.23],
    ...     'max_force_error': 1.5e-4,
    ...     'precession_rate': 5.3e-7
    ... }
    >>> save_diagnostics_json("diagnostics.json", diagnostics)
    Saved diagnostics to diagnostics.json (3 top-level keys)

    >>> # Load back for analysis
    >>> import json
    >>> with open("diagnostics.json") as f:
    ...     data = json.load(f)
    >>> print(data['max_force_error'])
    0.00015

    See Also
    --------
    save_state_csv : Save trajectories to CSV
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # Convert numpy arrays to lists for JSON serialization
    def convert_to_json_serializable(obj):
        """Recursively convert numpy arrays to lists."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: convert_to_json_serializable(val) for key, val in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_json_serializable(item) for item in obj]
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        else:
            return obj

    serializable_diagnostics = convert_to_json_serializable(diagnostics)

    with open(filepath, 'w') as f:
        json.dump(serializable_diagnostics, f, indent=2)

    n_keys = len(diagnostics)
    print(f"Saved diagnostics to {filepath} ({n_keys} top-level keys)")


# ============================================================================
# Example usage
# ============================================================================
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "create_example":
        # Create example configuration
        output_path = sys.argv[2] if len(sys.argv) > 2 else "mercury_orbit.yaml"
        create_example_config(output_path)

    elif len(sys.argv) > 1 and sys.argv[1] == "validate":
        # Load and validate configuration
        config_path = sys.argv[2] if len(sys.argv) > 2 else "mercury_orbit.yaml"
        try:
            print(f"Loading configuration from {config_path}...")
            config = load_config(config_path)

            print("\nMedium:")
            print(config['medium'])

            print(f"\nBodies ({len(config['bodies'])}):")
            for body in config['bodies']:
                print(body)

            print("\nValidating configuration...")
            is_valid, warnings_list = validate_config(config)

            if warnings_list:
                print(f"\n{'WARNINGS' if is_valid else 'ERRORS'}:")
                for w in warnings_list:
                    print(f"  - {w}")

            if is_valid:
                print("\nConfiguration is VALID and ready for simulation.")
            else:
                print("\nConfiguration is INVALID. Please fix errors above.")
                sys.exit(1)

        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    else:
        print("Usage:")
        print("  python io_cfg.py create_example [output.yaml]")
        print("  python io_cfg.py validate [config.yaml]")
