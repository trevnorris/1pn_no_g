"""Shared helpers for perihelion precession diagnostics.

These utilities centralise the setup of barycentric Sun--Mercury initial
conditions and the extraction of osculating orbital elements from simulation
trajectories.  They are imported by diagnostic scripts such as
``final_diagnosis.py`` and ``validate_1pn_precession.py`` so that every tool
uses the exact same conventions when reporting precession rates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import copy

import numpy as np

from slab.bodies import Body
from slab.dynamics import integrate_orbit
from slab.gr1pn import compute_orbit_elements
from slab.medium import Medium


@dataclass(frozen=True)
class TwoBodyParams:
    """Orbital parameters for a barycentric two-body configuration."""

    a: float
    e: float
    r_peri: float
    r_apo: float
    v_rel_peri: float
    T_orbit: float
    K: float
    M_primary: float
    M_secondary: float

    @property
    def M_total(self) -> float:
        return self.M_primary + self.M_secondary


def create_barycentric_mercury_config(
    eccentricity: float,
    *,
    rho0: float = 1.0,
    beta0: float = 0.045,
    cs: float = 63239.7263,
    semi_major_axis: float = 0.387,
    M_primary: float = 1.0,
    M_secondary: float = 3.3e-7,
) -> Tuple[Medium, List[Body], TwoBodyParams]:
    """Return a Mercury-like configuration in the barycentric frame.

    The initial conditions place the Sun and Mercury on opposite sides of the
    barycentre so that the centre-of-mass remains at rest throughout the
    integration.  This avoids spuriously attributing barycentric drift to
    periapsis motion when we compute osculating elements from the lab-frame
    trajectory.
    """

    medium = Medium(rho0=rho0, cs=cs, beta0=beta0, gamma_beta=0.0)
    K = medium.K

    a = semi_major_axis
    e = eccentricity
    r_peri = a * (1.0 - e)
    r_apo = a * (1.0 + e)

    M_total = M_primary + M_secondary
    mu = K * M_total
    v_rel_peri = np.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))

    x_primary = np.array([-M_secondary / M_total * r_peri, 0.0, 0.0])
    v_primary = np.array([0.0, -M_secondary / M_total * v_rel_peri, 0.0])
    x_secondary = np.array([M_primary / M_total * r_peri, 0.0, 0.0])
    v_secondary = np.array([0.0, M_primary / M_total * v_rel_peri, 0.0])

    primary = Body(
        name="Sun",
        M=M_primary,
        x=x_primary,
        v=v_primary,
        R=0.001,
        Q=M_primary / beta0,
    )

    secondary = Body(
        name="Mercury",
        M=M_secondary,
        x=x_secondary,
        v=v_secondary,
        R=0.0005,
        Q=M_secondary / beta0,
    )

    T_orbit = 2.0 * np.pi * np.sqrt(a**3 / (K * M_total))

    params = TwoBodyParams(
        a=a,
        e=e,
        r_peri=r_peri,
        r_apo=r_apo,
        v_rel_peri=v_rel_peri,
        T_orbit=T_orbit,
        K=K,
        M_primary=M_primary,
        M_secondary=M_secondary,
    )

    return medium, [primary, secondary], params


def compute_relative_state(traj: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return relative Sun→Mercury positions and velocities for a trajectory."""

    x_primary = traj["x"][:, 0, :]
    x_secondary = traj["x"][:, 1, :]
    v_primary = traj["v"][:, 0, :]
    v_secondary = traj["v"][:, 1, :]

    x_rel = x_secondary - x_primary
    v_rel = v_secondary - v_primary

    return traj["t"], x_rel, v_rel


def analyse_precession(
    traj: Dict[str, np.ndarray],
    params: TwoBodyParams,
) -> Dict[str, float]:
    """Extract multiple perihelion-precession estimators from a trajectory."""

    t, x_rel, v_rel = compute_relative_state(traj)

    omega = []
    elems_list = []
    for xr, vr in zip(x_rel, v_rel):
        elems = compute_orbit_elements(xr, vr, params.M_total, params.K)
        omega.append(elems["omega"])
        elems_list.append(elems)

    omega = np.asarray(omega)
    omega_unwrapped = np.unwrap(omega)

    if len(t) > 1:
        coeffs = np.polyfit(t, omega_unwrapped, deg=1)
        slope_per_orbit = coeffs[0] * params.T_orbit
    else:
        slope_per_orbit = np.nan

    r = np.linalg.norm(x_rel, axis=1)
    dr = np.diff(r)
    peri_mask = (dr[:-1] < 0) & (dr[1:] >= 0)
    peri_indices = np.where(peri_mask)[0] + 1

    if len(peri_indices) >= 2:
        peri_omega = np.unwrap(omega[peri_indices])
        peri_deltas = np.diff(peri_omega)
        peri_mean = peri_deltas.mean()
        peri_std = peri_deltas.std(ddof=1) if len(peri_deltas) > 1 else 0.0
    else:
        peri_mean = np.nan
        peri_std = np.nan

    return {
        "t": t,
        "omega": omega,
        "omega_unwrapped": omega_unwrapped,
        "elements": elems_list,
        "slope_per_orbit": slope_per_orbit,
        "peri_per_orbit": peri_mean,
        "peri_std": peri_std,
        "n_peri_cycles": max(0, len(peri_indices) - 1),
    }


def integrate_two_body(
    bodies: Iterable[Body],
    medium: Medium,
    *,
    dt: float,
    n_steps: int,
    save_every: int = 1,
    use_compressible: bool = True,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Convenience wrapper around :func:`integrate_orbit` for two bodies."""

    bodies_copy = [copy.deepcopy(b) for b in bodies]
    opts = {
        "use_compressible": use_compressible,
        "save_every": save_every,
        "use_quadrature": False,
        "verbose": False,
    }
    return integrate_orbit(bodies_copy, medium, dt, n_steps, opts)

