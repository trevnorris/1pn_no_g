"""
Velocity field calculations for superfluid hydrodynamics orbit simulator.

This module computes the velocity field from potential flow with point sinks.
Bodies create radial inflow patterns that superpose linearly.

Key equations (from plan_no_pde.md § 2):
    φ(x) = -Σ_b s_b/(4π|x - x_b|)
    v(x) = ∇φ = Σ_b (s_b/4π) r_b/r_b³
    s_b ≡ Q_b/ρ₀

Sign convention: sinks have Q > 0, producing inward velocity field.
For vector r_b = x - x_b (from body to field point), v points inward
when Q_b > 0.

Physics notes:
- Potential flow: irrotational, curl(v) = 0
- Incompressible away from sinks: div(v) = -Σ Q_b δ³(x - x_b)
- Each sink creates v ∝ 1/r² radial field
- Linear superposition holds in the ideal (ρ = ρ₀) limit
"""

import numpy as np
from typing import List, Tuple, Optional
from numpy.typing import NDArray


def v_self(
    x: NDArray[np.float64],
    x_body: NDArray[np.float64],
    Q: float,
    rho0: float,
    eps: float = 1e-12
) -> NDArray[np.float64]:
    """
    Velocity field from a single sink at x_body.

    Computes the radial inflow velocity at point x due to a single
    sink of strength Q located at x_body:

        v(x) = (Q/4πρ₀) * r/r³

    where r = x - x_body is the vector from the sink to the field point.

    Physical interpretation:
    - For Q > 0 (sink): velocity points inward toward x_body
    - Magnitude decays as 1/r²
    - Solid angle integral gives total flux Q

    Args:
        x: Field point position, shape (3,)
        x_body: Body (sink) position, shape (3,)
        Q: Volumetric intake rate [m³/s] (positive for sink)
        rho0: Background fluid density [kg/m³]
        eps: Regularization parameter for r → 0 (default 1e-12)

    Returns:
        Velocity vector at x, shape (3,). Returns zero if r < eps.

    Notes:
        - When x = x_body (r = 0), returns zero vector to avoid singularity
        - In practice, forces are computed on control surfaces at finite radius
        - The factor (Q/4πρ₀) has units [m³/s] = [velocity × area]

    Examples:
        >>> x = np.array([1.0, 0.0, 0.0])
        >>> x_body = np.array([0.0, 0.0, 0.0])
        >>> Q = 1.0
        >>> rho0 = 1.0
        >>> v = v_self(x, x_body, Q, rho0)
        >>> np.allclose(v, np.array([1.0/(4*np.pi), 0.0, 0.0]))
        True
    """
    # Vector from body to field point
    r_vec = x - x_body
    r = np.linalg.norm(r_vec)

    # Handle singularity at r = 0
    if r < eps:
        return np.zeros(3, dtype=np.float64)

    # v = (Q/4πρ₀) * r_vec/r³
    # Factor out 1/r³ = 1/(r² * r)
    prefactor = Q / (4.0 * np.pi * rho0)
    v = prefactor * r_vec / (r * r * r)

    return v


def v_ext_at(
    x_a: NDArray[np.float64],
    bodies: List,
    a_idx: int,
    rho0: float,
    eps: float = 1e-12
) -> NDArray[np.float64]:
    """
    External velocity at position x_a due to all bodies EXCEPT body a_idx.

    This is the key quantity for force calculations. By the control-surface
    lemma (plan_no_pde.md § 3.2, eq. 4):

        F_a^(inc) = (4/3) * (Q_a/4π) * v_ext(x_a)

    The external velocity is the superposition of all other sinks:

        v_ext(x_a) = Σ_{b≠a} (Q_b/4πρ₀) * r_ab/r_ab³

    where r_ab = x_a - x_b is the vector from body b to body a.

    Args:
        x_a: Position of body a, shape (3,)
        bodies: List of body objects with attributes .x (position) and .Q (intake)
        a_idx: Index of body a (to exclude from sum)
        rho0: Background fluid density [kg/m³]
        eps: Regularization for near-coincident bodies (default 1e-12)

    Returns:
        External velocity vector at x_a, shape (3,)

    Notes:
        - Self-field is excluded to compute forces properly
        - In two-body case: v_ext at body 1 comes only from body 2
        - Linear superposition of all external contributions

    Examples:
        >>> # Two bodies at distance d
        >>> from types import SimpleNamespace
        >>> b1 = SimpleNamespace(x=np.array([0., 0., 0.]), Q=1.0)
        >>> b2 = SimpleNamespace(x=np.array([1., 0., 0.]), Q=1.0)
        >>> bodies = [b1, b2]
        >>> v = v_ext_at(b1.x, bodies, 0, rho0=1.0)
        >>> # Should point from b2 toward b1 (in -x direction)
        >>> v[0] < 0  # True
    """
    v_ext = np.zeros(3, dtype=np.float64)

    for b_idx, body in enumerate(bodies):
        # Skip self
        if b_idx == a_idx:
            continue

        # Add contribution from body b
        v_b = v_self(x_a, body.x, body.Q, rho0, eps=eps)
        v_ext += v_b

    return v_ext


def v_total(
    x: NDArray[np.float64],
    bodies: List,
    rho0: float,
    eps: float = 1e-12
) -> NDArray[np.float64]:
    """
    Total velocity field at arbitrary point x from all bodies.

    Computes the full superfluid velocity field by summing contributions
    from all sinks:

        v(x) = Σ_b (Q_b/4πρ₀) * (x - x_b)/|x - x_b|³

    This is used for:
    - Field visualization
    - Computing velocity at control surface points for quadrature
    - Diagnostics and validation

    Args:
        x: Field point position, shape (3,)
        bodies: List of body objects with attributes .x and .Q
        rho0: Background fluid density [kg/m³]
        eps: Regularization parameter (default 1e-12)

    Returns:
        Total velocity vector at x, shape (3,)

    Notes:
        - At a body position, includes self-field (singular)
        - For force calculations, use v_ext_at instead
        - Far from all bodies, field decays as 1/r²

    Examples:
        >>> # Single body at origin
        >>> from types import SimpleNamespace
        >>> body = SimpleNamespace(x=np.array([0., 0., 0.]), Q=1.0)
        >>> x = np.array([1., 0., 0.])
        >>> v = v_total(x, [body], rho0=1.0)
        >>> # Should be radial inward
        >>> v[0] < 0  # True
    """
    v_sum = np.zeros(3, dtype=np.float64)

    for body in bodies:
        v_b = v_self(x, body.x, body.Q, rho0, eps=eps)
        v_sum += v_b

    return v_sum


def v_ext_vectorized(
    bodies: List,
    rho0: float,
    eps: float = 1e-12
) -> NDArray[np.float64]:
    """
    Compute external velocity for ALL bodies at once (at their own positions).

    This is the performance-critical function for the main integration loop.
    For N bodies, computes v_ext(x_a) for all a simultaneously using
    vectorized operations.

    Returns v_ext[a] = Σ_{b≠a} (Q_b/4πρ₀) * r_ab/r_ab³

    Performance strategy:
    - Build pairwise separation vectors (N×N×3 array)
    - Compute all pairwise velocities at once
    - Mask out self-interactions
    - Sum over source body index

    Args:
        bodies: List of N body objects with attributes .x and .Q
        rho0: Background fluid density [kg/m³]
        eps: Regularization for numerical stability (default 1e-12)

    Returns:
        Array of shape (N, 3) where v_ext[a] is the external velocity
        at body a's position due to all other bodies.

    Notes:
        - This is called every integration step
        - Vectorization provides ~10-100× speedup vs loops
        - Memory scales as O(N²) for separation vectors
        - For N > 1000, consider chunking or tree methods

    Algorithm:
        1. Extract positions into (N, 3) array
        2. Compute separations: r[a,b] = x[a] - x[b]
        3. Compute distances: d[a,b] = |r[a,b]|
        4. Compute velocities: v[a,b] = (Q[b]/4πρ₀) * r[a,b]/d[a,b]³
        5. Mask diagonal (self-interactions)
        6. Sum over source index b: v_ext[a] = Σ_b v[a,b]

    Examples:
        >>> from types import SimpleNamespace
        >>> b1 = SimpleNamespace(x=np.array([0., 0., 0.]), Q=1.0)
        >>> b2 = SimpleNamespace(x=np.array([1., 0., 0.]), Q=1.0)
        >>> v_ext = v_ext_vectorized([b1, b2], rho0=1.0)
        >>> v_ext.shape
        (2, 3)
    """
    N = len(bodies)

    if N == 0:
        return np.zeros((0, 3), dtype=np.float64)

    if N == 1:
        # Single body has no external field
        return np.zeros((1, 3), dtype=np.float64)

    # Extract positions and intakes
    # positions[a] = x_a, shape (N, 3)
    positions = np.array([body.x for body in bodies], dtype=np.float64)

    # intakes[b] = Q_b, shape (N,)
    intakes = np.array([body.Q for body in bodies], dtype=np.float64)

    # Compute all pairwise separations
    # r_ab[a, b] = x_a - x_b, shape (N, N, 3)
    # Broadcasting: positions[:, None, :] has shape (N, 1, 3)
    #               positions[None, :, :] has shape (1, N, 3)
    r_ab = positions[:, None, :] - positions[None, :, :]

    # Compute all pairwise distances
    # dist[a, b] = |x_a - x_b|, shape (N, N)
    dist = np.linalg.norm(r_ab, axis=2)

    # Regularize near-zero distances (diagonal and any coincident bodies)
    # We'll mask these out, but avoid divide-by-zero warnings
    dist = np.where(dist < eps, eps, dist)

    # Compute velocity contributions
    # v[a, b, :] = (Q_b/4πρ₀) * r_ab[a,b] / dist[a,b]³
    # Shape: (N, N, 3)
    prefactor = intakes / (4.0 * np.pi * rho0)  # shape (N,)

    # dist³ with shape (N, N)
    dist_cubed = dist * dist * dist

    # Broadcast prefactor[b] and divide by dist³[a,b]
    # prefactor[None, :] has shape (1, N)
    # dist_cubed has shape (N, N)
    # Result: shape (N, N)
    coeff = prefactor[None, :] / dist_cubed  # shape (N, N)

    # Multiply by separation vectors
    # coeff[:, :, None] has shape (N, N, 1)
    # r_ab has shape (N, N, 3)
    v_pairwise = coeff[:, :, None] * r_ab  # shape (N, N, 3)

    # Mask out self-interactions (diagonal)
    mask = ~np.eye(N, dtype=bool)  # False on diagonal
    v_pairwise_masked = np.where(mask[:, :, None], v_pairwise, 0.0)

    # Sum over source body index (axis 1)
    # v_ext[a] = Σ_b v_pairwise_masked[a, b]
    v_ext = np.sum(v_pairwise_masked, axis=1)  # shape (N, 3)

    return v_ext


# ============================================================================
# Utility functions for field analysis and validation
# ============================================================================

def v_magnitude(
    x: NDArray[np.float64],
    bodies: List,
    rho0: float
) -> float:
    """
    Compute magnitude of total velocity field at point x.

    Convenience function for scalar field visualization.

    Args:
        x: Field point, shape (3,)
        bodies: List of body objects
        rho0: Background density

    Returns:
        Speed |v(x)| as scalar float
    """
    v = v_total(x, bodies, rho0)
    return float(np.linalg.norm(v))


def potential(
    x: NDArray[np.float64],
    bodies: List,
    rho0: float,
    eps: float = 1e-12
) -> float:
    """
    Velocity potential φ(x) at point x.

    From plan_no_pde.md § 2:
        φ(x) = -Σ_b (s_b/4π) / |x - x_b|

    where s_b = Q_b/ρ₀.

    Note: Sign convention is chosen so that ∇φ gives inward velocity
    for sinks (Q > 0).

    Args:
        x: Field point, shape (3,)
        bodies: List of body objects
        rho0: Background density
        eps: Regularization parameter

    Returns:
        Scalar potential value

    Notes:
        - Potential is negative for sinks (Q > 0)
        - Related to field energy density: ½ ρ₀ |∇φ|²
    """
    phi = 0.0

    for body in bodies:
        r_vec = x - body.x
        r = np.linalg.norm(r_vec)

        if r < eps:
            # At body location, potential is singular
            # Return large negative value or raise warning
            continue

        s_b = body.Q / rho0
        phi -= s_b / (4.0 * np.pi * r)

    return phi


def check_divergence_free(
    x: NDArray[np.float64],
    bodies: List,
    rho0: float,
    delta: float = 1e-6
) -> float:
    """
    Numerical check that div(v) ≈ 0 away from bodies.

    Uses finite differences to estimate ∇·v at point x.
    Should be zero (within numerical error) if x is not near a body.

    Args:
        x: Field point, shape (3,)
        bodies: List of body objects
        rho0: Background density
        delta: Finite difference step size

    Returns:
        Estimated div(v) at x

    Notes:
        - Near a body (within ~delta), this will be large
        - Use for validation: check far-field is divergence-free
    """
    # Central differences in each direction
    divv = 0.0

    for i in range(3):
        x_plus = x.copy()
        x_minus = x.copy()
        x_plus[i] += delta
        x_minus[i] -= delta

        v_plus = v_total(x_plus, bodies, rho0)
        v_minus = v_total(x_minus, bodies, rho0)

        # ∂v_i/∂x_i
        dvi_dxi = (v_plus[i] - v_minus[i]) / (2.0 * delta)
        divv += dvi_dxi

    return divv


def check_curl_free(
    x: NDArray[np.float64],
    bodies: List,
    rho0: float,
    delta: float = 1e-6
) -> NDArray[np.float64]:
    """
    Numerical check that curl(v) ≈ 0 (irrotational flow).

    Uses finite differences to estimate ∇×v at point x.
    Should be zero vector everywhere for potential flow.

    Args:
        x: Field point, shape (3,)
        bodies: List of body objects
        rho0: Background density
        delta: Finite difference step size

    Returns:
        Estimated curl(v) at x, shape (3,)

    Notes:
        - Exact potential flow has curl(v) = 0 everywhere
        - Nonzero result indicates numerical error
        - Use for validation of field calculations
    """
    curl = np.zeros(3, dtype=np.float64)

    # curl(v) = (∂v_z/∂y - ∂v_y/∂z, ∂v_x/∂z - ∂v_z/∂x, ∂v_y/∂x - ∂v_x/∂y)

    # x-component: ∂v_z/∂y - ∂v_y/∂z
    x_plus_y = x.copy(); x_plus_y[1] += delta
    x_minus_y = x.copy(); x_minus_y[1] -= delta
    dvz_dy = (v_total(x_plus_y, bodies, rho0)[2] -
              v_total(x_minus_y, bodies, rho0)[2]) / (2.0 * delta)

    x_plus_z = x.copy(); x_plus_z[2] += delta
    x_minus_z = x.copy(); x_minus_z[2] -= delta
    dvy_dz = (v_total(x_plus_z, bodies, rho0)[1] -
              v_total(x_minus_z, bodies, rho0)[1]) / (2.0 * delta)

    curl[0] = dvz_dy - dvy_dz

    # y-component: ∂v_x/∂z - ∂v_z/∂x
    x_plus_x = x.copy(); x_plus_x[0] += delta
    x_minus_x = x.copy(); x_minus_x[0] -= delta
    dvx_dz = (v_total(x_plus_z, bodies, rho0)[0] -
              v_total(x_minus_z, bodies, rho0)[0]) / (2.0 * delta)
    dvz_dx = (v_total(x_plus_x, bodies, rho0)[2] -
              v_total(x_minus_x, bodies, rho0)[2]) / (2.0 * delta)

    curl[1] = dvx_dz - dvz_dx

    # z-component: ∂v_y/∂x - ∂v_x/∂y
    dvy_dx = (v_total(x_plus_x, bodies, rho0)[1] -
              v_total(x_minus_x, bodies, rho0)[1]) / (2.0 * delta)
    dvx_dy = (v_total(x_plus_y, bodies, rho0)[0] -
              v_total(x_minus_y, bodies, rho0)[0]) / (2.0 * delta)

    curl[2] = dvy_dx - dvx_dy

    return curl
