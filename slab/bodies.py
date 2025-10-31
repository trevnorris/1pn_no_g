"""Body dataclass for superfluid slab hydrodynamics.

This module defines bodies as localized fluid intakes ("mouths") that
remove superfluid at a rate Q_a. Each body has:
- Mass M_a and volumetric intake Q_a, related by M_a = beta * Q_a
- Position x_a and velocity v_a (3D vectors)
- Control surface radius R_a for force/flux calculations

Bodies do NOT generate gravitational fields. Instead, they perturb the
superfluid flow, and forces arise from momentum flux through control surfaces.
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class Body:
    """A body (fluid intake) in the superfluid slab.

    Bodies are "mouths" that remove superfluid at volumetric rate Q,
    perturbing the ambient flow. Forces on bodies arise from external
    momentum flux through small control surfaces of radius R.

    Attributes
    ----------
    name : str
        Identifier for this body (e.g., "Sun", "Mercury").
    M : float
        Inertial mass [kg in SI, or code units].
        Related to volumetric intake by M = beta * Q.
    x : np.ndarray
        Position vector [m in SI, or code units], shape (3,).
    v : np.ndarray
        Velocity vector [m/s in SI, or code units], shape (3,).
    R : float
        Control surface radius [m in SI, or code units].
        Must be small compared to body separations but large enough
        to enclose the mouth. Typical choice: R ~ 10^-3 * r_min.
    Q : float
        Volumetric intake [m³/s in SI, or code units].
        Positive for sinks (standard), negative for sources (exotic).
        Automatically synchronized with M via update methods.

    Notes
    -----
    **Mass-intake synchronization**:
    The relation M = beta * Q must be maintained consistently.
    Use the update methods:
    - `update_Q_from_M(medium)` after changing M
    - `update_M_from_Q(medium)` after changing Q

    **Control surface radius R**:
    - Must satisfy: mouth_size << R << separation to nearest body
    - Too small R: numerical issues in quadrature
    - Too large R: breakdown of near-field expansion
    - Typical: R ~ 10^-4 to 10^-3 in code units with separation ~ 1

    **Velocity field singularity**:
    At the mouth (r=0), v diverges as 1/r². The control surface R
    is chosen to avoid the singularity while staying in the linear
    near-field regime where v_ext is approximately constant.

    Examples
    --------
    >>> from slab.medium import Medium
    >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>>
    >>> # Create a body with mass, let Q be computed
    >>> sun = Body(
    ...     name="Sun",
    ...     M=1.0,
    ...     x=np.array([0.0, 0.0, 0.0]),
    ...     v=np.array([0.0, 0.0, 0.0]),
    ...     R=1e-3,
    ...     Q=0.0  # placeholder
    ... )
    >>> sun.update_Q_from_M(medium)
    >>> print(f"Sun: M={sun.M}, Q={sun.Q:.3e}")
    Sun: M=1.0, Q=1.000e-10
    >>>
    >>> # Create Mercury with semi-major axis a=0.387 AU (code units)
    >>> a = 0.387
    >>> v_circ = np.sqrt(medium.K * sun.M / a)
    >>> mercury = Body(
    ...     name="Mercury",
    ...     M=1.66e-7,
    ...     x=np.array([a, 0.0, 0.0]),
    ...     v=np.array([0.0, v_circ, 0.0]),
    ...     R=5e-4,
    ...     Q=0.0
    ... )
    >>> mercury.update_Q_from_M(medium)
    >>> print(mercury)
    Body 'Mercury': M=1.660e-07, Q=1.660e-17
      x = [3.870e-01, 0.000e+00, 0.000e+00]
      v = [0.000e+00, ..., 0.000e+00]
      R = 5.000e-04
      KE = ...
    """

    name: str
    M: float
    x: np.ndarray
    v: np.ndarray
    R: float
    Q: float = 0.0

    def __post_init__(self):
        """Validate body parameters and ensure arrays are proper numpy arrays."""
        # Convert to numpy arrays if needed
        if not isinstance(self.x, np.ndarray):
            self.x = np.array(self.x, dtype=float)
        if not isinstance(self.v, np.ndarray):
            self.v = np.array(self.v, dtype=float)

        # Validate shapes
        if self.x.shape != (3,):
            raise ValueError(f"Position x must have shape (3,), got {self.x.shape}")
        if self.v.shape != (3,):
            raise ValueError(f"Velocity v must have shape (3,), got {self.v.shape}")

        # Validate parameters
        if self.M <= 0:
            raise ValueError(f"Mass M must be positive, got {self.M}")
        if self.R <= 0:
            raise ValueError(f"Control surface radius R must be positive, got {self.R}")

        # Q can be positive (sink, standard) or negative (source, exotic)
        # but we warn if Q and M have inconsistent signs
        if self.Q != 0.0 and np.sign(self.Q) != np.sign(self.M):
            import warnings
            warnings.warn(
                f"Body '{self.name}': Q={self.Q} and M={self.M} have opposite signs. "
                "Standard bodies are sinks with Q>0, M>0.",
                UserWarning
            )

    def update_Q_from_M(self, medium) -> None:
        """Synchronize Q from current M using Q = M / beta0.

        Call this after updating M (e.g., via mass intake integration)
        to keep Q consistent.

        Parameters
        ----------
        medium : Medium
            The ambient medium containing beta0.

        Notes
        -----
        Uses constant beta0 from the medium. For density-dependent beta,
        this should be called with the appropriate local density value
        passed to medium.beta(rho).

        Examples
        --------
        >>> from slab.medium import Medium
        >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
        >>> body = Body("Test", M=1.0, x=[0,0,0], v=[0,0,0], R=1e-3, Q=0.0)
        >>> body.update_Q_from_M(medium)
        >>> body.Q
        1e-10
        """
        self.Q = self.M / medium.beta0

    def update_M_from_Q(self, medium) -> None:
        """Synchronize M from current Q using M = beta0 * Q.

        Call this after updating Q (e.g., via flux calculations)
        to keep M consistent.

        Parameters
        ----------
        medium : Medium
            The ambient medium containing beta0.

        Notes
        -----
        Uses constant beta0 from the medium. For density-dependent beta,
        this should be called with the appropriate local density value
        passed to medium.beta(rho).

        Examples
        --------
        >>> from slab.medium import Medium
        >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
        >>> body = Body("Test", M=1.0, x=[0,0,0], v=[0,0,0], R=1e-3, Q=1e-10)
        >>> body.Q = 2e-10  # Intake increased
        >>> body.update_M_from_Q(medium)
        >>> body.M
        2.0
        """
        self.M = medium.beta0 * self.Q

    @property
    def kinetic_energy(self) -> float:
        """Kinetic energy KE = (1/2) * M * v².

        Returns
        -------
        float
            Kinetic energy [kg·m²/s² in SI, or code units].

        Examples
        --------
        >>> body = Body("Test", M=1.0, x=[0,0,0], v=[1,0,0], R=1e-3, Q=0.0)
        >>> body.kinetic_energy
        0.5
        >>>
        >>> body.v = np.array([3.0, 4.0, 0.0])
        >>> body.kinetic_energy
        12.5
        """
        v_squared = np.dot(self.v, self.v)
        return 0.5 * self.M * v_squared

    def __str__(self) -> str:
        """Nice string representation for printing.

        Returns
        -------
        str
            Human-readable summary of body state.

        Examples
        --------
        >>> body = Body(
        ...     name="Earth",
        ...     M=3.0e-6,
        ...     x=np.array([1.0, 0.0, 0.0]),
        ...     v=np.array([0.0, 6.28, 0.0]),
        ...     R=1e-3,
        ...     Q=3e-16
        ... )
        >>> print(body)
        Body 'Earth': M=3.000e-06, Q=3.000e-16
          x = [1.000e+00, 0.000e+00, 0.000e+00]
          v = [0.000e+00, 6.280e+00, 0.000e+00]
          R = 1.000e-03
          KE = 5.918e-05
        """
        lines = [
            f"Body '{self.name}': M={self.M:.3e}, Q={self.Q:.3e}"
        ]
        lines.append(f"  x = [{self.x[0]:.3e}, {self.x[1]:.3e}, {self.x[2]:.3e}]")
        lines.append(f"  v = [{self.v[0]:.3e}, {self.v[1]:.3e}, {self.v[2]:.3e}]")
        lines.append(f"  R = {self.R:.3e}")
        lines.append(f"  KE = {self.kinetic_energy:.3e}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        """Unambiguous representation for debugging."""
        return (
            f"Body(name={self.name!r}, M={self.M!r}, "
            f"x={self.x!r}, v={self.v!r}, R={self.R!r}, Q={self.Q!r})"
        )
