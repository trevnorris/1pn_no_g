"""Medium dataclass for superfluid slab hydrodynamics.

This module defines the ambient medium properties for the superfluid slab
in which bodies (fluid intakes) reside. The medium is characterized by:
- Ambient density (rho0)
- Sound speed (cs), which controls compressible corrections
- Mass-intake factor (beta0), relating mass to volumetric intake Q
- Optional rarefaction exponent (gamma_beta) for density-dependent beta

No gravitational constant G appears. The orbital constant K = rho0/(4*pi*beta0^2)
emerges from the superfluid momentum flux and replaces G in Newtonian form.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class Medium:
    """Ambient superfluid medium properties.

    The medium defines the background state of the superfluid slab.
    Bodies (fluid intakes) perturb this background.

    Attributes
    ----------
    rho0 : float
        Ambient (background) density [kg/m³ in SI, or code units].
        The undisturbed density far from all bodies.
    cs : float
        Sound speed [m/s in SI, or code units].
        Controls the speed of acoustic disturbances and sets the scale
        for compressible (finite-Mach) corrections to forces.
        For incompressible limit, use cs → ∞.
    beta0 : float
        Mass-intake factor [kg·s/m³ in SI, or code units].
        Relates body mass M to volumetric intake Q via M = beta0 * Q.
        Equivalently, Q/M = 1/beta0 sets the universal flow-per-mass.
    gamma_beta : float, optional
        Rarefaction exponent for density-dependent beta (default: 0.0).
        When gamma_beta = 0, beta is constant (Solar System regime).
        When gamma_beta ≠ 0, beta(rho) = beta0 * (rho/rho0)^(-gamma_beta),
        modeling rarefaction effects in extreme environments.

    Notes
    -----
    The orbital constant K replaces Newton's G:
        K = rho0 / (4 * pi * beta0^2)  [m³/(kg·s²) in SI]

    This is the only combination of medium parameters that appears in
    the incompressible force law:
        F_a = sum_b (K * M_a * M_b / r_ab^2) * r_hat_ab

    Compressible corrections scale with cs^(-2) and generate perihelion
    precession proportional to (K*M)/(a*cs^2).

    Examples
    --------
    >>> # Incompressible medium (effectively cs → ∞)
    >>> medium_inc = Medium(rho0=1.0, cs=1e10, beta0=1e8)
    >>> print(medium_inc)
    Medium(rho0=1.000e+00, cs=1.000e+10, beta0=1.000e+08, gamma_beta=0.00)
      K = 7.958e-17 [orbital constant, replaces G]

    >>> # Compressible medium with finite sound speed
    >>> medium_comp = Medium(rho0=1.0, cs=1e4, beta0=1e10)
    >>> print(f"K = {medium_comp.K:.3e}")
    K = 7.958e-22

    >>> # Medium with rarefaction (density-dependent beta)
    >>> medium_rare = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=1.0)
    >>> rho_test = 0.8 * medium_rare.rho0
    >>> print(f"beta(rho=0.8*rho0) = {medium_rare.beta(rho_test):.3e}")
    beta(rho=0.8*rho0) = 1.250e+10
    """

    rho0: float
    cs: float
    beta0: float
    gamma_beta: float = 0.0

    def __post_init__(self):
        """Validate medium parameters."""
        if self.rho0 <= 0:
            raise ValueError(f"Ambient density rho0 must be positive, got {self.rho0}")
        if self.cs <= 0:
            raise ValueError(f"Sound speed cs must be positive, got {self.cs}")
        if self.beta0 <= 0:
            raise ValueError(f"Mass-intake factor beta0 must be positive, got {self.beta0}")
        if self.gamma_beta < 0:
            raise ValueError(f"Rarefaction exponent gamma_beta must be non-negative, got {self.gamma_beta}")

    @property
    def K(self) -> float:
        """Orbital constant K = rho0/(4*pi*beta0^2).

        This is the fundamental constant that replaces Newton's G in the
        superfluid formulation. It has the same dimensions as G and appears
        in the force law in the same way:

            F_a = sum_b (K * M_a * M_b / r_ab^2) * r_hat_ab

        Returns
        -------
        float
            Orbital constant [m³/(kg·s²) in SI, or code units].

        Examples
        --------
        >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
        >>> K = medium.K
        >>> # For two equal masses M at separation r:
        >>> # F = K * M^2 / r^2
        """
        return self.rho0 / (4.0 * np.pi * self.beta0**2)

    def beta(self, rho: float) -> float:
        """Mass-intake factor at density rho, accounting for rarefaction.

        Implements the rarefaction model:
            beta(rho) = beta0 * (rho/rho0)^(-gamma_beta)

        When gamma_beta = 0 (default), this returns beta0 (constant beta).
        When gamma_beta > 0, beta increases as density decreases,
        modeling rarefaction effects.

        Parameters
        ----------
        rho : float
            Local density [same units as rho0].

        Returns
        -------
        float
            Effective mass-intake factor at density rho.

        Notes
        -----
        For typical applications (Solar System), gamma_beta = 0 and beta
        is constant. Non-zero gamma_beta may be relevant for:
        - Extreme density gradients near compact objects
        - Capacity saturation regimes (Ma ~ 1)
        - Testing sensitivity to equation-of-state variations

        Examples
        --------
        >>> # Constant beta (default)
        >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=0.0)
        >>> medium.beta(0.5)  # beta independent of rho
        10000000000.0

        >>> # Rarefaction: beta increases as rho decreases
        >>> medium_rare = Medium(rho0=1.0, cs=1e4, beta0=1e10, gamma_beta=1.0)
        >>> medium_rare.beta(0.5)  # beta doubles when rho halves
        20000000000.0
        """
        if self.gamma_beta == 0.0:
            return self.beta0
        else:
            return self.beta0 * (rho / self.rho0)**(-self.gamma_beta)

    def __str__(self) -> str:
        """Nice string representation for printing.

        Returns
        -------
        str
            Human-readable summary of medium properties.

        Examples
        --------
        >>> medium = Medium(rho0=1.0, cs=1e4, beta0=1e10)
        >>> print(medium)
        Medium(rho0=1.000e+00, cs=1.000e+04, beta0=1.000e+10, gamma_beta=0.00)
          K = 7.958e-22 [orbital constant, replaces G]
        """
        lines = [
            f"Medium(rho0={self.rho0:.3e}, cs={self.cs:.3e}, "
            f"beta0={self.beta0:.3e}, gamma_beta={self.gamma_beta:.2f})"
        ]
        lines.append(f"  K = {self.K:.3e} [orbital constant, replaces G]")
        if self.gamma_beta != 0.0:
            lines.append(f"  Rarefaction: beta(rho) = beta0 * (rho/rho0)^(-{self.gamma_beta:.2f})")
        return "\n".join(lines)

    def __repr__(self) -> str:
        """Unambiguous representation for debugging."""
        return (f"Medium(rho0={self.rho0!r}, cs={self.cs!r}, "
                f"beta0={self.beta0!r}, gamma_beta={self.gamma_beta!r})")
