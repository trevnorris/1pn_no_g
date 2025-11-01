"""Visualization module for superfluid orbit simulator.

This module provides plotting functions to compare slab simulation trajectories
with GR 1PN predictions. It supports:
- Orbit comparison plots (2D trajectories + separation over time)
- Precession comparison (omega vs time with linear fits)
- 3D trajectory visualization

Design principles:
- Publication-ready plots (clear labels, legends, titles)
- Consistent color scheme (slab=blue, GR=red)
- High-DPI output (dpi=150 or 300)
- Grid lines and proper units
"""

from typing import Dict, List, Optional
import numpy as np
from pathlib import Path

# Import matplotlib with proper backend handling
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available. Visualization functions will not work.")


def _check_matplotlib():
    """Check if matplotlib is available, raise ImportError if not."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for visualization. "
            "Install with: pip install matplotlib"
        )


def plot_orbit_comparison(
    slab_trajectory: Dict[str, np.ndarray],
    gr_trajectory: Dict[str, np.ndarray],
    output_path: str,
    body_idx: int = 1,
    dpi: int = 150
) -> None:
    """Plot two-panel comparison of slab and GR trajectories.

    Creates a publication-quality figure with:
    - Top panel: Overlaid x-y trajectories (slab=blue solid, GR=red dashed)
    - Bottom panel: Separation |r_slab - r_gr| over time

    Parameters
    ----------
    slab_trajectory : Dict[str, np.ndarray]
        Slab trajectory from dynamics.integrate_orbit().
        Must contain 't', 'x', 'v' arrays.
    gr_trajectory : Dict[str, np.ndarray]
        GR trajectory from gr1pn.integrate_gr1pn_orbit().
        Must contain 't', 'x', 'v' arrays.
    output_path : str
        Output file path (e.g., "output/orbit_comparison.png")
    body_idx : int, optional
        Index of body to plot (default: 1 for orbiting body)
    dpi : int, optional
        Output resolution (default: 150)

    Notes
    -----
    The trajectories should have the same initial conditions and be sampled
    at similar time intervals for meaningful comparison. The separation plot
    shows how the two solutions diverge over time, which should be small if
    the slab correctly reproduces GR effects.

    Examples
    --------
    >>> # After running both simulations
    >>> plot_orbit_comparison(slab_traj, gr_traj, "output/orbit_comparison.png")
    Saved orbit comparison plot to output/orbit_comparison.png
    """
    _check_matplotlib()

    # Extract trajectories for the specified body
    t_slab = slab_trajectory['t']
    x_slab = slab_trajectory['x'][:, body_idx, :]

    t_gr = gr_trajectory['t']
    x_gr = gr_trajectory['x'][:, body_idx, :]

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Top panel: Overlaid x-y trajectories
    ax1.plot(x_slab[:, 0], x_slab[:, 1], 'b-', linewidth=1.5,
             label='Slab simulation', alpha=0.8)
    ax1.plot(x_gr[:, 0], x_gr[:, 1], 'r--', linewidth=1.5,
             label='GR 1PN', alpha=0.8)
    ax1.plot(x_slab[0, 0], x_slab[0, 1], 'go', markersize=8,
             label='Initial position')
    ax1.plot(0, 0, 'k*', markersize=15, label='Central body')

    ax1.set_xlabel('x [AU]', fontsize=12)
    ax1.set_ylabel('y [AU]', fontsize=12)
    ax1.set_title('Orbit Comparison: Slab vs GR 1PN', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')

    # Bottom panel: Separation over time
    # Interpolate GR trajectory to slab time points if needed
    if len(t_slab) != len(t_gr) or not np.allclose(t_slab, t_gr):
        # Linear interpolation of GR trajectory to slab time points
        x_gr_interp = np.zeros_like(x_slab)
        for dim in range(3):
            x_gr_interp[:, dim] = np.interp(t_slab, t_gr, x_gr[:, dim])
        separation = np.linalg.norm(x_slab - x_gr_interp, axis=1)
        t_sep = t_slab
    else:
        separation = np.linalg.norm(x_slab - x_gr, axis=1)
        t_sep = t_slab

    ax2.plot(t_sep, separation, 'k-', linewidth=1.5)
    ax2.set_xlabel('Time [yr]', fontsize=12)
    ax2.set_ylabel('Separation [AU]', fontsize=12)
    ax2.set_title('Position Difference: |r_slab - r_gr|', fontsize=12)
    ax2.grid(True, alpha=0.3)

    # Add statistics
    max_sep = np.max(separation)
    mean_sep = np.mean(separation)
    textstr = f'Max: {max_sep:.3e} AU\nMean: {mean_sep:.3e} AU'
    ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round',
             facecolor='wheat', alpha=0.5), fontsize=10)

    plt.tight_layout()

    # Save figure
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"Saved orbit comparison plot to {output_path}")


def plot_precession_comparison(
    slab_trajectory: Dict[str, np.ndarray],
    gr_trajectory: Dict[str, np.ndarray],
    body_idx: int,
    M_central: float,
    K_slab: float,
    G_gr: float,
    output_path: str,
    dpi: int = 150
) -> None:
    """Plot perihelion precession comparison between slab and GR.

    Computes omega(t) (argument of periapsis) for both trajectories, fits
    linear trends to extract precession rates, and plots comparison with
    agreement ratio.

    Parameters
    ----------
    slab_trajectory : Dict[str, np.ndarray]
        Slab trajectory from dynamics.integrate_orbit()
    gr_trajectory : Dict[str, np.ndarray]
        GR trajectory from gr1pn.integrate_gr1pn_orbit()
    body_idx : int
        Index of orbiting body
    M_central : float
        Mass of central body (same units for both slab and GR)
    K_slab : float
        Slab orbital constant (replaces G)
    G_gr : float
        Gravitational constant for GR calculation
    output_path : str
        Output file path (e.g., "output/precession_comparison.png")
    dpi : int, optional
        Output resolution (default: 150)

    Notes
    -----
    The precession rate is computed by:
    1. Computing osculating orbital elements at each time step
    2. Extracting omega (argument of periapsis)
    3. Unwrapping 2π discontinuities
    4. Fitting linear trend: omega(t) = omega_0 + (dω/dt) * t
    5. Converting to precession per orbit

    The agreement ratio should be close to 1.0 if the slab correctly
    reproduces GR 1PN precession.

    Examples
    --------
    >>> plot_precession_comparison(
    ...     slab_traj, gr_traj, body_idx=1,
    ...     M_central=1.0, K_slab=7.96e-22, G_gr=6.67e-11,
    ...     output_path="output/precession.png"
    ... )
    Saved precession comparison plot to output/precession.png
    """
    _check_matplotlib()

    from slab.gr1pn import compute_orbit_elements

    # Extract trajectories
    t_slab = slab_trajectory['t']
    x_slab = slab_trajectory['x'][:, body_idx, :]
    v_slab = slab_trajectory['v'][:, body_idx, :]

    t_gr = gr_trajectory['t']
    x_gr = gr_trajectory['x'][:, body_idx, :]
    v_gr = gr_trajectory['v'][:, body_idx, :]

    # Compute omega for slab trajectory
    omega_slab = []
    for i in range(len(t_slab)):
        elems = compute_orbit_elements(x_slab[i], v_slab[i], M_central, K_slab)
        omega_slab.append(elems['omega'])
    omega_slab = np.array(omega_slab)

    # Compute omega for GR trajectory
    omega_gr = []
    for i in range(len(t_gr)):
        elems = compute_orbit_elements(x_gr[i], v_gr[i], M_central, G_gr)
        omega_gr.append(elems['omega'])
    omega_gr = np.array(omega_gr)

    # Unwrap 2π discontinuities
    omega_slab_unwrapped = np.unwrap(omega_slab)
    omega_gr_unwrapped = np.unwrap(omega_gr)

    # Fit linear trends
    coeffs_slab = np.polyfit(t_slab, omega_slab_unwrapped, deg=1)
    domega_dt_slab = coeffs_slab[0]
    omega0_slab = coeffs_slab[1]
    fit_slab = omega0_slab + domega_dt_slab * t_slab

    coeffs_gr = np.polyfit(t_gr, omega_gr_unwrapped, deg=1)
    domega_dt_gr = coeffs_gr[0]
    omega0_gr = coeffs_gr[1]
    fit_gr = omega0_gr + domega_dt_gr * t_gr

    # Estimate orbital period
    a_list = []
    for i in range(len(t_slab)):
        elems = compute_orbit_elements(x_slab[i], v_slab[i], M_central, K_slab)
        if np.isfinite(elems['a']):
            a_list.append(elems['a'])
    a_mean = np.mean(a_list) if a_list else 1.0
    P_orb_slab = 2 * np.pi * np.sqrt(a_mean**3 / (K_slab * M_central))

    # Compute precession per orbit
    domega_per_orbit_slab = domega_dt_slab * P_orb_slab
    domega_per_orbit_gr = domega_dt_gr * P_orb_slab  # Use same period for comparison

    # Agreement ratio
    if domega_per_orbit_gr != 0:
        ratio = domega_per_orbit_slab / domega_per_orbit_gr
    else:
        ratio = np.nan

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot omega vs time
    ax.plot(t_slab, omega_slab_unwrapped, 'b.', markersize=3,
            label='Slab data', alpha=0.5)
    ax.plot(t_slab, fit_slab, 'b-', linewidth=2,
            label=f'Slab fit: dω/dt = {domega_dt_slab:.3e} rad/yr')

    ax.plot(t_gr, omega_gr_unwrapped, 'r.', markersize=3,
            label='GR data', alpha=0.5)
    ax.plot(t_gr, fit_gr, 'r-', linewidth=2,
            label=f'GR fit: dω/dt = {domega_dt_gr:.3e} rad/yr')

    ax.set_xlabel('Time [yr]', fontsize=12)
    ax.set_ylabel('Argument of Periapsis ω [rad]', fontsize=12)
    ax.set_title('Perihelion Precession: Slab vs GR 1PN', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add text box with precession info
    textstr = (
        f'Precession per orbit:\n'
        f'  Slab: {domega_per_orbit_slab:.3e} rad/orbit\n'
        f'  Slab: {domega_per_orbit_slab * 206265:.3f} arcsec/orbit\n'
        f'  GR:   {domega_per_orbit_gr:.3e} rad/orbit\n'
        f'  GR:   {domega_per_orbit_gr * 206265:.3f} arcsec/orbit\n'
        f'Agreement ratio: {ratio:.4f}'
    )
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round',
            facecolor='lightblue', alpha=0.7), fontsize=10,
            family='monospace')

    plt.tight_layout()

    # Save figure
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"Saved precession comparison plot to {output_path}")
    print(f"  Slab precession: {domega_per_orbit_slab * 206265:.3f} arcsec/orbit")
    print(f"  GR precession:   {domega_per_orbit_gr * 206265:.3f} arcsec/orbit")
    print(f"  Agreement ratio: {ratio:.4f}")


def plot_trajectory_3d(
    trajectory: Dict[str, np.ndarray],
    body_names: List[str],
    output_path: str,
    dpi: int = 150
) -> None:
    """Plot 3D trajectory visualization with color-coded bodies.

    Creates an interactive 3D plot showing orbital trajectories with:
    - Color-coded paths for each body
    - Initial positions marked with circles
    - Final positions marked with squares
    - Central body marked with star

    Parameters
    ----------
    trajectory : Dict[str, np.ndarray]
        Trajectory from integrate_orbit() or integrate_gr1pn_orbit()
        Must contain 't', 'x' arrays
    body_names : List[str]
        Names of bodies (for legend)
    output_path : str
        Output file path (e.g., "output/trajectory_3d.png")
    dpi : int, optional
        Output resolution (default: 150)

    Notes
    -----
    This function is useful for visualizing precession and orbital plane
    changes in 3D. For planar orbits (z=0), the 3D visualization may not
    add much over 2D plots.

    Examples
    --------
    >>> plot_trajectory_3d(
    ...     slab_traj,
    ...     body_names=["Sun", "Mercury"],
    ...     output_path="output/trajectory_3d.png"
    ... )
    Saved 3D trajectory plot to output/trajectory_3d.png
    """
    _check_matplotlib()

    # Extract data
    x = trajectory['x']  # Shape: (n_steps, n_bodies, 3)
    n_steps, n_bodies, _ = x.shape

    # Define colors for bodies
    colors = ['orange', 'blue', 'red', 'green', 'purple', 'brown', 'pink']

    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot trajectories
    for i in range(n_bodies):
        color = colors[i % len(colors)]
        label = body_names[i] if i < len(body_names) else f"Body {i}"

        # Plot trajectory
        ax.plot(x[:, i, 0], x[:, i, 1], x[:, i, 2],
                color=color, linewidth=1.5, label=label, alpha=0.7)

        # Mark initial position
        ax.scatter(x[0, i, 0], x[0, i, 1], x[0, i, 2],
                   color=color, marker='o', s=100, edgecolors='black',
                   linewidth=1.5)

        # Mark final position
        ax.scatter(x[-1, i, 0], x[-1, i, 1], x[-1, i, 2],
                   color=color, marker='s', s=100, edgecolors='black',
                   linewidth=1.5)

    # Mark central body (body 0) with star at origin
    if n_bodies > 0:
        ax.scatter(x[0, 0, 0], x[0, 0, 1], x[0, 0, 2],
                   color='gold', marker='*', s=500, edgecolors='black',
                   linewidth=2, zorder=10)

    ax.set_xlabel('x [AU]', fontsize=12)
    ax.set_ylabel('y [AU]', fontsize=12)
    ax.set_zlabel('z [AU]', fontsize=12)
    ax.set_title('3D Orbital Trajectories', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Equal aspect ratio for all axes
    max_range = np.array([
        x[:, :, 0].max() - x[:, :, 0].min(),
        x[:, :, 1].max() - x[:, :, 1].min(),
        x[:, :, 2].max() - x[:, :, 2].min()
    ]).max() / 2.0

    mid_x = (x[:, :, 0].max() + x[:, :, 0].min()) * 0.5
    mid_y = (x[:, :, 1].max() + x[:, :, 1].min()) * 0.5
    mid_z = (x[:, :, 2].max() + x[:, :, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.tight_layout()

    # Save figure
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"Saved 3D trajectory plot to {output_path}")
