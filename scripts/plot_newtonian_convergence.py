#!/usr/bin/env python3
"""
Plot Newtonian baseline convergence results.

This script generates publication-quality figures showing:
1. Spurious precession vs timestep (log-log plot showing dt² scaling)
2. Energy conservation vs timestep
3. Comparison of Newtonian vs Superfluid methods
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

# Data from validation study (8 orbits, e=0.1)
dt_values = np.array([0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, 0.000005])
steps_per_orbit = np.array([120.7, 241.3, 482.6, 1206.5, 2413.0, 4826.1, 12065.2, 24130.4, 48260.7])

# Spurious precession [arcsec/orbit]
precession = np.array([898.985, 226.788, 56.473, 9.008, 2.277, 0.567, 0.092, 0.023, 0.0057])
precession_std = np.array([69.805, 8.717, 1.155, 0.296, 0.062, 0.0, 0.0035, 0.0009, 0.0])

# Energy drift
energy_drift = np.array([2.5e-7, 1.42e-8, 1.87e-9, 1.8e-12, 2.15e-10, 7.49e-14, 8.56e-12, 2.03e-12, 2.42e-13])

# GR prediction
GR_precession = 0.0997  # arcsec/orbit

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# ============================================================================
# Plot 1: Spurious Precession vs Timestep
# ============================================================================

# Main data
ax1.loglog(dt_values, np.abs(precession), 'o-',
           color='#2E86AB', markersize=8, linewidth=2,
           label='Spurious precession', zorder=3)

# Error bars (only where std > 0)
mask = precession_std > 0
ax1.errorbar(dt_values[mask], np.abs(precession[mask]),
             yerr=precession_std[mask],
             fmt='none', ecolor='#2E86AB', alpha=0.3, capsize=4, zorder=2)

# dt² reference line
dt_ref = np.array([0.00001, 0.002])
prec_ref = 0.023 * (dt_ref / 0.00001)**2
ax1.loglog(dt_ref, prec_ref, '--',
           color='gray', linewidth=1.5, alpha=0.7,
           label='$\\propto \\mathrm{d}t^2$ scaling')

# GR prediction line
ax1.axhline(GR_precession, color='#A23B72', linestyle='--',
            linewidth=2, label='GR 1PN signal', zorder=1)

# 10% and 1% lines
ax1.axhline(0.1 * GR_precession, color='#F18F01', linestyle=':',
            linewidth=1.5, label='10% of GR', zorder=1)
ax1.axhline(0.01 * GR_precession, color='#C73E1D', linestyle=':',
            linewidth=1.5, label='1% of GR', zorder=1)

# Recommended region
dt_rec = 5e-6
ax1.axvline(dt_rec, color='green', linestyle='-.',
            linewidth=1.5, alpha=0.5, label='Recommended: $\\mathrm{d}t \\leq 5\\times10^{-6}$ yr')

ax1.set_xlabel('Timestep, $\\mathrm{d}t$ [yr]', fontsize=12)
ax1.set_ylabel('Spurious precession [arcsec/orbit]', fontsize=12)
ax1.set_title('(a) Spurious Precession vs Timestep', fontsize=13, fontweight='bold')
ax1.legend(loc='upper left', fontsize=9, framealpha=0.95)
ax1.grid(True, which='both', alpha=0.3, linestyle=':', linewidth=0.5)
ax1.set_xlim([3e-6, 3e-3])
ax1.set_ylim([0.003, 2000])

# ============================================================================
# Plot 2: Energy Conservation
# ============================================================================

ax2.loglog(dt_values, np.abs(energy_drift), 's-',
           color='#06A77D', markersize=8, linewidth=2,
           label='Energy drift $|\\Delta E/E|$', zorder=3)

# dt² reference line
dt_ref = np.array([0.000005, 0.002])
E_ref = 2.42e-13 * (dt_ref / 0.000005)**2
ax2.loglog(dt_ref, E_ref, '--',
           color='gray', linewidth=1.5, alpha=0.7,
           label='$\\propto \\mathrm{d}t^2$ scaling')

# Precision targets
ax2.axhline(1e-10, color='#F18F01', linestyle=':',
            linewidth=1.5, label='Good: $10^{-10}$', zorder=1)
ax2.axhline(1e-12, color='#06A77D', linestyle=':',
            linewidth=1.5, label='Excellent: $10^{-12}$', zorder=1)

# Recommended region
ax2.axvline(dt_rec, color='green', linestyle='-.',
            linewidth=1.5, alpha=0.5, label='Recommended: $\\mathrm{d}t \\leq 5\\times10^{-6}$ yr')

ax2.set_xlabel('Timestep, $\\mathrm{d}t$ [yr]', fontsize=12)
ax2.set_ylabel('Energy drift, $|\\Delta E/E|$', fontsize=12)
ax2.set_title('(b) Energy Conservation vs Timestep', fontsize=13, fontweight='bold')
ax2.legend(loc='lower right', fontsize=9, framealpha=0.95)
ax2.grid(True, which='both', alpha=0.3, linestyle=':', linewidth=0.5)
ax2.set_xlim([3e-6, 3e-3])
ax2.set_ylim([1e-14, 1e-6])

plt.tight_layout()
plt.savefig('newtonian_convergence.png', dpi=300, bbox_inches='tight')
plt.savefig('newtonian_convergence.pdf', bbox_inches='tight')
print("Saved: newtonian_convergence.png")
print("Saved: newtonian_convergence.pdf")

# ============================================================================
# Additional plot: Steps per orbit vs Precision
# ============================================================================

fig2, ax = plt.subplots(figsize=(8, 6))

# Precession vs steps per orbit
ax.loglog(steps_per_orbit, np.abs(precession), 'o-',
          color='#2E86AB', markersize=8, linewidth=2,
          label='Spurious precession', zorder=3)

# GR and target lines
ax.axhline(GR_precession, color='#A23B72', linestyle='--',
           linewidth=2, label='GR 1PN signal (0.1 arcsec/orbit)', zorder=1)
ax.axhline(0.1 * GR_precession, color='#F18F01', linestyle=':',
           linewidth=1.5, label='10% of GR (0.01 arcsec/orbit)', zorder=1)
ax.axhline(0.01 * GR_precession, color='#C73E1D', linestyle=':',
           linewidth=1.5, label='1% of GR (0.001 arcsec/orbit)', zorder=1)

# Recommended minimum
steps_rec = 48261
ax.axvline(steps_rec, color='green', linestyle='-.',
           linewidth=1.5, alpha=0.5,
           label=f'Recommended: $\\geq${steps_rec:.0f} steps/orbit')

ax.set_xlabel('Steps per orbit', fontsize=12)
ax.set_ylabel('Spurious precession [arcsec/orbit]', fontsize=12)
ax.set_title('Resolution Requirements for 1PN Studies', fontsize=13, fontweight='bold')
ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax.grid(True, which='both', alpha=0.3, linestyle=':', linewidth=0.5)
ax.set_xlim([100, 60000])
ax.set_ylim([0.003, 2000])

plt.tight_layout()
plt.savefig('steps_per_orbit_requirements.png', dpi=300, bbox_inches='tight')
plt.savefig('steps_per_orbit_requirements.pdf', bbox_inches='tight')
print("Saved: steps_per_orbit_requirements.png")
print("Saved: steps_per_orbit_requirements.pdf")

plt.show()

# ============================================================================
# Print summary table
# ============================================================================

print()
print("=" * 80)
print("NEWTONIAN BASELINE VALIDATION SUMMARY")
print("=" * 80)
print()
print(f"{'dt [yr]':<12} {'Steps/orbit':<12} {'Precession [arcsec]':<20} {'ΔE/E':<12}")
print("-" * 80)
for i in range(len(dt_values)):
    print(f"{dt_values[i]:<12.6f} {steps_per_orbit[i]:<12.1f} "
          f"{precession[i]:<20.6f} {energy_drift[i]:<12.2e}")
print("-" * 80)
print()
print(f"GR 1PN prediction: {GR_precession:.6f} arcsec/orbit")
print(f"Recommended minimum: dt ≤ {dt_rec:.6f} yr ({steps_rec:.0f} steps/orbit)")
print(f"  → Achieves {precession[-1]:.6f} arcsec/orbit spurious precession")
print(f"  → This is {precession[-1]/GR_precession*100:.1f}% of the GR signal")
print()
print("=" * 80)
