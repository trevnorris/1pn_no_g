#!/usr/bin/env python3
"""
Check if the precession measurement itself is the issue.

With only ~8 orbits and a precession of ~0.1 arcsec/orbit,
the total precession is ~0.8 arcsec ~ 4e-6 radians.

This is TINY compared to measurement noise from:
1. Integration errors
2. Non-exact initial conditions
3. Numerical precision in omega extraction

Let's check what we're actually measuring.
