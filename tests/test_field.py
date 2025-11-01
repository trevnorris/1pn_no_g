"""
Tests for velocity field calculations.

Validates:
1. Single-sink velocity field (1/r² decay)
2. External velocity computation
3. Vectorized performance function
4. Field properties (curl-free, divergence at sinks)
"""

import numpy as np
import pytest
from types import SimpleNamespace
import sys
sys.path.insert(0, '/var/projects/1pn_no_g')

from slab.field import (
    v_self,
    v_ext_at,
    v_total,
    v_ext_vectorized,
    v_magnitude,
    potential,
    check_curl_free,
)


class TestVSelf:
    """Tests for single-sink velocity field."""

    def test_radial_direction(self):
        """Velocity should point radially from body to field point."""
        x_body = np.array([0.0, 0.0, 0.0])
        x = np.array([1.0, 0.0, 0.0])
        Q = 1.0
        rho0 = 1.0

        v = v_self(x, x_body, Q, rho0)

        # Should point toward the sink (negative x direction)
        assert v[0] < 0
        assert abs(v[1]) < 1e-12
        assert abs(v[2]) < 1e-12

    def test_inverse_square_decay(self):
        """Velocity magnitude should decay as 1/r²."""
        x_body = np.array([0.0, 0.0, 0.0])
        Q = 1.0
        rho0 = 1.0

        r1 = 1.0
        r2 = 2.0

        x1 = np.array([r1, 0.0, 0.0])
        x2 = np.array([r2, 0.0, 0.0])

        v1 = v_self(x1, x_body, Q, rho0)
        v2 = v_self(x2, x_body, Q, rho0)

        mag1 = np.linalg.norm(v1)
        mag2 = np.linalg.norm(v2)

        # mag1/mag2 should equal (r2/r1)²
        ratio = mag1 / mag2
        expected_ratio = (r2 / r1) ** 2

        assert abs(ratio - expected_ratio) < 1e-10

    def test_zero_at_body(self):
        """Velocity should be zero at body position (singularity handled)."""
        x_body = np.array([1.0, 2.0, 3.0])
        Q = 1.0
        rho0 = 1.0

        v = v_self(x_body, x_body, Q, rho0)

        assert np.allclose(v, 0.0)

    def test_magnitude_formula(self):
        """Check exact magnitude: |v| = Q/(4πr²)."""
        x_body = np.array([0.0, 0.0, 0.0])
        x = np.array([2.0, 0.0, 0.0])
        Q = 3.0
        rho0 = 2.0

        v = v_self(x, x_body, Q, rho0)
        mag = np.linalg.norm(v)

        r = 2.0
        expected_mag = Q / (4.0 * np.pi * r * r)

        assert abs(mag - expected_mag) < 1e-10


class TestVExtAt:
    """Tests for external velocity at body positions."""

    def test_two_body_symmetry(self):
        """Two equal sinks should produce symmetric external velocities."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1.0, 0.0, 0.0]), Q=1.0)
        bodies = [b1, b2]
        rho0 = 1.0

        v_ext_1 = v_ext_at(b1.x, bodies, 0, rho0)
        v_ext_2 = v_ext_at(b2.x, bodies, 1, rho0)

        # Magnitudes should be equal
        assert abs(np.linalg.norm(v_ext_1) - np.linalg.norm(v_ext_2)) < 1e-10

        # Directions should be opposite
        v_ext_1_norm = v_ext_1 / np.linalg.norm(v_ext_1)
        v_ext_2_norm = v_ext_2 / np.linalg.norm(v_ext_2)

        assert np.allclose(v_ext_1_norm, -v_ext_2_norm)

    def test_excludes_self(self):
        """External velocity should not include self-field."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        bodies = [b1]
        rho0 = 1.0

        v_ext = v_ext_at(b1.x, bodies, 0, rho0)

        # Single body has no external field
        assert np.allclose(v_ext, 0.0)

    def test_three_body_superposition(self):
        """External velocity should be sum of individual contributions."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1.0, 0.0, 0.0]), Q=2.0)
        b3 = SimpleNamespace(x=np.array([0.0, 1.0, 0.0]), Q=3.0)
        bodies = [b1, b2, b3]
        rho0 = 1.0

        # External velocity at b1 should equal sum from b2 and b3
        v_ext = v_ext_at(b1.x, bodies, 0, rho0)

        v_from_b2 = v_self(b1.x, b2.x, b2.Q, rho0)
        v_from_b3 = v_self(b1.x, b3.x, b3.Q, rho0)
        v_expected = v_from_b2 + v_from_b3

        assert np.allclose(v_ext, v_expected)


class TestVTotal:
    """Tests for total velocity field."""

    def test_single_body_matches_v_self(self):
        """For single body, v_total should equal v_self."""
        body = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        bodies = [body]
        rho0 = 1.0

        x = np.array([1.0, 1.0, 0.0])

        v_tot = v_total(x, bodies, rho0)
        v_s = v_self(x, body.x, body.Q, rho0)

        assert np.allclose(v_tot, v_s)

    def test_superposition(self):
        """Total field should be sum of individual contributions."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([2.0, 0.0, 0.0]), Q=2.0)
        bodies = [b1, b2]
        rho0 = 1.0

        x = np.array([1.0, 1.0, 0.0])

        v_tot = v_total(x, bodies, rho0)
        v1 = v_self(x, b1.x, b1.Q, rho0)
        v2 = v_self(x, b2.x, b2.Q, rho0)

        assert np.allclose(v_tot, v1 + v2)


class TestVExtVectorized:
    """Tests for vectorized external velocity computation."""

    def test_matches_loop_version(self):
        """Vectorized version should match loop-based computation."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1.0, 0.0, 0.0]), Q=2.0)
        b3 = SimpleNamespace(x=np.array([0.0, 1.0, 0.0]), Q=3.0)
        bodies = [b1, b2, b3]
        rho0 = 1.0

        # Vectorized version
        v_ext_vec = v_ext_vectorized(bodies, rho0)

        # Loop version
        v_ext_loop = np.array([
            v_ext_at(bodies[i].x, bodies, i, rho0)
            for i in range(len(bodies))
        ])

        assert np.allclose(v_ext_vec, v_ext_loop)

    def test_single_body(self):
        """Single body should have zero external velocity."""
        body = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        bodies = [body]
        rho0 = 1.0

        v_ext = v_ext_vectorized(bodies, rho0)

        assert v_ext.shape == (1, 3)
        assert np.allclose(v_ext, 0.0)

    def test_two_body_opposite_directions(self):
        """Two bodies should have opposite external velocities."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1.0, 0.0, 0.0]), Q=1.0)
        bodies = [b1, b2]
        rho0 = 1.0

        v_ext = v_ext_vectorized(bodies, rho0)

        # Normalize and check opposite
        v1_norm = v_ext[0] / np.linalg.norm(v_ext[0])
        v2_norm = v_ext[1] / np.linalg.norm(v_ext[1])

        assert np.allclose(v1_norm, -v2_norm)

    def test_empty_bodies(self):
        """Empty body list should return empty array."""
        bodies = []
        rho0 = 1.0

        v_ext = v_ext_vectorized(bodies, rho0)

        assert v_ext.shape == (0, 3)


class TestFieldProperties:
    """Tests for physical properties of the field."""

    def test_curl_free_far_field(self):
        """Flow should be irrotational (curl = 0) away from bodies."""
        body = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        bodies = [body]
        rho0 = 1.0

        # Point far from body
        x = np.array([10.0, 10.0, 10.0])

        curl = check_curl_free(x, bodies, rho0, delta=1e-4)

        # Should be zero within numerical error
        assert np.linalg.norm(curl) < 1e-6

    def test_potential_gradient_matches_velocity(self):
        """∇φ should equal velocity field."""
        body = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        bodies = [body]
        rho0 = 1.0

        x = np.array([1.0, 0.5, 0.2])
        delta = 1e-6

        # Compute velocity
        v = v_total(x, bodies, rho0)

        # Compute ∇φ by finite differences
        grad_phi = np.zeros(3)
        for i in range(3):
            x_plus = x.copy()
            x_minus = x.copy()
            x_plus[i] += delta
            x_minus[i] -= delta

            phi_plus = potential(x_plus, bodies, rho0)
            phi_minus = potential(x_minus, bodies, rho0)

            grad_phi[i] = (phi_plus - phi_minus) / (2.0 * delta)

        # v = ∇φ
        assert np.allclose(v, grad_phi, atol=1e-6)


class TestEdgeCases:
    """Tests for edge cases and numerical stability."""

    def test_very_small_Q(self):
        """Should handle very small intake rates."""
        body = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1e-20)
        bodies = [body]
        rho0 = 1.0

        x = np.array([1.0, 0.0, 0.0])
        v = v_total(x, bodies, rho0)

        # Should be very small but finite
        assert np.isfinite(v).all()
        assert np.linalg.norm(v) < 1e-18

    def test_large_separation(self):
        """Should handle large separations correctly."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1e6, 0.0, 0.0]), Q=1.0)
        bodies = [b1, b2]
        rho0 = 1.0

        v_ext = v_ext_vectorized(bodies, rho0)

        # Should be very small but finite
        assert np.isfinite(v_ext).all()

    def test_near_coincident_bodies(self):
        """Should handle near-coincident bodies with regularization."""
        b1 = SimpleNamespace(x=np.array([0.0, 0.0, 0.0]), Q=1.0)
        b2 = SimpleNamespace(x=np.array([1e-15, 0.0, 0.0]), Q=1.0)
        bodies = [b1, b2]
        rho0 = 1.0

        # Should not crash or produce inf/nan
        v_ext = v_ext_vectorized(bodies, rho0, eps=1e-12)

        assert np.isfinite(v_ext).all()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
