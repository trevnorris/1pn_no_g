#!/usr/bin/env python3
"""
Test script to demonstrate mass drift diagnostic functionality.

This script tests the compute_mass_drift() function with both constant-mass
and variable-mass scenarios to validate the implementation.
"""

import numpy as np
from slab.diagnostics import compute_mass_drift


def test_constant_mass():
    """Test mass drift for constant-mass trajectory (Solar System regime)."""
    print("=" * 80)
    print("TEST 1: Constant Mass Trajectory (Solar System regime)")
    print("=" * 80)
    print()

    # Create a trajectory with constant masses
    trajectory = {
        't': np.array([0.0, 1.0, 2.0, 3.0, 4.0]),
        'M': np.array([
            [1.0, 3.3e-7],  # t=0: Sun, Mercury
            [1.0, 3.3e-7],  # t=1
            [1.0, 3.3e-7],  # t=2
            [1.0, 3.3e-7],  # t=3
            [1.0, 3.3e-7],  # t=4
        ]),
        'x': np.zeros((5, 2, 3)),
        'v': np.zeros((5, 2, 3)),
        'Q': np.zeros((5, 2)),
    }

    body_names = ["Sun", "Mercury"]

    drift = compute_mass_drift(trajectory, body_names)

    print("Bodies:", drift['bodies'])
    print()
    print("Initial masses:", drift['initial_mass'])
    print("Final masses:  ", drift['final_mass'])
    print()
    print("Absolute drift:", drift['abs_drift'])
    print("Relative drift:", drift['rel_drift'])
    print("Max rel. drift:", drift['max_rel_drift'])
    print()

    # Verify all drifts are zero
    assert np.all(drift['abs_drift'] == 0.0), "Expected zero absolute drift"
    assert np.all(drift['rel_drift'] == 0.0), "Expected zero relative drift"
    assert np.all(drift['max_rel_drift'] == 0.0), "Expected zero max drift"

    print("✓ All masses conserved (drift = 0.0)")
    print("✓ Test PASSED: Constant-mass regime validated")
    print()


def test_variable_mass():
    """Test mass drift for variable-mass trajectory (with intake)."""
    print("=" * 80)
    print("TEST 2: Variable Mass Trajectory (with mass intake)")
    print("=" * 80)
    print()

    # Create a trajectory where Mercury accretes mass
    # Sun: constant at 1.0
    # Mercury: grows from 3.3e-7 to 3.63e-7 (10% increase)
    trajectory = {
        't': np.array([0.0, 1.0, 2.0, 3.0, 4.0]),
        'M': np.array([
            [1.0, 3.30e-7],  # t=0
            [1.0, 3.40e-7],  # t=1: +3% for Mercury
            [1.0, 3.50e-7],  # t=2: +6%
            [1.0, 3.55e-7],  # t=3: +7.5%
            [1.0, 3.63e-7],  # t=4: +10%
        ]),
        'x': np.zeros((5, 2, 3)),
        'v': np.zeros((5, 2, 3)),
        'Q': np.zeros((5, 2)),
    }

    body_names = ["Sun", "Mercury"]

    drift = compute_mass_drift(trajectory, body_names)

    print("Bodies:", drift['bodies'])
    print()
    print("Initial masses:", drift['initial_mass'])
    print("Final masses:  ", drift['final_mass'])
    print()
    print("Absolute drift:", drift['abs_drift'])
    print("Relative drift:", drift['rel_drift'])
    print("Max rel. drift:", drift['max_rel_drift'])
    print()

    # Verify Sun is constant, Mercury has drifted
    assert drift['rel_drift'][0] == 0.0, "Sun should have zero drift"
    assert np.isclose(drift['rel_drift'][1], 0.1), "Mercury should have ~10% drift"
    assert np.isclose(drift['max_rel_drift'][1], 0.1), "Mercury max drift should be ~10%"

    print("✓ Sun mass conserved (drift = 0.0)")
    print(f"✓ Mercury mass increased by {drift['rel_drift'][1]*100:.1f}%")
    print("✓ Test PASSED: Variable-mass regime validated")
    print()


def test_edge_cases():
    """Test edge cases: zero mass, single timestep, empty trajectory."""
    print("=" * 80)
    print("TEST 3: Edge Cases")
    print("=" * 80)
    print()

    # Test 1: Zero initial mass
    print("Edge case 1: Zero initial mass")
    trajectory_zero = {
        't': np.array([0.0, 1.0]),
        'M': np.array([[0.0], [0.1]]),  # Massless to massive
        'x': np.zeros((2, 1, 3)),
        'v': np.zeros((2, 1, 3)),
        'Q': np.zeros((2, 1)),
    }
    drift = compute_mass_drift(trajectory_zero, ["TestBody"])
    assert np.isinf(drift['rel_drift'][0]), "Zero mass should give inf drift"
    print("  ✓ Zero initial mass handled correctly (rel_drift = inf)")
    print()

    # Test 2: Single timestep
    print("Edge case 2: Single timestep")
    trajectory_single = {
        't': np.array([0.0]),
        'M': np.array([[1.0]]),
        'x': np.zeros((1, 1, 3)),
        'v': np.zeros((1, 1, 3)),
        'Q': np.zeros((1, 1)),
    }
    drift = compute_mass_drift(trajectory_single, ["TestBody"])
    assert drift['rel_drift'][0] == 0.0, "Single timestep should have zero drift"
    print("  ✓ Single timestep handled correctly (drift = 0)")
    print()

    # Test 3: Empty trajectory
    print("Edge case 3: Empty trajectory")
    trajectory_empty = {
        't': np.array([]),
        'M': np.zeros((0, 0)),
        'x': np.zeros((0, 0, 3)),
        'v': np.zeros((0, 0, 3)),
        'Q': np.zeros((0, 0)),
    }
    drift = compute_mass_drift(trajectory_empty, [])
    assert len(drift['rel_drift']) == 0, "Empty trajectory should return empty arrays"
    print("  ✓ Empty trajectory handled correctly")
    print()

    print("✓ All edge cases PASSED")
    print()


def test_pass_fail_criteria():
    """Test pass/fail thresholds for mass drift."""
    print("=" * 80)
    print("TEST 4: Pass/Fail Criteria")
    print("=" * 80)
    print()

    # Test different drift levels
    test_cases = [
        (1e-15, "PASS", "Negligible (roundoff error)"),
        (1e-13, "PASS", "Excellent (< 1e-12)"),
        (1e-10, "WARN", "Good (< 1e-9)"),
        (1e-7, "FAIL", "WARNING (> 1e-9)"),
        (1e-3, "FAIL", "FAIL (> 1e-6)"),
    ]

    for drift_level, expected_status, description in test_cases:
        trajectory = {
            't': np.array([0.0, 1.0]),
            'M': np.array([[1.0], [1.0 + drift_level]]),
            'x': np.zeros((2, 1, 3)),
            'v': np.zeros((2, 1, 3)),
            'Q': np.zeros((2, 1)),
        }
        drift = compute_mass_drift(trajectory, ["TestBody"])
        rel_drift = drift['rel_drift'][0]

        # Check threshold
        if rel_drift < 1e-12:
            status = "PASS"
        elif rel_drift < 1e-9:
            status = "WARN"
        else:
            status = "FAIL"

        assert status == expected_status, f"Expected {expected_status}, got {status}"
        print(f"  ΔM/M = {rel_drift:.2e}  [{status}]  {description}")

    print()
    print("✓ Pass/fail thresholds validated")
    print()


if __name__ == "__main__":
    print()
    print("╔════════════════════════════════════════════════════════════════════════════╗")
    print("║                  MASS DRIFT DIAGNOSTIC TEST SUITE                          ║")
    print("╚════════════════════════════════════════════════════════════════════════════╝")
    print()

    test_constant_mass()
    test_variable_mass()
    test_edge_cases()
    test_pass_fail_criteria()

    print("=" * 80)
    print("ALL TESTS PASSED")
    print("=" * 80)
    print()
    print("The compute_mass_drift() diagnostic is working correctly!")
    print()
