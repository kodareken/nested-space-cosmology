#!/usr/bin/env python3
"""Physical metric and the linear trust-region step. No Dirac."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_gate_lm as LM
import derive_nsc_ks_gate_value as Value
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_physical_metric_is_positive_definite():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    metric = LM.physical_metric(8, family.center, samples=201)
    assert metric.shape == (16, 16)
    assert np.linalg.eigvalsh(metric).min() > 0


def test_trust_step_reduces_the_linear_model():
    residual = np.array([0.4, -0.2, 0.1])
    jacobian = np.array([[1.0, 0.0], [0.2, 1.0], [0.0, -0.5]])
    metric = np.eye(2)
    step, predicted = LM.levenberg_step(residual, jacobian, metric, radius=10.0)
    assert np.linalg.norm(predicted) < np.linalg.norm(residual)
    assert step @ metric @ step <= 10.0**2 + 1e-8


def test_trust_step_stays_inside_a_small_radius():
    residual = np.array([1.0, 1.0])
    jacobian = np.eye(2)
    metric = np.eye(2)
    step, _predicted = LM.levenberg_step(residual, jacobian, metric, radius=0.05)
    assert abs(np.sqrt(step @ step) - 0.05) < 1e-8


def test_saved_u0_gradient_matches_its_merit():
    gradient = Value.reconstruct_gradient('step-u0-h5e-4')
    assert gradient.shape == (47, 2)
    assert abs(float(np.max(np.abs(gradient))) - 0.0006395923256182039) < 1e-18
