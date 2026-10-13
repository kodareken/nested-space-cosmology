#!/usr/bin/env python3
"""Physical metric and the linear trust-region step. No Dirac."""
import sys
from pathlib import Path

import numpy as np
import pytest

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
    replay = Value.reconstruct_recorded_gradient(
        'step-u0-h5e-4',
        source_context='329ef809e466a86ed4a290ee32c19485f85add9d',
    )
    gradient = replay['gradient']
    provenance = replay['provenance']
    assert provenance['numeric_inputs_authenticated_from_current_local_bytes'] is True
    assert len(provenance['current_executable_math_original_AST_and_globals_verified']) == 9
    assert provenance['creating_producer_authenticated'] is False
    assert provenance['unverified_producer_bindings'] == [{
        'path': 'scripts/derive_nsc_ks_gate_value.py',
        'expected_sha256': '757a0192f3ca0073e62c2157bf13b98ef4d62136ae0e7fa320e071504e914f2f',
        'authenticated': False,
        'reason': 'original creating-CLI source bytes are unavailable in the preserved snapshot and bounded source history',
    }]
    assert provenance['physical_local_gate'] == 'OPEN'
    assert provenance['new_operator_solves'] == 0
    assert gradient.shape == (47, 2)
    # Preserve the original captured literal and its original allowance.
    # Independent current arithmetic has a separately reported finite-recipe
    # comparison envelope; neither allowance is a physical gate tolerance.
    diagnostic = provenance['numeric_replay']
    assert abs(diagnostic['recorded_constraint_maxima'][0] - 0.0006395923256182039) < 1e-18
    assert diagnostic['current_merit'] == float(np.max(np.abs(gradient)))
    assert abs(diagnostic['current_merit'] - diagnostic['recorded_merit']) <= diagnostic['merit_bound']
    assert diagnostic['physical_error_bound'] is None
    assert diagnostic['phase_quadrature_error_bound'] is None
    assert diagnostic['physical_constraint_status'] == 'OPEN'


def test_original_math_guard_rejects_algorithm_and_global_mutations():
    original = b"from hashlib import sha256\nscale = 2\ndef value(x):\n    return scale * x\n"
    Value._require_original_module_math(
        original, original + b"from hashlib import sha256\nnew_control = 3\n",
    )
    for altered in (
            original.replace(b"scale * x", b"scale + x"),
            original + b"scale = 3\n",
            original + b"del scale\n",
            original + b"from foreign_hash import sha256\n",
            original + b"scale.append(3)\n",
            original + b"def new_control(x=scale.append(3)): pass\n"):
        with pytest.raises(ValueError, match='mathematics|global bindings'):
            Value._require_original_module_math(original, altered)


def test_recorded_gradient_boundary_never_accepts_a_corrupt_numeric_hash():
    producer = 'scripts/derive_nsc_ks_gate_value.py'
    numeric = 'results/development/nsc-ks-gate-value-step-u0-h5e-4.json'
    report = {}
    with pytest.raises(ValueError, match='legacy recorded data hash mismatch'):
        Value._authenticate_legacy_replay_sources({
            producer: '757a0192f3ca0073e62c2157bf13b98ef4d62136ae0e7fa320e071504e914f2f',
            numeric: '0' * 64,
        }, source_context='329ef809e466a86ed4a290ee32c19485f85add9d',
            producer_report=report)
    assert report == {}


def test_recorded_register_integrity_ignores_later_path_history():
    relative='results/development/nsc-ks-gate-value-step-u0-h5e-4.json'
    raw=(ROOT/relative).read_bytes()
    assert Value._authenticate_recorded_register(
        relative,raw,'329ef809e466a86ed4a290ee32c19485f85add9d'
    )=='c468452bb75c7c2531f321eaf81487fa0d1ae1edc0640919d00d6da9fce12553'
    # Even whitespace changes must fail before parsing or numeric work.
    with pytest.raises(ValueError,match='local recorded register differs'):
        Value._authenticate_recorded_register(relative,raw+b' ',
            '329ef809e466a86ed4a290ee32c19485f85add9d')
    with pytest.raises(ValueError,match='predeclared original register'):
        Value._authenticate_recorded_register(relative,raw,'0'*40)
