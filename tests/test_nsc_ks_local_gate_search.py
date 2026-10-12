"""Admission, fresh acceptance, tangent control and authenticated resume."""
from hashlib import sha256
import json

import numpy as np
import pytest

from recursive_horizons.evidence_io import canonical_json_bytes
from recursive_horizons.nsc_ks_local_gate_search import (
    COMPONENTS,
    FullEvaluation,
    Gate2Admission,
    JacobianCheckError,
    SearchAdmissionError,
    SearchBinding,
    SearchBindingError,
    SearchConfig,
    ValueEvaluation,
    admit_gate2_budget,
    coefficients_sha256,
    run_search,
    verify_jacobian,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def digest(raw):
    return sha256(raw).hexdigest()


def write_budget(path, *, missing=False):
    components = {
        name: {"bound": ([None, None] if missing and name == COMPONENTS[0]
                          else [1e-13, 2e-13])}
        for name in reversed(COMPONENTS)
    }
    record = {
        "schema": "NSC-KS-GATE-ERROR-BUDGET-v5",
        "status": "PASS: nine components enclosed",
        "components": components,
        "missing_components": ([COMPONENTS[0]] if missing else []),
        "covered_positive_families": 60,
        "covered_signed_families": 120,
        "all_declared_source_paths_closed": True,
        "full_error_sum": [1e-12, 2e-12],
        "allocation_total": 2e-11,
        "residual_reserve": 1e-11,
        "physical_tolerance": 3e-11,
        "constraint_order": ["N", "beta"],
        "state_law": "C_Sigma[g]=U_g C_up U_g^dagger",
        "numerical_indicators_used_as_bounds": False,
    }
    path.write_text(json.dumps(record, sort_keys=True))
    return record


class LinearEvaluator:
    def __init__(self, target=1e-4, support_limit=None, wrong_jacobian=False):
        self.evaluator_identity = digest(b"synthetic-local-gate-evaluator")
        self.matrix = np.zeros((258, 64))
        self.matrix[:64] = np.eye(64)
        for row in range(64, 258):
            self.matrix[row, row % 64] = 0.1 / (1 + row)
        self.target = np.full(64, float(target))
        self.support_limit = support_limit
        self.wrong_jacobian = wrong_jacobian
        self.calls = 0

    def preparation_supported(self, family):
        if self.support_limit is None:
            return True
        return float(np.max(np.abs(np.asarray(family.coefficients)))) <= self.support_limit

    def _gradient(self, family):
        vector = np.asarray(family.coefficients, float).ravel()
        flat = self.matrix @ (vector - self.target)
        return np.column_stack((flat[:129], flat[129:]))

    def _receipt(self, family, kind):
        self.calls += 1
        coefficients = np.asarray(family.coefficients, float)
        coefficients_digest = coefficients_sha256(coefficients)
        identity = digest(
            kind.encode() + coefficients_digest.encode() + str(self.calls).encode())
        return coefficients, coefficients_digest, identity

    def value(self, family, *, fresh):
        coefficients, coefficients_digest, identity = self._receipt(family, "value")
        return ValueEvaluation(
            self._gradient(family), identity, coefficients_digest,
            self.evaluator_identity, 0.0, fresh)

    def full(self, family, *, fresh):
        coefficients, coefficients_digest, identity = self._receipt(family, "full")
        jacobian = np.empty((64, 129, 2))
        jacobian[:, :, 0] = self.matrix[:129].T
        jacobian[:, :, 1] = self.matrix[129:].T
        if self.wrong_jacobian:
            jacobian[0, 0, 0] += 1.0
        return FullEvaluation(
            self._gradient(family), identity, coefficients_digest,
            self.evaluator_identity, 0.0, fresh, jacobian)


def context(tmp_path, evaluator, *, support_budget=True):
    tmp_path.mkdir(parents=True, exist_ok=True)
    budget_path = tmp_path / "budget.json"
    write_budget(budget_path)
    admission = admit_gate2_budget(budget_path)
    initial = np.zeros((2, 32))
    identity = profile_identity(
        LocalIncomingFamily(initial), include_normal_window=True)
    dependencies = (("results/development/nsc-ks-gate-budget-v5.json",
                     admission.budget_sha256),)
    if not support_budget:
        dependencies = ((dependencies[0][0], digest(b"wrong budget")),)
    binding = SearchBinding(
        identity,
        digest(b"fixed upstream source"),
        evaluator.evaluator_identity,
        digest(b"evaluator implementation"),
        digest(b"solver settings"),
        dependencies,
    )
    config = SearchConfig(initial_trust_radius=1.0, maximum_trust_radius=2.0)
    return initial, np.eye(64), config, binding, admission


def test_gate2_admission_requires_all_nine_finite_components(tmp_path):
    good = tmp_path / "good.json"
    write_budget(good)
    admission = admit_gate2_budget(good)
    assert isinstance(admission, Gate2Admission)
    assert admission.component_sum == (1e-12, 2e-12)
    bad = tmp_path / "bad.json"
    write_budget(bad, missing=True)
    with pytest.raises(SearchAdmissionError, match="missing"):
        admit_gate2_budget(bad)


def test_linear_rectangular_search_accepts_only_fresh_value_and_resumes(tmp_path):
    evaluator = LinearEvaluator()
    initial, metric, config, binding, admission = context(tmp_path, evaluator)
    journal = tmp_path / "journal"
    result = run_search(
        evaluator, initial, metric, config, binding, admission, journal,
        trial_limit=1)
    assert result.status.startswith("NODAL_TARGET_REACHED")
    assert result.accepted_trials == 1
    assert result.merit < 1e-15
    event = json.loads((journal / "event-000000.json").read_text())
    assert event["accepted"] is True
    assert event["fresh_evaluation_identity"]
    assert event["jacobian_check"]["selected_columns"][:6] == [0, 15, 31, 32, 47, 63]
    assert len(event["jacobian_check"]["selected_columns"]) == 8
    resumed = run_search(
        evaluator, initial, metric, config, binding, admission, journal)
    assert resumed.trials == 1
    assert resumed.accepted_trials == 1
    np.testing.assert_allclose(np.asarray(resumed.coefficients), np.full((2, 32), 1e-4),
                               rtol=0, atol=1e-15)


def test_support_rejection_is_recorded_without_candidate_evolution(tmp_path):
    evaluator = LinearEvaluator(target=1e-4, support_limit=1e-5)
    initial, metric, config, binding, admission = context(tmp_path, evaluator)
    journal = tmp_path / "journal"
    run_search(evaluator, initial, metric, config, binding, admission, journal,
               trial_limit=1)
    event = json.loads((journal / "event-000000.json").read_text())
    assert event["accepted"] is False
    assert event["preparation_supported"] is False
    assert event["fresh_evaluation_identity"] is None
    assert event["trust_radius_after"] < event["trust_radius_before"]


def test_resume_rejects_mutated_event_and_budget_binding(tmp_path):
    evaluator = LinearEvaluator()
    initial, metric, config, binding, admission = context(tmp_path, evaluator)
    journal = tmp_path / "journal"
    run_search(evaluator, initial, metric, config, binding, admission, journal,
               trial_limit=1)
    event_path = journal / "event-000000.json"
    event = json.loads(event_path.read_text())
    event["binding_sha256"] = digest(b"mutated")
    event_path.write_bytes(canonical_json_bytes(event))
    with pytest.raises(SearchBindingError, match="event body digest|event chain or binding"):
        run_search(evaluator, initial, metric, config, binding, admission, journal,
                   trial_limit=0)
    _initial, _metric, _config, wrong, _admission = context(
        tmp_path / "other", evaluator, support_budget=False)
    with pytest.raises(SearchBindingError, match="budget digest"):
        run_search(evaluator, initial, metric, config, wrong, admission,
                   tmp_path / "other-journal", trial_limit=0)


def test_centered_control_rejects_wrong_full_jacobian(tmp_path):
    evaluator = LinearEvaluator(wrong_jacobian=True)
    initial, _metric, config, binding, _admission = context(tmp_path, evaluator)
    family = LocalIncomingFamily(initial)
    full = evaluator.full(family, fresh=True)
    with pytest.raises(JacobianCheckError, match="column 0"):
        verify_jacobian(evaluator, initial, full, config, binding)
