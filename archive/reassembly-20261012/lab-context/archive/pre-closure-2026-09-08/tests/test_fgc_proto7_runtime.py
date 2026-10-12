from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    AcceptedStateSourceFailureEvidence,
    ProposalFailureEvidence,
    Proto7TerminalStop,
    attempt_proto7_step,
)


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(
    *,
    source_calls: set[int] | None = None,
    lapse_calls: set[int] | None = None,
    kinetic_calls: set[int] | None = None,
):
    calls = 0
    source = set() if source_calls is None else set(source_calls)
    lapse = set() if lapse_calls is None else set(lapse_calls)
    kinetic = set() if kinetic_calls is None else set(kinetic_calls)

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        nonlocal calls
        index = calls
        calls += 1
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero,
            zero,
            zero,
            {
                "source_residual_infinity": 2.0e-12 if index in source else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 2.0e10 if index in kinetic else 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if index in lapse else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


class FGCProto7RuntimeTests(unittest.TestCase):
    def transaction(self) -> GR0RuntimeStageTransaction:
        return GR0RuntimeStageTransaction(
            thresholds=GR0UniversalThresholds(),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
            grid_spacing=0.125,
            cfl_maximum=0.125,
        )

    def attempt(self, rhs, transaction, *, preaccept=None):
        return attempt_proto7_step(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.01,
            state=_state(),
            rhs=rhs,
            projector=None,
            transaction=transaction,
            previous_step_index=0,
            previous_transaction_serial=0,
            preaccept=preaccept,
        )

    def test_accepted_state_source_miss_is_terminal_with_complete_evidence(self) -> None:
        transaction = self.transaction()
        causal_before = transaction.causal_state
        with self.assertRaises(Proto7TerminalStop) as caught:
            self.attempt(_rhs(source_calls={0}), transaction)
        evidence = caught.exception.evidence
        self.assertIsInstance(evidence, AcceptedStateSourceFailureEvidence)
        self.assertEqual(caught.exception.reason, "newton_residual_limit")
        self.assertEqual(evidence.accepted_state_sha256, array_content_sha256(
            _state().u, _state().p, _state().q
        ))
        self.assertFalse(evidence.proposal_constructed)
        self.assertEqual(transaction.state.accepted_stage_count, 0)
        self.assertEqual(transaction.causal_state, causal_before)

    def test_internal_source_only_miss_retries_from_exact_boundary(self) -> None:
        transaction = self.transaction()
        state_before = transaction.state
        causal_before = transaction.causal_state
        preaccept_calls = 0

        def preaccept(_proposal):
            nonlocal preaccept_calls
            preaccept_calls += 1

        result = self.attempt(
            _rhs(source_calls={4}),
            transaction,
            preaccept=preaccept,
        )
        self.assertIsNone(result.accepted)
        self.assertIsInstance(result.retry, ProposalFailureEvidence)
        self.assertTrue(result.retry.retryable_source_only)
        self.assertEqual(result.retry.failed_evaluations[0].stage_name, "rk4_k4")
        self.assertEqual(result.retry.accepted_step_index, 0)
        self.assertEqual(result.retry.accepted_transaction_serial, 0)
        self.assertEqual(transaction.state, state_before)
        self.assertEqual(transaction.causal_state, causal_before)
        self.assertEqual(preaccept_calls, 0)

    def test_candidate_endpoint_source_only_miss_is_now_retryable(self) -> None:
        transaction = self.transaction()
        result = self.attempt(_rhs(source_calls={5}), transaction)
        self.assertIsNone(result.accepted)
        self.assertTrue(result.retry.retryable_source_only)
        self.assertEqual(
            [item.stage_name for item in result.retry.failed_evaluations],
            ["candidate_endpoint"],
        )
        self.assertIsNone(result.retry.selected_terminal_reason)
        self.assertFalse(result.retry.non_source_failure_vetoed_retry)

    def test_every_failed_source_call_is_serialized_across_stage_labels(self) -> None:
        transaction = self.transaction()
        result = self.attempt(_rhs(source_calls={2, 4, 5}), transaction)
        self.assertEqual(
            [item.stage_name for item in result.retry.failed_evaluations],
            ["rk4_k2", "rk4_k4", "candidate_endpoint"],
        )
        self.assertEqual(result.retry.complete_failure_set, ("newton_residual_limit",))
        self.assertTrue(all(
            item.source_residual_infinity == 2.0e-12
            and item.source_residual_maximum == 1.0e-12
            for item in result.retry.failed_evaluations
        ))

    def test_later_non_source_failure_vetoes_earlier_source_retry(self) -> None:
        transaction = self.transaction()
        causal_before = transaction.causal_state
        with self.assertRaises(Proto7TerminalStop) as caught:
            self.attempt(
                _rhs(source_calls={2}, lapse_calls={4}),
                transaction,
            )
        evidence = caught.exception.evidence
        self.assertEqual(caught.exception.reason, "newton_residual_limit")
        self.assertTrue(evidence.non_source_failure_vetoed_retry)
        self.assertFalse(evidence.retryable_source_only)
        self.assertEqual(
            evidence.complete_failure_set,
            ("nonpositive_lapse", "newton_residual_limit"),
        )
        self.assertEqual(
            [item.stage_name for item in evidence.failed_evaluations],
            ["rk4_k2", "rk4_k4"],
        )
        self.assertEqual(transaction.causal_state, causal_before)

    def test_simultaneous_non_source_failure_keeps_inherited_priority(self) -> None:
        transaction = self.transaction()
        with self.assertRaises(Proto7TerminalStop) as caught:
            self.attempt(
                _rhs(source_calls={4}, lapse_calls={4}, kinetic_calls={4}),
                transaction,
            )
        evidence = caught.exception.evidence
        self.assertEqual(caught.exception.reason, "nonpositive_lapse")
        self.assertEqual(evidence.selected_terminal_reason, "nonpositive_lapse")
        self.assertEqual(
            evidence.failed_evaluations[0].failures,
            (
                "nonpositive_lapse",
                "newton_residual_limit",
                "kinetic_condition_limit",
            ),
        )

    def test_failure_evidence_is_canonical_dataclass_payload(self) -> None:
        transaction = self.transaction()
        result = self.attempt(_rhs(source_calls={5}), transaction)
        payload = asdict(result.retry)
        self.assertEqual(len(payload["accepted_state_sha256"]), 64)
        self.assertEqual(payload["accepted_step_index"], 0)
        self.assertEqual(payload["accepted_transaction_serial"], 0)
        self.assertTrue(payload["fields_bitwise_preserved"])
        self.assertFalse(payload["time_advanced"])
        self.assertFalse(payload["accepted_stage_count_advanced"])
        self.assertFalse(payload["causal_debit_advanced"])
        self.assertFalse(payload["preaccept_called"])

    def test_healthy_proposal_previews_external_state_then_commits_once(self) -> None:
        transaction = self.transaction()
        ordering = []

        def preaccept(_proposal):
            ordering.append(transaction.state.accepted_stage_count)
            return "external-preview"

        result = self.attempt(_rhs(), transaction, preaccept=preaccept)
        self.assertIsNone(result.retry)
        self.assertEqual(result.preaccept_payload, "external-preview")
        self.assertEqual(ordering, [0])
        self.assertEqual(transaction.state.accepted_stage_count, 5)
        self.assertEqual(result.accepted.transaction_serial, 5)

    def test_prior_terminal_latch_remains_immutable(self) -> None:
        transaction = self.transaction()
        with self.assertRaises(Proto7TerminalStop):
            self.attempt(_rhs(source_calls={0}), transaction)
        with self.assertRaises(GR0RuntimeStop) as repeated:
            self.attempt(_rhs(), transaction)
        self.assertEqual(repeated.exception.reason, "newton_residual_limit")


if __name__ == "__main__":
    unittest.main()
