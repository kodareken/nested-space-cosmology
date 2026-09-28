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
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement_runtime import (  # noqa: E402
    TDG5RefinementPathStop,
    binary64_cubic_absolute_envelope,
    commit_tdg5_gr0_refinement_pair,
    prepare_tdg5_gr0_refinement_pair,
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


def _rhs(*, failed_time: float | None = None):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        failed = failed_time is not None and (
            np.float64(time).tobytes() == np.float64(failed_time).tobytes()
        )
        return EvolutionRHS(
            state.u,
            state.p,
            state.q,
            {
                "source_residual_infinity": 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )


class TDG5StageCompleteRefinementRuntimeTests(unittest.TestCase):
    def prepare(self, method: str, *, rhs=None, transaction=None):
        return prepare_tdg5_gr0_refinement_pair(
            method=method,
            time=0.0,
            step_size=0.125,
            state=_state(),
            rhs=_rhs() if rhs is None else rhs,
            projector=None,
            transaction=_transaction() if transaction is None else transaction,
            previous_step_index=0,
            previous_transaction_serial=0,
        )

    def test_binary64_extrema_cover_endpoints_and_every_real_root(self) -> None:
        bump = binary64_cubic_absolute_envelope(np.asarray((0.0, 4.0, -4.0, 0.0)))
        self.assertEqual(bump.polynomial_count, 1)
        self.assertEqual(bump.endpoint_candidate_count, 2)
        self.assertEqual(bump.real_interior_root_count, 1)
        self.assertEqual(bump.raw_candidate_maximum, 1.0)
        self.assertGreaterEqual(bump.certified_continuous_upper_bound, 1.0)
        self.assertGreater(bump.outward_arithmetic_debit, 0.0)

        cubic = binary64_cubic_absolute_envelope(
            np.asarray((0.0, 9.0 / 16.0, -1.5, 1.0))
        )
        self.assertEqual(cubic.real_interior_root_count, 2)
        self.assertEqual(cubic.ambiguous_discriminant_count, 0)

    def test_outward_envelope_dominates_dense_independent_controls(self) -> None:
        controls = np.asarray(
            (
                (0.0, 4.0, -4.0, 0.0),
                (0.0, 9.0 / 16.0, -1.5, 1.0),
                (0.25, -0.75, 0.125, 0.5),
                (-1.0, 2.0, -3.0, 1.25),
            ),
            dtype=np.float64,
        )
        evidence = binary64_cubic_absolute_envelope(controls)
        theta = np.linspace(0.0, 1.0, 131_073)
        dense = np.max(
            np.abs(
                controls[:, 0, None]
                + theta
                * (
                    controls[:, 1, None]
                    + theta
                    * (controls[:, 2, None] + theta * controls[:, 3, None])
                )
            )
        )
        self.assertLessEqual(float(dense), evidence.certified_continuous_upper_bound)
        self.assertGreaterEqual(evidence.raw_candidate_maximum, float(dense) - 1.0e-12)

    def test_exact_zero_remains_zero_without_invented_floor(self) -> None:
        evidence = binary64_cubic_absolute_envelope(np.zeros((3, 4)))
        self.assertEqual(evidence.raw_candidate_maximum, 0.0)
        self.assertEqual(evidence.coefficient_construction_debit, 0.0)
        self.assertEqual(evidence.outward_arithmetic_debit, 0.0)
        self.assertEqual(evidence.bernstein_certification_slack, 0.0)
        self.assertEqual(evidence.certified_continuous_upper_bound, 0.0)

    def test_one_bit_coefficient_mutation_changes_public_evidence(self) -> None:
        baseline = np.asarray((0.0, 4.0, -4.0, 0.0))
        mutation = baseline.copy()
        mutation[1] = np.nextafter(mutation[1], np.inf)
        original = binary64_cubic_absolute_envelope(baseline)
        changed = binary64_cubic_absolute_envelope(mutation)
        self.assertNotEqual(asdict(original), asdict(changed))

    def test_derivative_scaling_and_subnormal_debit_cover_extreme_inputs(self) -> None:
        maximum = np.finfo(np.float64).max / 4.0
        large = binary64_cubic_absolute_envelope(
            np.asarray((0.0, 0.0, 0.0, maximum))
        )
        self.assertTrue(np.isfinite(large.certified_continuous_upper_bound))
        self.assertEqual(large.raw_candidate_maximum, maximum)

        subnormal = np.nextafter(0.0, np.inf)
        tiny = binary64_cubic_absolute_envelope(
            np.asarray((0.0, subnormal, 0.0, 0.0))
        )
        self.assertGreater(tiny.outward_arithmetic_debit, 0.0)
        self.assertGreaterEqual(tiny.certified_continuous_upper_bound, subnormal)

    def test_both_methods_prepare_same_grid_complete_state_pairs(self) -> None:
        for method, label, order, denominator, stage_count in (
            (PRIMARY_METHOD, "RK4", 4, 15, 5),
            (COMPARATOR_METHOD, "SSPRK3", 3, 7, 4),
        ):
            with self.subTest(method=label):
                transaction = _transaction()
                monitor_before = transaction.state
                causal_before = transaction.causal_state
                prepared = self.prepare(method, transaction=transaction)
                evidence = prepared.continuous_difference
                self.assertEqual(evidence.method_label, label)
                self.assertEqual(evidence.formal_order, order)
                self.assertEqual(evidence.richardson_denominator, denominator)
                self.assertEqual(evidence.owned_row_count, 4)
                self.assertEqual(evidence.field_count, 6)
                self.assertEqual(
                    evidence.conditional_fine_path_richardson_debit,
                    evidence.certified_continuous_upper_bound / denominator,
                )
                self.assertEqual(transaction.state, monitor_before)
                self.assertEqual(transaction.causal_state, causal_before)
                self.assertEqual(len(prepared.coarse.proposal.stages), stage_count)
                self.assertFalse(evidence.declared_cubic_is_exact_PDE_history)
                self.assertFalse(evidence.trajectory_asymptotic_regime_proved)

    def test_only_two_half_step_path_commits_to_the_real_ledger(self) -> None:
        for method, expected_stages, expected_serial in (
            (PRIMARY_METHOD, 10, 10),
            (COMPARATOR_METHOD, 8, 8),
        ):
            with self.subTest(method=method):
                transaction = _transaction()
                initial = _state()
                prepared = prepare_tdg5_gr0_refinement_pair(
                    method=method,
                    time=0.0,
                    step_size=0.125,
                    state=initial,
                    rhs=_rhs(),
                    projector=None,
                    transaction=transaction,
                    previous_step_index=0,
                    previous_transaction_serial=0,
                )
                committed = commit_tdg5_gr0_refinement_pair(
                    prepared,
                    transaction=transaction,
                    current_time=0.0,
                    current_state=initial,
                    current_step_index=0,
                    current_transaction_serial=0,
                )
                self.assertEqual(committed.committed_path, "fine_two_half_steps")
                self.assertFalse(committed.coarse_path_committed)
                self.assertTrue(committed.fine_half_steps_committed_atomically)
                self.assertEqual(committed.time, 0.125)
                self.assertEqual(committed.step_index, 2)
                self.assertEqual(committed.transaction_serial, expected_serial)
                self.assertEqual(
                    transaction.state.accepted_stage_count, expected_stages
                )
                self.assertEqual(
                    array_content_sha256(
                        committed.state.u, committed.state.p, committed.state.q
                    ),
                    array_content_sha256(
                        prepared.fine_right.accepted.state.u,
                        prepared.fine_right.accepted.state.p,
                        prepared.fine_right.accepted.state.q,
                    ),
                )
                self.assertNotEqual(
                    array_content_sha256(
                        committed.state.u, committed.state.p, committed.state.q
                    ),
                    array_content_sha256(
                        prepared.coarse.accepted.state.u,
                        prepared.coarse.accepted.state.p,
                        prepared.coarse.accepted.state.q,
                    ),
                )

    def test_fine_right_failure_preserves_real_state_and_ledgers(self) -> None:
        transaction = _transaction()
        monitor_before = transaction.state
        causal_before = transaction.causal_state
        initial = _state()
        initial_hash = array_content_sha256(initial.u, initial.p, initial.q)
        with self.assertRaises(TDG5RefinementPathStop) as caught:
            prepare_tdg5_gr0_refinement_pair(
                method=PRIMARY_METHOD,
                time=0.0,
                step_size=0.125,
                state=initial,
                rhs=_rhs(failed_time=0.09375),
                projector=None,
                transaction=transaction,
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        self.assertEqual(caught.exception.path, "fine_right")
        self.assertIsInstance(caught.exception.cause, Proto7TerminalStop)
        self.assertEqual(transaction.state, monitor_before)
        self.assertEqual(transaction.causal_state, causal_before)
        self.assertEqual(initial_hash, array_content_sha256(initial.u, initial.p, initial.q))
        self.assertFalse(caught.exception.real_transaction_advanced)
        self.assertFalse(caught.exception.accepted_state_advanced)

    def test_commit_fails_closed_on_state_or_ledger_drift(self) -> None:
        initial = _state()
        transaction = _transaction()
        prepared = prepare_tdg5_gr0_refinement_pair(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.125,
            state=initial,
            rhs=_rhs(),
            projector=None,
            transaction=transaction,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        transaction.causal_state = CausalBudgetState(previous_speed_upper=1.25)
        drifted = transaction.causal_state
        with self.assertRaisesRegex(RuntimeError, "drifted"):
            commit_tdg5_gr0_refinement_pair(
                prepared,
                transaction=transaction,
                current_time=0.0,
                current_state=initial,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(transaction.causal_state, drifted)
        self.assertEqual(transaction.state.accepted_stage_count, 0)

    def test_invalid_shapes_times_and_coefficients_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            binary64_cubic_absolute_envelope(np.ones(3))
        with self.assertRaises(ValueError):
            binary64_cubic_absolute_envelope(
                np.ones(4), np.asarray((0.0, 0.0, -1.0, 0.0))
            )
        transaction = _transaction()
        transaction.causal_state = CausalBudgetState(
            accepted_time=1.0e16,
            previous_speed_upper=1.0,
        )
        with self.assertRaisesRegex(ValueError, "bitwise divisible"):
            prepare_tdg5_gr0_refinement_pair(
                method=PRIMARY_METHOD,
                time=1.0e16,
                step_size=2.0,
                state=_state(),
                rhs=_rhs(),
                projector=None,
                transaction=transaction,
                previous_step_index=0,
                previous_transaction_serial=0,
            )

    def test_evidence_payload_keeps_every_nonclaim_public(self) -> None:
        prepared = self.prepare(PRIMARY_METHOD)
        payload = asdict(prepared.continuous_difference)
        self.assertFalse(payload["declared_cubic_is_exact_PDE_history"])
        self.assertFalse(payload["trajectory_asymptotic_regime_proved"])
        self.assertFalse(payload["rigorous_global_PDE_error_bound_proved"])
        self.assertIn("bernstein_certification_slack", payload["left_half"])
        self.assertIn("outward_arithmetic_debit", payload["right_half"])


if __name__ == "__main__":
    unittest.main()
