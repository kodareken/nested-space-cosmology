from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration_v7 import (  # noqa: E402
    Proto7RunMember,
    Proto7SourceRetryExhausted,
)
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)


def _flat_state() -> EvolutionState:
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
    source_calls: set[int],
    lapse_calls: set[int] | None = None,
):
    calls = 0
    lapse = set() if lapse_calls is None else set(lapse_calls)

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
                "source_residual_infinity": 2.0e-12 if index in source_calls else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if index in lapse else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


class _Tracers:
    def __init__(self) -> None:
        self.positions = np.array([0.5])
        self.proper_times = np.array([0.0])
        self.event_proper_times = []
        self.event_fields = []
        self.preview_count = 0
        self.commit_count = 0

    def preview_advance(self, **_kwargs):
        self.preview_count += 1
        return np.array([1.0]), np.array([2.0])

    def commit_advance(self, positions, proper_times):
        self.commit_count += 1
        self.positions = positions
        self.proper_times = proper_times

    def append_common_event(self, *_args):
        self.event_proper_times.append(self.proper_times.copy())
        self.event_fields.append(np.zeros((1, 6)))


def _member(operator) -> Proto7RunMember:
    grid = UniformRadialGrid(0.0, 1.0, 9)
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=grid.spacing,
        cfl_maximum=0.125,
    )
    return Proto7RunMember(
        amplitude="5/2",
        method_label="RK4",
        integrator_id=PRIMARY_METHOD,
        spatial_order=4,
        point_count=9,
        input_hash="0" * 64,
        initial=SimpleNamespace(grid=grid),
        state=_flat_state(),
        operator=operator,
        projector=None,
        transaction=transaction,
        tracers=_Tracers(),
    )


class FGCGR0CampaignRunnerV7Tests(unittest.TestCase):
    def test_candidate_endpoint_retry_is_serialized_then_step_accepts(self) -> None:
        member = _member(_rhs(source_calls={5}))
        records = []
        member.advance_to(
            0.01,
            retry_factor=0.5,
            maximum_CFL_retries=32,
            maximum_source_retries=32,
            minimum_step_size=2.0**-30,
            rejected_trial_sink=records.append,
        )
        self.assertEqual(member.time, 0.01)
        self.assertEqual(member.step_index, 2)
        self.assertEqual(member.source_retry_count, 1)
        self.assertEqual(member.tracers.preview_count, 2)
        self.assertEqual(member.tracers.commit_count, 2)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["event_type"], "rejected_unaccepted_source_only_proposal")
        self.assertEqual(record["failed_evaluations"][0]["stage_name"], "candidate_endpoint")
        self.assertEqual(record["retry_count_for_accepted_step"], 1)
        self.assertTrue(record["external_transaction_state_preserved"])
        self.assertTrue(record["complete_failure_evidence_serialized_before_retry"])

    def test_retry_exhaustion_retains_complete_last_proposal(self) -> None:
        member = _member(_rhs(source_calls={5, 11}))
        state_before = member.state
        causal_before = member.transaction.causal_state
        records = []
        with self.assertRaises(Proto7SourceRetryExhausted) as caught:
            member.advance_to(
                0.01,
                retry_factor=0.5,
                maximum_CFL_retries=32,
                maximum_source_retries=1,
                minimum_step_size=2.0**-30,
                rejected_trial_sink=records.append,
            )
        evidence = caught.exception.evidence
        self.assertEqual(evidence["exhaustion_kind"], "maximum_source_retries")
        self.assertEqual(len(records), 2)
        self.assertEqual(
            evidence["last_rejected_proposal"]["failed_evaluations"][0]["stage_name"],
            "candidate_endpoint",
        )
        self.assertIs(member.state, state_before)
        self.assertEqual(member.time, 0.0)
        self.assertEqual(member.step_index, 0)
        self.assertEqual(member.transaction.causal_state, causal_before)

    def test_non_source_veto_propagates_complete_terminal_evidence(self) -> None:
        member = _member(_rhs(source_calls={5}, lapse_calls={5}))
        with self.assertRaises(Exception) as caught:
            member.advance_to(
                0.01,
                retry_factor=0.5,
                maximum_CFL_retries=32,
                maximum_source_retries=32,
                minimum_step_size=2.0**-30,
                rejected_trial_sink=lambda _record: self.fail("terminal stop was serialized as retry"),
            )
        self.assertEqual(type(caught.exception).__name__, "Proto7TerminalStop")
        self.assertEqual(caught.exception.reason, "nonpositive_lapse")
        self.assertTrue(caught.exception.evidence.non_source_failure_vetoed_retry)

    def test_snapshot_restore_recovers_external_member_boundary(self) -> None:
        member = _member(_rhs(source_calls=set()))
        snapshot = member.snapshot()
        member.time = 1.0
        member.step_index = 7
        member.source_retry_count = 4
        member.tracers.positions[:] = 9.0
        member.restore(snapshot)
        self.assertEqual(member.time, 0.0)
        self.assertEqual(member.step_index, 0)
        self.assertEqual(member.source_retry_count, 0)
        np.testing.assert_array_equal(member.tracers.positions, np.array([0.5]))


if __name__ == "__main__":
    unittest.main()
