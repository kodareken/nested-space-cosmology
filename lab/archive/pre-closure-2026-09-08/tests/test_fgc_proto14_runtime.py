from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration import NormalFlowTracers  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution import proto14_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg6_temporal_admission_runtime as tdg6,
)


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(time: float, state: EvolutionState) -> EvolutionRHS:
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
            "minimum_lapse": 1.0,
            "minimum_radial_metric": 1.0,
            "minimum_areal_radius_away_from_center": 0.125,
        },
    )


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )


def _member() -> runtime.Proto14RunMember:
    state = _state()
    return runtime.Proto14RunMember(
        amplitude="3",
        method_label="RK4",
        integrator_id=PRIMARY_METHOD,
        spatial_order=4,
        point_count=9,
        input_hash="0" * 64,
        initial=SimpleNamespace(grid=SimpleNamespace(spacing=1.0, coordinates=np.linspace(0.0, 1.0, 9))),
        state=state,
        operator=_rhs,
        projector=None,
        transaction=_transaction(),
        tracers=NormalFlowTracers.create(
            minimum=0.25,
            maximum=0.75,
            spacing=0.25,
            state=state,
            coordinates=np.linspace(0.0, 1.0, 9),
            cutoff=1.0,
            outer_radius=128.0,
        ),
    )


def _failed_admission(method: str):
    def interval(value: float, count: int):
        values = np.tile(np.asarray((value, 0.0, 0.0, 0.0)), (count, 1))
        return tdg6._magnitude_interval(
            values,
            np.zeros_like(values),
            subinterval_count=count,
            owned_row_count=1,
        )

    outer = {name: interval(16.0, 2) for name in TDG6_COMPLETE_STATE_CHANNELS}
    fine = {name: interval(1.0, 4) for name in TDG6_COMPLETE_STATE_CHANNELS}
    outer[TDG6_COMPLETE_STATE_CHANNELS[0]] = interval(1.0e-12, 2)
    fine[TDG6_COMPLETE_STATE_CHANNELS[0]] = interval(9.0e-13, 4)
    return tdg6.classify_tdg6_runtime_intervals(
        method=method, outer_intervals=outer, finest_intervals=fine
    )


class Proto14RuntimeTests(unittest.TestCase):
    def _advance(self, member: runtime.Proto14RunMember, target: float) -> None:
        member.advance_to(
            target,
            retry_factor=0.5,
            maximum_CFL_retries=4,
            maximum_source_retries=4,
            minimum_step_size=2.0 ** -30,
            rejected_trial_sink=lambda record: None,
            temporal_rejection_sink=lambda record: None,
        )

    def test_fine_only_admission_advances_the_inherited_member(self) -> None:
        member = _member()
        self._advance(member, 0.125)
        self.assertEqual(member.time, 0.125)
        self.assertEqual(member.step_index, 4)
        self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 1)
        self.assertEqual(member.temporal_ledger.cumulative_temporal_retry_count, 0)
        self.assertEqual(member.transaction.state.accepted_stage_count, 20)

    def test_snapshot_round_trip_preserves_the_temporal_ledger(self) -> None:
        member = _member()
        self._advance(member, 0.125)
        snapshot = member.snapshot()
        self._advance(member, 0.25)
        member.restore(snapshot)
        self.assertEqual(member.time, 0.125)
        self.assertEqual(member.temporal_ledger, snapshot["temporal_ledger"])
        self.assertEqual(member.step_index, 4)

    def test_accepted_boundary_check_includes_complete_retry_and_tracer_history(self) -> None:
        member = _member()
        snapshot = member.snapshot()
        self.assertTrue(member._accepted_boundary_preserved(snapshot))
        member.CFL_retry_count += 1
        self.assertFalse(member._accepted_boundary_preserved(snapshot))
        member.CFL_retry_count -= 1
        member.tracers.event_fields[0][0, 0] += 1.0
        self.assertFalse(member._accepted_boundary_preserved(snapshot))

    def test_checkpoint_extension_refuses_a_time_mismatch(self) -> None:
        member = _member()
        self._advance(member, 0.125)
        metadata, debit = runtime.proto14_checkpoint_extension(member)
        restored = _member()
        with self.assertRaisesRegex(ValueError, "time differs"):
            runtime.restore_proto14_checkpoint_extension(restored, metadata, debit)
        restored.time = 0.125
        restored.temporal_ledger = tdg6.TDG6TemporalLedger.zero(initial_time=0.125)
        runtime.restore_proto14_checkpoint_extension(restored, metadata, debit)
        self.assertEqual(restored.temporal_ledger, member.temporal_ledger)

    def test_temporal_rejection_is_durable_before_ledger_adoption_and_retry(self) -> None:
        member = _member()
        original_prepare = runtime.prepare_tdg6_gr0_compositor
        calls = 0
        records: list[dict[str, object]] = []

        def prepare_once_failed(**kwargs):
            nonlocal calls
            calls += 1
            prepared = original_prepare(**kwargs)
            if calls == 1:
                return replace(
                    prepared,
                    continuous_admission=_failed_admission(member.integrator_id),
                )
            return prepared

        with patch.object(runtime, "prepare_tdg6_gr0_compositor", side_effect=prepare_once_failed):
            member.advance_to(
                0.125,
                retry_factor=0.5,
                maximum_CFL_retries=4,
                maximum_source_retries=4,
                minimum_step_size=2.0 ** -30,
                rejected_trial_sink=lambda record: None,
                temporal_rejection_sink=lambda record: records.append(dict(record)),
            )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["event_type"], "rejected_TDG6_temporal_admission")
        self.assertEqual(member.temporal_ledger.cumulative_temporal_retry_count, 1)
        self.assertEqual(len(member.temporal_ledger.serialized_temporal_rejections), 1)
        self.assertGreater(member.temporal_ledger.accepted_macro_step_count, 0)

    def test_failed_temporal_sink_leaves_the_ledger_and_accepted_boundary_unchanged(self) -> None:
        member = _member()
        before = member.snapshot()
        original_prepare = runtime.prepare_tdg6_gr0_compositor

        def failed_prepare(**kwargs):
            prepared = original_prepare(**kwargs)
            return replace(
                prepared,
                continuous_admission=_failed_admission(member.integrator_id),
            )

        with patch.object(runtime, "prepare_tdg6_gr0_compositor", side_effect=failed_prepare):
            with self.assertRaisesRegex(RuntimeError, "sink failed"):
                member.advance_to(
                    0.125,
                    retry_factor=0.5,
                    maximum_CFL_retries=4,
                    maximum_source_retries=4,
                    minimum_step_size=2.0 ** -30,
                    rejected_trial_sink=lambda record: None,
                    temporal_rejection_sink=lambda record: (_ for _ in ()).throw(RuntimeError("sink failed")),
                )
        self.assertTrue(member._accepted_boundary_preserved(before))

    def test_wrapped_cfl_refusal_remains_an_inherited_cfl_retry(self) -> None:
        member = _member()
        original_prepare = runtime.prepare_tdg6_gr0_compositor
        calls = 0

        def cfl_then_prepare(**kwargs):
            nonlocal calls
            calls += 1
            if calls == 1:
                raise tdg6.TDG6RefinementPathStop(
                    path="outer_0",
                    cause=CFLRetryRequired(observed_ratio=2.0, maximum_ratio=1.0),
                )
            return original_prepare(**kwargs)

        with patch.object(runtime, "prepare_tdg6_gr0_compositor", side_effect=cfl_then_prepare):
            self._advance(member, 0.125)
        self.assertEqual(member.CFL_retry_count, 1)
        self.assertEqual(member.source_retry_count, 0)
        self.assertEqual(member.temporal_ledger.cumulative_temporal_retry_count, 0)


if __name__ == "__main__":
    unittest.main()
