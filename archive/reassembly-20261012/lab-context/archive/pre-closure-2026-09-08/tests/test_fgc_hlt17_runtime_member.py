"""Focused in-memory HLT17 carrier controls; no campaign or store access."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    HLT17RuntimeMemberError,
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
    PINNED_PRIVATE_DEPENDENCIES,
    RUNTIME_PROTOCOL,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto14_runtime import Proto14RunMember  # noqa: E402
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    inherited_tdg6_ledger,
    seed_imp1_ledger,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)
from recursive_horizons.fgc.evolution.tdg11_c1r1_runtime import (  # noqa: E402
    C1R1TemporalRetryRequired,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6TemporalLedger,
)
from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState  # noqa: E402


MEMBER_SOURCE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_runtime_member.py"


@dataclass
class ProvidedGrid:
    coordinates: np.ndarray

    @property
    def point_count(self) -> int:
        return int(self.coordinates.size)

    @property
    def minimum(self) -> float:
        return float(self.coordinates[0])

    @property
    def maximum(self) -> float:
        return float(self.coordinates[-1])

    @property
    def spacing(self) -> float:
        return (self.maximum - self.minimum) / (self.point_count - 1)


@dataclass
class ProvidedInitial:
    grid: ProvidedGrid


def _identity(method: str, point_count: int) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_SYNTHETIC,
        method_label=label,
        point_count=point_count,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id="HLT17-SYNTHETIC-TEST",
        amplitude="synthetic",
    )


def _make_member(
    method: str = PRIMARY_METHOD,
    *,
    rhs=synthetic_exponential_rhs,
    signed_zeros: bool = False,
) -> HLT17RuntimeMember:
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    if signed_zeros:
        u = np.array(state.u, dtype="<f8", copy=True)
        p = np.array(state.p, dtype="<f8", copy=True)
        q = np.array(state.q, dtype="<f8", copy=True)
        u[0, 5] = -0.0
        p[2, 0] = -0.0
        q[1, 1] = -0.0
        state = EvolutionState(u, p, q)
        ledger = seed_imp1_ledger(
            inherited_tdg6_ledger(ledger),
            method=method,
            state_sha256=array_content_sha256(state.u, state.p, state.q),
            step_index=0,
            transaction_serial=0,
            origin_receipt_sha256=ledger.origin_receipt_sha256,
        )
        tracers.proper_times = np.array(tracers.proper_times, dtype="<f8", copy=True)
        tracers.proper_times[0] = -0.0
        transaction.causal_state = CausalBudgetState(
            accepted_time=transaction.causal_state.accepted_time,
            accumulated_characteristic_distance=-0.0,
            previous_speed_upper=transaction.causal_state.previous_speed_upper,
        )
    return HLT17RuntimeMember(
        identity=_identity(method, int(state.shape[0])),
        initial=ProvidedInitial(ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))),
        operator=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        state=state,
        temporal_ledger=ledger,
        source_operator_closure=HLT17ExternalReference(
            kind="source_operator_closure", sha256="a" * 64
        ),
        origin_reference=HLT17ExternalReference(
            kind="origin_receipt", sha256=ledger.origin_receipt_sha256
        ),
        time=float(ledger.last_accepted_time),
        step_index=int(ledger.current_step_index),
        transaction_serial=int(ledger.current_transaction_serial),
    )


def _quarter_impulse(time, current):
    du = np.zeros_like(current.u)
    if time == START + 1.0 / 128:
        du[:, 4] = 1.0
    return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), synthetic_diagnostics())


def _bits(member: HLT17RuntimeMember) -> dict[str, bytes]:
    return {
        "u": np.array(member.state.u, dtype="<f8", copy=True).tobytes(),
        "p": np.array(member.state.p, dtype="<f8", copy=True).tobytes(),
        "q": np.array(member.state.q, dtype="<f8", copy=True).tobytes(),
        "labels": np.asarray(member.tracers.labels).tobytes(),
        "positions": np.asarray(member.tracers.positions).tobytes(),
        "proper": np.asarray(member.tracers.proper_times).tobytes(),
        "events_t": b"".join(np.asarray(row).tobytes() for row in member.tracers.event_proper_times),
        "events_f": b"".join(np.asarray(row).tobytes() for row in member.tracers.event_fields),
        "ledger": member.temporal_ledger.inherited_snapshot.encode("ascii"),
        "cutoff": np.float64(member.tracers.cutoff).tobytes(),
        "outer": np.float64(member.tracers.outer_radius).tobytes(),
        "time": np.float64(member.time).tobytes(),
        "causal": np.float64(member.transaction.causal_state.accumulated_characteristic_distance).tobytes(),
    }


class HLT17RuntimeMemberTests(unittest.TestCase):
    def test_carrier_is_not_a_historical_subclass_or_relabel(self) -> None:
        self.assertFalse(issubclass(HLT17RuntimeMember, Proto14RunMember))
        source = MEMBER_SOURCE.read_text(encoding="utf-8")
        self.assertNotIn(".advance_to(", source)
        self.assertNotIn("hlt16_member_codec", source)
        self.assertNotIn("from .proto14_runtime", source)
        for name in PINNED_PRIVATE_DEPENDENCIES:
            self.assertIn(name.rsplit(".", 1)[-1], source)
        member = _make_member()
        self.assertEqual(member.key, "RK4-9")
        self.assertEqual(member.identity.scope, IDENTITY_SCOPE_SYNTHETIC)
        self.assertFalse(member.temporal_ledger.historical_origin_authenticated)
        self.assertFalse(member.origin_reference.authenticated)
        self.assertFalse(member.source_operator_closure.authenticated)
        self.assertEqual(RUNTIME_PROTOCOL, "FGC-2-SF1-PROTO19")

    def test_tdg6_ledger_and_false_live_keys_are_rejected(self) -> None:
        with self.assertRaises(HLT17RuntimeMemberError):
            HLT17MemberIdentity(
                scope=IDENTITY_SCOPE_LIVE,
                method_label="RK4",
                point_count=9,
                integrator_id=PRIMARY_METHOD,
                spatial_order=4,
                campaign_id="HLT17-SYNTHETIC-TEST",
                amplitude="3",
            )
        with self.assertRaises(HLT17RuntimeMemberError):
            HLT17MemberIdentity(
                scope=IDENTITY_SCOPE_SYNTHETIC,
                method_label="RK4",
                point_count=2049,
                integrator_id=PRIMARY_METHOD,
                spatial_order=4,
                campaign_id="HLT17-SYNTHETIC-TEST",
                amplitude="synthetic",
            )
        state, transaction, tracers, ledger, coordinates = make_synthetic_fixture()
        with self.assertRaises(TypeError):
            HLT17RuntimeMember(
                identity=_identity(PRIMARY_METHOD, 9),
                initial=ProvidedInitial(ProvidedGrid(coordinates)),
                operator=synthetic_exponential_rhs,
                projector=None,
                transaction=transaction,
                tracers=tracers,
                state=state,
                temporal_ledger=TDG6TemporalLedger.zero(initial_time=START),
                source_operator_closure=HLT17ExternalReference(
                    kind="source_operator_closure", sha256="a" * 64
                ),
                origin_reference=HLT17ExternalReference(
                    kind="origin_receipt", sha256="1" * 64
                ),
                time=START,
            )

    def test_both_methods_prepare_commit_restore_and_next_step(self) -> None:
        for method, serial in ((PRIMARY_METHOD, 20), (COMPARATOR_METHOD, 16)):
            with self.subTest(method=method):
                member = _make_member(method)
                inherited = member.temporal_ledger.inherited_snapshot
                inherited_debit = member.inherited_tdg6_sibling().accumulated_debit_vector
                operator = member.operator
                projector = member.projector
                prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
                self.assertTrue(prepared.assessment.admission_passed)
                committed = member.commit_prepared(prepared)
                self.assertEqual(committed.step_index, 4)
                self.assertEqual(committed.transaction_serial, serial)
                self.assertEqual(member.step_index, 4)
                self.assertIs(member.operator, operator)
                self.assertIs(member.projector, projector)
                self.assertEqual(member.temporal_ledger.inherited_snapshot, inherited)
                self.assertEqual(
                    member.inherited_tdg6_sibling().accumulated_debit_vector,
                    inherited_debit,
                )
                self.assertIsNotNone(member.executed_plan)
                before = _bits(member)
                capture = member.capture()
                member.restore(capture)
                self.assertEqual(_bits(member), before)
                nxt = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
                self.assertEqual(nxt.initial_time, member.time)
                self.assertTrue(nxt.assessment.admission_passed)

    def test_retry_ledger_preserves_physical_state_and_inherited_sibling(self) -> None:
        member = _make_member(rhs=_quarter_impulse)
        physical = member.physical_state_sha256()
        inherited = member.temporal_ledger.inherited_snapshot
        debit = member.temporal_ledger.accumulated_debit_vector
        labels = np.asarray(member.tracers.labels).tobytes()
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        self.assertFalse(prepared.assessment.admission_passed)
        written: list[object] = []
        with self.assertRaises(C1R1TemporalRetryRequired) as caught:
            member.require_admission(prepared, durable_rejection_sink=written.append)
        self.assertEqual(len(written), 1)
        member.adopt_retry_ledger(caught.exception)
        self.assertEqual(member.physical_state_sha256(), physical)
        self.assertEqual(member.temporal_ledger.inherited_snapshot, inherited)
        self.assertEqual(member.temporal_ledger.accumulated_debit_vector, debit)
        self.assertEqual(member.temporal_ledger.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(np.asarray(member.tracers.labels).tobytes(), labels)
        self.assertEqual(member.executed_plan, None)

    def test_signed_zeros_and_tracer_history_survive_capture_restore(self) -> None:
        member = _make_member(signed_zeros=True)
        self.assertEqual(np.float64(member.state.u[0, 5]).tobytes(), np.float64(-0.0).tobytes())
        self.assertEqual(np.float64(member.state.q[1, 1]).tobytes(), np.float64(-0.0).tobytes())
        self.assertEqual(
            np.float64(member.tracers.proper_times[0]).tobytes(),
            np.float64(-0.0).tobytes(),
        )
        self.assertEqual(
            np.float64(member.transaction.causal_state.accumulated_characteristic_distance).tobytes(),
            np.float64(-0.0).tobytes(),
        )
        original = _bits(member)
        capture = member.capture()
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        member.commit_prepared(prepared)
        self.assertNotEqual(_bits(member)["u"], original["u"])
        member.restore(capture)
        self.assertEqual(_bits(member), original)

    def test_stale_coordinates_refuse_before_mutation(self) -> None:
        member = _make_member()
        before = _bits(member)
        capture = member.capture()
        member.initial.grid.coordinates = np.linspace(0.0, 2.0, member.point_count)
        with self.assertRaises(HLT17RuntimeMemberError):
            member.restore(capture)
        member.initial.grid.coordinates = np.array(member.coordinates, dtype="<f8", copy=True)
        self.assertEqual(_bits(member), before)

    def test_late_restore_and_adoption_failures_roll_back_completely(self) -> None:
        member = _make_member()
        origin = member.capture()
        origin_bits = _bits(member)

        def boom_restore(name: str) -> None:
            if name == "after_restore_mutation":
                raise RuntimeError("injected late restore failure")

        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        member.commit_prepared(prepared)
        advanced = member.capture()
        member.restore(origin)
        with self.assertRaisesRegex(RuntimeError, "injected late restore failure"):
            member.restore(advanced, fault_hook=boom_restore)
        self.assertEqual(_bits(member), origin_bits)
        self.assertEqual(member.step_index, 0)
        self.assertIsNone(member.executed_plan)

        def boom_adopt(name: str) -> None:
            if name == "after_member_install":
                raise RuntimeError("injected late adoption failure")

        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        with self.assertRaisesRegex(RuntimeError, "injected late adoption failure"):
            member.commit_prepared(prepared, fault_hook=boom_adopt)
        self.assertEqual(_bits(member), origin_bits)
        self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 0)
        self.assertEqual(member.transaction.state.accepted_stage_count, 0)

    def test_adopt_committed_does_not_call_the_kernel_and_rolls_back(self) -> None:
        member = _make_member()
        origin = member.capture()
        origin_bits = _bits(member)
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        committed = member.commit_prepared(prepared)
        member.restore(origin)
        from recursive_horizons.fgc.evolution import tdg11_c1r1_runtime as c1r1
        from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as imp1

        with patch.object(
            c1r1, "commit_c1r1_runtime", side_effect=AssertionError("second kernel commit")
        ), patch.object(
            imp1, "commit_imp1_runtime", side_effect=AssertionError("imp1 kernel commit")
        ):
            with self.assertRaises(Exception):
                member.adopt_committed(committed, executed_plan=prepared.plan)
        self.assertEqual(_bits(member), origin_bits)
        self.assertEqual(member.step_index, 0)
        self.assertIsNone(member.executed_plan)


if __name__ == "__main__":
    unittest.main()
