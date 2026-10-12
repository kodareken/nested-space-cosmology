"""Joint in-memory HLT17 checkpoint controls; no store, origin, or campaign I/O."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes, load_canonical_json  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_c1r1_runtime as c1r1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as imp1  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    prepare_hlt17_next,
    replay_hlt17_rejection,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17_CFL_RETRY_CAP,
    HLT17_MEMBER_KEYS,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SOURCE_RETRY_CAP,
    HLT17_SYNTHETIC_MEMBER_KEYS,
    cursor_sha256,
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    ACCEPTED_FINE,
    FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM,
    HLT17CheckpointError,
    HLT17GenerationBundle,
    HLT17InMemoryCheckpoint,
    IDENTITY_LABELS,
    PARKED_AT_EVENT_TARGET,
    SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN,
    origin_predecessor_identity,
    synthetic_origin_predecessor_identity,
)
from recursive_horizons.fgc.evolution.hlt17_member_codec import (  # noqa: E402
    HLT17MemberCodecError,
    HLT17MemberEncoding,
    ORIGIN_PROTOCOL,
    RUNTIME_PROTOCOL,
    SNAPSHOT_ROLE_ACCEPTED,
)
from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
    TDG6_MINIMUM_OBSERVED_ORDER,
)
from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (  # noqa: E402
    plan_forward_proto14_subdivision,
)


CHECKPOINT_SOURCE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_member_checkpoint.py"
REQUESTED_CAP = 1.0 / 32.0


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


class ControllableRHS:
    def __init__(self) -> None:
        self.mode = "ok"
        self.calls = 0
        self.impulse_from = START

    def hlt17_synthetic_callable_binding(self) -> dict[str, object]:
        code = type(self).__call__.__code__
        return {
            "kind": "synthetic_failure_injection_double",
            "type": f"{type(self).__module__}.{type(self).__qualname__}",
            "call_code_sha256": sha256(code.co_code).hexdigest(),
            "mutable_injection_mode_excluded": True,
            "physical_source_authenticated": False,
        }

    def __call__(self, time: float, current) -> EvolutionRHS:
        self.calls += 1
        diagnostics = synthetic_diagnostics()
        zeros = np.zeros_like(current.u)
        if self.mode == "source":
            if time > self.impulse_from:
                diagnostics["source_residual_infinity"] = 1.0
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "cfl":
            if time > self.impulse_from:
                diagnostics["coordinate_speed_upper"] = 1.0e9
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "impulse":
            du = np.zeros_like(current.u)
            if time == self.impulse_from + 1.0 / 128.0:
                du[:, 4] = 1.0
            return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), diagnostics)
        return synthetic_exponential_rhs(time, current)


def _identity(
    method: str, point_count: int, *, campaign_id: str = "HLT17-SYNTHETIC-TEST"
) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_SYNTHETIC,
        method_label=label,
        point_count=point_count,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id=campaign_id,
        amplitude="synthetic",
    )


def _make_member(
    method: str = PRIMARY_METHOD,
    *,
    rhs=synthetic_exponential_rhs,
    closure: str = "a" * 64,
    campaign_id: str = "HLT17-SYNTHETIC-TEST",
) -> HLT17RuntimeMember:
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    return HLT17RuntimeMember(
        identity=_identity(method, int(state.shape[0]), campaign_id=campaign_id),
        initial=ProvidedInitial(ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))),
        operator=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        state=state,
        temporal_ledger=ledger,
        source_operator_closure=HLT17ExternalReference(
            kind="source_operator_closure", sha256=closure
        ),
        origin_reference=HLT17ExternalReference(
            kind="origin_receipt", sha256=ledger.origin_receipt_sha256
        ),
        time=float(ledger.last_accepted_time),
        step_index=int(ledger.current_step_index),
        transaction_serial=int(ledger.current_transaction_serial),
    )


def _bits(member: HLT17RuntimeMember) -> dict[str, bytes]:
    return {
        "u": np.array(member.state.u, dtype="<f8", copy=True).tobytes(),
        "p": np.array(member.state.p, dtype="<f8", copy=True).tobytes(),
        "q": np.array(member.state.q, dtype="<f8", copy=True).tobytes(),
        "labels": np.asarray(member.tracers.labels).tobytes(),
        "positions": np.asarray(member.tracers.positions).tobytes(),
        "proper": np.asarray(member.tracers.proper_times).tobytes(),
        "events_t": b"".join(
            np.asarray(row).tobytes() for row in member.tracers.event_proper_times
        ),
        "events_f": b"".join(
            np.asarray(row).tobytes() for row in member.tracers.event_fields
        ),
        "ledger": member.temporal_ledger.inherited_snapshot.encode("ascii"),
        "time": np.float64(member.time).tobytes(),
        "plan": (
            b""
            if member.executed_plan is None
            else member.executed_plan.macro_width.hex().encode("ascii")
        ),
        "source": str(member.source_retry_count).encode("ascii"),
        "cfl": str(member.CFL_retry_count).encode("ascii"),
        "serial": str(member.transaction_serial).encode("ascii"),
    }


class HLT17MemberCheckpointTests(unittest.TestCase):
    def test_bundle_rejects_internally_valid_cursor_with_changed_accepted_debit(self) -> None:
        checkpoint = HLT17InMemoryCheckpoint.seed(_make_member(), event_target=TARGET)
        checkpoint.advance_once(requested_cap=1.0 / 32)
        cursor = checkpoint.cursor
        changed_ledger = replace(
            cursor.ledger,
            accumulated_debit_vector=tuple(value + 1 for value in cursor.ledger.accumulated_debit_vector),
        )
        changed_cursor = replace(cursor, ledger=changed_ledger)
        raw_cursor = canonical_json_bytes(encode_hlt17_cursor(changed_cursor))
        with self.assertRaisesRegex(HLT17CheckpointError, "debit"):
            HLT17GenerationBundle(checkpoint.accepted_encoding.descriptor,
                                  checkpoint.accepted_encoding.payload, raw_cursor)

    def test_encoding_metadata_projection_does_not_mutate_its_byte_source(self) -> None:
        checkpoint = HLT17InMemoryCheckpoint.seed(_make_member(), event_target=TARGET)
        encoding = checkpoint.accepted_encoding
        projection = encoding.metadata
        projection["runtime_identity"]["accepted_generation"] = 123
        self.assertEqual(encoding.metadata["runtime_identity"]["accepted_generation"], 0)
        self.assertEqual(checkpoint.generation_bundle().accepted_generation, 0)

    def test_layer_is_in_memory_and_returns_publication_seams(self) -> None:
        source = CHECKPOINT_SOURCE.read_text(encoding="utf-8")
        self.assertNotIn("hlt16_state_store", source)
        self.assertNotIn("publish_exclusive", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertIn("does not write a store", source)
        self.assertEqual(
            IDENTITY_LABELS,
            (
                "physical_u_p_q",
                "accepted_descriptor",
                "full_kernel_cursor",
                "kernel_transition_digest",
                "future_store_journal",
            ),
        )
        self.assertEqual(HLT17_SYNTHETIC_MEMBER_KEYS, ("RK4-9", "SSPRK3-9"))
        self.assertEqual(HLT17_MEMBER_KEYS[-1], "SSPRK3-16385")
        self.assertEqual(HLT17_PRODUCTION_MEMBERS[-1].fallback_d12, 262080)
        self.assertEqual(HLT17_SOURCE_RETRY_CAP, 32)
        self.assertEqual(HLT17_CFL_RETRY_CAP, 32)
        self.assertEqual(len(IMP1_COMPLETE_STATE_CHANNELS), 18)
        self.assertEqual(TDG6_MINIMUM_OBSERVED_ORDER.numerator, 3)
        self.assertEqual(TDG6_MINIMUM_OBSERVED_ORDER.denominator, 2)
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        bundle = checkpoint.generation_bundle()
        self.assertIsNone(bundle.store_journal_parent_sha256)
        self.assertIsNone(bundle.store_journal_tip_sha256)
        self.assertIsNone(bundle.store_checkpoint_revision)
        self.assertFalse(bundle.durable_store_bound)
        self.assertFalse(bundle.campaign_execution_authorized)
        self.assertFalse(bundle.historical_origin_authenticated)
        self.assertFalse(bundle.publication_authorized)
        self.assertEqual(bundle.atomic_publication_seam, FUTURE_ATOMIC_GENERATION_BUNDLE_SEAM)
        self.assertNotEqual(bundle.kernel_transition_digest, bundle.cursor_sha256)
        self.assertNotEqual(bundle.kernel_transition_digest, bundle.descriptor_sha256)
        self.assertNotEqual(
            bundle.predecessor_cursor_identity_sha256, bundle.predecessor_cursor_sha256
        )
        seams = checkpoint.remaining_seams()
        self.assertIn("content_addressed_store", seams)
        self.assertIn("authenticated_origin", seams)
        self.assertIn("atomic_generation_bundle", seams)

    def test_initial_encoding_precedes_cursor_seed_with_actual_descriptor(self) -> None:
        member = _make_member()
        origin = origin_predecessor_identity(member)
        self.assertEqual(
            origin.cursor_identity["protocol_artifact_id"], ORIGIN_PROTOCOL
        )
        checkpoint = HLT17InMemoryCheckpoint.seed(
            member, event_target=TARGET, origin_predecessor=origin
        )
        self.assertEqual(
            checkpoint.cursor.descriptor_sha256,
            checkpoint.accepted_encoding.descriptor_sha256,
        )
        self.assertEqual(checkpoint.accepted_generation, 0)
        self.assertEqual(checkpoint.cursor.source_total, 0)
        metadata = load_canonical_json(checkpoint.accepted_encoding.descriptor)
        self.assertEqual(
            metadata["predecessor_cursor_identity"]["protocol_artifact_id"],
            ORIGIN_PROTOCOL,
        )
        inherited = _make_member()
        inherited.source_retry_count = 3
        inherited.CFL_retry_count = 1
        inherited.__post_init__()
        seeded = HLT17InMemoryCheckpoint.seed(inherited, event_target=TARGET)
        self.assertEqual(seeded.cursor.source_total, 3)
        self.assertEqual(seeded.cursor.cfl_total, 1)
        self.assertEqual(seeded.cursor.source_current, 0)
        self.assertEqual(seeded.member.source_retry_count, 3)

    def test_both_methods_prepare_admit_fine_encode_restore_next_step(self) -> None:
        for method, serial in ((PRIMARY_METHOD, 20), (COMPARATOR_METHOD, 16)):
            with self.subTest(method=method):
                member = _make_member(method)
                checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
                old_cursor = checkpoint.cursor
                old_full = cursor_sha256(old_cursor)
                commits = {"count": 0}
                real_commit = c1r1.commit_c1r1_runtime

                def counted(*args, **kwargs):
                    commits["count"] += 1
                    return real_commit(*args, **kwargs)

                with patch.object(c1r1, "commit_c1r1_runtime", side_effect=counted):
                    result = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
                self.assertEqual(result.disposition, ACCEPTED_FINE)
                self.assertEqual(commits["count"], 1)
                self.assertEqual(member.transaction_serial, serial)
                self.assertEqual(checkpoint.accepted_generation, 1)
                self.assertIsNotNone(checkpoint.last_accepted_executed_plan)
                encoding = result.encoding
                metadata = load_canonical_json(encoding.descriptor)
                self.assertEqual(metadata["snapshot_role"], SNAPSHOT_ROLE_ACCEPTED)
                self.assertEqual(
                    metadata["predecessor_cursor_identity"]["protocol_artifact_id"],
                    RUNTIME_PROTOCOL,
                )
                self.assertEqual(metadata["predecessor_cursor_sha256"], old_full)
                self.assertNotEqual(
                    metadata["predecessor_cursor_identity_sha256"],
                    metadata["predecessor_cursor_sha256"],
                )
                self.assertEqual(
                    checkpoint.cursor.descriptor_sha256, encoding.descriptor_sha256
                )
                self.assertNotEqual(
                    checkpoint.cursor.descriptor_sha256, old_cursor.descriptor_sha256
                )
                new_full = cursor_sha256(checkpoint.cursor)
                self.assertNotIn(new_full, encoding.descriptor.decode("ascii"))
                self.assertNotEqual(new_full, encoding.descriptor_sha256)
                self.assertNotEqual(
                    checkpoint.cursor.kernel_transition_digest, new_full
                )
                origin_member = _make_member(method)
                origin_checkpoint = HLT17InMemoryCheckpoint.seed(
                    origin_member, event_target=TARGET
                )
                origin_checkpoint.restore(
                    encoding, encode_hlt17_cursor(checkpoint.cursor)
                )
                self.assertEqual(origin_member.transaction_serial, serial)
                self.assertEqual(
                    origin_member.executed_plan, checkpoint.last_accepted_executed_plan
                )
                nxt = origin_checkpoint.prepare(requested_cap=REQUESTED_CAP)
                self.assertEqual(nxt.disposition, "prepared_fresh")
                self.assertTrue(nxt.prepared.assessment.admission_passed)
                with patch.object(
                    c1r1, "commit_c1r1_runtime", side_effect=counted
                ):
                    origin_checkpoint.adopt_fine(nxt.prepared)
                self.assertEqual(commits["count"], 2)

    def test_fine_then_temporal_source_cfl_restart_preserves_plan_and_totals(
        self,
    ) -> None:
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            with self.subTest(method=method):
                rhs = ControllableRHS()
                member = _make_member(method, rhs=rhs)
                checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
                rhs.mode = "ok"
                first = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
                self.assertEqual(first.disposition, ACCEPTED_FINE)
                accepted_encoding = checkpoint.accepted_encoding
                accepted_plan = checkpoint.last_accepted_executed_plan
                inherited = member.temporal_ledger.inherited_snapshot
                accepted_time = member.time.hex()
                rhs.impulse_from = member.time
                rhs.mode = "impulse"
                temporal = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
                self.assertEqual(temporal.disposition, "temporal_retry_required")
                self.assertEqual(member.executed_plan, accepted_plan)
                self.assertEqual(member.time.hex(), accepted_time)
                self.assertEqual(member.temporal_ledger.inherited_snapshot, inherited)
                checkpoint.restore(
                    accepted_encoding, encode_hlt17_cursor(checkpoint.cursor)
                )
                self.assertEqual(member.executed_plan, accepted_plan)
                self.assertEqual(checkpoint.member.source_retry_count, 0)
                rhs.impulse_from = member.time
                rhs.mode = "impulse"
                replay = replay_hlt17_rejection(
                    checkpoint.cursor, checkpoint.live_boundary()
                )
                rhs.mode = "source"
                source_prep = prepare_hlt17_next(
                    checkpoint.cursor, checkpoint.live_boundary(), replay=replay
                )
                source = checkpoint.absorb_nonfine(source_prep)
                self.assertEqual(source.disposition, "source_retry_required")
                self.assertEqual(member.source_retry_count, checkpoint.cursor.source_total)
                self.assertEqual(member.source_retry_count, 1)
                self.assertEqual(member.executed_plan, accepted_plan)
                checkpoint.restore(
                    accepted_encoding, encode_hlt17_cursor(checkpoint.cursor)
                )
                rhs.impulse_from = member.time
                rhs.mode = "impulse"
                replay = replay_hlt17_rejection(
                    checkpoint.cursor, checkpoint.live_boundary()
                )
                rhs.mode = "cfl"
                cfl_prep = prepare_hlt17_next(
                    checkpoint.cursor, checkpoint.live_boundary(), replay=replay
                )
                cfl = checkpoint.absorb_nonfine(cfl_prep)
                self.assertEqual(cfl.disposition, "cfl_retry_required")
                self.assertEqual(member.CFL_retry_count, 1)
                self.assertEqual(member.source_retry_count, 1)
                self.assertEqual(member.executed_plan, accepted_plan)
                self.assertEqual(member.temporal_ledger.inherited_snapshot, inherited)
                retry_cursor = encode_hlt17_cursor(checkpoint.cursor)
                origin_member = _make_member(method, rhs=ControllableRHS())
                restarted = HLT17InMemoryCheckpoint.seed(
                    origin_member, event_target=TARGET
                )
                restarted.restore(accepted_encoding, retry_cursor)
                self.assertEqual(origin_member.executed_plan, accepted_plan)
                self.assertEqual(origin_member.source_retry_count, 1)
                self.assertEqual(origin_member.CFL_retry_count, 1)
                self.assertEqual(origin_member.time.hex(), accepted_time)
                rhs_ok = origin_member.operator
                assert isinstance(rhs_ok, ControllableRHS)
                rhs_ok.impulse_from = origin_member.time
                rhs_ok.mode = "impulse"
                replay = replay_hlt17_rejection(
                    restarted.cursor, restarted.live_boundary()
                )
                rhs_ok.mode = "ok"
                nxt = prepare_hlt17_next(
                    restarted.cursor, restarted.live_boundary(), replay=replay
                )
                self.assertEqual(nxt.disposition, "prepared_overlay")
                admitted = restarted.require_admission(nxt.prepared)
                retry = restarted.adopt_fine(admitted.prepared)
                self.assertEqual(retry.disposition, ACCEPTED_FINE)
                self.assertEqual(restarted.cursor.source_current, 0)
                self.assertEqual(restarted.cursor.source_total, 1)
                self.assertEqual(restarted.cursor.cfl_total, 1)
                self.assertEqual(origin_member.source_retry_count, 1)
                self.assertEqual(
                    origin_member.temporal_ledger.inherited_snapshot, inherited
                )

    def test_late_failures_roll_back_the_real_member_boundary(self) -> None:
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        prepared = checkpoint.prepare(requested_cap=REQUESTED_CAP)
        admitted = checkpoint.require_admission(prepared.prepared)
        family = admitted.prepared
        assert family is not None
        before = _bits(member)
        before_cursor = encode_hlt17_cursor(checkpoint.cursor)
        before_descriptor = checkpoint.accepted_encoding.descriptor
        before_generation = checkpoint.accepted_generation

        for hook in (
            "after_kernel_adoption",
            "after_provisional_encoding",
            "after_cursor_finalization",
            "after_result_construction",
        ):
            with self.subTest(hook=hook):
                def boom(name: str, expected=hook) -> None:
                    if name == expected:
                        raise RuntimeError(f"injected {expected}")

                with self.assertRaisesRegex(RuntimeError, "injected"):
                    checkpoint.adopt_fine(family, fault_hook=boom)
                self.assertEqual(_bits(member), before)
                self.assertEqual(encode_hlt17_cursor(checkpoint.cursor), before_cursor)
                self.assertEqual(
                    checkpoint.accepted_encoding.descriptor, before_descriptor
                )
                self.assertEqual(checkpoint.accepted_generation, before_generation)
                self.assertEqual(member.transaction.state.accepted_stage_count, 0)
                self.assertIsNone(member.executed_plan)

        with patch.object(
            c1r1, "commit_c1r1_runtime", side_effect=RuntimeError("injected kernel")
        ):
            with self.assertRaisesRegex(RuntimeError, "injected kernel"):
                checkpoint.adopt_fine(family)
        self.assertEqual(_bits(member), before)
        self.assertEqual(member.transaction.state.accepted_stage_count, 0)

    def test_foreign_descriptor_cursor_ledger_and_source_fail_before_mutation(
        self,
    ) -> None:
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        result = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        encoding = result.encoding
        cursor_mapping = encode_hlt17_cursor(checkpoint.cursor)
        other = _make_member(closure="c" * 64, campaign_id="HLT17-FOREIGN-TEST")
        other_checkpoint = HLT17InMemoryCheckpoint.seed(other, event_target=TARGET)
        before = _bits(other)
        with self.assertRaises(Exception):
            other_checkpoint.restore(encoding, cursor_mapping)
        self.assertEqual(_bits(other), before)
        advanced = _bits(member)
        with self.assertRaises(TypeError):
            HLT17MemberEncoding(
                encoding.descriptor,
                encoding.payload,
                "d" * 64,
                encoding.physical_state_sha256,
                encoding.payload_sha256,
            )
        with self.assertRaises(Exception):
            checkpoint.restore(
                HLT17MemberEncoding(encoding.descriptor, encoding.payload[:-8]),
                cursor_mapping,
            )
        self.assertEqual(_bits(member), advanced)
        mapping = dict(cursor_mapping)
        mapping["descriptor_sha256"] = encoding.physical_state_sha256
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        self.assertEqual(other_checkpoint.accepted_generation, 0)
        self.assertEqual(_bits(other)["u"], before["u"])

    def test_park_at_event_target_does_not_call_source_or_retarget(self) -> None:
        rhs = ControllableRHS()
        member = _make_member(rhs=rhs)
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        first = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(first.disposition, ACCEPTED_FINE)
        second = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(second.disposition, ACCEPTED_FINE)
        self.assertEqual(member.time.hex(), TARGET.hex())
        calls_before = rhs.calls
        with patch.object(
            c1r1, "prepare_c1r1_initial_runtime", side_effect=AssertionError("source")
        ), patch.object(
            c1r1, "_prepare", side_effect=AssertionError("source")
        ), patch.object(
            imp1, "prepare_imp1_initial_runtime", side_effect=AssertionError("imp1 source")
        ), patch.object(
            imp1, "_prepare", side_effect=AssertionError("imp1 source")
        ):
            parked = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(parked.disposition, PARKED_AT_EVENT_TARGET)
        self.assertTrue(parked.parked)
        self.assertFalse(parked.accepted_state_advanced)
        self.assertFalse(parked.evidence["source_invoked"])
        self.assertFalse(parked.evidence["retarget_authorized"])
        self.assertEqual(rhs.calls, calls_before)
        self.assertEqual(member.time.hex(), TARGET.hex())

    def test_export_and_agree_reject_forged_bytes_hashes_and_live_mutation(self) -> None:
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        bundle = checkpoint.generation_bundle()
        with self.assertRaises(HLT17CheckpointError):
            replace(bundle, payload=b"bad")
        with self.assertRaises(TypeError):
            replace(bundle, payload_sha256="e" * 64)
        with self.assertRaises(TypeError):
            replace(bundle, campaign_execution_authorized=True)
        with self.assertRaises(TypeError):
            replace(bundle, accepted_generation=-1)
        with self.assertRaises(HLT17MemberCodecError):
            replace(checkpoint.accepted_encoding, payload=b"bad")
        with self.assertRaises((HLT17CheckpointError, HLT17MemberCodecError, TypeError)):
            checkpoint.accepted_encoding = replace(
                checkpoint.accepted_encoding, payload=b"bad"
            )
        checkpoint.agree()
        exported = checkpoint.generation_bundle()
        self.assertNotEqual(exported.payload, b"bad")
        self.assertEqual(exported.payload, bundle.payload)

        before = _bits(member)
        member.tracers.proper_times += 1.0
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.agree()
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.generation_bundle()
        member.tracers.proper_times -= 1.0
        checkpoint.agree()
        self.assertEqual(_bits(member)["proper"], before["proper"])

        original_monitor = member.transaction.state
        member.transaction.state = replace(
            original_monitor, accepted_stage_count=7
        )
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.agree()
        member.transaction.state = original_monitor
        original_causal = member.transaction.causal_state
        member.transaction.causal_state = replace(
            original_causal, previous_speed_upper=original_causal.previous_speed_upper * 2.0
        )
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.agree()
        member.transaction.causal_state = original_causal
        original_cfl = member.transaction.cfl_maximum
        member.transaction.cfl_maximum = original_cfl / 2.0
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.agree()
        member.transaction.cfl_maximum = original_cfl
        checkpoint.agree()

        mapping = encode_hlt17_cursor(checkpoint.cursor)
        mapping["kernel_transition_digest"] = "f" * 64
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(checkpoint.cursor)
        mapping["campaign_execution_authorized"] = True
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(checkpoint.cursor)
        mapping["cursor_generation"] = checkpoint.cursor.cursor_generation + 3
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(checkpoint.cursor)
        mapping["monitor_sha256"] = "a" * 64
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(checkpoint.cursor)
        mapping["overlay"] = {
            "owner": "source",
            "attempted_plan": None,
            "next_plan": None,
            "evidence": {},
            "owner_current_count": 1,
            "exhausted": False,
        }
        with self.assertRaises(Exception):
            restore_hlt17_cursor(mapping)
        with self.assertRaises(HLT17CheckpointError):
            HLT17GenerationBundle(
                descriptor=bundle.descriptor,
                payload=bundle.payload[:-8],
                cursor_bytes=bundle.cursor_bytes,
            )
        with self.assertRaises(HLT17CheckpointError):
            HLT17GenerationBundle(
                descriptor=bundle.descriptor[:-1],
                payload=bundle.payload,
                cursor_bytes=bundle.cursor_bytes,
            )
        with self.assertRaises(HLT17CheckpointError):
            HLT17GenerationBundle(
                descriptor=bundle.descriptor,
                payload=bundle.payload,
                cursor_bytes=bundle.cursor_bytes + b"\x00",
            )

    def test_synthetic_origin_placeholder_is_not_a_real_proto17_cursor(self) -> None:
        member = _make_member()
        predecessor = origin_predecessor_identity(member)
        identity = dict(predecessor.cursor_identity)
        expected = sha256(
            SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN + canonical_json_bytes(identity)
        ).hexdigest()
        self.assertEqual(predecessor.cursor_sha256, expected)
        self.assertNotEqual(predecessor.cursor_sha256, predecessor.identity_sha256)
        self.assertIn(b"synthetic-origin-cursor-placeholder", SYNTHETIC_ORIGIN_CURSOR_PLACEHOLDER_DOMAIN)
        self.assertIs(
            type(synthetic_origin_predecessor_identity(member)), type(predecessor)
        )

    def test_live_seeding_requires_explicit_origin_reference(self) -> None:
        state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(
            PRIMARY_METHOD, point_count=2049
        )
        live = HLT17RuntimeMember(
            identity=HLT17MemberIdentity(
                scope=IDENTITY_SCOPE_LIVE,
                method_label="RK4",
                point_count=2049,
                integrator_id=PRIMARY_METHOD,
                spatial_order=4,
                campaign_id="HLT17-LIVE-TEST",
                amplitude="3",
            ),
            initial=ProvidedInitial(
                ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))
            ),
            operator=synthetic_exponential_rhs,
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
        )
        with self.assertRaises(HLT17CheckpointError):
            synthetic_origin_predecessor_identity(live)
        with self.assertRaises(HLT17CheckpointError):
            origin_predecessor_identity(live)
        with self.assertRaises(HLT17CheckpointError):
            HLT17InMemoryCheckpoint.seed(live, event_target=TARGET)
        supplied = origin_predecessor_identity(live, origin_cursor_sha256="b" * 64)
        self.assertEqual(supplied.cursor_sha256, "b" * 64)
        self.assertFalse(supplied.cursor_authenticated)

    def test_restore_rejects_foreign_last_accepted_plan(self) -> None:
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        result = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        encoding = result.encoding
        cursor_mapping = encode_hlt17_cursor(checkpoint.cursor)
        before = _bits(member)
        wrong = plan_forward_proto14_subdivision(
            1.5, 1.5625, 1.0 / 32, minimum_width=float(TDG6_MINIMUM_MACRO_STEP)
        )
        self.assertNotEqual(wrong, encoding.executed_plan)
        with self.assertRaises(HLT17CheckpointError):
            checkpoint.restore(
                encoding,
                cursor_mapping,
                last_accepted_executed_plan=wrong,
            )
        self.assertEqual(_bits(member), before)
        checkpoint.restore(encoding, cursor_mapping, last_accepted_executed_plan=encoding.executed_plan)
        self.assertEqual(member.executed_plan, encoding.executed_plan)


if __name__ == "__main__":
    unittest.main()
