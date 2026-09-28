"""C1R1 in-memory HLT17 integration: provenance, both methods, and refusals."""

from __future__ import annotations

from dataclasses import dataclass
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
from recursive_horizons.fgc.evolution.hlt17_admission_runtime import (  # noqa: E402
    HLT17_IMPLEMENTATION_ID,
    hlt17_implementation_identity,
    require_hlt17_prepared,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    HLT17LiveBoundary,
    prepare_hlt17_fresh,
    prepare_hlt17_next,
    replay_hlt17_rejection,
    require_hlt17_admission,
    seed_cursor_from_boundary,
    with_ledger,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17IMP1Error,
    HLT17_CURSOR_SCHEMA,
    HLT17_LEGACY_IMP1_CURSOR_SCHEMA,
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    ACCEPTED_FINE,
    HLT17InMemoryCheckpoint,
)
from recursive_horizons.fgc.evolution.hlt17_member_codec import (  # noqa: E402
    CODEC_LEGACY_IMP1_SCHEMA,
    CODEC_SCHEMA,
    HLT17MemberCodecError,
    HLT17MemberEncoding,
)
from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_SYNTHETIC,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure import (  # noqa: E402
    C1R1_IMPLEMENTATION_ID,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)


REQUESTED_CAP = 1.0 / 8.0
BRIDGE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_imp1_bridge.py"
MEMBER = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_runtime_member.py"
CHECKPOINT = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_member_checkpoint.py"


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
        self.mode = "impulse"

    def hlt17_synthetic_callable_binding(self) -> dict[str, object]:
        code = type(self).__call__.__code__
        return {
            "kind": "synthetic_failure_injection_double",
            "type": f"{type(self).__module__}.{type(self).__qualname__}",
            "call_code_sha256": sha256(code.co_code).hexdigest(),
            "call_names": list(code.co_names),
            "mutable_injection_mode_excluded": True,
            "physical_source_authenticated": False,
        }

    def __call__(self, time: float, current) -> EvolutionRHS:
        diagnostics = synthetic_diagnostics()
        zeros = np.zeros_like(current.u)
        if self.mode == "source":
            if time > START:
                diagnostics["source_residual_infinity"] = 1.0
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "cfl":
            if time > START:
                diagnostics["coordinate_speed_upper"] = 1.0e9
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "ok":
            return synthetic_exponential_rhs(time, current)
        du = np.zeros_like(current.u)
        if time == START + 1.0 / 128.0:
            du[:, 4] = 1.0
        return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), diagnostics)


def _descriptor(physical: str) -> str:
    digest = sha256(b"hlt17-c1r1-descriptor:" + physical.encode("ascii")).hexdigest()
    if digest == physical:
        digest = sha256(b"hlt17-c1r1-descriptor-alt").hexdigest()
    return digest


def _identity(method: str, point_count: int) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_SYNTHETIC,
        method_label=label,
        point_count=point_count,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id="HLT17-C1R1-INTEGRATION",
        amplitude="synthetic",
    )


def _make_member(method: str = PRIMARY_METHOD, *, rhs=synthetic_exponential_rhs) -> HLT17RuntimeMember:
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
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


def _boundary(rhs: ControllableRHS | None = None, *, method: str = PRIMARY_METHOD):
    from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (
        imp1_enclosure_limits_for_member,
    )

    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    rhs = ControllableRHS() if rhs is None else rhs
    physical = array_content_sha256(state.u, state.p, state.q)
    member_key = "RK4-9" if method == PRIMARY_METHOD else "SSPRK3-9"
    return HLT17LiveBoundary(
        member_key=member_key,
        state=state,
        rhs=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
        descriptor_sha256=_descriptor(physical),
        event_target=TARGET,
        limits=imp1_enclosure_limits_for_member(member_key),
    )


class HLT17C1R1IntegrationTests(unittest.TestCase):
    def test_identity_is_the_named_c1r1_object_not_imp1_execution(self) -> None:
        identity = hlt17_implementation_identity()
        self.assertEqual(identity["implementation_id"], C1R1_IMPLEMENTATION_ID)
        self.assertEqual(identity["implementation_id"], HLT17_IMPLEMENTATION_ID)
        self.assertEqual(identity["reference_wire_artifact_id"], imp1.ARTIFACT_ID)
        self.assertNotEqual(identity["implementation_id"], imp1.ARTIFACT_ID)
        self.assertEqual(HLT17_CURSOR_SCHEMA, "FGC-1-HLT17-C1R1-cursor-v1")
        self.assertEqual(CODEC_SCHEMA, "FGC-1-HLT17-C1R1-member-state-v1")

    def test_owned_prepare_paths_do_not_call_sealed_imp1_prepare(self) -> None:
        for path in (BRIDGE, MEMBER, CHECKPOINT):
            source = path.read_text(encoding="utf-8")
            self.assertNotIn("prepare_imp1_initial_runtime(", source)
            self.assertNotIn("prepare_imp1_retry_runtime(", source)
            self.assertNotIn("commit_imp1_runtime(", source)
            self.assertNotIn("require_imp1_admission(", source)
            self.assertNotIn("imp1._prepare(", source)

    def test_both_methods_seed_encode_restore_and_fine_commit(self) -> None:
        for method, serial in ((PRIMARY_METHOD, 20), (COMPARATOR_METHOD, 16)):
            with self.subTest(method=method):
                member = _make_member(method)
                inherited = member.temporal_ledger.inherited_snapshot
                inherited_debit = member.inherited_tdg6_sibling().accumulated_debit_vector
                with patch.object(
                    imp1, "_prepare", side_effect=AssertionError("imp1 prepare")
                ), patch.object(
                    imp1,
                    "prepare_imp1_initial_runtime",
                    side_effect=AssertionError("imp1 prepare"),
                ), patch.object(
                    imp1,
                    "commit_imp1_runtime",
                    side_effect=AssertionError("imp1 commit"),
                ):
                    checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
                    mapping = encode_hlt17_cursor(checkpoint.cursor)
                    restored = restore_hlt17_cursor(mapping)
                    self.assertEqual(encode_hlt17_cursor(restored), mapping)
                    self.assertEqual(
                        mapping["implementation_id"], C1R1_IMPLEMENTATION_ID
                    )
                    metadata = load_canonical_json(
                        checkpoint.accepted_encoding.descriptor
                    )
                    self.assertEqual(metadata["schema"], CODEC_SCHEMA)
                    self.assertEqual(
                        metadata["runtime_identity"]["implementation_id"],
                        C1R1_IMPLEMENTATION_ID,
                    )
                    self.assertFalse(metadata["historical_origin_authenticated"])
                    self.assertFalse(metadata["campaign_execution_authorized"])
                    self.assertFalse(metadata["durable_store_bound"])
                    result = checkpoint.advance_once(requested_cap=1.0 / 32)
                self.assertEqual(result.disposition, ACCEPTED_FINE)
                self.assertIsInstance(result.prepared, c1r1.TDG11C1R1PreparedRuntime)
                self.assertIsInstance(result.committed, c1r1.TDG11C1R1CommittedRuntime)
                self.assertFalse(result.committed.corrected_endpoint_adopted)
                self.assertEqual(member.transaction_serial, serial)
                self.assertEqual(member.temporal_ledger.inherited_snapshot, inherited)
                self.assertEqual(
                    member.inherited_tdg6_sibling().accumulated_debit_vector,
                    inherited_debit,
                )
                self.assertEqual(
                    result.cursor.implementation_id, C1R1_IMPLEMENTATION_ID
                )

    def test_rejected_temporal_replay_source_overlay_retry_and_counters(self) -> None:
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            with self.subTest(method=method):
                rhs = ControllableRHS()
                boundary = _boundary(rhs, method=method)
                cursor = seed_cursor_from_boundary(boundary)
                inherited = cursor.ledger.inherited_snapshot
                debit = cursor.ledger.accumulated_debit_vector
                rhs.mode = "impulse"
                prepared = prepare_hlt17_fresh(
                    cursor, boundary, requested_cap=REQUESTED_CAP
                )
                self.assertFalse(prepared.prepared.assessment.admission_passed)
                self.assertEqual(
                    prepared.prepared.implementation_id, C1R1_IMPLEMENTATION_ID
                )
                closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
                self.assertEqual(closed.disposition, "temporal_retry_required")
                self.assertEqual(closed.sink_writes, 1)
                self.assertEqual(closed.cursor.ledger.inherited_snapshot, inherited)
                self.assertEqual(closed.cursor.ledger.accumulated_debit_vector, debit)
                restored = restore_hlt17_cursor(encode_hlt17_cursor(closed.cursor))
                live = with_ledger(boundary, restored.ledger)
                rhs.mode = "impulse"
                replay = replay_hlt17_rejection(restored, live)
                self.assertEqual(replay.sink_writes, 0)
                self.assertIsInstance(
                    replay.reconstructed, c1r1.TDG11C1R1PreparedRuntime
                )
                rhs.mode = "source"
                source = prepare_hlt17_next(restored, live, replay=replay)
                self.assertEqual(source.disposition, "source_retry_required")
                self.assertEqual(source.cursor.source_current, 1)
                self.assertEqual(source.cursor.source_total, 1)
                self.assertEqual(
                    source.cursor.ledger.current_macro_step_temporal_retry_count, 1
                )
                self.assertEqual(source.cursor.ledger.inherited_snapshot, inherited)
                self.assertEqual(source.cursor.ledger.accumulated_debit_vector, debit)
                live = with_ledger(boundary, source.cursor.ledger)
                rhs.mode = "impulse"
                replay = replay_hlt17_rejection(source.cursor, live)
                rhs.mode = "ok"
                nxt = prepare_hlt17_next(source.cursor, live, replay=replay)
                self.assertEqual(nxt.disposition, "prepared_overlay")
                self.assertTrue(nxt.prepared.assessment.admission_passed)

    def test_imp1_prepared_and_schema_swaps_are_refused(self) -> None:
        rhs = ControllableRHS()
        rhs.mode = "ok"
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        foreign = imp1.prepare_imp1_initial_runtime(
            method=cursor.method,
            time=START,
            event_target=TARGET,
            requested_cap=REQUESTED_CAP,
            state=boundary.state,
            rhs=boundary.rhs,
            projector=None,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
            limits=boundary.resolved_limits(),
        )
        self.assertTrue(foreign.assessment.admission_passed)
        with self.assertRaises(TypeError):
            require_hlt17_prepared(foreign)
        with self.assertRaises(TypeError):
            require_hlt17_admission(cursor, boundary, foreign)

        mapping = encode_hlt17_cursor(cursor)
        mapping["schema"] = HLT17_LEGACY_IMP1_CURSOR_SCHEMA
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(cursor)
        mapping["implementation_id"] = "not-the-fixed-c1r1-object"
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        metadata = dict(load_canonical_json(checkpoint.accepted_encoding.descriptor))
        metadata["schema"] = CODEC_LEGACY_IMP1_SCHEMA
        with self.assertRaises(HLT17MemberCodecError):
            HLT17MemberEncoding(
                canonical_json_bytes(metadata),
                checkpoint.accepted_encoding.payload,
            )

    def test_fine_rollback_does_not_substitute_endpoint_or_rerun_rejected_write(self) -> None:
        member = _make_member()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        origin = member.capture()
        origin_time = member.time.hex()
        prepared = checkpoint.prepare(requested_cap=1.0 / 32)
        admitted = checkpoint.require_admission(prepared.prepared)
        self.assertTrue(admitted.evidence["admission_passed"])
        with patch.object(
            c1r1, "commit_c1r1_runtime", side_effect=RuntimeError("injected kernel")
        ):
            with self.assertRaises(RuntimeError):
                checkpoint.adopt_fine(admitted.prepared)
        self.assertEqual(member.time.hex(), origin_time)
        self.assertEqual(member.step_index, origin.step_index)
        self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 0)
        self.assertIsNone(member.executed_plan)
        result = checkpoint.adopt_fine(admitted.prepared)
        self.assertEqual(result.disposition, ACCEPTED_FINE)
        self.assertFalse(result.committed.corrected_endpoint_adopted)
        self.assertFalse(result.committed.outer_path_adopted)
        self.assertFalse(result.committed.medium_path_adopted)


if __name__ == "__main__":
    unittest.main()
