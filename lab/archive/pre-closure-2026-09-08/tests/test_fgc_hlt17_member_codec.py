"""Focused lossless HLT17 codec controls; no store writer or campaign I/O."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes, load_canonical_json  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_member_codec import (  # noqa: E402
    ARRAY_NAMES,
    CODEC_SCHEMA,
    HLT17DecodedSnapshot,
    HLT17MemberCodecError,
    HLT17MemberCodecResourceError,
    HLT17MemberEncoding,
    HLT17PredecessorIdentity,
    MAX_PAYLOAD_BYTES,
    ORIGIN_PROTOCOL,
    PREDECESSOR_CURSOR_IDENTITY_KEYS,
    RUNTIME_PROTOCOL,
    SNAPSHOT_ROLE_ACCEPTED,
    SNAPSHOT_ROLE_RETRY,
    decode_member,
    encode_member,
    encode_plan,
    encode_retry_overlay,
    restore_member,
    validate_descriptor,
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
    EvolutionState,
    array_content_sha256,
)
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
    replan_c1r1_retry,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (  # noqa: E402
    plan_forward_proto14_subdivision,
)
from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState  # noqa: E402


CODEC_SOURCE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_member_codec.py"


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
        tracers.event_fields[0] = np.array(tracers.event_fields[0], dtype="<f8", copy=True)
        tracers.event_fields[0][0, 0] = -0.0
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


def _cursor_identity(
    member: HLT17RuntimeMember,
    *,
    protocol: str,
    accepted_generation: int,
) -> dict[str, object]:
    return {
        "protocol_artifact_id": protocol,
        "campaign_id": member.identity.campaign_id,
        "member_key": member.key,
        "method": member.identity.method_label,
        "point_count": member.identity.point_count,
        "accepted_boundary_time_hex": member.time.hex(),
        "accepted_state_sha256": member.physical_state_sha256(),
        "previous_step_index": member.step_index,
        "previous_transaction_serial": member.transaction_serial,
        "accepted_generation": accepted_generation,
    }


def _full_cursor_digest(identity: dict[str, object]) -> str:
    return sha256(
        b"hlt17-full-cursor-reference\n" + canonical_json_bytes(identity)
    ).hexdigest()


def _predecessor(
    member: HLT17RuntimeMember,
    *,
    protocol: str,
    accepted_generation: int,
    descriptor_sha256: str,
) -> HLT17PredecessorIdentity:
    identity = _cursor_identity(
        member, protocol=protocol, accepted_generation=accepted_generation
    )
    return HLT17PredecessorIdentity(
        descriptor_sha256=descriptor_sha256,
        cursor_identity=identity,
        cursor_sha256=_full_cursor_digest(identity),
    )


def _origin_predecessor(member: HLT17RuntimeMember) -> HLT17PredecessorIdentity:
    return _predecessor(
        member,
        protocol=ORIGIN_PROTOCOL,
        accepted_generation=0,
        descriptor_sha256=member.origin_reference.sha256,
    )


def _evolved_predecessor(
    member: HLT17RuntimeMember, origin_encoding, origin_member: HLT17RuntimeMember
) -> HLT17PredecessorIdentity:
    return _predecessor(
        origin_member,
        protocol=RUNTIME_PROTOCOL,
        accepted_generation=0,
        descriptor_sha256=origin_encoding.descriptor_sha256,
    )


def _same_boundary_predecessor(
    member: HLT17RuntimeMember, accepted_encoding: HLT17MemberEncoding, *, accepted_generation: int
) -> HLT17PredecessorIdentity:
    return _predecessor(
        member,
        protocol=RUNTIME_PROTOCOL,
        accepted_generation=accepted_generation,
        descriptor_sha256=accepted_encoding.descriptor_sha256,
    )


def _copy_decoded(decoded: HLT17DecodedSnapshot, **changes: object) -> HLT17DecodedSnapshot:
    payload = {
        "arrays": decoded.arrays,
        "metadata": decoded.metadata,
        "ledger": decoded.ledger,
        "executed_plan": decoded.executed_plan,
        "descriptor": decoded.descriptor,
        "payload": decoded.payload,
        "descriptor_sha256": decoded.descriptor_sha256,
        "physical_state_sha256": decoded.physical_state_sha256,
        "payload_sha256": decoded.payload_sha256,
    }
    payload.update(changes)
    return HLT17DecodedSnapshot(**payload)


def _quarter_impulse(time, current):
    du = np.zeros_like(current.u)
    if time == START + 1.0 / 128:
        du[:, 4] = 1.0
    return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), synthetic_diagnostics())


def _payload_map(encoding) -> dict[str, bytes]:
    decoded = decode_member(encoding.descriptor, encoding.payload)
    return {name: decoded.arrays[name].tobytes() for name in ARRAY_NAMES}


class HLT17MemberCodecTests(unittest.TestCase):
    def test_source_stays_pure_and_does_not_guess_the_cursor(self) -> None:
        source = CODEC_SOURCE.read_text(encoding="utf-8")
        self.assertNotIn("hlt16_member_codec", source)
        self.assertNotIn("hlt16_state_store", source)
        self.assertNotIn("publish_exclusive", source)
        self.assertNotIn(".advance_to(", source)
        self.assertIn("predecessor_cursor_identity", source)
        self.assertIn("predecessor_cursor_identity_sha256", source)
        self.assertIn("Not an HLT17 cursor object", source)
        self.assertIn(ORIGIN_PROTOCOL, source)
        self.assertNotIn("FGC-1-HLT15-origin-reference", source)
        self.assertEqual(tuple(PREDECESSOR_CURSOR_IDENTITY_KEYS)[-1], "protocol_artifact_id")
        self.assertEqual(ARRAY_NAMES[0], "u")
        self.assertEqual(len(ARRAY_NAMES), 9)
        self.assertEqual(SNAPSHOT_ROLE_ACCEPTED, "accepted_physical_state")
        self.assertEqual(SNAPSHOT_ROLE_RETRY, "retry_ledger_overlay")

    def test_prepare_commit_encode_decode_restore_next_step(self) -> None:
        for method, serial in ((PRIMARY_METHOD, 20), (COMPARATOR_METHOD, 16)):
            with self.subTest(method=method):
                member = _make_member(method, signed_zeros=True)
                origin = _origin_predecessor(member)
                origin_encoding = encode_member(member, origin, generation=0)
                origin_meta = load_canonical_json(origin_encoding.descriptor)
                self.assertEqual(origin_meta["snapshot_role"], SNAPSHOT_ROLE_ACCEPTED)
                self.assertEqual(
                    origin_meta["predecessor_cursor_identity"]["protocol_artifact_id"],
                    ORIGIN_PROTOCOL,
                )
                self.assertNotEqual(
                    origin_meta["predecessor_cursor_identity_sha256"],
                    origin_meta["predecessor_cursor_sha256"],
                )
                equal = dict(origin_meta)
                equal["predecessor_cursor_sha256"] = origin_meta[
                    "predecessor_cursor_identity_sha256"
                ]
                with self.assertRaises(HLT17MemberCodecError):
                    validate_descriptor(equal)
                self.assertNotEqual(
                    origin_encoding.descriptor_sha256,
                    origin_encoding.physical_state_sha256,
                )
                origin_payload = origin_encoding.payload
                origin_descriptor = origin_encoding.descriptor
                origin_bits = _payload_map(origin_encoding)
                self.assertEqual(
                    np.float64.fromhex(
                        load_canonical_json(origin_descriptor)["causal"][
                            "accumulated_characteristic_distance_hex"
                        ]
                    ).tobytes(),
                    np.float64(-0.0).tobytes(),
                )
                self.assertIn(np.float64(-0.0).tobytes(), origin_bits["u"])
                self.assertIn(np.float64(-0.0).tobytes(), origin_bits["event_fields"])
                prepared = member.prepare_fresh(
                    event_target=TARGET, requested_cap=1.0 / 32
                )
                self.assertTrue(prepared.assessment.admission_passed)
                prior = _make_member(method, signed_zeros=True)
                member.commit_prepared(prepared)
                self.assertEqual(member.transaction_serial, serial)
                evolved = encode_member(
                    member, _evolved_predecessor(member, origin_encoding, prior), generation=1
                )
                self.assertNotEqual(evolved.descriptor, origin_descriptor)
                self.assertIsNotNone(load_canonical_json(evolved.descriptor)["executed_plan"])
                twin = _make_member(method, signed_zeros=True)
                restore_member(twin, evolved)
                again = encode_member(
                    twin, _evolved_predecessor(twin, origin_encoding, prior), generation=1
                )
                self.assertEqual(again.descriptor, evolved.descriptor)
                self.assertEqual(again.payload, evolved.payload)
                self.assertEqual(again.physical_state_sha256, evolved.physical_state_sha256)
                nxt = twin.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
                self.assertTrue(nxt.assessment.admission_passed)
                self.assertEqual(twin.temporal_ledger.inherited_snapshot, member.temporal_ledger.inherited_snapshot)
                restore_member(twin, origin_encoding)
                restored_origin = encode_member(twin, origin, generation=0)
                self.assertEqual(restored_origin.descriptor, origin_descriptor)
                self.assertEqual(restored_origin.payload, origin_payload)

    def test_retry_ledger_round_trip_preserves_rejection_and_inherited_bytes(self) -> None:
        member = _make_member(rhs=_quarter_impulse)
        origin = _origin_predecessor(member)
        accepted = encode_member(member, origin, generation=0)
        inherited = member.temporal_ledger.inherited_snapshot
        debit = member.temporal_ledger.accumulated_debit_vector
        physical = member.physical_state_sha256()
        accepted_time = member.time.hex()
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        written: list[object] = []
        with self.assertRaises(C1R1TemporalRetryRequired) as caught:
            member.require_admission(prepared, durable_rejection_sink=written.append)
        member.adopt_retry_ledger(caught.exception)
        with self.assertRaises(HLT17MemberCodecError):
            encode_member(member, origin, generation=0)
        encoding = encode_retry_overlay(
            member,
            _same_boundary_predecessor(member, accepted, accepted_generation=0),
            accepted_generation=0,
        )
        metadata = load_canonical_json(encoding.descriptor)
        checkpoint = metadata["imp1_checkpoint"]
        self.assertEqual(metadata["snapshot_role"], SNAPSHOT_ROLE_RETRY)
        self.assertEqual(metadata["runtime_identity"]["accepted_generation"], 0)
        self.assertEqual(checkpoint["current_macro_step_temporal_retry_count"], 1)
        self.assertEqual(checkpoint["inherited_snapshot"], inherited)
        self.assertEqual(len(checkpoint["serialized_temporal_rejections"]), 1)
        self.assertEqual(member.physical_state_sha256(), physical)
        self.assertEqual(member.time.hex(), accepted_time)
        twin = _make_member(rhs=_quarter_impulse)
        restore_member(twin, encoding)
        self.assertEqual(twin.temporal_ledger.inherited_snapshot, inherited)
        self.assertEqual(twin.temporal_ledger.accumulated_debit_vector, debit)
        self.assertEqual(twin.temporal_ledger.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(twin.physical_state_sha256(), physical)
        self.assertEqual(twin.time.hex(), accepted_time)
        again = encode_retry_overlay(
            twin,
            _same_boundary_predecessor(twin, accepted, accepted_generation=0),
            accepted_generation=0,
        )
        self.assertEqual(again.descriptor, encoding.descriptor)
        self.assertEqual(again.payload, encoding.payload)
        successor = replan_c1r1_retry(prepared, caught.exception)
        nxt = twin.prepare_retry(successor)
        self.assertTrue(nxt.assessment.admission_passed)

    def test_stale_template_closure_and_coordinates_fail_before_mutation(self) -> None:
        member = _make_member()
        encoding = encode_member(member, _origin_predecessor(member), generation=0)
        before_hash = member.physical_state_sha256()
        metadata = load_canonical_json(encoding.descriptor)
        metadata["runtime_template_sha256"] = "f" * 64
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(
                member,
                type(encoding)(
                    canonical_json_bytes(metadata),
                    encoding.payload,
                ),
            )
        self.assertEqual(member.physical_state_sha256(), before_hash)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["source_operator_closure"]["sha256"] = "b" * 64
        mutated_descriptor = canonical_json_bytes(metadata)
        mutated = type(encoding)(
            mutated_descriptor,
            encoding.payload,
        )
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(member, mutated)
        other = _make_member()
        other.source_operator_closure = HLT17ExternalReference(
            kind="source_operator_closure", sha256="c" * 64
        )
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(other, encoding)
        self.assertEqual(other.step_index, 0)
        member.initial.grid.coordinates = np.linspace(0.0, 2.0, member.point_count)
        with self.assertRaises(Exception):
            encode_member(member, _origin_predecessor(member), generation=0)

    def test_descriptor_and_physical_hash_swaps_are_rejected(self) -> None:
        member = _make_member()
        encoding = encode_member(member, _origin_predecessor(member), generation=0)
        metadata = load_canonical_json(encoding.descriptor)
        swapped = dict(metadata)
        swapped["physical_state_sha256"] = encoding.descriptor_sha256
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(swapped), encoding.payload)
        swapped = dict(metadata)
        swapped["payload_sha256"] = encoding.physical_state_sha256
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(swapped), encoding.payload)
        swapped = dict(metadata)
        swapped["predecessor_descriptor_sha256"] = encoding.physical_state_sha256
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(swapped), encoding.payload)
        self.assertNotEqual(encoding.descriptor_sha256, encoding.physical_state_sha256)
        self.assertNotEqual(encoding.payload_sha256, encoding.physical_state_sha256)

    def test_malformed_payloads_fail_before_member_mutation(self) -> None:
        member = _make_member()
        encoding = encode_member(member, _origin_predecessor(member), generation=0)
        before = member.physical_state_sha256()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(encoding.descriptor, encoding.payload[:-8])
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(encoding.descriptor, encoding.payload + b"\x00")
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(encoding.descriptor[:-1], encoding.payload)
        pretty = json.dumps(load_canonical_json(encoding.descriptor), indent=2).encode("ascii")
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(pretty, encoding.payload)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["campaign_execution_authorized"] = 0
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), encoding.payload)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["runtime_identity"]["accepted_generation"] = 0.0
        with self.assertRaises(HLT17MemberCodecError):
            validate_descriptor(metadata)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["successor_cursor"] = {"sha256": "d" * 64}
        with self.assertRaises(HLT17MemberCodecError):
            validate_descriptor(metadata)
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(
                member,
                type(encoding)(encoding.descriptor, encoding.payload[:-8]),
            )
        self.assertEqual(member.physical_state_sha256(), before)

    def test_oversized_and_bool_aliases_are_resource_or_schema_stops(self) -> None:
        member = _make_member()
        encoding = encode_member(member, _origin_predecessor(member), generation=0)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["payload_byte_count"] = MAX_PAYLOAD_BYTES + 1
        with self.assertRaises(HLT17MemberCodecResourceError):
            validate_descriptor(metadata)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["historical_origin_authenticated"] = True
        with self.assertRaises(HLT17MemberCodecError):
            validate_descriptor(metadata)
        metadata = load_canonical_json(encoding.descriptor)
        metadata["source_operator_closure_authenticated"] = True
        with self.assertRaises(HLT17MemberCodecError):
            validate_descriptor(metadata)
        metadata = load_canonical_json(encoding.descriptor)
        self.assertEqual(metadata["schema"], CODEC_SCHEMA)
        from recursive_horizons.fgc.evolution.hlt17_admission_runtime import (
            HLT17_IMPLEMENTATION_ID,
            hlt17_implementation_identity,
        )

        identity = hlt17_implementation_identity()
        runtime = metadata["runtime_identity"]
        self.assertEqual(runtime["implementation_id"], HLT17_IMPLEMENTATION_ID)
        self.assertEqual(runtime["admission_runtime_id"], identity["admission_runtime_id"])
        self.assertEqual(
            runtime["reference_wire_artifact_id"],
            identity["reference_wire_artifact_id"],
        )

    def test_late_codec_restore_failure_rolls_back(self) -> None:
        member = _make_member()
        origin_encoding = encode_member(member, _origin_predecessor(member), generation=0)
        origin_hash = member.physical_state_sha256()
        prior = _make_member()
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        member.commit_prepared(prepared)
        evolved = encode_member(
            member, _evolved_predecessor(member, origin_encoding, prior), generation=1
        )
        restore_member(member, origin_encoding)
        self.assertEqual(member.physical_state_sha256(), origin_hash)

        def boom(name: str) -> None:
            if name == "after_restore_mutation":
                raise RuntimeError("injected late codec restore failure")

        with self.assertRaisesRegex(RuntimeError, "injected late codec restore failure"):
            restore_member(member, evolved, fault_hook=boom)
        self.assertEqual(member.physical_state_sha256(), origin_hash)
        self.assertEqual(member.step_index, 0)
        self.assertIsNone(member.executed_plan)
        self.assertEqual(
            encode_member(member, _origin_predecessor(member), generation=0).payload,
            origin_encoding.payload,
        )

    def test_forged_decoded_snapshot_and_encoding_hashes_fail_before_mutation(self) -> None:
        member = _make_member()
        encoding = encode_member(member, _origin_predecessor(member), generation=0)
        decoded = decode_member(encoding.descriptor, encoding.payload)
        before = member.physical_state_sha256()
        labels = np.asarray(member.tracers.proper_times).tobytes()
        forged_arrays = {
            name: np.array(decoded.arrays[name], dtype="<f8", copy=True)
            for name in ARRAY_NAMES
        }
        forged_arrays["tracer_proper_times"] = forged_arrays["tracer_proper_times"] + 1.0
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(member, _copy_decoded(decoded, arrays=forged_arrays))
        self.assertEqual(member.physical_state_sha256(), before)
        self.assertEqual(np.asarray(member.tracers.proper_times).tobytes(), labels)
        for name in ARRAY_NAMES:
            arrays = {
                key: np.array(decoded.arrays[key], dtype="<f8", copy=True)
                for key in ARRAY_NAMES
            }
            arrays[name] = arrays[name] + 1.0
            with self.assertRaises(HLT17MemberCodecError):
                restore_member(member, _copy_decoded(decoded, arrays=arrays))
        self.assertEqual(member.physical_state_sha256(), before)
        mutated_metadata = dict(decoded.metadata)
        mutated_metadata["counters"] = dict(mutated_metadata["counters"])
        mutated_metadata["counters"]["source_retry_count"] = 7
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(member, _copy_decoded(decoded, metadata=mutated_metadata))
        advanced = _make_member()
        advanced.commit_prepared(
            advanced.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        )
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(member, _copy_decoded(decoded, ledger=advanced.temporal_ledger))
        with self.assertRaises(TypeError):
            HLT17MemberEncoding(
                encoding.descriptor,
                encoding.payload,
                encoding.descriptor_sha256,
                encoding.physical_state_sha256,
                "e" * 64,
            )
        with self.assertRaises(HLT17MemberCodecError):
            HLT17MemberEncoding(encoding.descriptor, encoding.payload[:-8])
        with self.assertRaises(HLT17MemberCodecError):
            replace(encoding, payload=encoding.payload + b"\x00")
        self.assertEqual(member.physical_state_sha256(), before)
        self.assertEqual(np.asarray(member.tracers.proper_times).tobytes(), labels)

    def test_accepted_plan_and_monitor_must_match_the_fine_successor(self) -> None:
        member = _make_member()
        origin = _origin_predecessor(member)
        origin_encoding = encode_member(member, origin, generation=0)
        prior = _make_member()
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        member.commit_prepared(prepared)
        self.assertEqual(member.time, 1.46875)
        accepted = encode_member(
            member, _evolved_predecessor(member, origin_encoding, prior), generation=1
        )
        metadata = load_canonical_json(accepted.descriptor)
        self.assertEqual(metadata["snapshot_role"], SNAPSHOT_ROLE_ACCEPTED)
        self.assertEqual(metadata["runtime_identity"]["accepted_generation"], 1)
        self.assertEqual(metadata["monitor"]["last_accepted_time_hex"], member.time.hex())
        before = member.physical_state_sha256()
        wrong_plan = plan_forward_proto14_subdivision(
            1.5, 1.5625, 1.0 / 32, minimum_width=float(TDG6_MINIMUM_MACRO_STEP)
        )
        self.assertEqual(wrong_plan.boundaries[4], 1.53125)
        metadata["executed_plan"] = encode_plan(wrong_plan)
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), accepted.payload)
        metadata = load_canonical_json(accepted.descriptor)
        metadata["executed_plan"]["current_hex"] = (1.5).hex()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), accepted.payload)
        metadata = load_canonical_json(accepted.descriptor)
        pred = dict(metadata["predecessor_cursor_identity"])
        pred["previous_step_index"] = 1
        metadata["predecessor_cursor_identity"] = pred
        metadata["predecessor_cursor_identity_sha256"] = sha256(
            canonical_json_bytes(pred)
        ).hexdigest()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), accepted.payload)
        metadata = load_canonical_json(accepted.descriptor)
        pred = dict(metadata["predecessor_cursor_identity"])
        pred["previous_transaction_serial"] = 1
        metadata["predecessor_cursor_identity"] = pred
        metadata["predecessor_cursor_identity_sha256"] = sha256(
            canonical_json_bytes(pred)
        ).hexdigest()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), accepted.payload)
        metadata = load_canonical_json(accepted.descriptor)
        metadata["monitor"]["last_accepted_time_hex"] = (0.0).hex()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(metadata), accepted.payload)
        decoded = decode_member(accepted.descriptor, accepted.payload)
        with self.assertRaises(HLT17MemberCodecError):
            restore_member(
                member,
                _copy_decoded(decoded, executed_plan=wrong_plan),
            )
        self.assertEqual(member.physical_state_sha256(), before)
        ssp = _make_member(COMPARATOR_METHOD)
        ssp_origin = encode_member(ssp, _origin_predecessor(ssp), generation=0)
        ssp_prior = _make_member(COMPARATOR_METHOD)
        ssp.commit_prepared(
            ssp.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        )
        ssp_accepted = encode_member(
            ssp, _evolved_predecessor(ssp, ssp_origin, ssp_prior), generation=1
        )
        ssp_metadata = load_canonical_json(ssp_accepted.descriptor)
        self.assertEqual(ssp_metadata["counters"]["transaction_serial"], 16)
        pred = dict(ssp_metadata["predecessor_cursor_identity"])
        pred["previous_transaction_serial"] = 1
        ssp_metadata["predecessor_cursor_identity"] = pred
        ssp_metadata["predecessor_cursor_identity_sha256"] = sha256(
            canonical_json_bytes(pred)
        ).hexdigest()
        with self.assertRaises(HLT17MemberCodecError):
            decode_member(canonical_json_bytes(ssp_metadata), ssp_accepted.payload)

    def test_accepted_then_retry_overlay_preserves_boundary_and_revisions(self) -> None:
        fail_next = [False]
        boundary = [START]

        def gated(time, current):
            du = np.zeros_like(current.u)
            if fail_next[0] and time == boundary[0] + 1.0 / 128:
                du[:, 4] = 1.0
            return EvolutionRHS(
                du, np.zeros_like(du), np.zeros_like(du), synthetic_diagnostics()
            )

        member = _make_member(rhs=gated)
        origin = encode_member(member, _origin_predecessor(member), generation=0)
        prior = _make_member(rhs=gated)
        prepared = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        self.assertTrue(prepared.assessment.admission_passed)
        member.commit_prepared(prepared)
        boundary[0] = member.time
        fail_next[0] = True
        accepted_time = member.time.hex()
        accepted_hash = member.physical_state_sha256()
        accepted_step = member.step_index
        accepted_serial = member.transaction_serial
        accepted = encode_member(
            member, _evolved_predecessor(member, origin, prior), generation=1
        )
        fail_next[0] = True
        rejected = member.prepare_fresh(event_target=TARGET, requested_cap=1.0 / 32)
        written: list[object] = []
        with self.assertRaises(C1R1TemporalRetryRequired) as caught:
            member.require_admission(rejected, durable_rejection_sink=written.append)
        member.adopt_retry_ledger(caught.exception)
        self.assertEqual(member.time.hex(), accepted_time)
        self.assertEqual(member.physical_state_sha256(), accepted_hash)
        self.assertEqual(member.step_index, accepted_step)
        self.assertEqual(member.transaction_serial, accepted_serial)
        overlay = encode_retry_overlay(
            member,
            _same_boundary_predecessor(member, accepted, accepted_generation=1),
            accepted_generation=1,
        )
        metadata = load_canonical_json(overlay.descriptor)
        self.assertEqual(metadata["snapshot_role"], SNAPSHOT_ROLE_RETRY)
        self.assertEqual(metadata["runtime_identity"]["accepted_generation"], 1)
        self.assertEqual(metadata["accepted_time_hex"], accepted_time)
        self.assertEqual(metadata["physical_state_sha256"], accepted_hash)
        self.assertEqual(
            metadata["imp1_checkpoint"]["current_macro_step_temporal_retry_count"], 1
        )
        self.assertEqual(metadata["imp1_checkpoint"]["accepted_macro_step_count"], 1)
        self.assertIsNone(metadata["executed_plan"])
        twin = _make_member(rhs=gated)
        restore_member(twin, overlay)
        self.assertEqual(twin.time.hex(), accepted_time)
        self.assertEqual(twin.physical_state_sha256(), accepted_hash)
        self.assertEqual(twin.step_index, accepted_step)
        self.assertEqual(twin.temporal_ledger.current_macro_step_temporal_retry_count, 1)
        again = encode_retry_overlay(
            twin,
            _same_boundary_predecessor(twin, accepted, accepted_generation=1),
            accepted_generation=1,
        )
        self.assertEqual(again.descriptor, overlay.descriptor)
        self.assertEqual(again.payload, overlay.payload)
        successor = replan_c1r1_retry(rejected, caught.exception)
        nxt = twin.prepare_retry(successor)
        self.assertTrue(nxt.assessment.admission_passed)


if __name__ == "__main__":
    unittest.main()
