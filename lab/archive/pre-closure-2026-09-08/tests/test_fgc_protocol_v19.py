"""Pure PROTO19 first-event journal/checkpoint schema controls. No store I/O."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes, load_canonical_json  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_admission_runtime import (  # noqa: E402
    hlt17_implementation_identity,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    prepare_hlt17_next,
    replay_hlt17_rejection,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17_MEMBER_KEYS,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SYNTHETIC_MEMBER_KEYS,
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    ACCEPTED_FINE,
    HLT17InMemoryCheckpoint,
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
)
from recursive_horizons.fgc.evolution.protocol_v19 import (  # noqa: E402
    ACCEPTED_SUBSTEPS,
    ATTEMPT_PAYLOAD_KEYS,
    CHECKPOINT_KEYS,
    EVENT_ORIGIN_HEX,
    EVENT_TARGET_HEX,
    FORBIDDEN_JOURNAL_KEYS,
    GENESIS_DIGEST,
    HLT17ProtocolError,
    JOURNAL_ENVELOPE_KEYS,
    MEMBER_STATE_KEYS,
    PRODUCTION_MEMBER_IDENTITIES,
    PROTOCOL_ARTIFACT_ID,
    RECORD_KINDS,
    RK4_STAGE_RECORDS,
    ROOT_MANIFEST_KEYS,
    SSPRK3_STAGE_RECORDS,
    TERMINAL_KINDS,
    authenticate_generation_chain,
    build_root_manifest,
    build_seed_generation,
    build_successor_generation,
    member_state_from_bundle,
    parse_validated_generation,
    production_member_identities,
    revalidate_bundle,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_COMPLETE_STATE_CHANNELS,
    IMP1_REJECTION_EVENT_TYPE,
    append_imp1_rejection,
    imp1_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)


REQUESTED_CAP = 1.0 / 32.0
CAMPAIGN_ID = "HLT17-PROTO19-SYNTHETIC"


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


def _identity(method: str) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_SYNTHETIC,
        method_label=label,
        point_count=9,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id=CAMPAIGN_ID,
        amplitude="synthetic",
    )


def _make_member(method: str = PRIMARY_METHOD, *, rhs=synthetic_exponential_rhs) -> HLT17RuntimeMember:
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    return HLT17RuntimeMember(
        identity=_identity(method),
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


def _root():
    return build_root_manifest(
        campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
    )


def _minimum_step_exhaustion_evidence(bundle) -> dict[str, object]:
    cursor = restore_hlt17_cursor(load_canonical_json(bundle.cursor_bytes))
    current = cursor.ledger
    attempted = float(TDG6_MINIMUM_MACRO_STEP) * 1.5
    record = {
        "event_type": IMP1_REJECTION_EVENT_TYPE,
        "method": current.method,
        "time_hex": current.last_accepted_time.hex(),
        "event_target_hex": cursor.event_target.hex(),
        "attempted_width_hex": attempted.hex(),
        "next_cap_hex": (attempted / 2.0).hex(),
        "accepted_state_sha256": current.current_state_sha256,
        "accepted_step_index": current.current_step_index,
        "accepted_transaction_serial": current.current_transaction_serial,
        "assessment_sha256": "11" * 32,
        "preparation_sha256": "22" * 32,
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "failed_channels": [IMP1_COMPLETE_STATE_CHANNELS[0]],
        "retry_count_for_current_macro_step": current.current_macro_step_temporal_retry_count + 1,
        "cumulative_temporal_retry_count": current.cumulative_temporal_retry_count + 1,
    }
    updated = append_imp1_rejection(current, record)
    mapping = imp1_checkpoint_extension(updated)
    return {
        "kind": "temporal_retry_exhausted",
        "reason": "minimum_macro_step",
        "updated_imp1_ledger": mapping,
        "updated_ledger_sha256": sha256(canonical_json_bytes(mapping)).hexdigest(),
        "executable_retry_cursor": False,
    }


def _seed_live(inherited_source: int = 0, inherited_cfl: int = 0, *, controllable: bool = False):
    members = {}
    checkpoints = {}
    for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
        rhs = ControllableRHS() if controllable else synthetic_exponential_rhs
        member = _make_member(method, rhs=rhs)
        member.source_retry_count = inherited_source
        member.CFL_retry_count = inherited_cfl
        member.__post_init__()
        checkpoint = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
        checkpoints[member.key] = checkpoint
        members[member.key] = checkpoint.generation_bundle()
    root = _root()
    generation = build_seed_generation(root=root, bundles=members)
    return root, generation, checkpoints


def _rebind(files: dict[str, bytes], **replacements: bytes) -> dict[str, bytes]:
    updated = {path: bytes(payload) for path, payload in files.items()}
    updated.update(replacements)
    leaves = {
        path: sha256(payload).hexdigest()
        for path, payload in updated.items()
        if path != "generation_manifest.json"
    }
    manifest = load_canonical_json(updated["generation_manifest.json"])
    manifest["leaves"] = leaves
    if "checkpoint.json" in replacements:
        manifest["checkpoint_sha256"] = sha256(updated["checkpoint.json"]).hexdigest()
    if "journal.json" in replacements:
        digest = sha256(updated["journal.json"]).hexdigest()
        manifest["journal_record_sha256"] = digest
        checkpoint = load_canonical_json(updated["checkpoint.json"])
        checkpoint["journal_record_sha256"] = digest
        updated["checkpoint.json"] = canonical_json_bytes(checkpoint)
        manifest["checkpoint_sha256"] = sha256(updated["checkpoint.json"]).hexdigest()
        manifest["leaves"]["checkpoint.json"] = manifest["checkpoint_sha256"]
        manifest["leaves"]["journal.json"] = digest
    updated["generation_manifest.json"] = canonical_json_bytes(manifest)
    return updated


class ProtocolV19Tests(unittest.TestCase):
    def test_production_metadata_identities_are_exact_and_not_nine_point(self) -> None:
        identities = production_member_identities()
        self.assertEqual(identities, PRODUCTION_MEMBER_IDENTITIES)
        self.assertEqual(
            identities,
            (
                ("RK4-2049", 2049, 2044, 16352, 32704),
                ("RK4-4097", 4097, 4092, 32736, 65472),
                ("RK4-8193", 8193, 8188, 65504, 131008),
                ("SSPRK3-4097", 4097, 4092, 32736, 65472),
                ("SSPRK3-8193", 8193, 8188, 65504, 131008),
                ("SSPRK3-16385", 16385, 16380, 131040, 262080),
            ),
        )
        self.assertEqual(tuple(item.member_key for item in HLT17_PRODUCTION_MEMBERS), HLT17_MEMBER_KEYS)
        self.assertEqual(HLT17_SYNTHETIC_MEMBER_KEYS, ("RK4-9", "SSPRK3-9"))
        self.assertEqual(ACCEPTED_SUBSTEPS, 4)
        self.assertEqual(RK4_STAGE_RECORDS, 20)
        self.assertEqual(SSPRK3_STAGE_RECORDS, 16)
        root = build_root_manifest(campaign_id=CAMPAIGN_ID, member_keys=HLT17_MEMBER_KEYS)
        self.assertEqual(tuple(root["member_keys"]), HLT17_MEMBER_KEYS)
        self.assertFalse(root["production_write_authorized"])
        self.assertFalse(root["historical_origin_authenticated"])
        with self.assertRaisesRegex(HLT17ProtocolError, "synthetic RK4-9/SSPRK3-9"):
            build_seed_generation(root=root, bundles={})

    def test_seed_binds_distinct_hashes_and_genesis_parents(self) -> None:
        root, generation, checkpoints = _seed_live(inherited_source=3, inherited_cfl=1)
        self.assertEqual(generation.kind, "seed")
        self.assertEqual(generation.store_generation, 0)
        self.assertEqual(generation.journal_sequence, 0)
        journal = generation.journal_mapping
        checkpoint = generation.checkpoint_mapping
        self.assertEqual(journal["journal_parent_sha256"], GENESIS_DIGEST)
        self.assertEqual(checkpoint["predecessor_checkpoint_sha256"], GENESIS_DIGEST)
        self.assertEqual(checkpoint["journal_record_sha256"], generation.journal_sha256)
        self.assertNotIn("checkpoint_sha256", journal)
        for forbidden in FORBIDDEN_JOURNAL_KEYS:
            self.assertNotIn(forbidden, journal)
        self.assertNotEqual(generation.checkpoint_sha256, generation.journal_sha256)
        self.assertEqual(checkpoint["event_origin_time_hex"], EVENT_ORIGIN_HEX)
        self.assertEqual(checkpoint["event_target_time_hex"], EVENT_TARGET_HEX)
        self.assertFalse(generation.terminal)
        self.assertFalse(generation.first_event_complete)
        for key, checkpoint_obj in checkpoints.items():
            state = member_state_from_bundle(generation.bundles[key])
            for name, value in hlt17_implementation_identity().items():
                self.assertIn(name, MEMBER_STATE_KEYS)
                self.assertEqual(state[name], value)
            self.assertEqual(state["source_total"], 3)
            self.assertEqual(state["cfl_total"], 1)
            self.assertEqual(state["source_current"], 0)
            self.assertEqual(state["accepted_generation"], 0)
            self.assertNotEqual(state["cursor_sha256"], generation.checkpoint_sha256)
            self.assertNotEqual(state["cursor_sha256"], generation.journal_sha256)
            self.assertNotEqual(state["kernel_transition_digest"], generation.journal_sha256)
            self.assertNotEqual(state["kernel_transition_digest"], state["cursor_sha256"])
            self.assertEqual(state["accepted_time_hex"], START.hex())
            self.assertEqual(checkpoint_obj.cursor.source_total, 3)

    def test_accepted_fine_increments_generation_and_method_records(self) -> None:
        root, seed, checkpoints = _seed_live()
        for method, key, serial in (
            (PRIMARY_METHOD, "RK4-9", 20),
            (COMPARATOR_METHOD, "SSPRK3-9", 16),
        ):
            with self.subTest(member=key):
                result = checkpoints[key].advance_once(requested_cap=REQUESTED_CAP)
                self.assertEqual(result.disposition, ACCEPTED_FINE)
                successor = build_successor_generation(
                    root=root,
                    predecessor=seed,
                    predecessor_bundles=seed.bundles,
                    kind="accepted_fine",
                    member_key=key,
                    successor_bundle=result.bundle,
                )
                self.assertEqual(successor.store_generation, seed.store_generation + 1)
                self.assertEqual(successor.journal_sequence, seed.journal_sequence + 1)
                self.assertEqual(successor.journal_mapping["journal_parent_sha256"], seed.journal_sha256)
                self.assertEqual(
                    successor.checkpoint_mapping["predecessor_checkpoint_sha256"],
                    seed.checkpoint_sha256,
                )
                previous = member_state_from_bundle(seed.bundles[key])
                current = member_state_from_bundle(successor.bundles[key])
                self.assertEqual(current["accepted_generation"], previous["accepted_generation"] + 1)
                self.assertEqual(
                    current["kernel_checkpoint_generation"],
                    previous["kernel_checkpoint_generation"] + 1,
                )
                self.assertEqual(current["step_index"], previous["step_index"] + 4)
                self.assertEqual(current["transaction_serial"], previous["transaction_serial"] + serial)
                self.assertEqual(current["source_total"], previous["source_total"])
                self.assertNotEqual(current["descriptor_sha256"], previous["descriptor_sha256"])
                self.assertNotEqual(current["cursor_sha256"], previous["cursor_sha256"])
                self.assertNotEqual(
                    current["kernel_transition_digest"], successor.journal_sha256
                )
                other = "SSPRK3-9" if key == "RK4-9" else "RK4-9"
                self.assertEqual(
                    successor.bundles[other].cursor_bytes, seed.bundles[other].cursor_bytes
                )
                del method

    def test_retries_preserve_accepted_state_and_tdg6_and_advance_own_counters(self) -> None:
        root, seed, checkpoints = _seed_live(controllable=True)
        key = "RK4-9"
        checkpoint = checkpoints[key]
        rhs = checkpoint.member.operator
        assert isinstance(rhs, ControllableRHS)
        first = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(first.disposition, ACCEPTED_FINE)
        accepted = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="accepted_fine",
            member_key=key,
            successor_bundle=first.bundle,
        )
        accepted_state = member_state_from_bundle(accepted.bundles[key])
        tdg6 = accepted_state["tdg6_snapshot_sha256"]
        rhs.impulse_from = checkpoint.member.time
        rhs.mode = "impulse"
        temporal = checkpoint.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(temporal.disposition, "temporal_retry_required")
        pending = build_successor_generation(
            root=root,
            predecessor=accepted,
            predecessor_bundles=accepted.bundles,
            kind="temporal_retry",
            member_key=key,
            successor_bundle=temporal.bundle,
        )
        temporal_state = member_state_from_bundle(pending.bundles[key])
        self.assertEqual(temporal_state["accepted_generation"], accepted_state["accepted_generation"])
        self.assertEqual(temporal_state["descriptor_sha256"], accepted_state["descriptor_sha256"])
        self.assertEqual(temporal_state["physical_state_sha256"], accepted_state["physical_state_sha256"])
        self.assertEqual(temporal_state["accepted_time_hex"], accepted_state["accepted_time_hex"])
        self.assertEqual(temporal_state["tdg6_snapshot_sha256"], tdg6)
        self.assertEqual(temporal_state["cursor_generation"], accepted_state["cursor_generation"] + 1)
        self.assertTrue(pending.journal_mapping["payload"]["accepted_state_preserved"])
        self.assertTrue(pending.journal_mapping["payload"]["executable_retry_cursor"])
        self.assertNotEqual(pending.store_generation, temporal_state["kernel_checkpoint_generation"])
        checkpoint.restore(checkpoint.accepted_encoding, encode_hlt17_cursor(checkpoint.cursor))
        rhs.impulse_from = checkpoint.member.time
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(checkpoint.cursor, checkpoint.live_boundary())
        rhs.mode = "source"
        source = checkpoint.absorb_nonfine(
            prepare_hlt17_next(checkpoint.cursor, checkpoint.live_boundary(), replay=replay)
        )
        self.assertEqual(source.disposition, "source_retry_required")
        sourced = build_successor_generation(
            root=root,
            predecessor=pending,
            predecessor_bundles=pending.bundles,
            kind="source_retry",
            member_key=key,
            successor_bundle=source.bundle,
        )
        source_state = member_state_from_bundle(sourced.bundles[key])
        self.assertEqual(source_state["source_current"], temporal_state["source_current"] + 1)
        self.assertEqual(source_state["source_total"], temporal_state["source_total"] + 1)
        self.assertEqual(source_state["cfl_total"], temporal_state["cfl_total"])
        self.assertEqual(source_state["descriptor_sha256"], accepted_state["descriptor_sha256"])
        self.assertEqual(source_state["tdg6_snapshot_sha256"], tdg6)
        checkpoint.restore(checkpoint.accepted_encoding, encode_hlt17_cursor(checkpoint.cursor))
        rhs.impulse_from = checkpoint.member.time
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(checkpoint.cursor, checkpoint.live_boundary())
        rhs.mode = "cfl"
        cfl = checkpoint.absorb_nonfine(
            prepare_hlt17_next(checkpoint.cursor, checkpoint.live_boundary(), replay=replay)
        )
        self.assertEqual(cfl.disposition, "cfl_retry_required")
        cfl_gen = build_successor_generation(
            root=root,
            predecessor=sourced,
            predecessor_bundles=sourced.bundles,
            kind="cfl_retry",
            member_key=key,
            successor_bundle=cfl.bundle,
        )
        cfl_state = member_state_from_bundle(cfl_gen.bundles[key])
        self.assertEqual(cfl_state["cfl_current"], source_state["cfl_current"] + 1)
        self.assertEqual(cfl_state["source_total"], source_state["source_total"])
        self.assertEqual(cfl_state["tdg6_snapshot_sha256"], tdg6)
        self.assertEqual(cfl_state["physical_state_sha256"], accepted_state["physical_state_sha256"])

    def test_typed_terminals_preserve_history_and_forbid_fabricated_retry_cursors(self) -> None:
        root, seed, _checkpoints = _seed_live()
        key = "RK4-9"
        descriptor = load_canonical_json(seed.bundles[key].descriptor)
        ledger = descriptor["imp1_checkpoint"]
        premature = {
            "kind": "temporal_retry_exhausted",
            "reason": "maximum_temporal_retries",
            "updated_imp1_ledger": ledger,
            "updated_ledger_sha256": sha256(canonical_json_bytes(ledger)).hexdigest(),
            "executable_retry_cursor": False,
        }
        with self.assertRaisesRegex(HLT17ProtocolError, "exact one-rejection successor|premature"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
                kind="temporal_retry_exhausted",
                member_key=key,
                successor_bundle=seed.bundles[key],
                terminal_evidence=premature,
            )
        foreign = _minimum_step_exhaustion_evidence(seed.bundles["SSPRK3-9"])
        with self.assertRaisesRegex(HLT17ProtocolError, "method differs|origin|successor"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
                kind="temporal_retry_exhausted",
                member_key=key,
                successor_bundle=seed.bundles[key],
                terminal_evidence=foreign,
            )
        evidence = _minimum_step_exhaustion_evidence(seed.bundles[key])
        terminal = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="temporal_retry_exhausted",
            member_key=key,
            successor_bundle=seed.bundles[key],
            terminal_evidence=evidence,
        )
        self.assertTrue(terminal.terminal)
        self.assertEqual(terminal.disposition, "temporal_retry_exhausted")
        self.assertFalse(terminal.first_event_complete)
        self.assertEqual(
            terminal.bundles[key].cursor_bytes, seed.bundles[key].cursor_bytes
        )
        self.assertFalse(terminal.journal_mapping["payload"]["executable_retry_cursor"])
        self.assertIsNotNone(terminal.journal_mapping["payload"]["terminal_evidence"])
        with self.assertRaisesRegex(HLT17ProtocolError, "forbids further appends"):
            build_successor_generation(
                root=root,
                predecessor=terminal,
                predecessor_bundles=terminal.bundles,
                kind="accepted_fine",
                member_key=key,
                successor_bundle=seed.bundles[key],
            )
        invalid = {
            "kind": "invalid_premise",
            "reason": "synthetic invalid premise",
            "spends_retry": False,
            "executable_retry_cursor": False,
        }
        root2, seed2, _c2 = _seed_live()
        built = build_successor_generation(
            root=root2,
            predecessor=seed2,
            predecessor_bundles=seed2.bundles,
            kind="invalid_premise",
            member_key=key,
            successor_bundle=seed2.bundles[key],
            terminal_evidence=invalid,
        )
        self.assertTrue(built.terminal)
        resource = {
            "kind": "resource_exhausted",
            "reason": "synthetic resource stop without a production execution",
            "executable_retry_cursor": False,
        }
        root3, seed3, _c3 = _seed_live()
        exhausted = build_successor_generation(
            root=root3,
            predecessor=seed3,
            predecessor_bundles=seed3.bundles,
            kind="resource_exhausted",
            member_key=key,
            successor_bundle=seed3.bundles[key],
            terminal_evidence=resource,
        )
        self.assertTrue(exhausted.terminal)
        self.assertFalse(exhausted.first_event_complete)

    def test_first_event_completion_requires_every_member_fresh_at_target(self) -> None:
        root, generation, checkpoints = _seed_live()
        order = ("RK4-9", "SSPRK3-9", "RK4-9", "SSPRK3-9")
        for index, key in enumerate(order):
            result = checkpoints[key].advance_once(requested_cap=REQUESTED_CAP)
            self.assertEqual(result.disposition, ACCEPTED_FINE)
            generation = build_successor_generation(
                root=root,
                predecessor=generation,
                predecessor_bundles=generation.bundles,
                kind="accepted_fine",
                member_key=key,
                successor_bundle=result.bundle,
            )
            if index < 3:
                self.assertFalse(generation.first_event_complete)
                self.assertEqual(generation.disposition, "nonterminal")
        self.assertTrue(generation.first_event_complete)
        self.assertTrue(generation.terminal)
        self.assertEqual(generation.disposition, "first_event_complete")
        for key in HLT17_SYNTHETIC_MEMBER_KEYS:
            state = member_state_from_bundle(generation.bundles[key])
            self.assertEqual(state["accepted_time_hex"], TARGET.hex())
            self.assertTrue(state["parked"])
            self.assertTrue(state["public_fresh"])
        with self.assertRaisesRegex(HLT17ProtocolError, "forbids further appends|final boundary"):
            build_successor_generation(
                root=root,
                predecessor=generation,
                predecessor_bundles=generation.bundles,
                kind="accepted_fine",
                member_key="RK4-9",
                successor_bundle=generation.bundles["RK4-9"],
            )
        with self.assertRaisesRegex(HLT17ProtocolError, "disagrees with the immutable"):
            replace(generation, terminal=False)
        with self.assertRaisesRegex(HLT17ProtocolError, "final boundary|forbids further"):
            build_successor_generation(
                root=root,
                predecessor=generation,
                predecessor_bundles=generation.bundles,
                kind="resource_exhausted",
                member_key="RK4-9",
                successor_bundle=generation.bundles["RK4-9"],
                terminal_evidence={
                    "kind": "resource_exhausted",
                    "reason": "cannot overwrite first-event completion",
                    "executable_retry_cursor": False,
                },
            )

    def test_schema_mutations_and_circular_successor_hashes_fail(self) -> None:
        root, seed, checkpoints = _seed_live()
        result = checkpoints["RK4-9"].advance_once(requested_cap=REQUESTED_CAP)
        successor = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="accepted_fine",
            member_key="RK4-9",
            successor_bundle=result.bundle,
        )
        chain = authenticate_generation_chain(
            root=root, generations=(dict(seed.files), dict(successor.files))
        )
        self.assertEqual(len(chain), 2)
        checkpoint = load_canonical_json(successor.checkpoint_raw)
        self.assertEqual(set(checkpoint), set(CHECKPOINT_KEYS))
        journal = load_canonical_json(successor.journal_raw)
        self.assertEqual(set(journal), set(JOURNAL_ENVELOPE_KEYS))
        self.assertEqual(set(journal["payload"]), set(ATTEMPT_PAYLOAD_KEYS))
        self.assertEqual(set(root), set(ROOT_MANIFEST_KEYS))
        mutated = dict(checkpoint)
        mutated["store_generation"] = 99
        files = _rebind(
            dict(successor.files),
            **{"checkpoint.json": canonical_json_bytes(mutated)},
        )
        with self.assertRaisesRegex(HLT17ProtocolError, "store generation"):
            parse_validated_generation(
                files, root=root, predecessor=seed, predecessor_bundles=seed.bundles
            )
        mutated = dict(journal)
        mutated["journal_parent_sha256"] = successor.journal_sha256
        files = _rebind(
            dict(successor.files),
            **{"journal.json": canonical_json_bytes(mutated)},
        )
        with self.assertRaisesRegex(HLT17ProtocolError, "journal parent"):
            parse_validated_generation(
                files, root=root, predecessor=seed, predecessor_bundles=seed.bundles
            )
        mutated = dict(journal)
        mutated["checkpoint_sha256"] = successor.checkpoint_sha256
        with self.assertRaisesRegex(HLT17ProtocolError, "successor-checkpoint|keys differ"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"journal.json": canonical_json_bytes(mutated)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        mutated = dict(checkpoint)
        mutated["first_event_complete"] = True
        mutated["disposition"] = "first_event_complete"
        mutated["terminal"] = True
        with self.assertRaisesRegex(HLT17ProtocolError, "first-event completion"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"checkpoint.json": canonical_json_bytes(mutated)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        mutated = dict(checkpoint)
        mutated["campaign_execution_authorized"] = True
        with self.assertRaisesRegex(HLT17ProtocolError, "campaign_execution_authorized"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"checkpoint.json": canonical_json_bytes(mutated)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        mutated = dict(checkpoint)
        members = dict(mutated["members"])
        other = dict(members["SSPRK3-9"])
        other["accepted_generation"] = 3
        members["SSPRK3-9"] = other
        mutated["members"] = members
        with self.assertRaisesRegex(HLT17ProtocolError, "unchanged member|disagree"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"checkpoint.json": canonical_json_bytes(mutated)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        mutated = dict(journal)
        mutated["journal_sequence"] = True
        with self.assertRaisesRegex(HLT17ProtocolError, "journal_sequence"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"journal.json": canonical_json_bytes(mutated)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        hole = authenticate_generation_chain
        with self.assertRaisesRegex(HLT17ProtocolError, "hole|predecessor|seed cannot"):
            hole(root=root, generations=(dict(successor.files),))
        swapped = dict(checkpoint)
        swapped["journal_record_sha256"] = successor.checkpoint_sha256
        with self.assertRaisesRegex(HLT17ProtocolError, "journal"):
            parse_validated_generation(
                _rebind(
                    dict(successor.files),
                    **{"checkpoint.json": canonical_json_bytes(swapped)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        for kind in RECORD_KINDS:
            self.assertIsInstance(kind, str)
        for kind in TERMINAL_KINDS:
            self.assertIn(kind, RECORD_KINDS)
        rebuilt = revalidate_bundle(successor.bundles["RK4-9"])
        self.assertEqual(rebuilt.descriptor, successor.bundles["RK4-9"].descriptor)
        for key in CHECKPOINT_KEYS:
            mutated = dict(checkpoint)
            original = mutated[key]
            if original is False:
                mutated[key] = True
            elif original is True:
                mutated[key] = False
            elif type(original) is int:
                mutated[key] = original + 11
            elif type(original) is str and len(original) == 64:
                mutated[key] = "cd" * 32
            elif key == "kind":
                mutated[key] = "seed"
            elif key == "disposition":
                mutated[key] = "temporal_retry_exhausted"
            elif key == "changed_member_key":
                mutated[key] = "SSPRK3-9"
            elif key in {"member_keys", "members"}:
                continue
            else:
                continue
            with self.subTest(checkpoint_field=key):
                with self.assertRaises(HLT17ProtocolError):
                    parse_validated_generation(
                        _rebind(
                            dict(successor.files),
                            **{"checkpoint.json": canonical_json_bytes(mutated)},
                        ),
                        root=root,
                        predecessor=seed,
                        predecessor_bundles=seed.bundles,
                    )
        for key in JOURNAL_ENVELOPE_KEYS:
            mutated = dict(journal)
            original = mutated[key]
            if original is False:
                mutated[key] = True
            elif type(original) is int:
                mutated[key] = original + 5
            elif type(original) is str and len(original) == 64:
                mutated[key] = "ef" * 32
            elif key == "kind":
                mutated[key] = "seed"
            elif key == "payload":
                mutated[key] = {}
            else:
                continue
            with self.subTest(journal_field=key):
                with self.assertRaises(HLT17ProtocolError):
                    parse_validated_generation(
                        _rebind(
                            dict(successor.files),
                            **{"journal.json": canonical_json_bytes(mutated)},
                        ),
                        root=root,
                        predecessor=seed,
                        predecessor_bundles=seed.bundles,
                    )
        member_state = dict(checkpoint["members"]["RK4-9"])
        self.assertEqual(set(member_state), set(MEMBER_STATE_KEYS))
        for key in MEMBER_STATE_KEYS:
            mutated_state = dict(member_state)
            original = mutated_state[key]
            if original is False:
                mutated_state[key] = True
            elif original is True:
                mutated_state[key] = False
            elif type(original) is int:
                mutated_state[key] = original + 3
            elif type(original) is str and len(original) == 64:
                mutated_state[key] = "aa" * 32
            elif key in {"mode", "latest_owner", "accepted_time_hex"}:
                mutated_state[key] = "FRESH_READY" if original != "FRESH_READY" else "RETRY_PENDING"
            else:
                continue
            mutated = dict(checkpoint)
            members = dict(mutated["members"])
            members["RK4-9"] = mutated_state
            mutated["members"] = members
            with self.subTest(member_field=key):
                with self.assertRaises(HLT17ProtocolError):
                    parse_validated_generation(
                        _rebind(
                            dict(successor.files),
                            **{"checkpoint.json": canonical_json_bytes(mutated)},
                        ),
                        root=root,
                        predecessor=seed,
                        predecessor_bundles=seed.bundles,
                    )

    def test_unknown_keys_and_two_member_changes_are_rejected(self) -> None:
        root, seed, checkpoints = _seed_live()
        first = checkpoints["RK4-9"].advance_once(requested_cap=REQUESTED_CAP)
        second = checkpoints["SSPRK3-9"].advance_once(requested_cap=REQUESTED_CAP)
        valid = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="accepted_fine",
            member_key="RK4-9",
            successor_bundle=first.bundle,
        )
        extra = dict(valid.checkpoint_mapping)
        extra["latest"] = valid.checkpoint_sha256
        with self.assertRaisesRegex(HLT17ProtocolError, "latest|keys differ"):
            parse_validated_generation(
                _rebind(
                    dict(valid.files),
                    **{"checkpoint.json": canonical_json_bytes(extra)},
                ),
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
            )
        mixed = dict(seed.bundles)
        mixed["RK4-9"] = first.bundle
        mixed["SSPRK3-9"] = second.bundle
        with self.assertRaisesRegex(HLT17ProtocolError, "predecessor bundle bytes"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=mixed,
                kind="accepted_fine",
                member_key="RK4-9",
                successor_bundle=first.bundle,
            )
        with self.assertRaisesRegex(HLT17ProtocolError, "exactly one"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
                kind="accepted_fine",
                member_key="RK4-9",
                successor_bundle=seed.bundles["RK4-9"],
            )
        root_text = canonical_json_bytes(root).decode("ascii")
        self.assertIn(PROTOCOL_ARTIFACT_ID, root_text)
        self.assertNotIn("calibration claimed", root_text.lower())
        self.assertFalse(root["calibration_claimed"])
        self.assertFalse(root["eligibility_claimed"])


if __name__ == "__main__":
    unittest.main()
