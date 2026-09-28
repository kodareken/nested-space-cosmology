"""Focused codec and pure-transition controls for the HLT17 IMP1 cursor."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from hashlib import sha256
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    CURSOR_KEYS,
    FRESH_READY,
    HLT17IMP1Cursor,
    HLT17IMP1Error,
    HLT17OwnerRetryExhausted,
    HLT17_CFL_RETRY_CAP,
    HLT17_CURSOR_SCHEMA,
    HLT17_IMPLEMENTATION_ID,
    HLT17_LEGACY_IMP1_CURSOR_SCHEMA,
    HLT17_MEMBER_KEYS,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SYNTHETIC_MEMBER_KEYS,
    IDENTITY_SCOPE_SYNTHETIC,
    HLT17_SOURCE_RETRY_CAP,
    RETRY_PENDING,
    SOURCE_OWNER,
    close_fine_adoption,
    close_source_cfl_overlay,
    close_temporal_rejection,
    cursor_sha256,
    decode_tdg7_plan,
    encode_hlt17_cursor,
    encode_tdg7_plan,
    imp1_enclosure_limits_for_member,
    ledger_sha256,
    plan_imp1_lattice,
    restore_hlt17_cursor,
    seed_hlt17_cursor,
    strictly_smaller_tdg7_plan,
    temporal_successor_plan,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_COMPLETE_STATE_CHANNELS,
    IMP1_REJECTION_EVENT_TYPE,
    accept_imp1_step,
    append_imp1_rejection,
    imp1_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START,
    SYNTHETIC_TARGET,
    make_synthetic_fixture,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (  # noqa: E402
    cal11_event24_witness,
)


MODULE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_imp1_cursor.py"


def _sha(tag: str) -> str:
    return sha256(tag.encode("ascii")).hexdigest()


def _identities(physical: str) -> dict[str, str]:
    names = (
        "descriptor",
        "coordinates",
        "source",
        "monitor",
        "causal",
        "tracer",
        "guard",
    )
    values = {name: _sha(name) for name in names}
    if physical in values.values():
        values["descriptor"] = _sha("descriptor-distinct")
    return values


def _seed_cursor(ledger=None):
    if ledger is None:
        ledger = make_synthetic_fixture()[3]
    identities = _identities(ledger.current_state_sha256)
    cursor = seed_hlt17_cursor(
        member_key="RK4-9",
        event_target=SYNTHETIC_TARGET,
        ledger=ledger,
        descriptor_sha256=identities["descriptor"],
        physical_state_sha256=ledger.current_state_sha256,
        coordinates_sha256=identities["coordinates"],
        synthetic_callable_binding_sha256=identities["source"],
        monitor_sha256=identities["monitor"],
        causal_state_sha256=identities["causal"],
        tracer_sha256=identities["tracer"],
        guard_configuration_sha256=identities["guard"],
    )
    return cursor, ledger, identities


def _plan(cap: float = 1.0 / 8.0):
    return plan_imp1_lattice(SYNTHETIC_START, SYNTHETIC_TARGET, cap)


def _rejection(ledger, plan, *, preparation: str | None = None):
    preparation_hash = _sha("preparation") if preparation is None else preparation
    assessment_hash = _sha("assessment")
    return {
        "event_type": IMP1_REJECTION_EVENT_TYPE,
        "method": ledger.method,
        "time_hex": ledger.last_accepted_time.hex(),
        "event_target_hex": plan.event_target.hex(),
        "attempted_width_hex": plan.macro_width.hex(),
        "next_cap_hex": (plan.macro_width / 2.0).hex(),
        "accepted_state_sha256": ledger.current_state_sha256,
        "accepted_step_index": ledger.current_step_index,
        "accepted_transaction_serial": ledger.current_transaction_serial,
        "assessment_sha256": assessment_hash,
        "preparation_sha256": preparation_hash,
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "failed_channels": ["u:alpha"],
        "retry_count_for_current_macro_step": (
            ledger.current_macro_step_temporal_retry_count + 1
        ),
        "cumulative_temporal_retry_count": ledger.cumulative_temporal_retry_count + 1,
    }


def _pending_cursor():
    cursor, ledger, identities = _seed_cursor()
    plan = _plan()
    record = _rejection(ledger, plan)
    updated = append_imp1_rejection(ledger, record)
    pending = close_temporal_rejection(
        cursor,
        predecessor_plan=plan,
        predecessor_ledger=ledger,
        updated_ledger=updated,
        rejection=record,
        preparation_sha256=record["preparation_sha256"],
        assessment_sha256=record["assessment_sha256"],
    )
    return pending, ledger, updated, plan, identities, record


def _source_evidence(plan) -> dict[str, object]:
    return {
        "kind": "source_retry",
        "path": "outer_0",
        "attempted_requested_cap_hex": plan.requested_cap.hex(),
        "attempted_macro_width_hex": plan.macro_width.hex(),
        "source_retry": {"synthetic": True, "retryable_source_only": True},
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }


def _cfl_evidence(plan) -> dict[str, object]:
    return {
        "kind": "cfl_retry",
        "path": "outer_0",
        "attempted_requested_cap_hex": plan.requested_cap.hex(),
        "attempted_macro_width_hex": plan.macro_width.hex(),
        "observed_ratio_hex": (2.0).hex(),
        "maximum_ratio_hex": (1.0).hex(),
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }


class HLT17CursorTests(unittest.TestCase):
    def test_six_production_members_use_exact_n_and_8n_16n_ceilings(self) -> None:
        expected = (
            ("RK4-2049", 2049, 2044, 16352, 32704),
            ("RK4-4097", 4097, 4092, 32736, 65472),
            ("RK4-8193", 8193, 8188, 65504, 131008),
            ("SSPRK3-4097", 4097, 4092, 32736, 65472),
            ("SSPRK3-8193", 8193, 8188, 65504, 131008),
            ("SSPRK3-16385", 16385, 16380, 131040, 262080),
        )
        self.assertEqual(HLT17_MEMBER_KEYS, tuple(item[0] for item in expected))
        self.assertEqual(len(HLT17_PRODUCTION_MEMBERS), 6)
        for member, item in zip(HLT17_PRODUCTION_MEMBERS, expected, strict=True):
            self.assertEqual(member.member_key, item[0])
            self.assertEqual(member.point_count, item[1])
            self.assertEqual(member.owned_rows, item[1] - 5)
            self.assertEqual(member.owned_rows, item[2])
            self.assertEqual(member.fallback_d01, 8 * member.owned_rows)
            self.assertEqual(member.fallback_d12, 16 * member.owned_rows)
            self.assertEqual(member.fallback_d01, item[3])
            self.assertEqual(member.fallback_d12, item[4])
            limits = imp1_enclosure_limits_for_member(member.member_key)
            self.assertEqual(limits.maximum_owned_rows, member.owned_rows)
            self.assertEqual(limits.fallback_maximum_candidates_D01, member.fallback_d01)
            self.assertEqual(limits.fallback_maximum_candidates_D12, member.fallback_d12)
            self.assertEqual(limits.maximum_rational_bits, 32768)
        self.assertEqual(
            (
                HLT17_PRODUCTION_MEMBERS[-1].owned_rows,
                HLT17_PRODUCTION_MEMBERS[-1].fallback_d01,
                HLT17_PRODUCTION_MEMBERS[-1].fallback_d12,
            ),
            (16380, 131040, 262080),
        )
        self.assertEqual(HLT17_SYNTHETIC_MEMBER_KEYS, ("RK4-9", "SSPRK3-9"))
        rk4_9 = imp1_enclosure_limits_for_member("RK4-9")
        self.assertEqual(rk4_9.maximum_owned_rows, 4)
        self.assertEqual(rk4_9.fallback_maximum_candidates_D01, 32)
        self.assertEqual(rk4_9.fallback_maximum_candidates_D12, 64)

    def test_requested_cap_is_not_macro_width_for_cal11_and_target_limited_plans(
        self,
    ) -> None:
        witness = cal11_event24_witness()
        self.assertNotEqual(witness.requested_cap.hex(), witness.macro_width.hex())
        restored = decode_tdg7_plan(encode_tdg7_plan(witness))
        self.assertEqual(restored, witness)
        self.assertEqual(restored.requested_cap.hex(), witness.requested_cap.hex())
        self.assertEqual(len(restored.boundaries), 5)

        plan = _plan(1.0 / 8.0)
        self.assertEqual(plan.requested_cap, 1.0 / 8.0)
        self.assertEqual(plan.macro_width, 1.0 / 16.0)
        self.assertTrue(plan.target_limited)
        width_only = plan_imp1_lattice(
            plan.current, plan.event_target, plan.macro_width
        )
        self.assertNotEqual(encode_tdg7_plan(width_only), encode_tdg7_plan(plan))
        self.assertEqual(width_only.macro_width, plan.macro_width)
        self.assertEqual(width_only.requested_cap, plan.macro_width)

    def test_fresh_cursor_round_trip_and_refuses_origin_authentication(self) -> None:
        cursor, ledger, identities = _seed_cursor()
        mapping = encode_hlt17_cursor(cursor)
        restored = restore_hlt17_cursor(mapping)
        self.assertEqual(encode_hlt17_cursor(restored), mapping)
        self.assertEqual(cursor.mode, FRESH_READY)
        self.assertIsNone(cursor.temporal)
        self.assertIsNone(cursor.overlay)
        self.assertEqual(mapping["schema"], HLT17_CURSOR_SCHEMA)
        self.assertEqual(mapping["implementation_id"], HLT17_IMPLEMENTATION_ID)
        self.assertNotEqual(mapping["schema"], HLT17_LEGACY_IMP1_CURSOR_SCHEMA)
        self.assertFalse(mapping["historical_origin_authenticated"])
        self.assertFalse(mapping["origin_authenticated_by_this_cursor"])
        self.assertFalse(mapping["source_authenticated_by_this_cursor"])
        self.assertFalse(mapping["store_authenticated_by_this_cursor"])
        self.assertFalse(mapping["campaign_execution_authorized"])
        self.assertNotIn("pending_owner", mapping)
        self.assertEqual(set(mapping), set(CURSOR_KEYS))
        self.assertEqual(mapping["identity_scope"], IDENTITY_SCOPE_SYNTHETIC)
        self.assertEqual(mapping["member_key"], "RK4-9")
        self.assertEqual(mapping["point_count"], 9)
        self.assertEqual(mapping["descriptor_sha256"], identities["descriptor"])
        self.assertEqual(mapping["physical_state_sha256"], ledger.current_state_sha256)
        self.assertNotEqual(
            mapping["descriptor_sha256"], mapping["physical_state_sha256"]
        )
        self.assertEqual(mapping["imp1_ledger_sha256"], ledger_sha256(ledger))
        self.assertEqual(
            mapping["inherited_snapshot_sha256"],
            imp1_checkpoint_extension(ledger)["inherited_snapshot_sha256"],
        )
        self.assertTrue(cursor.public_fresh_permitted())

    def test_swapped_descriptor_and_physical_hashes_are_refused(self) -> None:
        cursor, _, _ = _seed_cursor()
        mapping = encode_hlt17_cursor(cursor)
        mapping["descriptor_sha256"], mapping["physical_state_sha256"] = (
            mapping["physical_state_sha256"],
            mapping["descriptor_sha256"],
        )
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)
        with self.assertRaises(HLT17IMP1Error):
            seed_hlt17_cursor(
                member_key="RK4-9",
                event_target=SYNTHETIC_TARGET,
                ledger=cursor.ledger,
                descriptor_sha256=cursor.physical_state_sha256,
                physical_state_sha256=cursor.physical_state_sha256,
                coordinates_sha256=cursor.coordinates_sha256,
                synthetic_callable_binding_sha256=cursor.synthetic_callable_binding_sha256,
                monitor_sha256=cursor.monitor_sha256,
                causal_state_sha256=cursor.causal_state_sha256,
                tracer_sha256=cursor.tracer_sha256,
                guard_configuration_sha256=cursor.guard_configuration_sha256,
            )

    def test_unknown_keys_bool_int_aliases_and_noncanonical_numbers_fail(self) -> None:
        cursor, _, _ = _seed_cursor()
        extra = encode_hlt17_cursor(cursor)
        extra["pending_owner"] = "temporal"
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(extra)

        boolean_step = encode_hlt17_cursor(cursor)
        boolean_step["accepted_step_index"] = True
        with self.assertRaises((HLT17IMP1Error, TypeError)):
            restore_hlt17_cursor(boolean_step)

        alias = encode_hlt17_cursor(cursor)
        alias["historical_origin_authenticated"] = 0
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(alias)

        short_hex = encode_hlt17_cursor(cursor)
        short_hex["event_target_hex"] = "0x1.8p+0"
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(short_hex)

        origin = encode_hlt17_cursor(cursor)
        origin["historical_origin_authenticated"] = True
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(origin)

    def test_temporal_pending_round_trip_keeps_requested_cap_and_half_cap(self) -> None:
        pending, ledger, updated, plan, _, record = _pending_cursor()
        self.assertEqual(pending.mode, RETRY_PENDING)
        self.assertFalse(pending.public_fresh_permitted())
        assert pending.temporal is not None
        self.assertNotEqual(plan.requested_cap.hex(), plan.macro_width.hex())
        self.assertEqual(
            pending.temporal.predecessor_plan.requested_cap.hex(), plan.requested_cap.hex()
        )
        self.assertEqual(
            pending.temporal.predecessor_plan.macro_width.hex(), plan.macro_width.hex()
        )
        successor = temporal_successor_plan(plan)
        self.assertEqual(pending.temporal.successor_plan, successor)
        self.assertEqual(
            pending.temporal.imp1_half_cap.hex(), (plan.macro_width / 2.0).hex()
        )
        self.assertEqual(pending.ledger.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(pending.ledger.inherited_snapshot, ledger.inherited_snapshot)
        self.assertEqual(pending.ledger.accumulated_debit_vector, ledger.accumulated_debit_vector)
        self.assertEqual(ledger_sha256(pending.ledger), ledger_sha256(updated))
        mapping = encode_hlt17_cursor(pending)
        restored = restore_hlt17_cursor(mapping)
        self.assertEqual(encode_hlt17_cursor(restored), mapping)
        self.assertEqual(restored.temporal.rejection.as_mapping(), record)
        self.assertNotIn("pending_owner", mapping)

    def test_source_and_cfl_overlay_do_not_overwrite_temporal_identity(self) -> None:
        pending, ledger, _, plan, _, _ = _pending_cursor()
        successor = pending.next_attempt_plan()
        assert successor is not None
        source = close_source_cfl_overlay(
            pending,
            owner=SOURCE_OWNER,
            attempted_plan=successor,
            evidence=_source_evidence(successor),
        )
        self.assertEqual(source.mode, RETRY_PENDING)
        assert source.temporal is not None and source.overlay is not None
        self.assertEqual(source.temporal.successor_plan, pending.temporal.successor_plan)
        self.assertEqual(source.overlay.attempted_plan, successor)
        self.assertNotEqual(source.overlay.next_plan, source.temporal.successor_plan)
        self.assertEqual(source.ledger.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(source.ledger.inherited_snapshot, ledger.inherited_snapshot)
        self.assertEqual(source.source_current, 1)
        self.assertEqual(source.source_total, 1)
        self.assertEqual(source.cfl_current, 0)
        self.assertEqual(source.latest_owner, SOURCE_OWNER)

        cfl_plan = source.next_attempt_plan()
        assert cfl_plan is not None
        mixed = close_source_cfl_overlay(
            source,
            owner="cfl",
            attempted_plan=cfl_plan,
            evidence=_cfl_evidence(cfl_plan),
        )
        self.assertEqual(mixed.mode, RETRY_PENDING)
        assert mixed.temporal is not None and mixed.overlay is not None
        self.assertEqual(mixed.temporal.successor_plan, pending.temporal.successor_plan)
        self.assertEqual(mixed.source_current, 1)
        self.assertEqual(mixed.cfl_current, 1)
        self.assertEqual(mixed.source_total, 1)
        self.assertEqual(mixed.cfl_total, 1)
        self.assertEqual(mixed.latest_owner, "cfl")
        self.assertEqual(mixed.ledger.current_macro_step_temporal_retry_count, 1)
        restored = restore_hlt17_cursor(encode_hlt17_cursor(mixed))
        self.assertEqual(restored.temporal.successor_plan, pending.temporal.successor_plan)
        self.assertEqual(restored.overlay.owner, "cfl")
        self.assertNotIn("pending_owner", encode_hlt17_cursor(restored))

    def test_new_temporal_after_overlay_replaces_temporal_keeps_source_cfl_counts(
        self,
    ) -> None:
        pending, _, _, _, _, _ = _pending_cursor()
        successor = pending.next_attempt_plan()
        assert successor is not None
        overlaid = close_source_cfl_overlay(
            pending,
            owner=SOURCE_OWNER,
            attempted_plan=successor,
            evidence=_source_evidence(successor),
        )
        next_plan = overlaid.next_attempt_plan()
        assert next_plan is not None
        record = _rejection(overlaid.ledger, next_plan, preparation=_sha("prep-2"))
        updated = append_imp1_rejection(overlaid.ledger, record)
        replaced = close_temporal_rejection(
            overlaid,
            predecessor_plan=next_plan,
            predecessor_ledger=overlaid.ledger,
            updated_ledger=updated,
            rejection=record,
            preparation_sha256=record["preparation_sha256"],
            assessment_sha256=record["assessment_sha256"],
        )
        self.assertEqual(replaced.mode, RETRY_PENDING)
        self.assertIsNone(replaced.overlay)
        assert replaced.temporal is not None
        self.assertEqual(replaced.temporal.predecessor_plan, next_plan)
        self.assertNotEqual(
            replaced.temporal.successor_plan, pending.temporal.successor_plan
        )
        self.assertEqual(replaced.source_current, 1)
        self.assertEqual(replaced.cfl_current, 0)
        self.assertEqual(replaced.ledger.current_macro_step_temporal_retry_count, 2)

    def test_fine_adoption_resets_per_macro_counts_and_retains_history(self) -> None:
        pending, ledger, _, plan, identities, _ = _pending_cursor()
        successor = pending.next_attempt_plan()
        assert successor is not None
        overlaid = close_source_cfl_overlay(
            pending,
            owner=SOURCE_OWNER,
            attempted_plan=successor,
            evidence=_source_evidence(successor),
        )
        debit = overlaid.ledger.accumulated_debit_vector
        committed = accept_imp1_step(
            overlaid.ledger,
            final_time=plan.endpoint,
            state_sha256=_sha("accepted-state"),
            step_index=overlaid.ledger.current_step_index + 4,
            transaction_serial=overlaid.ledger.current_transaction_serial + 20,
            debit_vector=debit,
            assessment_sha256=_sha("commit-assessment"),
        )
        new_descriptor = _sha("new-accepted-descriptor")
        adopted = close_fine_adoption(
            overlaid,
            committed_ledger=committed,
            descriptor_sha256=new_descriptor,
            physical_state_sha256=committed.current_state_sha256,
            coordinates_sha256=identities["coordinates"],
            monitor_sha256=_sha("monitor-after"),
            causal_state_sha256=_sha("causal-after"),
            tracer_sha256=_sha("tracer-after"),
            guard_configuration_sha256=identities["guard"],
        )
        self.assertEqual(adopted.descriptor_sha256, new_descriptor)
        self.assertNotEqual(adopted.descriptor_sha256, overlaid.descriptor_sha256)
        self.assertEqual(adopted.mode, FRESH_READY)
        self.assertIsNone(adopted.temporal)
        self.assertIsNone(adopted.overlay)
        self.assertEqual(adopted.source_current, 0)
        self.assertEqual(adopted.cfl_current, 0)
        self.assertEqual(adopted.source_total, 1)
        self.assertEqual(adopted.cfl_total, 0)
        self.assertEqual(adopted.ledger.current_macro_step_temporal_retry_count, 0)
        self.assertEqual(adopted.ledger.cumulative_temporal_retry_count, 1)
        self.assertEqual(
            adopted.ledger.serialized_temporal_rejections,
            overlaid.ledger.serialized_temporal_rejections,
        )
        self.assertEqual(adopted.ledger.inherited_snapshot, ledger.inherited_snapshot)
        self.assertEqual(adopted.checkpoint_generation, 1)
        self.assertEqual(
            restore_hlt17_cursor(encode_hlt17_cursor(adopted)).kernel_transition_digest,
            adopted.kernel_transition_digest,
        )
        mapping = encode_hlt17_cursor(adopted)
        self.assertIn("kernel_transition_digest", mapping)
        self.assertNotIn("journal_tip_sha256", mapping)
        self.assertNotIn("journal_parent_sha256", mapping)

    def test_fine_adoption_refuses_to_retain_the_old_descriptor(self) -> None:
        pending, _, _, plan, identities, _ = _pending_cursor()
        committed = accept_imp1_step(
            pending.ledger,
            final_time=plan.endpoint,
            state_sha256=_sha("accepted-state"),
            step_index=pending.ledger.current_step_index + 4,
            transaction_serial=pending.ledger.current_transaction_serial + 20,
            debit_vector=pending.ledger.accumulated_debit_vector,
            assessment_sha256=_sha("commit-assessment"),
        )
        with self.assertRaises(HLT17IMP1Error):
            close_fine_adoption(
                pending,
                committed_ledger=committed,
                descriptor_sha256=pending.descriptor_sha256,
                physical_state_sha256=committed.current_state_sha256,
                coordinates_sha256=identities["coordinates"],
                monitor_sha256=_sha("monitor-after"),
                causal_state_sha256=_sha("causal-after"),
                tracer_sha256=_sha("tracer-after"),
                guard_configuration_sha256=identities["guard"],
            )

    def test_inherited_source_cfl_baselines_are_not_zeroed_at_seed(self) -> None:
        ledger = make_synthetic_fixture()[3]
        identities = _identities(ledger.current_state_sha256)
        cursor = seed_hlt17_cursor(
            member_key="RK4-9",
            event_target=SYNTHETIC_TARGET,
            ledger=ledger,
            descriptor_sha256=identities["descriptor"],
            physical_state_sha256=ledger.current_state_sha256,
            coordinates_sha256=identities["coordinates"],
            synthetic_callable_binding_sha256=identities["source"],
            monitor_sha256=identities["monitor"],
            causal_state_sha256=identities["causal"],
            tracer_sha256=identities["tracer"],
            guard_configuration_sha256=identities["guard"],
            source_total=4,
            cfl_total=2,
        )
        self.assertEqual(cursor.source_current, 0)
        self.assertEqual(cursor.cfl_current, 0)
        self.assertEqual(cursor.source_total, 4)
        self.assertEqual(cursor.cfl_total, 2)
        restored = restore_hlt17_cursor(encode_hlt17_cursor(cursor))
        self.assertEqual(restored.source_total, 4)
        self.assertEqual(restored.cfl_total, 2)

    def test_fresh_source_overlay_does_not_fabricate_temporal_evidence(self) -> None:
        cursor, _, _ = _seed_cursor()
        plan = _plan()
        overlaid = close_source_cfl_overlay(
            cursor,
            owner=SOURCE_OWNER,
            attempted_plan=plan,
            evidence=_source_evidence(plan),
        )
        self.assertEqual(overlaid.mode, FRESH_READY)
        self.assertIsNone(overlaid.temporal)
        self.assertIsNotNone(overlaid.overlay)
        self.assertFalse(overlaid.public_fresh_permitted())
        self.assertEqual(overlaid.ledger.current_macro_step_temporal_retry_count, 0)

    def test_source_lattice_exhaustion_belongs_to_source_owner(self) -> None:
        tiny = plan_imp1_lattice(
            SYNTHETIC_START, SYNTHETIC_TARGET, float(TDG6_MINIMUM_MACRO_STEP)
        )
        with self.assertRaises(HLT17OwnerRetryExhausted) as caught:
            strictly_smaller_tdg7_plan(tiny, owner=SOURCE_OWNER)
        self.assertEqual(caught.exception.owner, SOURCE_OWNER)
        cursor, _, _ = _seed_cursor()
        exhausted = close_source_cfl_overlay(
            cursor,
            owner=SOURCE_OWNER,
            attempted_plan=tiny,
            evidence=_source_evidence(tiny),
        )
        assert exhausted.overlay is not None
        self.assertTrue(exhausted.overlay.exhausted)
        self.assertIsNone(exhausted.overlay.next_plan)
        self.assertEqual(exhausted.ledger.current_macro_step_temporal_retry_count, 0)
        self.assertEqual(exhausted.source_current, 1)

    def test_altered_plan_closure_generation_and_ledger_fail_restore(self) -> None:
        pending, _, _, _, _, _ = _pending_cursor()
        mapping = encode_hlt17_cursor(pending)
        temporal = dict(mapping["temporal"])
        plan = dict(temporal["successor_plan"])
        plan["target_limited"] = 1
        temporal["successor_plan"] = plan
        mapping["temporal"] = temporal
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        mapping = encode_hlt17_cursor(pending)
        mapping["guard_configuration_sha256"] = _sha("tampered-guard")
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        mapping = encode_hlt17_cursor(pending)
        mapping["cursor_generation"] = pending.cursor_generation + 5
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        mapping = encode_hlt17_cursor(pending)
        ledger = dict(mapping["imp1_ledger"])
        ledger["current_macro_step_temporal_retry_count"] = 0
        mapping["imp1_ledger"] = ledger
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

    def test_cursor_is_immutable_and_does_not_import_sealed_owners(self) -> None:
        cursor, _, _ = _seed_cursor()
        with self.assertRaises(FrozenInstanceError):
            cursor.mode = RETRY_PENDING  # type: ignore[misc]
        source = MODULE.read_text(encoding="utf-8")
        self.assertNotIn("from .hlt16", source)
        self.assertNotIn("import hlt16", source)
        self.assertNotIn("from .proto15_runtime", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertIn("pending_owner", source)
        digest = cursor_sha256(cursor)
        self.assertEqual(len(digest), 64)
        self.assertEqual(HLT17_SOURCE_RETRY_CAP, 32)
        self.assertEqual(HLT17_CFL_RETRY_CAP, 32)
        self.assertIsInstance(cursor, HLT17IMP1Cursor)
        self.assertEqual(len(IMP1_COMPLETE_STATE_CHANNELS), 18)
        self.assertEqual(IMP1_COMPLETE_STATE_CHANNELS[0], "u:alpha")
        self.assertEqual(IMP1_COMPLETE_STATE_CHANNELS[-1], "q:chi")

    def test_legacy_schema_and_implementation_identity_swaps_are_refused(self) -> None:
        cursor, _, _ = _seed_cursor()
        mapping = encode_hlt17_cursor(cursor)
        mapping["schema"] = HLT17_LEGACY_IMP1_CURSOR_SCHEMA
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(cursor)
        mapping["implementation_id"] = "tdg11_imp1_exact_bernstein_with_dual_fallback_v1"
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)
        mapping = encode_hlt17_cursor(cursor)
        mapping["admission_runtime_id"] = "silent-backend-change"
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)


if __name__ == "__main__":
    unittest.main()
