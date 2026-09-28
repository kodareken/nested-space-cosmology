"""Focused tests for the REC1 exact-recovery event wrapper."""
from __future__ import annotations

from pathlib import Path
import shutil
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from scripts import run_fgc_tdg8_rcv3_rec1_event as runner
from recursive_horizons.fgc.evolution import tdg8_rcv3_rec1_authority as authority
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore


ROOT = Path(__file__).resolve().parents[1]
LIVE_STORE = ROOT / authority.DESTINATION_PATH
AUTHORITY_COMMIT = "a" * 40


def _receipt() -> authority.TDG8RCV3REC1Authority:
    return authority.TDG8RCV3REC1Authority(
        authority_commit=AUTHORITY_COMMIT,
        prior_auth1_commit=authority.PRIOR_AUTH1_COMMIT,
        normalization_fix_commit=authority.NORMALIZATION_FIX_COMMIT,
        projection_id=authority.PROJECTION_ID,
        original_execution_commit=authority.ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=authority.ORIGINAL_PLAN_SHA256,
        campaign_id=authority.CAMPAIGN_ID,
        destination_path=authority.DESTINATION_PATH,
        entry_checkpoint_sha256=authority.ENTRY_CHECKPOINT_SHA256,
        predicted_checkpoint_sha256=authority.PREDICTED_CHECKPOINT_SHA256,
        config_sha256="b" * 64, result_sha256="c" * 64,
        implementation_inventory=(), environment={},
    )


def _plan() -> SimpleNamespace:
    return SimpleNamespace(
        sha256=authority.ORIGINAL_PLAN_SHA256,
        authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
        protocol_artifact_id=authority.TARGET_PROTOCOL,
        campaign_id=authority.CAMPAIGN_ID,
        branch="GR-0", amplitude="3", common_event_index=authority.EVENT,
        start_time=SimpleNamespace(rational="23/16"),
        target_time=SimpleNamespace(rational="3/2"),
        members=tuple(SimpleNamespace(member_key=key) for key in runner.MEMBER_KEYS),
    )


def _boundary(phase: authority.REC1Phase, **changes) -> authority.REC1Boundary:
    generation = 9 if phase in {
        authority.REC1Phase.EXACT_PRE_RECOVERY,
        authority.REC1Phase.GEN9_RECONCILE_PENDING,
    } else 10
    values = {
        "phase": phase, "checkpoint_generation": generation,
        "checkpoint_sha256": (
            authority.ENTRY_CHECKPOINT_SHA256 if generation == 9
            else authority.PREDICTED_CHECKPOINT_SHA256
        ),
        "journal_sequence": 10 if generation == 9 else 12,
        "journal_tip_sha256": authority.ENTRY_JOURNAL_SHA256 if generation == 9 else authority.SEQ12_SHA256,
        "event": authority.EVENT, "target": dict(authority.TARGET),
        "disposition": "nonterminal",
        "suffix_kinds": ("tdg6_rejection", "cursor_transition") if generation == 9 else (),
        "active_write": phase in {authority.REC1Phase.EXACT_PRE_RECOVERY,
                                  authority.REC1Phase.GEN9_RECONCILE_PENDING,
                                  authority.REC1Phase.GEN10_HANDOFF_PENDING},
        "writer_state": "stale_verified_lock" if phase is not authority.REC1Phase.EXACT_CLEAN_GEN10 else None,
        "writer_checkpoint_sha256": authority.ENTRY_CHECKPOINT_SHA256,
        "retired_writer_count": 0 if phase is authority.REC1Phase.EXACT_PRE_RECOVERY else 2,
        "exact_generation10_ancestor": generation >= 10,
        "physical_state_unchanged_at_recovery": True,
        "safe_to_continue": phase in {authority.REC1Phase.EXACT_PRE_RECOVERY,
                                      authority.REC1Phase.EXACT_CLEAN_GEN10},
        "terminal_lock_present": False,
    }
    values.update(changes)
    return authority.REC1Boundary(**values)


def _completed() -> dict[str, object]:
    return {
        "state": "event_one_complete",
        "authorization_commit": authority.ORIGINAL_EXECUTION_COMMIT,
        "campaign_id": authority.CAMPAIGN_ID,
        "plan_sha256": authority.ORIGINAL_PLAN_SHA256,
        "event": authority.SUCCESSOR_EVENT,
        "target": dict(authority.SUCCESSOR_TARGET),
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


class TDG8RCV3REC1RunnerTests(unittest.TestCase):
    def test_status_is_read_only(self) -> None:
        boundary = _boundary(authority.REC1Phase.EXACT_PRE_RECOVERY)
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_boundary", return_value=boundary),
            patch.object(runner, "HLT16CampaignStore") as store,
            patch.object(runner.event1, "_advance_authenticated_event") as engine,
        ):
            result = runner.status(ROOT, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(result["state"], authority.REC1Phase.EXACT_PRE_RECOVERY.value)
        store.assert_not_called()
        engine.assert_not_called()

    def test_real_copied_store_reconciles_to_clean_exact_gen10_without_PDE(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "calibration"
            shutil.copytree(LIVE_STORE, destination)
            store = HLT16CampaignStore(destination)
            plan = runner._original_plan(ROOT)
            before = authority.inspect_store(destination, _receipt())
            with patch.object(runner.event1, "_advance_authenticated_event") as engine:
                clean = runner._reconcile_to_clean_generation10(
                    plan, store, destination, _receipt(), before,
                )
            engine.assert_not_called()
            self.assertEqual(clean.phase, authority.REC1Phase.EXACT_CLEAN_GEN10)
            self.assertEqual(clean.checkpoint_sha256, authority.PREDICTED_CHECKPOINT_SHA256)
            self.assertFalse(clean.active_write)

    def test_engine_is_unreachable_when_exact_barrier_fails(self) -> None:
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_boundary", return_value=_boundary(authority.REC1Phase.EXACT_PRE_RECOVERY)),
            patch.object(runner, "_destination", return_value=ROOT / "copy"),
            patch.object(runner, "HLT16CampaignStore", return_value=SimpleNamespace()),
            patch.object(
                runner, "_reconcile_to_clean_generation10",
                side_effect=runner.TDG8RCV3REC1RunnerError("recovery", "barrier", "no gen10"),
            ),
            patch.object(runner.event1, "_advance_authenticated_event") as engine,
        ):
            with self.assertRaises(runner.TDG8RCV3REC1RunnerError):
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        engine.assert_not_called()

    def test_run_orders_clean_barrier_before_existing_engine(self) -> None:
        order: list[str] = []
        clean = _boundary(authority.REC1Phase.EXACT_CLEAN_GEN10)
        checkpoint = SimpleNamespace(
            generation=10, event=authority.EVENT,
            authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
            plan_sha256=authority.ORIGINAL_PLAN_SHA256,
        )
        templates = {key: object() for key in runner.MEMBER_KEYS}

        def reconcile(*_args, **_kwargs):
            order.append("barrier")
            return clean

        def advance(*_args, **_kwargs):
            order.append("engine")
            return _completed()

        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_boundary", return_value=_boundary(authority.REC1Phase.EXACT_PRE_RECOVERY)),
            patch.object(runner, "_destination", return_value=ROOT / "copy"),
            patch.object(runner, "HLT16CampaignStore", return_value=SimpleNamespace()),
            patch.object(runner, "_reconcile_to_clean_generation10", side_effect=reconcile),
            patch.object(runner, "_checkpoint", return_value=checkpoint),
            patch.object(runner, "build_static_gr0_shells", return_value=templates),
            patch.object(runner.event1, "_advance_authenticated_event", side_effect=advance),
        ):
            result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(order, ["barrier", "engine"])
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["candidate_branch_opened"])

    def test_event24_and_terminal_reentry_never_open_engine(self) -> None:
        completed = _boundary(
            authority.REC1Phase.EVENT24_COMPLETE,
            checkpoint_generation=20, event=24, target=dict(authority.SUCCESSOR_TARGET),
            disposition="event_complete", exact_generation10_ancestor=True,
            active_write=False, writer_state=None, writer_checkpoint_sha256=None,
        )
        for boundary in (completed, _boundary(
            authority.REC1Phase.TYPED_TERMINAL,
            disposition="invalid_terminal", terminal_lock_present=True,
            active_write=False, writer_state=None,
        )):
            with self.subTest(phase=boundary.phase):
                with (
                    patch.object(runner, "_fixed_root", return_value=ROOT),
                    patch.object(runner, "_execution_authority", return_value=_receipt()),
                    patch.object(runner, "_original_plan", return_value=_plan()),
                    patch.object(runner, "_boundary", return_value=boundary),
                    patch.object(runner.event1, "_advance_authenticated_event") as engine,
                ):
                    result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
                engine.assert_not_called()
                self.assertFalse(result["output_created"])

    def test_promoted_or_extra_event_result_is_rejected(self) -> None:
        for change in (
            {"physical_result_earned": True},
            {"event": 25},
            {"candidate_branch_opened": True},
        ):
            result = _completed()
            result.update(change)
            with self.subTest(change=change):
                with self.assertRaises(runner.TDG8RCV3REC1RunnerError):
                    runner._validate_event_result(result)


if __name__ == "__main__":
    unittest.main()
