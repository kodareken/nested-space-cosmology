from __future__ import annotations

"""Read-only status tests using byte copies of the authenticated GEN0 store."""

import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution import hlt16_lifecycle as lifecycle
from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import MemberCampaignState
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore
from recursive_horizons.fgc.evolution.hlt16_member_codec import GEN0_PROTOCOL, RUNTIME_PROTOCOL, encode_member
from recursive_horizons.fgc.evolution.hlt16_state_store import authorize_gen0_import
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS
from tests.fgc_gen0_fixture import copy_sealed_gen0_store, reconstruct_sealed_gen0_members


def _load_status_module():
    spec = importlib.util.spec_from_file_location(
        "fgc_pro19_event1_status_test_module",
        ROOT / "scripts" / "status_fgc_pro19_event1.py",
    )
    if spec is None or spec.loader is None:
        raise AssertionError("status command cannot be imported")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


status = _load_status_module()


class Proto19Event1StatusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.reconstructed = reconstruct_sealed_gen0_members(ROOT)
        cls.cursors = cls.reconstructed.evidence.checkpoint["cursors"]
        cls.campaign_id = str(cls.cursors[MEMBER_KEYS[0]]["campaign_id"])
        cls.snapshots = {
            key: encode_member(
                cls.reconstructed.members[key], cls.cursors[key],
                protocol_artifact_id=RUNTIME_PROTOCOL,
                campaign_id=cls.campaign_id,
                generation=0,
                provenance_source_protocol=GEN0_PROTOCOL,
            )
            for key in MEMBER_KEYS
        }
        cls.authorities = {key: authorize_gen0_import(cls.snapshots[key]) for key in MEMBER_KEYS}

    def copied_legacy(self) -> tuple[TemporaryDirectory, Path]:
        temporary = TemporaryDirectory()
        copied = Path(temporary.name) / "calibration"
        copy_sealed_gen0_store(ROOT, copied)
        return temporary, copied

    def _factory(self, key: str, descriptor: str) -> MemberCampaignState:
        ledger = self.reconstructed.members[key].temporal_ledger
        cursor = lifecycle.bridge_gen0(self.cursors[key], ledger, descriptor)
        return MemberCampaignState(
            descriptor_sha256=descriptor,
            cursor=dict(cursor.payload),
            ledger=p15._ledger_mapping(ledger),
        )

    def _bridge(self, copied: Path) -> HLT16CampaignStore:
        store = HLT16CampaignStore(copied)
        store.adopt_generation_zero(
            plan_sha256="1" * 64,
            authorization_commit="2" * 40,
            campaign_id=self.campaign_id,
            snapshots=self.snapshots,
            authorities=self.authorities,
            member_factory=self._factory,
            target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        )
        return store

    def test_legacy_generation_zero_is_safe_to_start_but_not_restart(self) -> None:
        temporary, copied = self.copied_legacy()
        with temporary:
            result = status.collect_status(ROOT, store_root=copied, process_lines=())
        self.assertEqual(result["state"], "verified_generation_zero_safe_to_start")
        self.assertEqual(result["checkpoint"], {
            "generation": 0,
            "sha256": "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d",
        })
        self.assertTrue(result["safe_to_start"])
        self.assertFalse(result["safe_to_restart"])
        self.assertEqual(set(result["members"]), set(MEMBER_KEYS))
        self.assertEqual({entry["mode"] for entry in result["members"].values()}, {"FRESH_READY"})
        self.assertEqual(result["last_complete_journal"]["sequence"], -1)

    def test_conservative_runner_detection_blocks_legacy_start(self) -> None:
        temporary, copied = self.copied_legacy()
        with temporary:
            result = status.collect_status(
                ROOT,
                store_root=copied,
                process_lines=("999 /usr/bin/python scripts/run_fgc_pro19_event1.py",),
            )
        self.assertEqual(result["state"], "verified_generation_zero_not_startable")
        self.assertFalse(result["safe_to_start"])
        self.assertTrue(result["runner_process"]["active"])
        self.assertEqual(result["runner_process"]["matches"][0]["pid"], 999)

    def test_piped_sid3_bootstrap_marker_is_detected(self) -> None:
        observed = status._runner_status((
            "777 /usr/bin/python -I -B - --process-marker "
            "FGC-PRO19-SID3-EVENT1 --resume",
        ))
        self.assertTrue(observed["active"])
        self.assertEqual(observed["matches"][0]["pid"], 777)

    def test_bridged_checkpoint_reports_restartable_complete_six_member_view(self) -> None:
        temporary, copied = self.copied_legacy()
        with temporary:
            bridged = self._bridge(copied)
            result = status.collect_status(ROOT, store_root=copied, process_lines=())
            checkpoint = bridged.authenticated_snapshot().checkpoint
        self.assertEqual(result["state"], "clean_checkpoint")
        self.assertEqual(result["checkpoint"], {"generation": 1, "sha256": checkpoint.sha256})
        self.assertTrue(result["safe_to_restart"])
        self.assertFalse(result["safe_to_start"])
        self.assertFalse(result["terminal_lock"])
        self.assertFalse(result["active_write_lock"])
        self.assertEqual(result["active_target"]["rational"], "3/2")
        self.assertEqual(set(result["members"]), set(MEMBER_KEYS))

    def test_legacy_runtime_suffix_is_not_misreported_as_safe_start(self) -> None:
        temporary, copied = self.copied_legacy()
        with temporary:
            (copied / "payloads").mkdir()
            (copied / "payloads" / ("a" * 64 + ".npz")).write_bytes(b"foreign")
            with self.assertRaises(status.Proto19StatusError):
                status.collect_status(ROOT, store_root=copied, process_lines=())


if __name__ == "__main__":
    unittest.main()
