from copy import deepcopy
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto17_hlt14_runtime import (
    Proto17HLT14Error,
    Proto17HLT14RecoveryStop,
    Proto17HLT14Runtime,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import (
    build_genesis,
    canonical,
    construct_common_event,
    cursor,
    digest,
    checkpoint,
    synthetic_common_event_fixture,
    synthetic_genesis_fixture,
)


class Proto17HLT14RuntimeTests(unittest.TestCase):
    def setUp(self):
        self.spec = synthetic_genesis_fixture()
        self.runtime = Proto17HLT14Runtime(self.spec)

    def _root(self, temporary: str) -> Path:
        return Path(temporary) / "hlt14-store"

    def test_materialize_replay_common_event_and_reopen(self):
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            genesis = self.runtime.materialize_generation_zero(root)
            self.assertEqual(genesis, self.runtime.genesis)
            self.assertEqual(genesis["checkpoint_sha256"], "01e8214dcbbef025fe7883670043659d3ec66bb47e9b3e9dd4b8089b8ba5f8de")
            prefix = self.runtime.replay_six_semantic_accepted_records(root)
            fixture = synthetic_common_event_fixture()
            self.assertEqual(prefix, fixture["predecessor"])
            self.assertEqual(self.runtime._records(root), fixture["prior_journal_records"])
            successor = self.runtime.commit_common_event(root)
            self.assertEqual(successor, construct_common_event(prefix, prior_journal_records=fixture["prior_journal_records"])["checkpoint"])
            self.assertEqual(self.runtime._records(root)[-1]["record_sha256"], "17f925e8489341392ac5df4fba1cc47b66a27153e64eca4455c3f701f37cbc17")
            self.assertEqual(successor["checkpoint_sha256"], "c4c31f53a6a97632addafd54d6b3ec08bdcffe4aa4321269f8c8bca957e5380b")
            self.assertEqual(self.runtime.open(root), successor)
            self.assertEqual(len(list((root / "states").glob("*.json"))), 6)
            self.assertEqual(len(list((root / "journal").glob("*.journal"))), 7)

    def test_lone_common_receipt_recovers_and_mixed_suffix_rejects(self):
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            self.runtime.replay_six_semantic_accepted_records(root)

            def interrupt(phase, target):
                if target == "common_event_receipt" and phase == "after_dir_fsync":
                    raise RuntimeError("injected")

            with self.assertRaises(RuntimeError):
                self.runtime.commit_common_event(root, fault_hook=interrupt)
            recovered = self.runtime.recover(root)
            self.assertEqual(recovered["committed_common_event_index"], 24)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            self.runtime.replay_six_semantic_accepted_records(root)
            self.runtime.commit_common_event(root)
            journal = root / "journal"
            copied = sorted(journal.glob("*.journal"))[-1]
            extra = journal / "00000000000000000008-bad.journal"
            extra.write_bytes(copied.read_bytes())
            with self.assertRaises((Proto17HLT14Error, Proto17HLT14RecoveryStop)):
                self.runtime.open(root)

    def test_semantic_and_authority_store_drift_reject(self):
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            self.runtime.replay_six_semantic_accepted_records(root)
            target = sorted((root / "journal").glob("*.journal"))[0]
            payload = target.read_bytes().replace(b'"synthetic":true', b'"synthetic":false')
            target.write_bytes(payload)
            with self.assertRaises(Proto17HLT14Error):
                self.runtime.open(root)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            other_input = deepcopy(self.spec)
            other_input["campaign_id"] = "foreign-campaign"
            other = build_genesis(other_input, _validate_derived=False)["genesis_spec"]
            foreign_runtime = Proto17HLT14Runtime(other)
            foreign_root = Path(temporary) / "foreign-store"
            foreign_runtime.materialize_generation_zero(foreign_root)
            with self.assertRaises(Proto17HLT14Error):
                self.runtime.open(foreign_root)

    def test_production_paths_are_refused(self):
        with self.assertRaises(Proto17HLT14Error):
            self.runtime.materialize_generation_zero(
                ROOT / "runs/fgc-2-sf1/proto17/calibration"
            )

    def test_generation_zero_requires_an_absent_root_and_incomplete_prefix_is_invalid(self):
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            root.mkdir()
            with self.assertRaises(Proto17HLT14Error):
                self.runtime.materialize_generation_zero(root)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)

            def interrupt(phase, target):
                if target == "accepted_journal" and phase == "after_dir_fsync":
                    raise RuntimeError("injected")

            with self.assertRaises(RuntimeError):
                self.runtime.replay_six_semantic_accepted_records(root, fault_hook=interrupt)
            with self.assertRaises(Proto17HLT14RecoveryStop):
                self.runtime.open(root)
            with self.assertRaises(Proto17HLT14RecoveryStop):
                self.runtime.recover(root)

    def test_generation_zero_staging_faults_never_adopt_partial_namespace(self):
        cases = (
            ("state_object", "after_fsync"),
            ("genesis_checkpoint", "after_dir_fsync"),
            ("namespace", "before_namespace_rename"),
            ("namespace", "after_namespace_rename"),
        )
        for target, phase in cases:
            with self.subTest(target=target, phase=phase), TemporaryDirectory() as temporary:
                root = self._root(temporary)

                def interrupt(observed_phase, observed_target):
                    if observed_target == target and observed_phase == phase:
                        raise RuntimeError("injected")

                with self.assertRaises(RuntimeError):
                    self.runtime.materialize_generation_zero(root, fault_hook=interrupt)
                self.assertFalse(root.exists())
                self.assertEqual(list(root.parent.glob(f".{root.name}.hlt14-stage-*")), [])

    def test_fully_rehashed_semantic_forgeries_reject(self):
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            predecessor = self.runtime.replay_six_semantic_accepted_records(root)
            records = self.runtime._records(root)

            # Rehash a structurally contiguous six-record chain after changing
            # only the semantic executed-plan value.  The semantic replayer,
            # not a stale digest, must reject it.
            forged = [__import__('copy').deepcopy(record) for record in records]
            forged[0]["payload"]["executed_plan"] = {"synthetic": False}
            parent = "0" * 64
            journal = root / "journal"
            for number, record in enumerate(forged, 1):
                record["journal_parent_sha256"] = parent
                bare = {key: value for key, value in record.items() if key != "record_sha256"}
                record["record_sha256"] = digest(bare)
                old = records[number - 1]
                old_path = journal / f"{number:020d}-{old['record_sha256']}.journal"
                old_path.unlink()
                encoded = canonical(record)
                (journal / f"{number:020d}-{record['record_sha256']}.journal").write_bytes(f"{len(encoded):08x} ".encode() + encoded + b"\n")
                parent = record["record_sha256"]
            # Rehash the checkpoint's journal tip too, keeping it internally
            # canonical but semantically inconsistent with the frozen prefix.
            altered = __import__('copy').deepcopy(predecessor)
            altered["journal_tip_sha256"] = parent
            bare = {key: value for key, value in altered.items() if key != "checkpoint_sha256"}
            altered["checkpoint_sha256"] = digest(bare)
            old_checkpoint = root / "checkpoints" / f"00000000000000000001-{predecessor['checkpoint_sha256']}.json"
            old_checkpoint.unlink()
            (root / "checkpoints" / f"00000000000000000001-{altered['checkpoint_sha256']}.json").write_bytes(canonical(altered))
            with self.assertRaises(Proto17HLT14Error):
                self.runtime.open(root)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            self.runtime.materialize_generation_zero(root)
            predecessor = self.runtime.replay_six_semantic_accepted_records(root)
            altered = __import__('copy').deepcopy(predecessor)
            old = altered["cursors"]["RK4-2049"]
            old["attempt_serial"] = 2
            old["attempt_id"] = f"{altered['campaign_id']}:RK4-2049:23:2"
            old.pop("cursor_chain_sha256")
            altered["cursors"]["RK4-2049"] = cursor(old)
            altered["attempt_high_water"]["RK4-2049"] = 2
            altered["cursor_generation_high_water"]["RK4-2049"] = 1
            forged = checkpoint(protocol_artifact_id=altered["protocol_artifact_id"], campaign_id=altered["campaign_id"], campaign_generation=1, committed_common_event_index=23, active_event_target_time=altered["active_event_target_time"], cursors=altered["cursors"], ledgers=altered["ledgers"], states=altered["states"], journal_tip_sha256=altered["journal_tip_sha256"], attempt_high_water=altered["attempt_high_water"], cursor_generation_high_water=altered["cursor_generation_high_water"], parent_checkpoint_sha256=altered["parent_checkpoint_sha256"])
            old_path = root / "checkpoints" / f"00000000000000000001-{predecessor['checkpoint_sha256']}.json"
            old_path.unlink()
            (root / "checkpoints" / f"00000000000000000001-{forged['checkpoint_sha256']}.json").write_bytes(canonical(forged))
            with self.assertRaises(Proto17HLT14Error):
                self.runtime.open(root)

    def test_all_common_event_fault_windows_recover_or_preserve_prior(self):
        phases = (
            "before_write", "after_write", "before_flush", "after_flush",
            "before_fsync", "after_fsync", "before_replace", "after_replace",
            "before_dir_fsync", "after_dir_fsync",
        )
        for target in ("common_event_receipt", "checkpoint"):
            for phase in phases:
                with self.subTest(target=target, phase=phase), TemporaryDirectory() as temporary:
                    root = self._root(temporary)
                    self.runtime.materialize_generation_zero(root)
                    self.runtime.replay_six_semantic_accepted_records(root)

                    def interrupt(observed_phase, observed_target):
                        if observed_target == target and observed_phase == phase:
                            raise RuntimeError("injected")

                    with self.assertRaises(RuntimeError):
                        self.runtime.commit_common_event(root, fault_hook=interrupt)
                    recovered = self.runtime.recover(root)
                    self.assertIn(recovered["committed_common_event_index"], {23, 24})
                    # A visible receipt always recovers to event 24; absent
                    # pre-replace receipt leaves the exact event-23 checkpoint.
                    self.assertEqual(self.runtime.open(root), recovered)


if __name__ == "__main__":
    unittest.main()
