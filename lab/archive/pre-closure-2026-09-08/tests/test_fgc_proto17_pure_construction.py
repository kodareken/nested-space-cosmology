from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto17_pure_construction import (
    Proto17ConstructionError,
    Proto17RecoveryInvalidStop,
    build_genesis,
    checkpoint,
    construct_common_event,
    recover_lone_common_event,
    run_pure_construction_qualification,
    synthetic_common_event_predecessor,
    synthetic_common_event_fixture,
    synthetic_genesis_fixture,
    validate_common_event_successor,
)


class Proto17PureConstructionTests(unittest.TestCase):
    def test_complete_genesis_binds_distinct_physical_and_state_object_hashes(self):
        built = build_genesis(synthetic_genesis_fixture())
        checkpoint = built["checkpoint"]
        self.assertEqual(checkpoint["campaign_generation"], 0)
        self.assertEqual(checkpoint["journal_tip_sha256"], "0" * 64)
        self.assertEqual(len(checkpoint["states"]), 6)
        self.assertEqual(built["store_plan"]["checkpoint"], checkpoint)
        self.assertTrue(built["store_plan"]["checkpoint_path"].startswith("checkpoints/00000000000000000000-"))
        descriptor = built["genesis_spec"]["member_descriptors"][0]
        self.assertNotEqual(
            descriptor["state_object"]["physical_state_sha256"], descriptor["state_object_canonical_sha256"]
        )
        self.assertEqual(checkpoint["states"][descriptor["member_key"]], descriptor["state_object_canonical_sha256"])

    def test_common_event_is_acyclic_and_recoverable(self):
        common = synthetic_common_event_fixture()
        predecessor, prior_journal = common["predecessor"], common["prior_journal_records"]
        edge = construct_common_event(predecessor, prior_journal_records=prior_journal)
        receipt = edge["receipt"]
        self.assertNotIn("checkpoint_sha256", receipt["payload"])
        self.assertNotIn("successor_cursor_set_sha256", receipt["payload"])
        self.assertEqual(receipt["payload"]["sequence"], len(prior_journal) + 1)
        recovered = recover_lone_common_event(
            predecessor, [receipt], prior_journal_records=prior_journal
        )
        self.assertEqual(recovered, edge)
        self.assertEqual(
            edge["checkpoint"]["ledger_set_sha256"], predecessor["ledger_set_sha256"]
        )
        self.assertEqual(
            edge["checkpoint"]["accepted_state_set_sha256"], predecessor["accepted_state_set_sha256"]
        )
        self.assertEqual(
            validate_common_event_successor(
                predecessor, prior_journal_records=prior_journal, receipt=receipt,
                successor_checkpoint=edge["checkpoint"],
            ), edge,
        )

    def test_mutations_and_noncontiguous_evidence_fail_closed(self):
        base = synthetic_genesis_fixture()
        for mutate in (
            lambda x: x["member_descriptors"][0]["state_object"].__setitem__("physical_state_sha256", "0" * 64),
            lambda x: x["member_descriptors"][0].__setitem__("state_object_canonical_sha256", "0" * 64),
            lambda x: x.__setitem__("physical_input_manifest_sha256", "0" * 64),
            lambda x: x["member_descriptors"][0]["state_object"]["physical_arrays"][0].__setitem__("npz_storage_key", "other"),
            lambda x: x["member_descriptors"][0]["initial_TDG6_ledger"].__setitem__("current_macro_step_temporal_retry_count", 1),
            lambda x: x.__setitem__("common_event_index", 24),
            lambda x: x.__setitem__("restart_coordinate_time", {"rational": "3/2", "binary64_hex": float(1.5).hex()}),
            lambda x: x.__setitem__("protocol_artifact_id", "foreign"),
            lambda x: x["member_descriptors"][0]["state_object"]["physical_arrays"][1].__setitem__("npz_storage_key", "arr_u"),
            lambda x: x["physical_input_manifest"].__setitem__("relative_input_root", "../escape"),
            lambda x: x["physical_input_manifest"]["bundle_members"][0].__setitem__("relative_input_path", "/absolute.npz"),
        ):
            altered = deepcopy(base)
            mutate(altered)
            with self.assertRaises(Proto17ConstructionError):
                build_genesis(altered)
        common = synthetic_common_event_fixture()
        predecessor, prior_journal = common["predecessor"], common["prior_journal_records"]
        with self.assertRaises(Proto17ConstructionError):
            construct_common_event(predecessor, prior_journal_records=[{"record_kind": "X"}])
        edge = construct_common_event(predecessor, prior_journal_records=prior_journal)
        malformed = deepcopy(edge["receipt"])
        malformed["payload"]["completed_common_event_index"] += 1
        with self.assertRaises(Proto17RecoveryInvalidStop):
            recover_lone_common_event(predecessor, [malformed], prior_journal_records=prior_journal)
        with self.assertRaises(Proto17RecoveryInvalidStop):
            recover_lone_common_event(predecessor, [edge["receipt"], edge["receipt"]], prior_journal_records=prior_journal)

    def test_nonzero_generation_cannot_restart_at_root_journal(self):
        common = synthetic_common_event_fixture()
        predecessor = common["predecessor"]
        with self.assertRaises(Proto17ConstructionError):
            checkpoint(
                protocol_artifact_id=predecessor["protocol_artifact_id"],
                campaign_id=predecessor["campaign_id"], campaign_generation=1,
                committed_common_event_index=predecessor["committed_common_event_index"],
                active_event_target_time=predecessor["active_event_target_time"],
                cursors=predecessor["cursors"], ledgers=predecessor["ledgers"],
                states=predecessor["states"], journal_tip_sha256="0" * 64,
                attempt_high_water=predecessor["attempt_high_water"],
                cursor_generation_high_water=predecessor["cursor_generation_high_water"],
                parent_checkpoint_sha256=predecessor["parent_checkpoint_sha256"],
            )

    def test_deterministic_qualification(self):
        evidence = run_pure_construction_qualification()
        self.assertTrue(evidence["lone_suffix_reconstructs_byte_identically"])
        self.assertTrue(evidence["all_qualified_mutations_rejected"])
        self.assertEqual(evidence["qualified_mutation_count"], len(evidence["qualified_mutation_names"]))
        self.assertIn("nonzero_generation_root_journal", evidence["qualified_mutation_names"])
        self.assertTrue(evidence["side_effect_free_reference_oracle_only"])
