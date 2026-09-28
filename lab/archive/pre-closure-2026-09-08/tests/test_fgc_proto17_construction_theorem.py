import ast
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto17_construction_theorem import (
    Proto17TheoremError, construct_successor, generation_zero_fixture,
    parse_frozen_contract, qualification, recover_only_lone_successor,
    sealed_synthetic_genesis, synchronized_predecessor,
)


PROTOCOL = (ROOT / "configs/fgc/fgc-2-sf1-protocol-v17.toml").read_bytes()
FREEZE = (ROOT / "configs/fgc/fgc-1-pro17-frz1.toml").read_bytes()


class Proto17ConstructionTheoremTests(unittest.TestCase):
    def test_parses_live_frozen_contract_without_oracle_import(self):
        contract = parse_frozen_contract(PROTOCOL, FREEZE)
        self.assertTrue(contract["historical_record_payload_semantics_deferred"])
        self.assertTrue(contract["raw_bytes_deferred"])

    def test_ast_proves_no_proto17_oracle_import_and_schema_bytes_are_live(self):
        source = ROOT / "src/recursive_horizons/fgc/evolution/proto17_construction_theorem.py"
        tree = ast.parse(source.read_text(encoding="utf-8"))
        modules = [node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertNotIn("recursive_horizons.fgc.evolution.protocol_v17", modules)
        self.assertNotIn("recursive_horizons.fgc.evolution.proto17_pure_construction", modules)
        with self.assertRaises(Proto17TheoremError):
            parse_frozen_contract(PROTOCOL.replace(b'"genesis_algorithm_id"', b'"bad_genesis_field"', 1), FREEZE)
        with self.assertRaises(Proto17TheoremError):
            parse_frozen_contract(PROTOCOL.replace(b'"source_bundle_key"', b'"bad_descriptor_field"', 1), FREEZE)
        with self.assertRaises(Proto17TheoremError):
            parse_frozen_contract(PROTOCOL.replace(b'"physical_state_sha256"', b'"bad_state_hash"', 1), FREEZE)
        with self.assertRaises(Proto17TheoremError):
            parse_frozen_contract(PROTOCOL.replace(b'"little_endian_c_bytes_sha256"', b'"bad_array_hash"', 1), FREEZE)

    def test_independent_generation_zero_and_successor(self):
        genesis = generation_zero_fixture()
        self.assertEqual(genesis["campaign_generation"], 0)
        self.assertEqual(genesis["journal_tip_sha256"], "0" * 64)
        predecessor, history = synchronized_predecessor()
        edge = construct_successor(predecessor, history)
        self.assertNotIn("checkpoint_sha256", edge["receipt"]["payload"])
        self.assertEqual(edge["checkpoint"]["ledger_set_sha256"], predecessor["ledger_set_sha256"])
        self.assertEqual(edge["checkpoint"]["accepted_state_set_sha256"], predecessor["accepted_state_set_sha256"])
        self.assertEqual(recover_only_lone_successor(predecessor, history, [edge["receipt"]]), edge)
        sealed = sealed_synthetic_genesis()
        self.assertEqual(len(sealed["genesis_spec"]["member_descriptors"]), 6)
        self.assertEqual(len(sealed["store_plan"]["state_objects"]), 6)
        self.assertEqual(sealed["store_plan"]["journal_records"], [])

    def test_generation_zero_helper_uses_frozen_state_schema_and_event_shapes(self):
        sealed = sealed_synthetic_genesis()
        self.assertEqual(generation_zero_fixture(), sealed["checkpoint"])
        for descriptor in sealed["genesis_spec"]["member_descriptors"]:
            state = descriptor["state_object"]
            self.assertEqual(state["schema_id"], "FGC-2-SF1-PROTO17-state-object-v1")
            shapes = {item["logical_name"]: item["shape"] for item in state["physical_arrays"]}
            self.assertEqual(shapes["event_proper_times"], [24, 48])
            self.assertEqual(shapes["event_fields"], [24, 48, 6])

    def test_root_reset_and_recovery_tails_fail_closed(self):
        predecessor, history = synchronized_predecessor()
        forged = deepcopy(predecessor); forged["journal_tip_sha256"] = "0" * 64
        with self.assertRaises(Proto17TheoremError):
            construct_successor(forged, [])
        edge = construct_successor(predecessor, history)
        with self.assertRaises(Proto17TheoremError):
            recover_only_lone_successor(predecessor, history, [])
        with self.assertRaises(Proto17TheoremError):
            recover_only_lone_successor(predecessor, history, [edge["receipt"], edge["receipt"]])

    def test_independent_mutation_corpus(self):
        record = qualification(PROTOCOL, FREEZE)
        self.assertTrue(record["lone_suffix_reconstructs_exactly"])
        self.assertTrue(record["all_mutations_rejected"])
        self.assertTrue(record["sealed_hashes_reproduced"])
        self.assertEqual(record["qualified_mutation_count"], 14)
        self.assertEqual(len(record["mutations_rejected"]), 14)
        self.assertTrue(record["structural_only"])
        self.assertFalse(record["runtime_or_physical_claim"])


if __name__ == "__main__":
    unittest.main()
