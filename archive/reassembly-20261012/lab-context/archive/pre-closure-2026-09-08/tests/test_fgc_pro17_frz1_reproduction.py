from copy import deepcopy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_pro17_frz1 as freeze


class PRO17FreezeTests(unittest.TestCase):
    def test_reproduces_canonical_freeze(self) -> None:
        freeze.verify()
        result = freeze.record(require_absent_namespace=False)
        self.assertTrue(result["gate_status"]["PROTO17_frozen"])
        self.assertFalse(result["gate_status"]["PROTO17_successor_runtime_implemented"])
        self.assertFalse(result["gate_status"]["PROTO17_output_namespace_created"])
        qualification = result["artifact_payload"]["synthetic_construction_qualification"]
        self.assertTrue(qualification["lone_suffix_reconstructs_byte_identically"])
        self.assertTrue(qualification["all_qualified_mutations_rejected"])
        self.assertEqual(qualification["qualified_mutation_count"], 14)
        self.assertEqual(len(qualification["qualified_mutation_names"]), 14)

    def test_rejects_genesis_schema_mutation(self) -> None:
        config = freeze._load_toml(freeze.CONFIG)
        protocol = freeze._load_toml(freeze.PROTOCOL)
        mutated = deepcopy(protocol)
        mutated["trusted_genesis"]["genesis_spec_required_fields"].remove(
            "physical_input_manifest_sha256"
        )
        with self.assertRaises(ValueError):
            freeze.validate_contract(config, mutated)

    def test_rejects_false_successor_claim(self) -> None:
        config = freeze._load_toml(freeze.CONFIG)
        protocol = freeze._load_toml(freeze.PROTOCOL)
        mutated = deepcopy(config)
        mutated["claims"]["PROTO17_successor_runtime_implemented"] = True
        with self.assertRaises(ValueError):
            freeze.validate_contract(mutated, protocol)
