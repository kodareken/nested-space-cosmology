from copy import deepcopy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_pro17_pref25 as binder


class PRO17PREF25ReproductionTests(unittest.TestCase):
    def test_reproduces_canonical_binder(self) -> None:
        binder.verify()
        record = binder.record()
        self.assertTrue(record["gate_status"]["HLT14_implementation_and_synthetic_qualification_authorized"])
        self.assertFalse(record["gate_status"]["HLT14_runtime_implemented"])
        audit = record["artifact_payload"]["qualified_controls_audit"]
        self.assertEqual(audit["sealed_control_count"], 14)
        self.assertTrue(audit["independent_equivalent_control_corpus_passed"])
        self.assertTrue(all(audit["binder_specific_mutations"].values()))

    def test_deferred_requirements_and_promotion_fail_closed(self) -> None:
        sealed = binder._sealed_json(binder._git_show("results/fgc-1-pro17-frz1.json"))
        protocol = binder.tomllib.loads(binder._git_show("configs/fgc/fgc-2-sf1-protocol-v17.toml").decode("utf-8"))
        altered = deepcopy(sealed["artifact_payload"]["mutation_contract"]["future_HLT14_runtime_requirements"])
        altered["semantic_historical_journal_payload_replay"]["completed"] = True
        with self.assertRaises(ValueError):
            binder._future(altered)
        promoted = deepcopy(sealed["gate_status"])
        promoted["PROTO17_successor_runtime_implemented"] = True
        self.assertNotEqual(promoted, sealed["gate_status"])
        self.assertTrue(all(binder._binder_mutations(sealed, protocol).values()))

    def test_config_lineage_and_claim_mutations_fail_closed(self) -> None:
        config = binder._load_config()
        mutated = deepcopy(config)
        mutated["immutable_lineage"]["checkpoint_commit"] = "0" * 40
        with self.assertRaises(ValueError):
            binder._validate_config(mutated)
        promoted = deepcopy(config)
        promoted["claims"]["HLT14_runtime_implemented"] = True
        with self.assertRaises(ValueError):
            binder._validate_config(promoted)

    def test_forbidden_import_audit_fails_closed(self) -> None:
        source = (ROOT / "src/recursive_horizons/fgc/evolution/proto17_construction_theorem.py").read_bytes()
        self.assertTrue(binder._no_oracle_import(source))
        self.assertFalse(binder._no_oracle_import(source + b"\nfrom recursive_horizons.fgc.evolution import protocol_v17 as hidden\n"))
        self.assertFalse(binder._no_oracle_import(source + b"\nimport recursive_horizons.fgc.evolution.proto17_pure_construction as hidden\n"))
        self.assertFalse(binder._no_oracle_import(source + b"\nfrom . import protocol_v17 as hidden\n"))

    def test_duplicate_and_schema_controls_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            binder._sealed_json(b'{"a":1,"a":2}')
        protocol = binder.tomllib.loads(binder._git_show("configs/fgc/fgc-2-sf1-protocol-v17.toml").decode("utf-8"))
        protocol["restart_inputs"]["member_keys_in_canonical_order"] = list(reversed(binder.MEMBERS))
        with self.assertRaises(ValueError):
            binder._schema(protocol)


if __name__ == "__main__":
    unittest.main()
