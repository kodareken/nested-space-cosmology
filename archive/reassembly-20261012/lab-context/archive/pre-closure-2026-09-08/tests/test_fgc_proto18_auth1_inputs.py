from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto18_auth1_inputs import (
    Proto18Auth1InputError, build_calibration_genesis,
    source_identities_from_pref26_result, validate_pref26_auth1_evidence,
)


class Auth1InputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((ROOT / "results/fgc-1-pro18-pref26.json").read_text())
        cls.evidence = cls.record["artifact_payload"]["AUTH1_input_evidence"]
        cls.hashes = {name: "a" * 64 for name in (
            "protocol_config_sha256", "protocol_freeze_result_sha256", "hlt13_result_sha256",
            "run_plan_sha256", "runtime_module_sha256", "adapter_module_sha256", "runner_sha256",
        )}

    def test_compact_evidence_derives_calibration_only_genesis_without_paths(self):
        output = build_calibration_genesis(
            self.evidence, campaign_id="proto17-gr0-calibration",
            namespace="runs/fgc-2-sf1/proto17/calibration", source_hashes=self.hashes,
            numerical_environment={"python_implementation": "CPython", "test": True},
        )
        self.assertEqual(output.genesis_spec["namespace"], "runs/fgc-2-sf1/proto17/calibration")
        self.assertEqual(output.genesis_checkpoint["campaign_generation"], 0)
        self.assertEqual(output.genesis_spec["member_keys_in_canonical_order"], [
            "RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
        ])

    def test_compact_evidence_and_holdout_mutations_fail(self):
        bad = deepcopy(self.evidence)
        bad["member_descriptors"][0]["step_index"] += 1
        with self.assertRaises(Proto18Auth1InputError):
            validate_pref26_auth1_evidence(bad)
        with self.assertRaises(Proto18Auth1InputError):
            build_calibration_genesis(
                self.evidence, campaign_id="x", namespace="runs/fgc-2-sf1/proto17/holdout",
                source_hashes=self.hashes, numerical_environment={"x": True},
            )

    def test_pref26_result_supplies_exact_predecompression_source_identities(self):
        identities = source_identities_from_pref26_result(self.record, self.evidence)
        self.assertEqual(tuple(identities), ("PROTO12", "RSP2"))
        self.assertEqual(
            identities["PROTO12"]["checkpoint_sha256"],
            "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
        )
        self.assertEqual(identities["RSP2"]["event_log_byte_count"], 15060)

    def test_pref26_source_manifest_mutation_fails(self):
        bad = deepcopy(self.record)
        bad["artifact_payload"]["source_archive_manifest"]["sources"][0]["checkpoint_byte_count"] += 1
        with self.assertRaises(Proto18Auth1InputError):
            source_identities_from_pref26_result(bad, self.evidence)


if __name__ == "__main__":
    unittest.main()
