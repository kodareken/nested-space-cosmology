from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_hlt14_mon14 as reproduction


class HLT14ReproductionTests(unittest.TestCase):
    def test_executed_qualification_matches_frozen_contract(self) -> None:
        result = reproduction.record()
        summary = result["artifact_payload"]["synthetic_qualification_summary"]
        self.assertEqual(summary["case_count"], 72)
        self.assertEqual(
            summary["case_digest"],
            "445765193eddba247b7452635ca382388950a55a3a2bfb128e0b125f530b6baf",
        )
        self.assertEqual(summary["sealed_PREF25_hash_vector"], reproduction.SEALED_HASHES)
        self.assertTrue(summary["sealed_PREF25_hash_vector_reproduced"])
        self.assertTrue(summary["temporary_directories_only"])
        self.assertTrue(summary["in_process_exception_faults_not_power_loss_proof"])
        self.assertFalse(summary["real_production_instance_requirements_completed"])
        self.assertEqual(len(summary["case_names"]), len(set(summary["case_names"])))

    def test_every_forbidden_scope_authorization_and_claim_promotion_fails(self) -> None:
        with reproduction.CONFIG.open("rb") as handle:
            base = tomllib.load(handle)
        mutations = (
            ("scope", "real_production_state_or_raw_bundle_opened"),
            ("scope", "real_historical_journal_payload_semantics_replayed"),
            ("scope", "external_production_launch_authority_instance_validated"),
            ("scope", "production_namespace_materialized"),
            ("scope", "pretrajectory_operation"),
            ("scope", "trajectory_read"),
            ("scope", "mechanism_question_answered"),
            ("authorization", "pretrajectory_operation_authorized"),
            ("authorization", "production_namespace_creation_authorized"),
            ("authorization", "fresh_GR0_calibration_authorized"),
            ("authorization", "candidate_execution_authorized"),
            ("authorization", "physical_claim_authorized"),
            ("claims", "HLT14_real_production_state_loaded"),
            ("claims", "HLT14_external_production_launch_authority_validated"),
            ("claims", "HLT14_pretrajectory_operation_authorized"),
            ("claims", "SGBL_execution_authorized"),
            ("claims", "FGCQR_holdout_execution_authorized"),
            ("claims", "physical_transition_claim_authorized"),
        )
        for section, key in mutations:
            altered = deepcopy(base)
            altered[section][key] = True
            with self.subTest(section=section, key=key):
                with self.assertRaises(ValueError):
                    reproduction.validate_config_mapping(altered)

    def test_result_nonclaim_promotion_and_schema_extension_fail(self) -> None:
        expected = reproduction.record()
        promoted = deepcopy(expected)
        promoted["gate_status"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaises(ValueError):
            reproduction.validate_result_mapping(promoted)
        extended = deepcopy(expected)
        extended["unexpected_claim"] = True
        with self.assertRaises(ValueError):
            reproduction.validate_result_mapping(extended)

    def test_cli_writes_and_verifies_one_canonical_result(self) -> None:
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "hlt14.json"
            subprocess.run(
                [sys.executable, str(ROOT / "scripts/reproduce_fgc_hlt14_mon14.py"),
                 "--output", str(output)],
                cwd=ROOT, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            reproduction.verify(output)
            source = output.read_text(encoding="utf-8")
            observed = json.loads(source)
            self.assertEqual(source, reproduction._canonical(observed))
            output.write_text('{"artifact_id":"FGC-1-HLT14-MON14","artifact_id":"forged"}\n')
            with self.assertRaises(ValueError):
                reproduction.load_canonical_result(output)


if __name__ == "__main__":
    unittest.main()
