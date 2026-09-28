from __future__ import annotations

from copy import deepcopy
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts import reproduce_fgc_pro18_pref26 as reproduction
from recursive_horizons.fgc.evolution.proto18_pref26_binder import Proto18Pref26BinderError


class Pref26ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = reproduction.reproduce()

    def test_canonical_result_round_trip_and_all_nonclaims(self) -> None:
        encoded = reproduction._canonical(self.result)
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "pref26.json"
            path.write_bytes(encoded)
            self.assertEqual(reproduction._load_result(path), self.result)
        self.assertTrue(self.result["gate_status"]["PREF26_completed"])
        self.assertTrue(self.result["gate_status"]["actual_raw_bundle_byte_and_hash_recomputed"])
        self.assertTrue(self.result["gate_status"]["actual_legacy_history_semantically_replayed"])
        self.assertTrue(self.result["gate_status"]["AUTH1_input_evidence_derived"])
        self.assertTrue(all(value is False for value in self.result["nonclaims"].values()))
        for forbidden in (
            "AUTH1_authorized", "AUTH1_committed", "HLT15_GEN1_authorized",
            "future_namespace_created", "candidate_execution_authorized",
            "FGCQR_holdout_execution_authorized", "DEF1_execution_authorized",
            "physical_transition_claim_authorized", "dark_sector_mechanism_derived",
        ):
            self.assertFalse(self.result["gate_status"][forbidden], forbidden)

    def test_stored_result_tamper_is_detectable_against_reproduction(self) -> None:
        tampered = deepcopy(self.result)
        tampered["gate_status"]["AUTH1_authorized"] = True
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "pref26.json"
            path.write_bytes(reproduction._canonical(tampered))
            stored = reproduction._load_result(path)
            self.assertNotEqual(stored, reproduction.reproduce())
        # Noncanonical byte tampering is separately rejected at load time.
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "pref26.json"
            path.write_bytes(reproduction._canonical(self.result) + b" ")
            with self.assertRaises(Proto18Pref26BinderError):
                reproduction._load_result(path)

    def test_sealed_lineage_contract_drift_fails_before_source_execution(self) -> None:
        config = reproduction._load_toml(reproduction.CONFIG)
        altered = deepcopy(config)
        altered["predecessor"]["sealed_commit"] = "0" * 40
        with self.assertRaises(Proto18Pref26BinderError) as caught:
            reproduction._verify_sealed_lineage(altered)
        self.assertEqual(caught.exception.stop_id, "PREF26_PREDECESSOR_LINEAGE_DRIFT")

    def test_cli_has_no_caller_selected_output_path(self) -> None:
        # PREF26 must never be redirectable into a future PROTO17 namespace.
        # Argparse rejects the option before reproduction or filesystem work.
        with redirect_stderr(StringIO()):
            with patch.object(
                sys,
                "argv",
                [
                    "reproduce_fgc_pro18_pref26.py",
                    "--write",
                    "--output",
                    "forbidden-result.json",
                ],
            ):
                with self.assertRaises(SystemExit) as caught:
                    reproduction.main()
        self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
