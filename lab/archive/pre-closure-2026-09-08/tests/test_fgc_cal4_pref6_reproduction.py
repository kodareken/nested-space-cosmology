from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal4_pref6 as cal4  # noqa: E402


class FGCCAL4PREF6ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.loaded = cal4.load_config()
        cls.result = cal4.load_canonical_result()

    def test_stored_result_and_raw_campaign_validate_without_replay(self) -> None:
        cal4.verify_canonical(replay=False)
        payload = self.result["artifact_payload"]
        diagnosis = payload["common_event_ownership_diagnosis"]
        self.assertTrue(
            diagnosis["amplitude_5over2"]["both_methods_prospectively_admitted"]
        )
        self.assertFalse(
            diagnosis["amplitude_3"]["both_methods_prospectively_admitted"]
        )
        self.assertFalse(payload["epistemic_boundary"]["mechanism_question_answered"])

    def test_claim_and_diagnosis_mutations_fail_closed(self) -> None:
        mutations = (
            lambda value: value["gate_status"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
            lambda value: value["artifact_payload"][
                "common_event_ownership_diagnosis"
            ]["amplitude_5over2"].__setitem__(
                "both_methods_prospectively_admitted", False
            ),
            lambda value: value["artifact_payload"][
                "common_event_ownership_diagnosis"
            ]["amplitude_3"].__setitem__(
                "both_methods_prospectively_admitted", True
            ),
            lambda value: value["artifact_payload"][
                "common_event_ownership_diagnosis"
            ]["amplitude_5over2"]["methods"]["RK4"].__setitem__(
                "old_full_domain_reduction_failures_boundary_localized", False
            ),
            lambda value: value["artifact_payload"]["decision"].__setitem__(
                "new_protocol_required", "FGC-2-SF1-PROTO7"
            ),
        )
        for mutation in mutations:
            changed = deepcopy(self.result)
            mutation(changed)
            with self.subTest(mutation=mutation):
                with self.assertRaises(ValueError):
                    cal4._validate_stored(changed, self.loaded)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(
                json.dumps(self.result, sort_keys=False) + "\n", encoding="utf-8"
            )
            with self.assertRaises(ValueError):
                cal4.load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"artifact_id":"x","artifact_id":"y"}\n')
            with self.assertRaises(ValueError):
                cal4.load_canonical_result(duplicate)

    def test_config_ownership_change_is_rejected_before_interpretation(self) -> None:
        source = cal4.DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as temporary:
            changed = Path(temporary) / "changed.toml"
            changed.write_text(
                source.replace(
                    "finest_pair_individual_budgets_required = true",
                    "finest_pair_individual_budgets_required = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                cal4.load_config(changed)


if __name__ == "__main__":
    unittest.main()
