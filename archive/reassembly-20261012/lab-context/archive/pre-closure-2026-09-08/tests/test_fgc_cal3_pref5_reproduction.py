from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal3_pref5 as cal3  # noqa: E402


class FGCCAL3PREF5ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.loaded = cal3.load_config()
        cls.result = cal3.load_canonical_result()

    def test_stored_result_and_raw_campaign_validate_without_replay(self) -> None:
        cal3.verify_canonical(replay=False)
        payload = self.result["artifact_payload"]
        replay = payload["deterministic_terminal_replay"]
        self.assertTrue(replay["last_accepted_state_passes_raw_gate"])
        self.assertEqual(
            replay["failed_terminal_stages"],
            ["rk4_k4", "candidate_endpoint"],
        )
        self.assertTrue(replay["one_further_halving_control"]["accepted"])
        self.assertFalse(payload["epistemic_boundary"]["mechanism_question_answered"])

    def test_claim_promotion_and_transaction_mutations_fail_closed(self) -> None:
        for mutation in (
            lambda value: value["gate_status"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
            lambda value: value["artifact_payload"][
                "deterministic_terminal_replay"
            ].__setitem__("last_accepted_state_passes_raw_gate", False),
            lambda value: value["artifact_payload"][
                "deterministic_terminal_replay"
            ].__setitem__("failed_terminal_stages", ["candidate_endpoint"]),
            lambda value: value["artifact_payload"][
                "deterministic_terminal_replay"
            ]["one_further_halving_control"].__setitem__("accepted", False),
            lambda value: value["artifact_payload"]["decision"].__setitem__(
                "new_protocol_required", "FGC-2-SF1-PROTO6"
            ),
        ):
            changed = deepcopy(self.result)
            mutation(changed)
            with self.assertRaises(ValueError):
                cal3._validate_stored(changed, self.loaded)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(
                json.dumps(self.result, sort_keys=False) + "\n", encoding="utf-8"
            )
            with self.assertRaises(ValueError):
                cal3.load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"artifact_id":"x","artifact_id":"y"}\n')
            with self.assertRaises(ValueError):
                cal3.load_canonical_result(duplicate)

    def test_config_change_is_rejected_before_interpretation(self) -> None:
        source = cal3.DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as temporary:
            changed = Path(temporary) / "changed.toml"
            changed.write_text(
                source.replace(
                    "every_source_only_failure_in_an_unaccepted_proposal_is_retryable = true",
                    "every_source_only_failure_in_an_unaccepted_proposal_is_retryable = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                cal3.load_config(changed)


if __name__ == "__main__":
    unittest.main()
