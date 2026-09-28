from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v7 import (  # noqa: E402
    validate_sf1_protocol_v7,
)
from scripts.reproduce_fgc_pro7_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCProtocolV7Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol_path = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v7.toml"
        cls.protocol = tomllib.loads(cls.protocol_path.read_text(encoding="utf-8"))

    def test_overlay_validates_and_remains_fail_closed(self) -> None:
        certificate = validate_sf1_protocol_v7(self.protocol)
        self.assertEqual(certificate["artifact_id"], "FGC-2-SF1-PROTO7")
        self.assertEqual(certificate["protocol_version"], 7)
        self.assertTrue(certificate["outcome_neutral_contract_validated"])
        self.assertTrue(certificate["accepted_state_source_failure_terminal"])
        self.assertTrue(certificate["all_source_only_unaccepted_proposal_failures_retryable"])
        self.assertTrue(certificate["candidate_endpoint_in_unaccepted_proposal"])
        self.assertTrue(certificate["non_source_failure_vetoes_retry"])
        self.assertTrue(certificate["complete_terminal_observability_required"])
        self.assertTrue(certificate["raw_source_residual_tolerance_unchanged"])
        self.assertTrue(certificate["all_non_source_stops_unchanged"])
        self.assertEqual(certificate["retry_factor"], "1/2")
        self.assertEqual(certificate["maximum_retry_count"], 32)
        self.assertFalse(certificate["resolved_holdout_manifest_present"])
        self.assertTrue(all(value is False for value in certificate["claims"].values()))

    def test_overlay_mutations_fail_closed(self) -> None:
        mutations = (
            (("replacement", "transaction_stop_ownership", "accepted_state_source_gate_failure_is_terminal"), False, "transaction-stop-ownership replacement"),
            (("replacement", "transaction_stop_ownership", "candidate_endpoint_is_unaccepted_until_transaction_commit"), False, "transaction-stop-ownership replacement"),
            (("replacement", "transaction_stop_ownership", "any_non_source_failure_in_unaccepted_proposal_vetoes_retry"), False, "transaction-stop-ownership replacement"),
            (("replacement", "transaction_stop_ownership", "retry_factor"), "3/4", "transaction-stop-ownership replacement"),
            (("replacement", "provenance", "holdout_output_root"), "runs/fgc-2-sf1/proto6/holdout", "provenance replacement"),
            (("claims", "FGCQR_holdout_execution_authorized"), True, "claims"),
        )
        for path, value, message in mutations:
            with self.subTest(path=path):
                mutated = deepcopy(self.protocol)
                target = mutated
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_sf1_protocol_v7(mutated)

    def test_freeze_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "every_source_only_unaccepted_proposal_failure_may_retry = true",
                    "every_source_only_unaccepted_proposal_failure_may_retry = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof contract"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "PROTO7_fresh_GR0_dynamic_calibration_authorized = false",
                    "PROTO7_fresh_GR0_dynamic_calibration_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)

    def test_canonical_freeze_reproduction(self) -> None:
        payload = load_canonical_result(DEFAULT_OUTPUT)
        self.assertEqual(payload, record(DEFAULT_CONFIG))
        gates = payload["gate_status"]
        self.assertTrue(gates["PROTO7_outcome_neutral_protocol_frozen"])
        self.assertTrue(gates["PROTO7_immutable_lineage_verified"])
        self.assertTrue(gates["PROTO7_unaccepted_proposal_source_retry_contract_frozen"])
        self.assertTrue(gates["PROTO7_non_source_retry_veto_contract_frozen"])
        self.assertTrue(gates["PROTO7_complete_failure_observability_contract_frozen"])
        self.assertFalse(gates["PROTO7_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        boundary = payload["artifact_payload"]["epistemic_boundary"]
        self.assertTrue(boundary["PROTO6_GR0_calibration_trajectory_was_seen"])
        self.assertFalse(boundary["first_nonzero_common_event_completed"])
        self.assertFalse(boundary["PROTO7_trajectory_read"])
        self.assertFalse(boundary["FGCQR_outcome_read"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        payload = record(DEFAULT_CONFIG)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(_canonical(payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
