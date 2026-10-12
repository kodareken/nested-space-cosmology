from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt1_mon1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCHLT1MON1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_every_stop_and_transaction_contract_pass(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "classical_health_and_typed_stop_monitoring_verified"
            ]
        )
        payload = self.payload["artifact_payload"]
        for key in (
            "canonical_dimensionless_norms_executable",
            "all_protocol_thresholds_implemented_strictly",
            "every_accepted_stage_monitored",
            "nonmonotone_internal_stage_abscissae_supported",
            "strict_global_stage_transaction_order_enforced",
            "stop_occurs_before_threshold_crossing",
            "first_failed_premise_is_preserved",
            "scale_stops_not_promoted_to_EFT_validity",
        ):
            self.assertTrue(payload[key])
        injected = payload["quantitative_evidence"]["typed_stop_injections"]
        self.assertEqual(len(injected["records"]), 26)
        self.assertTrue(injected["every_reason_independently_stopped"])
        self.assertEqual(
            injected["transaction_control"]["simultaneous_first_reason"],
            "nonpositive_lapse",
        )
        self.assertTrue(
            injected["transaction_control"][
                "first_failure_immutable_on_later_attempt"
            ]
        )
        self.assertTrue(
            injected["transaction_control"][
                "SSPRK3_nonmonotone_stage_abscissae_supported"
            ]
        )
        self.assertTrue(
            injected["transaction_control"][
                "strict_global_stage_transaction_order_enforced"
            ]
        )
        boundary = payload["quantitative_evidence"]["epistemic_boundary"]
        self.assertTrue(boundary["scale_thresholds_are_classical_stop_proxies_not_EFT_remainder_control"])
        self.assertFalse(boundary["time_evolution_performed"])
        self.assertFalse(boundary["FGCQR_holdout_outcome_inspected"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    'normalized_constraint_infinity_max = "1/10000000"',
                    'normalized_constraint_infinity_max = "1/1000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "thresholds differ"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "retained_EFT_remainder_control_supplied = false",
                    "retained_EFT_remainder_control_supplied = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "scope differs"):
                load_config(promoted)

            extra = root / "extra.toml"
            extra.write_text(
                source.replace(
                    "further_acceptance_after_stop_forbidden = true\n\n[proof_contract]",
                    "further_acceptance_after_stop_forbidden = true\nunexpected_transaction_flag = true\n\n[proof_contract]",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "transaction keys differ"):
                load_config(extra)

            noncanonical = root / "result.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)


if __name__ == "__main__":
    unittest.main()
