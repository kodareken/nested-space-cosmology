from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    validate_sf1_protocol,
)
from scripts.reproduce_fgc_pro4_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCProtocolV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol_path = (
            REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v4.toml"
        )
        cls.protocol_source = cls.protocol_path.read_text(encoding="utf-8")
        cls.protocol = tomllib.loads(cls.protocol_source)

    def test_overlay_validates_and_remains_fail_closed(self) -> None:
        certificate = validate_sf1_protocol(self.protocol)
        self.assertEqual(certificate["artifact_id"], "FGC-2-SF1-PROTO4")
        self.assertEqual(certificate["protocol_version"], 4)
        self.assertTrue(certificate["outcome_neutral_contract_validated"])
        self.assertTrue(certificate["branch_specific_stop_partition_validated"])
        self.assertTrue(certificate["weighted_spectral_admission_validated"])
        self.assertTrue(certificate["common_event_constraint_admission_validated"])
        self.assertFalse(certificate["resolved_holdout_manifest_present"])
        self.assertTrue(all(value is False for value in certificate["claims"].values()))

    def test_overlay_mutations_fail_closed(self) -> None:
        mutated = deepcopy(self.protocol)
        mutated["numerical_admission"]["spectrum"][
            "maximum_top_band_field_power_fraction"
        ] = "1/1024"
        with self.assertRaisesRegex(ValueError, "spectral admission"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        moved = mutated["stop_applicability"]["universal_runtime_stop_ids"].pop()
        mutated["stop_applicability"]["candidate_branch_only_stop_ids"].append(moved)
        with self.assertRaisesRegex(ValueError, "duplicates|stop-applicability partition"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        mutated["stop_applicability"]["universal_runtime_stop_ids"].append(
            mutated["stop_applicability"]["universal_runtime_stop_ids"][0]
        )
        with self.assertRaisesRegex(ValueError, "duplicates|stop-applicability partition"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        universal = mutated["stop_applicability"]["universal_runtime_stop_ids"]
        candidate = mutated["stop_applicability"]["candidate_branch_only_stop_ids"]
        universal[0], candidate[0] = candidate[0], universal[0]
        with self.assertRaisesRegex(ValueError, "stop-applicability partition"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        mutated["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "claims must remain fail-closed"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        mutated["amendment"]["predecessor_protocol_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "immutable CAL0 provenance"):
            validate_sf1_protocol(mutated)

        mutated = deepcopy(self.protocol)
        mutated["replacement"]["provenance"]["holdout_output_root"] = (
            "runs/fgc-2-sf1/proto3/holdout"
        )
        with self.assertRaisesRegex(ValueError, "output namespaces"):
            validate_sf1_protocol(mutated)

    def test_freeze_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "immutable_PROTO3_must_validate_unchanged = true",
                    "immutable_PROTO3_must_validate_unchanged = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_holdout_execution_authorized = false",
                    "FGCQR_holdout_execution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims must remain fail-closed"):
                load_config(promoted)

    def test_canonical_freeze_reproduction(self) -> None:
        payload = load_canonical_result(DEFAULT_OUTPUT)
        self.assertEqual(payload, record(DEFAULT_CONFIG))
        gates = payload["gate_status"]
        self.assertTrue(gates["PROTO4_outcome_neutral_protocol_frozen"])
        self.assertTrue(gates["PROTO4_immutable_lineage_verified"])
        self.assertFalse(gates["PROTO4_resolved_holdout_manifest_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        boundary = payload["artifact_payload"]["epistemic_boundary"]
        self.assertTrue(boundary["GR0_exploratory_outcomes_were_seen_before_freeze"])
        self.assertFalse(boundary["GR0_exploratory_runs_admissible_as_calibration_evidence"])
        self.assertFalse(boundary["FGCQR_evolution_outcome_inspected"])
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
