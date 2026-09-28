from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt4_mon4 as inherited  # noqa: E402
from scripts.reproduce_fgc_rsp1_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _validate_config_mapping,
    load_config,
    reproduce,
    validate_run_plan,
)


class FGCRSP1FreezeReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        observed = inherited._load_json_bytes(
            DEFAULT_OUTPUT.read_bytes(), "RSP1 authorization"
        )
        cls.observed = observed
        cls.reproduced = reproduce(
            DEFAULT_CONFIG,
            historical_namespace_evidence=observed["artifact_payload"][
                "namespace_precondition"
            ],
        )

    def test_stored_canonical_authorization_reproduces(self) -> None:
        self.assertEqual(self.observed, inherited._serial(self.reproduced))

    def test_authorization_exposes_t0_failures_without_reading_endpoint(self) -> None:
        payload = self.reproduced["artifact_payload"]
        self.assertEqual(len(payload["frozen_run_inputs"]), 6)
        self.assertEqual(len(payload["accepted_state_source_prechecks"]), 6)
        self.assertEqual(len(payload["initial_premise_events"]), 2)
        for event in payload["initial_premise_events"]:
            self.assertEqual(event["target_direct_pair_passed"], [False, False])
            self.assertFalse(event["study_admission_passed"])
            self.assertFalse(event["t0_target_failure_is_the_evolved_endpoint_outcome"])
            self.assertFalse(event["trajectory_advanced"])
        self.assertTrue(self.reproduced["gate_status"]["RSP1_execution_authorized"])
        self.assertFalse(
            self.reproduced["gate_status"]["amplitude_three_spectral_veto_cleared"]
        )
        self.assertFalse(
            self.reproduced["gate_status"]["FGCQR_holdout_execution_authorized"]
        )

    def test_claim_promotion_fails_closed(self) -> None:
        attacked = deepcopy(load_config(DEFAULT_CONFIG))
        attacked["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "scope, lineage, or claims"):
            _validate_config_mapping(attacked)

    def test_ratio_ceiling_mutation_fails_closed(self) -> None:
        source = (REPOSITORY / "configs/fgc/fgc-1-rsp1-run1.toml").read_text(
            encoding="utf-8"
        )
        mutated = source.replace(
            'unchanged_strict_nested_ratio_ceiling = "1/4"',
            'unchanged_strict_nested_ratio_ceiling = "1/3"',
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "frozen question"):
                validate_run_plan(path)


if __name__ == "__main__":
    unittest.main()
