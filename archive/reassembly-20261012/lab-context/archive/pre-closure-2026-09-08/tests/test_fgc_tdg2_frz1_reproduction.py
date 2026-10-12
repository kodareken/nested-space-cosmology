from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_tdg2_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _verify_optional_raw_hash,
    load_config,
    record,
    validate_config_data,
    verify_canonical,
)


class TDG2FreezeReproductionTests(unittest.TestCase):
    def test_canonical_freeze_reproduces(self) -> None:
        observed = verify_canonical()
        stored = json.loads(DEFAULT_OUTPUT.read_text())
        self.assertEqual(observed, stored)
        self.assertEqual(
            observed["classification"],
            "prospective_TDG2_absolute_tail_discriminator_freeze",
        )

    def test_synthetic_result_is_complete_but_pre_history(self) -> None:
        observed = record()
        payload = observed["artifact_payload"]
        self.assertTrue(
            payload["synthetic_preflight"]["all_control_classes_separated"]
        )
        self.assertFalse(
            payload["immutable_lineage"]["source_checkpoint_history_arrays_loaded"]
        )
        self.assertFalse(
            payload["claim_boundary"]["actual_terminal_histories_diagnosed"]
        )
        self.assertFalse(payload["claim_boundary"]["new_trajectory_authorized"])
        self.assertFalse(observed["gate_status"]["GR0_case_eligible"])
        self.assertTrue(
            observed["gate_status"]["TDG2_absolute_tail_discriminator_frozen"]
        )

    def test_threshold_lineage_and_claim_mutations_fail_closed(self) -> None:
        raw = load_config()
        attacks = []
        value = deepcopy(raw)
        value["frozen_discriminator"]["minimum_zero_power_order"] = "2"
        attacks.append(value)
        value = deepcopy(raw)
        value["frozen_discriminator"]["primary_point_counts"] = [1025, 2049, 4097]
        attacks.append(value)
        value = deepcopy(raw)
        value["immutable_lineage"]["source_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(raw)
        value["synthetic_controls"]["control_names"].pop()
        attacks.append(value)
        value = deepcopy(raw)
        value["claims"]["replacement_temporal_admission_defined"] = True
        attacks.append(value)
        value = deepcopy(raw)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                validate_config_data(attacked)

    def test_config_and_result_are_canonical_inputs(self) -> None:
        self.assertTrue(DEFAULT_CONFIG.is_file())
        raw = DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        expected = (
            json.dumps(
                stored,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode()
        self.assertEqual(raw, expected)

    def test_ignored_checkpoint_is_optional_but_fail_closed_when_present(self) -> None:
        payload = b"immutable-tdg2-checkpoint-control"
        expected = sha256(payload).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "checkpoint.npz"
            _verify_optional_raw_hash(checkpoint, expected)
            checkpoint.write_bytes(payload)
            _verify_optional_raw_hash(checkpoint, expected)
            checkpoint.write_bytes(payload + b"-mutated")
            with self.assertRaisesRegex(ValueError, "checkpoint hash differs"):
                _verify_optional_raw_hash(checkpoint, expected)


if __name__ == "__main__":
    unittest.main()
