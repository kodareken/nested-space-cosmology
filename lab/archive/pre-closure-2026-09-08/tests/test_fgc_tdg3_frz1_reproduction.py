from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_tdg3_frz1 as freeze


class FGCTDG3FRZ1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = freeze.verify_canonical()

    def test_canonical_prospective_freeze_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], freeze.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "prospective_TDG3_native_grid_discriminator_freeze",
        )
        self.assertTrue(
            self.record["gate_status"][
                "TDG3_interpolation_ownership_discriminator_frozen"
            ]
        )
        self.assertFalse(
            self.record["gate_status"]["TDG3_actual_histories_diagnosed"]
        )

    def test_freeze_consumes_no_raw_history_and_authorizes_no_gate(self) -> None:
        payload = self.record["artifact_payload"]
        lineage = payload["immutable_lineage"]
        boundary = payload["claim_boundary"]
        self.assertFalse(lineage["source_checkpoint_history_arrays_loaded"])
        self.assertFalse(lineage["source_checkpoint_resumed"])
        self.assertFalse(lineage["state_advanced"])
        self.assertFalse(boundary["actual_terminal_histories_diagnosed"])
        self.assertFalse(boundary["common_grid_resampling_performed"])
        self.assertFalse(boundary["continuum_surrogate_enclosure_claimed"])
        self.assertFalse(boundary["replacement_temporal_admission_defined"])
        self.assertFalse(boundary["new_trajectory_authorized"])
        self.assertFalse(boundary["candidate_branch_opened"])

    def test_synthetic_preflight_and_mutation_controls_are_complete(self) -> None:
        payload = self.record["artifact_payload"]
        synthetic = payload["synthetic_preflight"]
        self.assertTrue(synthetic["all_control_classes_separated"])
        self.assertTrue(
            synthetic["common_grid_interpolation_owner_control"][
                "TDG2_identifies_interpolation_owner"
            ]
        )
        self.assertTrue(
            synthetic["common_grid_interpolation_owner_control"][
                "TDG3_native_estimators_agree_nonzero"
            ]
        )
        self.assertEqual(set(payload["mutation_controls"].values()), {True})

    def test_tracked_result_is_canonical(self) -> None:
        raw = freeze.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, freeze._canonical_bytes(stored))

    def test_config_mutations_fail_closed(self) -> None:
        config = freeze.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["source_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["frozen_discriminator"]["top_bins"] = [27, 28, 29, 30, 31]
        attacks.append(value)
        value = deepcopy(config)
        value["frozen_discriminator"]["common_grid_resampling_allowed"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"][
            "piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure"
        ] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["replacement_temporal_admission_defined"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                freeze.validate_config_data(attacked)

    def test_optional_checkpoint_is_hashed_but_never_loaded(self) -> None:
        expected = freeze._expected_lineage()["source_checkpoint_sha256"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.npz"
            self.assertFalse(path.exists())
            path.write_bytes(b"not-the-checkpoint")
            self.assertNotEqual(freeze._sha(path), expected)


if __name__ == "__main__":
    unittest.main()
