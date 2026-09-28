from __future__ import annotations

import copy
from pathlib import Path
import tempfile
import unittest
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY))

from scripts import reproduce_fgc_pro12_frz1 as reproduction  # noqa: E402
from scripts import reproduce_fgc_hlt10_mon10 as successor  # noqa: E402


class FGCPro12Frz1ReproductionTests(unittest.TestCase):
    def test_canonical_freeze_is_proved_from_immutable_successor(self) -> None:
        # PRO12's absent-namespace observation is temporal pre-launch evidence.
        # Once the authorized campaign has populated that namespace, HLT10 is
        # the post-launch-safe verifier: it reads the exact freeze from its
        # immutable checkpoint and hash-matches the current artifact and every
        # inherited implementation premise.
        lineage = successor._immutable_lineage(successor.load_config())
        result = reproduction.load_canonical_result()
        self.assertEqual(lineage["freeze"]["artifact_id"], "FGC-1-PRO12-FRZ1")
        self.assertEqual(
            lineage["freeze"]["result_sha256"],
            successor.EXPECTED_LINEAGE["protocol_freeze_result_sha256"],
        )
        self.assertEqual(
            reproduction._sha(reproduction.DEFAULT_OUTPUT),
            lineage["freeze"]["result_sha256"],
        )
        gates = result["gate_status"]
        conditioning = result["artifact_payload"][
            "RSP1_pairwise_conditioning_recomputation"
        ]
        self.assertTrue(gates["PROTO12_frozen"])
        self.assertTrue(gates["PROTO12_pairwise_spectral_classifier_derived"])
        self.assertTrue(gates["RSP1_pairwise_conditioning_recomputed"])
        self.assertFalse(
            gates["RSP1_generic_all_field_raw_spectral_admission_passed"]
        )
        self.assertTrue(conditioning["all_pairwise_conditioning_admissions_passed"])
        self.assertFalse(conditioning["is_calibration_or_mechanism_evidence"])
        for method in ("RK4", "SSPRK3"):
            record = conditioning["methods"][method]
            self.assertFalse(record["raw_generic_all_field_admission_passed"])
            self.assertTrue(record["pairwise_admission"]["admission_passed"])
            self.assertTrue(record["pairwise_admission"]["saturation_used"])
            self.assertFalse(record["round_trip_is_a_continuum_error_bound"])

    def test_raw_ratios_and_all_map_witnesses_remain_serialized(self) -> None:
        result = reproduction.load_canonical_result()
        methods = result["artifact_payload"][
            "RSP1_pairwise_conditioning_recomputation"
        ]["methods"]
        for record in methods.values():
            raw = record["raw_spatial_spectral_admission"]
            guarded = record["pairwise_admission"]
            self.assertEqual(
                raw["field_power_tail_ratios"],
                guarded["field_power_tail_ratios"],
            )
            self.assertEqual(
                raw["derivative_power_tail_ratios"],
                guarded["derivative_power_tail_ratios"],
            )
            self.assertEqual(len(record["spectral_budgets"]), 3)
            self.assertEqual(len(record["spectral_tail_sensitivities"]), 3)
            for grid in record["spectral_tail_sensitivities"]:
                for witness in grid.values():
                    self.assertIn("top_band_erasure_perturbation_infinity", witness)
                    self.assertIn("round_trip_interpolation_infinity", witness)
                    self.assertIn("erasure_over_round_trip", witness)

    def test_removing_a_pairwise_guard_fails_closed(self) -> None:
        text = reproduction.DEFAULT_CONFIG.read_text()
        attacked = text.replace(
            "pairwise_saturation_requires_both_pair_budgets_and_map_witnesses = true",
            "pairwise_saturation_requires_both_pair_budgets_and_map_witnesses = false",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            path.write_text(attacked)
            with self.assertRaises(ValueError):
                reproduction.load_config(path)

    def test_checkpoint_hash_mismatch_fails_before_reclassification(self) -> None:
        config = reproduction.load_config()
        original_sha = reproduction._sha
        checkpoint = (
            reproduction.REPOSITORY
            / config["immutable_lineage"]["resolution_checkpoint"]
        ).resolve()

        def attacked_sha(path: Path) -> str:
            if path.resolve() == checkpoint:
                return "0" * 64
            return original_sha(path)

        with mock.patch.object(reproduction, "_sha", side_effect=attacked_sha):
            with self.assertRaises(ValueError):
                reproduction._recompute_rsp1_conditioning(config)

    def test_result_promotion_breaks_immutable_successor_binding(self) -> None:
        attacked = copy.deepcopy(reproduction.load_canonical_result())
        attacked["gate_status"]["PROTO12_fresh_GR0_dynamic_calibration_authorized"] = (
            True
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.json"
            path.write_text(reproduction._canonical(attacked))
            self.assertNotEqual(
                reproduction._sha(path),
                successor.EXPECTED_LINEAGE["protocol_freeze_result_sha256"],
            )


if __name__ == "__main__":
    unittest.main()
