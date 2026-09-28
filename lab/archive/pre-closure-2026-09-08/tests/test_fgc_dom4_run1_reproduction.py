from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_dom4_run1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCDOM4RUN1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_nonzero_domain_branch_witnesses_and_nonclaims(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "nonzero_classical_spherical_run_envelope_passed"
            ]
        )
        payload = self.payload["artifact_payload"]
        for key in (
            "nonzero_parameter_volume",
            "FGCQR_target_branch_only",
            "all_protocol_cases_contained",
            "branch_continuity_quantified",
            "classical_scope_only",
        ):
            self.assertTrue(payload[key])
        quantitative = payload["quantitative_evidence"]
        self.assertEqual(
            quantitative["parameter_container"]["dimensionless_coordinate_volume"],
            "9/4",
        )
        continuation = quantitative["seed_branch_continuation"]
        self.assertEqual(
            [item["phi_seed_factor"] for item in continuation["records"]],
            ["0", "1", "2"],
        )
        self.assertTrue(continuation["strict_nonsingular_branch_margins"])
        stress = quantitative["stress_initial_slice"]
        self.assertTrue(stress["regular_center"])
        self.assertTrue(stress["finite_mass"])
        self.assertTrue(stress["no_initial_trapped_sphere"])
        self.assertTrue(stress["inside_protocol_compactness_window"])
        boundary = quantitative["epistemic_boundary"]
        self.assertTrue(boundary["local_open_state_neighborhoods_are_not_a_global_run_tube"])
        self.assertFalse(boundary["time_evolution_performed"])
        self.assertFalse(boundary["FGCQR_holdout_outcome_inspected"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            volume = root / "volume.toml"
            volume.write_text(
                source.replace(
                    'dimensionless_coordinate_volume = "9/4"',
                    'dimensionless_coordinate_volume = "3"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "parameter envelope"):
                load_config(volume)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "retained_EFT_evolution_authorized = false",
                    "retained_EFT_evolution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "fail-closed"):
                load_config(promoted)

            noncanonical = root / "result.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)


if __name__ == "__main__":
    unittest.main()
