from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_dom1_qift1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_canonical_result,
    load_config,
    record,
)


class FGCHYP1DOM1QIFT1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = record()

    def test_frozen_config_canonical_result_and_readable_exact_bounds(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-DOM1-QIFT1")
        self.assertEqual(config.parameter_half_width, Q(1, 65536))
        self.assertEqual(config.acceleration_half_width, Q(1, 128))
        committed = load_canonical_result()
        self.assertEqual(self.payload, committed)

        certificate = self.payload["quantified_branch_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertEqual(certificate["formulation"]["parameter_dimension"], 30)
        self.assertEqual(certificate["formulation"]["acceleration_dimension"], 6)
        contraction = certificate["krawczyk_contraction_certificate"]
        self.assertLess(
            Q(contraction["contraction_infinity_norm_upper_bound"]), Q(1, 250)
        )
        self.assertGreater(
            Q(contraction["minimum_strict_interior_inclusion_margin"]), Q(3, 512)
        )
        for enclosure in contraction["krawczyk_image_box"]:
            self.assertLess(
                max(abs(Q(enclosure["lower"])), abs(Q(enclosure["upper"]))),
                Q(1, 512),
            )

        margins = certificate["regular_domain_ledger"]["strict_positive_margins"]
        readable_bounds = {
            "areal_radius_positive_lower_margin": Q(3),
            "base_metric_determinant_negative_lower_margin": Q(99, 100),
            "coordinate_time_inverse_metric_negative_lower_margin": Q(99, 100),
            "effective_planck_coefficient_positive_lower_margin": Q(3),
            "tilde_auxiliary_determinant_negative_lower_margin": Q(1, 65),
            "hat_auxiliary_determinant_negative_lower_margin": Q(1, 29),
        }
        for name, lower_bound in readable_bounds.items():
            self.assertGreater(Q(margins[name]), lower_bound)
        self.assertEqual(Q(margins["reference_annulus_coordinate_margin"]), Q(7, 2))
        self.assertFalse(self.payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_unknown_key_path_order_width_and_promoted_gate_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            unknown = root / "unknown.toml"
            unknown.write_text(source + "\nunknown_key = true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(unknown)

            traversal = root / "traversal.toml"
            traversal.write_text(
                source.replace(
                    'first_order_config = "configs/fgc/fgc-1-hyp1-fo1-rc1.toml"',
                    'first_order_config = "../fgc-1-hyp1-fo1-rc1.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|canonical"):
                load_config(traversal)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(
                source.replace('uniform_half_width = "1/65536"', 'uniform_half_width = "2/131072"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)

            reordered = root / "reordered.toml"
            reordered.write_text(
                source.replace(
                    'parameter_order = ["u.h_tt", "u.h_tr"',
                    'parameter_order = ["u.h_tr", "u.h_tt"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen ordered string list"):
                load_config(reordered)

            zero_width = root / "zero-width.toml"
            zero_width.write_text(
                source.replace('uniform_half_width = "1/65536"', 'uniform_half_width = "0"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "must be positive"):
                load_config(zero_width)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "uniform_open_domain_symmetrizer_proven = false",
                    "uniform_open_domain_symmetrizer_proven = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted)

            canonical_result = DEFAULT_OUTPUT.read_text(encoding="utf-8")
            duplicate_top = root / "duplicate-top.json"
            duplicate_top.write_text(
                canonical_result.replace(
                    "{\n  \"acceleration_box\"",
                    "{\n  \"artifact_id\": \"FGC-1-HYP1-DOM1-QIFT1\",\n  \"acceleration_box\"",
                    1,
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate_top)

            duplicate_nested = root / "duplicate-nested.json"
            duplicate_nested.write_text(
                canonical_result.replace(
                    '  "acceleration_box": {\n    "center_is_flat_root_zero_acceleration": true,',
                    '  "acceleration_box": {\n    "dimension": 6,\n    "center_is_flat_root_zero_acceleration": true,',
                    1,
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate_nested)

            noncanonical = root / "noncanonical-result.json"
            noncanonical.write_text(canonical_result.rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)


if __name__ == "__main__":
    unittest.main()
