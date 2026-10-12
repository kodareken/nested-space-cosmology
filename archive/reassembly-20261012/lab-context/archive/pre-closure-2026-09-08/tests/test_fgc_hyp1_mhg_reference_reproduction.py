from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_mhg_reference import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1MHGReferenceReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-MHG2-REF1")
        self.assertEqual(config.radial_domain_minimum.numerator, 1)
        self.assertEqual(config.radial_domain_minimum.denominator, 2)
        payload = record()
        committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(payload, committed)
        certificate = payload["reference_gauge_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            payload["gate_status"][
                "complete_reference_gauge_residual_passed"
            ]
        )
        self.assertFalse(
            payload["gate_status"]["implicit_second_time_derivative_branch_passed"]
        )
        self.assertFalse(payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_unknown_key_traversal_promoted_gate_and_center_fail_closed(self) -> None:
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
                    'action_config = "configs/fgc/fgc-1-action-gate.toml"',
                    'action_config = "../fgc-1-action-gate.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|does not exist"):
                load_config(traversal)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace("evolution_authorized = false", "evolution_authorized = true"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted)

            center = root / "center.toml"
            center.write_text(
                source.replace("center_included = false", "center_included = true"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "must not claim"):
                load_config(center)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(
                source.replace('radial_domain_minimum = "1/2"', 'radial_domain_minimum = "2/4"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)


if __name__ == "__main__":
    unittest.main()
