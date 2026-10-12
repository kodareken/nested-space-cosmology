from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_mhg_implicit import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1MHGImplicitReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-MHG3-IMP1")
        self.assertEqual(config.fixture["model_id"], "FGC-QR")
        payload = record()
        committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(payload, committed)
        certificate = payload["implicit_acceleration_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            certificate[
                "local_smooth_implicit_coordinate_time_acceleration_branch_proven"
            ]
        )
        self.assertTrue(
            payload["gate_status"][
                "flat_fgcqr_local_implicit_acceleration_branch_passed"
            ]
        )
        self.assertFalse(payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_unknown_key_traversal_mutations_and_promoted_gate_fail_closed(self) -> None:
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
                    'reference_config = "configs/fgc/fgc-1-hyp1-mhg-reference.toml"',
                    'reference_config = "../fgc-1-hyp1-mhg-reference.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|does not exist"):
                load_config(traversal)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(
                source.replace(
                    'radial_domain_minimum = "1/2"',
                    'radial_domain_minimum = "2/4"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)

            wrong_model = root / "wrong-model.toml"
            wrong_model.write_text(
                source.replace('model_id = "FGC-QR"', 'model_id = "GR-0"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "model_id is not frozen"):
                load_config(wrong_model)

            root_mutation = root / "root-mutation.toml"
            root_mutation.write_text(
                source.replace(
                    'alpha = { value = "1" }',
                    'alpha = { value = "2" }',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "flat vacuum root"):
                load_config(root_mutation)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "evolution_authorized = false",
                    "evolution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted)


if __name__ == "__main__":
    unittest.main()
