from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_modified_harmonic import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1ModifiedHarmonicReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-MHG1")
        self.assertEqual(config.tilde_normal_factor, 4)
        self.assertEqual(config.hat_normal_factor, 9)
        payload = record()
        committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(payload, committed)
        self.assertTrue(
            payload["modified_harmonic_certificate"][
                "all_declared_exact_checks_pass"
            ]
        )
        self.assertFalse(payload["gate_status"]["full_fgc1_hyp1_health_gate_passed"])
        self.assertFalse(payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_unknown_key_cone_order_and_path_traversal_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            unknown = root / "unknown.toml"
            unknown.write_text(source + "\nunknown_key = true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(unknown)

            reversed_cones = root / "reversed.toml"
            reversed_cones.write_text(
                source.replace(
                    "tilde_normal_factor = 4\nhat_normal_factor = 9",
                    "tilde_normal_factor = 9\nhat_normal_factor = 4",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "1 < tilde"):
                load_config(reversed_cones)

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


if __name__ == "__main__":
    unittest.main()
