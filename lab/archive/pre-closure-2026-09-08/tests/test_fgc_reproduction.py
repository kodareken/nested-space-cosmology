from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.reproduce_fgc_action import (
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCReproductionTests(unittest.TestCase):
    def test_frozen_configuration_and_record(self) -> None:
        config = load_config(DEFAULT_CONFIG)
        self.assertEqual(config.artifact_id, "FGC-1-ACT1")
        self.assertEqual([model.model_id.value for model in config.models], [
            "GR-0",
            "SGB-L",
            "FGC-QR",
        ])
        payload = record(DEFAULT_CONFIG)
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["project_version"], "0.11.0")
        self.assertEqual(
            payload["classification"],
            "action_and_background_algebra_not_principal_symbol_or_collapse_result",
        )
        self.assertTrue(all(payload["verified_algebraic_checks"].values()))
        self.assertTrue(
            payload["gate_status"][
                "declared_action_and_scalar_algebra_certificate_passed"
            ]
        )
        self.assertFalse(payload["gate_status"]["full_fgc1_health_gate_passed"])
        self.assertFalse(
            payload["gate_status"]["publication_local_defocusing_gate_passed"]
        )
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        self.assertEqual(rendered, DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_configuration_fails_closed_on_unknown_key(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.toml"
            path.write_text(source + "\nunknown = 1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(path)


if __name__ == "__main__":
    unittest.main()
