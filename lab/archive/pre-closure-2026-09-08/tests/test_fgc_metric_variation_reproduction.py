from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.reproduce_fgc_metric_variation import (
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCMetricVariationReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config(DEFAULT_CONFIG)
        self.assertEqual(config.artifact_id, "FGC-1-VAR1")
        self.assertEqual(config.dimension, 4)
        self.assertEqual(config.algorithm, "lcg64_v1")
        self.assertEqual(config.fixture_count, 8)
        payload = record(DEFAULT_CONFIG)
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(
            payload["classification"],
            "covariant_metric_equation_derivation_and_crosschecks_not_principal_symbol_or_collapse",
        )
        self.assertTrue(
            payload["gate_status"]
            ["var1_covariant_metric_equation_derived_and_cross_checked"]
        )
        self.assertFalse(payload["gate_status"]["full_fgc1_health_gate_passed"])
        self.assertFalse(
            payload["gate_status"]["publication_local_defocusing_gate_passed"]
        )
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        self.assertEqual(rendered, DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_unknown_key_wrong_dimension_and_path_traversal_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutations = (
            (source + "\nunknown = 1\n", "keys differ"),
            (source.replace("dimension = 4", "dimension = 5"), "must equal 4"),
            (
                source.replace(
                    'action_config = "configs/fgc/fgc-1-action-gate.toml"',
                    'action_config = "../fgc-1-action-gate.toml"',
                ),
                "inside the repository",
            ),
        )
        for mutated, expected in mutations:
            with self.subTest(expected=expected):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "bad.toml"
                    path.write_text(mutated, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, expected):
                        load_config(path)


if __name__ == "__main__":
    unittest.main()
