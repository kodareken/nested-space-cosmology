from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.reproduce_fgc_hyp1_reduction import (
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGC_HYP1ReductionReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config(DEFAULT_CONFIG)
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-RED1")
        self.assertEqual(
            config.configuration["formulation"]["gauge_id"],
            "generalized_radial_adm_dynamic_lambda_unfixed_areal_radius_v1",
        )
        self.assertEqual(
            [fixture["model_id"] for fixture in config.configuration["fixtures"]],
            ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"],
        )
        payload = record(DEFAULT_CONFIG)
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(
            payload["classification"],
            "exact_spherical_reduction_and_principal_part_preflight_not_hyperbolicity_or_evolution",
        )
        self.assertTrue(
            all(payload["reduction_certificate"]["verified_exact_checks"].values())
        )
        self.assertFalse(payload["gate_status"]["full_fgc1_health_gate_passed"])
        self.assertFalse(
            payload["gate_status"]["publication_local_defocusing_gate_passed"]
        )
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        self.assertEqual(rendered, DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_unknown_key_path_traversal_and_invalid_regular_state_fail_closed(
        self,
    ) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutations = (
            (source + "\nunknown = 1\n", "keys differ"),
            (
                source.replace(
                    'action_config = "configs/fgc/fgc-1-action-gate.toml"',
                    'action_config = "../fgc-1-action-gate.toml"',
                ),
                "inside the repository",
            ),
            (
                source.replace('areal_radius = "4"', 'areal_radius = "0"', 1),
                "areal radius must be positive",
            ),
            (
                source.replace('alpha = "1"', 'alpha = "1/0"', 1),
                "canonical exact rational string",
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
