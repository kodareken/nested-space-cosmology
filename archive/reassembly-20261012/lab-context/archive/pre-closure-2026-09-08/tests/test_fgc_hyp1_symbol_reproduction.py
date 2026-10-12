from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.reproduce_fgc_hyp1_symbol import (
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1SymbolReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config(DEFAULT_CONFIG)
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-SYM1")
        self.assertEqual(config.configuration["symbol"]["covector"], "xi_A=(-c,1)")
        self.assertEqual(
            config.configuration["symbol"]["alternate_metric_patch"],
            "metric_tt_row_covers_c_equals_minus_shift",
        )
        payload = record(DEFAULT_CONFIG)
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(
            payload["classification"],
            "exact_pointwise_covariant_scalar_quotient_and_adm_formulation_preflight_not_full_hyperbolicity_or_evolution",
        )
        self.assertTrue(
            payload["symbol_certificate"]["all_declared_exact_checks_pass"]
        )
        self.assertTrue(
            all(
                payload["symbol_certificate"]["verified_exact_checks"].values()
            )
        )
        self.assertFalse(payload["gate_status"]["full_fgc1_hyp1_health_gate_passed"])
        self.assertFalse(payload["gate_status"]["evolution_authorized"])
        self.assertFalse(
            payload["gate_status"]["publication_local_defocusing_gate_passed"]
        )
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        self.assertEqual(rendered, DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_unknown_key_path_traversal_and_patch_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutations = (
            (source + "\nunknown = 1\n", "keys differ"),
            (
                source.replace(
                    'reduction_config = "configs/fgc/fgc-1-hyp1-reduction.toml"',
                    'reduction_config = "../fgc-1-hyp1-reduction.toml"',
                ),
                "inside the repository",
            ),
            (
                source.replace(
                    'alternate_metric_patch = "metric_tt_row_covers_c_equals_minus_shift"',
                    'alternate_metric_patch = "none"',
                ),
                "alternate_metric_patch is not frozen",
            ),
            (
                source.replace(
                    "require_naive_gr_comparator_rejected = true",
                    "require_naive_gr_comparator_rejected = false",
                ),
                "must be true",
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
