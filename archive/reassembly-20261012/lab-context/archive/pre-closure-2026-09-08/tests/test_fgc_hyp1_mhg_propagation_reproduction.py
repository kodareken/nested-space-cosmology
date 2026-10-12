from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_mhg_propagation import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1MHGPropagationReproductionTests(unittest.TestCase):
    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-MHG4-PROP1")
        self.assertEqual(
            config.activated_fixture["fixture_id"], "FGCQR_activated_generic"
        )
        payload = record()
        committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(payload, committed)
        certificate = payload["gauge_propagation_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            certificate["conditional_lower_order_reference_gauge_operator_derived"]
        )
        self.assertEqual(
            certificate["formulation"]["output_order"],
            ["nu=t", "nu=r", "nu=theta", "nu=phi"],
        )
        self.assertTrue(
            payload["gate_status"]
            ["conditional_lower_order_reference_gauge_operator_passed"]
        )
        self.assertFalse(payload["gate_status"]["complete_gauge_propagation_passed"])
        self.assertFalse(payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_unknown_traversal_sign_order_and_promoted_claims_fail_closed(self) -> None:
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
                    'implicit_config = "configs/fgc/fgc-1-hyp1-mhg-implicit.toml"',
                    'implicit_config = "../fgc-1-hyp1-mhg-implicit.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|does not exist"):
                load_config(traversal)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(
                source.replace('value = "2/7"', 'value = "4/14"', 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)

            wrong_sign = root / "wrong-sign.toml"
            wrong_sign.write_text(
                source.replace(
                    "E_phi*nabla^nu_phi+E_chi*nabla^nu_chi=0",
                    "E_phi*nabla^nu_phi-E_chi*nabla^nu_chi=0",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "noether_identity is not frozen"):
                load_config(wrong_sign)

            truncated_output = root / "truncated-output.toml"
            truncated_output.write_text(
                source.replace(
                    'output_order = ["nu=t", "nu=r", "nu=theta", "nu=phi"]',
                    'output_order = ["nu=t", "nu=r"]',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen ordered string list"):
                load_config(truncated_output)

            promoted_direct = root / "promoted-direct.toml"
            promoted_direct.write_text(
                source.replace(
                    "direct_full_metric_residual_divergence_evaluated = false",
                    "direct_full_metric_residual_divergence_evaluated = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "cannot promote"):
                load_config(promoted_direct)

            promoted_gate = root / "promoted-gate.toml"
            promoted_gate.write_text(
                source.replace(
                    "evolution_authorized = false", "evolution_authorized = true"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted_gate)


if __name__ == "__main__":
    unittest.main()
