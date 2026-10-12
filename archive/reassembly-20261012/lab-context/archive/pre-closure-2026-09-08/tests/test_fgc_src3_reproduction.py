from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_src3 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCSRC3ReproductionTests(unittest.TestCase):
    def test_stored_canonical_result_reproduces(self) -> None:
        expected = _canonical(reproduce(DEFAULT_CONFIG)).encode("utf-8")
        self.assertEqual(DEFAULT_OUTPUT.read_bytes(), expected)

    def test_raw_threshold_weakening_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            'raw_complete_residual_maximum = "1e-12"',
            'raw_complete_residual_maximum = "2e-12"',
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "scope, evaluator, proof"):
                load_config(path)

    def test_claim_promotion_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "scope, evaluator, proof"):
                load_config(path)

    def test_clean_clone_path_reproduces_without_the_optional_raw_bundle(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        replacements = {
            "runs/fgc-2-sf1/src2/affine-wall-capture/manifest.json":
                "runs/fgc-2-sf1/src3-not-present/manifest.json",
            "runs/fgc-2-sf1/src2/affine-wall-capture/replay-result.json":
                "runs/fgc-2-sf1/src3-not-present/replay-result.json",
            "runs/fgc-2-sf1/src2/affine-wall-capture/affine-wall-fixture.npz":
                "runs/fgc-2-sf1/src3-not-present/affine-wall-fixture.npz",
            "runs/fgc-2-sf1/src2/affine-wall-capture/arithmetic-result.json":
                "runs/fgc-2-sf1/src3-not-present/arithmetic-result.json",
        }
        mutated = source
        for old, new in replacements.items():
            mutated = mutated.replace(old, new)
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory(dir=REPOSITORY) as temporary:
            path = Path(temporary) / "clean-clone.toml"
            path.write_text(mutated, encoding="utf-8")
            record = reproduce(path)
        evidence = record["artifact_payload"]["captured_full_grid_development_evidence"]
        self.assertTrue(evidence["raw_bundle_optional_for_clean_clone_reproduction"])
        self.assertTrue(evidence["stable_raw_gate_passed"])

    def test_partial_raw_bundle_presence_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            "runs/fgc-2-sf1/src2/affine-wall-capture/manifest.json",
            "runs/fgc-2-sf1/src3-not-present/manifest.json",
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory(dir=REPOSITORY) as temporary:
            path = Path(temporary) / "partial.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "only partially present"):
                reproduce(path)

    def test_rebound_duplicate_control_identifier_fails_closed(self) -> None:
        control_path = REPOSITORY / "configs/fgc/fgc-1-src3-controls.json"
        control = json.loads(control_path.read_text(encoding="utf-8"))
        control["independent_nontrivial_controls"][1]["control_id"] = "DYADIC-A"
        payload = _canonical(control).encode("utf-8")
        digest = hashlib.sha256(payload).hexdigest()
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory(dir=REPOSITORY) as temporary:
            temporary_path = Path(temporary)
            rebound = temporary_path / "controls.json"
            rebound.write_bytes(payload)
            relative = rebound.relative_to(REPOSITORY).as_posix()
            mutated = source.replace(
                'src3_control_fixture = "configs/fgc/fgc-1-src3-controls.json"',
                f'src3_control_fixture = "{relative}"',
            ).replace(
                'src3_control_fixture_sha256 = "2987770d68fd922aa56ef61b3a0a63064ea9e3461bf670eca765c217800c6e5e"',
                f'src3_control_fixture_sha256 = "{digest}"',
            )
            config_path = temporary_path / "mutated.toml"
            config_path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid or repeated"):
                reproduce(config_path)


if __name__ == "__main__":
    unittest.main()
