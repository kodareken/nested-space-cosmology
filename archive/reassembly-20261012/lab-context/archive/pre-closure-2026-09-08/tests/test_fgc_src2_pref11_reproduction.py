from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_src2_pref11 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCSRC2PREF11ReproductionTests(unittest.TestCase):
    def test_stored_canonical_result_reproduces(self) -> None:
        expected = _canonical(reproduce(DEFAULT_CONFIG)).encode("utf-8")
        self.assertEqual(DEFAULT_OUTPUT.read_bytes(), expected)

    def test_threshold_weakening_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            'exact_rounded_root_residual_maximum = "1e-26"',
            'exact_rounded_root_residual_maximum = "1e-12"',
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "scope, proof, successor, or claim"):
                load_config(path)

    def test_claim_promotion_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            "PROTO12_frozen = false",
            "PROTO12_frozen = true",
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "scope, proof, successor, or claim"):
                load_config(path)

    def test_compact_point_bit_mutation_cannot_escape_the_raw_binding(self) -> None:
        fixture_path = (
            REPOSITORY / "configs/fgc/fgc-1-src2-pref11-point0.json"
        )
        fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
        fixture["lower_jet"]["u"][1] = "0x1.5af0e4cec54d9p-57"
        payload = _canonical(fixture).encode("utf-8")
        digest = hashlib.sha256(payload).hexdigest()
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory(dir=REPOSITORY) as temporary:
            temporary_path = Path(temporary)
            compact_path = temporary_path / "point0.json"
            compact_path.write_bytes(payload)
            relative = compact_path.relative_to(REPOSITORY).as_posix()
            mutated = source.replace(
                'compact_point_fixture = "configs/fgc/fgc-1-src2-pref11-point0.json"',
                f'compact_point_fixture = "{relative}"',
            ).replace(
                'compact_point_fixture_sha256 = "d3eef30ce918cd03c44968656d9d2d01d4a7c41c91ece0681a88405d176f9354"',
                f'compact_point_fixture_sha256 = "{digest}"',
            )
            config_path = temporary_path / "mutated.toml"
            config_path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "compact point differs"):
                reproduce(config_path)

    def test_partial_raw_bundle_presence_fails_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            'cap1_manifest = "runs/fgc-2-sf1/src2/affine-wall-capture/manifest.json"',
            'cap1_manifest = "runs/fgc-2-sf1/src2/affine-wall-capture/not-present.json"',
        )
        self.assertNotEqual(source, mutated)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "only partially present"):
                reproduce(path)


if __name__ == "__main__":
    unittest.main()
