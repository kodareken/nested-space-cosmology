from __future__ import annotations

import ast
from inspect import signature
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.pro20_source_factory import (  # noqa: E402
    CAMPAIGN_ID,
    EVENT_TARGET,
    PROTOCOL,
    Pro20RuntimeOrigin,
    Pro20SourceFactoryError,
    build_pro20_runtime_origin,
)


class Pro20SourceFactoryContractTests(unittest.TestCase):
    def test_public_builder_has_no_execution_or_synthetic_success_switch(self) -> None:
        self.assertEqual(
            tuple(signature(build_pro20_runtime_origin).parameters),
            ("origin", "static_input_bytes", "source_closure_sha256"),
        )
        self.assertEqual(PROTOCOL, "FGC-2-SF1-PROTO19")
        self.assertEqual(CAMPAIGN_ID, "FGC-2-SF1-PRO20-EVENT1")
        self.assertEqual(EVENT_TARGET.hex(), (3.0 / 2.0).hex())

    def test_module_has_no_old_runner_store_or_endpoint_path(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_source_factory.py"
        tree = ast.parse(path.read_text("utf-8"))
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        joined = " ".join(imported)
        self.assertNotIn("run_fgc_gr0_calibration", joined)
        self.assertNotIn("progression_inputs", joined)
        self.assertNotIn("campaign_store", joined)
        source = path.read_text("utf-8")
        self.assertNotIn("write_bytes", source)
        self.assertNotIn("publish_", source)

    def test_invalid_origin_refuses_before_static_physical_construction(self) -> None:
        with patch(
            "recursive_horizons.fgc.evolution.pro20_source_factory.build_static_gr0_shells",
            side_effect=AssertionError("physical construction called"),
        ), self.assertRaises((TypeError, ValueError)):
            build_pro20_runtime_origin(
                object(),  # type: ignore[arg-type]
                static_input_bytes={},
                source_closure_sha256="a" * 64,
            )

    def test_runtime_origin_cannot_be_empty_or_take_a_fake_digest(self) -> None:
        for digest in ("a", "A" * 64):
            with self.subTest(digest=digest), self.assertRaises(Pro20SourceFactoryError):
                Pro20RuntimeOrigin(digest, "b" * 64, {})
        with self.assertRaisesRegex(Pro20SourceFactoryError, "member order"):
            Pro20RuntimeOrigin("a" * 64, "b" * 64, {})


if __name__ == "__main__":
    unittest.main()
