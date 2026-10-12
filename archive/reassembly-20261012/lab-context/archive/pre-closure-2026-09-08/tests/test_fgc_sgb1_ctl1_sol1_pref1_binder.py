from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from recursive_horizons.fgc import sgb1_ctl1_sol1_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _script():
    path = ROOT / "scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py"
    spec = importlib.util.spec_from_file_location("sol1_pref1_script", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SOL1PREF1BinderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = (ROOT / binder.CONFIG_PATH).read_bytes()
        self.result = (ROOT / binder.RESULT_PATH).read_bytes()

    def test_tracked_bundle_is_canonical_and_pinned(self) -> None:
        result = binder.verify_compact(ROOT)
        self.assertEqual(_sha(self.config), binder.CONFIG_SHA256)
        self.assertEqual(_sha(self.result), binder.RESULT_SHA256)
        self.assertEqual(result["classification"], binder.CLASSIFICATION)

    def test_live_reconstruction_matches_tracked_bytes(self) -> None:
        config, result = binder.compose_canonical_artifacts(ROOT)
        self.assertEqual(config, self.config)
        self.assertEqual(result, self.result)

    def test_nominal_outcome_and_claims_are_exact(self) -> None:
        result = binder.verify_compact(ROOT)
        self.assertEqual(result["claims"], binder.expected_claims())
        self.assertEqual(result["nominal"]["classification"], "interval_inconclusive")
        self.assertEqual(
            result["nominal"]["obstruction"],
            "prefix_endpoint_not_inside_declared_tube",
        )
        self.assertEqual(result["nominal"]["unique_tubes"], 0)
        self.assertEqual(
            result["nominal"]["contract_sha256"], binder.NOMINAL_CONTRACT_SHA256
        )
        self.assertEqual(
            result["control"]["contract_sha256"], binder.CONTROL_CONTRACT_SHA256
        )

    def test_hash_reclassification_and_promotion_attacks_fail(self) -> None:
        with self.assertRaises(binder.SGBLSOL1PREF1Error):
            binder.validate_compact_result(self.config + b" ", self.result)
        parsed = json.loads(self.result)
        for operation in ("class", "obstruction", "unique", "claim", "policy"):
            changed = deepcopy(parsed)
            if operation == "class":
                changed["nominal"]["classification"] = "unique_local_affine_orbit"
            elif operation == "obstruction":
                changed["nominal"]["obstruction"] = "jacobian_diagonal_contains_zero"
            elif operation == "unique":
                changed["nominal"]["unique_tubes"] = 1
            elif operation == "claim":
                changed["claims"]["SGBL_branch_owned_and_healthy"] = True
            else:
                changed["policy"]["tube_radius"] = "1/16"
            raw = binder.canonical_result(changed)
            with patch.object(binder, "RESULT_SHA256", _sha(raw)):
                with self.assertRaises(binder.SGBLSOL1PREF1Error):
                    binder.validate_compact_result(self.config, raw)

    def test_freeze_commit_parent_delta_and_hash_attacks_fail(self) -> None:
        for name, value in (
            ("FREEZE_COMMIT", "0" * 40),
            ("FREEZE_PARENT", "1" * 40),
            ("FREEZE_CONFIG_SHA256", "2" * 64),
            ("FREEZE_RESULT_SHA256", "3" * 64),
        ):
            with patch.object(binder, name, value):
                with self.assertRaises(binder.SGBLSOL1PREF1Error):
                    binder.compose_canonical_artifacts(ROOT)
        changed_delta = binder.FREEZE_DELTA[:-1]
        with patch.object(binder, "FREEZE_DELTA", changed_delta):
            with self.assertRaises(binder.SGBLSOL1PREF1Error):
                binder.compose_canonical_artifacts(ROOT)

    def test_owner_and_contract_hash_attacks_fail(self) -> None:
        with patch.object(binder, "OWNER_PATHS", binder.OWNER_PATHS[:-1]):
            with self.assertRaises(binder.SGBLSOL1PREF1Error):
                binder.compose_canonical_artifacts(ROOT)
        with patch.object(binder, "NOMINAL_CONTRACT_SHA256", "4" * 64):
            with self.assertRaises(binder.SGBLSOL1PREF1Error):
                binder.compose_canonical_artifacts(ROOT)

    def test_compact_check_is_git_source_and_reconstruction_blind(self) -> None:
        with (
            patch.object(binder, "compose_canonical_artifacts", side_effect=AssertionError()),
            patch.object(binder, "git_read", side_effect=AssertionError()),
        ):
            binder.verify_compact(ROOT)

    def test_compact_bundle_passes_without_git_source_or_runs(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            for relative, raw in ((binder.CONFIG_PATH, self.config), (binder.RESULT_PATH, self.result)):
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(raw)
            binder.verify_compact(root)
            self.assertFalse((root / ".git").exists())
            self.assertFalse((root / "src").exists())
            self.assertFalse((root / "runs").exists())

    def test_default_script_never_reconstructs(self) -> None:
        module = _script()
        with patch.object(binder, "compose_canonical_artifacts") as compose:
            self.assertEqual(module.main([]), 0)
        compose.assert_not_called()

    def test_top_level_imports_exclude_frz1_and_sol1_decision_code(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_sol1_pref1_binder.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in tree.body:
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        self.assertFalse(any("sol1_frz1" in item or "sgb1_ctl1_sol1" in item for item in imports))
        self.assertFalse(any("campaign" in item or "run_fgc" in item for item in imports))

    def test_binder_module_contains_no_write_or_publish_call(self) -> None:
        source = (
            ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_sol1_pref1_binder.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("publish_exclusive_file", source)
        self.assertNotIn("write_bytes", source)


if __name__ == "__main__":
    unittest.main()
