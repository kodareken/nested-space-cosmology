from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


def checker_module():
    spec = importlib.util.spec_from_file_location("pref27_check_repo_test", ROOT / "scripts/check_repo.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class Pref27CompactCheckerWalkerTests(unittest.TestCase):
    def test_live_store_observation_is_not_part_of_compact_default(self) -> None:
        checker = checker_module()
        with patch.object(
            checker,
            "_check_pro18_pref27_live_store",
            side_effect=AssertionError("ordinary compact verification reopened temporal evidence"),
        ):
            failures = checker._check_pro18_pref27_result()
        self.assertEqual(failures, [])

    def test_explicit_live_store_observation_rejects_successor_leaves(self) -> None:
        checker = checker_module()
        with TemporaryDirectory() as directory:
            root = Path(directory) / "repository"
            expected_root = "runs/fgc-2-sf1/proto17/calibration"
            store = root / expected_root
            (store / "states").mkdir(parents=True)
            expected_leaf = store / "states" / "generation-zero.json"
            expected_leaf.write_text("{}\n", encoding="utf-8")
            successor_leaf = store / "journal" / "000001-successor.journal"
            successor_leaf.parent.mkdir()
            successor_leaf.write_text("{}\n", encoding="utf-8")
            expected = {expected_leaf.relative_to(root).as_posix()}
            with patch.object(checker, "REPOSITORY", root):
                failures = checker._check_pro18_pref27_live_store(expected_root, expected)
        self.assertEqual(
            failures,
            ["FGC-1-PRO18-PREF27 ignored raw store is partially present or structurally drifted"],
        )

    def test_nested_directory_symlink_is_rejected_without_scanning_its_target(self) -> None:
        checker = checker_module()
        with TemporaryDirectory() as directory:
            root = Path(directory) / "repository"
            store = root / "runs/fgc-2-sf1/proto17/calibration"
            outside = root / "outside"
            (store / "states").mkdir(parents=True)
            outside.mkdir(parents=True)
            (outside / "must-not-be-enumerated").write_text("sentinel", encoding="utf-8")
            os.symlink(outside, store / "states" / "foreign-directory")
            calls: list[Path] = []
            original = checker.os.scandir

            def recording_scandir(path):
                calls.append(Path(path))
                return original(path)

            with patch.object(checker, "REPOSITORY", root), patch.object(
                checker.os, "scandir", side_effect=recording_scandir,
            ):
                leaves, unsafe = checker._pref27_nonfollowing_store_walk(store)
            self.assertTrue(unsafe)
            self.assertEqual(leaves, set())
            self.assertNotIn(outside, calls)

    def test_regular_tree_is_walked_lexically_without_reading_leaf_contents(self) -> None:
        checker = checker_module()
        with TemporaryDirectory() as directory:
            root = Path(directory) / "repository"
            store = root / "runs/fgc-2-sf1/proto17/calibration"
            (store / "zeta").mkdir(parents=True)
            (store / "alpha").mkdir()
            (store / "zeta" / "z.json").write_text("z", encoding="utf-8")
            (store / "alpha" / "a.json").write_text("a", encoding="utf-8")
            with patch.object(checker, "REPOSITORY", root):
                leaves, unsafe = checker._pref27_nonfollowing_store_walk(store)
            self.assertFalse(unsafe)
            self.assertEqual(leaves, {
                "runs/fgc-2-sf1/proto17/calibration/alpha/a.json",
                "runs/fgc-2-sf1/proto17/calibration/zeta/z.json",
            })


if __name__ == "__main__":
    unittest.main()
