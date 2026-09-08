"""Focused importer and nested-path checks for the v0.4.0 development records."""
from __future__ import annotations

import hashlib
import json
import unittest

from test_publication import ROOT, load_script


class ReleaseImportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.importer = load_script("stage_release_import.py")
        cls.release = json.loads((ROOT / "results/release-spec.json").read_text())

    def test_nested_development_paths_are_explicitly_allowed(self):
        for relative in (
            "results/development/compact-interaction.json",
            "results/development/torsion-uv-map.json",
            "results/development/flow-compatibility.json",
            "results/development/charged-sector.json",
        ):
            self.assertEqual(relative, self.importer.safe_relative(relative))

    def test_v040_import_files_match_lab_pin_and_local_bytes(self):
        self.assertEqual("0.4.0", self.release["release_version"])
        self.assertEqual(
            "95b96be312feb667377cdbc3bbfe453697a458dd",
            self.release["source_commit"],
        )
        self.assertEqual(17, len(self.release["import_files"]))
        self.assertEqual(83, len(self.release["preserved_81_scientific_files"]))
        for entry in self.release["import_files"]:
            path = ROOT / entry["path"]
            self.assertTrue(path.is_file(), entry["path"])
            self.assertEqual(
                entry["sha256"],
                hashlib.sha256(path.read_bytes()).hexdigest(),
                entry["path"],
            )

    def test_compact_mass_checker_is_pinned_but_not_a_result(self):
        paths = {entry["path"] for entry in self.release["import_files"]}
        self.assertIn("scripts/check_nsc_compact_mass_map.py", paths)
        self.assertIn("docs/nsc-compact-mass-map.md", paths)
        outputs = {row["output"] for row in self.release["scoped_follow_ups"]}
        self.assertNotIn("results/development/compact-mass-map.json", outputs)
        self.assertTrue(
            all("compact-mass" not in row["output"] for row in self.release["scoped_follow_ups"])
        )


if __name__ == "__main__":
    unittest.main()
