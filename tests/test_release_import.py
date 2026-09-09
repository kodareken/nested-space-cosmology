"""Focused importer and nested-path checks for the v0.5.0 development records."""
from __future__ import annotations

import hashlib
import json
import subprocess
import unittest

from test_publication import ROOT, load_script

V040_SOURCE = "95b96be312feb667377cdbc3bbfe453697a458dd"
V050_SOURCE = "6eeff9bfab26a18fcd029a59ec908643b63386b5"
NEW_IMPORTS = (
    "docs/nsc-compact-boundary-action.md",
    "docs/nsc-compact-casimir.md",
    "docs/nsc-vacuum-charge-matching.md",
    "results/development/compact-boundary-action.json",
    "results/development/compact-casimir.json",
    "results/development/vacuum-charge-matching.json",
    "scripts/check_nsc_compact_boundary_action.py",
    "scripts/check_nsc_compact_casimir.py",
    "scripts/check_nsc_vacuum_charge_matching.py",
    "src/recursive_horizons/nsc_compact_casimir.py",
)


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
            "results/development/vacuum-charge-matching.json",
            "results/development/compact-boundary-action.json",
            "results/development/compact-casimir.json",
        ):
            self.assertEqual(relative, self.importer.safe_relative(relative))

    def test_v050_import_files_match_lab_pin_and_local_bytes(self):
        self.assertEqual("0.5.0", self.release["release_version"])
        self.assertEqual(V050_SOURCE, self.release["source_commit"])
        self.assertEqual(
            "d8f89545c1476c4b5cc862ee2e4271eb89d244a7",
            self.release["checkpoint_source_commit"],
        )
        self.assertEqual(10, len(self.release["import_files"]))
        self.assertEqual(17, len(self.release["preserved_85_scientific_files"]))
        self.assertEqual(83, len(self.release["preserved_81_scientific_files"]))
        self.assertEqual({entry["path"] for entry in self.release["import_files"]}, set(NEW_IMPORTS))
        for entry in self.release["import_files"] + self.release["preserved_85_scientific_files"]:
            path = ROOT / entry["path"]
            self.assertTrue(path.is_file(), entry["path"])
            self.assertEqual(
                entry["sha256"],
                hashlib.sha256(path.read_bytes()).hexdigest(),
                entry["path"],
            )

    def test_v040_nested_records_are_preserved_not_reimported(self):
        preserved = {entry["path"] for entry in self.release["preserved_85_scientific_files"]}
        imported = {entry["path"] for entry in self.release["import_files"]}
        for relative in (
            "results/development/charged-sector.json",
            "results/development/compact-interaction.json",
            "results/development/flow-compatibility.json",
            "results/development/torsion-uv-map.json",
        ):
            self.assertIn(relative, preserved)
            self.assertNotIn(relative, imported)

    def test_v040_scientific_bytes_match_tagged_snapshot(self):
        original = json.loads(
            subprocess.check_output(
                ["git", "show", "v0.4.0:results/manifest.json"],
                cwd=ROOT,
            )
        )
        self.assertEqual(85, original["result_count"])
        self.assertEqual(V040_SOURCE, original["release_source_commit"])
        previous = {step["output"]: step for step in original["steps"]}
        current = {step["output"]: step for step in json.loads(
            (ROOT / "results/manifest.json").read_text()
        )["steps"]}
        self.assertEqual(set(previous), set(current) - {
            "results/development/vacuum-charge-matching.json",
            "results/development/compact-boundary-action.json",
            "results/development/compact-casimir.json",
        })
        for output, step in previous.items():
            now = current[output]
            self.assertEqual(step["output_sha256"], now["output_sha256"], output)
            self.assertEqual(step["generator_sha256"], now["generator_sha256"], output)
            self.assertEqual(step["comparison_policy"], now["comparison_policy"], output)
            self.assertEqual(step.get("source_dependency_hashes"), now.get("source_dependency_hashes"), output)

    def test_new_record_dependency_closure_has_forty_one_existing_paths(self):
        seen = set()
        files = set()

        def walk(path):
            if path in seen:
                return
            seen.add(path)
            files.add(path)
            value = json.loads((ROOT / path).read_text())
            for field in ("source_hashes", "input_hashes"):
                for relative in value.get(field, {}):
                    files.add(relative)
                    if relative.endswith(".json") and relative.startswith("results/"):
                        walk(relative)

        for path in (
            "results/development/vacuum-charge-matching.json",
            "results/development/compact-boundary-action.json",
            "results/development/compact-casimir.json",
        ):
            walk(path)
        self.assertEqual(41, len(files))
        for relative in sorted(files):
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            record = json.loads(path.read_text()) if relative.endswith(".json") else None
            if record and relative in seen:
                for field in ("source_hashes", "input_hashes"):
                    for nested, expected in record.get(field, {}).items():
                        self.assertEqual(expected, hashlib.sha256((ROOT / nested).read_bytes()).hexdigest(), nested)

    def test_compact_mass_checker_is_pinned_but_not_a_result(self):
        paths = {entry["path"] for entry in self.release["preserved_85_scientific_files"]}
        self.assertIn("scripts/check_nsc_compact_mass_map.py", paths)
        self.assertIn("docs/nsc-compact-mass-map.md", paths)
        outputs = {row["output"] for row in self.release["scoped_follow_ups"]}
        self.assertNotIn("results/development/compact-mass-map.json", outputs)
        self.assertTrue(
            all("compact-mass" not in row["output"] for row in self.release["scoped_follow_ups"])
        )


if __name__ == "__main__":
    unittest.main()
