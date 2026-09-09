"""v0.6.0 identity, count, policy and metadata agreement checks."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import tomllib
import unittest

from test_publication import ROOT, load_script


V050_SOURCE = "6eeff9bfab26a18fcd029a59ec908643b63386b5"
V060_SOURCE = "eea43512f63e61d8d7686261a683bc36e99efa19"
NEW_OUTPUTS = (
    "results/development/horizon-source.json",
    "results/development/warped-source.json",
    "results/development/compact-matching.json",
)
NEW_STEPS = {
    "NSC-19-HORIZON-SOURCE": {
        "output": "results/development/horizon-source.json",
        "schema": "NSC-HORIZON-SOURCE-v1",
        "claim": "derived-horizon-source",
        "relative_tolerance": 3e-13,
        "absolute_tolerance": 3e-13,
    },
    "NSC-20-WARPED-SOURCE": {
        "output": "results/development/warped-source.json",
        "schema": "NSC-WARPED-SOURCE-v1",
        "claim": "derived-warped-source",
        "relative_tolerance": 3e-08,
        "absolute_tolerance": 3e-09,
    },
    "NSC-21-COMPACT-MATCHING": {
        "output": "results/development/compact-matching.json",
        "schema": "NSC-COMPACT-MATCHING-v1",
        "claim": "derived-compact-matching",
        "relative_tolerance": 3e-08,
        "absolute_tolerance": 3e-09,
    },
}


class V060ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer = load_script("reproduce_public_results.py")
        cls.manifest = json.loads((ROOT / "results/manifest.json").read_text())
        cls.release = json.loads((ROOT / "results/release-spec.json").read_text())
        cls.steps = {row["artifact_id"]: row for row in cls.manifest["steps"]}
        cls.paper = json.loads((ROOT / "paper/metadata.json").read_text())
        cls.project = tomllib.loads((ROOT / "pyproject.toml").read_text())
        cls.citation = (ROOT / "CITATION.cff").read_text()

    def test_package_citation_paper_and_release_versions_agree(self):
        self.assertEqual("0.6.0", self.project["project"]["version"])
        self.assertEqual("0.6.0", self.paper["version"])
        self.assertEqual("0.6.0", self.release["release_version"])
        self.assertEqual("0.6.0", self.manifest["release_version"])
        self.assertEqual("9 September 2026", self.paper["date"])
        self.assertEqual("2026-09-09", self.paper["date_iso"])
        self.assertEqual(1788912000, self.paper["source_date_epoch"])
        self.assertRegex(self.citation, r"(?m)^date-released: 2026-09-09$")
        match = re.search(r"(?m)^version:\s*(.+)$", self.citation)
        self.assertIsNotNone(match)
        self.assertEqual("0.6.0", match.group(1).strip().strip("'\""))
        self.assertNotIn('version: 0.6.0"', self.citation)
        self.assertNotIn('version: 0.5.0"', self.citation)

    def test_result_count_is_fifty_eight_historical_plus_thirty_three_follow_ups(self):
        self.assertEqual(58, self.manifest["historical_result_count"])
        self.assertEqual(33, len(self.release["scoped_follow_ups"]))
        self.assertEqual(91, self.manifest["result_count"])
        self.assertEqual(91, len(self.manifest["steps"]))
        self.assertEqual(58 + 33, len(self.manifest["steps"]))
        self.assertEqual(V060_SOURCE, self.release["source_commit"])
        self.assertEqual(V060_SOURCE, self.manifest["release_source_commit"])

    def test_previous_eighty_eight_manifest_steps_are_byte_identical(self):
        original = json.loads(
            subprocess.check_output(["git", "show", "v0.5.0:results/manifest.json"], cwd=ROOT)
        )
        self.assertEqual(88, original["result_count"])
        self.assertEqual(V050_SOURCE, original["release_source_commit"])
        current = self.manifest["steps"]
        self.assertEqual(original["steps"], current[:88])
        self.assertEqual(NEW_OUTPUTS, tuple(row["output"] for row in current[88:]))

    def test_new_records_use_schema_status_all_fields_and_source_closure(self):
        for artifact, expected in NEW_STEPS.items():
            step = self.steps[artifact]
            value = json.loads((ROOT / step["output"]).read_text())
            self.assertEqual(expected["output"], step["output"])
            self.assertEqual("schema_status", step["identity_policy"]["kind"])
            self.assertEqual(expected["schema"], step["identity_policy"]["schema"])
            self.assertEqual(expected["schema"], value["schema"])
            self.assertEqual(step["identity_policy"]["status"], value["status"])
            self.assertNotIn("artifact_id", value)
            self.assertNotIn("terminal", value)
            self.assertNotIn("exceptions", step["comparison_policy"])
            self.assertEqual("all_fields", step["comparison_policy"]["kind"])
            self.assertEqual(
                (expected["relative_tolerance"], expected["absolute_tolerance"]),
                (
                    step["comparison_policy"]["relative_tolerance"],
                    step["comparison_policy"]["absolute_tolerance"],
                ),
            )
            self.assertEqual([expected["claim"]], step["paper_claim_ids"])
            self.assertEqual(V060_SOURCE, step["follow_up_source_commit"])
            self.assertEqual(set(step["source_dependencies"]), set(step["source_dependency_hashes"]))
            self.assertTrue(step["source_dependency_hashes"])
            self.reproducer.validate_identity(value, step)
            locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
            self.assertTrue(all(location[0] == "input_hashes" for location in locations))
            for relative, digest in step["source_dependency_hashes"].items():
                self.assertEqual(digest, hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), relative)
            for relative, digest in value["source_hashes"].items():
                self.assertEqual(digest, hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), relative)
            for relative, digest in value["input_hashes"].items():
                self.assertEqual(digest, hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), relative)
                self.assertIn(relative, step["dependencies"])

    def test_six_physical_targets_remain_open_with_scoped_support(self):
        self.assertEqual(6, len(self.release["current_targets"]))
        by_id = {row["claim_id"]: row for row in self.release["current_targets"]}
        for claim in (
            "target-constants",
            "target-continuation",
            "target-dark-sector",
            "target-measurement",
            "target-antimatter",
            "target-recursion",
        ):
            self.assertEqual("open_target", by_id[claim]["status"])
            self.assertLessEqual(set(by_id[claim]["supporting_artifact_ids"]), set(self.steps))
        self.assertIn("NSC-19-HORIZON-SOURCE", by_id["target-continuation"]["supporting_artifact_ids"])
        self.assertIn("NSC-20-WARPED-SOURCE", by_id["target-constants"]["supporting_artifact_ids"])
        self.assertIn("NSC-21-COMPACT-MATCHING", by_id["target-dark-sector"]["supporting_artifact_ids"])
        self.assertNotIn("closed", json.dumps(self.release["current_targets"]))


if __name__ == "__main__":
    unittest.main()
