"""Focused identity, hash and comparison checks for nested development records."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

from test_publication import ROOT, load_script


DEVELOPMENT = (
    "NSC-12-COMPACT-INTERACTION",
    "NSC-13-TORSION-UV-MAP",
    "NSC-14-FLOW-COMPATIBILITY",
    "NSC-15-CHARGED-SECTOR",
)
NEW_DEVELOPMENT = (
    "NSC-16-VACUUM-CHARGE-MATCHING",
    "NSC-17-COMPACT-BOUNDARY-ACTION",
    "NSC-18-COMPACT-CASIMIR",
)
V060_DEVELOPMENT = (
    "NSC-19-HORIZON-SOURCE",
    "NSC-20-WARPED-SOURCE",
    "NSC-21-COMPACT-MATCHING",
)


class ReproducePublicResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer = load_script("reproduce_public_results.py")
        cls.manifest = cls.reproducer.load_manifest()
        cls.steps = {row["artifact_id"]: row for row in cls.manifest["steps"]}

    def test_schema_status_assigns_manifest_ids_without_json_identity(self):
        for artifact in DEVELOPMENT + NEW_DEVELOPMENT + V060_DEVELOPMENT:
            step = self.steps[artifact]
            value = json.loads((ROOT / step["output"]).read_text())
            self.assertNotIn("artifact_id", value)
            self.assertNotIn("terminal", value)
            self.assertEqual("schema_status", step["identity_policy"]["kind"])
            self.reproducer.validate_identity(value, step)
            self.assertEqual(step["artifact_id"], artifact)
            self.assertTrue(step["output"].startswith("results/development/"))
            wrong = copy.deepcopy(value)
            wrong["status"] = value["status"] + " altered"
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_identity(wrong, step)
            injected = copy.deepcopy(value)
            injected["artifact_id"] = artifact
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_identity(injected, step)

    def test_source_and_input_hashes_are_not_digest_normalized(self):
        for artifact in DEVELOPMENT:
            step = self.steps[artifact]
            value = json.loads((ROOT / step["output"]).read_text())
            locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
            self.assertEqual(set(), locations)
            normalized = self.reproducer.normalize_dynamic_hashes(value, locations)
            self.assertEqual(value["source_hashes"], normalized["source_hashes"])
            if "input_hashes" in value:
                self.assertEqual(value["input_hashes"], normalized["input_hashes"])
            tampered = copy.deepcopy(value)
            key = next(iter(tampered["source_hashes"]))
            tampered["source_hashes"][key] = "0" * 64
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_authenticated_inputs(ROOT, tampered, step)

    def test_new_records_normalize_only_result_input_digests(self):
        expected_counts = {
            "NSC-16-VACUUM-CHARGE-MATCHING": 3,
            "NSC-17-COMPACT-BOUNDARY-ACTION": 2,
            "NSC-18-COMPACT-CASIMIR": 2,
            "NSC-19-HORIZON-SOURCE": 3,
            "NSC-20-WARPED-SOURCE": 3,
            "NSC-21-COMPACT-MATCHING": 4,
        }
        for artifact, count in expected_counts.items():
            step = self.steps[artifact]
            value = json.loads((ROOT / step["output"]).read_text())
            locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
            self.assertEqual(count, len(locations))
            self.assertTrue(all(location[0] == "input_hashes" for location in locations))
            normalized = self.reproducer.normalize_dynamic_hashes(value, locations)
            self.assertEqual(value["source_hashes"], normalized["source_hashes"])
            for path, digest in value["input_hashes"].items():
                actual = normalized["input_hashes"][path]
                self.assertEqual(
                    "<validated-generated-result>" if path.startswith("results/") else digest,
                    actual,
                )
            tampered = copy.deepcopy(value)
            key = next(iter(tampered["source_hashes"]))
            tampered["source_hashes"][key] = "0" * 64
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_authenticated_inputs(ROOT, tampered, step)
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.normalize_dynamic_hashes(
                    value, {("source_hashes", next(iter(value["source_hashes"])))}
                )

    def test_recorded_all_field_policies_are_preserved(self):
        compact = self.steps["NSC-12-COMPACT-INTERACTION"]
        self.assertEqual(
            (3e-12, 3e-12),
            (
                compact["comparison_policy"]["relative_tolerance"],
                compact["comparison_policy"]["absolute_tolerance"],
            ),
        )
        for artifact in (
            "NSC-13-TORSION-UV-MAP",
            "NSC-14-FLOW-COMPATIBILITY",
            "NSC-15-CHARGED-SECTOR",
        ):
            policy = self.steps[artifact]["comparison_policy"]
            self.assertEqual("all_fields", policy["kind"])
            self.assertEqual(0.0, policy["relative_tolerance"])
            self.assertEqual(0.0, policy["absolute_tolerance"])
        for artifact in (
            "NSC-16-VACUUM-CHARGE-MATCHING",
            "NSC-17-COMPACT-BOUNDARY-ACTION",
        ):
            policy = self.steps[artifact]["comparison_policy"]
            self.assertEqual("all_fields", policy["kind"])
            self.assertEqual((3e-13, 3e-13), (policy["relative_tolerance"], policy["absolute_tolerance"]))
            self.assertNotIn("exceptions", policy)
        casimir = self.steps["NSC-18-COMPACT-CASIMIR"]["comparison_policy"]
        self.assertEqual("all_fields", casimir["kind"])
        self.assertEqual((3e-08, 3e-09), (casimir["relative_tolerance"], casimir["absolute_tolerance"]))
        self.assertNotIn("exceptions", casimir)
        horizon = self.steps["NSC-19-HORIZON-SOURCE"]["comparison_policy"]
        self.assertEqual("all_fields", horizon["kind"])
        self.assertEqual((3e-13, 3e-13), (horizon["relative_tolerance"], horizon["absolute_tolerance"]))
        self.assertNotIn("exceptions", horizon)
        for artifact in ("NSC-20-WARPED-SOURCE", "NSC-21-COMPACT-MATCHING"):
            policy = self.steps[artifact]["comparison_policy"]
            self.assertEqual("all_fields", policy["kind"])
            self.assertEqual((3e-08, 3e-09), (policy["relative_tolerance"], policy["absolute_tolerance"]))
            self.assertNotIn("exceptions", policy)

    def test_v030_status_strings_remain_historically_true(self):
        compact = json.loads(
            (ROOT / self.steps["NSC-12-COMPACT-INTERACTION"]["output"]).read_text()
        )
        self.assertEqual(
            "development calculation; not part of the v0.3.0 release collection",
            compact["status"],
        )
        charged = json.loads(
            (ROOT / self.steps["NSC-15-CHARGED-SECTOR"]["output"]).read_text()
        )
        self.assertEqual(
            ["(I8-tau3 tensor gamma5)/2", "(I8+tau3 tensor gamma5)/2"],
            charged["allowed_zero_projectors"],
        )
        self.assertFalse(charged["gauge_heat_coefficient"]["total_renormalized_gauge_coupling_fixed"])

    def test_selected_only_parser_accepts_nested_outputs(self):
        selected = self.reproducer.selected_steps(
            self.manifest,
            "results/development/charged-sector.json,NSC-12-COMPACT-INTERACTION,NSC-18-COMPACT-CASIMIR,NSC-21-COMPACT-MATCHING",
        )
        self.assertEqual(
            {
                "results/development/compact-interaction.json",
                "results/development/charged-sector.json",
                "results/development/compact-casimir.json",
                "results/development/compact-matching.json",
            },
            {row["output"] for row in selected},
        )
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.selected_steps(self.manifest, "results/missing.json")

    def test_unselected_outputs_are_seeded_for_recursive_authentication(self):
        selected = {"results/development/compact-matching.json"}
        with tempfile.TemporaryDirectory() as workspace:
            expected = Path(workspace) / "expected"
            work = Path(workspace) / "work"
            for step in self.manifest["steps"]:
                self.reproducer.copy_relative(ROOT, expected, step["output"])
            complete = self.reproducer.seed_unselected_outputs(
                work, expected, self.manifest, selected
            )
            self.assertNotIn("results/development/compact-matching.json", complete)
            for relative in (
                "results/development/warped-source.json",
                "results/development/horizon-source.json",
                "results/development/compact-casimir.json",
                "results/development/compact-boundary-action.json",
                "results/development/vacuum-charge-matching.json",
                "results/development/charged-sector.json",
                "results/nsc-9-covariant-source.json",
                "results/nsc-8-finite-terms.json",
                "results/nsc-4-covariant-measure.json",
            ):
                self.assertIn(relative, complete)
                self.assertEqual(
                    (ROOT / relative).read_bytes(),
                    (work / relative).read_bytes(),
                )

    def test_empty_explicit_selection_never_triggers_full_reproduction(self):
        for selection in ("", " ", ",,", " , "):
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.selected_steps(self.manifest, selection)


if __name__ == "__main__":
    unittest.main()
