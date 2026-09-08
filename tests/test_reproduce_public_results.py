"""Focused identity, hash and comparison checks for the four nested records."""
from __future__ import annotations

import copy
import json
import unittest

from test_publication import ROOT, load_script


DEVELOPMENT = (
    "NSC-12-COMPACT-INTERACTION",
    "NSC-13-TORSION-UV-MAP",
    "NSC-14-FLOW-COMPATIBILITY",
    "NSC-15-CHARGED-SECTOR",
)


class ReproducePublicResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer = load_script("reproduce_public_results.py")
        cls.manifest = cls.reproducer.load_manifest()
        cls.steps = {row["artifact_id"]: row for row in cls.manifest["steps"]}

    def test_schema_status_assigns_manifest_ids_without_json_identity(self):
        for artifact in DEVELOPMENT:
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
            "results/development/charged-sector.json,NSC-12-COMPACT-INTERACTION",
        )
        self.assertEqual(
            {
                "results/development/compact-interaction.json",
                "results/development/charged-sector.json",
            },
            {row["output"] for row in selected},
        )
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.selected_steps(self.manifest, "results/missing.json")

    def test_empty_explicit_selection_never_triggers_full_reproduction(self):
        for selection in ("", " ", ",,", " , "):
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.selected_steps(self.manifest, selection)


if __name__ == "__main__":
    unittest.main()
