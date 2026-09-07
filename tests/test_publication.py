from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "results" / "manifest.json"


def load_script(name: str):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name.removesuffix(".py"), path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicationManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        cls.steps = cls.manifest["steps"]

    def test_exact_frontier_and_count(self) -> None:
        self.assertEqual(58, self.manifest["result_count"])
        self.assertEqual(58, len(self.steps))
        self.assertEqual(
            "NSC-2-ZETA1-RECURSION-MAP", self.manifest["frontier_artifact_id"]
        )
        current = [
            step["artifact_id"]
            for step in self.steps
            if step["category"] == "current_frontier"
        ]
        self.assertEqual(["NSC-2-ZETA1-RECURSION-MAP"], current)

    def test_manifest_is_topological(self) -> None:
        complete: set[str] = set()
        for step in self.steps:
            self.assertLessEqual(set(step["dependencies"]), complete)
            complete.add(step["output"])

    def test_declared_hashes_match_checkout(self) -> None:
        for step in self.steps:
            for path_key, hash_key in (
                ("output", "output_sha256"),
                ("generator", "generator_sha256"),
            ):
                path = ROOT / step[path_key]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(
                    step[hash_key], hashlib.sha256(path.read_bytes()).hexdigest()
                )

    def test_every_result_is_terminal_and_identified(self) -> None:
        for step in self.steps:
            result = json.loads((ROOT / step["output"]).read_text(encoding="utf-8"))
            self.assertIs(True, result["terminal"])
            self.assertEqual(step["artifact_id"], result["artifact_id"])

    def test_determinant_scale_roots_are_not_frontier(self) -> None:
        superseded = {
            step["artifact_id"]
            for step in self.steps
            if step["category"] == "superseded_candidate"
        }
        self.assertTrue(
            {
                "NSC-2-ZETA1-LOWEST-MODE",
                "NSC-2-ZETA1-REGULATED-DETERMINANT",
                "NSC-2-ZETA1-ORBIFOLD-PARITY",
                "NSC-2-ZETA1-WARPED-Y",
            }
            <= superseded
        )

    def test_known_black_universe_is_imported(self) -> None:
        categories = {step["artifact_id"]: step["category"] for step in self.steps}
        self.assertEqual(
            "imported_benchmark_reproduction",
            categories["NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING"],
        )

    def test_frontier_preserves_open_scale_claims(self) -> None:
        frontier = json.loads(
            (ROOT / self.manifest["frontier_output"]).read_text(encoding="utf-8")
        )
        self.assertIs(False, frontier["gate"]["mode_resolved_tail_solved"])
        self.assertIs(False, frontier["gate"]["zeta_derived"])
        self.assertIs(False, frontier["nonclaims"]["Omega_value_selected"])
        self.assertIs(False, frontier["nonclaims"]["physical_zeta_promoted"])

    def test_manifest_builder_round_trip(self) -> None:
        builder = load_script("build_result_manifest.py")
        self.assertEqual(self.manifest, builder.build())

    def test_numeric_comparison_policy_propagates_through_dynamic_steps(self) -> None:
        policies = {
            step["artifact_id"]: step["comparison_policy"]["kind"]
            for step in self.steps
        }
        self.assertEqual(
            "portable_numeric",
            policies["NSC-1-S-ONE-FULL-EXPONENTIAL-RELATIVE-KERNEL"],
        )
        self.assertEqual(
            "portable_numeric",
            policies["NSC-1-S-ONE-GAP-TRANSITION-CLOSURE"],
        )
        self.assertEqual(
            "portable_numeric",
            policies["NSC-2-ZETA1-ANOMALY-DECOMPOSITION"],
        )

    def test_reproduction_validator_accepts_frozen_checkout(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        reproducer.validate_checkout(self.manifest)

    def test_all_authenticated_field_variants_are_validated_and_normalized(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        for name in (
            "nsc-1-s-one-spectral-room-bind.json",
            "nsc-1-s-one-shadow-lensing-observation.json",
            "nsc-1-s-one-resolution-channel-closure.json",
            "nsc-1-s-one-flat-spectral-poles.json",
        ):
            value = json.loads((ROOT / "results" / name).read_text(encoding="utf-8"))
            reproducer.validate_authenticated_inputs(ROOT, value)
            normalized = reproducer.normalize_dynamic_hashes(value)
            entries = list(reproducer.authenticated_entries(normalized))
            self.assertTrue(entries)
            self.assertTrue(
                all(
                    entry["sha256"] == "<validated-generated-result>"
                    for entry in entries
                    if entry["path"].startswith("results/")
                )
            )

    def test_portable_comparison_ignores_only_undeclared_numerics(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        expected = {
            "classification": "same",
            "gate": {"passes": True},
            "scan": 5.0,
            "argmin_index": 2,
        }
        actual = {
            "classification": "same",
            "gate": {"passes": True},
            "scan": 9.0,
            "argmin_index": 3,
        }
        reproducer.compare_portable(
            expected,
            actual,
            path="",
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
            compare_numbers=False,
        )
        with self.assertRaises(reproducer.ReproductionError):
            reproducer.compare_portable(
                expected,
                {
                    "classification": "same",
                    "gate": {"passes": False},
                    "scan": 5.0,
                    "argmin_index": 2,
                },
                path="",
                relative_tolerance=1e-8,
                absolute_tolerance=1e-10,
                compare_numbers=False,
            )

    def test_headline_pointer_is_resolved(self) -> None:
        reproducer = load_script("reproduce_public_results.py")
        value = {"outer": {"rows": [{"zeta": 4.75}]}}
        self.assertEqual(4.75, reproducer.resolve_pointer(value, "/outer/rows/0/zeta"))

    def test_headline_observables_are_granular_numeric_paths(self) -> None:
        for step in self.steps:
            for observable in step["headline_observables"]:
                self.assertEqual(
                    {"pointer", "relative_tolerance", "absolute_tolerance"},
                    set(observable),
                )
                value = json.loads(
                    (ROOT / step["output"]).read_text(encoding="utf-8")
                )
                pointer_value = load_script("reproduce_public_results.py").resolve_pointer(
                    value, observable["pointer"]
                )
                self.assertIsInstance(pointer_value, (int, float))


if __name__ == "__main__":
    unittest.main()
