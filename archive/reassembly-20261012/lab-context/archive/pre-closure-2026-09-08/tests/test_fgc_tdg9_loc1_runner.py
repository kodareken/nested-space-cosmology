"""Focused static and pure controls for the LOC1 read-only runner."""

from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction
import math
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_loc1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_loc1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema_independent as independent  # noqa: E402


class LOC1RunnerTests(unittest.TestCase):
    @staticmethod
    def _zero_fraction(value: object) -> bool:
        return value == {"numerator": "0", "denominator": "1"}

    def test_forbidden_imports_and_store_writes_are_absent(self) -> None:
        source = (ROOT / "scripts/run_fgc_tdg9_loc1.py").read_text()
        imported = {
            alias.name
            for node in ast.walk(ast.parse(source))
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("run_fgc_tdg9_ar1" in name for name in imported))
        self.assertFalse(any("tdg8_persisted_retry_replay" in name for name in imported))
        self.assertFalse(any("tdg9_ar1_pref1_binder" in name for name in imported))
        for forbidden in ("acquire_writer", "publish_checkpoint", "persist_snapshot", "SGB-L", "FGC-QR"):
            self.assertNotIn(forbidden, source)

    def test_exact_cubic_population_and_row_regions(self) -> None:
        segment = (0.0, 0.0, 1.0, 0.0, 0.25)
        rows = tuple((segment, (segment, segment), (segment,) * 4) for _ in range(authority.OWNED_ROW_COUNT))
        d01, d01_independent = runner._cubics(rows, "D01")  # type: ignore[attr-defined]
        d12, d12_independent = runner._cubics(rows, "D12")  # type: ignore[attr-defined]
        self.assertEqual((len(d01), len(d01_independent)), (4_088, 4_088))
        self.assertEqual((len(d12), len(d12_independent)), (8_176, 8_176))
        self.assertEqual(dict(d01[0].metadata)["global_index"], 1)
        self.assertEqual(dict(d01[0].metadata)["row_region"], "centre_parity_derivative_row")
        self.assertEqual(dict(d12[-1].metadata)["row_region"], "outer_projector_influenced")

    def test_publication_contains_only_manifest_and_terminal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})  # type: ignore[attr-defined]
            leaves = sorted((root / authority.OUTPUT_NAMESPACE).iterdir())
            self.assertEqual([item.name for item in leaves], ["manifest.json", "terminal.json"])

    def test_publication_rejects_symlinked_output_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            os.symlink(foreign, root / "runs" / "fgc-2-sf1")
            with self.assertRaisesRegex(runner.LOC1RunnerError, "unsafe_parent_component"):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})  # type: ignore[attr-defined]
            self.assertEqual(tuple(Path(foreign).iterdir()), ())

    def test_ulp_record_is_explanatory(self) -> None:
        record = runner._ulp_record(1.0)  # type: ignore[attr-defined]
        self.assertEqual(record["hex"], "0x1.0000000000000p+0")
        self.assertIn("upward", record)
        self.assertIn("downward", record)

    def test_route_agreement_rejects_candidate_count_and_classification_drift(self) -> None:
        coefficients = (Fraction(1, 24), Fraction(5, 12), Fraction(-3, 2), Fraction(1))
        first = primary.localize_absolute_maximum(
            [primary.LocalCubic(coefficients, {"ordinal": 0})],
            maximum_candidates=4,
            refinement_depth=96,
        )
        second = independent.localize_absolute_maximum_independently(
            [independent.IndependentLocalCubic(coefficients, {"ordinal": 0})],
            maximum_candidates=4,
            refinement_bits=128,
        )
        runner._require_route_agreement(first, second, "control")  # type: ignore[attr-defined]
        for changed in (
            replace(second, candidate_count=second.candidate_count - 1),
            replace(second, classification="unique_maximum"),
        ):
            with self.assertRaisesRegex(runner.LOC1RunnerError, "independent_disagreement"):
                runner._require_route_agreement(first, changed, "mutation")  # type: ignore[attr-defined]

    def test_right_half_linear_decomposition_uses_restricted_parent(self) -> None:
        outer = (0.0, 1.0, 1.0, 1.0, 1.0)
        medium = (
            (0.0, 1.0, 0.5, 1.0, 0.5),
            (0.5, 1.0, 1.0, 1.0, 0.5),
        )
        fine = tuple((index / 4, 1.0, (index + 1) / 4, 1.0, 0.25) for index in range(4))
        rows = ((outer, medium, fine),)
        candidate = SimpleNamespace(metadata=(("owned_row", 0), ("subinterval", 1)))
        evidence = runner._hermite_inputs(rows, "D01", candidate)  # type: ignore[attr-defined]
        self.assertTrue(all(self._zero_fraction(value) for value in evidence["endpoint_state_difference"].values()))
        self.assertTrue(all(self._zero_fraction(value) for value in evidence["width_scaled_endpoint_RHS_difference"].values()))
        self.assertTrue(all(self._zero_fraction(value) for value in evidence["complete_difference_polynomial"]))
        self.assertEqual(evidence["updates"]["restricted_parent"], evidence["updates"]["child"])

    def test_all_quarter_linear_decompositions_use_restricted_medium_parent(self) -> None:
        outer = (0.0, 1.0, 1.0, 1.0, 1.0)
        medium = (
            (0.0, 1.0, 0.5, 1.0, 0.5),
            (0.5, 1.0, 1.0, 1.0, 0.5),
        )
        fine = tuple((index / 4, 1.0, (index + 1) / 4, 1.0, 0.25) for index in range(4))
        rows = ((outer, medium, fine),)
        for subinterval in range(4):
            with self.subTest(subinterval=subinterval):
                candidate = SimpleNamespace(metadata=(("owned_row", 0), ("subinterval", subinterval)))
                evidence = runner._hermite_inputs(rows, "D12", candidate)  # type: ignore[attr-defined]
                self.assertTrue(all(self._zero_fraction(value) for value in evidence["endpoint_state_difference"].values()))
                self.assertTrue(all(self._zero_fraction(value) for value in evidence["width_scaled_endpoint_RHS_difference"].values()))
                self.assertTrue(all(self._zero_fraction(value) for value in evidence["complete_difference_polynomial"]))
                self.assertEqual(evidence["updates"]["restricted_parent"], evidence["updates"]["child"])

    def test_value_only_and_slope_only_components_are_exact(self) -> None:
        value_only = runner._component_coefficients(  # type: ignore[attr-defined]
            (0.0, 1.0, 1.0, 1.0, 1.0),
            (0.0, 1.0, 0.0, 1.0, 0.5),
            0,
        )
        self.assertEqual(value_only["slope_S"], (Fraction(0),) * 4)
        self.assertNotEqual(value_only["value_V"], (Fraction(0),) * 4)
        self.assertEqual(
            tuple(a + b for a, b in zip(value_only["value_V"], value_only["slope_S"], strict=True)),
            value_only["complete_C"],
        )
        slope_only = runner._component_coefficients(  # type: ignore[attr-defined]
            (0.0, 1.0, 1.0, 1.0, 1.0),
            (0.0, 0.0, 0.5, 0.0, 0.5),
            0,
        )
        self.assertEqual(slope_only["value_V"], (Fraction(0),) * 4)
        self.assertNotEqual(slope_only["slope_S"], (Fraction(0),) * 4)
        self.assertEqual(
            tuple(a + b for a, b in zip(slope_only["value_V"], slope_only["slope_S"], strict=True)),
            slope_only["complete_C"],
        )

    @staticmethod
    def _maximum(value: Fraction) -> SimpleNamespace:
        return SimpleNamespace(global_absolute_lower=value, global_absolute_upper=value)

    def test_contraction_boundary_and_one_bit_controls_apply_to_all_components(self) -> None:
        zero = self._maximum(Fraction(0))
        below = self._maximum(Fraction(*math.nextafter(1 / math.sqrt(8), 0.0).as_integer_ratio()))
        above = self._maximum(Fraction(*math.nextafter(1 / math.sqrt(8), math.inf).as_integer_ratio()))
        outer = self._maximum(Fraction(1))
        for component in ("value_V", "slope_S", "complete_C"):
            with self.subTest(component=component, case="exact_nonzero_equality"):
                self.assertEqual(
                    runner._classify_contraction_bounds(  # type: ignore[attr-defined]
                        Fraction(1), Fraction(1), Fraction(1), Fraction(1), exact_zero=False
                    ),
                    "sufficient_contraction_pass",
                )
            with self.subTest(component=component, case="rational_below"):
                self.assertEqual(
                    runner._classify_contraction_bounds(  # type: ignore[attr-defined]
                        Fraction(1), Fraction(1), Fraction(999, 1000), Fraction(999, 1000), exact_zero=False
                    ),
                    "sufficient_contraction_pass",
                )
            with self.subTest(component=component, case="rational_above"):
                self.assertEqual(
                    runner._classify_contraction_bounds(  # type: ignore[attr-defined]
                        Fraction(1), Fraction(1), Fraction(1001, 1000), Fraction(1001, 1000), exact_zero=False
                    ),
                    "sufficient_contraction_failure",
                )
            with self.subTest(component=component, case="exact_zero_precedence"):
                self.assertEqual(runner._contraction_assessment(zero, zero)["classification"], "exact_zero")  # type: ignore[attr-defined]
            with self.subTest(component=component, case="one_bit_below"):
                self.assertEqual(runner._contraction_assessment(outer, below)["classification"], "sufficient_contraction_pass")  # type: ignore[attr-defined]
            with self.subTest(component=component, case="one_bit_above"):
                self.assertEqual(runner._contraction_assessment(outer, above)["classification"], "sufficient_contraction_failure")  # type: ignore[attr-defined]

    def test_mixed_only_and_inconclusive_sole_ownership_are_preserved(self) -> None:
        maximum = self._maximum(Fraction(1))

        def component(classification: str) -> dict[str, object]:
            return {"D01": maximum, "D12": maximum, "contraction": {"classification": classification}}

        mixed = runner._component_ownership({  # type: ignore[attr-defined]
            "value_V": component("sufficient_contraction_pass"),
            "slope_S": component("sufficient_contraction_pass"),
            "complete_C": component("sufficient_contraction_failure"),
        })
        self.assertEqual(mixed["classification"], "mixed_only_failure")
        self.assertTrue(mixed["mixed_only_failure"])
        cancellation = runner._component_ownership({  # type: ignore[attr-defined]
            "value_V": {"D01": self._maximum(Fraction(2)), "D12": self._maximum(Fraction(2)), "contraction": {"classification": "threshold_inconclusive"}},
            "slope_S": {"D01": self._maximum(Fraction(2)), "D12": self._maximum(Fraction(2)), "contraction": {"classification": "threshold_inconclusive"}},
            "complete_C": {"D01": self._maximum(Fraction(1)), "D12": self._maximum(Fraction(1)), "contraction": {"classification": "threshold_inconclusive"}},
        })
        self.assertEqual(cancellation["cancellation_certified_by_level"], {"D01": True, "D12": True})
        for value_class, slope_class in (
            ("sufficient_contraction_failure", "threshold_inconclusive"),
            ("threshold_inconclusive", "sufficient_contraction_failure"),
        ):
            with self.subTest(value=value_class, slope=slope_class):
                result = runner._component_ownership({  # type: ignore[attr-defined]
                    "value_V": component(value_class),
                    "slope_S": component(slope_class),
                    "complete_C": component("sufficient_contraction_failure"),
                })
                self.assertEqual(result["classification"], "component_ownership_inconclusive")
                self.assertFalse(result["value_only_failure"])
                self.assertFalse(result["slope_only_failure"])


if __name__ == "__main__":
    unittest.main()
