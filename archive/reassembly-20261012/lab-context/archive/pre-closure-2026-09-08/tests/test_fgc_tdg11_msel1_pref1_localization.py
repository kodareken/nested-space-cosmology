"""Focused synthetic controls for the independent TDG11-MSEL1-PREF1 complete-C adapter.

These tests do not call the production complete-C adapter, reconstruction
owner, or runner decision logic.  Runtime execution is withheld while the
authorized FRZ1 diagnostic occupies the process gate.
"""

from __future__ import annotations

import ast
from dataclasses import fields, replace
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
    three_halves_order_passes_squared,
)
from recursive_horizons.fgc.evolution.tdg9_local_extrema import (  # noqa: E402
    localize_absolute_maximum as primary_localize,
)
from recursive_horizons.fgc.evolution.tdg9_local_extrema_independent_v2 import (  # noqa: E402
    IndependentLocalCubicV2,
    IndependentLocalizationEvidenceV2,
    RootIsolationInconclusive,
    localize_absolute_maximum_independently_v2 as independent_localize,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization import (  # noqa: E402
    IndependentCompleteCAdmissionEvidence,
    IndependentCompleteCClosedEvidence,
    IndependentCompleteCDifferenceEvidence,
    IndependentCompleteCResourceExhausted,
    IndependentCompleteCRouteDisagreement,
    PREF1_COMPLETE_C_EVALUATOR_ID,
    assess_complete_c_independently,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg11_msel1_pref1_localization as pref1_core,
)


Q = Fraction
MODULE_PATH = (
    ROOT / "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_localization.py"
)
PRODUCTION_ADAPTER_PATH = (
    ROOT / "src/recursive_horizons/fgc/evolution/tdg11_rational_complete_c.py"
)
ROW_HASH_DOMAIN = b"TDG11-RATIONAL-COMPLETE-C-ROW-STREAM-v1\n"
COEFFICIENT_HASH_DOMAIN = b"TDG11-RATIONAL-COMPLETE-C-CUBIC-STREAM-v1\n"
COMBINED_HASH_DOMAIN = b"TDG11-RATIONAL-COMPLETE-C-COMBINED-v1\n"
SURVIVOR_KEY_HASH_DOMAIN = b"TDG11-RATIONAL-COMPLETE-C-SURVIVOR-KEY-STREAM-v1\n"
STATIONARY_COUNT_DOMAIN = b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n"
TDG10_ROW_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-ROW-STREAM-v1\n"
TDG10_COEFFICIENT_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-CUBIC-STREAM-v1\n"
FORBIDDEN_IMPORT_MARKERS = (
    "tdg5_",
    "tdg6_temporal_admission_runtime",
    "tdg6_temporal_admission_theorem",
    "tdg7_",
    "tdg8_",
    "tdg9_ac1",
    "tdg9_ar1",
    "tdg9_loc1",
    "tdg9_loc2",
    "tdg9_ti1",
    "tdg9_ti2",
    "tdg9_ur1",
    "tdg10_exact_complete_c_runtime",
    "tdg10_qa",
    "tdg11_compensated",
    "tdg11_msel1_reconstruction",
    "tdg11_msel1_runtime",
    "tdg11_msel1_authority",
    "tdg11_msel1_contract",
    "tdg11_rational_complete_c",
    "assess_rational_complete_c_rows",
    "RationalCompleteCAdmissionEvidence",
    "binder",
    "scripts",
    "exact_hermite_coefficients",
    "numpy",
    "proto",
)


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module)
            names.update(alias.name for alias in node.names)
    return names


def _smoothstep(t: Fraction) -> Fraction:
    return 3 * t**2 - 2 * t**3


def _smoothstep_deriv(t: Fraction) -> Fraction:
    return 6 * t - 6 * t**2


def _smoothstep_segment(
    scale: Fraction, t0: Fraction, t1: Fraction, h: Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
    width = (t1 - t0) * h
    return (
        scale * _smoothstep(t0),
        scale * _smoothstep_deriv(t0) / h,
        scale * _smoothstep(t1),
        scale * _smoothstep_deriv(t1) / h,
        width,
    )


def _smoothstep_row(
    a: object, b: object, *, h: object = 1
) -> tuple[object, object, object]:
    amplitude_01 = Q(a)
    amplitude_12 = Q(b)
    width = Q(h)
    outer = _smoothstep_segment(Q(0), Q(0), Q(1), width)
    medium = (
        _smoothstep_segment(-amplitude_01, Q(0), Q(1, 2), width),
        _smoothstep_segment(-amplitude_01, Q(1, 2), Q(1), width),
    )
    fine = tuple(
        _smoothstep_segment(
            -(amplitude_01 + amplitude_12), Q(index, 4), Q(index + 1, 4), width
        )
        for index in range(4)
    )
    return (outer, medium, fine)


def _zero_row(*, h: object = 1) -> tuple[object, object, object]:
    return _smoothstep_row(0, 0, h=h)


def _hidden_interior_row() -> tuple[object, object, object]:
    """Global cubic T^2-T^3: its interior maximum exceeds all D01 endpoints."""

    def value(time: Fraction) -> Fraction:
        return time**2 - time**3

    def deriv(time: Fraction) -> Fraction:
        return 2 * time - 3 * time**2

    def sample(t0: Fraction, t1: Fraction) -> tuple[Fraction, ...]:
        return (-value(t0), -deriv(t0), -value(t1), -deriv(t1), t1 - t0)

    outer = (Q(0), Q(0), Q(0), Q(0), Q(1))
    medium = (sample(Q(0), Q(1, 2)), sample(Q(1, 2), Q(1)))
    fine = tuple(sample(Q(index, 4), Q(index + 1, 4)) for index in range(4))
    return (outer, medium, fine)


def _assess(
    row: object,
    *,
    expected: int = 1,
    d01: int = 64,
    d12: int = 64,
    depth: int = 32,
):
    return assess_complete_c_independently(
        (row,),
        expected_row_count=expected,
        maximum_candidates_D01=d01,
        maximum_candidates_D12=d12,
        refinement_depth=depth,
    )


def _row_digest(rows: tuple[object, ...]) -> str:
    hasher = sha256(ROW_HASH_DOMAIN)
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        values = [Q(value) for segment in (outer, *medium, *fine) for value in segment]
        hasher.update(
            (
                f"{ordinal}|"
                + "|".join(f"{item.numerator}/{item.denominator}" for item in values)
                + "\n"
            ).encode("ascii")
        )
    return hasher.hexdigest()


def _coefficient_digest(coefficients: tuple[tuple[Fraction, ...], ...]) -> str:
    hasher = sha256(COEFFICIENT_HASH_DOMAIN)
    for ordinal, cubic in enumerate(coefficients):
        hasher.update(
            (
                f"{ordinal}|"
                + "|".join(f"{item.numerator}/{item.denominator}" for item in cubic)
                + "\n"
            ).encode("ascii")
        )
    return hasher.hexdigest()


def _combined_digest(outer: str, finest: str) -> str:
    return sha256(
        COMBINED_HASH_DOMAIN + f"{outer}\n{finest}\n".encode("ascii")
    ).hexdigest()


def _mutate_segment(segment: object, index: int, value: object) -> tuple[object, ...]:
    items = list(segment)  # type: ignore[arg-type]
    items[index] = value
    return tuple(items)


def _exact_flags(d01_lower: Fraction, d01_upper: Fraction, d12_lower: Fraction, d12_upper: Fraction) -> tuple[bool, bool, bool]:
    left = 8 * d12_upper**2
    right = d01_lower**2
    exact_zero = d01_upper == 0 and d12_upper == 0
    if exact_zero:
        return False, False, False
    passed = left <= right
    failed = (not passed) and (8 * d12_lower**2 > d01_upper**2)
    return passed, failed, (not passed) and (not failed)


class TDG11MSEL1PREF1LocalizationTests(unittest.TestCase):
    def test_module_has_no_production_adapter_reconstruction_or_runner_imports(self) -> None:
        imported = _imported_modules(MODULE_PATH)
        source = MODULE_PATH.read_text()
        self.assertIn("tdg6_temporal_admission_design", imported)
        self.assertIn("classify_tdg6_channel", imported)
        self.assertIn("three_halves_order_passes_squared", imported)
        self.assertIn("localize_absolute_maximum", imported)
        self.assertIn("localize_absolute_maximum_independently_v2", imported)
        self.assertNotIn("tdg11_rational_complete_c", imported)
        self.assertNotIn("tdg11_msel1_reconstruction", imported)
        self.assertNotIn("assess_rational_complete_c_rows", imported)
        self.assertNotIn("exact_hermite_coefficients", imported)
        self.assertNotIn("restrict_exact_cubic_to_half", imported)
        self.assertNotIn("subtract_exact_cubics", imported)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )
        self.assertNotIn("import scripts", source)
        self.assertNotIn("from scripts", source)
        self.assertIn("assess_complete_c_independently", source)
        self.assertIn("TDG11-RATIONAL-COMPLETE-C-ROW-STREAM-v1", source)
        production_source = PRODUCTION_ADAPTER_PATH.read_text()
        for domain in (
            "TDG11-RATIONAL-COMPLETE-C-ROW-STREAM-v1",
            "TDG11-RATIONAL-COMPLETE-C-CUBIC-STREAM-v1",
            "TDG11-RATIONAL-COMPLETE-C-COMBINED-v1",
            "TDG11-RATIONAL-COMPLETE-C-SURVIVOR-KEY-STREAM-v1",
            "TDG9-LOC2-STATIONARY-COUNT-STREAM-v1",
        ):
            self.assertIn(domain, production_source)
            self.assertIn(domain, source)
        test_imported = _imported_modules(Path(__file__))
        self.assertFalse(
            any("binder" in name or "scripts" in name for name in test_imported)
        )
        self.assertNotIn("tdg11_rational_complete_c", test_imported)
        self.assertNotIn("tdg11_msel1_reconstruction", test_imported)
        self.assertNotIn("tdg11_msel1_runtime", test_imported)
        self.assertNotIn("assess_rational_complete_c_rows", test_imported)
        self.assertEqual(
            assess_complete_c_independently.__module__,
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization",
        )
        self.assertEqual(
            IndependentCompleteCAdmissionEvidence.__module__,
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization",
        )

    def test_both_tdg9_localizers_are_rerun(self) -> None:
        row = _zero_row()
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum",
            wraps=primary_localize,
        ) as primary_mock, patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            wraps=independent_localize,
        ) as independent_mock:
            evidence = _assess(row)
        self.assertEqual(primary_mock.call_count, 2)
        self.assertEqual(independent_mock.call_count, 2)
        self.assertEqual(
            evidence.d01.primary_evaluator_id,
            "tdg9_loc1_derivative_monotone_bisection_v1",
        )
        self.assertEqual(
            evidence.d01.independent_evaluator_id,
            "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2",
        )
        self.assertTrue(evidence.d01.routes_agree)
        self.assertTrue(evidence.d12.routes_agree)

    def test_hermite_and_half_restriction_match_basis_expansion(self) -> None:
        segment = (Q(1), Q(0), Q(2), Q(0), Q(1))
        cubic = pref1_core._hermite_cubic_on_unit_interval(segment)
        self.assertEqual(cubic, (Q(1), Q(0), Q(3), Q(-2)))
        left = pref1_core._restrict_cubic_to_half(cubic, half=0)
        self.assertEqual(left, (Q(1), Q(0), Q(3, 4), Q(-1, 4)))
        right = pref1_core._restrict_cubic_to_half(cubic, half=1)
        self.assertEqual(right[0], Q(3, 2))
        self.assertEqual(right[0], cubic[0] + cubic[1] / 2 + cubic[2] / 4 + cubic[3] / 8)
        self.assertEqual(right[3], cubic[3] / 8)
        self.assertEqual(
            pref1_core._subtract_cubics(cubic, cubic),
            (Q(0), Q(0), Q(0), Q(0)),
        )
        with self.assertRaisesRegex(ValueError, "half must be zero or one"):
            pref1_core._restrict_cubic_to_half(cubic, half=True)  # type: ignore[arg-type]

    def test_float_bool_and_malformed_inputs_fail_closed(self) -> None:
        row = _zero_row()
        with self.assertRaisesRegex(ValueError, "exactly 4"):
            _assess((row[0], row[1], tuple(row[2][:-1])))
        with self.assertRaisesRegex(ValueError, "exactly 2"):
            _assess((row[0], tuple(row[1][:-1]), row[2]))
        with self.assertRaisesRegex(ValueError, "exactly 3"):
            _assess((row[0], row[1]))
        with self.assertRaisesRegex(ValueError, "exactly 5"):
            _assess((row[0][:-1], row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            _assess((_mutate_segment(row[0], 0, True), row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            _assess((_mutate_segment(row[0], 0, 0.0), row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            _assess((_mutate_segment(row[0], 4, 1.0), row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            _assess((_mutate_segment(row[0], 0, Decimal("0")), row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "must be an integer"):
            assess_complete_c_independently(
                (row,),
                expected_row_count=True,  # type: ignore[arg-type]
                maximum_candidates_D01=8,
                maximum_candidates_D12=8,
            )
        with self.assertRaisesRegex(TypeError, "must be an integer"):
            _assess(row, d01=True)  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "must be an integer"):
            _assess(row, depth=True)  # type: ignore[arg-type]
        integer_row = (
            (0, 0, 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        evidence = _assess(integer_row)
        self.assertEqual(evidence.decision.classification, "exact_zero")

    def test_width_and_shared_boundary_slope_rejection(self) -> None:
        row = _zero_row()
        medium = list(row[1])
        medium[0] = _mutate_segment(medium[0], 4, Q(1, 3))
        with self.assertRaisesRegex(ValueError, "half the outer width"):
            _assess((row[0], tuple(medium), row[2]))
        fine = list(row[2])
        fine[0] = _mutate_segment(fine[0], 4, Q(1, 2))
        with self.assertRaisesRegex(ValueError, "one quarter of the outer width"):
            _assess((row[0], row[1], tuple(fine)))
        with self.assertRaisesRegex(ValueError, "width must be positive"):
            _assess((_mutate_segment(row[0], 4, 0), row[1], row[2]))

        passing = _smoothstep_row(3, 1)
        medium = list(passing[1])
        medium[1] = _mutate_segment(medium[1], 0, Q(1))
        with self.assertRaisesRegex(ValueError, "medium endpoint values"):
            _assess((passing[0], tuple(medium), passing[2]))
        medium = list(passing[1])
        medium[1] = _mutate_segment(medium[1], 1, Q(1))
        with self.assertRaisesRegex(ValueError, "medium endpoint RHS"):
            _assess((passing[0], tuple(medium), passing[2]))
        fine = list(passing[2])
        fine[1] = _mutate_segment(fine[1], 0, Q(1))
        with self.assertRaisesRegex(ValueError, "fine endpoint values"):
            _assess((passing[0], passing[1], tuple(fine)))
        fine = list(passing[2])
        fine[1] = _mutate_segment(fine[1], 1, Q(1))
        with self.assertRaisesRegex(ValueError, "fine endpoint RHS"):
            _assess((passing[0], passing[1], tuple(fine)))
        fine = list(passing[2])
        fine[0] = _mutate_segment(fine[0], 0, Q(1))
        with self.assertRaisesRegex(ValueError, "initial values"):
            _assess((passing[0], passing[1], tuple(fine)))
        medium = list(passing[1])
        medium[0] = _mutate_segment(medium[0], 1, Q(1))
        with self.assertRaisesRegex(ValueError, "initial RHS"):
            _assess((passing[0], tuple(medium), passing[2]))

    def test_row_underflow_and_overflow_are_bounded(self) -> None:
        row = _zero_row()
        with self.assertRaisesRegex(ValueError, "exactly expected_row_count"):
            assess_complete_c_independently(
                (),
                expected_row_count=1,
                maximum_candidates_D01=8,
                maximum_candidates_D12=8,
            )
        with self.assertRaisesRegex(ValueError, "exactly expected_row_count"):
            _assess(row, expected=2)
        with self.assertRaisesRegex(ValueError, "exceed expected_row_count"):
            assess_complete_c_independently(
                (row, row),
                expected_row_count=1,
                maximum_candidates_D01=8,
                maximum_candidates_D12=8,
            )
        evidence = assess_complete_c_independently(
            (row, row),
            expected_row_count=2,
            maximum_candidates_D01=16,
            maximum_candidates_D12=16,
            refinement_depth=32,
        )
        self.assertEqual(evidence.row_count, 2)
        self.assertEqual(evidence.expected_row_count, 2)
        self.assertEqual(evidence.d01.polynomial_count, 4)
        self.assertEqual(evidence.d12.polynomial_count, 8)

    def test_exact_zero_is_admitted_without_inventing_order(self) -> None:
        evidence = _assess(_zero_row())
        self.assertEqual(evidence.decision.classification, "exact_zero")
        self.assertTrue(evidence.decision.admission_passed)
        self.assertIsNone(evidence.decision.order_threshold_passed)
        self.assertFalse(evidence.decision.order_threshold_resolved)
        self.assertEqual(evidence.d01.lower, 0)
        self.assertEqual(evidence.d01.upper, 0)
        self.assertEqual(evidence.d12.lower, 0)
        self.assertEqual(evidence.d12.upper, 0)
        self.assertEqual(evidence.decision.finest_pair_debit, 0)
        self.assertFalse(evidence.sufficient_contraction_pass)
        self.assertFalse(evidence.sufficient_contraction_failure)
        self.assertFalse(evidence.threshold_inconclusive)
        self.assertEqual(evidence.d01.polynomial_count, 2)
        self.assertEqual(evidence.d12.polynomial_count, 4)
        self.assertEqual(evidence.d01.survivor_count, 4)
        self.assertEqual(evidence.d12.survivor_count, 8)
        self.assertEqual(_exact_flags(Q(0), Q(0), Q(0), Q(0)), (False, False, False))
        self.assertEqual(
            evidence.decision,
            classify_tdg6_channel(
                CertifiedMagnitudeInterval(evidence.d01.lower, evidence.d01.upper),
                CertifiedMagnitudeInterval(evidence.d12.lower, evidence.d12.upper),
            ),
        )
        self.assertEqual(evidence.evaluator_id, PREF1_COMPLETE_C_EVALUATOR_ID)

    def test_smoothstep_family_passes_and_fails_contraction(self) -> None:
        passing = _assess(_smoothstep_row(3, 1))
        self.assertEqual(passing.d01.lower, 3)
        self.assertEqual(passing.d01.upper, 3)
        self.assertEqual(passing.d12.lower, 1)
        self.assertEqual(passing.d12.upper, 1)
        self.assertEqual(passing.sufficient_pass_left, 8)
        self.assertEqual(passing.sufficient_pass_right, 9)
        self.assertTrue(passing.sufficient_contraction_pass)
        self.assertFalse(passing.sufficient_contraction_failure)
        self.assertFalse(passing.threshold_inconclusive)
        self.assertEqual(passing.decision.classification, "resolved_order_pass")
        self.assertTrue(passing.decision.admission_passed)
        self.assertTrue(passing.decision.order_threshold_passed)
        self.assertEqual(passing.decision.minimum_observed_order, Q(3, 2))
        self.assertTrue(passing.d01.routes_agree)
        self.assertTrue(passing.d12.routes_agree)
        self.assertEqual(passing.d01.survivor_count, 1)
        self.assertEqual(passing.d12.survivor_count, 1)
        self.assertEqual(
            passing.d01.primary_stationary_count_stream_sha256,
            passing.d01.independent_stationary_count_stream_sha256,
        )
        self.assertTrue(8 * passing.d12.upper**2 <= passing.d01.lower**2)
        self.assertTrue(
            three_halves_order_passes_squared(
                outer_lower_squared=passing.d01.lower**2,
                finest_upper_squared=passing.d12.upper**2,
            )
        )

        failing = _assess(_smoothstep_row(2, 1))
        self.assertEqual(failing.d01.upper, 2)
        self.assertEqual(failing.d12.lower, 1)
        self.assertGreater(8 * failing.d12.lower**2, failing.d01.upper**2)
        self.assertTrue(failing.sufficient_contraction_failure)
        self.assertFalse(failing.sufficient_contraction_pass)
        self.assertFalse(failing.threshold_inconclusive)
        self.assertFalse(failing.decision.admission_passed)
        self.assertEqual(failing.decision.classification, "resolved_order_failure")
        self.assertIs(failing.decision.order_threshold_passed, False)
        self.assertEqual(
            _exact_flags(failing.d01.lower, failing.d01.upper, failing.d12.lower, failing.d12.upper),
            (False, True, False),
        )

    def test_three_halves_equality_uses_squared_predicate(self) -> None:
        self.assertTrue(
            three_halves_order_passes_squared(
                outer_lower_squared=Q(8),
                finest_upper_squared=Q(1),
            )
        )
        self.assertFalse(
            three_halves_order_passes_squared(
                outer_lower_squared=Q(8) - Q(1, 2**52),
                finest_upper_squared=Q(1),
            )
        )
        left, right, passed, failed, inconclusive = pref1_core._contraction_flags(
            outer_lower=Q(8),
            outer_upper=Q(8),
            finest_lower=Q(1),
            finest_upper=Q(1),
        )
        self.assertEqual(left, 8)
        self.assertEqual(right, 64)
        self.assertTrue(passed)
        self.assertFalse(failed)
        self.assertFalse(inconclusive)
        self.assertTrue(left <= right)

    def test_threshold_straddling_is_not_sufficient_failure(self) -> None:
        left, right, passed, failed, inconclusive = pref1_core._contraction_flags(
            outer_lower=Q(2),
            outer_upper=Q(3),
            finest_lower=Q(1),
            finest_upper=Q(1),
        )
        self.assertEqual(left, 8)
        self.assertEqual(right, 4)
        self.assertFalse(passed)
        self.assertFalse(failed)
        self.assertTrue(inconclusive)
        self.assertFalse(8 * Q(1) ** 2 > Q(3) ** 2)

        passing = _assess(_smoothstep_row(3, 1))
        straddle_d01 = replace(passing.d01, lower=Q(2), upper=Q(3))
        decision = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Q(2), Q(3)),
            CertifiedMagnitudeInterval(Q(1), Q(1)),
        )
        self.assertEqual(decision.classification, "resolved_order_failure")
        self.assertFalse(decision.admission_passed)
        straddle = replace(
            passing,
            d01=straddle_d01,
            decision=decision,
            sufficient_pass_left=8 * passing.d12.upper**2,
            sufficient_pass_right=Q(2) ** 2,
            sufficient_contraction_pass=False,
            sufficient_contraction_failure=False,
            threshold_inconclusive=True,
        )
        self.assertTrue(straddle.threshold_inconclusive)
        self.assertFalse(straddle.sufficient_contraction_failure)
        self.assertFalse(straddle.sufficient_contraction_pass)
        self.assertEqual(straddle.decision.classification, "resolved_order_failure")
        self.assertEqual(
            _exact_flags(straddle.d01.lower, straddle.d01.upper, straddle.d12.lower, straddle.d12.upper),
            (False, False, True),
        )
        with self.assertRaisesRegex(ValueError, "sufficient-failure flag"):
            replace(straddle, sufficient_contraction_failure=True)
        with self.assertRaisesRegex(ValueError, "sufficient-pass flag"):
            replace(straddle, sufficient_contraction_pass=True)

    def test_interior_extrema_hidden_by_endpoints_remain_candidates(self) -> None:
        evidence = _assess(_hidden_interior_row())
        self.assertEqual(evidence.d01.lower, Q(4, 27))
        self.assertEqual(evidence.d01.upper, Q(4, 27))
        self.assertGreater(evidence.d01.lower, Q(1, 8))
        self.assertEqual(evidence.d12.lower, 0)
        self.assertEqual(evidence.d12.upper, 0)
        self.assertEqual(evidence.d01.polynomial_count, 2)
        self.assertEqual(evidence.d01.candidate_count, 5)
        self.assertEqual(evidence.d01.survivor_count, 1)
        self.assertEqual(evidence.d01.localization_classification, "unique_maximum")
        expected_stationary = sha256(
            STATIONARY_COUNT_DOMAIN + b"0|0\n1|1\n"
        ).hexdigest()
        self.assertEqual(
            evidence.d01.primary_stationary_count_stream_sha256,
            expected_stationary,
        )
        self.assertEqual(
            evidence.d01.independent_stationary_count_stream_sha256,
            expected_stationary,
        )
        self.assertTrue(evidence.sufficient_contraction_pass)
        self.assertEqual(evidence.decision.classification, "resolved_order_pass")

    def test_maximum_of_row_ordering_does_not_sum_or_use_last_row(self) -> None:
        smaller = _smoothstep_row(3, 1)
        larger = _smoothstep_row(5, 2)
        forward = assess_complete_c_independently(
            (smaller, larger),
            expected_row_count=2,
            maximum_candidates_D01=64,
            maximum_candidates_D12=64,
            refinement_depth=32,
        )
        reversed_rows = assess_complete_c_independently(
            (larger, smaller),
            expected_row_count=2,
            maximum_candidates_D01=64,
            maximum_candidates_D12=64,
            refinement_depth=32,
        )
        self.assertEqual(forward.d01.lower, 5)
        self.assertEqual(forward.d01.upper, 5)
        self.assertEqual(forward.d12.lower, 2)
        self.assertEqual(forward.d12.upper, 2)
        self.assertEqual(forward.d01.lower, reversed_rows.d01.lower)
        self.assertEqual(forward.d12.upper, reversed_rows.d12.upper)
        self.assertNotEqual(forward.d01.lower, 3)
        self.assertNotEqual(forward.d01.lower, 8)
        self.assertNotEqual(forward.d12.lower, 1)
        self.assertEqual(forward.d01.polynomial_count, 4)
        self.assertEqual(forward.d12.polynomial_count, 8)
        self.assertNotEqual(forward.row_stream_sha256, reversed_rows.row_stream_sha256)
        self.assertNotEqual(
            forward.d01.coefficient_stream_sha256,
            reversed_rows.d01.coefficient_stream_sha256,
        )
        self.assertNotEqual(
            forward.d01.survivor_key_stream_sha256,
            reversed_rows.d01.survivor_key_stream_sha256,
        )

    def test_hashes_use_tdg11_domains_and_not_tdg10(self) -> None:
        rows = (_zero_row(),)
        first = _assess(rows[0])
        second = _assess(rows[0])
        zero = (Q(0), Q(0), Q(0), Q(0))
        expected_row = _row_digest(rows)
        expected_d01 = _coefficient_digest((zero, zero))
        expected_d12 = _coefficient_digest((zero, zero, zero, zero))
        expected_combined = _combined_digest(expected_d01, expected_d12)
        for evidence in (first, second):
            self.assertEqual(evidence.row_stream_sha256, expected_row)
            self.assertEqual(evidence.d01.coefficient_stream_sha256, expected_d01)
            self.assertEqual(evidence.d12.coefficient_stream_sha256, expected_d12)
            self.assertEqual(
                evidence.combined_coefficient_stream_sha256, expected_combined
            )
        self.assertEqual(first.row_stream_sha256, second.row_stream_sha256)
        payload = b"0|0/1\n"
        self.assertNotEqual(
            sha256(ROW_HASH_DOMAIN + payload).hexdigest(),
            sha256(TDG10_ROW_HASH_DOMAIN + payload).hexdigest(),
        )
        self.assertNotEqual(ROW_HASH_DOMAIN, TDG10_ROW_HASH_DOMAIN)
        self.assertNotEqual(COEFFICIENT_HASH_DOMAIN, TDG10_COEFFICIENT_HASH_DOMAIN)
        mutated = _assess(_smoothstep_row(0, 1))
        self.assertNotEqual(
            mutated.d12.coefficient_stream_sha256,
            expected_d12,
        )
        empty_survivor = sha256(SURVIVOR_KEY_HASH_DOMAIN).hexdigest()
        self.assertNotEqual(first.d01.survivor_key_stream_sha256, empty_survivor)

    def test_changed_metadata_width_and_continuity_rehash_or_reject(self) -> None:
        baseline = _assess(_zero_row())
        wider = _assess(_zero_row(h=2))
        self.assertNotEqual(baseline.row_stream_sha256, wider.row_stream_sha256)
        self.assertEqual(baseline.d01.lower, wider.d01.lower)
        two_zero = assess_complete_c_independently(
            (_zero_row(), _zero_row()),
            expected_row_count=2,
            maximum_candidates_D01=16,
            maximum_candidates_D12=16,
            refinement_depth=32,
        )
        self.assertNotEqual(baseline.row_stream_sha256, two_zero.row_stream_sha256)
        self.assertNotEqual(
            baseline.d01.survivor_key_stream_sha256,
            two_zero.d01.survivor_key_stream_sha256,
        )
        row = _zero_row()
        with self.assertRaisesRegex(ValueError, "half the outer width"):
            _assess((row[0], (_mutate_segment(row[1][0], 4, Q(1, 3)), row[1][1]), row[2]))
        passing = _smoothstep_row(3, 1)
        broken = list(passing[1])
        broken[0] = _mutate_segment(broken[0], 2, Q(1))
        with self.assertRaisesRegex(ValueError, "medium endpoint values"):
            _assess((passing[0], tuple(broken), passing[2]))

    def test_flags_nonclaims_and_tampered_fields_fail_closed(self) -> None:
        evidence = _assess(_zero_row())
        self.assertFalse(evidence.absolute_tolerance_used)
        self.assertFalse(evidence.physical_signal_used_for_normalization)
        self.assertFalse(evidence.declared_cubic_is_exact_PDE_history)
        self.assertFalse(evidence.production_authority)
        self.assertFalse(evidence.decision.absolute_state_tolerance_used)
        self.assertFalse(evidence.decision.physical_signal_used_for_normalization)
        self.assertNotIn("candidates", evidence.__dataclass_fields__)
        self.assertNotIn("cubics", evidence.__dataclass_fields__)
        self.assertNotIn("candidates", {item.name for item in fields(evidence.d01)})
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, absolute_tolerance_used=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, physical_signal_used_for_normalization=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, declared_cubic_is_exact_PDE_history=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, production_authority=True)
        for field in (
            "absolute_tolerance_used",
            "physical_signal_used_for_normalization",
            "declared_cubic_is_exact_PDE_history",
            "production_authority",
        ):
            for malformed in (None, 0, "", []):
                with self.subTest(field=field, malformed=malformed):
                    with self.assertRaisesRegex(ValueError, "crossed its scope"):
                        replace(evidence, **{field: malformed})
        for malformed in (True, 1.0):
            with self.assertRaisesRegex(ValueError, "crossed its scope"):
                replace(evidence, schema_version=malformed)
        with self.assertRaisesRegex(ValueError, "incompatible with cubic extrema"):
            replace(evidence.d01, polynomial_count=3)
        with self.assertRaisesRegex(ValueError, "row_count must equal"):
            replace(evidence, row_count=2)
        with self.assertRaisesRegex(ValueError, "combined coefficient hash"):
            replace(evidence, combined_coefficient_stream_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "sufficient-failure flag"):
            replace(evidence, sufficient_contraction_failure=True)
        passing = _assess(_smoothstep_row(3, 1))
        with self.assertRaisesRegex(ValueError, "sufficient-pass flag"):
            replace(passing, sufficient_contraction_pass=False)
        with self.assertRaisesRegex(ValueError, "cannot record disagreement"):
            replace(passing.d01, routes_agree=False)
        with self.assertRaisesRegex(ValueError, "unique maximum must retain one"):
            replace(passing.d01, survivor_count=2)
        with self.assertRaisesRegex(ValueError, "stationary-count digests disagree"):
            replace(passing.d01, independent_stationary_count_stream_sha256="0" * 64)

    def test_saved_tdg6_label_is_not_evidence(self) -> None:
        zero = _assess(_zero_row())
        passing = _assess(_smoothstep_row(3, 1))
        with self.assertRaisesRegex(ValueError, "channel decision differs"):
            replace(passing, decision=zero.decision)
        fresh = classify_tdg6_channel(
            CertifiedMagnitudeInterval(passing.d01.lower, passing.d01.upper),
            CertifiedMagnitudeInterval(passing.d12.lower, passing.d12.upper),
        )
        self.assertEqual(passing.decision, fresh)
        self.assertNotEqual(passing.decision, zero.decision)

    def test_resource_ceilings_and_route_disagreement_are_typed_stops(self) -> None:
        row = _zero_row()
        with self.assertRaises(IndependentCompleteCResourceExhausted) as exhausted:
            _assess(row, d01=3, d12=64)
        self.assertEqual(
            exhausted.exception.evidence.reason, "candidate_ceiling_exhausted"
        )
        self.assertEqual(exhausted.exception.evidence.level, "D01")
        self.assertEqual(exhausted.exception.evidence.maximum_candidates_D01, 3)
        self.assertEqual(exhausted.exception.evidence.maximum_candidates_D12, 64)
        self.assertIn("exhausted:candidate_ceiling_exhausted", str(exhausted.exception))
        self.assertNotEqual(
            exhausted.exception.evidence.reason, "resolved_order_failure"
        )

        with self.assertRaises(IndependentCompleteCResourceExhausted) as exhausted_d12:
            _assess(row, d01=64, d12=7)
        self.assertEqual(
            exhausted_d12.exception.evidence.reason, "candidate_ceiling_exhausted"
        )
        self.assertEqual(exhausted_d12.exception.evidence.level, "D12")
        self.assertEqual(exhausted_d12.exception.evidence.maximum_candidates_D12, 7)

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            side_effect=RootIsolationInconclusive(
                "root_location_inconclusive", "synthetic"
            ),
        ):
            with self.assertRaises(IndependentCompleteCResourceExhausted) as isolated:
                _assess(row)
        self.assertEqual(
            isolated.exception.evidence.reason, "root_isolation_inconclusive"
        )
        self.assertIn("root_location_inconclusive", isolated.exception.evidence.detail)

        real = independent_localize(
            [
                IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    {"level": "D01", "row_index": 0, "subinterval": 0},
                ),
                IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    {"level": "D01", "row_index": 0, "subinterval": 1},
                ),
            ],
            maximum_candidates=64,
        )
        shifted = IndependentLocalizationEvidenceV2(
            real.classification,
            real.polynomial_count,
            real.candidate_count,
            real.candidates,
            real.global_absolute_lower + 1,
            real.global_absolute_upper + 1,
            real.maximum_candidates,
            real.initial_refinement_bits,
            real.refinement_schedule,
            real.global_proof_bit_ceiling,
            real.stationary_count_stream_sha256,
        )
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            return_value=shifted,
        ):
            with self.assertRaises(IndependentCompleteCRouteDisagreement) as disagreed:
                _assess(row)
        self.assertEqual(disagreed.exception.evidence.reason, "route_disagreement")
        self.assertIn(
            "global_intervals_overlap=false", disagreed.exception.evidence.detail
        )
        self.assertNotEqual(
            disagreed.exception.evidence.reason, "resolved_order_failure"
        )
        count_shifted = IndependentLocalizationEvidenceV2(
            real.classification,
            real.polynomial_count,
            real.candidate_count + 1,
            real.candidates,
            real.global_absolute_lower,
            real.global_absolute_upper,
            real.maximum_candidates,
            real.initial_refinement_bits,
            real.refinement_schedule,
            real.global_proof_bit_ceiling,
            real.stationary_count_stream_sha256,
        )
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            return_value=count_shifted,
        ):
            with self.assertRaises(IndependentCompleteCRouteDisagreement) as count_case:
                _assess(row)
        self.assertIn("candidate_count_agrees=false", count_case.exception.evidence.detail)
        digest_shifted = IndependentLocalizationEvidenceV2(
            real.classification,
            real.polynomial_count,
            real.candidate_count,
            real.candidates,
            real.global_absolute_lower,
            real.global_absolute_upper,
            real.maximum_candidates,
            real.initial_refinement_bits,
            real.refinement_schedule,
            real.global_proof_bit_ceiling,
            "0" * 64,
        )
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            return_value=digest_shifted,
        ):
            with self.assertRaises(IndependentCompleteCRouteDisagreement) as digest_case:
                _assess(row)
        self.assertIn(
            "stationary_count_digest_agrees=false",
            digest_case.exception.evidence.detail,
        )
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            return_value=replace(real, classification="unrecognized"),
        ):
            with self.assertRaises(IndependentCompleteCRouteDisagreement):
                _assess(row)

    def test_route_bounds_are_reconciled_by_exact_intersection(self) -> None:
        row = _zero_row()
        real = independent_localize(
            [
                IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    (("level", "D01"), ("row_index", 0), ("subinterval", 0)),
                ),
                IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    (("level", "D01"), ("row_index", 0), ("subinterval", 1)),
                ),
            ],
            maximum_candidates=64,
        )

        def looser(cubics, *, maximum_candidates):  # type: ignore[no-untyped-def]
            actual = independent_localize(cubics, maximum_candidates=maximum_candidates)
            if actual.polynomial_count != 2:
                return actual
            return IndependentLocalizationEvidenceV2(
                actual.classification,
                actual.polynomial_count,
                actual.candidate_count,
                actual.candidates,
                actual.global_absolute_lower,
                actual.global_absolute_upper + 1,
                actual.maximum_candidates,
                actual.initial_refinement_bits,
                actual.refinement_schedule,
                actual.global_proof_bit_ceiling,
                actual.stationary_count_stream_sha256,
            )

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization."
            "localize_absolute_maximum_independently_v2",
            side_effect=looser,
        ):
            evidence = _assess(row)
        self.assertEqual(evidence.d01.lower, 0)
        self.assertEqual(evidence.d01.upper, 0)
        self.assertNotEqual(evidence.d01.upper, real.global_absolute_upper + 1)

    def test_falsely_promoted_or_malformed_dataclasses_fail_closed(self) -> None:
        evidence = _assess(_zero_row())
        closed = IndependentCompleteCClosedEvidence(
            reason="candidate_ceiling_exhausted",
            level="D01",
            row_count=1,
            maximum_candidates_D01=3,
            maximum_candidates_D12=64,
            refinement_depth=32,
            detail="primary_localizer",
        )
        with self.assertRaisesRegex(TypeError, "d01 must be"):
            replace(evidence, d01=closed)  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "decision must be"):
            replace(evidence, decision="resolved_order_pass")  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "unknown rational complete-C closed reason"):
            replace(closed, reason="resolved_order_failure")
        with self.assertRaisesRegex(ValueError, "resource evidence must use a resource reason"):
            IndependentCompleteCResourceExhausted(
                replace(closed, reason="route_disagreement")
            )
        with self.assertRaisesRegex(ValueError, "disagreement evidence must use"):
            IndependentCompleteCRouteDisagreement(closed)
        with self.assertRaisesRegex(TypeError, "d01 must be"):
            IndependentCompleteCAdmissionEvidence(
                d01={"lower": 0, "upper": 0},  # type: ignore[arg-type]
                d12=evidence.d12,
                decision=evidence.decision,
                row_count=1,
                expected_row_count=1,
                row_stream_sha256=evidence.row_stream_sha256,
                combined_coefficient_stream_sha256=evidence.combined_coefficient_stream_sha256,
                sufficient_pass_left=0,
                sufficient_pass_right=0,
                sufficient_contraction_pass=False,
                sufficient_contraction_failure=False,
                threshold_inconclusive=False,
                maximum_candidates_D01=64,
                maximum_candidates_D12=64,
                refinement_depth=32,
            )
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, evaluator_id="tdg11_rational_complete_c_dual_route_v1")
        self.assertIsInstance(evidence.d01, IndependentCompleteCDifferenceEvidence)
        self.assertNotIn("candidates", evidence.d01.__dataclass_fields__)


if __name__ == "__main__":
    unittest.main()
