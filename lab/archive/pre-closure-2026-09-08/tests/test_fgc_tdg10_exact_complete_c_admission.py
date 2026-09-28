"""Focused controls for the pure TDG10 exact complete-C admission core."""

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
    three_halves_order_passes_squared,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (  # noqa: E402
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    assess_exact_complete_c_rows,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg9_local_extrema_independent_v2 as independent_v2,
    tdg10_exact_complete_c_admission as tdg10_core,
)


Q = Fraction
Cubic = tuple[Fraction, Fraction, Fraction, Fraction]
MODULE_PATH = (
    ROOT
    / "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_admission.py"
)
ROW_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-ROW-STREAM-v1\n"
COEFFICIENT_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-CUBIC-STREAM-v1\n"
COMBINED_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-COMBINED-v1\n"
SURVIVOR_KEY_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-SURVIVOR-KEY-STREAM-v1\n"
STATIONARY_COUNT_DOMAIN = b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n"
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
    "binder",
    "scripts",
    "tdg9_exact_temporal_arithmetic_independent",
)


def _evaluate(coefficients: Cubic, point: Fraction) -> Fraction:
    a0, a1, a2, a3 = coefficients
    return a0 + point * (a1 + point * (a2 + point * a3))


def _derivative(coefficients: Cubic, point: Fraction) -> Fraction:
    return coefficients[1] + point * (2 * coefficients[2] + point * 3 * coefficients[3])


def _segment(coefficients: Cubic, width: Fraction) -> tuple[float, ...]:
    values = (
        _evaluate(coefficients, Q(0)),
        _derivative(coefficients, Q(0)) / width,
        _evaluate(coefficients, Q(1)),
        _derivative(coefficients, Q(1)) / width,
        width,
    )
    answer = tuple(float(item) for item in values)
    if any(Q(*item.as_integer_ratio()) != exact for item, exact in zip(answer, values)):
        raise AssertionError("test fixture is not exactly binary64 representable")
    return answer


def _restrict(coefficients: Cubic, half: int) -> Cubic:
    a0, a1, a2, a3 = coefficients
    if half == 0:
        return a0, a1 / 2, a2 / 4, a3 / 8
    return (
        a0 + a1 / 2 + a2 / 4 + a3 / 8,
        a1 / 2 + a2 / 2 + 3 * a3 / 8,
        a2 / 4 + 3 * a3 / 8,
        a3 / 8,
    )


def _subtract(left: Cubic, right: Cubic) -> Cubic:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _scaled(coefficients: Cubic, factor: Fraction) -> Cubic:
    return tuple(factor * item for item in coefficients)  # type: ignore[return-value]


def _row(d01: Cubic, d12: Cubic) -> tuple[object, object, object]:
    zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
    medium_polynomials = tuple(_scaled(d01, Q(-1)) for _ in range(2))
    fine_polynomials = tuple(
        _subtract(
            _restrict(medium_polynomials[medium_index], half),
            d12,
        )
        for medium_index in range(2)
        for half in (0, 1)
    )
    return (
        _segment(zero, Q(1)),
        tuple(_segment(item, Q(1, 2)) for item in medium_polynomials),
        tuple(_segment(item, Q(1, 4)) for item in fine_polynomials),
    )


def _assess(row: object, *, d01: int = 64, d12: int = 64, depth: int = 32):
    return assess_exact_complete_c_rows(
        (row,),
        maximum_candidates_D01=d01,
        maximum_candidates_D12=d12,
        refinement_depth=depth,
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


def _row_digest(rows: tuple[object, ...]) -> str:
    hasher = sha256(ROW_HASH_DOMAIN)
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        values = [value for segment in (outer, *medium, *fine) for value in segment]
        hasher.update(
            (f"{ordinal}|" + "|".join(item.hex() for item in values) + "\n").encode(
                "ascii"
            )
        )
    return hasher.hexdigest()


def _coefficient_digest(coefficients: tuple[Cubic, ...]) -> str:
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
    return sha256(COMBINED_HASH_DOMAIN + f"{outer}\n{finest}\n".encode("ascii")).hexdigest()


class TDG10ExactCompleteCAdmissionTests(unittest.TestCase):
    def test_module_has_no_historical_binder_runtime_or_script_imports(self) -> None:
        imported = _imported_modules(MODULE_PATH)
        source = MODULE_PATH.read_text()
        self.assertIn("tdg6_temporal_admission_design", imported)
        self.assertIn("classify_tdg6_channel", imported)
        self.assertIn("exact_hermite_coefficients", imported)
        self.assertIn("restrict_exact_cubic_to_half", imported)
        self.assertIn("subtract_exact_cubics", imported)
        self.assertIn("localize_absolute_maximum", imported)
        self.assertIn("localize_absolute_maximum_independently_v2", imported)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )
        self.assertNotIn("tdg6_temporal_admission_runtime", imported)
        self.assertNotIn("tdg9_exact_temporal_arithmetic_independent", imported)
        self.assertNotIn("import scripts", source)
        self.assertNotIn("from scripts", source)
        test_imported = _imported_modules(Path(__file__))
        self.assertFalse(
            any("binder" in name or "scripts" in name for name in test_imported)
        )
        self.assertNotIn("tdg6_temporal_admission_runtime", test_imported)

    def test_malformed_nonfinite_and_nonbuiltin_inputs_fail_closed(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        broken_fine = (row[0], row[1], tuple(row[2][:-1]))
        with self.assertRaisesRegex(ValueError, "exactly 4"):
            _assess(broken_fine)
        with self.assertRaisesRegex(ValueError, "exactly 2"):
            _assess((row[0], tuple(row[1][:-1]), row[2]))
        with self.assertRaisesRegex(ValueError, "exactly 3"):
            _assess((row[0], row[1]))
        with self.assertRaisesRegex(ValueError, "exactly 5"):
            _assess((row[0][:-1], row[1], row[2]))
        with self.assertRaisesRegex(ValueError, "at least one"):
            assess_exact_complete_c_rows(
                (),
                maximum_candidates_D01=8,
                maximum_candidates_D12=8,
            )
        with self.assertRaisesRegex(TypeError, "built-in binary64"):
            poisoned = (True, 0.0, 0.0, 0.0, 1.0)
            _assess((poisoned, row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "built-in binary64"):
            poisoned = (1, 0.0, 0.0, 0.0, 1.0)
            _assess((poisoned, row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "built-in binary64"):
            poisoned = (Decimal("0"), 0.0, 0.0, 0.0, 1.0)
            _assess((poisoned, row[1], row[2]))
        with self.assertRaisesRegex(ValueError, "must be finite"):
            poisoned = (float("nan"), 0.0, 0.0, 0.0, 1.0)
            _assess((poisoned, row[1], row[2]))
        with self.assertRaisesRegex(ValueError, "must be finite"):
            poisoned = (float("inf"), 0.0, 0.0, 0.0, 1.0)
            _assess((poisoned, row[1], row[2]))
        with self.assertRaisesRegex(TypeError, "must be an integer"):
            assess_exact_complete_c_rows(
                (row,),
                maximum_candidates_D01=True,  # type: ignore[arg-type]
                maximum_candidates_D12=8,
            )

    def test_width_mismatch_fails_before_classification(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        medium = list(row[1])
        left = list(medium[0])
        left[4] = 0.25
        medium[0] = tuple(left)
        with self.assertRaisesRegex(ValueError, "half the outer width"):
            _assess((row[0], tuple(medium), row[2]))
        fine = list(row[2])
        first = list(fine[0])
        first[4] = 0.5
        fine[0] = tuple(first)
        with self.assertRaisesRegex(ValueError, "one quarter of the outer width"):
            _assess((row[0], row[1], tuple(fine)))

    def test_exact_zero_is_admitted_without_inventing_order(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        evidence = _assess(_row(zero, zero))
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

    def test_exact_three_halves_equality_is_admitted(self) -> None:
        self.assertTrue(
            three_halves_order_passes_squared(
                outer_lower_squared=Q(8),
                finest_upper_squared=Q(1),
            )
        )
        # Rational cubics cannot place positive D01 and D12 exactly on
        # sqrt(8):1.  The SSPRK3 freeze control (8, 1) is an exact
        # p >= 3/2 admission, and the squared helper equality is the
        # p = 3/2 boundary itself.
        constant_eight: Cubic = (Q(8), Q(0), Q(0), Q(0))
        constant_one: Cubic = (Q(1), Q(0), Q(0), Q(0))
        evidence = _assess(_row(constant_eight, constant_one))
        self.assertEqual(evidence.d01.lower, 8)
        self.assertEqual(evidence.d01.upper, 8)
        self.assertEqual(evidence.d12.lower, 1)
        self.assertEqual(evidence.d12.upper, 1)
        self.assertEqual(evidence.sufficient_pass_left, 8)
        self.assertEqual(evidence.sufficient_pass_right, 64)
        self.assertTrue(evidence.sufficient_contraction_pass)
        self.assertEqual(evidence.decision.classification, "resolved_order_pass")
        self.assertTrue(evidence.decision.admission_passed)
        self.assertTrue(evidence.decision.order_threshold_passed)
        self.assertEqual(evidence.decision.minimum_observed_order, Q(3, 2))

    def test_one_rational_perturbation_below_threshold_is_not_admitted(self) -> None:
        self.assertFalse(
            three_halves_order_passes_squared(
                outer_lower_squared=Q(8) - Q(1, 2**52),
                finest_upper_squared=Q(1),
            )
        )
        constant_eight: Cubic = (Q(8), Q(0), Q(0), Q(0))
        constant_three: Cubic = (Q(3), Q(0), Q(0), Q(0))
        evidence = _assess(_row(constant_eight, constant_three))
        self.assertEqual(evidence.d01.lower, 8)
        self.assertEqual(evidence.d12.upper, 3)
        self.assertEqual(evidence.sufficient_pass_left, 72)
        self.assertEqual(evidence.sufficient_pass_right, 64)
        self.assertFalse(evidence.decision.admission_passed)
        self.assertEqual(evidence.decision.classification, "resolved_order_failure")
        self.assertIs(evidence.decision.order_threshold_passed, False)
        self.assertTrue(evidence.decision.temporal_retry_permitted)

    def test_algebraic_threshold_inconclusive_does_not_become_a_pass(self) -> None:
        # max |4t^3-6t| on [0,1] is 2*sqrt(2), so the true squared comparison
        # is 8 == 8.  Rational isolating intervals must not round that into
        # either a sufficient pass or a sufficient failure.
        boundary: Cubic = (Q(0), Q(-6), Q(0), Q(4))
        constant: Cubic = (Q(1), Q(0), Q(0), Q(0))
        evidence = _assess(_row(boundary, constant), depth=8)
        self.assertTrue(evidence.threshold_inconclusive)
        self.assertFalse(evidence.sufficient_contraction_pass)
        self.assertFalse(evidence.sufficient_contraction_failure)
        self.assertGreater(evidence.d01.lower, 0)
        self.assertEqual(evidence.d12.lower, 1)
        self.assertEqual(evidence.d12.upper, 1)
        self.assertGreater(evidence.sufficient_pass_left, evidence.sufficient_pass_right)
        self.assertLessEqual(
            8 * evidence.d12.lower**2,
            evidence.d01.upper**2,
        )
        self.assertFalse(evidence.decision.admission_passed)
        self.assertEqual(evidence.decision.classification, "resolved_order_failure")
        self.assertTrue(evidence.d01.routes_agree)
        self.assertTrue(evidence.d12.routes_agree)

    def test_strict_sufficient_failure_is_certified(self) -> None:
        constant_two: Cubic = (Q(2), Q(0), Q(0), Q(0))
        constant_one: Cubic = (Q(1), Q(0), Q(0), Q(0))
        evidence = _assess(_row(constant_two, constant_one))
        self.assertEqual(evidence.d01.upper, 2)
        self.assertEqual(evidence.d12.lower, 1)
        self.assertGreater(
            8 * evidence.d12.lower**2,
            evidence.d01.upper**2,
        )
        self.assertTrue(evidence.sufficient_contraction_failure)
        self.assertFalse(evidence.sufficient_contraction_pass)
        self.assertFalse(evidence.threshold_inconclusive)
        self.assertFalse(evidence.decision.admission_passed)
        self.assertEqual(evidence.decision.classification, "resolved_order_failure")

    def test_routes_agree_on_counts_digests_survivors_and_global_overlap(self) -> None:
        linear: Cubic = (Q(0), Q(1), Q(0), Q(0))
        evidence = _assess(_row(_scaled(linear, Q(16)), linear))
        self.assertTrue(evidence.d01.routes_agree)
        self.assertTrue(evidence.d12.routes_agree)
        self.assertEqual(
            evidence.d01.primary_stationary_count_stream_sha256,
            evidence.d01.independent_stationary_count_stream_sha256,
        )
        self.assertEqual(
            evidence.d12.primary_stationary_count_stream_sha256,
            evidence.d12.independent_stationary_count_stream_sha256,
        )
        self.assertEqual(evidence.d01.polynomial_count, 2)
        self.assertEqual(evidence.d12.polynomial_count, 4)
        self.assertGreaterEqual(evidence.d01.candidate_count, evidence.d01.co_maximizer_count)
        self.assertGreaterEqual(evidence.d12.candidate_count, evidence.d12.co_maximizer_count)
        expected_stationary = sha256(STATIONARY_COUNT_DOMAIN + b"0|0\n1|0\n").hexdigest()
        self.assertEqual(
            evidence.d01.primary_stationary_count_stream_sha256,
            expected_stationary,
        )
        expected_d12_stationary = sha256(
            STATIONARY_COUNT_DOMAIN + b"0|0\n1|0\n2|0\n3|0\n"
        ).hexdigest()
        self.assertEqual(
            evidence.d12.independent_stationary_count_stream_sha256,
            expected_d12_stationary,
        )
        self.assertEqual(evidence.decision.classification, "resolved_order_pass")
        field_names = {item.name for item in fields(evidence.d01)}
        self.assertNotIn("candidates", field_names)
        self.assertNotIn("primary_candidates", field_names)

    def test_deterministic_coefficient_and_row_hashes_use_tdg10_domains(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        rows = (_row(zero, zero),)
        first = _assess(rows[0])
        second = _assess(rows[0])
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
        self.assertEqual(
            first.d01.survivor_key_stream_sha256,
            second.d01.survivor_key_stream_sha256,
        )

        constant_one: Cubic = (Q(1), Q(0), Q(0), Q(0))
        mutated = _assess(_row(zero, constant_one))
        self.assertNotEqual(
            mutated.d12.coefficient_stream_sha256,
            expected_d12,
        )
        self.assertNotEqual(
            mutated.combined_coefficient_stream_sha256,
            expected_combined,
        )
        hasher = sha256(SURVIVOR_KEY_HASH_DOMAIN)
        self.assertTrue(first.d01.survivor_key_stream_sha256)
        self.assertNotEqual(first.d01.survivor_key_stream_sha256, hasher.hexdigest())

    def test_flags_and_nonclaims_are_frozen(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        evidence = _assess(_row(zero, zero))
        self.assertTrue(evidence.radius_free_complete_C)
        self.assertFalse(evidence.absolute_tolerance_used)
        self.assertFalse(evidence.physical_signal_used_for_normalization)
        self.assertFalse(evidence.declared_cubic_is_exact_PDE_history)
        self.assertFalse(evidence.pde_trajectory_order_certified)
        self.assertFalse(evidence.decision.absolute_state_tolerance_used)
        self.assertFalse(evidence.decision.physical_signal_used_for_normalization)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, radius_free_complete_C=False)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, absolute_tolerance_used=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, physical_signal_used_for_normalization=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, declared_cubic_is_exact_PDE_history=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(evidence, pde_trajectory_order_certified=True)
        self.assertNotIn("candidates", evidence.__dataclass_fields__)

    def test_candidate_ceiling_and_route_disagreement_fail_closed(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        with self.assertRaises(ExactCompleteCResourceExhausted) as exhausted:
            _assess(row, d01=3, d12=64)
        self.assertEqual(
            exhausted.exception.evidence.reason, "candidate_ceiling_exhausted"
        )
        self.assertEqual(exhausted.exception.evidence.level, "D01")

        real = independent_v2.localize_absolute_maximum_independently_v2(
            [
                independent_v2.IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    {"level": "D01", "row_index": 0, "subinterval": 0},
                ),
                independent_v2.IndependentLocalCubicV2(
                    (Q(0), Q(0), Q(0), Q(0)),
                    {"level": "D01", "row_index": 0, "subinterval": 1},
                ),
            ],
            maximum_candidates=64,
        )
        shifted = independent_v2.IndependentLocalizationEvidenceV2(
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
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission."
            "localize_absolute_maximum_independently_v2",
            return_value=shifted,
        ):
            with self.assertRaises(ExactCompleteCRouteDisagreement) as disagreed:
                _assess(row)
        self.assertEqual(disagreed.exception.evidence.reason, "route_disagreement")
        self.assertIn(
            "global_intervals_overlap=false", disagreed.exception.evidence.detail
        )

    def test_tolerance_or_boundary_clipping_flags_fail_closed(self) -> None:
        constant: Cubic = (Q(1), Q(0), Q(0), Q(0))
        row = _row(constant, constant)
        real_primary = tdg10_core.localize_absolute_maximum
        real_independent = tdg10_core.localize_absolute_maximum_independently_v2

        def tolerant_primary(cubics, **kwargs):
            return replace(real_primary(cubics, **kwargs), tolerance_used=True)

        with patch.object(
            tdg10_core, "localize_absolute_maximum", side_effect=tolerant_primary
        ):
            with self.assertRaisesRegex(ValueError, "absolute tolerance"):
                _assess(row)

        def clipped_independent(cubics, **kwargs):
            return replace(
                real_independent(cubics, **kwargs), boundary_clipping_used=True
            )

        with patch.object(
            tdg10_core,
            "localize_absolute_maximum_independently_v2",
            side_effect=clipped_independent,
        ):
            with self.assertRaisesRegex(ValueError, "boundary clipping"):
                _assess(row)


if __name__ == "__main__":
    unittest.main()
