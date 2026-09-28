"""Focused independent controls for the TDG11-IMP1 Bernstein enclosure."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from decimal import Decimal
from fractions import Fraction
from inspect import signature
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    propose_step,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    classify_tdg6_channel,
)
from recursive_horizons.fgc.evolution import tdg11_imp1_enclosure as imp1_core  # noqa: E402
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization import (  # noqa: E402
    assess_complete_c_independently,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_reconstruction import (  # noqa: E402
    reconstruct_channel_independently,
    validate_recorded_family_independently,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_reconstruction import (  # noqa: E402
    ORIGINAL_ARITHMETIC_ID,
    RecordedFamilyError,
    reconstruct_channel,
    validate_recorded_family,
)
from recursive_horizons.fgc.evolution.tdg11_rational_complete_c import (  # noqa: E402
    RationalCompleteCClosedEvidence,
    RationalCompleteCResourceExhausted,
    RationalCompleteCRouteDisagreement,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_enclosure import (  # noqa: E402
    DEFAULT_IMP1_ENCLOSURE_LIMITS,
    IMP1_ENCLOSURE_EVALUATOR_ID,
    IMP1_ENCLOSURE_SCHEMA_VERSION,
    IMP1ChannelAssessment,
    IMP1EnclosureDisagreement,
    IMP1EnclosureFallbackStop,
    IMP1EnclosureLimits,
    IMP1EnclosureResourceStop,
    IMP1FamilyAssessment,
    IMP1PathAssessment,
    assess_imp1_channel,
    assess_imp1_family,
    assess_imp1_rows,
    bound_exact_cubics,
)


Q = Fraction
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_imp1_enclosure.py"
OWNED = slice(1, -4)
BLOCKS = ("u", "p", "q")
FIELDS = ("alpha", "v", "lambda", "R", "phi", "chi")
POINT_COUNT = 9
FIELD_COUNT = 6
OWNED_ROWS = POINT_COUNT - 5
LARGE_BASE = 2.0**52
MACRO_WIDTH = 0.25
EXPONENTIAL_START = 23 / 16
EXPONENTIAL_WIDTHS = (1 / 16, 1 / 32)
EXACT_BUMP_MAX = Q(4, 27)
EXACT_SIX_FIVE_MAX = Q(32, 25)
FORBIDDEN_IMPORT_MARKERS = (
    "tdg11_msel1_pref1",
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
    "tdg10_",
    "hlt16",
    "binder",
    "scripts",
    "proto",
    "campaign",
    "store",
    "numpy",
)


def _imported_modules(path: Path) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
            names.update(alias.name for alias in node.names)
    return names


def _odd_part(denominator: int) -> int:
    odd = denominator
    while odd % 2 == 0:
        odd //= 2
    return odd


def _in_direct_ring(value: Fraction) -> bool:
    return _odd_part(value.denominator) in (1, 3)


def _tiny_limits(**changes: int) -> IMP1EnclosureLimits:
    values = {
        "maximum_owned_rows": 8,
        "maximum_rational_bits": 256,
        "maximum_subdivision_depth": 0,
        "maximum_subdivisions_per_pair": 0,
        "fallback_maximum_candidates_D01": 64,
        "fallback_maximum_candidates_D12": 128,
        "fallback_refinement_depth": 32,
    }
    values.update(changes)
    return IMP1EnclosureLimits(**values)


def _zero_row(*, h: object = 1) -> tuple[object, object, object]:
    width = Q(h)
    outer = (0, 0, 0, 0, width)
    medium = ((0, 0, 0, 0, width / 2), (0, 0, 0, 0, width / 2))
    fine = tuple((0, 0, 0, 0, width / 4) for _ in range(4))
    return (outer, medium, fine)


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


def _bump(t: Fraction) -> Fraction:
    return t**2 - t**3


def _bump_deriv(t: Fraction) -> Fraction:
    return 2 * t - 3 * t**2


def _bump_row(
    outer_scale: object, medium_scale: object, fine_scale: object, *, h: object = 1
) -> tuple[object, object, object]:
    width = Q(h)

    def sample(scale: Fraction, t0: Fraction, t1: Fraction) -> tuple[Fraction, ...]:
        return (
            scale * _bump(t0),
            scale * _bump_deriv(t0) / width,
            scale * _bump(t1),
            scale * _bump_deriv(t1) / width,
            (t1 - t0) * width,
        )

    outer = sample(Q(outer_scale), Q(0), Q(1))
    medium = (sample(Q(medium_scale), Q(0), Q(1, 2)), sample(Q(medium_scale), Q(1, 2), Q(1)))
    fine = tuple(
        sample(Q(fine_scale), Q(index, 4), Q(index + 1, 4)) for index in range(4)
    )
    return (outer, medium, fine)


def _d01_only_row() -> tuple[object, object, object]:
    return _bump_row(0, -1, -1)


def _d12_only_row() -> tuple[object, object, object]:
    return _bump_row(0, 0, -1)


def _straddle_row() -> tuple[object, object, object]:
    # D01 is the unit bump and D12 is (3/8) of the same bump.  Initial
    # Bernstein hulls straddle 8 U12^2 ? L01^2, while the exact maxima
    # 4/27 and 1/18 are a sufficient failure.
    return _bump_row(0, -1, Q(-11, 8))


def _independent(rows: tuple[object, ...], *, row_count: int | None = None):
    return assess_complete_c_independently(
        rows,
        expected_row_count=len(rows) if row_count is None else row_count,
        maximum_candidates_D01=64,
        maximum_candidates_D12=128,
        refinement_depth=32,
    )


def _state(*, u=1.0, p=0.0, q=0.0, points=POINT_COUNT, fields=FIELD_COUNT):
    return EvolutionState(
        np.full((points, fields), u, dtype=np.float64),
        np.full((points, fields), p, dtype=np.float64),
        np.full((points, fields), q, dtype=np.float64),
    )


def _rhs_from_functions(du, dp, dq):
    def rhs(_time, state):
        return EvolutionRHS(du(state), dp(state), dq(state), {})

    return rhs


def _constant_rhs(du=0.0, dp=0.0, dq=0.0):
    return _rhs_from_functions(
        lambda state: np.full_like(state.u, du),
        lambda state: np.full_like(state.p, dp),
        lambda state: np.full_like(state.q, dq),
    )


def _exponential_rhs():
    return _rhs_from_functions(
        lambda state: np.array(state.u, copy=True),
        lambda state: np.array(state.p, copy=True),
        lambda state: np.array(state.q, copy=True),
    )


def _has_factor_three(value: Fraction) -> bool:
    return _odd_part(value.denominator) == 3


def _path(method, state, rhs, start, width, steps):
    proposals = []
    current_state, current_time, step = state, start, width / steps
    for _ in range(steps):
        proposal = propose_step(
            method=method,
            time=current_time,
            step_size=step,
            state=current_state,
            rhs=rhs,
        )
        proposals.append(proposal)
        current_state = proposal.candidate_state
        current_time = proposal.final_time
    return tuple(proposals)


def _paths(method, state, rhs, *, width=MACRO_WIDTH, start=0.0):
    return (
        _path(method, state, rhs, start, width, 1),
        _path(method, state, rhs, start, width, 2),
        _path(method, state, rhs, start, width, 4),
    )


def _family(method, state, rhs, **kwargs):
    paths = _paths(method, state, rhs, **kwargs)
    return paths, validate_recorded_family(
        paths,
        method=method,
        arithmetic_id=ORIGINAL_ARITHMETIC_ID,
        expected_owned_row_count=OWNED_ROWS,
    )


def _closed_evidence(reason: str) -> RationalCompleteCClosedEvidence:
    return RationalCompleteCClosedEvidence(
        reason=reason,
        level="D01",
        row_count=1,
        maximum_candidates_D01=64,
        maximum_candidates_D12=128,
        refinement_depth=32,
        detail="synthetic",
    )


class ImportBoundaryAndSignatureTests(unittest.TestCase):
    def test_producer_stays_inside_the_sealed_instrument_boundary(self) -> None:
        imported = _imported_modules(MODULE_PATH)
        source = MODULE_PATH.read_text()
        self.assertIn("tdg11_msel1_reconstruction", imported)
        self.assertIn("tdg11_rational_complete_c", imported)
        self.assertIn("tdg6_temporal_admission_design", imported)
        self.assertIn("validate_recorded_family", source)
        self.assertIn("reconstruct_channel", source)
        self.assertIn("assess_rational_complete_c_rows", source)
        self.assertIn(IMP1_ENCLOSURE_EVALUATOR_ID, source)
        self.assertNotIn("tdg11_msel1_pref1", source)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )

    def test_public_signatures_match_the_runtime_contract(self) -> None:
        self.assertEqual(
            IMP1_ENCLOSURE_EVALUATOR_ID,
            "tdg11_imp1_exact_bernstein_with_dual_fallback_v1",
        )
        self.assertEqual(IMP1_ENCLOSURE_SCHEMA_VERSION, 1)
        self.assertEqual(
            tuple(signature(IMP1EnclosureLimits).parameters),
            (
                "maximum_owned_rows",
                "maximum_rational_bits",
                "maximum_subdivision_depth",
                "maximum_subdivisions_per_pair",
                "fallback_maximum_candidates_D01",
                "fallback_maximum_candidates_D12",
                "fallback_refinement_depth",
            ),
        )
        self.assertEqual(DEFAULT_IMP1_ENCLOSURE_LIMITS.maximum_owned_rows, 2044)
        self.assertEqual(DEFAULT_IMP1_ENCLOSURE_LIMITS.maximum_rational_bits, 32768)
        self.assertEqual(DEFAULT_IMP1_ENCLOSURE_LIMITS.maximum_subdivision_depth, 12)
        self.assertEqual(
            DEFAULT_IMP1_ENCLOSURE_LIMITS.maximum_subdivisions_per_pair, 32768
        )
        self.assertEqual(
            DEFAULT_IMP1_ENCLOSURE_LIMITS.fallback_maximum_candidates_D01, 16352
        )
        self.assertEqual(
            DEFAULT_IMP1_ENCLOSURE_LIMITS.fallback_maximum_candidates_D12, 32704
        )
        self.assertEqual(DEFAULT_IMP1_ENCLOSURE_LIMITS.fallback_refinement_depth, 160)
        self.assertEqual(
            tuple(signature(assess_imp1_channel).parameters),
            ("family", "channel", "limits"),
        )
        self.assertEqual(
            tuple(signature(assess_imp1_family).parameters),
            ("family", "limits"),
        )
        self.assertEqual(
            tuple(signature(assess_imp1_rows).parameters),
            (
                "rows",
                "expected_row_count",
                "limits",
                "accumulation_bounds",
                "refine",
                "allow_fallback",
            ),
        )
        self.assertEqual(
            tuple(signature(bound_exact_cubics).parameters),
            ("cubics", "limits", "refine"),
        )


class BernsteinCubicControlTests(unittest.TestCase):
    def test_six_t_squared_minus_five_t_cubed_keeps_direct_ring_debit(self) -> None:
        controls = (Q(0), Q(0), Q(2), Q(1))
        initial = bound_exact_cubics((controls,))
        self.assertEqual(initial.lower, Q(81, 64))
        self.assertEqual(initial.upper, 2)
        self.assertLessEqual(initial.lower, EXACT_SIX_FIVE_MAX)
        self.assertLessEqual(EXACT_SIX_FIVE_MAX, initial.upper)
        self.assertEqual(EXACT_SIX_FIVE_MAX.denominator, 25)
        self.assertTrue(_in_direct_ring(initial.lower))
        self.assertTrue(_in_direct_ring(initial.upper))
        self.assertNotEqual(initial.upper, EXACT_SIX_FIVE_MAX)
        refined = bound_exact_cubics((controls,), refine=True)
        self.assertGreaterEqual(refined.lower, initial.lower)
        self.assertLessEqual(refined.upper, initial.upper)
        self.assertLessEqual(refined.lower, EXACT_SIX_FIVE_MAX)
        self.assertLessEqual(EXACT_SIX_FIVE_MAX, refined.upper)
        self.assertTrue(_in_direct_ring(refined.upper))
        self.assertNotEqual(_odd_part(EXACT_SIX_FIVE_MAX.denominator), 1)

    def test_interior_negative_and_repeated_root_bump_enclose_four_over_27(self) -> None:
        bump = (Q(0), Q(0), Q(1, 3), Q(0))
        negative = (Q(0), Q(0), Q(-1, 3), Q(0))
        for controls in (bump, negative):
            enclosure = bound_exact_cubics((controls,))
            self.assertEqual(enclosure.lower, Q(9, 64))
            self.assertEqual(enclosure.upper, Q(1, 3))
            self.assertLessEqual(enclosure.lower, EXACT_BUMP_MAX)
            self.assertLessEqual(EXACT_BUMP_MAX, enclosure.upper)
            self.assertFalse(enclosure.identically_zero)
        zero = bound_exact_cubics(((0, 0, 0, 0),))
        self.assertTrue(zero.identically_zero)
        self.assertEqual(zero.lower, 0)
        self.assertEqual(zero.upper, 0)

    def test_nonzero_cubic_cannot_vanish_all_five_samples(self) -> None:
        # Zero at 0, 1/2, and 1, but the quarter samples remain nonzero.
        interior = bound_exact_cubics(((0, 1, -1, 0),))
        self.assertGreater(interior.lower, 0)
        self.assertFalse(interior.identically_zero)
        mixed = bound_exact_cubics(((0, 0, 0, 0), (0, 0, 1, 0)))
        self.assertEqual(mixed.polynomial_count, 2)
        self.assertGreater(mixed.lower, 0)
        self.assertFalse(mixed.identically_zero)


class RowClassificationTests(unittest.TestCase):
    def test_exact_zero_d01_only_and_d12_only(self) -> None:
        zero = assess_imp1_rows((_zero_row(),), expected_row_count=1)
        self.assertTrue(zero.exact_zero)
        self.assertEqual(zero.decision.classification, "exact_zero")
        self.assertTrue(zero.decision.admission_passed)
        self.assertFalse(zero.sufficient_contraction_pass)
        self.assertFalse(zero.sufficient_contraction_failure)
        self.assertFalse(zero.fallback_called)
        self.assertEqual(zero.fallback_calls, 0)
        self.assertEqual(zero.public_fine_debit, 0)
        self.assertEqual(zero.extra_enclosure_debit, 0)
        independent_zero = _independent((_zero_row(),))
        self.assertEqual(independent_zero.d01.lower, 0)
        self.assertEqual(independent_zero.d01.upper, 0)
        self.assertEqual(independent_zero.d12.lower, 0)
        self.assertEqual(independent_zero.d12.upper, 0)

        d01_only = assess_imp1_rows((_d01_only_row(),), expected_row_count=1)
        independent_d01 = _independent((_d01_only_row(),))
        self.assertEqual(independent_d01.d12.upper, 0)
        self.assertEqual(independent_d01.d01.lower, EXACT_BUMP_MAX)
        self.assertGreater(d01_only.bernstein.d01_lower, 0)
        self.assertEqual(d01_only.bernstein.d12_upper, 0)
        self.assertLessEqual(d01_only.bernstein.d01_lower, independent_d01.d01.lower)
        self.assertLessEqual(independent_d01.d01.upper, d01_only.bernstein.d01_upper)
        self.assertTrue(d01_only.sufficient_contraction_pass)
        self.assertEqual(d01_only.decision.classification, "resolved_order_pass")
        self.assertFalse(d01_only.fallback_called)

        d12_only = assess_imp1_rows((_d12_only_row(),), expected_row_count=1)
        independent_d12 = _independent((_d12_only_row(),))
        self.assertEqual(independent_d12.d01.upper, 0)
        self.assertEqual(independent_d12.d12.lower, EXACT_BUMP_MAX)
        self.assertEqual(d12_only.bernstein.d01_upper, 0)
        self.assertGreater(d12_only.bernstein.d12_lower, 0)
        self.assertLessEqual(d12_only.bernstein.d12_lower, independent_d12.d12.lower)
        self.assertLessEqual(independent_d12.d12.upper, d12_only.bernstein.d12_upper)
        self.assertTrue(d12_only.sufficient_contraction_failure)
        self.assertEqual(d12_only.decision.classification, "order_inconclusive")
        self.assertFalse(d12_only.decision.admission_passed)

    def test_smoothstep_pass_and_failure_match_independent_exact_maxima(self) -> None:
        passing = assess_imp1_rows((_smoothstep_row(3, 1),), expected_row_count=1)
        independent = _independent((_smoothstep_row(3, 1),))
        self.assertEqual(passing.bernstein.d01_lower, independent.d01.lower)
        self.assertEqual(passing.bernstein.d01_upper, independent.d01.upper)
        self.assertEqual(passing.bernstein.d12_lower, independent.d12.lower)
        self.assertEqual(passing.bernstein.d12_upper, independent.d12.upper)
        self.assertEqual(passing.decision.classification, "resolved_order_pass")
        self.assertTrue(passing.sufficient_contraction_pass)
        self.assertFalse(passing.fallback_called)
        self.assertEqual(passing.public_fine_debit, 1)
        self.assertTrue(_in_direct_ring(passing.public_fine_debit))

        failing = assess_imp1_rows((_smoothstep_row(2, 1),), expected_row_count=1)
        self.assertTrue(failing.sufficient_contraction_failure)
        self.assertEqual(failing.decision.classification, "resolved_order_failure")
        self.assertFalse(failing.decision.admission_passed)
        self.assertFalse(failing.fallback_called)


class ThresholdFallbackAndResourceTests(unittest.TestCase):
    def test_threshold_straddle_forces_fallback_then_classifies_from_intersection(
        self,
    ) -> None:
        limits = _tiny_limits()
        path = assess_imp1_rows(
            (_straddle_row(),),
            expected_row_count=1,
            limits=limits,
        )
        independent = _independent((_straddle_row(),))
        exact_d12 = Q(1, 18)
        self.assertTrue(path.bernstein.threshold_inconclusive)
        self.assertTrue(path.fallback_called)
        self.assertEqual(path.fallback_calls, 1)
        self.assertEqual(independent.d01.lower, EXACT_BUMP_MAX)
        self.assertEqual(independent.d01.upper, EXACT_BUMP_MAX)
        self.assertEqual(independent.d12.lower, exact_d12)
        self.assertEqual(independent.d12.upper, exact_d12)
        self.assertLessEqual(path.bernstein.d01_lower, EXACT_BUMP_MAX)
        self.assertLessEqual(EXACT_BUMP_MAX, path.bernstein.d01_upper)
        self.assertLessEqual(path.bernstein.d12_lower, exact_d12)
        self.assertLessEqual(exact_d12, path.bernstein.d12_upper)
        self.assertEqual(path.gate_d01_lower, EXACT_BUMP_MAX)
        self.assertEqual(path.gate_d12_upper, exact_d12)
        self.assertTrue(path.sufficient_contraction_failure)
        self.assertEqual(path.decision.classification, "resolved_order_failure")
        self.assertEqual(path.public_fine_debit, path.bernstein.d12_upper)
        self.assertEqual(path.gate_debit, exact_d12)
        self.assertEqual(
            path.extra_enclosure_debit, path.bernstein.d12_upper - exact_d12
        )
        self.assertGreater(path.extra_enclosure_debit, 0)
        self.assertTrue(_in_direct_ring(path.public_fine_debit))
        self.assertEqual(_odd_part(EXACT_BUMP_MAX.denominator), 27)
        self.assertNotIn("survivor", path.as_mapping())
        self.assertNotIn("stationary", path.as_mapping())

    def test_empty_intersection_and_fallback_errors_are_typed_stops(self) -> None:
        limits = _tiny_limits()

        class _Interval:
            def __init__(self, lower: Fraction, upper: Fraction) -> None:
                self.lower = lower
                self.upper = upper

        class _Evidence:
            d01 = _Interval(Q(10), Q(11))
            d12 = _Interval(Q(10), Q(11))

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_imp1_enclosure."
            "assess_rational_complete_c_rows",
            return_value=_Evidence(),
        ):
            with self.assertRaises(IMP1EnclosureDisagreement) as empty:
                assess_imp1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(empty.exception.level, "D01")

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_imp1_enclosure."
            "assess_rational_complete_c_rows",
            side_effect=RationalCompleteCRouteDisagreement(
                _closed_evidence("route_disagreement")
            ),
        ):
            with self.assertRaises(IMP1EnclosureFallbackStop) as disagreed:
                assess_imp1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(disagreed.exception.reason, "route_disagreement")

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_imp1_enclosure."
            "assess_rational_complete_c_rows",
            side_effect=RationalCompleteCResourceExhausted(
                _closed_evidence("candidate_ceiling_exhausted")
            ),
        ):
            with self.assertRaises(IMP1EnclosureFallbackStop) as exhausted:
                assess_imp1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(exhausted.exception.reason, "candidate_ceiling_exhausted")

    def test_hard_caps_are_resource_stops_and_work_cap_is_not_a_nonpass(self) -> None:
        with self.assertRaises(IMP1EnclosureResourceStop) as rows:
            assess_imp1_rows(
                (_zero_row(), _zero_row()),
                expected_row_count=2,
                limits=_tiny_limits(maximum_owned_rows=1),
            )
        self.assertEqual(rows.exception.reason, "maximum_owned_rows")
        with self.assertRaises(IMP1EnclosureResourceStop) as bits:
            assess_imp1_rows(
                (_zero_row(),),
                expected_row_count=1,
                limits=_tiny_limits(maximum_rational_bits=2),
            )
        self.assertEqual(bits.exception.reason, "maximum_rational_bits")
        unresolved = assess_imp1_rows(
            (_straddle_row(),),
            expected_row_count=1,
            limits=_tiny_limits(),
            refine=True,
            allow_fallback=False,
        )
        self.assertTrue(unresolved.bernstein.threshold_inconclusive)
        self.assertTrue(unresolved.threshold_inconclusive)
        self.assertFalse(unresolved.sufficient_contraction_pass)
        self.assertFalse(unresolved.sufficient_contraction_failure)
        self.assertFalse(unresolved.fallback_called)
        self.assertEqual(unresolved.fallback_calls, 0)


class DirectRingAndMalformedInputTests(unittest.TestCase):
    def test_dyadic_scaling_and_illegal_generic_rationals(self) -> None:
        base = assess_imp1_rows((_d01_only_row(),), expected_row_count=1)
        scaled = assess_imp1_rows((_bump_row(0, Q(-1, 4), Q(-1, 4)),), expected_row_count=1)
        self.assertEqual(scaled.bernstein.d01_lower, base.bernstein.d01_lower / 4)
        self.assertEqual(scaled.bernstein.d01_upper, base.bernstein.d01_upper / 4)
        self.assertTrue(_in_direct_ring(scaled.bernstein.d01_upper))
        illegal_endpoint = (
            (Q(1, 5), 0, 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            assess_imp1_rows((illegal_endpoint,), expected_row_count=1)
        illegal_rhs = (
            (0, Q(1, 5), 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        with self.assertRaisesRegex(ValueError, "dyadic"):
            assess_imp1_rows((illegal_rhs,), expected_row_count=1)
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            bound_exact_cubics(((Q(1, 5), 0, 0, 0),))
        third = assess_imp1_rows(
            (_zero_row(),),
            expected_row_count=1,
            accumulation_bounds=(Q(1, 3), Q(1, 3), Q(1, 3)),
        )
        self.assertEqual(third.public_fine_debit, Q(1, 3))
        self.assertTrue(_in_direct_ring(third.public_fine_debit))
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_imp1_rows(
                (((0.0, 0, 0, 0, 1), _zero_row()[1], _zero_row()[2]),),
                expected_row_count=1,
            )
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_imp1_rows(
                (((Decimal("0"), 0, 0, 0, 1), _zero_row()[1], _zero_row()[2]),),
                expected_row_count=1,
            )
        integer_row = (
            (0, 0, 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        evidence = assess_imp1_rows((integer_row,), expected_row_count=1)
        self.assertTrue(evidence.exact_zero)


class RecordedFamilyTests(unittest.TestCase):
    def test_both_methods_all_eighteen_channels_and_large_base_debit(self) -> None:
        elapsed = {}
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            _, family = _family(method, _state(), _constant_rhs(du=1.0))
            started = time.perf_counter()
            assessment = assess_imp1_family(family)
            elapsed[method] = time.perf_counter() - started
            self.assertIsInstance(assessment, IMP1FamilyAssessment)
            self.assertEqual(assessment.method, method)
            self.assertEqual(assessment.family_sha256, family.family_sha256)
            self.assertEqual(len(assessment.channels), 18)
            self.assertEqual(
                tuple(record.channel for record in assessment.channels),
                TDG6_COMPLETE_STATE_CHANNELS,
            )
            self.assertTrue(assessment.admission_passed)
            self.assertEqual(assessment.evaluator_id, IMP1_ENCLOSURE_EVALUATOR_ID)
            self.assertEqual(len(assessment.public_fine_debits), 18)
            mapping = assessment.as_mapping()
            self.assertEqual(mapping["assessment_sha256"], assessment.assessment_sha256)
            self.assertEqual(len(mapping["channels"]), 18)
            for record, debit in zip(
                assessment.channels, assessment.public_fine_debits, strict=True
            ):
                independent = _independent(
                    reconstruct_channel(family, record.channel).corrected_rows,
                    row_count=OWNED_ROWS,
                )
                self.assertEqual(record.corrected.decision.classification, "exact_zero")
                self.assertEqual(independent.decision.classification, "exact_zero")
                self.assertEqual(record.corrected.public_fine_debit, 0)
                self.assertEqual(debit, 0)
                self.assertFalse(record.corrected.fallback_called)
                self.assertEqual(record.corrected.fallback_calls, 0)
                self.assertFalse(record.raw_unresolved)
                self.assertNotIn("fallback", record.corrected.as_mapping())
                self.assertTrue(_in_direct_ring(record.corrected.public_fine_debit))
        self.assertLess(elapsed[PRIMARY_METHOD], 5.0)
        self.assertLess(elapsed[COMPARATOR_METHOD], 5.0)

        paths, large = _family(
            PRIMARY_METHOD, _state(u=LARGE_BASE), _constant_rhs(du=1.0)
        )
        self.assertEqual(LARGE_BASE + MACRO_WIDTH, LARGE_BASE)
        rebuilt = reconstruct_channel(large, "u:alpha")
        channel = assess_imp1_channel(large, "u:alpha")
        independent_raw = _independent(rebuilt.raw_rows, row_count=OWNED_ROWS)
        independent_corr = _independent(rebuilt.corrected_rows, row_count=OWNED_ROWS)
        self.assertEqual(independent_raw.d01.lower, Q(1, 36))
        self.assertEqual(independent_raw.d01.upper, Q(1, 36))
        self.assertEqual(independent_raw.d12.lower, Q(1, 72))
        self.assertEqual(independent_raw.d12.upper, Q(1, 72))
        self.assertLessEqual(channel.raw.d01_lower, Q(1, 36))
        self.assertLessEqual(Q(1, 36), channel.raw.d01_upper)
        self.assertLessEqual(channel.raw.d12_lower, Q(1, 72))
        self.assertLessEqual(Q(1, 72), channel.raw.d12_upper)
        self.assertEqual(independent_corr.decision.classification, "exact_zero")
        self.assertEqual(channel.corrected.decision.classification, "exact_zero")
        self.assertEqual(channel.corrected.public_fine_debit, Q(1, 4))
        self.assertEqual(channel.corrected.accumulation_bounds[2], Q(1, 4))
        self.assertEqual(channel.corrected.extra_enclosure_debit, 0)
        self.assertFalse(channel.corrected.fallback_called)
        self.assertTrue(_in_direct_ring(channel.corrected.public_fine_debit))
        rest = assess_imp1_channel(large, "p:alpha")
        self.assertEqual(rest.corrected.public_fine_debit, 0)
        self.assertEqual(len(paths[0]) + len(paths[1]) + len(paths[2]), 7)

    def test_missing_and_altered_stage_records_fail_closed(self) -> None:
        paths, family = _family(
            PRIMARY_METHOD, _state(), _constant_rhs(du=1.0)
        )
        truncated = (paths[0], paths[1], paths[2][:-1])
        with self.assertRaises(RecordedFamilyError):
            validate_recorded_family(
                truncated,
                method=PRIMARY_METHOD,
                arithmetic_id=ORIGINAL_ARITHMETIC_ID,
                expected_owned_row_count=OWNED_ROWS,
            )
        array = family.paths[0][0].stages[1].rhs.du
        array.setflags(write=True)
        array[OWNED][0, 0] += 1.0
        array.setflags(write=False)
        with self.assertRaises(RecordedFamilyError):
            assess_imp1_channel(family, "u:alpha")

    def test_exponential_y_prime_y_families_accept_factor_three_corrected_debt(self) -> None:
        self.assertEqual(Q.from_float(EXPONENTIAL_START), Q(23, 16))
        self.assertEqual(Q.from_float(EXPONENTIAL_WIDTHS[0]), Q(1, 16))
        self.assertEqual(Q.from_float(EXPONENTIAL_WIDTHS[1]), Q(1, 32))
        localized = False
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            for width in EXPONENTIAL_WIDTHS:
                _, family = _family(
                    method,
                    _state(u=1.0),
                    _exponential_rhs(),
                    start=EXPONENTIAL_START,
                    width=width,
                )
                assessment = assess_imp1_family(family)
                self.assertEqual(len(assessment.channels), 18)
                self.assertEqual(
                    tuple(record.channel for record in assessment.channels),
                    TDG6_COMPLETE_STATE_CHANNELS,
                )
                independent_family = validate_recorded_family_independently(
                    family.paths,
                    method=family.method,
                    arithmetic_id=ORIGINAL_ARITHMETIC_ID,
                    owned_row_count=OWNED_ROWS,
                )
                independent = reconstruct_channel_independently(
                    independent_family, "u:alpha"
                )
                produced = reconstruct_channel(family, "u:alpha")
                self.assertEqual(
                    independent.accumulation_bounds, produced.accumulation_bounds
                )
                factor_three = False
                for row in independent.corrected_rows:
                    for segment in (row[0], *row[1], *row[2]):
                        self.assertTrue(_in_direct_ring(segment[0]))
                        self.assertTrue(_in_direct_ring(segment[2]))
                        self.assertEqual(_odd_part(segment[1].denominator), 1)
                        self.assertEqual(_odd_part(segment[3].denominator), 1)
                        self.assertEqual(_odd_part(segment[4].denominator), 1)
                        factor_three = factor_three or _has_factor_three(segment[0])
                        factor_three = factor_three or _has_factor_three(segment[2])
                factor_three = factor_three or any(
                    _has_factor_three(item)
                    for item in independent.accumulation_bounds
                )
                self.assertTrue(factor_three)
                channel = assess_imp1_channel(family, "u:alpha")
                self.assertTrue(_in_direct_ring(channel.corrected.public_fine_debit))
                if not localized:
                    independent_loc = _independent(
                        independent.corrected_rows, row_count=OWNED_ROWS
                    )
                    self.assertLessEqual(
                        channel.corrected.bernstein.d01_lower,
                        independent_loc.d01.lower,
                    )
                    self.assertLessEqual(
                        independent_loc.d01.upper,
                        channel.corrected.bernstein.d01_upper,
                    )
                    self.assertLessEqual(
                        channel.corrected.bernstein.d12_lower,
                        independent_loc.d12.lower,
                    )
                    self.assertLessEqual(
                        independent_loc.d12.upper,
                        channel.corrected.bernstein.d12_upper,
                    )
                    localized = True
        self.assertTrue(localized)


class LargeIntegerAndDerivedOverflowTests(unittest.TestCase):
    def test_chunked_decimal_formatting_keeps_global_digit_limit(self) -> None:
        before = sys.get_int_max_str_digits()
        huge = Q(1, 1 << 20000)
        with self.assertRaises(ValueError):
            str(1 << 20000)
        text = imp1_core._canonical_text(huge)
        self.assertTrue(text.startswith("1/"))
        self.assertEqual(text[-1], "6")
        enclosure = bound_exact_cubics(((huge, huge, huge, huge),))
        self.assertEqual(enclosure.lower, huge)
        self.assertEqual(enclosure.upper, huge)
        self.assertEqual(sys.get_int_max_str_digits(), before)
        with self.assertRaises(IMP1EnclosureResourceStop) as oversized:
            bound_exact_cubics(((Q(1, 1 << 40000), 0, 0, 0),))
        self.assertEqual(oversized.exception.reason, "maximum_rational_bits")
        self.assertEqual(sys.get_int_max_str_digits(), before)

    def test_squared_contraction_witness_overflow_is_a_resource_stop(self) -> None:
        before = sys.get_int_max_str_digits()
        with self.assertRaises(IMP1EnclosureResourceStop) as overflow:
            assess_imp1_rows(
                (_smoothstep_row(0, 1 << 20000),),
                expected_row_count=1,
            )
        self.assertEqual(overflow.exception.reason, "maximum_rational_bits")
        self.assertIn("sufficient", overflow.exception.detail)
        self.assertEqual(sys.get_int_max_str_digits(), before)


class ConstructorAndMappingTests(unittest.TestCase):
    def test_changed_classifications_bool_aliases_and_malformed_fields_fail(self) -> None:
        path = assess_imp1_rows((_smoothstep_row(3, 1),), expected_row_count=1)
        with self.assertRaisesRegex(ValueError, "sufficient-pass flag"):
            replace(path, sufficient_contraction_pass=False)
        with self.assertRaises(TypeError):
            replace(path, sufficient_contraction_pass=1)  # type: ignore[arg-type]
        with self.assertRaises(TypeError):
            replace(path, fallback_called=1)  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "public fine debit"):
            replace(path, public_fine_debit=Q(0))
        with self.assertRaisesRegex(ValueError, "extra enclosure debit"):
            replace(path, extra_enclosure_debit=Q(1))
        with self.assertRaisesRegex(ValueError, "lower bound exceeds"):
            replace(path.bernstein, d01_lower=Q(4), d01_upper=Q(1))
        with self.assertRaisesRegex(ValueError, "must be nonnegative"):
            replace(path.bernstein, subdivisions_used=-1)
        with self.assertRaisesRegex(ValueError, "fallback bounds"):
            replace(
                assess_imp1_rows((_zero_row(),), expected_row_count=1),
                fallback_called=True,
                fallback_calls=1,
            )
        _, family = _family(PRIMARY_METHOD, _state(), _constant_rhs(du=1.0))
        assessment = assess_imp1_family(family)
        with self.assertRaises(FrozenInstanceError):
            assessment.admission_passed = False  # type: ignore[misc]
        with self.assertRaisesRegex(ValueError, "TDG6 complete-state order"):
            replace(assessment, channels=tuple(reversed(assessment.channels)))
        with self.assertRaisesRegex(ValueError, "exactly 18"):
            replace(assessment, channels=assessment.channels[:-1])
        with self.assertRaisesRegex(ValueError, "public_fine_debits differ"):
            replace(
                assessment,
                public_fine_debits=(Q(1),) + assessment.public_fine_debits[1:],
            )
        with self.assertRaises(TypeError):
            replace(assessment, admission_passed=1)  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "family assessment hash"):
            replace(assessment, assessment_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(assessment, global_PDE_error_certified=True)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(assessment, production_authority=True)
        channel = assessment.channels[0]
        self.assertIsInstance(channel, IMP1ChannelAssessment)
        with self.assertRaisesRegex(ValueError, "crossed its scope"):
            replace(channel, raw_medium_difference_bound_certified=True)
        with self.assertRaisesRegex(ValueError, "channel assessment hash"):
            replace(channel, channel_sha256="0" * 64)
        mapping = assessment.as_mapping()
        self.assertEqual(mapping["evaluator_id"], IMP1_ENCLOSURE_EVALUATOR_ID)
        self.assertFalse(mapping["global_PDE_error_certified"])
        self.assertEqual(len(mapping["public_fine_debits"]), 18)

    def test_forged_fallback_gate_is_rejected_even_with_coherent_flags(self) -> None:
        path = assess_imp1_rows(
            (_straddle_row(),),
            expected_row_count=1,
            limits=_tiny_limits(),
        )
        self.assertIsInstance(path, IMP1PathAssessment)
        self.assertTrue(path.fallback_called)
        forged_d01 = (path.bernstein.d01_upper, path.bernstein.d01_upper)
        forged_d12 = (path.gate_d12_lower, path.gate_d12_upper)
        left, right, passed, failed, inconclusive, exact_zero = (
            imp1_core._contraction_flags(
                outer_lower=forged_d01[0],
                outer_upper=forged_d01[1],
                finest_lower=forged_d12[0],
                finest_upper=forged_d12[1],
            )
        )
        decision = classify_tdg6_channel(
            CertifiedMagnitudeInterval(forged_d01[0], forged_d01[1]),
            CertifiedMagnitudeInterval(forged_d12[0], forged_d12[1]),
        )
        r_fine = path.accumulation_bounds[2]
        with self.assertRaisesRegex(
            ValueError, "gate does not equal Bernstein/fallback intersection"
        ):
            replace(
                path,
                gate_d01_lower=forged_d01[0],
                gate_d01_upper=forged_d01[1],
                gate_d12_lower=forged_d12[0],
                gate_d12_upper=forged_d12[1],
                decision=decision,
                sufficient_pass_left=left,
                sufficient_pass_right=right,
                sufficient_contraction_pass=passed,
                sufficient_contraction_failure=failed,
                threshold_inconclusive=inconclusive,
                exact_zero=exact_zero,
                public_fine_debit=path.bernstein.d12_upper + r_fine,
                gate_debit=forged_d12[1] + r_fine,
                extra_enclosure_debit=path.bernstein.d12_upper - forged_d12[1],
            )
        with self.assertRaisesRegex(ValueError, "empty Bernstein/fallback intersection"):
            replace(
                path,
                fallback_d01_lower=Q(10),
                fallback_d01_upper=Q(11),
                fallback_d12_lower=Q(10),
                fallback_d12_upper=Q(11),
            )


class IndependentOracleAndPerformanceTests(unittest.TestCase):
    def test_independent_pref1_encloses_and_cheap_family_assessment_is_timed(self) -> None:
        rows = (_smoothstep_row(3, 1), _d01_only_row())
        produced = assess_imp1_rows(rows, expected_row_count=2)
        independent = _independent(rows, row_count=2)
        self.assertLessEqual(produced.bernstein.d01_lower, independent.d01.lower)
        self.assertLessEqual(independent.d01.upper, produced.bernstein.d01_upper)
        self.assertLessEqual(produced.bernstein.d12_lower, independent.d12.lower)
        self.assertLessEqual(independent.d12.upper, produced.bernstein.d12_upper)
        self.assertEqual(
            produced.decision.classification, independent.decision.classification
        )
        self.assertFalse(produced.fallback_called)
        _, family = _family(PRIMARY_METHOD, _state(), _constant_rhs(du=1.0))
        started = time.perf_counter()
        assessment = assess_imp1_family(family)
        elapsed = time.perf_counter() - started
        self.assertTrue(assessment.admission_passed)
        self.assertLess(elapsed, 5.0)
        self.assertGreaterEqual(elapsed, 0.0)


if __name__ == "__main__":
    unittest.main()
