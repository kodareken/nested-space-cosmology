"""C1R1 enclosure equivalence against the sealed IMP1 Bernstein corpus."""

from __future__ import annotations

import ast
from decimal import Decimal
from fractions import Fraction
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
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure import (  # noqa: E402
    C1R1_IMPLEMENTATION_ID,
    C1R1_MATHEMATICAL_OBJECT,
    C1R1_REFERENCE_WIRE_EVALUATOR_ID,
    assess_c1r1_channel,
    assess_c1r1_family,
    assess_c1r1_rows,
    bound_c1r1_cubics,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_enclosure import (  # noqa: E402
    IMP1_ENCLOSURE_EVALUATOR_ID,
    IMP1EnclosureDisagreement,
    IMP1EnclosureFallbackStop,
    IMP1EnclosureLimits,
    IMP1EnclosureResourceStop,
    assess_imp1_channel,
    assess_imp1_family,
    assess_imp1_rows,
    bound_exact_cubics,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_reconstruction import (  # noqa: E402
    ORIGINAL_ARITHMETIC_ID,
    reconstruct_channel,
    validate_recorded_family,
)
from recursive_horizons.fgc.evolution.tdg11_rational_complete_c import (  # noqa: E402
    RationalCompleteCClosedEvidence,
    RationalCompleteCResourceExhausted,
    RationalCompleteCRouteDisagreement,
)


Q = Fraction
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_c1r1_enclosure.py"
OWNED = slice(1, -4)
LARGE_BASE = 2.0**52
MACRO_WIDTH = 0.25
EXPONENTIAL_START = 23 / 16
EXPONENTIAL_WIDTHS = (1 / 16, 1 / 32)
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


def _width_row(h: object) -> tuple[object, object, object]:
    width = Q(h)
    outer = (0, 0, 0, 0, width)
    medium = ((0, 0, 0, 0, width / 2), (0, 0, 0, 0, width / 2))
    fine = tuple((0, 0, 0, 0, width / 4) for _ in range(4))
    return (outer, medium, fine)


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
    return _bump_row(0, -1, Q(-11, 8))


def _state(*, u=1.0, p=0.0, q=0.0, points=9, fields=6):
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


def _inhomogeneous_state(*, points=9):
    u = np.zeros((points, 6), dtype=np.float64)
    for row in range(points):
        u[row, 0] = 1.0 + row / 16.0
        u[row, 4] = row / 32.0
    return EvolutionState(u, np.zeros_like(u), np.zeros_like(u))


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
        expected_owned_row_count=state.shape[0] - 5,
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


class ImportBoundaryTests(unittest.TestCase):
    def test_c1r1_stays_inside_its_sealed_instrument_boundary(self) -> None:
        imported = _imported_modules(MODULE_PATH)
        source = MODULE_PATH.read_text()
        self.assertIn("tdg11_c1r1_ring", imported)
        self.assertIn("tdg11_imp1_enclosure", imported)
        self.assertIn("tdg11_msel1_reconstruction", imported)
        self.assertIn("tdg11_rational_complete_c", imported)
        self.assertIn("assess_rational_complete_c_rows", source)
        self.assertIn("C1R1_IMPLEMENTATION_ID", source)
        self.assertNotEqual(C1R1_IMPLEMENTATION_ID, IMP1_ENCLOSURE_EVALUATOR_ID)
        self.assertEqual(C1R1_REFERENCE_WIRE_EVALUATOR_ID, IMP1_ENCLOSURE_EVALUATOR_ID)
        self.assertEqual(
            C1R1_MATHEMATICAL_OBJECT,
            "exact_accumulation_reconstruction_with_debit",
        )
        self.assertNotIn("tdg11_msel1_pref1", source)
        self.assertNotIn("numpy", imported)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )


class ExactCubicEquivalenceTests(unittest.TestCase):
    def test_independent_cubics_match_imp1_bounds_and_hashes(self) -> None:
        cases = (
            ((Q(0), Q(0), Q(2), Q(1)),),
            ((Q(0), Q(0), Q(1, 3), Q(0)),),
            ((Q(0), Q(0), Q(-1, 3), Q(0)),),
            ((0, 0, 0, 0),),
            ((0, 1, -1, 0),),
            ((0, 0, 0, 0), (0, 0, 1, 0)),
        )
        for cubics in cases:
            for refine in (False, True):
                with self.subTest(cubics=cubics, refine=refine):
                    produced = bound_c1r1_cubics(cubics, refine=refine)
                    reference = bound_exact_cubics(cubics, refine=refine)
                    self.assertEqual(produced.lower, reference.lower)
                    self.assertEqual(produced.upper, reference.upper)
                    self.assertEqual(
                        produced.coefficient_stream_sha256,
                        reference.coefficient_stream_sha256,
                    )
                    self.assertEqual(produced.subdivisions_used, reference.subdivisions_used)
                    self.assertEqual(produced.identically_zero, reference.identically_zero)

    def test_non_ring_cubics_follow_the_imp1_reference_path(self) -> None:
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            bound_c1r1_cubics(((Q(1, 5), 0, 0, 0),))


class RowCorpusEquivalenceTests(unittest.TestCase):
    def test_zero_d01_d12_smoothstep_straddle_and_inhomogeneous_rows(self) -> None:
        cases = (
            ((_zero_row(),), {"expected_row_count": 1}),
            ((_d01_only_row(),), {"expected_row_count": 1}),
            ((_d12_only_row(),), {"expected_row_count": 1}),
            ((_smoothstep_row(3, 1),), {"expected_row_count": 1}),
            ((_smoothstep_row(2, 1),), {"expected_row_count": 1}),
            (
                (_zero_row(), _d01_only_row(), _smoothstep_row(3, 1)),
                {"expected_row_count": 3},
            ),
            (
                (_straddle_row(),),
                {"expected_row_count": 1, "limits": _tiny_limits()},
            ),
            (
                (_straddle_row(),),
                {
                    "expected_row_count": 1,
                    "limits": _tiny_limits(),
                    "refine": True,
                    "allow_fallback": False,
                },
            ),
            (
                (_zero_row(),),
                {
                    "expected_row_count": 1,
                    "accumulation_bounds": (Q(1, 3), Q(1, 3), Q(1, 3)),
                },
            ),
            ((_bump_row(0, Q(-1, 4), Q(-1, 4)),), {"expected_row_count": 1}),
        )
        for rows, kwargs in cases:
            with self.subTest(rows=rows, kwargs=kwargs):
                produced = assess_c1r1_rows(rows, **kwargs)
                reference = assess_imp1_rows(rows, **kwargs)
                self.assertEqual(produced.as_mapping(), reference.as_mapping())
                self.assertEqual(
                    produced.decision.classification, reference.decision.classification
                )
                self.assertEqual(
                    produced.public_fine_debit, reference.public_fine_debit
                )
                self.assertEqual(
                    produced.sufficient_pass_left, reference.sufficient_pass_left
                )
                self.assertEqual(
                    produced.sufficient_pass_right, reference.sufficient_pass_right
                )

    def test_non_ring_malformed_and_resource_failures_match_imp1(self) -> None:
        illegal_endpoint = (
            (Q(1, 5), 0, 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            assess_c1r1_rows((illegal_endpoint,), expected_row_count=1)
        illegal_rhs = (
            (0, Q(1, 5), 0, 0, 1),
            ((0, 0, 0, 0, Q(1, 2)), (0, 0, 0, 0, Q(1, 2))),
            tuple((0, 0, 0, 0, Q(1, 4)) for _ in range(4)),
        )
        with self.assertRaisesRegex(ValueError, "dyadic"):
            assess_c1r1_rows((illegal_rhs,), expected_row_count=1)
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_c1r1_rows(
                (((0.0, 0, 0, 0, 1), _zero_row()[1], _zero_row()[2]),),
                expected_row_count=1,
            )
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_c1r1_rows(
                (((Decimal("0"), 0, 0, 0, 1), _zero_row()[1], _zero_row()[2]),),
                expected_row_count=1,
            )
        with self.assertRaises(IMP1EnclosureResourceStop) as rows:
            assess_c1r1_rows(
                (_zero_row(), _zero_row()),
                expected_row_count=2,
                limits=_tiny_limits(maximum_owned_rows=1),
            )
        self.assertEqual(rows.exception.reason, "maximum_owned_rows")
        with self.assertRaises(IMP1EnclosureResourceStop) as bits:
            assess_c1r1_rows(
                (_zero_row(),),
                expected_row_count=1,
                limits=_tiny_limits(maximum_rational_bits=2),
            )
        self.assertEqual(bits.exception.reason, "maximum_rational_bits")
        with self.assertRaises(IMP1EnclosureResourceStop) as overflow:
            assess_c1r1_rows(
                (_smoothstep_row(0, 1 << 20000),),
                expected_row_count=1,
            )
        self.assertEqual(overflow.exception.reason, "maximum_rational_bits")
        self.assertIn("sufficient", overflow.exception.detail)

    def test_forced_fallback_errors_are_typed_stops(self) -> None:
        limits = _tiny_limits()

        class _Interval:
            def __init__(self, lower: Fraction, upper: Fraction) -> None:
                self.lower = lower
                self.upper = upper

        class _Evidence:
            d01 = _Interval(Q(10), Q(11))
            d12 = _Interval(Q(10), Q(11))

        with patch(
            "recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure."
            "assess_rational_complete_c_rows",
            return_value=_Evidence(),
        ):
            with self.assertRaises(IMP1EnclosureDisagreement) as empty:
                assess_c1r1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(empty.exception.level, "D01")
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure."
            "assess_rational_complete_c_rows",
            side_effect=RationalCompleteCRouteDisagreement(
                _closed_evidence("route_disagreement")
            ),
        ):
            with self.assertRaises(IMP1EnclosureFallbackStop) as disagreed:
                assess_c1r1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(disagreed.exception.reason, "route_disagreement")
        with patch(
            "recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure."
            "assess_rational_complete_c_rows",
            side_effect=RationalCompleteCResourceExhausted(
                _closed_evidence("candidate_ceiling_exhausted")
            ),
        ):
            with self.assertRaises(IMP1EnclosureFallbackStop) as exhausted:
                assess_c1r1_rows(
                    (_straddle_row(),),
                    expected_row_count=1,
                    limits=limits,
                )
        self.assertEqual(exhausted.exception.reason, "candidate_ceiling_exhausted")

    def test_non_ring_prefix_iterator_refuses_like_imp1_and_does_not_skip(self) -> None:
        def _raise_c1r1() -> None:
            assess_c1r1_rows(
                iter((_width_row(Q(1, 3)), _width_row(Q(1)))),
                expected_row_count=1,
            )

        def _raise_imp1() -> None:
            assess_imp1_rows(
                iter((_width_row(Q(1, 3)), _width_row(Q(1)))),
                expected_row_count=1,
            )

        with self.assertRaisesRegex(
            ValueError, r"rows\[0\]\.outer\[4\] is not an exact dyadic rational"
        ) as produced:
            _raise_c1r1()
        with self.assertRaisesRegex(
            ValueError, r"rows\[0\]\.outer\[4\] is not an exact dyadic rational"
        ) as reference:
            _raise_imp1()
        self.assertEqual(type(produced.exception), type(reference.exception))
        self.assertEqual(str(produced.exception), str(reference.exception))
        self.assertNotIn("admission_passed", str(produced.exception))

    def test_late_extra_malformed_generators_and_list_parity(self) -> None:
        valid = (_zero_row(), _d01_only_row())
        from_list = assess_c1r1_rows(valid, expected_row_count=2)
        from_tuple_gen = assess_c1r1_rows((row for row in valid), expected_row_count=2)
        from_imp1 = assess_imp1_rows(valid, expected_row_count=2)
        self.assertEqual(from_list.as_mapping(), from_imp1.as_mapping())
        self.assertEqual(from_tuple_gen.as_mapping(), from_imp1.as_mapping())

        pulls = {"n": 0}

        def _late_non_ring() -> object:
            for item in (_width_row(Q(1)), _width_row(Q(1, 3))):
                pulls["n"] += 1
                yield item

        with self.assertRaisesRegex(
            ValueError, r"rows\[1\]\.outer\[4\] is not an exact dyadic rational"
        ):
            assess_c1r1_rows(_late_non_ring(), expected_row_count=2)
        self.assertEqual(pulls["n"], 2)

        extra_pulls = {"n": 0}

        def _extra() -> object:
            for item in (_zero_row(), _zero_row(), _zero_row()):
                extra_pulls["n"] += 1
                yield item

        with self.assertRaisesRegex(ValueError, "rows exceed expected_row_count"):
            assess_c1r1_rows(_extra(), expected_row_count=1)
        with self.assertRaisesRegex(ValueError, "rows exceed expected_row_count"):
            assess_imp1_rows(
                iter((_zero_row(), _zero_row(), _zero_row())),
                expected_row_count=1,
            )
        self.assertEqual(extra_pulls["n"], 2)

        def _malformed() -> object:
            yield _zero_row()
            yield (
                (0.0, 0, 0, 0, 1),
                _zero_row()[1],
                _zero_row()[2],
            )

        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_c1r1_rows(_malformed(), expected_row_count=2)
        with self.assertRaisesRegex(TypeError, "exact Fraction or built-in int"):
            assess_imp1_rows(
                iter(
                    (
                        _zero_row(),
                        ((0.0, 0, 0, 0, 1), _zero_row()[1], _zero_row()[2]),
                    )
                ),
                expected_row_count=2,
            )

        unused = {"n": 0}

        def _unused() -> object:
            unused["n"] += 1
            yield _zero_row()

        with self.assertRaises(IMP1EnclosureResourceStop) as produced:
            assess_c1r1_rows(
                _unused(),
                expected_row_count=2,
                limits=_tiny_limits(maximum_owned_rows=1),
            )
        with self.assertRaises(IMP1EnclosureResourceStop) as reference:
            assess_imp1_rows(
                iter((_zero_row(), _zero_row())),
                expected_row_count=2,
                limits=_tiny_limits(maximum_owned_rows=1),
            )
        self.assertEqual(produced.exception.reason, reference.exception.reason)
        self.assertEqual(unused["n"], 0)


class DeterministicDifferentialCorpusTests(unittest.TestCase):
    def test_cubics_and_rows_match_imp1_mappings_decisions_hashes_and_stops(
        self,
    ) -> None:
        cubic_cases = (
            ((Q(1, 2), 0, 0, 0),),
            ((Q(1, 3), 0, 0, 0),),
            ((Q(0), Q(0), Q(-1, 3), Q(0)),),
            ((Q(-1, 4), Q(1, 8), Q(-1, 16), Q(1, 32)),),
            ((1, -1, 1, -1),),
            ((0, 0, 0, 0), (Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2))),
            ((Q(0), Q(0), Q(2), Q(1)), (2, 2, 2, 2)),
            ((Q(1, 3), Q(0), Q(-1, 3), Q(0)), (0, 0, Q(1, 3), 0)),
        )
        for cubics in cubic_cases:
            for refine in (False, True):
                with self.subTest(cubics=cubics, refine=refine):
                    produced = bound_c1r1_cubics(cubics, refine=refine)
                    reference = bound_exact_cubics(cubics, refine=refine)
                    self.assertEqual(produced.lower, reference.lower)
                    self.assertEqual(produced.upper, reference.upper)
                    self.assertEqual(
                        produced.coefficient_stream_sha256,
                        reference.coefficient_stream_sha256,
                    )
                    self.assertEqual(
                        produced.subdivisions_used, reference.subdivisions_used
                    )
                    self.assertEqual(
                        produced.identically_zero, reference.identically_zero
                    )

        row_cases = (
            ((_width_row(Q(1, 4)),), {"expected_row_count": 1}),
            ((_d01_only_row(), _d12_only_row()), {"expected_row_count": 2}),
            ((_bump_row(0, -1, Q(-1, 8)),), {"expected_row_count": 1}),
            (
                (_smoothstep_row(3, 1), _smoothstep_row(3, 1)),
                {"expected_row_count": 2},
            ),
            (
                (_straddle_row(),),
                {"expected_row_count": 1, "limits": _tiny_limits()},
            ),
        )
        for rows, kwargs in row_cases:
            with self.subTest(rows=rows, kwargs=kwargs):
                produced = assess_c1r1_rows(rows, **kwargs)
                reference = assess_imp1_rows(rows, **kwargs)
                self.assertEqual(produced.as_mapping(), reference.as_mapping())
                self.assertEqual(
                    produced.decision.classification,
                    reference.decision.classification,
                )
                self.assertEqual(
                    produced.bernstein.row_stream_sha256,
                    reference.bernstein.row_stream_sha256,
                )
                self.assertEqual(
                    produced.bernstein.combined_coefficient_stream_sha256,
                    reference.bernstein.combined_coefficient_stream_sha256,
                )

        tiny = _tiny_limits(maximum_rational_bits=2)
        with self.assertRaises(IMP1EnclosureResourceStop) as produced_bits:
            bound_c1r1_cubics(((Q(1, 8), 0, 0, 0),), limits=tiny)
        with self.assertRaises(IMP1EnclosureResourceStop) as reference_bits:
            bound_exact_cubics(((Q(1, 8), 0, 0, 0),), limits=tiny)
        self.assertEqual(produced_bits.exception.reason, reference_bits.exception.reason)
        self.assertEqual(str(produced_bits.exception), str(reference_bits.exception))
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            bound_c1r1_cubics(((Q(1, 5), 0, 0, 0), (0, 0, 1, 0)))
        with self.assertRaisesRegex(ValueError, "direct Bernstein ring"):
            bound_exact_cubics(((Q(1, 5), 0, 0, 0), (0, 0, 1, 0)))


class RecordedFamilyEquivalenceTests(unittest.TestCase):
    def test_both_methods_nine_thirty_three_one_twenty_nine_and_inhomogeneous(
        self,
    ) -> None:
        elapsed: dict[tuple[str, int, str], float] = {}
        sizes = (9, 33, 129)
        for points in sizes:
            for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
                _, family = _family(
                    method, _state(points=points), _constant_rhs(du=1.0)
                )
                started = time.perf_counter()
                produced = assess_c1r1_family(family)
                elapsed[(method, points, "c1r1")] = time.perf_counter() - started
                started = time.perf_counter()
                reference = assess_imp1_family(family)
                elapsed[(method, points, "imp1")] = time.perf_counter() - started
                self.assertEqual(produced.as_mapping(), reference.as_mapping())
                self.assertEqual(produced.assessment_sha256, reference.assessment_sha256)
                self.assertEqual(produced.evaluator_id, IMP1_ENCLOSURE_EVALUATOR_ID)
                self.assertEqual(len(produced.channels), 18)
                self.assertEqual(
                    tuple(record.channel for record in produced.channels),
                    TDG6_COMPLETE_STATE_CHANNELS,
                )
                self.assertTrue(produced.admission_passed)

        _, inhomogeneous = _family(
            PRIMARY_METHOD, _inhomogeneous_state(), _exponential_rhs()
        )
        produced = assess_c1r1_family(inhomogeneous)
        reference = assess_imp1_family(inhomogeneous)
        self.assertEqual(produced.as_mapping(), reference.as_mapping())
        hashes = {
            record.corrected.bernstein.row_stream_sha256
            for record in produced.channels
        }
        self.assertGreater(len(hashes), 1)

        paths, large = _family(
            PRIMARY_METHOD, _state(u=LARGE_BASE), _constant_rhs(du=1.0)
        )
        self.assertEqual(LARGE_BASE + MACRO_WIDTH, LARGE_BASE)
        produced = assess_c1r1_channel(large, "u:alpha")
        reference = assess_imp1_channel(large, "u:alpha")
        self.assertEqual(produced.as_mapping(), reference.as_mapping())
        self.assertEqual(produced.corrected.public_fine_debit, Q(1, 4))
        self.assertEqual(len(paths[0]) + len(paths[1]) + len(paths[2]), 7)

        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            for width in EXPONENTIAL_WIDTHS:
                _, family = _family(
                    method,
                    _state(u=1.0),
                    _exponential_rhs(),
                    start=EXPONENTIAL_START,
                    width=width,
                )
                produced = assess_c1r1_family(family)
                reference = assess_imp1_family(family)
                self.assertEqual(produced.as_mapping(), reference.as_mapping())
                rebuilt = reconstruct_channel(family, "u:alpha")
                channel = assess_c1r1_channel(family, "u:alpha")
                self.assertEqual(
                    channel.corrected.accumulation_bounds, rebuilt.accumulation_bounds
                )

        # Small-profile evidence only; a speedup is not assumed or required.
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            for points in sizes:
                self.assertGreaterEqual(elapsed[(method, points, "c1r1")], 0.0)
                self.assertGreaterEqual(elapsed[(method, points, "imp1")], 0.0)
        self.assertLess(elapsed[(PRIMARY_METHOD, 9, "c1r1")], 5.0)
        self.assertLess(elapsed[(COMPARATOR_METHOD, 9, "c1r1")], 5.0)


if __name__ == "__main__":
    unittest.main()
