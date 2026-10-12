"""Focused controls for the pure TDG10 exact complete-C runtime adapter."""

from __future__ import annotations

import ast
from copy import copy
from dataclasses import fields, replace
from fractions import Fraction
from hashlib import sha256
import inspect
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6PreparedGR0Compositor,
    TDG6TemporalLedger,
    prepare_tdg6_gr0_compositor,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (  # noqa: E402
    ExactCompleteCClosedEvidence,
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    assess_exact_complete_c_rows,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime import (  # noqa: E402
    TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01,
    TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12,
    TDG10_PRODUCTION_OWNED_ROW_COUNT,
    TDG10_PRODUCTION_REFINEMENT_DEPTH,
    TDG10ExactCompleteCRuntimeClosed,
    TDG10ExactCompleteCRuntimeEvidence,
    assess_tdg10_exact_complete_c,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg10_exact_complete_c_runtime as tdg10_runtime,
)


Q = Fraction
Cubic = tuple[Fraction, Fraction, Fraction, Fraction]
MODULE_PATH = (
    ROOT / "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_runtime.py"
)
POINT_COUNT = 9
OWNED_ROW_COUNT = 4
OWNED_INNER_ROW = 1
PROJECTOR_OWNED_OUTER_ROWS = 4
CEILING_D01 = 64
CEILING_D12 = 64
DEPTH = 8
_FIELD_NAMES = ("alpha", "v", "lambda", "R", "phi", "chi")
_BLOCK_NAMES = ("u", "p", "q")
FORBIDDEN_IMPORT_MARKERS = (
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
)


def _marker(block: int, field: int, row: int) -> float:
    return float(1000 * block + 100 * field + row + 1)


def _state() -> EvolutionState:
    arrays = []
    for block in range(3):
        values = np.zeros((POINT_COUNT, 6), dtype=np.float64)
        for field in range(6):
            for row in range(POINT_COUNT):
                values[row, field] = _marker(block, field, row)
        arrays.append(values)
    return EvolutionState(*arrays)


def _rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
    zero = np.zeros_like(state.u)
    return EvolutionRHS(
        zero,
        zero,
        zero,
        {
            "source_residual_infinity": 0.0,
            "source_refinement_iterations": 0,
            "source_residual_decreased_monotonically": True,
            "kinetic_condition_infinity": 1.0,
            "coordinate_speed_upper": 1.0,
            "minimum_lapse": 1.0,
            "minimum_radial_metric": 1.0,
            "minimum_areal_radius_away_from_center": 0.125,
        },
    )


def _poisoned_rhs(rhs: EvolutionRHS) -> EvolutionRHS:
    du = np.asarray(rhs.du, dtype=np.float64).copy()
    du[0, 0] = np.nextafter(du[0, 0], np.inf)
    return EvolutionRHS(du, rhs.dp, rhs.dq, rhs.diagnostics)


class _Tracers:
    def __init__(self) -> None:
        self.labels = np.asarray((0.0, 1.0, 2.0), dtype=np.float64)
        self.positions = np.asarray((0.25, 0.50, 0.75), dtype=np.float64)
        self.proper_times = np.zeros(3, dtype=np.float64)
        self.event_proper_times: list[np.ndarray] = []
        self.event_fields: list[np.ndarray] = []
        self.cutoff = 1.0
        self.outer_radius = 128.0

    def preview_advance(self, *, old_state, new_state, coordinates, step_size):
        return (
            np.asarray(self.positions, dtype=np.float64).copy(),
            np.asarray(self.proper_times, dtype=np.float64).copy() + float(step_size),
        )

    def commit_advance(self, positions, proper_times) -> None:
        self.positions = np.asarray(positions, dtype=np.float64).copy()
        self.proper_times = np.asarray(proper_times, dtype=np.float64).copy()


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )


def _prepare(method: str = COMPARATOR_METHOD):
    state = _state()
    transaction = _transaction()
    tracers = _Tracers()
    ledger = TDG6TemporalLedger.zero(initial_time=0.0)
    prepared = prepare_tdg6_gr0_compositor(
        method=method,
        time=0.0,
        step_size=0.125,
        state=state,
        rhs=_rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=np.linspace(0.0, 1.0, POINT_COUNT),
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger


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


def _compositor_digest(prepared: TDG6PreparedGR0Compositor) -> str:
    parts = [
        prepared.method,
        prepared.initial_state_sha256,
        float(prepared.initial_time).hex(),
        float(prepared.final_time).hex(),
        str(prepared.previous_step_index),
        str(prepared.previous_transaction_serial),
    ]
    for path in (prepared.outer, prepared.medium, prepared.fine):
        parts.append(path.initial_state_sha256)
        for attempt in path.attempts:
            proposal = attempt.proposal
            parts.append(
                array_content_sha256(
                    proposal.initial_state.u,
                    proposal.initial_state.p,
                    proposal.initial_state.q,
                    proposal.candidate_state.u,
                    proposal.candidate_state.p,
                    proposal.candidate_state.q,
                )
            )
            for record in proposal.stages:
                parts.append(
                    array_content_sha256(record.rhs.du, record.rhs.dp, record.rhs.dq)
                )
    return sha256("\n".join(parts).encode("ascii")).hexdigest()


def _external_snapshot(state, transaction, tracers, ledger) -> dict[str, object]:
    return {
        "state": array_content_sha256(state.u, state.p, state.q),
        "monitor": transaction.state,
        "causal": transaction.causal_state,
        "tracer_positions": array_content_sha256(tracers.positions),
        "tracer_proper_times": array_content_sha256(tracers.proper_times),
        "ledger": ledger,
    }


def _closed_evidence(*, reason: str) -> ExactCompleteCClosedEvidence:
    return ExactCompleteCClosedEvidence(
        reason=reason,
        level="D01",
        row_count=OWNED_ROW_COUNT,
        maximum_candidates_D01=CEILING_D01,
        maximum_candidates_D12=CEILING_D12,
        refinement_depth=DEPTH,
        detail="test",
    )


class TDG10ExactCompleteCRuntimeTests(unittest.TestCase):
    passing_evidence = None
    failing_evidence = None

    @classmethod
    def setUpClass(cls) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        rows = (_row(zero, zero),) * OWNED_ROW_COUNT
        cls.passing_evidence = assess_exact_complete_c_rows(
            rows,
            maximum_candidates_D01=CEILING_D01,
            maximum_candidates_D12=CEILING_D12,
            refinement_depth=DEPTH,
        )
        cls.failing_evidence = assess_exact_complete_c_rows(
            (_row((Q(2), Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(0), Q(0))),)
            * OWNED_ROW_COUNT,
            maximum_candidates_D01=CEILING_D01,
            maximum_candidates_D12=CEILING_D12,
            refinement_depth=DEPTH,
        )

    def assess(self, prepared, fake):
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime."
            "assess_exact_complete_c_rows",
            side_effect=fake,
        ) as mocked:
            result = assess_tdg10_exact_complete_c(
                prepared,
                maximum_candidates_D01=CEILING_D01,
                maximum_candidates_D12=CEILING_D12,
                refinement_depth=DEPTH,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )
        return result, mocked

    def test_module_is_pure_adapter_and_isolates_private_tdg5_helpers(self) -> None:
        source = MODULE_PATH.read_text()
        imported = _imported_modules(MODULE_PATH)
        self.assertIn("tdg6_temporal_admission_runtime", imported)
        self.assertIn("TDG6PreparedGR0Compositor", imported)
        self.assertIn("tdg10_exact_complete_c_admission", imported)
        self.assertIn("assess_exact_complete_c_rows", imported)
        self.assertIn("tdg5_stage_complete_refinement_runtime", imported)
        self.assertIn("tdg5._owned_state", source)
        self.assertIn("tdg5._owned_rhs", source)
        self.assertIn("tdg5._proposal_endpoint_records", source)
        self.assertIn("Historical private helpers", source)
        self.assertNotIn("prepare_tdg6_gr0_compositor(", source)
        self.assertNotIn("commit_tdg6_gr0_compositor(", source)
        self.assertNotIn("attempt_proto7_step", source)
        self.assertNotIn(".continuous_admission", source)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )
        test_imported = _imported_modules(Path(__file__))
        self.assertFalse(
            any("binder" in name or "scripts" in name for name in test_imported)
        )

    def test_channel_order_and_one_two_four_segmentation(self) -> None:
        prepared, *_ = _prepare()
        captured: list[tuple[str, tuple[object, ...], dict[str, int]]] = []

        def fake(rows, **kwargs):
            materialized = tuple(rows)
            captured.append(
                (
                    TDG6_COMPLETE_STATE_CHANNELS[len(captured)],
                    materialized,
                    {
                        "maximum_candidates_D01": kwargs["maximum_candidates_D01"],
                        "maximum_candidates_D12": kwargs["maximum_candidates_D12"],
                        "refinement_depth": kwargs["refinement_depth"],
                    },
                )
            )
            return self.passing_evidence

        result, mocked = self.assess(prepared, fake)
        self.assertEqual(mocked.call_count, 18)
        self.assertEqual(
            tuple(channel for channel, _, _ in captured), TDG6_COMPLETE_STATE_CHANNELS
        )
        self.assertEqual(
            tuple(item.channel for item in result.channel_evidence),
            TDG6_COMPLETE_STATE_CHANNELS,
        )
        outer_width = prepared.final_time - prepared.initial_time
        for channel, rows, request in captured:
            self.assertEqual(len(rows), OWNED_ROW_COUNT)
            self.assertEqual(request["maximum_candidates_D01"], CEILING_D01)
            self.assertEqual(request["maximum_candidates_D12"], CEILING_D12)
            self.assertEqual(request["refinement_depth"], DEPTH)
            for row in rows:
                outer, medium, fine = row
                self.assertEqual(len(outer), 5)
                self.assertEqual(len(medium), 2)
                self.assertEqual(len(fine), 4)
                self.assertTrue(all(len(item) == 5 for item in medium))
                self.assertTrue(all(len(item) == 5 for item in fine))
                self.assertTrue(
                    all(type(value) is float for value in (outer + sum(medium, ()) + sum(fine, ())))
                )
                self.assertEqual(outer[4], outer_width)
                self.assertTrue(all(item[4] == outer_width / 2.0 for item in medium))
                self.assertTrue(all(item[4] == outer_width / 4.0 for item in fine))
            _ = channel

    def test_owned_row_slicing_excludes_centre_and_projector_rows(self) -> None:
        prepared, *_ = _prepare()
        captured: dict[str, tuple[object, ...]] = {}

        def fake(rows, **kwargs):
            captured[TDG6_COMPLETE_STATE_CHANNELS[len(captured)]] = tuple(rows)
            return self.passing_evidence

        self.assess(prepared, fake)
        owned_rows = range(
            OWNED_INNER_ROW, POINT_COUNT - PROJECTOR_OWNED_OUTER_ROWS
        )
        excluded = (0, *range(POINT_COUNT - PROJECTOR_OWNED_OUTER_ROWS, POINT_COUNT))
        self.assertEqual(tuple(owned_rows), (1, 2, 3, 4))
        self.assertEqual(excluded, (0, 5, 6, 7, 8))
        for block_index, block in enumerate(_BLOCK_NAMES):
            for field_index, field in enumerate(_FIELD_NAMES):
                channel = f"{block}:{field}"
                left_values = tuple(
                    row[0][0] for row in captured[channel]
                )
                expected = tuple(
                    _marker(block_index, field_index, row) for row in owned_rows
                )
                self.assertEqual(left_values, expected)
                forbidden = {
                    _marker(block_index, field_index, row) for row in excluded
                }
                self.assertTrue(forbidden.isdisjoint(left_values))
                for row in captured[channel]:
                    outer, medium, fine = row
                    self.assertEqual(outer[0], outer[2])
                    self.assertEqual(outer[1], 0.0)
                    self.assertEqual(outer[3], 0.0)
                    self.assertTrue(all(item[0] == outer[0] and item[2] == outer[0] for item in medium))
                    self.assertTrue(all(item[0] == outer[0] and item[2] == outer[0] for item in fine))

    def test_all_of_pass(self) -> None:
        prepared, *_ = _prepare()
        result, mocked = self.assess(prepared, lambda rows, **kwargs: self.passing_evidence)
        self.assertEqual(mocked.call_count, 18)
        self.assertTrue(result.complete_admission_passed)
        self.assertEqual(result.failed_channels, ())
        self.assertEqual(result.owned_row_count, OWNED_ROW_COUNT)
        self.assertEqual(
            result.finest_pair_exact_upper_vector, (Q(0),) * 18
        )
        self.assertTrue(
            all(item.evidence.decision.admission_passed for item in result.channel_evidence)
        )
        self.assertFalse(result.state_advance_authorized)
        self.assertFalse(result.production_method_earned)

    def test_all_of_fail_records_failed_channels(self) -> None:
        prepared, *_ = _prepare()
        failed_channel = "u:R"
        self.assertEqual(TDG6_COMPLETE_STATE_CHANNELS[3], failed_channel)

        def fake(rows, **kwargs):
            if len(fake.calls) == 3:
                fake.calls.append(failed_channel)
                return self.failing_evidence
            fake.calls.append(TDG6_COMPLETE_STATE_CHANNELS[len(fake.calls)])
            return self.passing_evidence

        fake.calls = []
        result, mocked = self.assess(prepared, fake)
        self.assertEqual(mocked.call_count, 18)
        self.assertFalse(result.complete_admission_passed)
        self.assertEqual(result.failed_channels, (failed_channel,))
        self.assertFalse(result.channel_evidence[3].evidence.decision.admission_passed)
        self.assertTrue(
            all(
                item.evidence.decision.admission_passed
                for item in result.channel_evidence
                if item.channel != failed_channel
            )
        )
        self.assertEqual(result.finest_pair_exact_upper_vector[3], Q(1))
        self.assertTrue(
            all(
                item == Q(0)
                for index, item in enumerate(result.finest_pair_exact_upper_vector)
                if index != 3
            )
        )
        self.assertFalse(result.state_advance_authorized)

    def test_real_small_exact_core_exact_zero_integration(self) -> None:
        prepared, state, transaction, tracers, ledger = _prepare()
        before = (
            _compositor_digest(prepared),
            _external_snapshot(state, transaction, tracers, ledger),
        )
        result = assess_tdg10_exact_complete_c(
            prepared,
            maximum_candidates_D01=CEILING_D01,
            maximum_candidates_D12=CEILING_D12,
            refinement_depth=DEPTH,
            expected_owned_row_count=OWNED_ROW_COUNT,
        )
        self.assertTrue(result.complete_admission_passed)
        self.assertEqual(result.failed_channels, ())
        self.assertEqual(len(result.channel_evidence), 18)
        for item in result.channel_evidence:
            self.assertEqual(item.evidence.decision.classification, "exact_zero")
            self.assertTrue(item.evidence.decision.admission_passed)
            self.assertEqual(item.evidence.d01.upper, 0)
            self.assertEqual(item.evidence.d12.upper, 0)
            self.assertEqual(item.evidence.row_count, OWNED_ROW_COUNT)
            self.assertNotIn("candidates", item.evidence.__dataclass_fields__)
        self.assertEqual(result.finest_pair_exact_upper_vector, (Q(0),) * 18)
        self.assertEqual(_compositor_digest(prepared), before[0])
        self.assertEqual(
            _external_snapshot(state, transaction, tracers, ledger), before[1]
        )

    def test_no_mutation_of_prepared_or_external_ledgers(self) -> None:
        prepared, state, transaction, tracers, ledger = _prepare()
        before_prepared = _compositor_digest(prepared)
        before_external = _external_snapshot(state, transaction, tracers, ledger)
        result, _ = self.assess(
            prepared, lambda rows, **kwargs: self.passing_evidence
        )
        self.assertIsInstance(result, TDG10ExactCompleteCRuntimeEvidence)
        self.assertEqual(_compositor_digest(prepared), before_prepared)
        self.assertEqual(
            _external_snapshot(state, transaction, tracers, ledger), before_external
        )
        self.assertIs(state.u.flags.writeable, False)
        self.assertEqual(transaction.state.accepted_stage_count, 0)
        self.assertEqual(ledger.accepted_macro_step_count, 0)

    def test_malformed_paths_and_component_mismatch_fail_closed(self) -> None:
        prepared, *_ = _prepare()
        with self.assertRaisesRegex(TypeError, "TDG6PreparedGR0Compositor"):
            assess_tdg10_exact_complete_c(
                object(),
                maximum_candidates_D01=CEILING_D01,
                maximum_candidates_D12=CEILING_D12,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )

        rk4, *_ = _prepare(PRIMARY_METHOD)
        with self.assertRaisesRegex(ValueError, "requires SSPRK3"):
            assess_tdg10_exact_complete_c(
                rk4,
                maximum_candidates_D01=CEILING_D01,
                maximum_candidates_D12=CEILING_D12,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )

        relabeled = copy(prepared.outer)
        object.__setattr__(relabeled, "level", "medium")
        swapped = copy(prepared)
        object.__setattr__(swapped, "outer", relabeled)
        with self.assertRaisesRegex(ValueError, "wrong level name"):
            assess_tdg10_exact_complete_c(
                swapped,
                maximum_candidates_D01=CEILING_D01,
                maximum_candidates_D12=CEILING_D12,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )

        truncated = copy(prepared.fine)
        object.__setattr__(truncated, "attempts", prepared.fine.attempts[:-1])
        broken_fine = copy(prepared)
        object.__setattr__(broken_fine, "fine", truncated)
        with self.assertRaisesRegex(ValueError, "1/2/4"):
            assess_tdg10_exact_complete_c(
                broken_fine,
                maximum_candidates_D01=CEILING_D01,
                maximum_candidates_D12=CEILING_D12,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )

        with self.assertRaisesRegex(ValueError, "component/order"):
            assess_tdg10_exact_complete_c(
                prepared,
                maximum_candidates_D01=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01,
                maximum_candidates_D12=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12,
            )

        poisoned = np.zeros((3, OWNED_ROW_COUNT, 5), dtype=np.float64)
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime._owned_state",
            return_value=poisoned,
        ):
            with self.assertRaisesRegex(ValueError, "component/order"):
                assess_tdg10_exact_complete_c(
                    prepared,
                    maximum_candidates_D01=CEILING_D01,
                    maximum_candidates_D12=CEILING_D12,
                    expected_owned_row_count=OWNED_ROW_COUNT,
                )

        poisoned_dtype = np.zeros(
            (3, OWNED_ROW_COUNT, 6), dtype=np.float32
        )
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime._owned_state",
            return_value=poisoned_dtype,
        ):
            with self.assertRaisesRegex(ValueError, "little-endian binary64"):
                assess_tdg10_exact_complete_c(
                    prepared,
                    maximum_candidates_D01=CEILING_D01,
                    maximum_candidates_D12=CEILING_D12,
                    expected_owned_row_count=OWNED_ROW_COUNT,
                )

        poisoned_layout = np.zeros(
            (3, OWNED_ROW_COUNT, 6), dtype=np.float64
        )[:, ::-1, :]
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime._owned_state",
            return_value=poisoned_layout,
        ):
            with self.assertRaisesRegex(ValueError, "C-contiguous"):
                assess_tdg10_exact_complete_c(
                    prepared,
                    maximum_candidates_D01=CEILING_D01,
                    maximum_candidates_D12=CEILING_D12,
                    expected_owned_row_count=OWNED_ROW_COUNT,
                )

        with self.assertRaisesRegex(TypeError, "must be an integer"):
            assess_tdg10_exact_complete_c(
                prepared,
                maximum_candidates_D01=True,  # type: ignore[arg-type]
                maximum_candidates_D12=CEILING_D12,
                expected_owned_row_count=OWNED_ROW_COUNT,
            )

    def test_shared_node_rhs_identity_and_production_limits_fail_closed(self) -> None:
        prepared, *_ = _prepare()
        real_endpoint_records = tdg10_runtime._proposal_endpoint_records

        first_medium = prepared.medium.attempts[0].proposal
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime."
            "_proposal_endpoint_records",
            side_effect=lambda proposal: (
                (_poisoned_rhs(real_endpoint_records(proposal)[0]), real_endpoint_records(proposal)[1])
                if proposal is first_medium
                else real_endpoint_records(proposal)
            ),
        ):
            with self.assertRaisesRegex(ValueError, "shared initial state produced different RHS"):
                self.assess(prepared, lambda rows, **kwargs: self.passing_evidence)

        second_medium = prepared.medium.attempts[1].proposal
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime."
            "_proposal_endpoint_records",
            side_effect=lambda proposal: (
                (_poisoned_rhs(real_endpoint_records(proposal)[0]), real_endpoint_records(proposal)[1])
                if proposal is second_medium
                else real_endpoint_records(proposal)
            ),
        ):
            with self.assertRaisesRegex(ValueError, "shared state produced different RHS"):
                self.assess(prepared, lambda rows, **kwargs: self.passing_evidence)

        with self.assertRaisesRegex(ValueError, "frozen D01/D12 candidate ceilings"):
            assess_tdg10_exact_complete_c(
                prepared,
                maximum_candidates_D01=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01 + 1,
                maximum_candidates_D12=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12,
                refinement_depth=TDG10_PRODUCTION_REFINEMENT_DEPTH,
                expected_owned_row_count=TDG10_PRODUCTION_OWNED_ROW_COUNT,
            )

        self.assertEqual(
            tdg10_runtime._builtin_finite_float(np.float64(1.0), name="test"),
            1.0,
        )
        with self.assertRaisesRegex(TypeError, "built-in binary64"):
            tdg10_runtime._builtin_finite_float(np.float32(1.0), name="test")

    def test_nonfinite_values_fail_closed(self) -> None:
        prepared, *_ = _prepare()
        poisoned = np.zeros((3, OWNED_ROW_COUNT, 6), dtype=np.float64)
        poisoned[0, 0, 0] = np.nan
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime._owned_state",
            return_value=poisoned,
        ):
            with self.assertRaisesRegex(ValueError, "must be finite"):
                assess_tdg10_exact_complete_c(
                    prepared,
                    maximum_candidates_D01=CEILING_D01,
                    maximum_candidates_D12=CEILING_D12,
                    expected_owned_row_count=OWNED_ROW_COUNT,
                )
        poisoned_inf = np.zeros((3, OWNED_ROW_COUNT, 6), dtype=np.float64)
        poisoned_inf[1, 2, 3] = np.inf
        with patch(
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime._owned_rhs",
            return_value=poisoned_inf,
        ):
            with self.assertRaisesRegex(ValueError, "must be finite"):
                assess_tdg10_exact_complete_c(
                    prepared,
                    maximum_candidates_D01=CEILING_D01,
                    maximum_candidates_D12=CEILING_D12,
                    expected_owned_row_count=OWNED_ROW_COUNT,
                )

    def test_exact_core_exceptions_are_not_reinterpreted_as_physics(self) -> None:
        prepared, *_ = _prepare()
        exhausted = ExactCompleteCResourceExhausted(
            _closed_evidence(reason="candidate_ceiling_exhausted")
        )
        with self.assertRaises(TDG10ExactCompleteCRuntimeClosed) as closed:
            self.assess(prepared, lambda rows, **kwargs: (_ for _ in ()).throw(exhausted))
        self.assertIs(closed.exception.physical_classification, False)
        self.assertEqual(closed.exception.channel, TDG6_COMPLETE_STATE_CHANNELS[0])
        self.assertIs(closed.exception.cause, exhausted)
        self.assertEqual(
            closed.exception.evidence.reason, "candidate_ceiling_exhausted"
        )

        disagreed = ExactCompleteCRouteDisagreement(
            _closed_evidence(reason="route_disagreement")
        )
        with self.assertRaises(TDG10ExactCompleteCRuntimeClosed) as closed:
            self.assess(prepared, lambda rows, **kwargs: (_ for _ in ()).throw(disagreed))
        self.assertIs(closed.exception.physical_classification, False)
        self.assertEqual(closed.exception.cause.evidence.reason, "route_disagreement")

    def test_frozen_nonclaims_and_compact_result(self) -> None:
        prepared, *_ = _prepare()
        result, _ = self.assess(
            prepared, lambda rows, **kwargs: self.passing_evidence
        )
        self.assertEqual(result.method, COMPARATOR_METHOD)
        self.assertTrue(result.admission_is_all_of)
        self.assertFalse(result.state_advance_authorized)
        self.assertFalse(result.production_method_earned)
        self.assertFalse(result.PDE_state_committed)
        self.assertFalse(result.GR0_calibration_completed)
        self.assertFalse(result.candidate_branch_opened)
        names = {item.name for item in fields(result)}
        self.assertNotIn("outer", names)
        self.assertNotIn("medium", names)
        self.assertNotIn("fine", names)
        self.assertNotIn("prepared", names)
        self.assertNotIn("candidates", names)
        self.assertNotIn("prepared_tdg6", names)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, state_advance_authorized=True)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, production_method_earned=True)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, PDE_state_committed=True)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, GR0_calibration_completed=True)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, candidate_branch_opened=True)
        with self.assertRaisesRegex(ValueError, "nonclaim"):
            replace(result, admission_is_all_of=False)
        self.assertEqual(TDG10_PRODUCTION_OWNED_ROW_COUNT, 2044)

    def test_ssprk3_required_and_default_owned_row_count_is_production(self) -> None:
        prepared, *_ = _prepare()
        with self.assertRaisesRegex(ValueError, "component/order"):
            assess_tdg10_exact_complete_c(
                prepared,
                maximum_candidates_D01=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01,
                maximum_candidates_D12=TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12,
                refinement_depth=TDG10_PRODUCTION_REFINEMENT_DEPTH,
            )
        default = inspect.signature(
            assess_tdg10_exact_complete_c
        ).parameters["expected_owned_row_count"].default
        self.assertEqual(default, 2044)
        self.assertEqual(default, TDG10_PRODUCTION_OWNED_ROW_COUNT)


if __name__ == "__main__":
    unittest.main()
