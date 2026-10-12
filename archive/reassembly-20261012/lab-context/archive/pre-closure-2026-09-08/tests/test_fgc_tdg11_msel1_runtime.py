"""Synthetic guarded SHADOW family controls; no campaign or runner inputs."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from recursive_horizons.fgc.evolution.boundary_domain import (
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (
    AcceptedStateSourceFailureEvidence,
    ProposalFailureEvidence,
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.tdg11_compensated_rk import (
    ARITHMETIC_ID as COMPENSATED_ARITHMETIC_ID,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_reconstruction import (
    ORIGINAL_ARITHMETIC_ID,
    ValidatedRecordedFamily,
)
from recursive_horizons.fgc.evolution import numerical_engine as numeric_engine
from recursive_horizons.fgc.evolution import proto7_runtime as proto7
from recursive_horizons.fgc.evolution import tdg11_compensated_rk as compensated_rk
from recursive_horizons.fgc.evolution import tdg11_msel1_runtime as runtime


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_msel1_runtime.py"
POINT_COUNT = 9
STEP_SIZE = 0.125
FINE_TWO_MIDPOINT = 0.078125
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
    "tdg10_",
    "hlt16",
    "binder",
    "scripts",
    "campaign",
    "store",
    "classify_tdg6",
)


def _imported_modules(path: Path) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module)
            names.update(alias.name for alias in node.names)
    return names


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, POINT_COUNT)
    u = np.zeros((POINT_COUNT, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(*, failed_time: float | None = None, source_only: bool = False):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        failed = failed_time is not None and (
            np.float64(time).tobytes() == np.float64(failed_time).tobytes()
        )
        return EvolutionRHS(
            state.u,
            state.p,
            state.q,
            {
                "source_residual_infinity": 2.0 if failed and source_only else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed and not source_only else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


class _CountingRHS:
    def __init__(self, inner):
        self.inner = inner
        self.calls = 0
        self.times: list[float] = []

    def __call__(self, time: float, state: EvolutionState) -> EvolutionRHS:
        self.calls += 1
        self.times.append(float(time))
        return self.inner(time, state)


class _Tracers:
    def __init__(self, *, fail_preview_at: int | None = None) -> None:
        self.labels = np.asarray((0.0, 1.0, 2.0), dtype=np.float64)
        self.positions = np.asarray((0.25, 0.50, 0.75), dtype=np.float64)
        self.proper_times = np.zeros(3, dtype=np.float64)
        self.event_proper_times: list[np.ndarray] = []
        self.event_fields: list[np.ndarray] = []
        self.cutoff = 1.0
        self.outer_radius = 128.0
        self.preview_calls = 0
        self.commit_calls = 0
        self.fail_preview_at = fail_preview_at

    def preview_advance(self, *, old_state, new_state, coordinates, step_size):
        self.preview_calls += 1
        if (
            self.fail_preview_at is not None
            and self.preview_calls == self.fail_preview_at
        ):
            raise RuntimeError("injected tracer preview failure")
        return (
            np.asarray(self.positions, dtype=np.float64).copy() + float(step_size),
            np.asarray(self.proper_times, dtype=np.float64).copy() + float(step_size),
        )

    def commit_advance(self, positions, proper_times) -> None:
        self.commit_calls += 1
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


def _coordinates() -> np.ndarray:
    return np.linspace(0.0, 1.0, POINT_COUNT)


def _array_fingerprint(array: np.ndarray) -> tuple[object, ...]:
    return (
        id(array),
        array.dtype.str,
        tuple(array.shape),
        tuple(array.strides),
        array.tobytes(),
        bool(array.flags.writeable),
        bool(array.flags.c_contiguous),
        bool(array.flags.owndata),
    )


def _fingerprint(state, transaction, tracers, coordinates) -> tuple[object, ...]:
    return (
        id(state),
        _array_fingerprint(state.u),
        _array_fingerprint(state.p),
        _array_fingerprint(state.q),
        array_content_sha256(state.u, state.p, state.q),
        id(transaction),
        id(transaction.state),
        id(transaction.causal_state),
        transaction.state,
        transaction.causal_state,
        transaction.grid_spacing,
        transaction.cfl_maximum,
        id(tracers),
        tracers.preview_calls,
        tracers.commit_calls,
        tracers.cutoff,
        tracers.outer_radius,
        _array_fingerprint(tracers.labels),
        _array_fingerprint(tracers.positions),
        _array_fingerprint(tracers.proper_times),
        tuple(_array_fingerprint(item) for item in tracers.event_proper_times),
        tuple(_array_fingerprint(item) for item in tracers.event_fields),
        _array_fingerprint(coordinates),
    )


def _build(
    *,
    method: str = PRIMARY_METHOD,
    arithmetic_id: str = ORIGINAL_ARITHMETIC_ID,
    rhs=None,
    state=None,
    transaction=None,
    tracers=None,
    coordinates=None,
    step_size: float = STEP_SIZE,
):
    initial = _state() if state is None else state
    return runtime.build_guarded_shadow_family(
        method=method,
        arithmetic_id=arithmetic_id,
        time=0.0,
        step_size=step_size,
        state=initial,
        rhs=_rhs() if rhs is None else rhs,
        projector=None,
        transaction=_transaction() if transaction is None else transaction,
        tracers=_Tracers() if tracers is None else tracers,
        coordinates=_coordinates() if coordinates is None else coordinates,
        previous_step_index=0,
        previous_transaction_serial=0,
    )


class ImportBoundaryTests(unittest.TestCase):
    def test_module_isolates_private_helpers_and_avoids_historical_runners(
        self,
    ) -> None:
        source = MODULE_PATH.read_text()
        imported = _imported_modules(MODULE_PATH)
        self.assertIn("tdg5_stage_complete_refinement_runtime", imported)
        self.assertIn("tdg6_temporal_admission_runtime", imported)
        self.assertIn("tdg11_msel1_reconstruction", imported)
        self.assertIn("tdg11_compensated_rk", imported)
        self.assertIn("validate_recorded_family", imported)
        self.assertIn("_latch_accepted_source_failure", imported)
        self.assertIn("_complete_proposal_failure_evidence", imported)
        self.assertIn("propose_compensated_step", imported)
        self.assertIn("tdg5._clone_gr0_transaction", source)
        self.assertIn("tdg6._clone_tracers", source)
        self.assertIn("numeric_engine.propose_step", source)
        self.assertIn("guards=(transaction,)", source)
        self.assertNotIn("attempt_proto7_step", source)
        self.assertNotIn("prepare_tdg6_gr0_compositor", source)
        self.assertNotIn("commit_tdg6_gr0_compositor", source)
        self.assertNotIn("TDG6TemporalLedger", source)
        self.assertNotIn("classify_tdg6", source)
        self.assertNotIn("proto7.propose_step =", source)
        self.assertNotIn("numeric_engine.propose_step =", source)
        self.assertNotIn("propose_step =", source)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )
        test_imported = _imported_modules(Path(__file__))
        self.assertFalse(
            any("binder" in name or "scripts" in name for name in test_imported)
        )
        self.assertNotIn("propose_step", runtime.__dict__)


class GuardedShadowFamilyTests(unittest.TestCase):
    def test_both_methods_and_arithmetic_ids_build_one_two_four_family(self) -> None:
        for method, stages, rhs_calls in (
            (PRIMARY_METHOD, 35, 42),
            (COMPARATOR_METHOD, 28, 35),
        ):
            for arithmetic_id in (ORIGINAL_ARITHMETIC_ID, COMPENSATED_ARITHMETIC_ID):
                with self.subTest(method=method, arithmetic_id=arithmetic_id):
                    counter = _CountingRHS(_rhs())
                    state = _state()
                    transaction = _transaction()
                    tracers = _Tracers()
                    coordinates = _coordinates()
                    before = _fingerprint(state, transaction, tracers, coordinates)
                    family = _build(
                        method=method,
                        arithmetic_id=arithmetic_id,
                        rhs=counter,
                        state=state,
                        transaction=transaction,
                        tracers=tracers,
                        coordinates=coordinates,
                    )
                    self.assertIsInstance(family, runtime.GuardedShadowFamily)
                    self.assertEqual(
                        tuple(path.level for path in family.paths),
                        ("outer", "medium", "fine"),
                    )
                    self.assertEqual(
                        tuple(len(path.attempts) for path in family.paths),
                        (1, 2, 4),
                    )
                    self.assertIsInstance(family.recorded, ValidatedRecordedFamily)
                    self.assertEqual(family.recorded.method, method)
                    self.assertEqual(family.recorded.arithmetic_id, arithmetic_id)
                    self.assertEqual(family.accepted_source_precheck_count, 7)
                    self.assertEqual(family.stage_record_count, stages)
                    self.assertEqual(family.rhs_call_count, rhs_calls)
                    self.assertEqual(counter.calls, rhs_calls)
                    self.assertEqual(
                        counter.calls,
                        family.accepted_source_precheck_count
                        + family.stage_record_count,
                    )
                    self.assertEqual(
                        _fingerprint(state, transaction, tracers, coordinates), before
                    )
                    self.assertEqual(tracers.preview_calls, 0)
                    self.assertEqual(tracers.commit_calls, 0)
                    with self.assertRaises(FrozenInstanceError):
                        family.rhs_call_count = 0  # type: ignore[misc]

    def test_stage_counts_use_fresh_endpoint_calls(self) -> None:
        for method, stages_per_proposal in (
            (PRIMARY_METHOD, 5),
            (COMPARATOR_METHOD, 4),
        ):
            with self.subTest(method=method):
                counter = _CountingRHS(_rhs())
                family = _build(method=method, rhs=counter)
                endpoint_times = []
                for path in family.paths:
                    for attempt in path.attempts:
                        stages = attempt.proposal.stages
                        self.assertEqual(len(stages), stages_per_proposal)
                        self.assertEqual(stages[-1].stage_name, "candidate_endpoint")
                        self.assertIsNot(stages[-1].rhs, stages[-2].rhs)
                        self.assertEqual(
                            np.float64(stages[-1].time).tobytes(),
                            np.float64(attempt.proposal.final_time).tobytes(),
                        )
                        endpoint_times.append(stages[-1].time)
                        self.assertIsNotNone(attempt.accepted)
                        self.assertIsNone(attempt.retry)
                self.assertEqual(len(endpoint_times), 7)
                self.assertGreaterEqual(
                    sum(
                        np.float64(time).tobytes() == np.float64(item).tobytes()
                        for time in counter.times
                        for item in endpoint_times
                    ),
                    7,
                )

    def test_accepted_source_failure_stops_family_and_preserves_inputs(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _Tracers()
        coordinates = _coordinates()
        counter = _CountingRHS(_rhs(failed_time=0.0, source_only=True))
        before = _fingerprint(state, transaction, tracers, coordinates)
        with self.assertRaises(runtime.TDG11ShadowPremiseStop) as caught:
            _build(
                rhs=counter,
                state=state,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
        stop = caught.exception
        self.assertEqual(stop.level, "outer")
        self.assertEqual(stop.proposal_index, 0)
        self.assertIsInstance(stop.cause, Proto7TerminalStop)
        self.assertIsInstance(stop.cause.evidence, AcceptedStateSourceFailureEvidence)
        self.assertFalse(stop.cause.evidence.proposal_constructed)
        self.assertIsNone(stop.source_retry)
        self.assertFalse(stop.scientific_nonpass)
        self.assertFalse(stop.informal_retry)
        self.assertFalse(stop.temporal_retry)
        self.assertEqual(stop.accepted_source_precheck_count, 1)
        self.assertEqual(stop.stage_record_count, 0)
        self.assertEqual(stop.rhs_call_count, 1)
        self.assertEqual(counter.calls, 1)
        self.assertEqual(_fingerprint(state, transaction, tracers, coordinates), before)
        self.assertEqual(transaction.state, GR0RuntimeMonitorState())

    def test_later_source_only_retry_stops_without_informal_retry(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _Tracers()
        coordinates = _coordinates()
        counter = _CountingRHS(_rhs(failed_time=FINE_TWO_MIDPOINT, source_only=True))
        before = _fingerprint(state, transaction, tracers, coordinates)
        with self.assertRaises(runtime.TDG11ShadowPremiseStop) as caught:
            _build(
                rhs=counter,
                state=state,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
        stop = caught.exception
        self.assertEqual(stop.level, "fine")
        self.assertEqual(stop.proposal_index, 2)
        self.assertIsNone(stop.cause)
        self.assertIsInstance(stop.source_retry, ProposalFailureEvidence)
        self.assertTrue(stop.source_retry.retryable_source_only)
        self.assertTrue(stop.source_retry_owned_by_PROTO7)
        self.assertFalse(stop.informal_retry)
        self.assertFalse(stop.temporal_retry)
        self.assertFalse(stop.scientific_nonpass)
        self.assertEqual(stop.accepted_source_precheck_count, 6)
        self.assertEqual(stop.stage_record_count, 30)
        self.assertEqual(stop.rhs_call_count, 36)
        self.assertEqual(counter.calls, 36)
        self.assertEqual(_fingerprint(state, transaction, tracers, coordinates), before)
        self.assertEqual(transaction.state.accepted_stage_count, 0)

    def test_non_source_veto_preserves_owner_and_inputs(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _Tracers()
        coordinates = _coordinates()
        counter = _CountingRHS(_rhs(failed_time=FINE_TWO_MIDPOINT, source_only=False))
        before = _fingerprint(state, transaction, tracers, coordinates)
        with self.assertRaises(runtime.TDG11ShadowPremiseStop) as caught:
            _build(
                rhs=counter,
                state=state,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
        stop = caught.exception
        self.assertEqual(stop.level, "fine")
        self.assertEqual(stop.proposal_index, 2)
        self.assertIsInstance(stop.cause, Proto7TerminalStop)
        self.assertIsInstance(stop.cause.evidence, ProposalFailureEvidence)
        self.assertTrue(stop.cause.evidence.non_source_failure_vetoed_retry)
        self.assertEqual(
            stop.cause.evidence.selected_terminal_reason, "nonpositive_lapse"
        )
        self.assertEqual(stop.cause.reason, "nonpositive_lapse")
        self.assertIsNone(stop.source_retry)
        self.assertFalse(stop.scientific_nonpass)
        self.assertEqual(stop.accepted_source_precheck_count, 6)
        self.assertEqual(stop.stage_record_count, 30)
        self.assertEqual(stop.rhs_call_count, 36)
        self.assertEqual(counter.calls, 36)
        self.assertEqual(_fingerprint(state, transaction, tracers, coordinates), before)

    def test_tracer_preview_failure_stops_family(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _Tracers(fail_preview_at=1)
        coordinates = _coordinates()
        counter = _CountingRHS(_rhs())
        before = _fingerprint(state, transaction, tracers, coordinates)
        with self.assertRaises(runtime.TDG11ShadowPremiseStop) as caught:
            _build(
                rhs=counter,
                state=state,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
        stop = caught.exception
        self.assertEqual(stop.level, "outer")
        self.assertEqual(stop.proposal_index, 0)
        self.assertIsInstance(stop.cause, RuntimeError)
        self.assertIn("tracer preview", str(stop.cause))
        self.assertIsNone(stop.source_retry)
        self.assertFalse(stop.scientific_nonpass)
        self.assertEqual(stop.accepted_source_precheck_count, 1)
        self.assertEqual(stop.stage_record_count, 5)
        self.assertEqual(stop.rhs_call_count, 6)
        self.assertEqual(counter.calls, 6)
        self.assertEqual(_fingerprint(state, transaction, tracers, coordinates), before)
        self.assertEqual(tracers.preview_calls, 0)
        self.assertEqual(tracers.commit_calls, 0)

    def test_no_historical_global_proposer_substitution(self) -> None:
        proto7_propose = proto7.propose_step
        engine_propose = numeric_engine.propose_step
        compensated_propose = compensated_rk.propose_compensated_step
        attempt = proto7.attempt_proto7_step
        with patch.object(
            proto7, "propose_step", side_effect=AssertionError("global proposer")
        ):
            with patch.object(
                proto7,
                "attempt_proto7_step",
                side_effect=AssertionError("historical attempt"),
            ):
                family = _build()
        self.assertEqual(family.rhs_call_count, 42)
        self.assertIs(proto7.propose_step, proto7_propose)
        self.assertIs(numeric_engine.propose_step, engine_propose)
        self.assertIs(compensated_rk.propose_compensated_step, compensated_propose)
        self.assertIs(proto7.attempt_proto7_step, attempt)
        with patch.object(
            numeric_engine,
            "propose_step",
            side_effect=AssertionError("legacy proposer"),
        ):
            compensated = _build(arithmetic_id=COMPENSATED_ARITHMETIC_ID)
        self.assertEqual(compensated.recorded.arithmetic_id, COMPENSATED_ARITHMETIC_ID)
        with patch.object(
            runtime,
            "propose_compensated_step",
            side_effect=AssertionError("compensated proposer"),
        ):
            legacy = _build()
        self.assertEqual(legacy.recorded.arithmetic_id, ORIGINAL_ARITHMETIC_ID)

    def test_input_mutation_raises_integrity_error(self) -> None:
        cases = ("state", "transaction", "tracers", "coordinates")
        for target in cases:
            with self.subTest(target=target):
                state = _state()
                transaction = _transaction()
                tracers = _Tracers()
                coordinates = _coordinates()
                inner = _rhs()
                mutated = {"done": False}

                def rhs(time: float, current: EvolutionState) -> EvolutionRHS:
                    if not mutated["done"]:
                        mutated["done"] = True
                        if target == "state":
                            array = state.u
                            array.setflags(write=True)
                            array[0, 0] += 1.0
                        elif target == "transaction":
                            transaction.state = GR0RuntimeMonitorState(
                                accepted_stage_count=99
                            )
                        elif target == "tracers":
                            tracers.positions[0] = np.nextafter(
                                tracers.positions[0], np.inf
                            )
                        else:
                            coordinates[0] = np.nextafter(coordinates[0], np.inf)
                    return inner(time, current)

                with self.assertRaises(runtime.TDG11ShadowIntegrityError) as caught:
                    _build(
                        rhs=rhs,
                        state=state,
                        transaction=transaction,
                        tracers=tracers,
                        coordinates=coordinates,
                    )
                self.assertNotIsInstance(
                    caught.exception, runtime.TDG11ShadowPremiseStop
                )
                self.assertFalse(isinstance(caught.exception, Proto7TerminalStop))

    def test_keyboard_interrupt_propagates_without_wrapping(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _Tracers()
        coordinates = _coordinates()
        before = _fingerprint(state, transaction, tracers, coordinates)

        def rhs(_time: float, _state: EvolutionState) -> EvolutionRHS:
            raise KeyboardInterrupt

        with self.assertRaises(KeyboardInterrupt):
            _build(
                rhs=rhs,
                state=state,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
        self.assertEqual(_fingerprint(state, transaction, tracers, coordinates), before)

    def test_unrecognized_ids_and_nonuniform_lattice_are_rejected(self) -> None:
        counter = _CountingRHS(_rhs())
        with self.assertRaises(ValueError):
            _build(method="RK4", rhs=counter)
        with self.assertRaises(ValueError):
            _build(arithmetic_id="not-an-arithmetic", rhs=counter)
        with self.assertRaises(ValueError):
            _build(step_size=0.1, rhs=counter)
        for coordinates in ([0.0] * POINT_COUNT, _coordinates().astype(np.float32)):
            with self.assertRaises(TypeError):
                _build(coordinates=coordinates, rhs=counter)
        self.assertEqual(counter.calls, 0)

    def test_frozen_metadata_alias_mutation_is_detected(self) -> None:
        transaction = _transaction()
        inner = _rhs()
        changed = False

        def rhs(time, state):
            nonlocal changed
            if not changed:
                changed = True
                object.__setattr__(transaction.state, "accepted_stage_count", 99)
            return inner(time, state)

        with self.assertRaises(runtime.TDG11ShadowIntegrityError):
            _build(rhs=rhs, transaction=transaction)

    def test_counter_alias_types_are_rejected(self) -> None:
        family = _build()
        for name in (
            "accepted_source_precheck_count",
            "stage_record_count",
            "rhs_call_count",
        ):
            with self.assertRaises(ValueError):
                replace(family, **{name: float(getattr(family, name))})
