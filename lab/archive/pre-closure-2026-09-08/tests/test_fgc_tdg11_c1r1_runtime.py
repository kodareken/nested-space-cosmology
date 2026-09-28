"""C1R1 runtime adapter: provenance, shadows, ledger and fine-only commit."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_c1r1_runtime as c1r1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as imp1  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.tdg11_c1r1_enclosure import (  # noqa: E402
    C1R1_IMPLEMENTATION_ID,
    C1R1_REFERENCE_WIRE_EVALUATOR_ID,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_enclosure import (  # noqa: E402
    IMP1EnclosureLimits,
    IMP1_ENCLOSURE_EVALUATOR_ID,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_REJECTION_EVENT_TYPE,
    imp1_checkpoint_extension,
    restore_imp1_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture as _objects,
    synthetic_diagnostics as _diagnostics,
    synthetic_exponential_rhs as _rhs,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)


MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_c1r1_runtime.py"


def _prepare_c1r1(
    method: str = PRIMARY_METHOD, *, rhs=_rhs, saturated: bool = False, target=TARGET
):
    state, transaction, tracers, ledger, coordinates = _objects(
        method, saturated=saturated
    )
    prepared = c1r1.prepare_c1r1_initial_runtime(
        method=method,
        time=START,
        event_target=target,
        requested_cap=1.0 / 32,
        state=state,
        rhs=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger, coordinates


def _prepare_imp1(
    method: str = PRIMARY_METHOD, *, rhs=_rhs, saturated: bool = False, target=TARGET
):
    state, transaction, tracers, ledger, coordinates = _objects(
        method, saturated=saturated
    )
    prepared = imp1.prepare_imp1_initial_runtime(
        method=method,
        time=START,
        event_target=target,
        requested_cap=1.0 / 32,
        state=state,
        rhs=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger, coordinates


def _current(prepared, state, transaction, tracers, ledger, coordinates):
    return dict(
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=ledger,
        current_time=prepared.initial_time,
        current_state=state,
        current_step_index=prepared.previous_step_index,
        current_transaction_serial=prepared.previous_transaction_serial,
        coordinates=coordinates,
    )


def _quarter_impulse(time, current):
    du = np.zeros_like(current.u)
    if time == START + 1.0 / 128:
        du[:, 4] = 1.0
    return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), _diagnostics())


def _prepare_source() -> str:
    source = MODULE_PATH.read_text()
    start = source.index("def _prepare(")
    end = source.index("def prepare_c1r1_initial_runtime(")
    return source[start:end]


class C1R1RuntimeProvenanceTests(unittest.TestCase):
    def test_prepare_uses_c1r1_assessor_and_does_not_wrap_imp1_prepare(self) -> None:
        body = _prepare_source()
        self.assertIn("assess_c1r1_family", body)
        self.assertNotIn("assess_imp1_family", body)
        self.assertNotIn("prepare_imp1_initial_runtime", body)
        self.assertNotIn("prepare_imp1_retry_runtime", body)
        source = MODULE_PATH.read_text()
        self.assertIn("import tdg11_imp1_runtime as _imp1_reference_runtime", source)
        self.assertIn("_imp1_reference_shadow_path", source)

    def test_receipt_distinguishes_implementation_from_reference_wire(self) -> None:
        c1r1_items = _prepare_c1r1()
        imp1_items = _prepare_imp1()
        prepared = c1r1_items[0]
        reference = imp1_items[0]
        self.assertIsInstance(prepared, c1r1.TDG11C1R1PreparedRuntime)
        self.assertNotIsInstance(prepared, imp1.TDG11IMP1PreparedRuntime)
        self.assertEqual(prepared.implementation_id, C1R1_IMPLEMENTATION_ID)
        self.assertEqual(
            prepared.reference_wire_evaluator_id, C1R1_REFERENCE_WIRE_EVALUATOR_ID
        )
        self.assertEqual(prepared.assessment.evaluator_id, IMP1_ENCLOSURE_EVALUATOR_ID)
        self.assertEqual(
            prepared.assessment.as_mapping(), reference.assessment.as_mapping()
        )
        self.assertEqual(
            prepared.recorded_family.family_sha256,
            reference.recorded_family.family_sha256,
        )
        self.assertEqual(
            prepared.assessment.assessment_sha256, reference.assessment.assessment_sha256
        )
        self.assertNotEqual(prepared.receipt_sha256, reference.receipt_sha256)
        mapping = c1r1._prepared_mapping(prepared)
        self.assertEqual(mapping["artifact_id"], "FGC-1-TDG11-C1R1")
        self.assertEqual(mapping["reference_wire_artifact_id"], "FGC-1-TDG11-IMP1")
        self.assertFalse(mapping["corrected_endpoint_committable"])
        self.assertFalse(mapping["production_authority"])


class C1R1RuntimeGuardTests(unittest.TestCase):
    def test_both_methods_commit_fine_only_and_preserve_ledger_wire(self) -> None:
        for method, records, serial in (
            (PRIMARY_METHOD, 35, 20),
            (COMPARATOR_METHOD, 28, 16),
        ):
            with self.subTest(method=method):
                items = _prepare_c1r1(method)
                prepared, state, tx, tracers, ledger, coordinates = items
                initial_state_hash = array_content_sha256(state.u, state.p, state.q)
                inherited = ledger.inherited_snapshot
                old_debit = ledger.accumulated_debit_vector
                self.assertEqual(prepared.stage_record_count, records)
                self.assertTrue(prepared.assessment.admission_passed)
                self.assertEqual(len(prepared.assessment.public_fine_debits), 18)
                committed = c1r1.commit_c1r1_runtime(prepared, **_current(*items))
                self.assertEqual(committed.step_index, 4)
                self.assertEqual(committed.transaction_serial, serial)
                self.assertIs(committed.state, prepared.fine.final_accepted.state)
                self.assertFalse(committed.corrected_endpoint_adopted)
                self.assertFalse(committed.outer_path_adopted)
                self.assertFalse(committed.medium_path_adopted)
                self.assertFalse(committed.campaign_execution_authorized)
                self.assertEqual(committed.implementation_id, C1R1_IMPLEMENTATION_ID)
                self.assertEqual(
                    committed.temporal_ledger.inherited_snapshot, inherited
                )
                self.assertEqual(
                    committed.temporal_ledger.accumulated_debit_vector,
                    tuple(
                        a + b
                        for a, b in zip(
                            old_debit,
                            prepared.assessment.public_fine_debits,
                            strict=True,
                        )
                    ),
                )
                self.assertEqual(
                    array_content_sha256(state.u, state.p, state.q), initial_state_hash
                )
                encoded = imp1_checkpoint_extension(committed.temporal_ledger)
                restored = restore_imp1_checkpoint_extension(encoded)
                self.assertEqual(
                    canonical_json_bytes(encoded),
                    canonical_json_bytes(imp1_checkpoint_extension(restored)),
                )

    def test_nonpass_retry_uses_imp1_ledger_event_and_preserves_boundary(self) -> None:
        items = _prepare_c1r1(rhs=_quarter_impulse)
        prepared, state, tx, tracers, ledger, coordinates = items
        self.assertFalse(prepared.assessment.admission_passed)
        before = (
            c1r1._imp1_reference_tracer_snapshot(tracers),
            tx.state,
            tx.causal_state,
        )
        written = []
        with self.assertRaises(c1r1.C1R1TemporalRetryRequired) as caught:
            c1r1.require_c1r1_admission(
                prepared,
                **_current(*items),
                durable_rejection_sink=lambda record: written.append(
                    canonical_json_bytes(record)
                ),
            )
        self.assertEqual(len(written), 1)
        self.assertEqual(
            before,
            (c1r1._imp1_reference_tracer_snapshot(tracers), tx.state, tx.causal_state),
        )
        self.assertEqual(ledger.current_macro_step_temporal_retry_count, 0)
        retry = caught.exception
        self.assertEqual(retry.evidence["event_type"], IMP1_REJECTION_EVENT_TYPE)
        self.assertEqual(
            retry.evidence["assessment_sha256"], prepared.assessment.assessment_sha256
        )
        self.assertEqual(retry.evidence["preparation_sha256"], prepared.receipt_sha256)
        successor = c1r1.replan_c1r1_retry(prepared, retry)
        self.assertEqual(successor.plan.macro_width, prepared.plan.macro_width / 2)
        next_prepared = c1r1.prepare_c1r1_retry_runtime(
            successor,
            state=state,
            rhs=_quarter_impulse,
            projector=None,
            transaction=tx,
            tracers=tracers,
            coordinates=coordinates,
        )
        self.assertTrue(next_prepared.assessment.admission_passed)

    def test_source_failure_and_stale_guard_preserve_real_state(self) -> None:
        state, tx, tracers, ledger, coordinates = _objects()
        before = (
            c1r1._imp1_reference_tracer_snapshot(tracers),
            tx.state,
            tx.causal_state,
            canonical_json_bytes(imp1_checkpoint_extension(ledger)),
        )

        def bad_rhs(time, current):
            if time > START:
                raise ArithmeticError("synthetic source failure")
            return _rhs(time, current)

        with self.assertRaises(c1r1.C1R1RefinementPathStop) as caught:
            c1r1.prepare_c1r1_initial_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=0.03125,
                state=state,
                rhs=bad_rhs,
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        self.assertFalse(caught.exception.temporal_retry)
        self.assertEqual(
            before,
            (
                c1r1._imp1_reference_tracer_snapshot(tracers),
                tx.state,
                tx.causal_state,
                canonical_json_bytes(imp1_checkpoint_extension(ledger)),
            ),
        )

        items = _prepare_c1r1()
        prepared, state, tx, tracers, ledger, coordinates = items
        tx.thresholds = replace(tx.thresholds, source_residual_maximum=1.0e-6)
        with self.assertRaises(ValueError):
            c1r1.commit_c1r1_runtime(prepared, **_current(*items))
        self.assertEqual(ledger.accepted_macro_step_count, 0)

    def test_lattice_and_row_stops_precede_source(self) -> None:
        class Poison:
            def __getattr__(self, name):
                raise AssertionError(f"premature runtime hook: {name}")

        with self.assertRaises(c1r1.C1R1CoordinateLatticeStop):
            c1r1.prepare_c1r1_initial_runtime(
                method=PRIMARY_METHOD,
                time=np.nextafter(2.0, -np.inf),
                event_target=2.0,
                requested_cap=0.125,
                state=Poison(),
                rhs=Poison(),
                projector=Poison(),
                transaction=Poison(),
                tracers=Poison(),
                coordinates=Poison(),
                temporal_ledger=Poison(),
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        state, tx, tracers, ledger, coordinates = _objects()
        calls = []
        with self.assertRaises(c1r1.C1R1RuntimeResourceStop):
            c1r1.prepare_c1r1_initial_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=0.03125,
                state=state,
                rhs=lambda *args: calls.append(args),
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
                limits=replace(IMP1EnclosureLimits(), maximum_owned_rows=1),
            )
        self.assertEqual(calls, [])

    def test_saturated_endpoint_is_not_replaced_by_corrected_path(self) -> None:
        state, tx, tracers, ledger, coordinates = _objects(saturated=True)

        def constant_rhs(time, current):
            du = np.zeros_like(current.u)
            du[:, 4] = 1.0
            return EvolutionRHS(
                du, np.zeros_like(du), np.zeros_like(du), _diagnostics()
            )

        prepared = c1r1.prepare_c1r1_initial_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=START + 0.25,
            requested_cap=0.25,
            state=state,
            rhs=constant_rhs,
            projector=None,
            transaction=tx,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        items = (prepared, state, tx, tracers, ledger, coordinates)
        index = TDG6_COMPLETE_STATE_CHANNELS.index("u:phi")
        self.assertEqual(prepared.assessment.public_fine_debits[index], Fraction(1, 4))
        committed = c1r1.commit_c1r1_runtime(prepared, **_current(*items))
        self.assertTrue(np.array_equal(committed.state.u[:, 4], state.u[:, 4]))
        self.assertFalse(committed.corrected_endpoint_adopted)


if __name__ == "__main__":
    unittest.main()
