"""Supplementary end-to-end retry-limit controls; no campaign or stored data."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    make_synthetic_fixture,
    synthetic_diagnostics,
    _current,
)


def _jump_rhs(time, state):
    """Intentionally nonsmooth fixed source: repeated scientific nonadmission."""

    du = np.zeros_like(state.u)
    du[:, 4] = 1.0 if time > START else 0.0
    return EvolutionRHS(
        du, np.zeros_like(du), np.zeros_like(du), synthetic_diagnostics()
    )


def _initial(cap):
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture()
    # The large-cap counter control is a vector ODE fixture, not physical CFL
    # qualification. Zero shift keeps synthetic tracers inside their domain.
    u = state.u.copy()
    u[:, 1] = 0.0
    state = EvolutionState(u, state.p, state.q)
    state_hash = array_content_sha256(state.u, state.p, state.q)
    ledger = replace(
        ledger, origin_state_sha256=state_hash, current_state_sha256=state_hash
    )
    transaction.grid_spacing = 32.0
    tracers = type(tracers).create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=coordinates,
        cutoff=1.0,
        outer_radius=128.0,
    )
    prepared = runtime.prepare_imp1_initial_runtime(
        method=PRIMARY_METHOD,
        time=START,
        event_target=32.0,
        requested_cap=cap,
        state=state,
        rhs=_jump_rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger, coordinates


class IMP1RetryLimitTests(unittest.TestCase):
    def test_exact_minimum_width_rejection_is_retained_without_retry(self):
        items = _initial(2.0**-30)
        written = []
        with self.assertRaises(runtime.IMP1TemporalRetryExhausted) as caught:
            runtime.require_imp1_admission(
                items[0], **_current(items), durable_rejection_sink=written.append
            )
        self.assertEqual(caught.exception.reason, "minimum_macro_step")
        self.assertEqual(len(written), 1)
        self.assertEqual(
            caught.exception.updated_ledger.current_macro_step_temporal_retry_count, 1
        )
        self.assertEqual(items[2].state.accepted_stage_count, 0)
        self.assertEqual(items[4].current_macro_step_temporal_retry_count, 0)

    def test_thirty_third_rejection_exhausts_count_before_minimum_width(self):
        items = _initial(8.0)
        written = []
        initial_snapshot = runtime._tracer_snapshot(items[3])
        for attempt in range(1, 34):
            prepared, state, transaction, tracers, ledger, coordinates = items
            self.assertFalse(prepared.assessment.admission_passed)
            if attempt == 33:
                with self.assertRaises(runtime.IMP1TemporalRetryExhausted) as caught:
                    runtime.require_imp1_admission(
                        prepared,
                        **_current(items),
                        durable_rejection_sink=written.append,
                    )
                self.assertEqual(caught.exception.reason, "maximum_temporal_retries")
                self.assertEqual(
                    caught.exception.updated_ledger.current_macro_step_temporal_retry_count,
                    33,
                )
                self.assertEqual(
                    len(caught.exception.updated_ledger.serialized_temporal_rejections),
                    33,
                )
                break
            with self.assertRaises(runtime.IMP1TemporalRetryRequired) as caught:
                runtime.require_imp1_admission(
                    prepared, **_current(items), durable_rejection_sink=written.append
                )
            successor = runtime.replan_imp1_retry(prepared, caught.exception)
            next_prepared = runtime.prepare_imp1_retry_runtime(
                successor,
                state=state,
                rhs=_jump_rhs,
                projector=None,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
            )
            items = (
                next_prepared,
                state,
                transaction,
                tracers,
                successor.updated_ledger,
                coordinates,
            )
        self.assertEqual(len(written), 33)
        self.assertEqual(items[2].state.accepted_stage_count, 0)
        self.assertEqual(items[2].causal_state.accepted_time, START)
        self.assertEqual(runtime._tracer_snapshot(items[3]), initial_snapshot)
        self.assertEqual(items[4].accepted_macro_step_count, 0)

    def test_injected_retry_lattice_stop_keeps_acknowledged_rejection(self):
        items = _initial(0.125)
        written = []
        with self.assertRaises(runtime.IMP1TemporalRetryRequired) as caught:
            runtime.require_imp1_admission(
                items[0], **_current(items), durable_rejection_sink=written.append
            )
        retry = caught.exception
        with patch.object(
            runtime,
            "_plan",
            side_effect=runtime.IMP1CoordinateLatticeStop("injected_no_aligned_width"),
        ):
            with self.assertRaises(runtime.IMP1TemporalRetryExhausted) as stopped:
                runtime.replan_imp1_retry(items[0], retry)
        self.assertEqual(stopped.exception.reason, "coordinate_lattice")
        self.assertIs(stopped.exception.updated_ledger, retry.updated_ledger)
        self.assertEqual(len(written), 1)
        self.assertEqual(items[2].state.accepted_stage_count, 0)


if __name__ == "__main__":
    unittest.main()
