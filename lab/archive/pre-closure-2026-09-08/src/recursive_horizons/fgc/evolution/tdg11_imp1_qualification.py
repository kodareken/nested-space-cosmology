"""No-trajectory synthetic qualification fixtures for the new IMP1 kernel.

This module is a qualification owner, never imported by the production
runtime. Synthetic health diagnostics do not certify any physical source.
Independent reconstruction/localization is used only on the reference side.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
from time import perf_counter
from typing import Callable

import numpy as np

from recursive_horizons.evidence_io import canonical_json_bytes

from .boundary_domain import BoundaryGeometry, CausalBudgetState
from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeStageTransaction, GR0UniversalThresholds
from .proto19_gr0_static_factory import NormalFlowTracers
from .tdg11_imp1_ledger import (
    imp1_checkpoint_extension,
    restore_imp1_checkpoint_extension,
    seed_imp1_ledger,
)
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS
from . import tdg6_temporal_admission_runtime as tdg6


SYNTHETIC_START = 23.0 / 16.0
SYNTHETIC_TARGET = 24.0 / 16.0


def synthetic_diagnostics() -> dict[str, object]:
    """Declared guard fixture; not GR-0/FGC source or constraint evidence."""

    return {
        "source_residual_infinity": 0.0,
        "source_refinement_iterations": 0,
        "source_residual_decreased_monotonically": True,
        "kinetic_condition_infinity": 1.0,
        "coordinate_speed_upper": 1.0,
        "minimum_lapse": 1.0,
        "minimum_radial_metric": 1.0,
        "minimum_areal_radius_away_from_center": 0.125,
    }


def synthetic_exponential_rhs(time: float, state: EvolutionState) -> EvolutionRHS:
    return EvolutionRHS(state.u, state.p, state.q, synthetic_diagnostics())


def synthetic_saturation_rhs(time: float, state: EvolutionState) -> EvolutionRHS:
    du = np.zeros_like(state.u)
    du[:, 4] = 1.0
    return EvolutionRHS(
        du, np.zeros_like(du), np.zeros_like(du), synthetic_diagnostics()
    )


def make_synthetic_fixture(
    method: str = PRIMARY_METHOD, *, saturated: bool = False, point_count: int = 9
):
    """Independent, explicit small fixture; opens no file or stored state."""

    if method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
        raise ValueError("unknown synthetic IMP1 method")
    if type(point_count) is not int or point_count not in (9, 17, 33, 129, 2049):
        raise ValueError("point count differs from the prospective synthetic controls")
    coordinates = np.linspace(0.0, 1.0, point_count)
    u = np.zeros((point_count, 6))
    u[:, 0], u[:, 1], u[:, 2], u[:, 3] = 1.0, 0.025, 1.0, coordinates
    u[:, 4] = 2.0**52 if saturated else 0.125
    u[:, 5] = 0.0625
    p, q = np.zeros_like(u), np.zeros_like(u)
    q[:, 3] = 1.0
    state = EvolutionState(u, p, q)
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(
            accepted_time=SYNTHETIC_START, previous_speed_upper=1.0
        ),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )
    tracers = NormalFlowTracers.create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=coordinates,
        cutoff=1.0,
        outer_radius=128.0,
    )
    inherited = replace(
        tdg6.TDG6TemporalLedger.zero(initial_time=SYNTHETIC_START),
        accumulated_debit_vector=tuple(float(index + 1) / 1024 for index in range(18)),
    )
    ledger = seed_imp1_ledger(
        inherited,
        method=method,
        state_sha256=array_content_sha256(state.u, state.p, state.q),
        step_index=0,
        transaction_serial=0,
        origin_receipt_sha256="1" * 64,
    )
    return state, transaction, tracers, ledger, coordinates


def rational_pair(value: Fraction) -> dict[str, str]:
    """Only for the bounded small synthetic report, never a permissive parser."""

    if (
        type(value) is not Fraction
        or max(value.numerator.bit_length(), value.denominator.bit_length()) > 4096
    ):
        raise ValueError("synthetic report rational is outside its small control bound")
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _prove(condition: bool, label: str) -> bool:
    if condition is not True:
        raise ValueError(f"IMP1 synthetic qualification did not prove {label}")
    return True


def _prepare_fixture(method: str, *, point_count: int = 9, saturated: bool = False):
    from .tdg11_imp1_runtime import prepare_imp1_initial_runtime

    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(
        method, point_count=point_count, saturated=saturated
    )
    prepared = prepare_imp1_initial_runtime(
        method=method,
        time=SYNTHETIC_START,
        event_target=SYNTHETIC_START + 0.25 if saturated else SYNTHETIC_TARGET,
        requested_cap=0.25 if saturated else 1.0 / 32,
        state=state,
        rhs=synthetic_saturation_rhs if saturated else synthetic_exponential_rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger, coordinates


def _current(items) -> dict[str, object]:
    prepared, state, transaction, tracers, ledger, coordinates = items
    return {
        "transaction": transaction,
        "tracers": tracers,
        "temporal_ledger": ledger,
        "current_time": prepared.initial_time,
        "current_state": state,
        "current_step_index": prepared.previous_step_index,
        "current_transaction_serial": prepared.previous_transaction_serial,
        "coordinates": coordinates,
    }


def _independent_family(prepared):
    from .tdg11_msel1_pref1_reconstruction import (
        ORIGINAL_ARITHMETIC_ID,
        validate_recorded_family_independently,
    )

    return validate_recorded_family_independently(
        prepared.recorded_family.paths,
        method=prepared.method,
        arithmetic_id=ORIGINAL_ARITHMETIC_ID,
        owned_row_count=prepared.recorded_family.owned_row_count,
    )


def _reference_bounds(rows):
    from .tdg11_msel1_pref1_localization import assess_complete_c_independently

    return assess_complete_c_independently(
        rows,
        expected_row_count=len(rows),
        maximum_candidates_D01=8 * len(rows),
        maximum_candidates_D12=16 * len(rows),
        refinement_depth=160,
    )


def _overlap(left_lower, left_upper, right_lower, right_upper) -> bool:
    return max(left_lower, right_lower) <= min(left_upper, right_upper)


def _independent_reference_check(prepared) -> int:
    from .tdg11_msel1_reconstruction import reconstruct_channel
    from .tdg11_msel1_pref1_reconstruction import reconstruct_channel_independently

    reference_family = _independent_family(prepared)
    _prove(
        tuple(record.channel for record in prepared.assessment.channels)
        == TDG6_COMPLETE_STATE_CHANNELS,
        "reference channel coverage",
    )
    for channel, record in zip(
        TDG6_COMPLETE_STATE_CHANNELS, prepared.assessment.channels, strict=True
    ):
        independent = reconstruct_channel_independently(reference_family, channel)
        produced = reconstruct_channel(prepared.recorded_family, channel)
        _prove(
            independent.raw_rows == produced.raw_rows
            and independent.corrected_rows == produced.corrected_rows
            and independent.accumulation_bounds == produced.accumulation_bounds,
            "independent raw/corrected reconstruction and accumulation",
        )
        corrected = _reference_bounds(independent.corrected_rows)
        raw = _reference_bounds(independent.raw_rows)
        actual = record.corrected
        _prove(
            actual.accumulation_bounds == independent.accumulation_bounds,
            "complete noncancelling debit",
        )
        _prove(
            actual.public_fine_debit
            == actual.bernstein.d12_upper + independent.accumulation_bounds[2]
            and actual.gate_debit
            == actual.gate_d12_upper + independent.accumulation_bounds[2]
            and actual.extra_enclosure_debit
            == actual.bernstein.d12_upper - actual.gate_d12_upper
            and actual.extra_enclosure_debit >= 0,
            "conservative public/gate/debit separation",
        )
        for bounds, reference in (
            ((actual.gate_d01_lower, actual.gate_d01_upper), corrected.d01),
            ((actual.gate_d12_lower, actual.gate_d12_upper), corrected.d12),
            ((record.raw.d01_lower, record.raw.d01_upper), raw.d01),
            ((record.raw.d12_lower, record.raw.d12_upper), raw.d12),
        ):
            _prove(
                _overlap(*bounds, reference.lower, reference.upper),
                "independent complete-C enclosure overlap",
            )
        _prove(
            actual.decision.admission_passed is True
            and corrected.decision.admission_passed is True,
            "independent synthetic order/zero agreement",
        )
    return len(TDG6_COMPLETE_STATE_CHANNELS)


def _method_control(method: str) -> dict[str, object]:
    from .tdg11_imp1_runtime import commit_imp1_runtime, prepare_imp1_initial_runtime

    items = _prepare_fixture(method)
    prepared, _, transaction, tracers, ledger, coordinates = items
    reference_count = _independent_reference_check(prepared)
    before_count = transaction.state.accepted_stage_count
    committed = commit_imp1_runtime(prepared, **_current(items))
    checkpoint = imp1_checkpoint_extension(committed.temporal_ledger)
    restored = restore_imp1_checkpoint_extension(checkpoint)
    next_prepared = prepare_imp1_initial_runtime(
        method=method,
        time=committed.time,
        event_target=SYNTHETIC_TARGET,
        requested_cap=1.0 / 32,
        state=committed.state,
        rhs=synthetic_exponential_rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=restored,
        previous_step_index=committed.step_index,
        previous_transaction_serial=committed.transaction_serial,
    )
    return {
        "method": method,
        "point_count": 9,
        "owned_row_count": prepared.recorded_family.owned_row_count,
        "channels": [record.channel for record in prepared.assessment.channels],
        "shadow_stage_records": prepared.stage_record_count,
        "next_family_shadow_stage_records": next_prepared.stage_record_count,
        "total_shadow_proposals": sum(
            len(path.attempts)
            for instance in (prepared, next_prepared)
            for path in (instance.outer, instance.medium, instance.fine)
        ),
        "committed_proposals": committed.step_index - prepared.previous_step_index,
        "committed_stage_records": transaction.state.accepted_stage_count
        - before_count,
        "independent_reference_channels": reference_count,
        "endpoint_sha256": array_content_sha256(
            committed.state.u, committed.state.p, committed.state.q
        ),
        "assessment_sha256": prepared.assessment.assessment_sha256,
        "checkpoint_sha256": sha256(canonical_json_bytes(checkpoint)).hexdigest(),
        "admission_passed": _prove(
            prepared.assessment.admission_passed is True, "all-channel admission"
        ),
        "inherited_snapshot_preserved": _prove(
            restored.inherited_snapshot == ledger.inherited_snapshot,
            "inherited snapshot",
        ),
        "exact_debit_added": _prove(
            restored.accumulated_debit_vector
            == tuple(
                a + b
                for a, b in zip(
                    ledger.accumulated_debit_vector,
                    prepared.assessment.public_fine_debits,
                    strict=True,
                )
            ),
            "exact debit addition",
        ),
        "fine_endpoint_adopted": _prove(
            committed.state is prepared.fine.final_accepted.state
            and committed.corrected_endpoint_adopted is False,
            "actual fine-only adoption",
        ),
        "checkpoint_roundtrip_identical": _prove(
            canonical_json_bytes(checkpoint)
            == canonical_json_bytes(imp1_checkpoint_extension(restored)),
            "checkpoint roundtrip",
        ),
        "next_macro_prepared": _prove(
            next_prepared.initial_time == committed.time
            and next_prepared.assessment.admission_passed is True,
            "next macro from restored ledger",
        ),
        "independent_reference_agreement": _prove(
            reference_count == 18, "all independent channels"
        ),
    }


def _saturation_control() -> dict[str, object]:
    from .tdg11_imp1_runtime import commit_imp1_runtime
    from .tdg11_msel1_pref1_reconstruction import reconstruct_channel_independently

    items = _prepare_fixture(PRIMARY_METHOD, saturated=True)
    prepared, state, _, _, _, _ = items
    reference = reconstruct_channel_independently(
        _independent_family(prepared), "u:phi"
    )
    raw = _reference_bounds(reference.raw_rows)
    corrected = _reference_bounds(reference.corrected_rows)
    for bound in (raw.d01, raw.d12, corrected.d01, corrected.d12):
        _prove(bound.lower == bound.upper, "exact synthetic saturation maximum")
    debit = prepared.assessment.public_fine_debits[
        TDG6_COMPLETE_STATE_CHANNELS.index("u:phi")
    ]
    _prove(debit == Fraction(1, 4), "retained saturation debit")
    committed = commit_imp1_runtime(prepared, **_current(items))
    return {
        "method": PRIMARY_METHOD,
        "channel": "u:phi",
        "initial_value_hex": float(state.u[1, 4]).hex(),
        "width_hex": prepared.plan.macro_width.hex(),
        "actual_endpoint_unchanged": _prove(
            bool(np.array_equal(committed.state.u[:, 4], state.u[:, 4])),
            "saturated actual endpoint",
        ),
        "raw_D01": rational_pair(raw.d01.upper),
        "raw_D12": rational_pair(raw.d12.upper),
        "corrected_D01": rational_pair(corrected.d01.upper),
        "corrected_D12": rational_pair(corrected.d12.upper),
        "public_fine_debit": rational_pair(debit),
        "raw_failure_not_relabeled": _prove(
            raw.sufficient_contraction_failure is True
            and corrected.decision.classification == "exact_zero",
            "separate raw and corrected semantics",
        ),
    }


def run_qualification(
    *, progress: Callable[[dict], None] | None = None
) -> dict[str, object]:
    """Run only the prospectively declared synthetic controls, serially.

    Timings go to optional progress output and are never a scientific gate
    or part of the deterministic compact commitment.
    """

    from .tdg11_imp1_enclosure import IMP1EnclosureLimits

    _prove(
        IMP1EnclosureLimits().as_mapping()
        == {
            "maximum_owned_rows": 2044,
            "maximum_rational_bits": 32768,
            "maximum_subdivision_depth": 12,
            "maximum_subdivisions_per_pair": 32768,
            "fallback_maximum_candidates_D01": 16352,
            "fallback_maximum_candidates_D12": 32704,
            "fallback_refinement_depth": 160,
        },
        "prospective enclosure resource settings",
    )
    method_controls = []
    for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
        started = perf_counter()
        method_controls.append(_method_control(method))
        if progress is not None:
            progress(
                {
                    "control": "independent_method_reference",
                    "method": method,
                    "elapsed_seconds": perf_counter() - started,
                }
            )
    saturation = _saturation_control()
    sizes = []
    for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
        for count in (129, 2049):
            started = perf_counter()
            prepared = _prepare_fixture(method, point_count=count)[0]
            sizes.append(
                {
                    "method": method,
                    "point_count": count,
                    "owned_row_count": prepared.recorded_family.owned_row_count,
                    "channel_count": len(prepared.assessment.channels),
                    "admission_passed": _prove(
                        prepared.assessment.admission_passed is True,
                        "all-channel size control",
                    ),
                    "assessment_sha256": prepared.assessment.assessment_sha256,
                    "shadow_stage_records": prepared.stage_record_count,
                    "physical_source_tested": False,
                }
            )
            if progress is not None:
                progress(
                    {
                        "control": "synthetic_size",
                        "method": method,
                        "point_count": count,
                        "elapsed_seconds": perf_counter() - started,
                        "physical_source_tested": False,
                    }
                )
    return {
        "method_controls": method_controls,
        "saturation_control": saturation,
        "size_controls": sizes,
    }


__all__ = [
    "SYNTHETIC_START",
    "SYNTHETIC_TARGET",
    "make_synthetic_fixture",
    "synthetic_diagnostics",
    "synthetic_exponential_rhs",
    "synthetic_saturation_rhs",
    "run_qualification",
]
