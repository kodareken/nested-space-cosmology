from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (
    OwnedConstraintAdmission,
    OwnedConstraintSnapshot,
)
from recursive_horizons.fgc.evolution.numerical_engine import (
    EvolutionState,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto4_admission import (
    PROTO4_CONSTRAINT_ORDER,
    ConstraintNorms,
)
from recursive_horizons.fgc.evolution.rsp2_constraint_runtime import (
    RSP2_POINT_COUNTS,
    rsp2_constraint_assessment,
)
from recursive_horizons.fgc.evolution.rsp2_checkpoint import (
    RSP2_CAL9_STAGE_COUNTS,
    RSP2_CAL9_STATE_HASHES,
)


TARGET = PROTO4_CONSTRAINT_ORDER.index("radial_momentum")


def _norms(point_count: int, target: float) -> ConstraintNorms:
    values = np.full(len(PROTO4_CONSTRAINT_ORDER), target / 4.0)
    values[TARGET] = target
    return ConstraintNorms(
        sample_count=point_count,
        component_names=PROTO4_CONSTRAINT_ORDER,
        component_infinity=values,
        global_infinity=float(np.max(values)),
        minimum_denominator=1.0,
        maximum_denominator=2.0,
    )


def _snapshot(point_count: int, raw: float, enclosure: float) -> OwnedConstraintSnapshot:
    owned_count = point_count - 5
    effective = max(raw - enclosure, 0.0)
    return OwnedConstraintSnapshot(
        point_count=point_count,
        coordinate_time=23.0 / 16.0,
        diagnostic_spatial_order=2,
        accepted_stage_count=10,
        fixed_outer_rows=4,
        derivative_stencil_reach_rows=1,
        excluded_outer_rows=5,
        last_owned_radius=127.0,
        full_domain_norms=_norms(point_count, raw),
        owned_domain_norms=_norms(owned_count, raw),
        effective_owned_norms_for_order_only=_norms(owned_count, effective),
        accumulated_roundoff_enclosure=enclosure,
        full_component_maximum_radii=np.ones(len(PROTO4_CONSTRAINT_ORDER)),
        owned_component_maximum_radii=np.ones(len(PROTO4_CONSTRAINT_ORDER)),
    )


def _admission(order: float, *, passed: bool = True) -> OwnedConstraintAdmission:
    orders = {name: 2.0 for name in PROTO4_CONSTRAINT_ORDER}
    orders["radial_momentum"] = order
    statuses = {name: "finite_order" for name in PROTO4_CONSTRAINT_ORDER}
    zeros = {name: (0.01, 0.003, 0.001) for name in PROTO4_CONSTRAINT_ORDER}
    radii = {name: (1.0, 1.0, 1.0) for name in PROTO4_CONSTRAINT_ORDER}
    return OwnedConstraintAdmission(
        method="SSPRK3",
        point_counts=RSP2_POINT_COUNTS,
        coordinate_time=23.0 / 16.0,
        excluded_outer_rows=(5, 5, 5),
        last_owned_radii=(127.0, 127.0, 127.0),
        accepted_stage_counts=(1, 2, 3),
        accumulated_roundoff_enclosures=(1e-12, 2e-12, 3e-12),
        full_domain_component_infinity=zeros,
        owned_domain_component_infinity=zeros,
        full_component_maximum_radii=radii,
        owned_component_maximum_radii=radii,
        owned_raw_global_norms=(0.01, 0.003, 0.001),
        component_finest_pair_orders=orders,
        component_status=statuses,
        minimum_finite_finest_pair_order=min(orders.values()),
        coarsest_guard_passed=True,
        finest_guard_passed=True,
        monotone_refinement_passed=True,
        finest_pair_order_passed=passed,
        common_event_alignment_passed=True,
        admission_passed=passed,
    )


class RSP2ConstraintRuntimeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.states = [
            EvolutionState(
                np.zeros((count, 6)),
                np.zeros((count, 6)),
                np.zeros((count, 6)),
            )
            for count in RSP2_POINT_COUNTS
        ]
        self.grids = [UniformRadialGrid(0.0, 128.0, count) for count in RSP2_POINT_COUNTS]

    def _evaluate(self, snapshots, admission):
        with patch(
            "recursive_horizons.fgc.evolution.rsp2_constraint_runtime.owned_semidiscrete_constraint_snapshot",
            side_effect=snapshots,
        ), patch(
            "recursive_horizons.fgc.evolution.rsp2_constraint_runtime.owned_constraint_admission",
            return_value=admission,
        ):
            return rsp2_constraint_assessment(
                self.states,
                self.grids,
                accepted_stage_counts=(1, 2, 3),
                coordinate_time=23.0 / 16.0,
                enforce_frozen_predecessors=False,
            )

    def test_exact_three_halves_passes(self) -> None:
        snapshots = (
            _snapshot(4097, 0.02, 1e-12),
            _snapshot(8193, 0.008, 2e-12),
            _snapshot(16385, 0.008 / (2.0**1.5), 3e-12),
        )
        result = self._evaluate(snapshots, _admission(1.5, passed=True))
        self.assertEqual(result.classification, "completed_target_and_complete_constraint_pass")
        self.assertTrue(result.target_order_passed)
        self.assertTrue(result.target_preasymptotic_on_tested_ladder)
        self.assertFalse(result.target_persistent_through_tested_ladder)

    def test_one_bit_below_three_halves_persists(self) -> None:
        below = np.nextafter(1.5, -np.inf)
        snapshots = (
            _snapshot(4097, 0.02, 1e-12),
            _snapshot(8193, 0.008, 2e-12),
            _snapshot(16385, 0.0029, 3e-12),
        )
        result = self._evaluate(snapshots, _admission(below, passed=False))
        self.assertEqual(result.classification, "completed_target_order_persistent_below_threshold")
        self.assertFalse(result.target_order_passed)
        self.assertTrue(result.target_persistent_through_tested_ladder)

    def test_target_clear_with_other_component_obstruction_is_separate(self) -> None:
        snapshots = (
            _snapshot(4097, 0.02, 1e-12),
            _snapshot(8193, 0.008, 2e-12),
            _snapshot(16385, 0.002, 3e-12),
        )
        admission = _admission(2.0, passed=False)
        orders = dict(admission.component_finest_pair_orders)
        orders["gauge_r"] = 1.0
        admission = replace(
            admission,
            component_finest_pair_orders=orders,
            minimum_finite_finest_pair_order=1.0,
        )
        result = self._evaluate(snapshots, admission)
        self.assertEqual(result.classification, "completed_target_cleared_other_constraint_obstruction")
        self.assertEqual(result.non_target_order_obstructions, ("gauge_r",))

    def test_nonmonotone_target_is_inconclusive_not_persistent(self) -> None:
        snapshots = (
            _snapshot(4097, 0.02, 1e-12),
            _snapshot(8193, 0.008, 2e-12),
            _snapshot(16385, 0.009, 3e-12),
        )
        result = self._evaluate(snapshots, _admission(1.0, passed=False))
        self.assertEqual(result.classification, "completed_target_inconclusive_premise_failure")
        self.assertFalse(result.target_persistent_through_tested_ladder)

    def test_wrong_grid_ladder_fails_closed(self) -> None:
        grids = list(self.grids)
        grids[-1] = UniformRadialGrid(0.0, 128.0, 32769)
        with self.assertRaisesRegex(ValueError, "4097 -> 8193 -> 16385"):
            rsp2_constraint_assessment(
                self.states,
                grids,
                accepted_stage_counts=(1, 2, 3),
                coordinate_time=23.0 / 16.0,
            )

    def test_wrong_or_nonfinite_time_fails_closed(self) -> None:
        for value in (1.5, float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                rsp2_constraint_assessment(
                    self.states,
                    self.grids,
                    accepted_stage_counts=(*RSP2_CAL9_STAGE_COUNTS, 3),
                    coordinate_time=value,
                )

    def test_predecessor_stage_count_mutation_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "stage counts differ"):
            rsp2_constraint_assessment(
                self.states,
                self.grids,
                accepted_stage_counts=(RSP2_CAL9_STAGE_COUNTS[0] + 1, 3948, 3),
                coordinate_time=23.0 / 16.0,
            )

    def test_predecessor_state_mutation_fails_closed(self) -> None:
        with patch(
            "recursive_horizons.fgc.evolution.rsp2_constraint_runtime.array_content_sha256",
            side_effect=(RSP2_CAL9_STATE_HASHES[0], "0" * 64),
        ), self.assertRaisesRegex(ValueError, "predecessor states differ"):
            rsp2_constraint_assessment(
                self.states,
                self.grids,
                accepted_stage_counts=(*RSP2_CAL9_STAGE_COUNTS, 3),
                coordinate_time=23.0 / 16.0,
            )

    def test_roundoff_subtraction_changes_only_order_evidence(self) -> None:
        snapshots = (
            _snapshot(4097, 0.02, 0.001),
            _snapshot(8193, 0.008, 0.002),
            _snapshot(16385, 0.003, 0.0025),
        )
        result = self._evaluate(snapshots, _admission(2.0, passed=True))
        self.assertEqual(result.target_raw_owned_norms, (0.02, 0.008, 0.003))
        self.assertEqual(
            result.target_accumulated_roundoff_enclosures,
            (0.001, 0.002, 0.0025),
        )
        self.assertEqual(
            result.target_effective_norms_for_order_only,
            (0.019, 0.006, 0.0005),
        )
        self.assertEqual(
            result.constraint_admission.owned_raw_global_norms,
            (0.01, 0.003, 0.001),
        )

    def test_wrong_method_or_component_record_fails_closed(self) -> None:
        snapshots = (
            _snapshot(4097, 0.02, 1e-12),
            _snapshot(8193, 0.008, 2e-12),
            _snapshot(16385, 0.002, 3e-12),
        )
        result = self._evaluate(snapshots, _admission(2.0, passed=True))
        with self.assertRaisesRegex(ValueError, "method differs"):
            replace(result, method="RK4")
        with self.assertRaisesRegex(ValueError, "target component differs"):
            replace(result, target_component="gauge_r")


if __name__ == "__main__":
    unittest.main()
