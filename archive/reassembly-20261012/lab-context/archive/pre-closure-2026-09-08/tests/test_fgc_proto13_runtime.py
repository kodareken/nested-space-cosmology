from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto13_runtime import (  # noqa: E402
    PROTO13_COMPARATOR_POINT_COUNTS,
    PROTO13_PRIMARY_POINT_COUNTS,
    proto13_method_owned_trapped_assessment,
)


@dataclass(frozen=True)
class Event:
    method: str
    point_counts: tuple[int, int, int]
    coordinate_time: float
    admission_passed: bool = True


def state(count: int, normal: float) -> EvolutionState:
    radii = np.linspace(0.0, 128.0, count)
    u = np.zeros((count, 6), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p[:, 3] = normal * radii
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


class FGCProto13RuntimeTests(unittest.TestCase):
    def assessment(self, *, comparator_shift: float = 0.0, admission: bool = True):
        primary_grids = [
            UniformRadialGrid(0.0, 128.0, count)
            for count in PROTO13_PRIMARY_POINT_COUNTS
        ]
        comparator_grids = [
            UniformRadialGrid(0.0, 128.0, count)
            for count in PROTO13_COMPARATOR_POINT_COUNTS
        ]
        return proto13_method_owned_trapped_assessment(
            primary_states=[state(count, -2.0) for count in PROTO13_PRIMARY_POINT_COUNTS],
            primary_grids=primary_grids,
            comparator_states=[
                state(count, -2.0 + comparator_shift)
                for count in PROTO13_COMPARATOR_POINT_COUNTS
            ],
            comparator_grids=comparator_grids,
            primary_common_event=Event(
                "RK4", PROTO13_PRIMARY_POINT_COUNTS, 1.5, admission
            ),
            comparator_common_event=Event(
                "SSPRK3", PROTO13_COMPARATOR_POINT_COUNTS, 1.5, admission
            ),
        )

    def test_exact_common_node_interval_pass(self) -> None:
        observed = self.assessment()
        self.assertEqual(observed.common_physical_point_count, 2049)
        self.assertEqual(len(set(observed.primary_trapped_scores)), 1)
        self.assertEqual(
            observed.primary_trapped_scores, observed.comparator_trapped_scores
        )
        self.assertGreater(observed.primary_trapped_scores[-1], 0.0)
        self.assertTrue(observed.cross_method_interval_overlap)
        self.assertTrue(observed.trapped_sign_passed)

    def test_disjoint_method_intervals_fail(self) -> None:
        observed = self.assessment(comparator_shift=0.25)
        self.assertGreater(observed.cross_method_interval_gap, 0.0)
        self.assertFalse(observed.cross_method_interval_overlap)
        self.assertFalse(observed.trapped_sign_passed)

    def test_common_event_admission_remains_a_veto(self) -> None:
        observed = self.assessment(admission=False)
        self.assertFalse(observed.both_method_owned_common_event_admissions_passed)
        self.assertFalse(observed.trapped_sign_passed)

    def test_wrong_ladder_and_domain_fail_closed(self) -> None:
        primary_grids = [
            UniformRadialGrid(0.0, 128.0, count)
            for count in PROTO13_PRIMARY_POINT_COUNTS
        ]
        comparator_grids = [
            UniformRadialGrid(0.0, 128.0, count)
            for count in PROTO13_COMPARATOR_POINT_COUNTS
        ]
        comparator_grids[-1] = UniformRadialGrid(0.0, 64.0, 16385)
        with self.assertRaisesRegex(ValueError, "radial domain"):
            proto13_method_owned_trapped_assessment(
                primary_states=[state(count, -2.0) for count in PROTO13_PRIMARY_POINT_COUNTS],
                primary_grids=primary_grids,
                comparator_states=[state(count, -2.0) for count in PROTO13_COMPARATOR_POINT_COUNTS],
                comparator_grids=comparator_grids,
                primary_common_event=Event("RK4", PROTO13_PRIMARY_POINT_COUNTS, 1.5),
                comparator_common_event=Event("SSPRK3", PROTO13_COMPARATOR_POINT_COUNTS, 1.5),
            )


if __name__ == "__main__":
    unittest.main()
