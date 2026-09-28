from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
    proto9_gr0_common_event,
    validate_proto9_common_event_assessment,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)


def _flat(count: int, order: int = 4) -> tuple[EvolutionState, UniformRadialGrid]:
    grid = UniformRadialGrid(0.0, 128.0, count)
    u = np.zeros((count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    return EvolutionState(u, p, q), grid


class FGCProto9RuntimeTests(unittest.TestCase):
    def test_exact_new_ladder_delegates_to_unchanged_PROTO8_compositor(self) -> None:
        pairs = [_flat(count) for count in PROTO9_POINT_COUNTS]
        assessed = proto9_gr0_common_event(
            [item[0] for item in pairs],
            [item[1] for item in pairs],
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
        self.assertIs(validate_proto9_common_event_assessment(assessed), assessed)
        self.assertEqual(assessed.point_counts, PROTO9_POINT_COUNTS)
        self.assertEqual(
            assessed.admission_passed,
            assessed.constraint_admission.admission_passed
            and assessed.spatial_spectral_admission.admission_passed,
        )

    def test_old_partial_and_reordered_ladders_fail_closed(self) -> None:
        for counts in (
            (1025, 2049, 4097),
            (2049, 4097, 4097),
            tuple(reversed(PROTO9_POINT_COUNTS)),
        ):
            pairs = [_flat(count) for count in counts]
            with self.assertRaisesRegex(ValueError, "exact point counts"):
                proto9_gr0_common_event(
                    [item[0] for item in pairs],
                    [item[1] for item in pairs],
                    accepted_stage_counts=(0, 0, 0),
                    method="RK4",
                    coordinate_time=0.0,
                )


if __name__ == "__main__":
    unittest.main()
