from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal2_source_diagnosis import (  # noqa: E402
    diagnose_gr0_grid_accelerations,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    ADM_CENTER_PARITIES,
    Q_CENTER_PARITIES,
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    SBPFirstDerivative,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


class FGCCAL2SourceDiagnosisTests(unittest.TestCase):
    def test_diagnostic_preserves_raw_root_evidence_on_a_real_slice(self) -> None:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(chi_amplitude=2.5),
            point_count=129,
            outer_radius=128.0,
            constraint_method="RK4",
            diagnostic_spatial_order=4,
        )
        derivative = SBPFirstDerivative(initial.grid, 4)
        state = initial.state
        p_r = derivative.differentiate(
            state.p, center_parities=ADM_CENTER_PARITIES
        )
        q_r = derivative.differentiate(
            state.q, center_parities=Q_CENTER_PARITIES
        )
        diagnosis = diagnose_gr0_grid_accelerations(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            initial.grid.coordinates[1:],
        )
        self.assertTrue(diagnosis.raw_gate_passed)
        self.assertLess(diagnosis.residual_infinity, diagnosis.raw_tolerance)
        self.assertLess(diagnosis.kinetic_condition_infinity_maximum, 1.0e10)
        self.assertGreaterEqual(len(diagnosis.iterations), 1)

    def test_invalid_limits_fail_closed(self) -> None:
        arrays = np.zeros((2, 6), dtype=np.float64)
        radii = np.array([1.0, 2.0], dtype=np.float64)
        with self.assertRaisesRegex(ValueError, "raw_tolerance"):
            diagnose_gr0_grid_accelerations(
                arrays,
                arrays,
                arrays,
                arrays,
                arrays,
                radii,
                raw_tolerance=0.0,
            )
        with self.assertRaisesRegex(ValueError, "maximum_refinement_iterations"):
            diagnose_gr0_grid_accelerations(
                arrays,
                arrays,
                arrays,
                arrays,
                arrays,
                radii,
                maximum_refinement_iterations=17,
            )


if __name__ == "__main__":
    unittest.main()
