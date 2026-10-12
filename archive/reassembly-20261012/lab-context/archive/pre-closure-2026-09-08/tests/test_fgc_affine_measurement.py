from __future__ import annotations

from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.affine_measurement import (  # noqa: E402
    analyze_affine_null,
    polynomial_flrw_history,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)


class FGCAffineMeasurementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.history = polynomial_flrw_history(
            time_count=65,
            radius_count=129,
            final_time=0.25,
            maximum_radius=4.0,
        )

    def test_both_integrators_recover_null_affine_raychaudhuri_identity(self) -> None:
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            with self.subTest(method=method):
                result = analyze_affine_null(
                    self.history,
                    method=method,
                    launch_time=0.0,
                    launch_radius=1.0,
                    final_time=0.25,
                    step_size=1.0 / 2048.0,
                )
                self.assertEqual(result.initial_normalization_residual, 0.0)
                self.assertLessEqual(result.maximum_null_residual, 1.0e-12)
                self.assertLessEqual(result.maximum_affine_residual, 1.0e-10)
                self.assertLessEqual(result.maximum_route_disagreement, 1.0e-8)
                self.assertTrue((result.k_t > 0.0).all())
                self.assertTrue((result.affine_parameters[1:] > result.affine_parameters[:-1]).all())

    def test_history_rejects_trajectory_outside_interpolation_domain(self) -> None:
        with self.assertRaisesRegex(ValueError, "left the stored history"):
            analyze_affine_null(
                self.history,
                method=PRIMARY_METHOD,
                launch_time=0.0,
                launch_radius=3.99,
                final_time=0.25,
                step_size=1.0 / 128.0,
            )


if __name__ == "__main__":
    unittest.main()
