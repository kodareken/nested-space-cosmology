from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_src2_affine_capture import (  # noqa: E402
    TargetedSourceWallCapture,
    _required_tables,
    check,
    load_config,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    Proto11GR0EvolutionOperator,
    exact_spherical_minkowski_reference,
)


class FGCSRC2CaptureRunnerTests(unittest.TestCase):
    def test_canonical_config_reaches_only_the_capture_boundary(self) -> None:
        result = check()
        self.assertTrue(result["authorized_to_capture"])
        self.assertEqual(result["source_only_rejection_count"], 164)
        self.assertEqual(
            result["first_common_event_state_sha256"],
            "212724c1b7615d894d3e0d9ad2488dfe60c8d11957f0fa74c2c61ec324c25e7d",
        )
        self.assertEqual(
            result["terminal_accepted_state_sha256"],
            "414e3805cd2693d73770f6f7222ea27cdaf568299c7f913e6029457320792bab",
        )
        self.assertFalse(result["output_created_by_check"])
        self.assertFalse(result["PROTO12_frozen"])
        self.assertFalse(result["FGCQR_holdout_execution_authorized"])

    def test_threshold_scope_and_claim_mutations_fail_closed(self) -> None:
        config = load_config()
        attacked = deepcopy(config)
        attacked["arithmetic"]["raw_complete_residual_maximum"] = "2/1000000000000"
        with self.assertRaisesRegex(ValueError, "arithmetic contract"):
            _required_tables(attacked)
        attacked = deepcopy(config)
        attacked["scope"]["FGCQR_trajectory_read"] = True
        with self.assertRaisesRegex(ValueError, "scope"):
            _required_tables(attacked)
        attacked = deepcopy(config)
        attacked["claims"]["PROTO12_frozen"] = True
        with self.assertRaisesRegex(ValueError, "claim"):
            _required_tables(attacked)

    def test_target_capture_requires_bitwise_time_and_residual_identity(self) -> None:
        grid = UniformRadialGrid(0.0, 8.0, 65)
        state = exact_spherical_minkowski_reference(grid)
        operator = Proto11GR0EvolutionOperator(
            grid,
            spatial_order=4,
            ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
        )
        residual = float(operator(0.0, state).diagnostics["source_residual_infinity"])
        capture = TargetedSourceWallCapture(
            operator,
            target_time=0.0,
            target_residual=residual,
        )
        capture(0.0, state)
        capture(float(np.nextafter(0.0, 1.0)), state)
        self.assertEqual(capture.total_calls, 2)
        self.assertEqual(capture.target_matches, 1)
        self.assertIsNotNone(capture.captured_state)
        self.assertIsNotNone(capture.captured_rhs)


if __name__ == "__main__":
    unittest.main()
