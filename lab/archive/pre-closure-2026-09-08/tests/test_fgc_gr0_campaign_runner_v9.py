from __future__ import annotations

from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts.run_fgc_gr0_calibration_v9 import (  # noqa: E402
    DEFAULT_AUTHORIZATION,
    DEFAULT_PLAN,
    RUNNER_ID,
    _INHERITED_EVENT_ASSESSMENT,
    _INHERITED_RUNNER_ID,
    _INHERITED_VALIDATE_AUTHORIZATION,
    _bound_inherited_engine,
    _validate_authorization,
)


class FGCGR0CampaignRunnerV9Tests(unittest.TestCase):
    def test_closed_HLT7_cannot_launch(self) -> None:
        with self.assertRaisesRegex(ValueError, "absent, incomplete, or over-broad"):
            _validate_authorization(
                DEFAULT_PLAN,
                DEFAULT_AUTHORIZATION,
                require_fresh_namespace=True,
            )

    def test_process_local_bindings_are_exact_and_exception_safe(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "synthetic interruption"):
            with _bound_inherited_engine():
                self.assertEqual(proto8_runner.RUNNER_ID, RUNNER_ID)
                self.assertIsNot(
                    proto8_runner._validate_authorization,
                    _INHERITED_VALIDATE_AUTHORIZATION,
                )
                self.assertIsNot(
                    proto8_runner._event_assessment,
                    _INHERITED_EVENT_ASSESSMENT,
                )
                raise RuntimeError("synthetic interruption")
        self.assertEqual(proto8_runner.RUNNER_ID, _INHERITED_RUNNER_ID)
        self.assertIs(
            proto8_runner._validate_authorization,
            _INHERITED_VALIDATE_AUTHORIZATION,
        )
        self.assertIs(proto8_runner._event_assessment, _INHERITED_EVENT_ASSESSMENT)

    def test_CAL6_plan_carries_only_the_frozen_ladder_at_runtime(self) -> None:
        with DEFAULT_PLAN.open("rb") as handle:
            plan = tomllib.load(handle)
        self.assertEqual(plan["numerics"]["resolutions"], [2049, 4097, 8193])
        self.assertEqual(plan["scope"]["branch"], "GR-0")
        self.assertTrue(
            plan["scope"]["candidate_action_fields_or_health_stops_forbidden"]
        )
        self.assertFalse(plan["claims"]["FGCQR_holdout_execution_authorized"])


if __name__ == "__main__":
    unittest.main()
