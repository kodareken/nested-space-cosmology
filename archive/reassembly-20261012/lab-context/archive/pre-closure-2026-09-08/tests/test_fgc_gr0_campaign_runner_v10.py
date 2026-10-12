from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts.reproduce_fgc_hlt4_mon4 import _canonical  # noqa: E402
from scripts.reproduce_fgc_hlt8_mon8 import load_canonical_result  # noqa: E402
from scripts.run_fgc_gr0_calibration_v10 import (  # noqa: E402
    DEFAULT_AUTHORIZATION,
    DEFAULT_PLAN,
    RUNNER_ID,
    _INHERITED_RUNNER_ID,
    _INHERITED_VALIDATE_AUTHORIZATION,
    _bound_inherited_engine,
    _event_assessment,
    _validate_authorization,
)
from scripts import run_fgc_gr0_calibration_v9 as proto9_runner  # noqa: E402


class FGCGR0CampaignRunnerV10Tests(unittest.TestCase):
    def test_HLT8_authorizes_only_the_frozen_GR0_plan(self) -> None:
        plan, authorization = _validate_authorization(
            DEFAULT_PLAN,
            DEFAULT_AUTHORIZATION,
            require_fresh_namespace=False,
        )
        self.assertEqual(plan["artifact_id"], "FGC-1-CAL7-RUN1-PLAN")
        gates = authorization["gate_status"]
        self.assertTrue(gates["PROTO10_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])

    def test_promoted_or_closed_authorization_fails(self) -> None:
        attacked = load_canonical_result(DEFAULT_AUTHORIZATION)
        attacked["gate_status"][
            "PROTO10_fresh_GR0_dynamic_calibration_authorized"
        ] = False
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "closed.json"
            path.write_text(_canonical(attacked), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "absent, incomplete, or over-broad"):
                _validate_authorization(
                    DEFAULT_PLAN,
                    path,
                    require_fresh_namespace=False,
                )

    def test_process_local_bindings_are_exact_and_exception_safe(self) -> None:
        original_event = proto9_runner._INHERITED_EVENT_ASSESSMENT
        with self.assertRaisesRegex(RuntimeError, "synthetic interruption"):
            with _bound_inherited_engine():
                self.assertEqual(proto8_runner.RUNNER_ID, RUNNER_ID)
                self.assertIsNot(
                    proto8_runner._validate_authorization,
                    _INHERITED_VALIDATE_AUTHORIZATION,
                )
                self.assertIs(proto8_runner._event_assessment, _event_assessment)
                raise RuntimeError("synthetic interruption")
        self.assertEqual(proto8_runner.RUNNER_ID, _INHERITED_RUNNER_ID)
        self.assertIs(
            proto8_runner._validate_authorization,
            _INHERITED_VALIDATE_AUTHORIZATION,
        )
        self.assertIs(proto8_runner._event_assessment, original_event)

    def test_CAL7_plan_preserves_physics_and_freezes_guard_semantics(self) -> None:
        with DEFAULT_PLAN.open("rb") as handle:
            plan = tomllib.load(handle)
        self.assertEqual(plan["numerics"]["resolutions"], [2049, 4097, 8193])
        self.assertEqual(plan["scope"]["branch"], "GR-0")
        self.assertTrue(
            plan["scope"]["candidate_action_fields_or_health_stops_forbidden"]
        )
        self.assertTrue(
            plan["numerics"][
                "common_event_spatial_coarse_to_medium_tail_ratios_required_directly"
            ]
        )
        self.assertTrue(
            plan["numerics"][
                "common_event_spatial_finest_pair_direct_or_conditioning_guarded_saturation_required"
            ]
        )
        self.assertFalse(plan["claims"]["FGCQR_holdout_execution_authorized"])


if __name__ == "__main__":
    unittest.main()
