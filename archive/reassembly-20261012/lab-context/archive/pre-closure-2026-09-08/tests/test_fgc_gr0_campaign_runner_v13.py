from __future__ import annotations

from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration_v13 as runner  # noqa: E402
from scripts.reproduce_fgc_hlt11_mon11 import (  # noqa: E402
    DEFAULT_OUTPUT,
    load_canonical_result,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)


class FGCGR0CampaignRunnerV13Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authorization = load_canonical_result(DEFAULT_OUTPUT)

    def test_complete_restart_restores_exactly(self) -> None:
        members = runner._restore_frozen_members(self.authorization)
        expected = {
            item["key"]: item
            for item in __import__("json").loads(
                (REPOSITORY / "results/fgc-1-pro13-frz1.json").read_text(
                    encoding="utf-8"
                )
            )["artifact_payload"]["restart_members"]
        }
        self.assertEqual(set(members), set(expected))
        for key, member in members.items():
            self.assertEqual(member.time, 23.0 / 16.0)
            self.assertEqual(len(member.tracers.event_fields), 24)
            self.assertEqual(
                array_content_sha256(member.state.u, member.state.p, member.state.q),
                expected[key]["state_sha256"],
            )

    def test_restart_event_uses_method_owned_ladders(self) -> None:
        members = runner._restore_frozen_members(self.authorization)
        plan = runner.hlt11.validate_run_plan(runner.DEFAULT_PLAN)
        event = runner._event_assessment(plan, members, 23.0 / 16.0)
        self.assertEqual(
            list(event["primary_common_event"].point_counts), [2049, 4097, 8193]
        )
        self.assertEqual(
            list(event["comparator_common_event"].point_counts),
            [4097, 8193, 16385],
        )
        self.assertFalse(event["temporal_spectral_admission"]["available"])
        self.assertTrue(event["method_owned_ladders_enforced"])

    def test_manifest_records_the_frozen_numerical_runtime(self) -> None:
        plan = runner.hlt11.validate_run_plan(runner.DEFAULT_PLAN)
        manifest = runner._manifest(
            runner.DEFAULT_PLAN, runner.DEFAULT_AUTHORIZATION, self.authorization
        )
        self.assertEqual(
            manifest["numerical_runtime_contract"], plan["numerical_runtime"]
        )
        self.assertEqual(manifest["runtime_environment"]["python_version"], "3.14.3")
        self.assertEqual(manifest["runtime_environment"]["numpy_version"], "2.5.1")
        self.assertEqual(manifest["runtime_environment"]["blas"]["name"], "accelerate")


if __name__ == "__main__":
    unittest.main()
