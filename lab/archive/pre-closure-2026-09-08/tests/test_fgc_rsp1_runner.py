from __future__ import annotations

from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_rsp1_resolution_study import (  # noqa: E402
    DEFAULT_AUTHORIZATION,
    DEFAULT_PLAN,
    _build_members,
    _validate_authorization,
)
from recursive_horizons.fgc.evolution.rsp1_resolution_runtime import (  # noqa: E402
    RSP1GR0EvolutionOperator,
    RSP1_POINT_COUNTS,
)


class FGCRSP1RunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.plan, cls.authorization = _validate_authorization(
            DEFAULT_PLAN,
            DEFAULT_AUTHORIZATION,
            require_fresh_namespace=False,
        )

    def test_authorization_is_narrow_and_pretrajectory(self) -> None:
        gates = self.authorization["gate_status"]
        self.assertTrue(gates["RSP1_runtime_implemented"])
        self.assertTrue(gates["RSP1_execution_authorized"])
        self.assertFalse(gates["amplitude_three_spectral_veto_cleared"])
        self.assertFalse(gates["classical_spherical_diagnostic_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])

    def test_runtime_builds_exactly_six_SRC3_members(self) -> None:
        members = _build_members(self.plan, self.authorization)
        self.assertEqual(len(members), 6)
        self.assertEqual(
            {(item.method_label, item.point_count) for item in members.values()},
            {
                (method, point_count)
                for method in ("RK4", "SSPRK3")
                for point_count in RSP1_POINT_COUNTS
            },
        )
        for member in members.values():
            self.assertIsInstance(member.operator, RSP1GR0EvolutionOperator)
            self.assertEqual(member.amplitude, "3")


if __name__ == "__main__":
    unittest.main()
