from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts.reproduce_fgc_hlt4_mon4 import _canonical  # noqa: E402
from scripts.reproduce_fgc_hlt9_mon9 import load_canonical_result  # noqa: E402
from scripts.run_fgc_gr0_calibration_v11 import (  # noqa: E402
    DEFAULT_AUTHORIZATION,
    DEFAULT_PLAN,
    RUNNER_ID,
    _INHERITED_BUILD_MEMBERS,
    _INHERITED_PROTO6_BUILD_MEMBERS,
    _INHERITED_RUNNER_ID,
    _INHERITED_VALIDATE_AUTHORIZATION,
    _bound_inherited_engine,
    _build_members,
    _event_assessment,
    _validate_authorization,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    Proto11GR0EvolutionOperator,
)


class FGCGR0CampaignRunnerV11Tests(unittest.TestCase):
    def test_HLT9_authorizes_only_the_frozen_GR0_plan(self) -> None:
        plan, authorization = _validate_authorization(
            DEFAULT_PLAN,
            DEFAULT_AUTHORIZATION,
            require_fresh_namespace=False,
        )
        self.assertEqual(plan["artifact_id"], "FGC-1-CAL8-RUN1-PLAN")
        gates = authorization["gate_status"]
        self.assertTrue(gates["PROTO11_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["PROTO11_resolved_holdout_manifest_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])

    def test_promoted_or_closed_authorization_fails(self) -> None:
        attacked = load_canonical_result(DEFAULT_AUTHORIZATION)
        attacked["gate_status"][
            "PROTO11_fresh_GR0_dynamic_calibration_authorized"
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
        with self.assertRaisesRegex(RuntimeError, "synthetic interruption"):
            with _bound_inherited_engine():
                self.assertEqual(proto8_runner.RUNNER_ID, RUNNER_ID)
                self.assertIsNot(
                    proto8_runner._validate_authorization,
                    _INHERITED_VALIDATE_AUTHORIZATION,
                )
                self.assertIs(proto8_runner._build_members, _build_members)
                self.assertIs(proto8_runner._event_assessment, _event_assessment)
                self.assertIs(proto6_runner._build_members, _build_members)
                raise RuntimeError("synthetic interruption")
        self.assertEqual(proto8_runner.RUNNER_ID, _INHERITED_RUNNER_ID)
        self.assertIs(
            proto8_runner._validate_authorization,
            _INHERITED_VALIDATE_AUTHORIZATION,
        )
        self.assertIs(proto8_runner._build_members, _INHERITED_BUILD_MEMBERS)
        self.assertIs(
            proto6_runner._build_members, _INHERITED_PROTO6_BUILD_MEMBERS
        )

    def test_runtime_members_use_frozen_PROTO11_states_and_operator(self) -> None:
        plan, authorization = _validate_authorization(
            DEFAULT_PLAN,
            DEFAULT_AUTHORIZATION,
            require_fresh_namespace=False,
        )
        members = _build_members(plan, authorization, "5/2")
        self.assertEqual(len(members), 6)
        frozen = {
            (item["method"], item["point_count"]): item
            for item in authorization["artifact_payload"]["frozen_run_inputs"]
            if item["amplitude"] == "5/2"
        }
        for member in members.values():
            self.assertIsInstance(member.operator, Proto11GR0EvolutionOperator)
            self.assertEqual(
                member.input_hash,
                frozen[(member.method_label, member.point_count)][
                    "expanded_run_config_sha256"
                ],
            )
            self.assertFalse(
                member.operator(0.0, member.state).diagnostics[
                    "PROTO11_interior_q_reprojected"
                ]
            )


if __name__ == "__main__":
    unittest.main()
