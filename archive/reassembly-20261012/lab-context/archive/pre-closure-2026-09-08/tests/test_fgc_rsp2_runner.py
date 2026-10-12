from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts import run_fgc_rsp2_constraint_study as runner
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256


class RSP2RunnerTests(unittest.TestCase):
    def test_authorization_and_frozen_input_bind(self) -> None:
        plan, authorization = runner._validate_authorization(
            runner.DEFAULT_PLAN,
            runner.DEFAULT_AUTHORIZATION,
            # RSP2 is an immutable completed predecessor. Its prospective
            # absence was proved at authorization; successor tests must not
            # replay that temporal observation after the namespace exists.
            require_fresh_namespace=False,
        )
        member = runner._build_member(plan, authorization)
        self.assertEqual(member.key, "SSPRK3-16385")
        self.assertEqual(member.time, 0.0)
        self.assertEqual(member.source_retry_count, 0)
        self.assertFalse(authorization["gate_status"]["FGCQR_holdout_execution_authorized"])

    def test_runtime_and_invalid_stops_cannot_be_persistence(self) -> None:
        typed = runner._typed_runtime_stop(
            reason="constraint_loss",
            member="SSPRK3-16385",
            evidence={"gate": "constraint"},
            failures=("radial_momentum",),
        )
        self.assertEqual(typed["classification"], "stopped_before_RSP2_endpoint")
        self.assertTrue(typed["runtime_stop_is_not_target_persistence"])
        self.assertNotIn("RSP2_target_order_cleared", typed)

        invalid = runner._invalid_runtime_stop(
            RuntimeError("synthetic implementation failure"),
            member="SSPRK3-16385",
        )
        self.assertEqual(
            invalid["classification"],
            "invalid_implementation_or_nonconverged_run",
        )
        self.assertTrue(invalid["invalid_run_is_not_target_persistence"])

    def test_initial_checkpoint_round_trip_is_state_exact(self) -> None:
        plan, authorization = runner._validate_authorization(
            runner.DEFAULT_PLAN,
            runner.DEFAULT_AUTHORIZATION,
            require_fresh_namespace=False,
        )
        member = runner._build_member(plan, authorization)
        manifest = runner._manifest(
            runner.DEFAULT_PLAN,
            runner.DEFAULT_AUTHORIZATION,
            authorization,
        )
        expected_hash = array_content_sha256(
            member.state.u, member.state.p, member.state.q
        )
        with TemporaryDirectory() as directory:
            root = Path(directory)
            event_log = root / "events.jsonl"
            event_log.write_bytes(b'{"event_type":"test_checkpoint"}\n')
            checkpoint_path = root / "latest-checkpoint.npz"
            metadata = runner._checkpoint_metadata(
                manifest=manifest,
                boundary_index=0,
                event_log_path=event_log,
                member=member,
                terminal=False,
                terminal_result=None,
            )
            runner.inherited._write_checkpoint(
                checkpoint_path,
                metadata=metadata,
                members={member.key: member},
                event_log_path=event_log,
            )
            restored_metadata, restored, restored_log = runner._restore_checkpoint(
                checkpoint_path,
                plan=plan,
                authorization=authorization,
                manifest=manifest,
            )
        self.assertEqual(restored_metadata, metadata)
        self.assertEqual(restored_log, b'{"event_type":"test_checkpoint"}\n')
        self.assertEqual(
            array_content_sha256(restored.state.u, restored.state.p, restored.state.q),
            expected_hash,
        )
        self.assertEqual(restored.transaction.state, member.transaction.state)
        self.assertEqual(restored.transaction.causal_state, member.transaction.causal_state)


if __name__ == "__main__":
    unittest.main()
