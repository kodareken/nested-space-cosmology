"""Focused orchestration tests for the one-event TDG8 successor runner."""
from __future__ import annotations

import ast
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import replace
from hashlib import sha256
import io
import inspect
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from scripts import run_fgc_tdg8_successor_event as runner


ROOT = Path(__file__).resolve().parents[1]
AUTHORITY = "2" * 40
PLAN_SHA256 = "3" * 64
CHECKPOINT_SHA256 = "4" * 64
CONFIG_RAW = b"frozen successor config\n"
RESULT_RAW = b"frozen successor result\n"


def _completed(returncode: int, stdout: bytes = b"") -> SimpleNamespace:
    return SimpleNamespace(returncode=returncode, stdout=stdout, stderr=b"")


def _plan(**changes):
    values = {
        "authorization_commit": AUTHORITY,
        "protocol_artifact_id": runner.EXPECTED_PROTOCOL,
        "campaign_id": runner.tdg8.CAMPAIGN_ID,
        "branch": "GR-0",
        "amplitude": "3",
        "common_event_index": 23,
        "start_time": SimpleNamespace(
            rational="23/16", binary64_hex=(23 / 16).hex()
        ),
        "target_time": SimpleNamespace(
            rational="3/2", binary64_hex=(3 / 2).hex()
        ),
        "members": tuple(
            SimpleNamespace(member_key=key) for key in runner.MEMBER_KEYS
        ),
        "sha256": PLAN_SHA256,
    }
    values.update(changes)
    return SimpleNamespace(**values)


def _source(**changes):
    values = {
        "branch": "GR-0",
        "amplitude": "3",
        "members": {key: object() for key in runner.MEMBER_KEYS},
    }
    values.update(changes)
    return SimpleNamespace(**values)


def _authority_receipt(**changes):
    inventory = tuple(
        runner.authority.ImplementationBinding(
            role=role,
            path=path,
            sha256=(f"{index + 1:x}" * 64)[:64],
        )
        for index, (role, path) in enumerate(
            runner.authority.REQUIRED_IMPLEMENTATION_INVENTORY
        )
    )
    values = {
        "authority_commit": AUTHORITY,
        "predecessor_repair_commit": runner.authority.PREDECESSOR_REPAIR_COMMIT,
        "config_sha256": sha256(CONFIG_RAW).hexdigest(),
        "result_sha256": sha256(RESULT_RAW).hexdigest(),
        "campaign_id": runner.tdg8.CAMPAIGN_ID,
        "target_protocol": runner.EXPECTED_PROTOCOL,
        "branch": "GR-0",
        "amplitude": "3",
        "event": 23,
        "start_rational": "23/16",
        "target_rational": "3/2",
        "event_count": 1,
        "source_campaign_path": runner.EXPECTED_SOURCE_RELATIVE.as_posix(),
        "destination_campaign_path": (
            runner.EXPECTED_DESTINATION_RELATIVE.as_posix()
        ),
        "source_terminal_read_only": True,
        "candidate_branches_forbidden": True,
        "immutable_bindings": runner.authority.IMMUTABLE_BINDINGS,
        "implementation_inventory": inventory,
        "environment": dict(runner.authority.PINNED_ENVIRONMENT),
    }
    values.update(changes)
    return runner.authority.TDG8SuccessorExecutionAuthority(**values)


def _authority_config(receipt) -> dict[str, object]:
    return {
        "implementation_inventory": [
            {"role": item.role, "path": item.path, "sha256": item.sha256}
            for item in receipt.implementation_inventory
        ]
    }


def _checkpoint(plan=None, *, event: int = 23, disposition: str = "nonterminal"):
    plan = _plan() if plan is None else plan
    target = (
        {"rational": "3/2", "binary64_hex": (3 / 2).hex()}
        if event == 23
        else {"rational": "25/16", "binary64_hex": (25 / 16).hex()}
    )
    return SimpleNamespace(
        authorization_commit=plan.authorization_commit,
        plan_sha256=plan.sha256,
        campaign_id=plan.campaign_id,
        protocol=plan.protocol_artifact_id,
        event=event,
        target=target,
        members={key: object() for key in runner.MEMBER_KEYS},
        generation=1,
        sha256=CHECKPOINT_SHA256,
        disposition=disposition,
    )


def _result(plan=None, *, event: int = 24, state: str = "event_one_complete"):
    plan = _plan() if plan is None else plan
    return {
        "authorization_commit": plan.authorization_commit,
        "plan_sha256": plan.sha256,
        "campaign_id": plan.campaign_id,
        "event": event,
        "target": {
            "rational": "25/16" if event == 24 else "3/2",
            "binary64_hex": (25 / 16).hex() if event == 24 else (3 / 2).hex(),
        },
        "state": state,
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


class TDG8SuccessorRunnerTests(unittest.TestCase):
    def test_src_package_precedes_untracked_repository_imports(self) -> None:
        self.assertLess(
            sys.path.index(str(ROOT / "src")),
            sys.path.index(str(ROOT)),
        )

    def test_git_authority_requires_exact_head_and_ignores_untracked_inventory(self) -> None:
        calls: list[tuple[str, ...]] = []

        def git(_root, *arguments):
            calls.append(arguments)
            if arguments[0] == "rev-parse":
                return _completed(0, f"{AUTHORITY}\n".encode("ascii"))
            return _completed(0)

        with patch.object(runner, "_git", side_effect=git):
            self.assertEqual(
                runner._require_git_authority(ROOT, AUTHORITY), AUTHORITY
            )
        self.assertEqual(calls[0], ("rev-parse", "--verify", "HEAD"))
        self.assertEqual(calls[1][0], "diff")
        self.assertNotIn("status", calls[1])
        self.assertNotIn("--others", calls[1])

    def test_git_queries_scrub_hostile_ambient_repository_state(self) -> None:
        completed = _completed(0, b"ok\n")
        with patch.dict(
            runner.os.environ,
            {
                "GIT_DIR": "/tmp/foreign.git",
                "GIT_WORK_TREE": "/tmp/foreign-tree",
                "GIT_CONFIG_COUNT": "1",
                "GIT_CONFIG_KEY_0": "core.fsmonitor",
                "GIT_CONFIG_VALUE_0": "hostile",
            },
            clear=False,
        ), patch.object(
            runner.subprocess, "run", return_value=completed
        ) as invoke:
            self.assertIs(runner._git(ROOT, "rev-parse", "HEAD"), completed)

        command = invoke.call_args.args[0]
        options = invoke.call_args.kwargs
        self.assertIn("--no-replace-objects", command)
        self.assertIn("--no-optional-locks", command)
        self.assertIs(options["stdin"], runner.subprocess.DEVNULL)
        self.assertEqual(options["env"]["GIT_OPTIONAL_LOCKS"], "0")
        for hostile in (
            "GIT_DIR",
            "GIT_WORK_TREE",
            "GIT_CONFIG_COUNT",
            "GIT_CONFIG_KEY_0",
            "GIT_CONFIG_VALUE_0",
        ):
            self.assertNotIn(hostile, options["env"])

    def test_wrong_commit_and_dirty_tracked_tree_fail_before_execution(self) -> None:
        with patch.object(
            runner,
            "_git",
            return_value=_completed(0, f"{'1' * 40}\n".encode("ascii")),
        ):
            with self.assertRaisesRegex(
                runner.TDG8SuccessorRunnerError, "differs from live HEAD"
            ):
                runner._require_git_authority(ROOT, AUTHORITY)

        with patch.object(
            runner,
            "_git",
            side_effect=(
                _completed(0, f"{AUTHORITY}\n".encode("ascii")),
                _completed(1),
            ),
        ):
            with self.assertRaisesRegex(
                runner.TDG8SuccessorRunnerError, "tracked tree differs"
            ):
                runner._require_git_authority(ROOT, AUTHORITY)

    def test_fresh_destination_authenticates_installs_and_initializes_once(self) -> None:
        plan = _plan()
        source = _source()
        checkpoint = _checkpoint(plan)
        store = object()
        with patch.object(runner, "_execution_authority", return_value=_authority_receipt()), \
             patch.object(runner, "_require_git_authority", return_value=AUTHORITY), \
             patch.object(runner, "_destination_state", return_value="absent"), \
             patch.object(runner.tdg8, "construct_successor_plan", return_value=plan) as construct, \
             patch.object(runner.tdg8, "authenticate_and_reconstruct", return_value=source) as authenticate, \
             patch.object(runner.tdg8, "install_generation_zero_prefix") as install, \
             patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner.event1, "_recoverable_bootstrap_prefix", return_value=True), \
             patch.object(runner.runtime, "initialize_or_recover_generation_one", return_value=checkpoint) as initialize, \
             patch.object(runner.runtime, "recover_persisted_generation_one") as persisted, \
             patch.object(runner, "build_static_gr0_shells") as static_shells, \
             patch.object(runner.event1, "_advance_authenticated_event", return_value=_result(plan)) as advance:
            result = runner.run_successor_event(ROOT, authorization_commit=AUTHORITY)

        authenticate.assert_called_once_with(ROOT.resolve())
        construct.assert_called_once_with(ROOT.resolve(), AUTHORITY)
        install.assert_called_once_with(
            ROOT.resolve(),
            ROOT.resolve() / runner.EXPECTED_DESTINATION_RELATIVE,
            source,
            plan.sha256,
        )
        initialize.assert_called_once_with(plan, store, source)
        persisted.assert_not_called()
        static_shells.assert_not_called()
        advance.assert_called_once_with(
            plan, store, checkpoint, source.members, attempt_runner=None
        )
        self.assertEqual(result["event"], 24)
        self.assertEqual(result["event_limit"], 1)
        self.assertFalse(result["candidate_branch_opened"])

    def test_generation_one_resume_is_persisted_only(self) -> None:
        plan = _plan()
        checkpoint = _checkpoint(plan)
        store = object()
        templates = {key: object() for key in runner.MEMBER_KEYS}
        with patch.object(runner, "_execution_authority", return_value=_authority_receipt()), \
             patch.object(runner, "_require_git_authority", return_value=AUTHORITY), \
             patch.object(runner, "_destination_state", return_value="directory"), \
             patch.object(runner.tdg8, "construct_successor_plan", return_value=plan), \
             patch.object(runner.tdg8, "authenticate_and_reconstruct") as authenticate, \
             patch.object(runner.tdg8, "install_generation_zero_prefix") as install, \
             patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner.event1, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner.runtime, "initialize_or_recover_generation_one") as initialize, \
             patch.object(runner.runtime, "recover_persisted_generation_one", return_value=checkpoint) as persisted, \
             patch.object(runner, "build_static_gr0_shells", return_value=templates) as static_shells, \
             patch.object(runner.event1, "_advance_authenticated_event", return_value=_result(plan)) as advance:
            result = runner.run_successor_event(ROOT, authorization_commit=AUTHORITY)

        authenticate.assert_not_called()
        install.assert_not_called()
        initialize.assert_not_called()
        persisted.assert_called_once_with(plan, store)
        static_shells.assert_called_once_with(ROOT.resolve())
        advance.assert_called_once_with(
            plan, store, checkpoint, templates, attempt_runner=None
        )
        self.assertEqual(result["state"], "event_one_complete")

    def test_mutated_typed_authority_cannot_reach_bootstrap_or_evolution(self) -> None:
        exact = _authority_receipt()
        changed_inventory = (
            replace(exact.implementation_inventory[0], sha256="f" * 64),
            *exact.implementation_inventory[1:],
        )
        mutations = (
            replace(exact, campaign_id="foreign-campaign"),
            replace(exact, source_campaign_path="runs/foreign/source"),
            replace(exact, destination_campaign_path="runs/foreign/destination"),
            replace(exact, candidate_branches_forbidden=False),
            replace(exact, implementation_inventory=changed_inventory),
        )
        for mutated in mutations:
            with self.subTest(mutated=mutated), patch.object(
                runner.authority,
                "_read_repository_leaf",
                side_effect=(CONFIG_RAW, RESULT_RAW),
            ), patch.object(
                runner.authority,
                "authorize_execution",
                return_value=mutated,
            ), patch.object(
                runner.authority,
                "_parse_config",
                return_value=_authority_config(exact),
            ), patch.object(
                runner, "_destination_state"
            ) as destination_state, patch.object(
                runner.tdg8, "construct_successor_plan"
            ) as construct, patch.object(
                runner.tdg8, "authenticate_and_reconstruct"
            ) as authenticate, patch.object(
                runner.tdg8, "install_generation_zero_prefix"
            ) as install, patch.object(
                runner.event1, "_advance_authenticated_event"
            ) as advance:
                with self.assertRaisesRegex(
                    runner.TDG8SuccessorRunnerError,
                    "execution authority scope differs",
                ):
                    runner.run_successor_event(
                        ROOT,
                        authorization_commit=AUTHORITY,
                    )
            destination_state.assert_not_called()
            construct.assert_not_called()
            authenticate.assert_not_called()
            install.assert_not_called()
            advance.assert_not_called()

    def test_failed_authority_validation_cannot_inspect_or_create_destination(self) -> None:
        with patch.object(
            runner.authority,
            "_read_repository_leaf",
            side_effect=(CONFIG_RAW, RESULT_RAW),
        ), patch.object(
            runner.authority,
            "authorize_execution",
            side_effect=runner.authority.TDG8SuccessorAuthorityError(
                "injected invalid authority"
            ),
        ), patch.object(
            runner, "_destination_state"
        ) as destination_state, patch.object(
            runner.tdg8, "construct_successor_plan"
        ) as construct, patch.object(
            runner.tdg8, "install_generation_zero_prefix"
        ) as install, patch.object(
            runner.event1, "_advance_authenticated_event"
        ) as advance:
            with self.assertRaisesRegex(
                runner.TDG8SuccessorRunnerError,
                "execution authority failed",
            ):
                runner.run_successor_event(ROOT, authorization_commit=AUTHORITY)
        destination_state.assert_not_called()
        construct.assert_not_called()
        install.assert_not_called()
        advance.assert_not_called()

    def test_plan_branch_campaign_and_destination_are_immutable(self) -> None:
        for changes in (
            {"branch": "candidate"},
            {"amplitude": "5/2"},
            {"campaign_id": "foreign"},
            {"common_event_index": 24},
        ):
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(
                    runner.TDG8SuccessorRunnerError, "GR-0 event scope"
                ):
                    runner._require_plan_scope(_plan(**changes), AUTHORITY)

        with patch.object(
            runner.tdg8, "DESTINATION_RELATIVE", Path("runs/foreign/calibration")
        ):
            with self.assertRaisesRegex(
                runner.TDG8SuccessorRunnerError, "destination authority differs"
            ):
                runner._fixed_destination(ROOT.resolve())

        tree = ast.parse(inspect.getsource(runner))
        imported = {
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertFalse(any("sgbl" in name.lower() for name in imported))
        self.assertFalse(any("fgcqr" in name.lower() for name in imported))
        main_source = inspect.getsource(runner.main)
        for forbidden_option in ("--root", "--destination", "--store-root"):
            self.assertNotIn(forbidden_option, main_source)

    def test_cli_requires_exactly_one_run_or_status_mode(self) -> None:
        with patch.object(
            sys,
            "argv",
            [
                "run_fgc_tdg8_successor_event.py",
                "--authority-commit",
                AUTHORITY,
                "--run",
                "--status",
            ],
        ):
            with redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    runner.main()
        self.assertEqual(caught.exception.code, 2)

    def test_result_cannot_advance_beyond_one_event_or_hide_a_nonterminal(self) -> None:
        plan = _plan()
        wrong_terminal_target = _result(
            plan, event=23, state="scientific_terminal"
        )
        wrong_terminal_target["target"] = {
            "rational": "25/16",
            "binary64_hex": (25 / 16).hex(),
        }
        for result in (
            _result(plan, event=25, state="event_one_complete"),
            _result(plan, event=23, state="nonterminal"),
            wrong_terminal_target,
            {**_result(plan), "target": "malformed"},
        ):
            with self.subTest(event=result["event"], state=result["state"]):
                with self.assertRaises(runner.TDG8SuccessorRunnerError):
                    runner._require_terminal_result(result, plan)

    def test_cli_serializes_a_reused_event_runner_failure_fail_closed(self) -> None:
        with patch.object(
            sys,
            "argv",
            [
                "run_fgc_tdg8_successor_event.py",
                "--authority-commit",
                AUTHORITY,
                "--run",
            ],
        ), patch.object(
            runner,
            "run_successor_event",
            side_effect=runner.event1.Proto19Event1RunnerError("injected failure"),
        ):
            output = io.StringIO()
            with redirect_stdout(output):
                code = runner.main()
        self.assertEqual(code, 2)
        self.assertIn("invalid_implementation_or_nonconverged_run", output.getvalue())
        self.assertIn('"candidate_branch_opened":false', output.getvalue())

    def test_status_absence_is_read_only_and_not_restartable(self) -> None:
        plan = _plan()
        with patch.object(runner, "_execution_authority", return_value=_authority_receipt()), \
             patch.object(runner, "_require_git_authority", return_value=AUTHORITY), \
             patch.object(runner, "_destination_state", return_value="absent"), \
             patch.object(runner.tdg8, "construct_successor_plan", return_value=plan), \
             patch.object(runner, "HLT16CampaignStore") as store, \
             patch.object(runner.tdg8, "authenticate_and_reconstruct") as authenticate, \
             patch.object(runner.tdg8, "install_generation_zero_prefix") as install:
            result = runner.collect_status(ROOT, authorization_commit=AUTHORITY)
        store.assert_not_called()
        authenticate.assert_not_called()
        install.assert_not_called()
        self.assertEqual(result["state"], "destination_absent")
        self.assertTrue(result["safe_to_start"])
        self.assertFalse(result["safe_to_restart"])

    def test_status_reports_authenticated_checkpoint_without_recovery_mutation(self) -> None:
        plan = _plan()
        checkpoint = _checkpoint(plan)
        recovery = SimpleNamespace(
            authorization_commit=AUTHORITY,
            plan_sha256=PLAN_SHA256,
            checkpoint_generation=1,
            checkpoint_sha256=CHECKPOINT_SHA256,
            member_times={},
            member_modes={},
            active_target=checkpoint.target,
            last_complete_journal_sequence=0,
            journal_tip_sha256="5" * 64,
            state="clean_checkpoint",
            safe_to_restart=True,
            terminal=False,
            active_write=False,
            writer_state=None,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=checkpoint)
            ),
            inspect_recovery=Mock(return_value=recovery),
        )
        with patch.object(runner, "_execution_authority", return_value=_authority_receipt()), \
             patch.object(runner, "_require_git_authority", return_value=AUTHORITY), \
             patch.object(runner, "_destination_state", return_value="directory"), \
             patch.object(runner.tdg8, "construct_successor_plan", return_value=plan), \
             patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner.event1, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner.runtime, "recover_persisted_generation_one") as recover, \
             patch.object(runner.tdg8, "authenticate_and_reconstruct") as authenticate, \
             patch.object(runner.tdg8, "install_generation_zero_prefix") as install:
            result = runner.collect_status(ROOT, authorization_commit=AUTHORITY)
        recover.assert_not_called()
        authenticate.assert_not_called()
        install.assert_not_called()
        store.authenticated_snapshot.assert_called_once_with()
        store.inspect_recovery.assert_called_once_with()
        self.assertEqual(result["state"], "clean_checkpoint")
        self.assertTrue(result["safe_to_restart"])


if __name__ == "__main__":
    unittest.main()
