"""Focused synthetic controls for the bounded TDG10-QA2 runner."""

from __future__ import annotations

import ast
import errno
from fractions import Fraction
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (
    ExactCompleteCClosedEvidence,
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime import (
    TDG10ExactCompleteCRuntimeClosed,
)

from scripts import run_fgc_tdg10_qa2 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg10_qa2_authority as authority  # noqa: E402


AUTHORITY_COMMIT = "a" * 40


def _fingerprint(retry: int) -> dict[str, object]:
    width = authority.width_spec(retry)
    return {
        "member_key": authority.MEMBER_KEY,
        "accepted_time_hex": authority.ACCEPTED_TIME_HEX,
        "state_sha256": authority.PHYSICAL_STATE_SHA256,
        "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": width["transaction_sha256"],
    }


def _member() -> SimpleNamespace:
    return SimpleNamespace(key=authority.MEMBER_KEY)


def _restored(retry: int) -> runner.ti2_runner.RestoredReplay:
    width = authority.width_spec(retry)
    replay = dict(authority.REPLAYS[0 if retry == 4 else 1])
    return runner.ti2_runner.RestoredReplay(
        replay,
        _member(),
        _fingerprint(retry),
        str(width["historical_journal_sha256"]),
    )


def _assessment(*, passing: bool = True) -> SimpleNamespace:
    channels = []
    failed: list[str] = []
    for name in authority.CHANNEL_ORDER:
        passed = passing or name != "u:R"
        if not passed:
            failed.append(name)
        channels.append(
            SimpleNamespace(
                channel=name,
                admission_passed=passed,
                classification=(
                    "resolved_order_pass" if passed else "resolved_order_failure"
                ),
                d12=SimpleNamespace(upper=Fraction(1, 8 if passed else 2)),
            )
        )
    return SimpleNamespace(
        channels=channels,
        complete_admission_passed=not failed,
        failed_channels=tuple(failed),
    )


def _receipt() -> authority.QA2Authority:
    return authority.QA2Authority(authority_commit=AUTHORITY_COMMIT)


def _sealed() -> tuple[int, str]:
    return (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )


def _closed(reason: str = "route_disagreement") -> ExactCompleteCClosedEvidence:
    return ExactCompleteCClosedEvidence(
        reason=reason,
        level="D01",
        row_count=1,
        maximum_candidates_D01=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        maximum_candidates_D12=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        detail="synthetic",
    )


class QA2RunnerTests(unittest.TestCase):
    def test_forbidden_production_apis_are_not_imported_or_called(self) -> None:
        source = (ROOT / authority.RUNNER_PATH).read_text()
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertNotIn("encode_member", imported)
        self.assertNotIn("encode_payload", imported)
        self.assertNotIn("commit_tdg6_gr0_compositor", imported)
        self.assertNotIn("require_tdg6_temporal_admission", imported)
        self.assertIn("_restore_replay", source)
        self.assertIn("_prepare_shadow", source)
        self.assertIn("EXACT_COMPLETE_C_RUNTIME_MODULE", source)
        self.assertIn("assess_tdg10_exact_complete_c", source)
        for forbidden in (
            "commit_tdg6_gr0_compositor(",
            "require_tdg6_temporal_admission(",
            "encode_member(",
            "encode_payload(",
            "serialize_diagnostic_endpoint(",
            "publish_checkpoint(",
            "acquire_writer(",
            "write_campaign(",
        ):
            self.assertNotIn(forbidden, source)
        self.assertEqual(
            runner.EXACT_COMPLETE_C_RUNTIME_MODULE,
            "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime",
        )
        self.assertEqual(
            runner.EXACT_COMPLETE_C_ASSESS_NAME,
            "assess_tdg10_exact_complete_c",
        )

    def test_restore_replay_refuses_generation_9_and_forbidden_successors(self) -> None:
        replay = dict(authority.REPLAYS[0])
        replay["predecessor_generation"] = 9
        with self.assertRaisesRegex(runner.QA2RunnerError, "must_not_restore_generation_9"):
            runner._restore_replay(ROOT, replay)
        replay = dict(authority.REPLAYS[0])
        replay["predecessor_generation"] = 11
        with self.assertRaisesRegex(
            runner.QA2RunnerError, "forbidden_predecessor_generation"
        ):
            runner._restore_replay(ROOT, replay)
        replay = dict(authority.REPLAYS[1])
        replay["predecessor_generation"] = 12
        with self.assertRaisesRegex(
            runner.QA2RunnerError, "forbidden_predecessor_generation"
        ):
            runner._restore_replay(ROOT, replay)

    def test_terminal_reduction_is_closed(self) -> None:
        def record(retry: int, width_class: str) -> dict[str, object]:
            return {"retry": retry, "width_class": width_class}

        self.assertEqual(
            runner.reduce_qa2_terminal(
                (
                    record(4, runner.WIDTH_PASS_CLASS),
                    record(5, runner.WIDTH_PASS_CLASS),
                )
            ),
            runner.BOTH_PASS_CLASS,
        )
        self.assertEqual(
            runner.reduce_qa2_terminal(
                (
                    record(4, runner.WIDTH_PASS_CLASS),
                    record(5, runner.WIDTH_NONPASS_CLASS),
                )
            ),
            runner.NONPASS_CLASS,
        )
        self.assertEqual(
            runner.reduce_qa2_terminal(
                (record(4, "inconclusive"), record(5, runner.WIDTH_PASS_CLASS))
            ),
            runner.INCONCLUSIVE_CLASS,
        )
        self.assertEqual(
            runner.reduce_qa2_terminal(
                (record(4, "premise_stop"), record(5, "inconclusive"))
            ),
            runner.PREMISE_STOP_CLASS,
        )
        with self.assertRaisesRegex(runner.QA2RunnerError, "unknown_width_class"):
            runner.reduce_qa2_terminal(
                (record(4, "unknown"), record(5, runner.WIDTH_PASS_CLASS))
            )

    def _run_mocked(
        self,
        root: Path,
        *,
        assess: object,
        prepare: object | None = None,
        assess_by_retry: dict[int, object] | None = None,
        prepare_by_retry: dict[int, object] | None = None,
    ) -> dict[str, object]:
        restored = {4: _restored(4), 5: _restored(5)}
        prepared = SimpleNamespace(fine=SimpleNamespace())

        def restore(_root: Path, replay: object) -> object:
            return restored[int(replay["retry"])]

        def prepare_one(item: object) -> object:
            retry = int(item.replay["retry"])
            if prepare_by_retry is not None and retry in prepare_by_retry:
                value = prepare_by_retry[retry]
                if isinstance(value, BaseException):
                    raise value
                return value
            if isinstance(prepare, BaseException):
                raise prepare
            return prepared if prepare is None else prepare

        def assess_one(_prepared: object) -> object:
            if assess_by_retry is not None:
                # assessments are requested in retry order 4 then 5
                retry = 4 + len(getattr(assess_one, "seen", ()))
                assess_one.seen = (*getattr(assess_one, "seen", ()), retry)
                value = assess_by_retry.get(retry, assess)
            else:
                value = assess
            if isinstance(value, BaseException):
                raise value
            return value

        def fingerprint(member: object, descriptor: object) -> object:
            for item in restored.values():
                if item.member is member:
                    return item.fingerprint
            return restored[4].fingerprint

        with (
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_snapshot_store", return_value=_sealed()),
            patch.object(runner, "_restore_replay", side_effect=restore),
            patch.object(runner, "_prepare_shadow", side_effect=prepare_one),
            patch.object(
                runner, "_assess_prepared_exact_complete_c", side_effect=assess_one
            ),
            patch.object(
                runner.ti2_runner,
                "_member_fingerprint",
                side_effect=fingerprint,
            ),
        ):
            return runner.run(root, authority_commit=AUTHORITY_COMMIT)

    def test_both_pass_creates_only_manifest_and_terminal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=_assessment(passing=True))
            namespace = root / authority.OUTPUT_NAMESPACE
            names = tuple(sorted(path.name for path in namespace.iterdir()))
            self.assertEqual(names, ("manifest.json", "terminal.json"))
            self.assertEqual(terminal["classification"], runner.BOTH_PASS_CLASS)
            self.assertTrue(terminal["width_robustness_passed"])
            self.assertFalse(terminal["diagnostic_fine_endpoint_serialized"])
            self.assertNotIn("diagnostic_fine_endpoint", terminal)
            self.assertFalse((namespace / "payloads").exists())
            self.assertFalse((namespace / "fine-endpoint").exists())
            self.assertFalse((namespace / "checkpoints").exists())
            self.assertFalse(terminal["state_advance_authorized"])
            self.assertFalse(terminal["production_method_earned"])
            self.assertTrue(terminal["store_unchanged"])
            self.assertEqual([item["retry"] for item in terminal["widths"]], [4, 5])
            observed = runner._inspect_output(
                root, authority_commit=AUTHORITY_COMMIT
            )
            self.assertEqual(observed["classification"], runner.BOTH_PASS_CLASS)
            self.assertTrue(observed["width_robustness_passed"])

    def test_one_width_nonpass_is_not_robustness(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(
                root,
                assess=_assessment(passing=True),
                assess_by_retry={
                    4: _assessment(passing=True),
                    5: _assessment(passing=False),
                },
            )
            namespace = root / authority.OUTPUT_NAMESPACE
            names = tuple(sorted(path.name for path in namespace.iterdir()))
            self.assertEqual(names, ("manifest.json", "terminal.json"))
            self.assertEqual(terminal["classification"], runner.NONPASS_CLASS)
            self.assertFalse(terminal["width_robustness_passed"])
            self.assertFalse((namespace / "payloads").exists())
            runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_nonpass_inconclusive_and_stop_create_no_payload(self) -> None:
        stop = runner.tdg6.TDG6RefinementPathStop(
            path="outer_0", cause=RuntimeError("synthetic premise")
        )
        disagreement = ExactCompleteCRouteDisagreement(_closed())
        cases = (
            ("nonpass", _assessment(passing=False), runner.NONPASS_CLASS),
            ("inconclusive", disagreement, runner.INCONCLUSIVE_CLASS),
            ("stop", "prepare-stop", runner.PREMISE_STOP_CLASS),
        )
        for label, assess, classification in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                if assess == "prepare-stop":
                    terminal = self._run_mocked(
                        root, assess=_assessment(passing=True), prepare=stop
                    )
                else:
                    terminal = self._run_mocked(root, assess=assess)
                namespace = root / authority.OUTPUT_NAMESPACE
                names = tuple(sorted(path.name for path in namespace.iterdir()))
                self.assertEqual(names, ("manifest.json", "terminal.json"))
                self.assertEqual(terminal["classification"], classification)
                self.assertFalse(terminal["width_robustness_passed"])
                self.assertFalse((namespace / "payloads").exists())
                self.assertNotIn("diagnostic_fine_endpoint", terminal)
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_member_mutation_is_rejected_on_every_shadow_exit(self) -> None:
        stop = runner.tdg6.TDG6RefinementPathStop(
            path="outer_0", cause=RuntimeError("synthetic premise")
        )
        disagreement = ExactCompleteCRouteDisagreement(_closed())
        cases = (
            ("completed", None, _assessment(passing=True)),
            ("premise", stop, _assessment(passing=True)),
            ("inconclusive", None, disagreement),
        )
        for label, prepare, assess in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                restored = {4: _restored(4), 5: _restored(5)}
                mutated = dict(_fingerprint(4))
                mutated["state_sha256"] = "f" * 64

                def fingerprint(member: object, descriptor: object) -> object:
                    del descriptor
                    return mutated if member is restored[4].member else restored[5].fingerprint

                with (
                    patch.object(runner, "_execution_authority", return_value=_receipt()),
                    patch.object(runner, "_snapshot_store", return_value=_sealed()),
                    patch.object(
                        runner,
                        "_restore_replay",
                        side_effect=lambda _root, replay: restored[int(replay["retry"])],
                    ),
                    patch.object(
                        runner,
                        "_prepare_shadow",
                        side_effect=prepare if prepare is not None else None,
                        return_value=SimpleNamespace(fine=SimpleNamespace()),
                    ),
                    patch.object(
                        runner,
                        "_assess_prepared_exact_complete_c",
                        side_effect=assess if isinstance(assess, BaseException) else None,
                        return_value=assess,
                    ),
                    patch.object(
                        runner.ti2_runner,
                        "_member_fingerprint",
                        side_effect=fingerprint,
                    ),
                    self.assertRaisesRegex(runner.QA2RunnerError, "member_mutated"),
                ):
                    runner.run(root, authority_commit=AUTHORITY_COMMIT)
                self.assertFalse((root / authority.OUTPUT_NAMESPACE).exists())

    def test_invalid_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(runner, "_execution_authority", return_value=_receipt()),
                patch.object(runner, "_snapshot_store", return_value=_sealed()),
                patch.object(
                    runner,
                    "_restore_replay",
                    side_effect=runner.QA2RunnerError("replay", "broken", "synthetic"),
                ),
                self.assertRaisesRegex(runner.QA2RunnerError, "broken"),
            ):
                runner.run(root, authority_commit=AUTHORITY_COMMIT)
            self.assertFalse((root / authority.OUTPUT_NAMESPACE).exists())
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            if parent.exists():
                self.assertFalse(
                    any(
                        item.name.startswith(authority.STAGING_PREFIX)
                        for item in parent.iterdir()
                    )
                )

    def test_output_publication_and_status_reject_symlink_or_partial_namespace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaises(OSError):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})
            self.assertEqual(tuple(Path(foreign).iterdir()), ())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            (target / "manifest.json").write_text("{}\n")
            with self.assertRaisesRegex(runner.QA2RunnerError, "partial_or_foreign"):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_publication_rejects_foreign_or_racing_namespace_without_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            marker = target / "foreign"
            marker.write_text("keep")
            with self.assertRaisesRegex(runner.QA2RunnerError, "namespace_exists"):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})
            self.assertEqual(marker.read_text(), "keep")

        class _RaceRename:
            argtypes = None
            restype = None

            def __call__(self, *_arguments: object) -> int:
                import ctypes

                ctypes.set_errno(errno.EEXIST)
                return -1

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_library = SimpleNamespace(renameatx_np=_RaceRename())
            with (
                patch.object(runner.ctypes, "CDLL", return_value=fake_library),
                self.assertRaisesRegex(runner.QA2RunnerError, "namespace_arrived"),
            ):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            self.assertFalse((root / authority.OUTPUT_NAMESPACE).exists())
            self.assertFalse(
                any(
                    item.name.startswith(authority.STAGING_PREFIX)
                    for item in parent.iterdir()
                )
            )

    def test_no_terminal_can_publish_after_store_drift(self) -> None:
        with (
            patch.object(runner, "_publish") as publish,
            self.assertRaisesRegex(runner.QA2RunnerError, "campaign_store_mutated"),
        ):
            runner._publish_terminal(
                ROOT,
                {"kind": "manifest"},
                {"classification": runner.NONPASS_CLASS},
                store_before=(115, "a" * 64),
                store_after=(116, "b" * 64),
            )
        publish.assert_not_called()

    def test_status_is_read_only_and_does_not_require_absence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._run_mocked(root, assess=_assessment(passing=False))
            with (
                patch.object(runner, "_execution_authority", return_value=_receipt()),
                patch.object(runner, "_snapshot_store", return_value=_sealed()),
                patch.object(
                    authority,
                    "require_output_absent",
                    side_effect=AssertionError("status required absence"),
                ) as absence,
            ):
                observed = runner.status(root, authority_commit=AUTHORITY_COMMIT)
            absence.assert_not_called()
            self.assertTrue(observed["status_read_only"])
            self.assertFalse(observed["safe_to_run"])
            self.assertFalse(observed["state_advance_authorized"])
            self.assertEqual(
                observed["output"]["classification"], runner.NONPASS_CLASS
            )

    def test_resource_exhaustion_is_inconclusive_not_invalid(self) -> None:
        exhausted = ExactCompleteCResourceExhausted(
            _closed("candidate_ceiling_exhausted")
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=exhausted)
        self.assertEqual(terminal["classification"], runner.INCONCLUSIVE_CLASS)
        self.assertEqual(
            terminal["typed_inconclusive"]["type"],
            "ExactCompleteCResourceExhausted",
        )
        self.assertFalse(terminal["width_robustness_passed"])

    def test_runtime_closed_wrapper_is_inconclusive(self) -> None:
        closed = TDG10ExactCompleteCRuntimeClosed(
            channel="u:R",
            cause=ExactCompleteCRouteDisagreement(_closed()),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=closed)
        self.assertEqual(terminal["classification"], runner.INCONCLUSIVE_CLASS)
        self.assertEqual(
            terminal["typed_inconclusive"]["type"],
            "TDG10ExactCompleteCRuntimeClosed",
        )
        self.assertFalse((root / authority.OUTPUT_NAMESPACE / "payloads").exists())

    def test_invalid_class_is_not_publishable(self) -> None:
        widths = [
            {
                "retry": 4,
                "width_class": "premise_stop",
                "predecessor_generation": 10,
                "attempted_width_hex": authority.WIDTHS[0]["attempted_width_hex"],
                "replay_receipt": authority.expected_raw_replay_receipt(4),
                "typed_stop": {"type": "TDG6RefinementPathStop", "detail": "x"},
            },
            {
                "retry": 5,
                "width_class": "premise_stop",
                "predecessor_generation": 11,
                "attempted_width_hex": authority.WIDTHS[1]["attempted_width_hex"],
                "replay_receipt": authority.expected_raw_replay_receipt(5),
                "typed_stop": {"type": "TDG6RefinementPathStop", "detail": "x"},
            },
        ]
        terminal = {
            **runner._terminal_base(
                AUTHORITY_COMMIT,
                runner.PREMISE_STOP_CLASS,
                _sealed(),
                _sealed(),
                width_robustness_passed=False,
            ),
            "typed_stop": {"type": "TDG6RefinementPathStop", "detail": "x", "retry": 4},
            "widths": widths,
        }
        terminal["classification"] = runner.INVALID_CLASS
        with self.assertRaisesRegex(
            runner.QA2RunnerError, "terminal_class_not_publishable"
        ):
            runner._validate_terminal(terminal, authority_commit=AUTHORITY_COMMIT)

    def test_width_robustness_requires_both_passes(self) -> None:
        with self.assertRaisesRegex(
            runner.QA2RunnerError, "robustness_without_both_pass"
        ):
            runner._terminal_base(
                AUTHORITY_COMMIT,
                runner.NONPASS_CLASS,
                _sealed(),
                _sealed(),
                width_robustness_passed=True,
            )


if __name__ == "__main__":
    unittest.main()
