"""Focused synthetic controls for the bounded TDG9 AC1 runner."""

from __future__ import annotations

import ast
from copy import deepcopy
import errno
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    classify_tdg6_channel,
)


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ac1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ac1_authority as authority  # noqa: E402


AUTHORITY_COMMIT = "a" * 40


def _exact(value: float) -> dict[str, str]:
    return runner._exact_binary64(value, label="synthetic")


def _channel(
    name: str,
    *,
    outer_lo: float,
    outer_hi: float,
    fine_lo: float,
    fine_hi: float,
) -> dict[str, object]:
    decision = classify_tdg6_channel(
        CertifiedMagnitudeInterval(
            Fraction.from_float(outer_lo), Fraction.from_float(outer_hi)
        ),
        CertifiedMagnitudeInterval(
            Fraction.from_float(fine_lo), Fraction.from_float(fine_hi)
        ),
    )
    return {
        "channel": name,
        "classification": decision.classification,
        "admission_passed": decision.admission_passed,
        "temporal_retry_permitted": decision.temporal_retry_permitted,
        "order_threshold_resolved": decision.order_threshold_resolved,
        "order_threshold_passed": decision.order_threshold_passed,
        "outer": {"lower": _exact(outer_lo), "upper": _exact(outer_hi)},
        "finest": {"lower": _exact(fine_lo), "upper": _exact(fine_hi)},
        "finest_pair_debit": _exact(fine_hi),
    }


def _admission(channels: list[dict[str, object]]) -> dict[str, object]:
    failed = [
        item["channel"] for item in channels if item["admission_passed"] is not True
    ]
    serialized = {
        "method": "second_order_diagonal_norm_SBP_plus_SSPRK3",
        "owned_row_count": authority.OWNED_ROW_COUNT,
        "channel_order": list(authority.CHANNEL_ORDER),
        "channels": channels,
        "complete_admission_passed": not failed,
        "failed_channels": failed,
        "finest_pair_debit_vector": [item["finest_pair_debit"] for item in channels],
        "admission_is_all_of": True,
    }
    runner.validate_all_of_identity(serialized)
    return serialized


def _passing_channels() -> list[dict[str, object]]:
    return [
        _channel(name, outer_lo=16.0, outer_hi=16.0, fine_lo=1.0, fine_hi=1.0)
        for name in TDG6_COMPLETE_STATE_CHANNELS
    ]


def _premise_terminal(commit: str = AUTHORITY_COMMIT) -> dict[str, object]:
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    return {
        **runner._terminal_base(
            commit, "shadow_proposal_premise_stop", sealed, sealed
        ),
        "typed_stop": {
            "type": "TDG6RefinementPathStop",
            "detail": "synthetic bounded premise stop",
        },
        "replay_receipt": runner.expected_replay_receipt(),
    }


def _completed_terminal(*, passing: bool = True) -> dict[str, object]:
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    channels = _passing_channels()
    if not passing:
        channels[1] = _channel(
            "u:v",
            outer_lo=1e-12,
            outer_hi=1e-12,
            fine_lo=9e-13,
            fine_hi=9e-13,
        )
    admission = _admission(channels)
    classification = runner.reduce_retry3_terminal(
        complete_admission_passed=bool(admission["complete_admission_passed"]),
        failed_channels=tuple(admission["failed_channels"]),
    )
    return {
        **runner._terminal_base(AUTHORITY_COMMIT, classification, sealed, sealed),
        "replay_receipt": runner.expected_replay_receipt(),
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "continuous_admission": admission,
    }


def _publish_premise(root: Path, *, terminal: dict[str, object] | None = None) -> None:
    runner._publish(
        root,
        runner._manifest(AUTHORITY_COMMIT),
        _premise_terminal() if terminal is None else terminal,
    )


class AC1RunnerTests(unittest.TestCase):
    def test_runner_contains_no_admission_commit_or_localization_path(self) -> None:
        source = (ROOT / authority.RUNNER_PATH).read_text()
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("candidate" in name.lower() for name in imported))
        self.assertNotIn("run_fgc_tdg9_loc1", imported)
        self.assertNotIn("run_fgc_tdg9_loc2", imported)
        for forbidden in (
            "require_tdg6_temporal_admission(",
            "commit_tdg6",
            "accept_step(",
            "acquire_writer",
            "publish_checkpoint",
            "localize_absolute_maximum",
            "FAILED_OCCURRENCES",
        ):
            self.assertNotIn(forbidden, source)
        self.assertIn("require_output_absent", source)
        self.assertIn("_restore_replay", source)
        self.assertIn("_prepare_shadow", source)

    def test_all_of_identity_holds_for_eighteen_passing_channels(self) -> None:
        serialized = _admission(_passing_channels())
        self.assertTrue(serialized["complete_admission_passed"])
        self.assertEqual(serialized["failed_channels"], [])
        self.assertEqual(
            runner.reduce_retry3_terminal(
                complete_admission_passed=True, failed_channels=()
            ),
            "completed_retry3_all_18_channels_pass",
        )

    def test_one_failed_channel_vetoes_complete_admission(self) -> None:
        channels = _passing_channels()
        channels[1] = _channel(
            "u:v",
            outer_lo=1e-12,
            outer_hi=1e-12,
            fine_lo=9e-13,
            fine_hi=9e-13,
        )
        serialized = _admission(channels)
        self.assertFalse(serialized["complete_admission_passed"])
        self.assertEqual(serialized["failed_channels"], ["u:v"])
        self.assertEqual(
            runner.reduce_retry3_terminal(
                complete_admission_passed=False, failed_channels=("u:v",)
            ),
            "completed_retry3_one_or_more_channels_fail",
        )

    def test_order_inconclusive_is_a_channel_failure_not_a_terminal(self) -> None:
        channels = _passing_channels()
        channels[2] = _channel(
            "u:lambda",
            outer_lo=0.0,
            outer_hi=1.0,
            fine_lo=0.25,
            fine_hi=0.5,
        )
        serialized = _admission(channels)
        self.assertEqual(channels[2]["classification"], "order_inconclusive")
        self.assertFalse(channels[2]["admission_passed"])
        self.assertEqual(serialized["failed_channels"], ["u:lambda"])
        self.assertNotEqual(
            runner.reduce_retry3_terminal(
                complete_admission_passed=False, failed_channels=("u:lambda",)
            ),
            "order_inconclusive",
        )

    def test_reordered_or_disagreed_flags_fail_closed(self) -> None:
        serialized = _admission(_passing_channels())
        reordered = deepcopy(serialized)
        reordered["channels"] = list(reversed(reordered["channels"]))
        with self.assertRaisesRegex(runner.AC1RunnerError, "channel_order"):
            runner.validate_all_of_identity(reordered)
        disagreed = deepcopy(serialized)
        disagreed["complete_admission_passed"] = False
        with self.assertRaisesRegex(runner.AC1RunnerError, "complete_pass_identity"):
            runner.validate_all_of_identity(disagreed)
        forged_fail = deepcopy(serialized)
        forged_fail["failed_channels"] = ["u:R"]
        with self.assertRaisesRegex(runner.AC1RunnerError, "failed_channels_identity"):
            runner.validate_all_of_identity(forged_fail)

    def test_reclassification_and_debit_identity_are_exact(self) -> None:
        serialized = _admission(_passing_channels())
        mutated = deepcopy(serialized)
        mutated["channels"][0]["classification"] = "resolved_order_failure"
        mutated["channels"][0]["admission_passed"] = False
        mutated["failed_channels"] = ["u:alpha"]
        mutated["complete_admission_passed"] = False
        with self.assertRaisesRegex(runner.AC1RunnerError, "reclassification"):
            runner.validate_all_of_identity(mutated)
        debit = deepcopy(serialized)
        debit["channels"][0]["finest_pair_debit"] = _exact(2.0)
        debit["finest_pair_debit_vector"][0] = _exact(2.0)
        with self.assertRaisesRegex(runner.AC1RunnerError, "reclassification"):
            runner.validate_all_of_identity(debit)

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
            with self.assertRaisesRegex(runner.AC1RunnerError, "partial_or_foreign"):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_status_authenticates_manifest_terminal_and_all_nonclaims(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _publish_premise(root)
            observed = runner._inspect_output(
                root, authority_commit=AUTHORITY_COMMIT
            )
            self.assertEqual(observed["state"], "terminal")
            self.assertEqual(
                observed["classification"], "shadow_proposal_premise_stop"
            )

        mutations: dict[str, tuple[dict[str, object], str]] = {}
        terminal = _premise_terminal()
        for name, value in (
            ("schema", "wrong-schema"),
            ("runner_id", "wrong-runner"),
            ("authority_commit", "b" * 40),
            ("production_SSPRK3_comparator", True),
            ("production_method_earned", True),
            ("store_unchanged", False),
        ):
            changed = deepcopy(terminal)
            changed[name] = value
            mutations[f"wrong-{name}"] = (changed, "terminal_identity")
        invalid = deepcopy(terminal)
        invalid["classification"] = "invalid_provenance_or_implementation"
        mutations["invalid-class"] = (invalid, "terminal_class_not_publishable")

        for label, (changed, error) in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _publish_premise(root, terminal=changed)
                with self.assertRaisesRegex(runner.AC1RunnerError, error):
                    runner._inspect_output(
                        root, authority_commit=AUTHORITY_COMMIT
                    )

    def test_status_authenticates_published_terminal_without_requiring_absence(self) -> None:
        sealed = (
            authority.SEALED_STORE_LEAF_COUNT,
            authority.SEALED_STORE_SNAPSHOT_SHA256,
        )
        receipt = authority.AC1Authority(authority_commit=AUTHORITY_COMMIT)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _publish_premise(root)
            with (
                patch.object(runner, "_execution_authority", return_value=receipt),
                patch.object(runner, "_snapshot_store", return_value=sealed),
            ):
                observed = runner.status(root, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(observed["output"]["state"], "terminal")
        self.assertEqual(
            observed["output"]["classification"], "shadow_proposal_premise_stop"
        )
        self.assertFalse(observed["safe_to_run"])
        self.assertFalse(observed["state_advance_authorized"])

    def test_replay_receipt_is_exact_on_premise_and_completed_terminals(self) -> None:
        runner._validate_terminal(
            _premise_terminal(), authority_commit=AUTHORITY_COMMIT
        )
        runner._validate_terminal(
            _completed_terminal(passing=True), authority_commit=AUTHORITY_COMMIT
        )
        runner._validate_terminal(
            _completed_terminal(passing=False), authority_commit=AUTHORITY_COMMIT
        )
        builders = (
            _premise_terminal,
            lambda: _completed_terminal(passing=True),
            lambda: _completed_terminal(passing=False),
        )
        for builder in builders:
            with self.subTest(builder=builder):
                empty = deepcopy(builder())
                empty["replay_receipt"] = {}
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "replay_receipt_schema"
                ):
                    runner._validate_terminal(
                        empty, authority_commit=AUTHORITY_COMMIT
                    )
                nested = deepcopy(builder())
                nested["replay_receipt"] = {
                    "predecessor": dict(runner.expected_replay_receipt()),
                    "historical_journal_sha256": authority.HISTORICAL_JOURNAL_SHA256,
                }
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "replay_receipt_schema"
                ):
                    runner._validate_terminal(
                        nested, authority_commit=AUTHORITY_COMMIT
                    )
                extra = deepcopy(builder())
                extra["replay_receipt"]["shadow_path_count"] = 7
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "replay_receipt_schema"
                ):
                    runner._validate_terminal(
                        extra, authority_commit=AUTHORITY_COMMIT
                    )
                wrong = deepcopy(builder())
                wrong["replay_receipt"]["member_key"] = "SSPRK3-2049"
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "replay_receipt_identity"
                ):
                    runner._validate_terminal(
                        wrong, authority_commit=AUTHORITY_COMMIT
                    )
                wrong_journal = deepcopy(builder())
                wrong_journal["replay_receipt"]["historical_journal_sha256"] = "0" * 64
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "replay_receipt_identity"
                ):
                    runner._validate_terminal(
                        wrong_journal, authority_commit=AUTHORITY_COMMIT
                    )

    def test_raw_output_rejects_nan_and_infinity(self) -> None:
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                target = root / authority.OUTPUT_NAMESPACE
                target.mkdir(parents=True)
                (target / "manifest.json").write_bytes(
                    runner._pretty(runner._manifest(AUTHORITY_COMMIT))
                )
                (target / "terminal.json").write_text(
                    '{\n  "classification": ' + token + "\n}\n"
                )
                with self.assertRaisesRegex(
                    runner.AC1RunnerError, "nonfinite_json_constant"
                ):
                    runner._inspect_output(
                        root, authority_commit=AUTHORITY_COMMIT
                    )

    def test_publication_rejects_foreign_or_racing_namespace_without_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            marker = target / "foreign"
            marker.write_text("keep")
            with self.assertRaisesRegex(runner.AC1RunnerError, "namespace_exists"):
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
                self.assertRaisesRegex(runner.AC1RunnerError, "namespace_arrived"),
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
            self.assertRaisesRegex(runner.AC1RunnerError, "campaign_store_mutated"),
        ):
            runner._publish_terminal(
                ROOT,
                {"kind": "manifest"},
                {"kind": "terminal"},
                store_before=(115, "a" * 64),
                store_after=(116, "b" * 64),
            )
        publish.assert_not_called()

    def test_shadow_premise_stop_is_not_rewrapped_as_implementation_failure(self) -> None:
        member = SimpleNamespace(
            temporal_ledger=object(),
            time=1.0,
            state=object(),
            operator=object(),
            projector=object(),
            transaction=object(),
            tracers=object(),
            initial=SimpleNamespace(grid=SimpleNamespace(coordinates=object())),
            step_index=4,
            transaction_serial=9,
        )
        fingerprint = {"descriptor_sha256": "a" * 64}
        restored = runner.ti2_runner.RestoredReplay(
            {"retry": 3, "attempted_width_hex": "0x1p-8"},
            member,
            fingerprint,
            "b" * 64,
        )
        stop = runner.tdg6.TDG6RefinementPathStop(
            path="outer_0", cause=RuntimeError("synthetic premise")
        )
        with (
            patch.object(
                runner.ti2_runner, "_member_fingerprint", return_value=fingerprint
            ),
            patch.object(
                runner.ti2_runner.tdg6,
                "prepare_tdg6_gr0_compositor",
                side_effect=stop,
            ),
            self.assertRaises(runner.tdg6.TDG6RefinementPathStop),
        ):
            runner.ti2_runner._prepare_shadow(restored)

    def test_premise_stop_after_restore_publishes_authenticated_receipt(self) -> None:
        fingerprint = {
            "member_key": authority.MEMBER_KEY,
            "method_label": "RK4",
            "accepted_time_hex": authority.ACCEPTED_TIME_HEX,
            "state_sha256": authority.PHYSICAL_STATE_SHA256,
            "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
            "transaction_sha256": authority.TRANSACTION_SHA256,
            "tracer_history_sha256": "c" * 64,
        }
        restored = runner.ti2_runner.RestoredReplay(
            {"retry": 3, "attempted_width_hex": authority.WIDTH_HEX},
            SimpleNamespace(),
            fingerprint,
            authority.HISTORICAL_JOURNAL_SHA256,
        )
        stop = runner.tdg6.TDG6RefinementPathStop(
            path="outer_0", cause=RuntimeError("synthetic premise")
        )
        sealed = (
            authority.SEALED_STORE_LEAF_COUNT,
            authority.SEALED_STORE_SNAPSHOT_SHA256,
        )
        published: dict[str, object] = {}

        def publish(
            _root: Path,
            _manifest: object,
            terminal: dict[str, object],
            **_kwargs: object,
        ) -> None:
            published["terminal"] = terminal

        with (
            patch.object(
                runner,
                "_execution_authority",
                return_value=authority.AC1Authority(
                    authority_commit=AUTHORITY_COMMIT
                ),
            ),
            patch.object(runner.authority, "require_output_absent"),
            patch.object(runner, "_snapshot_store", return_value=sealed),
            patch.object(runner, "_restore_replay", return_value=restored),
            patch.object(runner, "_prepare_shadow", side_effect=stop),
            patch.object(runner, "_publish_terminal", side_effect=publish),
        ):
            terminal = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        expected = runner.expected_replay_receipt()
        self.assertEqual(terminal["classification"], "shadow_proposal_premise_stop")
        self.assertEqual(terminal["replay_receipt"], expected)
        self.assertEqual(published["terminal"]["replay_receipt"], expected)
        self.assertNotIn("predecessor", terminal["replay_receipt"])


if __name__ == "__main__":
    unittest.main()
