"""Focused synthetic controls for the bounded TDG10-QA1 runner."""

from __future__ import annotations

import ast
import errno
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState
from recursive_horizons.fgc.evolution.hlt16_state_store import decode_payload
from recursive_horizons.fgc.evolution.proto5_runtime import GR0RuntimeMonitorState
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (
    ExactCompleteCAdmissionEvidence,
    ExactCompleteCClosedEvidence,
    ExactCompleteCDifferenceEvidence,
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    INDEPENDENT_LOCALIZER_ID,
    PRIMARY_LOCALIZER_ID,
    _combined_coefficient_hash,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime import (
    TDG10ExactCompleteCChannelEvidence,
    TDG10ExactCompleteCRuntimeClosed,
    TDG10ExactCompleteCRuntimeEvidence,
)


from scripts import run_fgc_tdg10_qa1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg10_qa1_authority as authority  # noqa: E402


AUTHORITY_COMMIT = "a" * 40


class FakeTracers:
    def __init__(self) -> None:
        self.labels = np.array([0.0, 1.0], dtype=np.float64)
        self.positions = np.array([0.25, 0.5], dtype=np.float64)
        self.proper_times = np.array([0.0, 0.0], dtype=np.float64)
        self.event_proper_times = [np.array([0.0], dtype=np.float64)]
        self.event_fields = [np.zeros((1, 2), dtype=np.float64)]
        self.cutoff = 1.0
        self.outer_radius = 128.0
        self.commits = 0

    def preview_advance(self, *_args: object, **_kwargs: object) -> tuple[object, ...]:
        return ("preview",)

    def commit_advance(self, *_payload: object) -> None:
        self.commits += 1
        self.positions = np.asarray(self.positions, dtype=np.float64) + 0.125
        self.proper_times = np.asarray(self.proper_times, dtype=np.float64) + 0.03125


class _Attempt:
    def __init__(self, final_time: float) -> None:
        self.preaccept_payload = (final_time,)
        self.proposal = SimpleNamespace(final_time=final_time)


def _fingerprint() -> dict[str, object]:
    return {
        "member_key": authority.MEMBER_KEY,
        "accepted_time_hex": authority.ACCEPTED_TIME_HEX,
        "state_sha256": authority.PHYSICAL_STATE_SHA256,
        "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": authority.TRANSACTION_SHA256,
    }


def _member() -> SimpleNamespace:
    return SimpleNamespace(
        tracers=FakeTracers(),
        initial=SimpleNamespace(
            grid=SimpleNamespace(coordinates=np.array([0.0, 1.0], dtype=np.float64))
        ),
    )


def _prepared(member: SimpleNamespace) -> SimpleNamespace:
    attempts = tuple(_Attempt(1.0 + 0.25 * (index + 1)) for index in range(4))
    clone = runner.tdg6._clone_tracers(member.tracers)
    for attempt in attempts:
        clone.commit_advance(*attempt.preaccept_payload)
    snapshot = runner.tdg6._tracer_snapshot(clone)
    fine = SimpleNamespace(
        attempts=attempts,
        initial_time=1.0,
        final_time=2.0,
        final_tracer_snapshot=snapshot,
        final_tracer_positions=np.asarray(clone.positions, dtype=np.float64).copy(),
        final_tracer_proper_times=np.asarray(clone.proper_times, dtype=np.float64).copy(),
        final_monitor_state=GR0RuntimeMonitorState(
            accepted_stage_count=4,
            last_transaction_serial=40,
            last_accepted_time=2.0,
        ),
        final_causal_state=CausalBudgetState(
            accepted_time=2.0,
            accumulated_characteristic_distance=0.5,
            previous_speed_upper=1.0,
        ),
        final_accepted=SimpleNamespace(
            time=2.0,
            step_index=12,
            transaction_serial=40,
            state=SimpleNamespace(
                u=np.array([1.0], dtype=np.float64),
                p=np.array([2.0], dtype=np.float64),
                q=np.array([3.0], dtype=np.float64),
            ),
        ),
    )
    return SimpleNamespace(fine=fine)


def _restored(member: SimpleNamespace) -> runner.ti2_runner.RestoredReplay:
    return runner.ti2_runner.RestoredReplay(
        dict(authority.REPLAY),
        member,
        _fingerprint(),
        authority.HISTORICAL_JOURNAL_SHA256,
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


def _receipt() -> authority.QA1Authority:
    return authority.QA1Authority(authority_commit=AUTHORITY_COMMIT)


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


class QA1RunnerTests(unittest.TestCase):
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
        self.assertNotIn("commit_tdg6_gr0_compositor", imported)
        self.assertNotIn("require_tdg6_temporal_admission", imported)
        self.assertIn("encode_payload", imported)
        self.assertIn("_restore_replay", source)
        self.assertIn("_prepare_shadow", source)
        self.assertIn("EXACT_COMPLETE_C_RUNTIME_MODULE", source)
        self.assertIn("assess_tdg10_exact_complete_c", source)
        for forbidden in (
            "commit_tdg6_gr0_compositor(",
            "require_tdg6_temporal_admission(",
            "encode_member(",
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

    def test_restore_replay_refuses_generation_10(self) -> None:
        replay = dict(authority.REPLAY)
        replay["predecessor_generation"] = 10
        with (
            patch.object(authority, "REPLAY", replay),
            self.assertRaisesRegex(
                runner.QA1RunnerError, "retry3_must_not_restore_generation_10"
            ),
        ):
            runner._restore_replay(ROOT)
        self.assertEqual(int(authority.REPLAY["predecessor_generation"]), 9)

    def test_adapter_channel_evidence_shape_serializes(self) -> None:
        channels = []
        for name in authority.CHANNEL_ORDER:
            channels.append(
                SimpleNamespace(
                    channel=name,
                    evidence=SimpleNamespace(
                        decision=SimpleNamespace(
                            admission_passed=True,
                            classification="resolved_order_pass",
                        ),
                        d12=SimpleNamespace(upper=Fraction(1, 8)),
                    ),
                )
            )
        serialized = runner.serialize_exact_assessment(
            SimpleNamespace(
                channel_evidence=channels,
                complete_admission_passed=True,
                failed_channels=(),
            )
        )
        self.assertTrue(serialized["complete_admission_passed"])
        self.assertEqual(serialized["failed_channels"], [])
        self.assertEqual(len(serialized["channels"]), 18)

    def test_terminal_reduction_is_closed(self) -> None:
        self.assertEqual(
            runner.reduce_qa1_terminal(
                complete_admission_passed=True, failed_channels=()
            ),
            runner.PASS_CLASS,
        )
        self.assertEqual(
            runner.reduce_qa1_terminal(
                complete_admission_passed=False, failed_channels=("u:R",)
            ),
            runner.NONPASS_CLASS,
        )

    def test_reconstructed_tracer_and_payload_hashes(self) -> None:
        member = _member()
        prepared = _prepared(member)
        restored = _restored(member)
        exact = runner.serialize_exact_assessment(_assessment(passing=True))
        original_positions = np.array(member.tracers.positions, copy=True)
        with patch.object(
            runner.ti2_runner, "_member_fingerprint", return_value=restored.fingerprint
        ):
            extras, record = runner.serialize_diagnostic_endpoint(
                restored, prepared, exact
            )
        self.assertEqual(member.tracers.commits, 0)
        self.assertTrue(np.array_equal(member.tracers.positions, original_positions))
        payload_relative = record["payload_relative"]
        descriptor_relative = record["descriptor_relative"]
        self.assertEqual(set(extras), {payload_relative, descriptor_relative})
        raw = extras[payload_relative]
        self.assertEqual(sha256(raw).hexdigest(), record["payload_raw_sha256"])
        self.assertTrue(payload_relative.startswith("payloads/"))
        self.assertTrue(payload_relative.endswith(".npz"))
        self.assertTrue(descriptor_relative.startswith("fine-endpoint/"))
        arrays, semantic, manifest = decode_payload(
            raw, expected_semantic_sha256=record["payload_semantic_sha256"]
        )
        self.assertEqual(semantic, record["payload_semantic_sha256"])
        self.assertEqual([item["name"] for item in manifest], list(arrays))
        self.assertTrue(np.array_equal(arrays["u"], np.array([1.0], dtype=np.float64)))
        descriptor = runner._parse_output_json(
            extras[descriptor_relative], name=descriptor_relative
        )
        self.assertEqual(descriptor["descriptor_sha256"], record["descriptor_sha256"])
        self.assertFalse(descriptor["accepted_state"])
        self.assertFalse(descriptor["campaign_descriptor"])
        self.assertFalse(descriptor["proto15_cursor"])
        self.assertFalse(descriptor["journal"])
        self.assertFalse(descriptor["checkpoint"])
        self.assertFalse(descriptor["ledger"])
        self.assertEqual(descriptor["accepted_step_tableau"], "SSPRK3")
        self.assertEqual(
            descriptor["actual_spatial_operator"], "inherited_RK4_2049_SBP4"
        )
        self.assertEqual(len(descriptor["exact_D12_upper_vector"]), 18)
        self.assertFalse(descriptor["state_advance_authorized"])

    def _run_mocked(
        self,
        root: Path,
        *,
        assess: object,
        prepare: object | None = None,
    ) -> dict[str, object]:
        member = _member()
        restored = _restored(member)
        prepared = _prepared(member)
        assess_kw: dict[str, object]
        if isinstance(assess, BaseException):
            assess_kw = {"side_effect": assess}
        else:
            assess_kw = {"return_value": assess}
        prepare_kw: dict[str, object]
        if prepare is None:
            prepare_kw = {"return_value": prepared}
        else:
            prepare_kw = {"side_effect": prepare}
        with (
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_snapshot_store", return_value=_sealed()),
            patch.object(runner, "_restore_replay", return_value=restored),
            patch.object(runner, "_prepare_shadow", **prepare_kw),
            patch.object(runner, "_assess_prepared_exact_complete_c", **assess_kw),
            patch.object(
                runner.ti2_runner,
                "_member_fingerprint",
                return_value=restored.fingerprint,
            ),
        ):
            return runner.run(root, authority_commit=AUTHORITY_COMMIT)

    def test_all_pass_creates_diagnostic_payload_but_no_campaign_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=_assessment(passing=True))
            namespace = root / authority.OUTPUT_NAMESPACE
            names = tuple(sorted(path.name for path in namespace.iterdir()))
            self.assertEqual(
                names,
                ("fine-endpoint", "manifest.json", "payloads", "terminal.json"),
            )
            self.assertEqual(terminal["classification"], runner.PASS_CLASS)
            self.assertTrue(terminal["diagnostic_fine_endpoint"]["present"])
            self.assertFalse(terminal["diagnostic_fine_endpoint"]["accepted_state"])
            payload = namespace / terminal["diagnostic_fine_endpoint"]["payload_relative"]
            descriptor = (
                namespace / terminal["diagnostic_fine_endpoint"]["descriptor_relative"]
            )
            self.assertTrue(payload.is_file())
            self.assertTrue(descriptor.is_file())
            self.assertFalse((namespace / "checkpoints").exists())
            self.assertFalse((namespace / "journal").exists())
            self.assertFalse((namespace / "members").exists())
            self.assertFalse(terminal["state_advance_authorized"])
            self.assertFalse(terminal["campaign_state_write_authorized"])
            self.assertFalse(terminal["fine_path_committed"])
            self.assertTrue(terminal["store_unchanged"])
            observed = runner._inspect_output(
                root, authority_commit=AUTHORITY_COMMIT
            )
            self.assertEqual(observed["classification"], runner.PASS_CLASS)
            self.assertEqual(
                observed["hashes"][terminal["diagnostic_fine_endpoint"]["payload_relative"]],
                terminal["diagnostic_fine_endpoint"]["payload_raw_sha256"],
            )

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
                self.assertFalse((namespace / "payloads").exists())
                self.assertFalse((namespace / "fine-endpoint").exists())
                self.assertNotIn("diagnostic_fine_endpoint", terminal)
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_invalid_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(runner, "_execution_authority", return_value=_receipt()),
                patch.object(runner, "_snapshot_store", return_value=_sealed()),
                patch.object(
                    runner,
                    "_restore_replay",
                    side_effect=runner.QA1RunnerError("replay", "broken", "synthetic"),
                ),
                self.assertRaisesRegex(runner.QA1RunnerError, "broken"),
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
            with self.assertRaisesRegex(runner.QA1RunnerError, "partial_or_foreign"):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_publication_rejects_foreign_or_racing_namespace_without_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            marker = target / "foreign"
            marker.write_text("keep")
            with self.assertRaisesRegex(runner.QA1RunnerError, "namespace_exists"):
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
                self.assertRaisesRegex(runner.QA1RunnerError, "namespace_arrived"),
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
            self.assertRaisesRegex(runner.QA1RunnerError, "campaign_store_mutated"),
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
        terminal = {
            **runner._terminal_base(
                AUTHORITY_COMMIT, runner.PREMISE_STOP_CLASS, _sealed(), _sealed()
            ),
            "typed_stop": {"type": "TDG6RefinementPathStop", "detail": "x"},
            "replay_receipt": runner.expected_replay_receipt(),
        }
        terminal["classification"] = runner.INVALID_CLASS
        with self.assertRaisesRegex(
            runner.QA1RunnerError, "terminal_class_not_publishable"
        ):
            runner._validate_terminal(terminal, authority_commit=AUTHORITY_COMMIT)

    def test_serialize_exact_assessment_consumes_runtime_evidence_nesting(self) -> None:
        digest_coeff = "a" * 64
        digest_keys = "b" * 64
        digest_stationary = "c" * 64
        d01 = ExactCompleteCDifferenceEvidence(
            lower=Fraction(0),
            upper=Fraction(0),
            polynomial_count=2,
            candidate_count=1,
            co_maximizer_count=1,
            localization_classification="unique_maximum",
            coefficient_stream_sha256=digest_coeff,
            survivor_key_stream_sha256=digest_keys,
            primary_stationary_count_stream_sha256=digest_stationary,
            independent_stationary_count_stream_sha256=digest_stationary,
            primary_evaluator_id=PRIMARY_LOCALIZER_ID,
            independent_evaluator_id=INDEPENDENT_LOCALIZER_ID,
            maximum_candidates=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        )
        d12 = ExactCompleteCDifferenceEvidence(
            lower=Fraction(0),
            upper=Fraction(0),
            polynomial_count=4,
            candidate_count=1,
            co_maximizer_count=1,
            localization_classification="unique_maximum",
            coefficient_stream_sha256=digest_coeff,
            survivor_key_stream_sha256=digest_keys,
            primary_stationary_count_stream_sha256=digest_stationary,
            independent_stationary_count_stream_sha256=digest_stationary,
            primary_evaluator_id=PRIMARY_LOCALIZER_ID,
            independent_evaluator_id=INDEPENDENT_LOCALIZER_ID,
            maximum_candidates=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        )
        zero = CertifiedMagnitudeInterval(Fraction(0), Fraction(0))
        decision = classify_tdg6_channel(zero, zero)
        self.assertEqual(decision.classification, "exact_zero")
        self.assertTrue(decision.admission_passed)
        evidence = ExactCompleteCAdmissionEvidence(
            d01=d01,
            d12=d12,
            decision=decision,
            row_count=1,
            row_stream_sha256="d" * 64,
            combined_coefficient_stream_sha256=_combined_coefficient_hash(
                d01.coefficient_stream_sha256, d12.coefficient_stream_sha256
            ),
            sufficient_pass_left=Fraction(0),
            sufficient_pass_right=Fraction(0),
            sufficient_contraction_pass=False,
            sufficient_contraction_failure=False,
            threshold_inconclusive=False,
            maximum_candidates_D01=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
            maximum_candidates_D12=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
            refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        )
        runtime = TDG10ExactCompleteCRuntimeEvidence(
            method=COMPARATOR_METHOD,
            owned_row_count=1,
            channel_evidence=tuple(
                TDG10ExactCompleteCChannelEvidence(channel=name, evidence=evidence)
                for name in authority.CHANNEL_ORDER
            ),
            complete_admission_passed=True,
            failed_channels=(),
            finest_pair_exact_upper_vector=tuple(
                Fraction(0) for _ in authority.CHANNEL_ORDER
            ),
            maximum_candidates_D01=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
            maximum_candidates_D12=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
            refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        )
        serialized = runner.serialize_exact_assessment(runtime)
        self.assertTrue(serialized["complete_admission_passed"])
        self.assertEqual(serialized["failed_channels"], [])
        self.assertEqual(serialized["channel_count"], 18)
        self.assertEqual(
            [item["channel"] for item in serialized["channels"]],
            list(authority.CHANNEL_ORDER),
        )
        self.assertTrue(all(item["admission_passed"] is True for item in serialized["channels"]))
        self.assertEqual(
            {item["classification"] for item in serialized["channels"]},
            {"exact_zero"},
        )
        self.assertEqual(
            serialized["channels"][0]["d12_upper"],
            {"numerator": "0", "denominator": "1"},
        )

    def test_status_rejects_rewritten_payload_or_descriptor_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._run_mocked(root, assess=_assessment(passing=True))
            terminal_path = root / authority.OUTPUT_NAMESPACE / "terminal.json"
            terminal = json.loads(terminal_path.read_text())
            del terminal["diagnostic_fine_endpoint"]["payload_relative"]
            terminal_path.write_bytes(runner._pretty(terminal))
            with self.assertRaisesRegex(
                runner.QA1RunnerError, "diagnostic fine endpoint"
            ):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=_assessment(passing=True))
            namespace = root / authority.OUTPUT_NAMESPACE
            payload_rel = str(terminal["diagnostic_fine_endpoint"]["payload_relative"])
            descriptor_rel = str(
                terminal["diagnostic_fine_endpoint"]["descriptor_relative"]
            )
            payload_path = namespace / payload_rel
            descriptor_path = namespace / descriptor_rel
            terminal_path = namespace / "terminal.json"

            tampered = b"not-a-valid-hlt16-payload"
            payload_path.write_bytes(tampered)
            rewritten = json.loads(terminal_path.read_text())
            rewritten["diagnostic_fine_endpoint"]["payload_raw_sha256"] = sha256(
                tampered
            ).hexdigest()
            terminal_path.write_bytes(runner._pretty(rewritten))
            with self.assertRaisesRegex(runner.QA1RunnerError, "payload_decode"):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            terminal = self._run_mocked(root, assess=_assessment(passing=True))
            namespace = root / authority.OUTPUT_NAMESPACE
            descriptor_rel = str(
                terminal["diagnostic_fine_endpoint"]["descriptor_relative"]
            )
            descriptor_path = namespace / descriptor_rel
            descriptor = json.loads(descriptor_path.read_text())
            descriptor["predecessor"]["generation"] = 10
            descriptor_path.write_bytes(runner._pretty(descriptor))
            with self.assertRaisesRegex(runner.QA1RunnerError, "descriptor_address"):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)


if __name__ == "__main__":
    unittest.main()
