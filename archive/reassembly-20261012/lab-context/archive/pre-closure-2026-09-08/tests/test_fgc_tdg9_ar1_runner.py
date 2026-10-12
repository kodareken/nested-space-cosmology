"""Focused safety and classification controls for the TDG9 AR1 runner."""

from __future__ import annotations

import ast
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ar1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402


Q = Fraction


def _interval(lower: int, upper: int) -> SimpleNamespace:
    return SimpleNamespace(lower=Q(lower), upper=Q(upper))


def _evidence(
    classification: str,
    *,
    evaluator_id: str,
    row_count: int = 1,
) -> SimpleNamespace:
    outer = SimpleNamespace(
        interval=_interval(4, 4),
        coefficient_stream_sha256="1" * 64,
        polynomial_count=2,
        nodes_visited=2,
        maximum_depth_reached=0,
        evaluator_id=evaluator_id,
    )
    finest = SimpleNamespace(
        interval=_interval(1, 1),
        coefficient_stream_sha256="2" * 64,
        polynomial_count=4,
        nodes_visited=4,
        maximum_depth_reached=0,
        evaluator_id=evaluator_id,
    )
    return SimpleNamespace(
        outer_difference=outer,
        finest_difference=finest,
        row_count=row_count,
        classification=classification,
        sufficient_pass_certified=classification == "sufficient_contraction_pass",
        contraction_threshold_resolved=classification != "exact_zero",
        contraction_threshold_passed=(
            True
            if classification == "sufficient_contraction_pass"
            else False
            if classification == "sufficient_contraction_failure"
            else None
        ),
        sufficient_pass_left=Q(8),
        sufficient_pass_right=Q(16),
        max_depth=64,
        max_nodes=4_194_304,
        nodes_visited=6,
        combined_coefficient_stream_sha256="3" * 64,
        schema_version=1,
        evaluator_id=evaluator_id,
        absolute_tolerance_used=False,
        physical_signal_used_for_normalization=False,
        replay_identity_certified=False,
        historical_proposal_certified=False,
        pde_trajectory_order_certified=False,
    )


def _rows() -> tuple[object, ...]:
    segment = (0.0, 0.0, 0.0, 0.0, 0.5)
    return ((segment, (segment, segment), (segment,) * 4),)


def _terminal_rows(classification: str) -> list[dict[str, object]]:
    return [
        {"channels": [{"classification": classification} for _ in ar1.CHANNEL_ORDER]}
        for _ in ar1.REPLAYS
    ]


def _selection_fixture() -> tuple[object, object, object]:
    target = {
        "binary64_hex": ar1.TARGET_BINARY64_HEX,
        "rational": ar1.TARGET_RATIONAL,
    }
    descriptor_sha256 = "d" * 64
    cursor = {
        "protocol_artifact_id": ar1.TARGET_PROTOCOL,
        "campaign_id": ar1.CAMPAIGN_ID,
        "member_key": ar1.MEMBER_KEY,
        "method": ar1.METHOD,
        "point_count": ar1.POINT_COUNT,
        "committed_common_event_index": ar1.EVENT,
        "event_target_time": target,
        "accepted_state_sha256": descriptor_sha256,
        "mode": "RETRY_PENDING",
    }
    member_state = SimpleNamespace(
        descriptor_sha256=descriptor_sha256,
        cursor=cursor,
    )
    checkpoint = SimpleNamespace(
        protocol=ar1.TARGET_PROTOCOL,
        campaign_id=ar1.CAMPAIGN_ID,
        plan_sha256=ar1.CAMPAIGN_PLAN_SHA256,
        authorization_commit=ar1.CAMPAIGN_AUTHORIZATION_COMMIT,
        event=ar1.EVENT,
        target=target,
        members={ar1.MEMBER_KEY: member_state},
        generation=9,
    )
    descriptor = {
        "metadata": {
            "runtime_identity": {
                "protocol_artifact_id": ar1.TARGET_PROTOCOL,
                "campaign_id": ar1.CAMPAIGN_ID,
                "branch": ar1.BRANCH,
                "amplitude": ar1.AMPLITUDE,
                "member_key": ar1.MEMBER_KEY,
                "method": ar1.METHOD,
                "point_count": ar1.POINT_COUNT,
            },
            "template": {
                "amplitude": ar1.AMPLITUDE,
                "method": ar1.METHOD,
                "point_count": ar1.POINT_COUNT,
            },
        }
    }
    store = SimpleNamespace(
        load_state=lambda _digest: SimpleNamespace(descriptor=descriptor)
    )
    return store, checkpoint, member_state


class TDG9AR1RunnerTests(unittest.TestCase):
    def test_runner_imports_no_continuation_or_candidate_path_and_has_no_store_writes(self) -> None:
        path = ROOT / "scripts/run_fgc_tdg9_ar1.py"
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("event1" in name or "candidate" in name for name in imported))
        for forbidden in (
            "_advance_authenticated_event",
            "acquire_writer",
            "persist_snapshot",
            "publish_journal",
            "publish_checkpoint",
            "publish_terminal_lock",
            "SGB-L",
            "FGC-QR",
        ):
            self.assertNotIn(forbidden, source)

    def test_terminal_reduction_has_only_three_scientific_outcomes(self) -> None:
        self.assertEqual(
            runner._terminal_classification(  # type: ignore[attr-defined]
                _terminal_rows("sufficient_contraction_pass")
            ),
            "legacy_enclosure_owned_nonadmission_on_all_frozen_samples",
        )
        rows = _terminal_rows("sufficient_contraction_pass")
        rows[1]["channels"][4]["classification"] = "sufficient_contraction_failure"  # type: ignore[index]
        self.assertEqual(
            runner._terminal_classification(rows),  # type: ignore[attr-defined]
            "legacy_enclosure_not_sole_owner_on_frozen_samples",
        )
        rows = _terminal_rows("exact_zero")
        rows[2]["channels"][7]["classification"] = "resource_inconclusive"  # type: ignore[index]
        self.assertEqual(
            runner._terminal_classification(rows),  # type: ignore[attr-defined]
            "arithmetic_discriminator_inconclusive",
        )
        with self.assertRaises(runner.TDG9AR1RunnerError):
            runner._terminal_classification(  # type: ignore[attr-defined]
                _terminal_rows("unresolved_without_typed_resource_stop")
            )

        self.assertEqual(
            runner._terminal_classification(_terminal_rows("exact_zero")),  # type: ignore[attr-defined]
            "legacy_enclosure_owned_nonadmission_on_all_frozen_samples",
        )
        rows = _terminal_rows("resource_inconclusive")
        rows[2]["channels"][17]["classification"] = "sufficient_contraction_failure"  # type: ignore[index]
        self.assertEqual(
            runner._terminal_classification(rows),  # type: ignore[attr-defined]
            "legacy_enclosure_not_sole_owner_on_frozen_samples",
        )
        rows = _terminal_rows("sufficient_contraction_pass")
        rows.pop()
        with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "classification_shape"):
            runner._terminal_classification(rows)  # type: ignore[attr-defined]

    def test_checkpoint_cursor_and_descriptor_selection_are_exact(self) -> None:
        store, checkpoint, member_state = _selection_fixture()
        identity = runner._require_checkpoint_selection(  # type: ignore[attr-defined]
            store, checkpoint, member_state
        )
        self.assertEqual(identity["branch"], "GR-0")
        self.assertEqual(identity["point_count"], 2049)

        cases = (
            ("campaign_id", "wrong"),
            ("plan_sha256", "0" * 64),
            ("authorization_commit", "0" * 40),
            ("event", 24),
            (
                "target",
                {"binary64_hex": "0x1.0000000000000p+0", "rational": "1"},
            ),
        )
        for name, value in cases:
            store, checkpoint, member_state = _selection_fixture()
            setattr(checkpoint, name, value)
            with self.subTest(checkpoint=name), self.assertRaisesRegex(
                runner.TDG9AR1RunnerError, "checkpoint_selection"
            ):
                runner._require_checkpoint_selection(  # type: ignore[attr-defined]
                    store, checkpoint, member_state
                )

        for name, value in (
            ("campaign_id", "wrong"),
            ("method", "SSPRK3"),
            ("point_count", 4097),
            ("committed_common_event_index", 24),
            ("accepted_state_sha256", "0" * 64),
            ("mode", "FRESH_READY"),
        ):
            store, checkpoint, member_state = _selection_fixture()
            member_state.cursor[name] = value
            with self.subTest(cursor=name), self.assertRaisesRegex(
                runner.TDG9AR1RunnerError, "cursor_selection"
            ):
                runner._require_checkpoint_selection(  # type: ignore[attr-defined]
                    store, checkpoint, member_state
                )

        store, checkpoint, member_state = _selection_fixture()
        persisted = store.load_state(member_state.descriptor_sha256)
        persisted.descriptor["metadata"]["runtime_identity"]["branch"] = "FGC-QR"
        store.load_state = lambda _digest: persisted
        with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "descriptor_selection"):
            runner._require_checkpoint_selection(  # type: ignore[attr-defined]
                store, checkpoint, member_state
            )

    def test_owned_Hermite_arrays_and_width_are_strict_binary64(self) -> None:
        shape = (3, ar1.OWNED_ROW_COUNT, 6)
        valid = np.zeros(shape, dtype="<f8", order="C")
        self.assertIs(
            runner._require_owned_array(valid, label="valid"),  # type: ignore[attr-defined]
            valid,
        )
        invalid = (
            np.zeros(shape, dtype=np.float32),
            np.zeros(shape, dtype=">f8"),
            np.asfortranarray(valid),
            np.zeros(shape, dtype=object),
            np.zeros((3, ar1.OWNED_ROW_COUNT - 1, 6), dtype="<f8"),
            valid.tolist(),
        )
        for value in invalid:
            with self.subTest(type=type(value).__name__), self.assertRaises(
                runner.TDG9AR1RunnerError
            ):
                runner._require_owned_array(value, label="mutated")  # type: ignore[attr-defined]
        for nonfinite in (float("nan"), float("inf"), -float("inf")):
            mutated = valid.copy()
            mutated[2, -1, -1] = nonfinite
            with self.subTest(nonfinite=nonfinite), self.assertRaisesRegex(
                runner.TDG9AR1RunnerError, "nonfinite"
            ):
                runner._require_owned_array(mutated, label="mutated")  # type: ignore[attr-defined]

        proposal = SimpleNamespace(
            initial_state=valid,
            candidate_state=valid,
            initial_time=0.0,
            final_time=0.5,
        )
        with (
            patch.object(runner.tdg5, "_proposal_endpoint_records", return_value=(valid, valid)),
            patch.object(runner.tdg5, "_owned_state", side_effect=lambda value: value),
            patch.object(runner.tdg5, "_owned_rhs", side_effect=lambda value: value),
        ):
            surface = runner._segment_surface(proposal)  # type: ignore[attr-defined]
            self.assertEqual(surface.width.hex(), "0x1.0000000000000p-1")
            proposal.final_time = np.float64(0.5)
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "Hermite_width"):
                runner._segment_surface(proposal)  # type: ignore[attr-defined]
            for invalid_width in (
                0.0,
                -0.5,
                float("nan"),
                float("inf"),
                -float("inf"),
            ):
                proposal.final_time = invalid_width
                with self.subTest(width=invalid_width), self.assertRaisesRegex(
                    runner.TDG9AR1RunnerError, "Hermite_width"
                ):
                    runner._segment_surface(proposal)  # type: ignore[attr-defined]

    def test_replay_mismatch_is_rejected_before_exact_classification(self) -> None:
        replay = ar1.REPLAYS[0]
        historical = {"retry_count_for_current_macro_step": 3, "nested": [1, 2]}
        normalized = runner._require_replayed_evidence(  # type: ignore[attr-defined]
            {"retry_count_for_current_macro_step": 3, "nested": (1, 2)},
            historical,
            replay,
            rejected_retry_count=3,
        )
        self.assertEqual(normalized, historical)
        with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "historical_evidence_mismatch"):
            runner._require_replayed_evidence(  # type: ignore[attr-defined]
                {"retry_count_for_current_macro_step": 3, "nested": (1, 3)},
                historical,
                replay,
                rejected_retry_count=3,
            )
        with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "historical_evidence_mismatch"):
            runner._require_replayed_evidence(  # type: ignore[attr-defined]
                historical, historical, replay, rejected_retry_count=4
            )

    def test_two_exact_evaluators_must_agree(self) -> None:
        primary = _evidence(
            "sufficient_contraction_pass",
            evaluator_id="tdg9_exact_bernstein_de_casteljau_v1",
        )
        independent = _evidence(
            "sufficient_contraction_pass",
            evaluator_id="tdg9_exact_derivative_isolation_interval_horner_v1",
        )
        with (
            patch.object(
                runner.exact_primary,
                "assess_exact_temporal_refinement",
                return_value=primary,
            ),
            patch.object(
                runner.exact_independent,
                "assess_exact_temporal_refinement_independently",
                return_value=independent,
            ),
            patch.object(runner, "_json_exact", return_value={"evidence": "exact"}),
        ):
            result = runner._assess_channel(lambda: iter(_rows()))  # type: ignore[attr-defined]
        self.assertEqual(result["classification"], "sufficient_contraction_pass")

        independent.classification = "sufficient_contraction_failure"
        with (
            patch.object(
                runner.exact_primary,
                "assess_exact_temporal_refinement",
                return_value=primary,
            ),
            patch.object(
                runner.exact_independent,
                "assess_exact_temporal_refinement_independently",
                return_value=independent,
            ),
            patch.object(runner, "_json_exact", return_value={"evidence": "exact"}),
            self.assertRaisesRegex(runner.TDG9AR1RunnerError, "evaluator_disagreement"),
        ):
            runner._assess_channel(lambda: iter(_rows()))  # type: ignore[attr-defined]

    def test_output_namespace_is_exclusive_canonical_and_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.assertEqual(runner._inspect_output(root)["state"], "absent")  # type: ignore[attr-defined]
            manifest = {"schema": runner.RAW_RESULT_SCHEMA, "kind": "manifest"}
            terminal = {
                "schema": runner.RAW_RESULT_SCHEMA,
                "classification": "arithmetic_discriminator_inconclusive",
            }
            runner._publish(root, manifest, terminal)  # type: ignore[attr-defined]
            status = runner._inspect_output(root)  # type: ignore[attr-defined]
            self.assertEqual(status["state"], "terminal")
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "namespace_preexisting"):
                runner._publish(root, manifest, terminal)  # type: ignore[attr-defined]

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = root / Path(ar1.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{ar1.STAGING_PREFIX}fault-cut").mkdir()
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "private_stage_present"):
                runner._inspect_output(root)  # type: ignore[attr-defined]

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / ar1.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            (target / "manifest.json").write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "partial_or_foreign"):
                runner._inspect_output(root)  # type: ignore[attr-defined]

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = root / Path(ar1.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            destination = root / "foreign"
            destination.mkdir()
            (root / ar1.OUTPUT_NAMESPACE).symlink_to(destination, target_is_directory=True)
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "namespace_symlink"):
                runner._inspect_output(root)  # type: ignore[attr-defined]

    def test_run_enforces_store_identity_and_exact_108_call_budget(self) -> None:
        receipt = SimpleNamespace(
            authority_commit="a" * 40,
            pref2_commit=ar1.PREF2_COMMIT,
            evaluator_commit=ar1.EVALUATOR_COMMIT,
            store_manifest_sha256=ar1.PREF2_STORE_MANIFEST_SHA256,
            output_namespace=ar1.OUTPUT_NAMESPACE,
            channel_order=ar1.CHANNEL_ORDER,
            max_depth=ar1.MAX_DEPTH,
            max_nodes_per_channel_per_evaluator=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
            max_evaluator_calls=ar1.MAX_EVALUATOR_CALLS,
            environment=tuple(ar1.ENVIRONMENT.items()),
            selection=tuple(ar1.SELECTION.items()),
            authority_delta_paths=ar1.AUTHORITY_DELTA_PATHS,
        )
        shells = {
            key: object()
            for key in (
                "RK4-2049", "RK4-4097", "RK4-8193",
                "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
            )
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {"classification": "sufficient_contraction_pass"}
            with (
                patch.object(runner, "_fixed_root", return_value=root),
                patch.object(runner, "_execution_authority", return_value=receipt),
                patch.object(runner, "_inspect_output", return_value={"state": "absent"}),
                patch.object(runner, "_authenticate_pref2", side_effect=["same", "same"]),
                patch.object(runner, "HLT16CampaignStore", return_value=object()),
                patch.object(runner, "build_static_gr0_shells", return_value=shells),
                patch.object(runner, "_prepare_replay", return_value=(object(), {})) as prepare,
                patch.object(runner, "_hermite_surface", return_value=object()),
                patch.object(runner, "_rows_factory", return_value=(2044, lambda: iter(_rows()))),
                patch.object(runner, "_assess_channel", return_value=result) as assess,
                patch.object(runner, "_publish") as publish,
            ):
                terminal = runner.run(ROOT, authority_commit="a" * 40)
            self.assertEqual(terminal["evaluator_calls"], 108)
            self.assertEqual(prepare.call_count, 3)
            self.assertEqual(assess.call_count, 54)
            publish.assert_called_once()

            with (
                patch.object(runner, "_fixed_root", return_value=root),
                patch.object(runner, "_execution_authority", return_value=receipt),
                patch.object(runner, "_inspect_output", return_value={"state": "absent"}),
                patch.object(runner, "_authenticate_pref2", side_effect=["before", "after"]),
                patch.object(runner, "HLT16CampaignStore", return_value=object()),
                patch.object(runner, "build_static_gr0_shells", return_value=shells),
                patch.object(runner, "_prepare_replay", return_value=(object(), {})),
                patch.object(runner, "_hermite_surface", return_value=object()),
                patch.object(runner, "_rows_factory", return_value=(2044, lambda: iter(_rows()))),
                patch.object(runner, "_assess_channel", return_value=result),
                patch.object(runner, "_publish") as publish,
                self.assertRaisesRegex(runner.TDG9AR1RunnerError, "store_changed_during_replay"),
            ):
                runner.run(ROOT, authority_commit="a" * 40)
            publish.assert_not_called()

    def test_raw_manifest_binds_environment_selection_and_committed_image_delta(self) -> None:
        receipt = SimpleNamespace(
            authority_commit="a" * 40,
            pref2_commit=ar1.PREF2_COMMIT,
            evaluator_commit=ar1.EVALUATOR_COMMIT,
            environment=tuple(ar1.ENVIRONMENT.items()),
            selection=tuple(ar1.SELECTION.items()),
            authority_delta_paths=ar1.AUTHORITY_DELTA_PATHS,
            store_manifest_sha256=ar1.PREF2_STORE_MANIFEST_SHA256,
            output_namespace=ar1.OUTPUT_NAMESPACE,
            channel_order=ar1.CHANNEL_ORDER,
            max_depth=ar1.MAX_DEPTH,
            max_nodes_per_channel_per_evaluator=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
            max_evaluator_calls=ar1.MAX_EVALUATOR_CALLS,
        )
        manifest = runner._manifest(receipt)  # type: ignore[attr-defined]
        self.assertEqual(manifest["environment"], ar1.ENVIRONMENT)
        self.assertEqual(manifest["selection"], ar1.SELECTION)
        self.assertEqual(
            manifest["authority_delta_paths"], list(ar1.AUTHORITY_DELTA_PATHS)
        )
        self.assertEqual(manifest["authority_commit"], "a" * 40)
        self.assertNotIn("authority_commit", ar1.COMMITTED_IMAGE)

    def test_user_supplied_repository_and_authority_inputs_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(runner.TDG9AR1RunnerError, "repository_scope"):
                runner._fixed_root(Path(temporary))  # type: ignore[attr-defined]
        with self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "malformed"):
            ar1.authorize_diagnostic(ROOT, b"", b"", "not-a-commit")


if __name__ == "__main__":
    unittest.main()
