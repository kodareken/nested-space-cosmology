"""Focused synthetic controls for the bounded TDG9 TI2 runner."""

from __future__ import annotations

import ast
from copy import deepcopy
import ctypes
from decimal import Decimal
import errno
from fractions import Fraction
import math
import os
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ti2 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ti2_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD  # noqa: E402


AUTHORITY_COMMIT = "a" * 40


def _occurrence(*, passed: bool = False, failed: bool = False, inconclusive: bool = False) -> dict[str, object]:
    return {
        "original_ownership_class": "both_components_independently_fail",
        "shadow_ownership_class": "complete_not_failure" if passed else "both_components_independently_fail",
        "complete_failure_cleared": passed,
        "shadow_complete_failure": failed,
        "shadow_complete_inconclusive": inconclusive,
    }


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
        "replay_receipts": {},
        "occurrences": [],
    }


def _publish_premise(root: Path, *, terminal: dict[str, object] | None = None) -> None:
    runner._publish(
        root,
        runner._manifest(AUTHORITY_COMMIT),
        _premise_terminal() if terminal is None else terminal,
    )


class TI2RunnerTests(unittest.TestCase):
    def test_runner_contains_no_admission_commit_or_candidate_path(self) -> None:
        source = (ROOT / authority.RUNNER_PATH).read_text()
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("candidate" in name.lower() for name in imported))
        for forbidden in (
            "require_tdg6_temporal_admission(",
            "commit_tdg6",
            "accept_step(",
            "acquire_writer",
            "publish_checkpoint",
        ):
            self.assertNotIn(forbidden, source)

    def test_prepare_shadow_changes_only_tableau_selector(self) -> None:
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
        stages = tuple(
            SimpleNamespace(stage_name=name)
            for name in ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint")
        )
        proposal = SimpleNamespace(method=COMPARATOR_METHOD, stages=stages)
        paths = tuple(
            SimpleNamespace(attempts=tuple(SimpleNamespace(proposal=proposal) for _ in range(count)))
            for count in (1, 2, 4)
        )
        prepared = SimpleNamespace(
            method=COMPARATOR_METHOD,
            initial_state_sha256=authority.PHYSICAL_STATE_SHA256,
            outer=paths[0],
            medium=paths[1],
            fine=paths[2],
        )
        restored = runner.RestoredReplay(
            replay={"retry": 3, "attempted_width_hex": "0x1.0000000000000p-8"},
            member=member,
            fingerprint=fingerprint,
            historical_journal_sha256="b" * 64,
        )
        with (
            patch.object(runner, "_member_fingerprint", return_value=fingerprint),
            patch.object(runner.tdg6, "prepare_tdg6_gr0_compositor", return_value=prepared) as prepare,
        ):
            self.assertIs(runner._prepare_shadow(restored), prepared)
        call = prepare.call_args.kwargs
        self.assertEqual(call["method"], COMPARATOR_METHOD)
        self.assertIs(call["state"], member.state)
        self.assertIs(call["rhs"], member.operator)
        self.assertIs(call["projector"], member.projector)
        self.assertIs(call["transaction"], member.transaction)
        self.assertIs(call["tracers"], member.tracers)
        self.assertIs(call["temporal_ledger"], member.temporal_ledger)

    def test_wrong_state_time_grid_operator_projector_tracer_or_ledger_fails(self) -> None:
        valid = {
            "member_key": authority.MEMBER_KEY,
            "method_label": "RK4",
            "source_integrator": runner.PRIMARY_METHOD,
            "point_count": authority.POINT_COUNT,
            "accepted_time_hex": authority.ACCEPTED_TIME_HEX,
            "state_sha256": authority.PHYSICAL_STATE_SHA256,
            "coordinates_sha256": authority.COORDINATES_SHA256,
            "grid_spacing_hex": authority.GRID_SPACING_HEX,
            "outer_radius_hex": authority.OUTER_RADIUS_HEX,
            "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
            "operator_class": "x.Proto12GR0EvolutionOperator",
            "projector_callable": (
                "x.make_gr0_center_boundary_projector.<locals>.projector"
            ),
            "transaction_class": "x.GR0RuntimeStageTransaction",
            "tracer_class": "x.NormalFlowTracers",
            "ledger_class": "x.TDG6TemporalLedger",
            "transaction_sha256": authority.TRANSACTION_SHA256_BY_RETRY[3],
            "tracer_history_sha256": "4" * 64,
        }
        runner._validate_one_variable_fingerprint(valid, retry=3)
        mutations = {
            "method_label": "SSPRK3",
            "accepted_time_hex": "0x1p+0",
            "state_sha256": "0" * 64,
            "grid_spacing_hex": "0x1p-3",
            "operator_class": "x.OtherOperator",
            "projector_callable": "x.other_projector",
            "transaction_class": "x.OtherTransaction",
            "tracer_class": "x.OtherTracers",
            "ledger_class": "x.OtherLedger",
        }
        for key, value in mutations.items():
            with self.subTest(key=key):
                changed = deepcopy(valid)
                changed[key] = value
                with self.assertRaisesRegex(
                    runner.TI1RunnerError, "one_variable_precondition"
                ):
                    runner._validate_one_variable_fingerprint(changed, retry=3)

    def test_fingerprint_encoder_is_exact_and_strictly_scoped(self) -> None:
        positive_zero = runner._fingerprint_exact(0.0)
        negative_zero = runner._fingerprint_exact(-0.0)
        adjacent = math.nextafter(1.0, 2.0)
        self.assertEqual(positive_zero, {"binary64_hex": "0x0.0p+0"})
        self.assertEqual(negative_zero, {"binary64_hex": "-0x0.0p+0"})
        self.assertNotEqual(positive_zero, negative_zero)
        self.assertEqual(
            runner._fingerprint_exact(adjacent),
            {"binary64_hex": adjacent.hex()},
        )
        self.assertEqual(runner._fingerprint_exact(True), True)
        self.assertEqual(runner._fingerprint_exact(1), 1)
        self.assertEqual(
            runner._fingerprint_exact(Fraction(3, 7)),
            {"numerator": "3", "denominator": "7"},
        )
        first = runner._canonical(runner._fingerprint_exact({"b": 2, "a": 1.0}))
        second = runner._canonical(runner._fingerprint_exact({"a": 1.0, "b": 2}))
        self.assertEqual(first, second)
        self.assertEqual(
            runner._fingerprint_exact((1, 2.0)),
            runner._fingerprint_exact([1, 2.0]),
        )
        for value in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(value=value), self.assertRaisesRegex(
                runner.TI1RunnerError, "nonfinite_fingerprint_binary64"
            ):
                runner._fingerprint_exact(value)
        for value in (np.float64(1.0), Decimal("1.0")):
            with self.subTest(kind=type(value)), self.assertRaisesRegex(
                runner.TI1RunnerError, "unsupported_fingerprint_value"
            ):
                runner._fingerprint_exact(value)
        with self.assertRaisesRegex(runner.TI1RunnerError, "nonexact_value"):
            runner._json_exact(1.0)

    def test_fingerprint_mutation_and_indices_fail_closed(self) -> None:
        left = runner._canonical(runner._fingerprint_exact({"value": 1.0}))
        right = runner._canonical(
            runner._fingerprint_exact({"value": math.nextafter(1.0, 2.0)})
        )
        self.assertNotEqual(left, right)
        for value in (-1, True, 1.0, np.int64(1)):
            with self.subTest(value=value), self.assertRaisesRegex(
                runner.TI1RunnerError,
                "negative_or_noninteger_fingerprint_index",
            ):
                runner._fingerprint_index("index", value)
        with self.assertRaisesRegex(
            runner.TI1RunnerError, "negative_or_noninteger_fingerprint_index"
        ):
            runner._validate_one_variable_fingerprint({}, retry=-1)

    def test_wrong_shadow_method_fails_one_variable_identity(self) -> None:
        member = SimpleNamespace(
            temporal_ledger=object(), time=1.0, state=object(), operator=object(),
            projector=object(), transaction=object(), tracers=object(),
            initial=SimpleNamespace(grid=SimpleNamespace(coordinates=object())),
            step_index=4, transaction_serial=9,
        )
        fingerprint = {"descriptor_sha256": "a" * 64}
        bad = SimpleNamespace(
            method="RK4", initial_state_sha256=authority.PHYSICAL_STATE_SHA256,
            outer=SimpleNamespace(attempts=()), medium=SimpleNamespace(attempts=()),
            fine=SimpleNamespace(attempts=()),
        )
        restored = runner.RestoredReplay(
            {"retry": 3, "attempted_width_hex": "0x1p-8"}, member, fingerprint, "b" * 64
        )
        with (
            patch.object(runner, "_member_fingerprint", return_value=fingerprint),
            patch.object(runner.tdg6, "prepare_tdg6_gr0_compositor", return_value=bad),
            self.assertRaisesRegex(runner.TI1RunnerError, "one_variable_identity"),
        ):
            runner._prepare_shadow(restored)

    def test_terminal_reduction_distinguishes_clear_persist_and_inconclusive(self) -> None:
        clear, _ = runner._reduce_completed([_occurrence(passed=True) for _ in range(10)])
        persist, _ = runner._reduce_completed(
            [_occurrence(passed=True) for _ in range(9)] + [_occurrence(failed=True)]
        )
        inconclusive, _ = runner._reduce_completed(
            [_occurrence(passed=True) for _ in range(9)]
            + [_occurrence(inconclusive=True)]
        )
        self.assertEqual(
            clear, "completed_all_ten_complete_failures_clear_under_tableau_shadow"
        )
        self.assertEqual(
            persist,
            "completed_one_or_more_complete_failures_persist_under_tableau_shadow",
        )
        self.assertEqual(inconclusive, "bounded_localization_or_resource_inconclusive")

    def test_changed_ownership_alone_does_not_clear_failure(self) -> None:
        items = [_occurrence(failed=True) for _ in range(10)]
        items[0]["shadow_ownership_class"] = "endpoint_state_owned_failure"
        classification, matrix = runner._reduce_completed(items)
        self.assertIn("one_or_more", classification)
        self.assertEqual(sum(matrix.values()), 10)

    def test_primary_v2_disagreement_is_invalid_not_a_tableau_result(self) -> None:
        evidence = SimpleNamespace(
            candidate_count=1,
            stationary_count_stream_sha256="a" * 64,
        )
        with (
            patch.object(
                runner.loc2,
                "_cubics",
                return_value=([object()] * 4088, [object()] * 4088),
            ),
            patch.object(
                runner.primary, "localize_absolute_maximum", return_value=evidence
            ),
            patch.object(
                runner.independent,
                "localize_absolute_maximum_independently_v2",
                return_value=evidence,
            ),
            patch.object(
                runner.loc2,
                "_primary_stationary_count_digest",
                return_value="a" * 64,
            ),
            patch.object(
                runner.loc2,
                "_require_route_agreement",
                side_effect=runner.loc2.LOC2RunnerError(
                    "localization", "independent_disagreement", "synthetic"
                ),
            ),
            self.assertRaisesRegex(
                runner.TI1RunnerError, "primary_v2_disagreement"
            ),
        ):
            runner._localize_shadow_occurrence((), 3, "u:alpha")

    def test_candidate_ceiling_is_bounded_inconclusive(self) -> None:
        with (
            patch.object(
                runner.loc2,
                "_cubics",
                return_value=([object()] * 4088, [object()] * 4088),
            ),
            patch.object(
                runner.primary,
                "localize_absolute_maximum",
                side_effect=RuntimeError("candidate_ceiling_exhausted"),
            ),
            self.assertRaisesRegex(
                runner.TI1BoundedInconclusive, "candidate_ceiling_exhausted"
            ),
        ):
            runner._localize_shadow_occurrence((), 3, "u:alpha")

    def test_shadow_premise_stop_is_not_rewrapped_as_implementation_failure(self) -> None:
        member = SimpleNamespace(
            temporal_ledger=object(), time=1.0, state=object(), operator=object(),
            projector=object(), transaction=object(), tracers=object(),
            initial=SimpleNamespace(grid=SimpleNamespace(coordinates=object())),
            step_index=4, transaction_serial=9,
        )
        fingerprint = {"descriptor_sha256": "a" * 64}
        restored = runner.RestoredReplay(
            {"retry": 3, "attempted_width_hex": "0x1p-8"}, member, fingerprint, "b" * 64
        )
        stop = runner.tdg6.TDG6RefinementPathStop(
            path="outer_0", cause=RuntimeError("synthetic premise")
        )
        with (
            patch.object(runner, "_member_fingerprint", return_value=fingerprint),
            patch.object(
                runner.tdg6, "prepare_tdg6_gr0_compositor", side_effect=stop
            ),
            self.assertRaises(runner.tdg6.TDG6RefinementPathStop),
        ):
            runner._prepare_shadow(restored)

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
            with self.assertRaisesRegex(runner.TI1RunnerError, "partial_or_foreign"):
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
            ("independent_method_agreement", True),
            ("store_unchanged", False),
            (
                "store_snapshot_after",
                {"leaf_count": 116, "sha256": "b" * 64},
            ),
        ):
            changed = deepcopy(terminal)
            changed[name] = value
            mutations[f"wrong-{name}"] = (changed, "terminal_identity")
        missing = deepcopy(terminal)
        del missing["candidate_branch_opened"]
        mutations["missing-nonclaim"] = (missing, "terminal_identity")
        invalid = deepcopy(terminal)
        invalid["classification"] = "invalid_provenance_or_implementation"
        mutations["invalid-class"] = (invalid, "terminal_class_not_publishable")

        for label, (changed, error) in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _publish_premise(root, terminal=changed)
                with self.assertRaisesRegex(runner.TI1RunnerError, error):
                    runner._inspect_output(
                        root, authority_commit=AUTHORITY_COMMIT
                    )

    def test_status_rejects_empty_or_wrong_manifest(self) -> None:
        for label, manifest in (
            ("empty", {}),
            (
                "wrong-authority",
                {**runner._manifest(AUTHORITY_COMMIT), "authority_commit": "b" * 40},
            ),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                runner._publish(root, manifest, _premise_terminal())
                with self.assertRaisesRegex(
                    runner.TI1RunnerError, "manifest_identity"
                ):
                    runner._inspect_output(
                        root, authority_commit=AUTHORITY_COMMIT
                    )

    def test_status_rejects_leaf_and_directory_substitution(self) -> None:
        original = runner._read_output_leaf
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _publish_premise(root)
            replaced = False

            def replace_leaf(
                directory_fd: int, name: str, before: os.stat_result
            ) -> bytes:
                nonlocal replaced
                raw = original(directory_fd, name, before)
                if name == "manifest.json" and not replaced:
                    replaced = True
                    terminal = runner._pretty(_premise_terminal())
                    os.unlink("terminal.json", dir_fd=directory_fd)
                    leaf = os.open(
                        "terminal.json",
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                        0o644,
                        dir_fd=directory_fd,
                    )
                    try:
                        os.write(leaf, terminal)
                    finally:
                        os.close(leaf)
                return raw

            with (
                patch.object(runner, "_read_output_leaf", side_effect=replace_leaf),
                self.assertRaisesRegex(
                    runner.TI1RunnerError,
                    "output_(leaf_substituted|directory_changed)",
                ),
            ):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _publish_premise(root)
            target = root / authority.OUTPUT_NAMESPACE
            backup = target.with_name(f"{target.name}.old")
            replaced = False

            def replace_directory(
                directory_fd: int, name: str, before: os.stat_result
            ) -> bytes:
                nonlocal replaced
                raw = original(directory_fd, name, before)
                if name == "manifest.json" and not replaced:
                    replaced = True
                    target.rename(backup)
                    target.mkdir()
                    for leaf_name in ("manifest.json", "terminal.json"):
                        (target / leaf_name).write_bytes((backup / leaf_name).read_bytes())
                return raw

            with (
                patch.object(
                    runner, "_read_output_leaf", side_effect=replace_directory
                ),
                self.assertRaisesRegex(
                    runner.TI1RunnerError, "output_directory_(changed|substituted)"
                ),
            ):
                runner._inspect_output(root, authority_commit=AUTHORITY_COMMIT)

    def test_publication_rejects_foreign_or_racing_namespace_without_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            marker = target / "foreign"
            marker.write_text("keep")
            with self.assertRaisesRegex(runner.TI1RunnerError, "namespace_exists"):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})
            self.assertEqual(marker.read_text(), "keep")

        class _RaceRename:
            argtypes = None
            restype = None

            def __call__(self, *_arguments: object) -> int:
                ctypes.set_errno(errno.EEXIST)
                return -1

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_library = SimpleNamespace(renameatx_np=_RaceRename())
            with (
                patch.object(runner.ctypes, "CDLL", return_value=fake_library),
                self.assertRaisesRegex(runner.TI1RunnerError, "namespace_arrived"),
            ):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            self.assertFalse((root / authority.OUTPUT_NAMESPACE).exists())
            self.assertFalse(
                any(item.name.startswith(authority.STAGING_PREFIX) for item in parent.iterdir())
            )

    def test_no_terminal_can_publish_after_store_drift(self) -> None:
        with (
            patch.object(runner, "_publish") as publish,
            self.assertRaisesRegex(runner.TI1RunnerError, "campaign_store_mutated"),
        ):
            runner._publish_terminal(
                ROOT,
                {"kind": "manifest"},
                {"kind": "terminal"},
                store_before=(115, "a" * 64),
                store_after=(116, "b" * 64),
            )
        publish.assert_not_called()

    def test_compact_LOC2_original_map_is_all_ten_failures(self) -> None:
        observed = runner._original_occurrences(ROOT)
        self.assertEqual(tuple(observed), authority.FAILED_OCCURRENCES)
        self.assertTrue(all(item["original_complete_failure"] for item in observed.values()))


if __name__ == "__main__":
    unittest.main()
