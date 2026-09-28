"""Raw/store/shadow/Git-blind controls for the QA2-PREF1 binder core."""

from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg10_qa2_pref1_binder as binder,
)


def _synthetic_exact(retry: int, failed: tuple[str, ...]) -> dict[str, object]:
    channels = []
    for ordinal, name in enumerate(binder.CHANNEL_ORDER, start=1):
        passed = name not in failed
        channels.append(
            {
                "admission_passed": passed,
                "channel": name,
                "classification": (
                    "resolved_order_pass" if passed else "resolved_order_failure"
                ),
                "d12_upper": binder._encode_rational(
                    Fraction(retry * 100 + ordinal, 1 << 60)
                ),
            }
        )
    return {
        "admission_is_all_of": True,
        "channel_count": 18,
        "channel_order": list(binder.CHANNEL_ORDER),
        "channels": channels,
        "complete_admission_passed": not failed,
        "failed_channels": list(failed),
        "independent_route_required": True,
        "interval_owner": binder.INTERVAL_OWNER,
    }


def _synthetic_width(width: dict[str, object]) -> dict[str, object]:
    failed = tuple(width["failed_channels"])
    return {
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "attempted_width_hex": width["attempted_width_hex"],
        "complete_admission_passed": not failed,
        "exact_complete_C": _synthetic_exact(int(width["retry"]), failed),
        "predecessor_generation": width["generation"],
        "replay_receipt": binder._expected_replay_receipt(width),
        "retry": width["retry"],
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "width_class": (
            "all_18_channel_pass" if not failed else "one_or_more_channel_nonpass"
        ),
    }


def _synthetic_terminal() -> dict[str, object]:
    return {
        **binder._expected_terminal_base(binder.RAW_CLASSIFICATION, robustness=False),
        "widths": [_synthetic_width(dict(width)) for width in binder.WIDTHS],
    }


def _synthetic_raw_bytes() -> tuple[bytes, bytes]:
    return (
        binder.canonical_result(binder._expected_manifest()),
        binder.canonical_result(_synthetic_terminal()),
    )


def _difference(
    *,
    level: str,
    lower: Fraction,
    upper: Fraction,
    coefficient: str,
    classification: str = "unique_maximum",
    co_maximizers: int = 1,
    stationary: str = "3" * 64,
) -> dict[str, object]:
    d01 = level == "D01"
    return {
        "candidate_count": max(1, co_maximizers),
        "co_maximizer_count": co_maximizers,
        "coefficient_stream_sha256": coefficient,
        "independent_evaluator_id": binder.INDEPENDENT_EVALUATOR_ID,
        "independent_stationary_count_stream_sha256": stationary,
        "localization_classification": classification,
        "lower": binder._encode_rational(lower),
        "maximum_candidates": (
            binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D01
            if d01
            else binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        ),
        "polynomial_count": (
            binder.D01_POLYNOMIAL_COUNT if d01 else binder.D12_POLYNOMIAL_COUNT
        ),
        "primary_evaluator_id": binder.PRIMARY_EVALUATOR_ID,
        "primary_stationary_count_stream_sha256": stationary,
        "routes_agree": True,
        "survivor_key_stream_sha256": "4" * 64,
        "upper": binder._encode_rational(upper),
    }


def _rich_channel(*, channel: str = "u:R", passing: bool = True) -> dict[str, object]:
    d01_lower = Fraction(4 if passing else 1)
    d01_upper = d01_lower
    d12_lower = Fraction(1)
    d12_upper = Fraction(1)
    d01_coefficient = "1" * 64
    d12_coefficient = "2" * 64
    decision = binder._reclassify_channel(
        channel=channel,
        d01_lower=d01_lower,
        d01_upper=d01_upper,
        d12_lower=d12_lower,
        d12_upper=d12_upper,
    )
    combined = sha256(
        binder._COMBINED_HASH_DOMAIN
        + f"{d01_coefficient}\n{d12_coefficient}\n".encode()
    ).hexdigest()
    classification = decision["decision"]
    return {
        "D01": _difference(
            level="D01",
            lower=d01_lower,
            upper=d01_upper,
            coefficient=d01_coefficient,
        ),
        "D12": _difference(
            level="D12",
            lower=d12_lower,
            upper=d12_upper,
            coefficient=d12_coefficient,
        ),
        "admission_passed": classification.admission_passed,
        "candidate_evidence_available": True,
        "channel": channel,
        "classification": classification.classification,
        "combined_coefficient_stream_sha256": combined,
        "evaluator_id": binder.TDG10_EVALUATOR_ID,
        "hash_evidence_available": True,
        "maximum_candidates_D01": binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "maximum_candidates_D12": binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "order_threshold_passed": classification.order_threshold_passed,
        "order_threshold_resolved": classification.order_threshold_resolved,
        "refinement_depth": binder.PRIMARY_REFINEMENT_DEPTH,
        "route_evidence_available": True,
        "row_count": binder.OWNED_ROW_COUNT,
        "row_stream_sha256": "5" * 64,
        "sufficient_condition": binder.SUFFICIENT_CONDITION,
        "sufficient_condition_holds": decision["sufficient_condition_holds"],
        "sufficient_contraction_failure": decision["sufficient_contraction_failure"],
        "sufficient_contraction_pass": decision["sufficient_contraction_pass"],
        "sufficient_pass_left": binder._encode_rational(decision["left"]),
        "sufficient_pass_right": binder._encode_rational(decision["right"]),
        "temporal_retry_permitted": classification.temporal_retry_permitted,
        "threshold_inconclusive": decision["threshold_inconclusive"],
    }


def _load_reproducer():
    path = ROOT / "scripts/reproduce_fgc_tdg10_qa2_pref1.py"
    spec = importlib.util.spec_from_file_location("qa2_pref1_reproducer", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TDG10QA2PREF1BinderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # Tracked compact/config input only. No ignored raw or historical store.
        cls.config_raw = (ROOT / binder.CONFIG_PATH).read_bytes()
        cls.result_raw = (ROOT / binder.RESULT_PATH).read_bytes()
        cls.result = json.loads(cls.result_raw)

    def assert_semantic_mutation_rejected(self, mutation) -> None:
        changed = deepcopy(self.result)
        mutation(changed["artifact_payload"])
        raw = binder.canonical_result(changed)
        with patch.object(binder, "COMPACT_RESULT_SHA256", None):
            try:
                binder.validate_compact_result(self.config_raw, raw)
            except binder.TDG10QA2PREF1Error as exc:
                self.assertNotEqual(exc.stop_id, "PREF1_COMPACT_HASH_DRIFT")
            else:
                self.fail("semantic compact mutation was accepted")

    def test_config_is_exact_two_width_environment_bound_contract(self) -> None:
        parsed = binder._config(self.config_raw)
        self.assertEqual(parsed, binder.expected_config())
        self.assertEqual(parsed["authority"]["blob_count"], 16)
        self.assertEqual(parsed["live_dependency_authority"]["blob_count"], 10)
        self.assertEqual(parsed["selection"]["selected_retries"], [4, 5])
        self.assertEqual(
            [item["generation"] for item in parsed["selection"]["widths"]],
            [10, 11],
        )
        self.assertEqual(parsed["work_budget"]["shadow_paths"], 14)
        self.assertEqual(
            parsed["work_budget"]["maximum_stage_and_endpoint_RHS_records"], 56
        )
        self.assertEqual(parsed["environment"], binder.ENVIRONMENT)
        self.assertTrue(
            parsed["claims"]["exact_remedy_rejected_on_tested_retry_neighborhood"]
        )
        self.assertFalse(parsed["claims"]["state_advance_authorized"])

    def test_environment_is_exact_and_mismatch_fails_closed(self) -> None:
        with (
            patch.object(
                binder.platform,
                "python_implementation",
                return_value=binder.ENVIRONMENT["python_implementation"],
            ),
            patch.object(
                binder.platform,
                "python_version",
                return_value=binder.ENVIRONMENT["python_version"],
            ) as python_version,
            patch.object(
                binder.np, "__version__", binder.ENVIRONMENT["numpy_version"]
            ),
            patch.object(
                binder.platform,
                "system",
                return_value=binder.ENVIRONMENT["system"],
            ),
            patch.object(
                binder.platform,
                "machine",
                return_value=binder.ENVIRONMENT["machine"],
            ),
            patch.object(binder.sys, "byteorder", binder.ENVIRONMENT["byteorder"]),
        ):
            self.assertEqual(binder._environment(), binder.ENVIRONMENT)
            python_version.return_value = "3.14.4"
            with self.assertRaisesRegex(
                binder.TDG10QA2PREF1Error, "PREF1_ENVIRONMENT_DRIFT"
            ):
                binder._environment()

    def test_live_environment_check_precedes_raw_or_scientific_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(
                    binder, "_static_bound_payload", return_value={"synthetic": True}
                ),
                patch.object(
                    binder,
                    "_environment",
                    side_effect=binder.TDG10QA2PREF1Error(
                        "PREF1_ENVIRONMENT_DRIFT", "synthetic"
                    ),
                ) as environment,
                patch.object(binder, "_raw_snapshot") as raw,
                patch.object(binder, "_restore_shadow") as shadow,
                self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_ENVIRONMENT_DRIFT"
                ),
            ):
                binder.bind_raw_result(self.config_raw, root)
            environment.assert_called_once_with()
            raw.assert_not_called()
            shadow.assert_not_called()

    def test_localization_classes_match_proven_exact_owners(self) -> None:
        self.assertEqual(
            binder._LOCALIZATION_CLASSES,
            frozenset({"unique_maximum", "nonunique_or_interval_inconclusive"}),
        )
        for obsolete in (
            "exact_zero",
            "single_maximizer",
            "coincident_co_maximizers",
        ):
            self.assertNotIn(obsolete, binder._LOCALIZATION_CLASSES)

    def test_forbidden_decision_imports_and_production_calls_are_absent(self) -> None:
        source = (
            ROOT / "src/recursive_horizons/fgc/evolution/tdg10_qa2_pref1_binder.py"
        ).read_text()
        tree = ast.parse(source)
        imports: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                imports.add(module)
                imports.update(f"{module}.{alias.name}" for alias in node.names)
        for fragment in (
            "tdg10_qa2_authority",
            "run_fgc_tdg10_qa2",
            "tdg10_exact_complete_c_runtime",
            "tdg10_exact_complete_c_admission",
            "tdg10_qa1_authority",
            "tdg10_qa1_pref1_binder",
            "tdg9_ti2_authority",
            "run_fgc_tdg9_ti2",
        ):
            self.assertFalse(any(fragment in item for item in imports), fragment)
        for forbidden_call in (
            "commit_tdg6_gr0_compositor(",
            "require_tdg6_temporal_admission(",
            "encode_member(",
            "encode_payload(",
            "serialize_diagnostic_endpoint(",
            "publish_checkpoint(",
            "write_campaign(",
        ):
            self.assertNotIn(forbidden_call, source)

    def test_synthetic_two_leaf_tree_hashes_and_terminal_shape(self) -> None:
        manifest_raw, terminal_raw = _synthetic_raw_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            namespace = root / "synthetic" / "raw"
            namespace.mkdir(parents=True)
            (namespace / "manifest.json").write_bytes(manifest_raw)
            (namespace / "terminal.json").write_bytes(terminal_raw)
            with (
                patch.object(binder, "RAW_NAMESPACE", "synthetic/raw"),
                patch.object(
                    binder,
                    "RAW_MANIFEST_SHA256",
                    sha256(manifest_raw).hexdigest(),
                ),
                patch.object(
                    binder,
                    "RAW_TERMINAL_SHA256",
                    sha256(terminal_raw).hexdigest(),
                ),
            ):
                blobs, manifest, terminal = binder._raw_snapshot(root)
        self.assertEqual(tuple(sorted(blobs)), ("manifest.json", "terminal.json"))
        self.assertEqual(manifest, binder._expected_manifest())
        widths = binder._validate_raw_terminal_shape(terminal)
        self.assertEqual([item["retry"] for item in widths], [4, 5])

    def test_authority_validation_uses_mocked_git_and_synthetic_blobs(self) -> None:
        commit = "c" * 40
        parent = "d" * 40
        blobs = tuple(
            (f"synthetic/path-{index}.txt", f"blob-{index}\n".encode())
            for index in range(16)
        )
        inventory = tuple((path, sha256(raw).hexdigest()) for path, raw in blobs)
        payloads = dict(blobs)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()

            def fake_git(_root: Path, *arguments: str) -> bytes:
                if arguments == ("rev-parse", "--show-toplevel"):
                    return f"{root}\n".encode()
                if arguments == ("cat-file", "-t", commit):
                    return b"commit\n"
                if arguments == ("show", "-s", "--format=%P", commit):
                    return f"{parent}\n".encode()
                if arguments == (
                    "for-each-ref",
                    "--format=%(refname)",
                    "refs/replace",
                ):
                    return b""
                if arguments == (
                    "diff-tree",
                    "--root",
                    "--no-commit-id",
                    "--name-only",
                    "-r",
                    commit,
                ):
                    return ("\n".join(sorted(payloads)) + "\n").encode()
                if arguments[:1] == ("show",) and len(arguments) == 2:
                    prefix = f"{commit}:"
                    return payloads[arguments[1][len(prefix) :]]
                raise AssertionError(arguments)

            with (
                patch.object(binder, "AUTHORITY_COMMIT", commit),
                patch.object(binder, "AUTHORITY_PARENT", parent),
                patch.object(binder, "_AUTHORITY_BLOBS", inventory),
                patch.object(binder, "_git", side_effect=fake_git),
            ):
                binder._authenticate_authority(root, {"authority_commit": commit})

    def test_live_dependency_authority_is_synthetic_and_rejects_drift(self) -> None:
        commit = "e" * 40
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            blobs = []
            for index in range(10):
                relative = f"deps/dependency-{index}.py"
                raw = f"VALUE = {index}\n".encode()
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(raw)
                blobs.append((relative, raw))
            inventory = tuple((path, sha256(raw).hexdigest()) for path, raw in blobs)
            payloads = dict(blobs)

            def fake_git(_root: Path, *arguments: str) -> bytes:
                if arguments == ("rev-parse", "--verify", "HEAD^{commit}"):
                    return f"{commit}\n".encode()
                if arguments[0] in {"diff", "ls-files"}:
                    return b""
                if arguments[:1] == ("show",):
                    prefix = f"{commit}:"
                    return payloads[arguments[1][len(prefix) :]]
                raise AssertionError(arguments)

            with (
                patch.object(binder, "AUTHORITY_COMMIT", commit),
                patch.object(binder, "_LIVE_DEPENDENCY_BLOBS", inventory),
                patch.object(binder, "_git", side_effect=fake_git),
            ):
                binder._authenticate_live_dependencies(root)
                (root / blobs[0][0]).write_bytes(b"VALUE = 'drift'\n")
                with self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_DEPENDENCY_BLOB_DRIFT"
                ):
                    binder._authenticate_live_dependencies(root)

    def test_live_working_boundary_uses_mocked_exact_pre_result_delta(self) -> None:
        commit = "f" * 40
        paths = tuple(sorted(f"planned/path-{index}.txt" for index in range(13)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for index, relative in enumerate(paths):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"planned {index}\n")
            state = {"extra": False, "staged": False}

            def stream(items: tuple[str, ...]) -> bytes:
                return b"" if not items else ("\0".join(items) + "\0").encode()

            def fake_git(_root: Path, *arguments: str) -> bytes:
                if arguments == ("rev-parse", "--verify", "HEAD^{commit}"):
                    return f"{commit}\n".encode()
                if arguments == ("diff", "--name-only", "-z", "--"):
                    return stream(paths[:7])
                if arguments == ("diff", "--cached", "--name-only", "-z", "--"):
                    return stream((paths[0],)) if state["staged"] else b""
                if arguments == (
                    "ls-files",
                    "--others",
                    "--exclude-standard",
                    "-z",
                    "--",
                ):
                    items = (*paths[7:], ".qdrant-initialized")
                    if state["extra"]:
                        items = (*items, "foreign.txt")
                    return stream(items)
                raise AssertionError(arguments)

            with (
                patch.object(binder, "AUTHORITY_COMMIT", commit),
                patch.object(binder, "_PREF1_PRE_RESULT_PATHS", paths),
                patch.object(
                    binder,
                    "_ALLOWED_ADDITIONAL_UNTRACKED",
                    (".qdrant-initialized",),
                ),
                patch.object(binder, "_git", side_effect=fake_git),
            ):
                binder._authenticate_live_working_boundary(root)
                state["extra"] = True
                with self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_WORKING_BOUNDARY_DRIFT"
                ):
                    binder._authenticate_live_working_boundary(root)
                state["extra"] = False
                state["staged"] = True
                with self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_WORKING_BOUNDARY_DRIFT"
                ):
                    binder._authenticate_live_working_boundary(root)

    def test_store_snapshot_algorithm_uses_only_synthetic_temp_store(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = root / "synthetic-store"
            (store / "nested").mkdir(parents=True)
            (store / "a.txt").write_bytes(b"alpha\n")
            (store / "nested" / "b.txt").write_bytes(b"beta\n")
            with patch.object(binder, "STORE_PATH", "synthetic-store"):
                before = binder._snapshot_store(root)
                after = binder._snapshot_store(root)
        self.assertEqual(before[0], 2)
        self.assertEqual(after, before)

    def test_terminal_reduction_covers_all_four_publishable_branches(self) -> None:
        def records(first: str, second: str) -> list[dict[str, object]]:
            return [
                {"retry": 4, "width_class": first},
                {"retry": 5, "width_class": second},
            ]

        self.assertEqual(
            binder.reduce_qa2_terminal(
                records("all_18_channel_pass", "all_18_channel_pass")
            ),
            binder.BOTH_PASS_CLASS,
        )
        self.assertEqual(
            binder.reduce_qa2_terminal(
                records("all_18_channel_pass", "one_or_more_channel_nonpass")
            ),
            binder.RAW_CLASSIFICATION,
        )
        self.assertEqual(
            binder.reduce_qa2_terminal(records("inconclusive", "all_18_channel_pass")),
            binder.INCONCLUSIVE_CLASS,
        )
        self.assertEqual(
            binder.reduce_qa2_terminal(records("premise_stop", "inconclusive")),
            binder.PREMISE_STOP_CLASS,
        )

    def test_synthetic_terminal_rejects_identity_schema_and_promotion_mutations(
        self,
    ) -> None:
        mutations = []
        changed = _synthetic_terminal()
        changed["widths"][0]["predecessor_generation"] = 11
        mutations.append(changed)
        changed = _synthetic_terminal()
        changed["widths"][1]["attempted_width_hex"] = "0x1.0p-13"
        mutations.append(changed)
        changed = _synthetic_terminal()
        changed["widths"][0].pop("exact_complete_C")
        mutations.append(changed)
        changed = _synthetic_terminal()
        changed["state_advance_authorized"] = True
        mutations.append(changed)
        for changed in mutations:
            with self.subTest(keys=sorted(changed)):
                with self.assertRaises(binder.TDG10QA2PREF1Error):
                    binder._validate_raw_terminal_shape(changed)

    def test_top_typed_premise_stop_must_equal_owning_width(self) -> None:
        widths = [_synthetic_width(dict(width)) for width in binder.WIDTHS]
        stop = {"type": "TDG6RefinementPathStop", "detail": "bounded"}
        widths[0] = {
            "retry": 4,
            "width_class": "premise_stop",
            "predecessor_generation": 10,
            "attempted_width_hex": binder.WIDTHS[0]["attempted_width_hex"],
            "replay_receipt": binder._expected_replay_receipt(binder.WIDTHS[0]),
            "typed_stop": stop,
        }
        terminal = {
            **binder._expected_terminal_base(
                binder.PREMISE_STOP_CLASS, robustness=False
            ),
            "typed_stop": {**stop, "retry": 4},
            "widths": widths,
        }
        binder._validate_raw_terminal_shape(terminal)
        terminal["typed_stop"]["detail"] = "different"
        with self.assertRaisesRegex(binder.TDG10QA2PREF1Error, "PREF1_STOP_DRIFT"):
            binder._validate_raw_terminal_shape(terminal)
        terminal["typed_stop"] = {**stop, "retry": 4, "extra": False}
        with self.assertRaises(binder.TDG10QA2PREF1Error):
            binder._validate_raw_terminal_shape(terminal)

    def test_top_typed_inconclusive_must_equal_owning_width(self) -> None:
        widths = [_synthetic_width(dict(width)) for width in binder.WIDTHS]
        closed = {
            "type": "ExactCompleteCResourceExhausted",
            "reason": "candidate_ceiling_exhausted",
            "detail": "bounded",
        }
        widths[0] = {
            "retry": 4,
            "width_class": "inconclusive",
            "predecessor_generation": 10,
            "attempted_width_hex": binder.WIDTHS[0]["attempted_width_hex"],
            "replay_receipt": binder._expected_replay_receipt(binder.WIDTHS[0]),
            "shadow_path_count": 7,
            "shadow_proposal_count": 7,
            "SSPRK3_stage_and_endpoint_record_count": 28,
            "typed_inconclusive": closed,
        }
        terminal = {
            **binder._expected_terminal_base(
                binder.INCONCLUSIVE_CLASS, robustness=False
            ),
            "typed_inconclusive": {**closed, "retry": 4},
            "widths": widths,
        }
        binder._validate_raw_terminal_shape(terminal)
        for field, replacement in (
            ("type", "UnknownClosure"),
            ("reason", "different"),
            ("detail", "different"),
        ):
            changed = deepcopy(terminal)
            changed["typed_inconclusive"][field] = replacement
            with (
                self.subTest(field=field),
                self.assertRaises(binder.TDG10QA2PREF1Error),
            ):
                binder._validate_raw_terminal_shape(changed)

    def test_compact_rich_channel_recomputes_all_semantics(self) -> None:
        valid = _rich_channel()
        binder._validate_rich_channel(valid, expected_channel="u:R")
        for field in (
            "sufficient_contraction_failure",
            "threshold_inconclusive",
        ):
            changed = deepcopy(valid)
            changed[field] = not changed[field]
            with (
                self.subTest(field=field),
                self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_RECLASSIFY_DRIFT"
                ),
            ):
                binder._validate_rich_channel(changed, expected_channel="u:R")
        changed = deepcopy(valid)
        changed["combined_coefficient_stream_sha256"] = "0" * 64
        with self.assertRaisesRegex(
            binder.TDG10QA2PREF1Error, "PREF1_RECLASSIFY_DRIFT"
        ):
            binder._validate_rich_channel(changed, expected_channel="u:R")

    def test_real_compact_result_has_exact_bound_byte_hash(self) -> None:
        self.assertEqual(
            binder.COMPACT_RESULT_SHA256,
            "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6",
        )
        self.assertEqual(
            sha256(self.result_raw).hexdigest(), binder.COMPACT_RESULT_SHA256
        )
        self.assertEqual(len(self.result_raw), 253881)
        binder.validate_compact_result(self.config_raw, self.result_raw)
        with self.assertRaisesRegex(
            binder.TDG10QA2PREF1Error, "PREF1_COMPACT_HASH_DRIFT"
        ):
            binder.validate_compact_result(self.config_raw, self.result_raw + b" ")

    def test_real_compact_rejects_rational_channel_and_all_of_mutations(self) -> None:
        def rational(payload: dict[str, object]) -> None:
            upper = payload["independent_replay"]["widths"][0]["channels"][0]["D12"][
                "upper"
            ]
            original = int(upper["numerator"])
            changed = original ^ 0b10
            self.assertNotEqual(changed, original)
            upper["numerator"] = str(changed)

        def channel_class(payload: dict[str, object]) -> None:
            payload["independent_replay"]["widths"][0]["channels"][0][
                "classification"
            ] = "resolved_order_pass"

        def channel_pass(payload: dict[str, object]) -> None:
            payload["independent_replay"]["widths"][0]["channels"][0][
                "admission_passed"
            ] = True

        def failed_channels(payload: dict[str, object]) -> None:
            payload["independent_replay"]["widths"][0]["failed_channels"] = [
                "u:lambda",
                "u:R",
            ]

        def all_of(payload: dict[str, object]) -> None:
            payload["independent_replay"]["all_of_identity_verified"] = False

        def robustness(payload: dict[str, object]) -> None:
            payload["independent_replay"]["width_robustness_passed"] = True

        for name, mutation in (
            ("rational", rational),
            ("channel_class", channel_class),
            ("channel_pass", channel_pass),
            ("failed_channels", failed_channels),
            ("all_of", all_of),
            ("robustness", robustness),
        ):
            with self.subTest(name=name):
                self.assert_semantic_mutation_rejected(mutation)

    def test_real_compact_rejects_deep_exact_route_mutations(self) -> None:
        def sufficient_failure(payload: dict[str, object]) -> None:
            channel = payload["independent_replay"]["widths"][0]["channels"][0]
            channel["sufficient_contraction_failure"] = not channel[
                "sufficient_contraction_failure"
            ]

        def threshold_inconclusive(payload: dict[str, object]) -> None:
            channel = payload["independent_replay"]["widths"][0]["channels"][0]
            channel["threshold_inconclusive"] = not channel["threshold_inconclusive"]

        def stationary(payload: dict[str, object]) -> None:
            difference = payload["independent_replay"]["widths"][0]["channels"][0][
                "D01"
            ]
            difference["independent_stationary_count_stream_sha256"] = "0" * 64

        def localization(payload: dict[str, object]) -> None:
            difference = payload["independent_replay"]["widths"][0]["channels"][0][
                "D01"
            ]
            difference["localization_classification"] = (
                "nonunique_or_interval_inconclusive"
            )

        def co_maximizer(payload: dict[str, object]) -> None:
            difference = payload["independent_replay"]["widths"][0]["channels"][0][
                "D01"
            ]
            difference["co_maximizer_count"] = 2

        def combined(payload: dict[str, object]) -> None:
            payload["independent_replay"]["widths"][0]["channels"][0][
                "combined_coefficient_stream_sha256"
            ] = "0" * 64

        for name, mutation in (
            ("sufficient_failure", sufficient_failure),
            ("threshold_inconclusive", threshold_inconclusive),
            ("stationary", stationary),
            ("localization", localization),
            ("co_maximizer", co_maximizer),
            ("combined", combined),
        ):
            with self.subTest(name=name):
                self.assert_semantic_mutation_rejected(mutation)

    def test_real_compact_rejects_counts_conclusion_nonclaim_and_hash_mutations(
        self,
    ) -> None:
        def replay_count(payload: dict[str, object]) -> None:
            payload["independent_replay"]["shadow_path_count"] = 13

        def width_count(payload: dict[str, object]) -> None:
            payload["independent_replay"]["widths"][0][
                "SSPRK3_stage_and_endpoint_record_count"
            ] = 27

        def conclusion(payload: dict[str, object]) -> None:
            payload["conclusion"]["production_method_earned"] = True

        def scope_nonclaim(payload: dict[str, object]) -> None:
            payload["scope"]["production_commit_called"] = True

        def claim(payload: dict[str, object]) -> None:
            payload["claims"]["state_advance_authorized"] = True

        def authority_hash(payload: dict[str, object]) -> None:
            payload["authority"]["blobs"][0]["sha256"] = "0" * 64

        def raw_hash(payload: dict[str, object]) -> None:
            payload["raw"]["terminal_sha256"] = "0" * 64

        for name, mutation in (
            ("replay_count", replay_count),
            ("width_count", width_count),
            ("conclusion", conclusion),
            ("scope_nonclaim", scope_nonclaim),
            ("claim", claim),
            ("authority_hash", authority_hash),
            ("raw_hash", raw_hash),
        ):
            with self.subTest(name=name):
                self.assert_semantic_mutation_rejected(mutation)

    def test_compact_difference_rejects_route_and_class_count_mutations(self) -> None:
        valid = _difference(
            level="D01",
            lower=Fraction(1),
            upper=Fraction(1),
            coefficient="1" * 64,
        )
        binder._validate_difference(valid, level="u:R.D01")
        changed = deepcopy(valid)
        changed["independent_stationary_count_stream_sha256"] = "9" * 64
        with self.assertRaisesRegex(
            binder.TDG10QA2PREF1Error, "PREF1_DIFFERENCE_CONTRACT_DRIFT"
        ):
            binder._validate_difference(changed, level="u:R.D01")
        changed = deepcopy(valid)
        changed["candidate_count"] = 2
        changed["co_maximizer_count"] = 2
        with self.assertRaisesRegex(
            binder.TDG10QA2PREF1Error, "PREF1_DIFFERENCE_CONTRACT_DRIFT"
        ):
            binder._validate_difference(changed, level="u:R.D01")
        changed = deepcopy(valid)
        changed["localization_classification"] = "nonunique_or_interval_inconclusive"
        with self.assertRaisesRegex(
            binder.TDG10QA2PREF1Error, "PREF1_DIFFERENCE_CONTRACT_DRIFT"
        ):
            binder._validate_difference(changed, level="u:R.D01")

    def test_two_leaf_reader_rejects_partial_and_symlinked_synthetic_tree(self) -> None:
        manifest_raw, _ = _synthetic_raw_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            namespace = root / "synthetic-raw"
            namespace.mkdir()
            (namespace / "manifest.json").write_bytes(manifest_raw)
            with patch.object(binder, "RAW_NAMESPACE", "synthetic-raw"):
                with self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_RAW_TREE_DRIFT"
                ):
                    binder._raw_tree(root)
        with (
            tempfile.TemporaryDirectory() as temporary,
            tempfile.TemporaryDirectory() as foreign,
        ):
            root = Path(temporary)
            (root / "synthetic-raw").symlink_to(foreign)
            with patch.object(binder, "RAW_NAMESPACE", "synthetic-raw"):
                with self.assertRaisesRegex(
                    binder.TDG10QA2PREF1Error, "PREF1_PATH_UNSAFE"
                ):
                    binder._raw_tree(root)

    def test_compact_parser_rejects_duplicate_nonfinite_and_partial_results(
        self,
    ) -> None:
        for raw in (
            b'{"artifact_id":"FGC-1-TDG10-QA2-PREF1"}\n',
            b'{"artifact_id":"x","artifact_id":"y","artifact_payload":{}}\n',
            b'{"artifact_id":"x","artifact_payload":{"value":NaN}}\n',
        ):
            with self.assertRaises(binder.TDG10QA2PREF1Error):
                binder.validate_compact_result(self.config_raw, raw)

    def test_reproducer_default_missing_result_is_clear_and_never_live(self) -> None:
        missing = "results/__qa2_pref1_intentionally_absent__.json"
        completed = subprocess.run(
            (
                sys.executable,
                "-B",
                "scripts/reproduce_fgc_tdg10_qa2_pref1.py",
                "--output",
                missing,
            ),
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("compact result is absent", completed.stderr)
        self.assertFalse((ROOT / missing).exists())

    def test_reproducer_reads_and_publication_reject_symlinks_and_clobber(self) -> None:
        module = _load_reproducer()
        with (
            tempfile.TemporaryDirectory() as temporary,
            tempfile.TemporaryDirectory() as foreign,
        ):
            root = Path(temporary)
            (root / "safe").mkdir()
            (root / "safe" / "leaf.json").write_bytes(b"first\n")
            self.assertEqual(
                module.read_repository_leaf(root, "safe/leaf.json"), b"first\n"
            )
            (root / "safe" / "link.json").symlink_to(root / "safe" / "leaf.json")
            with self.assertRaisesRegex(SystemExit, "unsafe PREF1 leaf"):
                module.read_repository_leaf(root, "safe/link.json")
            (root / "linked-parent").symlink_to(foreign)
            with self.assertRaisesRegex(SystemExit, "unsafe PREF1 parent"):
                module.read_repository_leaf(root, "linked-parent/leaf.json")
            module.write_canonical_result(root, "safe/result.json", b"first\n")
            with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                module.write_canonical_result(root, "safe/result.json", b"second\n")
            (root / "safe" / "symlink-result.json").symlink_to(
                root / "safe" / "leaf.json"
            )
            with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                module.write_canonical_result(
                    root, "safe/symlink-result.json", b"second\n"
                )
            self.assertEqual((root / "safe" / "result.json").read_bytes(), b"first\n")
            self.assertFalse(
                any(".tmp-" in path.name for path in (root / "safe").iterdir())
            )

    def test_ordinary_suite_source_has_no_live_evidence_access(self) -> None:
        source = Path(__file__).read_text()
        forbidden = (
            "binder." + "RAW_NAMESPACE",
            "_snapshot_store(" + "ROOT)",
            "ROOT / binder." + "STORE_PATH",
            "subprocess.run((" + '"git"',
            "git " + "show",
            "git " + "diff-tree",
        )
        for text in forbidden:
            self.assertNotIn(text, source)
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            function = node.func
            if isinstance(function, ast.Attribute) and function.attr in {
                "_raw_snapshot",
                "_raw_tree",
                "_snapshot_store",
            }:
                self.assertTrue(node.args)
                first = node.args[0]
                self.assertFalse(
                    isinstance(first, ast.Name) and first.id == "ROOT",
                    ast.unparse(node),
                )


if __name__ == "__main__":
    unittest.main()
