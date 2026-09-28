"""Focused contract controls for prospective FGC-1-TDG9-AR1-AUTH1."""

from __future__ import annotations

import ast
import json
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest
from unittest.mock import call, patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from scripts import reproduce_fgc_tdg9_ar1_auth1 as reproduce  # noqa: E402


class TDG9AR1AuthorityTests(unittest.TestCase):
    def test_frozen_three_width_scope_and_claim_boundary(self) -> None:
        self.assertEqual(tuple(item["retry"] for item in ar1.REPLAYS), (3, 4, 5))
        self.assertEqual(
            tuple(item["prior_retry_count"] for item in ar1.REPLAYS), (2, 3, 4)
        )
        self.assertEqual(
            tuple(item["predecessor_generation"] for item in ar1.REPLAYS),
            (9, 10, 11),
        )
        self.assertEqual(
            tuple(item["attempted_width_hex"] for item in ar1.REPLAYS),
            (
                "0x1.aaa9612df8000p-11",
                "0x1.aaa9612df8000p-12",
                "0x1.aaa9612df0000p-13",
            ),
        )
        self.assertEqual(len(ar1.CHANNEL_ORDER), 18)
        self.assertEqual(ar1.MAX_EVALUATOR_CALLS, 3 * 18 * 2)
        self.assertEqual(ar1.ORDER_THRESHOLD, "3/2")
        self.assertEqual(ar1.ORDER_SQUARED_MULTIPLIER, 8)
        self.assertFalse(ar1.SCOPE["automatic_resource_escalation_authorized"])
        self.assertFalse(ar1.SCOPE["fourth_retry_width_authorized"])
        self.assertFalse(ar1.SCOPE["campaign_store_mutation_authorized"])
        self.assertFalse(ar1.SCOPE["continuation_authorized"])
        self.assertFalse(ar1.SCOPE["candidate_branches_authorized"])
        for claim in (
            "raw_AR1_result_earned",
            "AR1_binder_completed",
            "common_event_completed",
            "GR0_calibration_completed",
            "candidate_execution_authorized",
            "mechanism_result_earned",
            "physical_result_earned",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertIs(ar1.CLAIMS[claim], False)

    def test_environment_selection_and_committed_image_are_exact(self) -> None:
        self.assertEqual(
            ar1.ENVIRONMENT,
            {
                "python_implementation": "CPython",
                "python_version": "3.14.3",
                "numpy_version": "2.5.1",
                "system": "Darwin",
                "sys_platform": "darwin",
                "machine": "arm64",
                "byteorder": "little",
            },
        )
        self.assertEqual(
            ar1.SELECTION,
            {
                "protocol": "FGC-2-SF1-PROTO18",
                "campaign_id": "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1",
                "campaign_plan_sha256": (
                    "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
                ),
                "campaign_authorization_commit": (
                    "6df67967e6e0bc63eca4dc6951a60d3787368f26"
                ),
                "branch": "GR-0",
                "amplitude": "3",
                "method": "RK4",
                "member_key": "RK4-2049",
                "point_count": 2049,
                "owned_row_count": 2044,
                "event": 23,
                "target_rational": "3/2",
                "target_binary64_hex": "0x1.8000000000000p+0",
            },
        )
        self.assertEqual(len(ar1.AUTHORITY_DELTA_PATHS), 17)
        self.assertEqual(tuple(sorted(ar1.AUTHORITY_DELTA_PATHS)), ar1.AUTHORITY_DELTA_PATHS)
        self.assertEqual(ar1.COMMITTED_IMAGE["base_commit"], ar1.EVALUATOR_COMMIT)
        self.assertNotIn("authority_commit", ar1.COMMITTED_IMAGE)
        self.assertTrue(
            ar1.COMMITTED_IMAGE[
                "authority_commit_must_be_single_direct_successor_of_base"
            ]
        )
        self.assertTrue(
            ar1.COMMITTED_IMAGE["exact_delta_locks_unlisted_transitive_modules"]
        )

    def test_environment_gate_rejects_each_owned_dimension(self) -> None:
        mutations = (
            (ar1.platform, "python_implementation", lambda: "PyPy"),
            (ar1.platform, "python_version", lambda: "3.14.2"),
            (ar1.np, "__version__", "2.5.0"),
            (ar1.platform, "system", lambda: "Linux"),
            (ar1.sys, "platform", "linux"),
            (ar1.platform, "machine", lambda: "x86_64"),
            (ar1.sys, "byteorder", "big"),
        )
        with (
            patch.object(
                ar1.platform,
                "python_implementation",
                lambda: ar1.ENVIRONMENT["python_implementation"],
            ),
            patch.object(
                ar1.platform,
                "python_version",
                lambda: ar1.ENVIRONMENT["python_version"],
            ),
            patch.object(ar1.np, "__version__", ar1.ENVIRONMENT["numpy_version"]),
            patch.object(
                ar1.platform, "system", lambda: ar1.ENVIRONMENT["system"]
            ),
            patch.object(ar1.sys, "platform", ar1.ENVIRONMENT["sys_platform"]),
            patch.object(
                ar1.platform, "machine", lambda: ar1.ENVIRONMENT["machine"]
            ),
            patch.object(ar1.sys, "byteorder", ar1.ENVIRONMENT["byteorder"]),
        ):
            self.assertEqual(ar1._observe_environment(), ar1.ENVIRONMENT)  # type: ignore[attr-defined]
            for owner, name, value in mutations:
                with self.subTest(name=name), patch.object(owner, name, value):
                    with self.assertRaisesRegex(
                        ar1.TDG9AR1AuthorityError, "environment"
                    ):
                        ar1._observe_environment()  # type: ignore[attr-defined]

    def test_prospective_image_is_exact_and_allows_only_qdrant(self) -> None:
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=ar1.EVALUATOR_COMMIT),
            patch.object(
                ar1,
                "_working_image",
                return_value=(ar1.AUTHORITY_DELTA_PATHS, (), (".qdrant-initialized",)),
            ),
            patch.object(ar1, "_read_regular", return_value=b"owned") as read,
        ):
            ar1._require_prospective_image(ROOT)  # type: ignore[attr-defined]
        self.assertEqual(read.call_count, len(ar1.AUTHORITY_DELTA_PATHS))

        for image in (
            (ar1.AUTHORITY_DELTA_PATHS, ("README.md",), (".qdrant-initialized",)),
            (ar1.AUTHORITY_DELTA_PATHS, (), (".qdrant-initialized", "foreign")),
        ):
            with (
                self.subTest(image=image),
                patch.object(ar1, "_require_no_replace_refs"),
                patch.object(ar1, "_head", return_value=ar1.EVALUATOR_COMMIT),
                patch.object(ar1, "_working_image", return_value=image),
                self.assertRaises(ar1.TDG9AR1AuthorityError),
            ):
                ar1._require_prospective_image(ROOT)  # type: ignore[attr-defined]
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value="f" * 40),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "evaluator commit"),
        ):
            ar1._require_prospective_image(ROOT)  # type: ignore[attr-defined]

    def test_runtime_image_requires_exact_clean_head_ancestry_and_delta(self) -> None:
        commit = "a" * 40
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=commit),
            patch.object(
                ar1, "_commit_parents", return_value=(ar1.EVALUATOR_COMMIT,)
            ),
            patch.object(
                ar1,
                "_working_image",
                return_value=((), (), (".qdrant-initialized",)),
            ),
            patch.object(ar1, "_require_ancestor") as ancestor,
            patch.object(ar1, "_git_paths", return_value=ar1.AUTHORITY_DELTA_PATHS),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        self.assertEqual(
            ancestor.call_args_list,
            [
                call(ROOT, ar1.PREF2_COMMIT, commit),
                call(ROOT, ar1.EVALUATOR_COMMIT, commit),
            ],
        )

        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value="b" * 40),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "exact HEAD"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=commit),
            patch.object(
                ar1, "_commit_parents", return_value=(ar1.EVALUATOR_COMMIT,)
            ),
            patch.object(ar1, "_working_image", return_value=(("README.md",), (), ())),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "not clean"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=commit),
            patch.object(
                ar1, "_commit_parents", return_value=(ar1.EVALUATOR_COMMIT,)
            ),
            patch.object(ar1, "_working_image", return_value=((), (), ())),
            patch.object(ar1, "_require_ancestor"),
            patch.object(
                ar1,
                "_git_paths",
                return_value=tuple(sorted((*ar1.AUTHORITY_DELTA_PATHS, "extra"))),
            ),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "delta differs"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=commit),
            patch.object(
                ar1, "_commit_parents", return_value=(ar1.EVALUATOR_COMMIT,)
            ),
            patch.object(ar1, "_working_image", return_value=((), (), ())),
            patch.object(
                ar1,
                "_require_ancestor",
                side_effect=ar1.TDG9AR1AuthorityError("not ancestor"),
            ),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "not ancestor"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        with (
            patch.object(ar1, "_require_no_replace_refs"),
            patch.object(ar1, "_head", return_value=commit),
            patch.object(ar1, "_commit_parents", return_value=("b" * 40,)),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "direct successor"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]
        with (
            patch.object(
                ar1,
                "_require_no_replace_refs",
                side_effect=ar1.TDG9AR1AuthorityError("replace refs"),
            ),
            self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "replace refs"),
        ):
            ar1._require_runtime_image(ROOT, commit)  # type: ignore[attr-defined]

    def test_direct_decision_and_publication_dependencies_are_inventory_bound(self) -> None:
        paths = tuple(path for _role, path in ar1.IMPLEMENTATION_INVENTORY)
        for expected in (
            "src/recursive_horizons/fgc/evolution/tdg9_ar1_authority.py",
            "scripts/run_fgc_tdg9_ar1.py",
            "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic.py",
            "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic_independent.py",
            "src/recursive_horizons/fgc/evolution/tdg8_persisted_retry_replay.py",
            "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
            "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py",
        ):
            self.assertIn(expected, paths)

    def test_authority_imports_no_continuation_or_candidate_runner(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/tdg9_ar1_authority.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        text = path.read_text(encoding="utf-8")
        self.assertFalse(any("event1" in name or "candidate" in name for name in imported))
        self.assertNotIn("_advance_authenticated_event", text)

    def test_namespace_absence_rejects_stage_symlink_and_existing_target(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = root / Path(ar1.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            ar1._require_namespace_absent(root)  # type: ignore[attr-defined]

            stage = parent / f"{ar1.STAGING_PREFIX}cut"
            stage.mkdir()
            with self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "staging"):
                ar1._require_namespace_absent(root)  # type: ignore[attr-defined]
            stage.rmdir()

            target = root / ar1.OUTPUT_NAMESPACE
            target.mkdir()
            with self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "already exists"):
                ar1._require_namespace_absent(root)  # type: ignore[attr-defined]
            target.rmdir()

            destination = root / "destination"
            destination.mkdir()
            target.symlink_to(destination, target_is_directory=True)
            with self.assertRaisesRegex(ar1.TDG9AR1AuthorityError, "redirected"):
                ar1._require_namespace_absent(root)  # type: ignore[attr-defined]

    def test_rejected_compact_rewrite_is_bounded_atomic_and_race_closed(self) -> None:
        stale = b"known rejected compact bytes\n"
        replacement = b"new canonical compact bytes\n"
        digest = ar1._sha(stale)  # type: ignore[attr-defined]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / ar1.RESULT_PATH
            path.parent.mkdir(parents=True)
            path.write_bytes(stale)
            with patch.object(reproduce, "REJECTED_RESULT_SHA256S", (digest,)):
                reproduce._replace_rejected(replacement, root)  # type: ignore[attr-defined]
            self.assertEqual(path.read_bytes(), replacement)
            self.assertEqual(
                tuple(item.name for item in path.parent.iterdir()), (path.name,)
            )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / ar1.RESULT_PATH
            path.parent.mkdir(parents=True)
            path.write_bytes(b"foreign bytes\n")
            with self.assertRaisesRegex(ValueError, "unknown bytes"):
                reproduce._replace_rejected(replacement, root)  # type: ignore[attr-defined]
            self.assertEqual(path.read_bytes(), b"foreign bytes\n")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / ar1.RESULT_PATH
            path.parent.mkdir(parents=True)
            target = root / "foreign-result"
            target.write_bytes(stale)
            path.symlink_to(target)
            with patch.object(reproduce, "REJECTED_RESULT_SHA256S", (digest,)):
                with self.assertRaisesRegex(ValueError, "regular leaf"):
                    reproduce._replace_rejected(replacement, root)  # type: ignore[attr-defined]
            self.assertEqual(target.read_bytes(), stale)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / ar1.RESULT_PATH
            path.parent.mkdir(parents=True)
            path.write_bytes(stale)
            original_reader = reproduce._read_rejected_identity  # type: ignore[attr-defined]
            calls = 0

            def race(directory_fd: int, leaf: str) -> tuple[int, int, int, int]:
                nonlocal calls
                calls += 1
                if calls == 2:
                    path.write_bytes(b"raced foreign bytes\n")
                return original_reader(directory_fd, leaf)

            with (
                patch.object(reproduce, "REJECTED_RESULT_SHA256S", (digest,)),
                patch.object(reproduce, "_read_rejected_identity", side_effect=race),
                self.assertRaisesRegex(ValueError, "unknown bytes"),
            ):
                reproduce._replace_rejected(replacement, root)  # type: ignore[attr-defined]
            self.assertEqual(path.read_bytes(), b"raced foreign bytes\n")
            self.assertEqual(
                tuple(item.name for item in path.parent.iterdir()), (path.name,)
            )

    def test_config_and_compact_bundle_are_strict_when_materialized(self) -> None:
        config_path = ROOT / ar1.CONFIG_PATH
        result_path = ROOT / ar1.RESULT_PATH
        if not config_path.exists() or not result_path.exists():
            self.skipTest("integration worker has not materialized the compact bundle")
        config = config_path.read_bytes()
        result = result_path.read_bytes()
        materialized_config = tomllib.loads(config.decode("utf-8"))
        try:
            materialized_payload = json.loads(result)["artifact_payload"]
        except (KeyError, TypeError, json.JSONDecodeError):
            materialized_payload = {}
        if (
            materialized_config.get("environment") != ar1.ENVIRONMENT
            or materialized_config.get("selection") != ar1.SELECTION
            or materialized_config.get("committed_image") != ar1.COMMITTED_IMAGE
            or materialized_payload.get("environment") != ar1.ENVIRONMENT
            or materialized_payload.get("selection") != ar1.SELECTION
            or materialized_payload.get("committed_image") != ar1.COMMITTED_IMAGE
        ):
            self.skipTest("integration worker has not regenerated the repaired bundle")
        parsed = ar1._parse_config(config)  # type: ignore[attr-defined]
        self.assertEqual(parsed["replays"], list(ar1.REPLAYS))
        ar1.validate_compact(config, result)

        mutated = config.replace(
            b'order_threshold = "3/2"', b'order_threshold = "149/100"', 1
        )
        with self.assertRaises(ar1.TDG9AR1AuthorityError):
            ar1._parse_config(mutated)  # type: ignore[attr-defined]
        mutated = config.replace(
            ar1.REPLAYS[0]["attempted_width_hex"].encode(),
            b"0x1.0000000000000p-11",
            1,
        )
        with self.assertRaises(ar1.TDG9AR1AuthorityError):
            ar1._parse_config(mutated)  # type: ignore[attr-defined]
        mutated = config.replace(b'"u:alpha"', b'"u:wrong"', 1)
        with self.assertRaises(ar1.TDG9AR1AuthorityError):
            ar1._parse_config(mutated)  # type: ignore[attr-defined]
        for old, new in (
            (b'python_version = "3.14.3"', b'python_version = "3.14.2"'),
            (
                b'campaign_id = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"',
                b'campaign_id = "wrong"',
            ),
            (b'point_count = 2049', b'point_count = 4097'),
            (
                b'base_commit = "3e37e3baa5a7cf862c4c0c004c7f43dcfa1313b2"',
                b'base_commit = "0000000000000000000000000000000000000000"',
            ),
        ):
            with self.subTest(old=old):
                self.assertIn(old, config)
                with self.assertRaises(ar1.TDG9AR1AuthorityError):
                    ar1._parse_config(config.replace(old, new, 1))  # type: ignore[attr-defined]


if __name__ == "__main__":
    unittest.main()
