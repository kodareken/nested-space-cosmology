from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for entry in (ROOT, SRC):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg10_qa1_pref1_binder as binder,
)
from scripts import reproduce_fgc_tdg10_qa1_pref1 as reproducer  # noqa: E402


CONFIG_RAW = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT = ROOT / binder.RESULT_PATH


def _rational(value: Fraction) -> dict[str, str]:
    return {
        "denominator": str(value.denominator),
        "numerator": str(value.numerator),
    }


def _difference(level: str, lower: Fraction, upper: Fraction) -> dict[str, object]:
    is_d01 = level == "D01"
    return {
        "candidate_count": 1,
        "co_maximizer_count": 1,
        "coefficient_stream_sha256": "1" * 64,
        "independent_evaluator_id": binder.INDEPENDENT_EVALUATOR_ID,
        "independent_stationary_count_stream_sha256": "2" * 64,
        "localization_classification": "unique_maximum",
        "lower": _rational(lower),
        "maximum_candidates": (
            binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D01
            if is_d01
            else binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        ),
        "polynomial_count": (
            binder.D01_POLYNOMIAL_COUNT if is_d01 else binder.D12_POLYNOMIAL_COUNT
        ),
        "primary_evaluator_id": binder.PRIMARY_EVALUATOR_ID,
        "primary_stationary_count_stream_sha256": "3" * 64,
        "routes_agree": True,
        "survivor_key_stream_sha256": "4" * 64,
        "upper": _rational(upper),
    }


def _channel(
    expected: tuple[str, str, str],
) -> dict[str, object]:
    channel, numerator, denominator = expected
    d12 = Fraction(int(numerator), int(denominator))
    d01 = 4 * d12
    left = binder.TDG6_ORDER_SQUARED_MULTIPLIER * d12**2
    right = d01**2
    return {
        "D01": _difference("D01", d01, d01),
        "D12": _difference("D12", d12, d12),
        "admission_passed": True,
        "candidate_evidence_available": True,
        "channel": channel,
        "classification": "resolved_order_pass",
        "combined_coefficient_stream_sha256": "5" * 64,
        "evaluator_id": binder.TDG10_EVALUATOR_ID,
        "hash_evidence_available": True,
        "maximum_candidates_D01": binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "maximum_candidates_D12": binder.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "order_threshold_passed": True,
        "order_threshold_resolved": True,
        "refinement_depth": binder.PRIMARY_REFINEMENT_DEPTH,
        "route_evidence_available": True,
        "row_count": binder.OWNED_ROW_COUNT,
        "row_stream_sha256": "6" * 64,
        "sufficient_condition": binder.SUFFICIENT_CONDITION,
        "sufficient_condition_holds": True,
        "sufficient_contraction_failure": False,
        "sufficient_contraction_pass": True,
        "sufficient_pass_left": _rational(left),
        "sufficient_pass_right": _rational(right),
        "temporal_retry_permitted": False,
        "threshold_inconclusive": False,
    }


def compact_fixture() -> dict[str, object]:
    """Build a compact-only fixture without opening raw/store/shadow inputs."""

    result = binder.expected_compact_result(CONFIG_RAW)
    result["artifact_payload"]["independent_replay"]["channels"] = [
        _channel(item) for item in binder.EXPECTED_D12_UPPERS
    ]
    return result


def _write_raw_fixture(root: Path) -> dict[str, bytes]:
    raw = root / binder.RAW_NAMESPACE
    (raw / "fine-endpoint").mkdir(parents=True)
    (raw / "payloads").mkdir()
    blobs = {
        "manifest.json": b'{"fixture":"manifest"}\n',
        "terminal.json": b'{"fixture":"terminal"}\n',
        binder.DESCRIPTOR_RELATIVE: b'{"fixture":"descriptor"}\n',
        binder.PAYLOAD_RELATIVE: b"fixture-payload\n",
    }
    for relative, value in blobs.items():
        path = raw / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)
    return blobs


def _raw_hash_patches(blobs: dict[str, bytes]):
    return patch.multiple(
        binder,
        RAW_MANIFEST_SHA256=sha256(blobs["manifest.json"]).hexdigest(),
        RAW_TERMINAL_SHA256=sha256(blobs["terminal.json"]).hexdigest(),
        DESCRIPTOR_RAW_SHA256=sha256(blobs[binder.DESCRIPTOR_RELATIVE]).hexdigest(),
        PAYLOAD_RAW_SHA256=sha256(blobs[binder.PAYLOAD_RELATIVE]).hexdigest(),
    )


class TDG10QA1PREF1CompactTests(unittest.TestCase):
    def test_checked_in_config_is_the_exact_typed_contract(self) -> None:
        config = tomllib.loads(CONFIG_RAW.decode("utf-8"))
        self.assertEqual(config, binder.expected_config())
        self.assertEqual(config["authority"]["blob_count"], 20)
        self.assertEqual(len(config["authority"]["blobs"]), 20)
        self.assertEqual(len(binder._AUTHORITY_BLOBS), 20)
        self.assertEqual(config["decision"]["sufficient_condition"], "8*U12^2 <= L01^2")
        self.assertTrue(config["decision"]["equality_passes"])

    def test_compact_validation_is_raw_store_shadow_and_git_blind(self) -> None:
        encoded = binder.canonical_result(compact_fixture())
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "COMPACT_RESULT_SHA256", sha256(encoded).hexdigest()),
            patch.object(binder, "_raw_tree", side_effect=forbidden),
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "_restore_shadow", side_effect=forbidden),
            patch.object(
                binder, "_assess_prepared_independently", side_effect=forbidden
            ),
            patch.object(binder, "bind_raw_result", side_effect=forbidden),
            patch.object(binder, "_git", side_effect=forbidden),
        ):
            result = binder.validate_compact_result(CONFIG_RAW, encoded)
        replay = result["artifact_payload"]["independent_replay"]
        self.assertTrue(replay["complete_admission_passed"])
        self.assertEqual(len(replay["channels"]), 18)
        self.assertEqual(replay["failed_channels"], [])

    def test_tracked_result_validates_compactly_when_present(self) -> None:
        if not RESULT.is_file():
            self.skipTest("compact QA1 PREF1 result has not been materialized")
        encoded = RESULT.read_bytes()
        self.assertEqual(sha256(encoded).hexdigest(), binder.COMPACT_RESULT_SHA256)
        forbidden = AssertionError("tracked compact verification opened live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "_restore_shadow", side_effect=forbidden),
            patch.object(binder, "_git", side_effect=forbidden),
        ):
            binder.validate_compact_result(CONFIG_RAW, encoded)

    def test_claim_one_bit_classifier_inequality_and_channel_mutations_fail_closed(
        self,
    ) -> None:
        mutations = (
            (
                "claim",
                lambda value: value["artifact_payload"]["claims"].__setitem__(
                    "accepted_state", True
                ),
            ),
            (
                "one_bit_interval",
                lambda value: value["artifact_payload"]["independent_replay"][
                    "channels"
                ][0]["D12"]["upper"].__setitem__(
                    "numerator",
                    str(
                        int(
                            value["artifact_payload"]["independent_replay"]["channels"][
                                0
                            ]["D12"]["upper"]["numerator"]
                        )
                        + 1
                    ),
                ),
            ),
            (
                "classifier",
                lambda value: value["artifact_payload"]["independent_replay"][
                    "channels"
                ][3].__setitem__("classification", "order_inconclusive"),
            ),
            (
                "inequality",
                lambda value: value["artifact_payload"]["independent_replay"][
                    "channels"
                ][0].__setitem__(
                    "D01",
                    deepcopy(
                        value["artifact_payload"]["independent_replay"]["channels"][0][
                            "D12"
                        ]
                    ),
                ),
            ),
            (
                "channel_order",
                lambda value: value["artifact_payload"]["independent_replay"][
                    "channels"
                ].reverse(),
            ),
            (
                "all_of",
                lambda value: value["artifact_payload"][
                    "independent_replay"
                ].__setitem__("failed_channels", ["u:R"]),
            ),
        )
        for name, mutate in mutations:
            with self.subTest(name=name):
                changed = compact_fixture()
                mutate(changed)
                encoded = binder.canonical_result(changed)
                with self.assertRaises(binder.TDG10QA1PREF1Error):
                    with patch.object(
                        binder,
                        "COMPACT_RESULT_SHA256",
                        sha256(encoded).hexdigest(),
                    ):
                        binder.validate_compact_result(CONFIG_RAW, encoded)

    def test_duplicate_noncanonical_and_nonfinite_json_are_rejected(self) -> None:
        for name, raw in (
            ("duplicate", b'{"a":1,"a":2}\n'),
            ("noncanonical", b'{"a": 1}\n'),
            ("nan", b'{"a":NaN}\n'),
            ("infinity", b'{"a":Infinity}\n'),
        ):
            with self.subTest(name=name), self.assertRaises(binder.TDG10QA1PREF1Error):
                binder._json(raw, name)
        with self.assertRaises(binder.TDG10QA1PREF1Error):
            binder.canonical_result({"nonfinite": math.nan})

    def test_common_mode_import_boundary_excludes_qa1_decision_owners(self) -> None:
        source = Path(binder.__file__).read_text(encoding="utf-8")
        imports: set[str] = set()
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    imports.add(node.module)
                imports.update(alias.name for alias in node.names)
        forbidden = (
            "tdg10_qa1_authority",
            "tdg10_exact_complete_c_runtime",
            "tdg10_exact_complete_c_admission",
            "tdg9_ti2_authority",
            "tdg9_ti2_pref1_binder",
            "run_fgc_tdg10_qa1",
        )
        self.assertFalse(
            {name for name in imports if any(needle in name for needle in forbidden)}
        )
        self.assertNotIn("serialize_continuous_admission", imports)
        self.assertNotIn("validate_all_of_identity", imports)


class TDG10QA1PREF1GitIsolationTests(unittest.TestCase):
    def test_git_subprocess_strips_hostile_ambient_state(self) -> None:
        completed = subprocess.CompletedProcess(
            args=(), returncode=0, stdout=b"answer\n", stderr=b""
        )
        hostile = {
            "GIT_DIR": "/attacker/repository",
            "GIT_WORK_TREE": "/attacker/worktree",
            "GIT_CONFIG_GLOBAL": "/attacker/config",
            "GIT_CONFIG_SYSTEM": "/attacker/system-config",
            "GIT_CONFIG_COUNT": "1",
            "GIT_CONFIG_KEY_0": "include.path",
            "GIT_CONFIG_VALUE_0": "/attacker/redirect",
            "GIT_OBJECT_DIRECTORY": "/attacker/objects",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": "/attacker/alternate",
            "GIT_REPLACE_REF_BASE": "refs/attacker/replace/",
            "GIT_OPTIONAL_LOCKS": "1",
        }
        with (
            patch.dict(os.environ, hostile),
            patch.object(binder.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(
                binder._git(Path("/authorized"), "rev-parse", "HEAD"),
                b"answer\n",
            )
        command = run.call_args.args[0]
        self.assertEqual(
            command[:8],
            (
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                "rev-parse",
            ),
        )
        self.assertEqual(run.call_args.kwargs["cwd"], Path("/authorized"))
        self.assertIs(run.call_args.kwargs["stdin"], subprocess.DEVNULL)
        environment = run.call_args.kwargs["env"]
        self.assertEqual(environment["GIT_OPTIONAL_LOCKS"], "0")
        self.assertEqual(environment["LC_ALL"], "C")
        self.assertFalse(set(environment) & (set(hostile) - {"GIT_OPTIONAL_LOCKS"}))

    @staticmethod
    def _make_repository(root: Path, content: str) -> str:
        environment = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("GIT_")
        }
        environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
        subprocess.run(
            ["git", "init", "-q", str(root)],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )
        (root / "marker").write_text(content, encoding="utf-8")
        subprocess.run(
            ["git", "-C", str(root), "add", "marker"],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )
        subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "-c",
                "user.name=QA1 PREF1 test",
                "-c",
                "user.email=qa1-pref1@example.invalid",
                "commit",
                "-qm",
                "fixture",
            ],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )
        return (
            subprocess.run(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                check=True,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=environment,
            )
            .stdout.decode("ascii")
            .strip()
        )

    def test_hostile_git_paths_cannot_redirect_real_lookup(self) -> None:
        with TemporaryDirectory(prefix="qa1-pref1-hostile-git-") as directory:
            base = Path(directory)
            authorized = base / "authorized"
            attacker = base / "attacker"
            authorized_commit = self._make_repository(authorized, "authorized\n")
            self._make_repository(attacker, "attacker\n")
            attacker_config = base / "attacker.gitconfig"
            attacker_config.write_text("not valid git config\n", encoding="utf-8")
            hostile = {
                "GIT_DIR": str(attacker / ".git"),
                "GIT_WORK_TREE": str(attacker),
                "GIT_CONFIG_GLOBAL": str(attacker_config),
                "GIT_CONFIG_COUNT": "1",
                "GIT_CONFIG_KEY_0": "include.path",
                "GIT_CONFIG_VALUE_0": str(attacker_config),
                "GIT_OBJECT_DIRECTORY": str(attacker / ".git" / "objects"),
                "GIT_ALTERNATE_OBJECT_DIRECTORIES": str(attacker / ".git" / "objects"),
                "GIT_REPLACE_REF_BASE": "refs/attacker/replace/",
            }
            with patch.dict(os.environ, hostile):
                self.assertEqual(
                    binder._git(authorized, "rev-parse", "HEAD"),
                    f"{authorized_commit}\n".encode("ascii"),
                )
                self.assertEqual(
                    binder._git(authorized, "show", "HEAD:marker"), b"authorized\n"
                )


class TDG10QA1PREF1RawSafetyTests(unittest.TestCase):
    def test_canonical_temporary_four_leaf_tree_passes_hash_binding(self) -> None:
        with TemporaryDirectory(prefix="qa1-pref1-raw-control-") as directory:
            root = Path(directory)
            blobs = _write_raw_fixture(root)
            with (
                _raw_hash_patches(blobs),
                patch.object(
                    binder,
                    "_restore_shadow",
                    side_effect=AssertionError("raw fixture constructed shadows"),
                ),
            ):
                self.assertEqual(binder._raw_tree(root), blobs)

    def test_nofollow_hardlink_special_extra_missing_and_partial_attacks_fail_closed(
        self,
    ) -> None:
        attacks = (
            "root_symlink",
            "inner_directory_symlink",
            "leaf_symlink",
            "hardlink",
            "extra",
            "missing",
            "partial_nested",
            "special",
        )
        for attack in attacks:
            if attack == "special" and not hasattr(os, "mkfifo"):
                continue
            with (
                self.subTest(attack=attack),
                TemporaryDirectory(prefix=f"qa1-pref1-{attack}-") as directory,
            ):
                root = Path(directory)
                blobs = _write_raw_fixture(root)
                raw = root / binder.RAW_NAMESPACE
                if attack == "root_symlink":
                    real = root / "real-raw"
                    raw.rename(real)
                    raw.symlink_to(real, target_is_directory=True)
                elif attack == "inner_directory_symlink":
                    nested = raw / "fine-endpoint"
                    real = root / "real-fine-endpoint"
                    nested.rename(real)
                    nested.symlink_to(real, target_is_directory=True)
                elif attack == "leaf_symlink":
                    path = raw / "manifest.json"
                    target = root / "manifest-target.json"
                    path.rename(target)
                    path.symlink_to(target)
                elif attack == "hardlink":
                    path = raw / "terminal.json"
                    target = root / "terminal-target.json"
                    path.rename(target)
                    os.link(target, path)
                elif attack == "extra":
                    (raw / "foreign").write_bytes(b"foreign\n")
                elif attack == "missing":
                    (raw / "terminal.json").unlink()
                elif attack == "partial_nested":
                    (raw / binder.DESCRIPTOR_RELATIVE).unlink()
                elif attack == "special":
                    path = raw / "terminal.json"
                    path.unlink()
                    try:
                        os.mkfifo(path)
                    except OSError as exc:
                        self.skipTest(
                            f"FIFO unsupported by temporary filesystem: {exc}"
                        )
                with (
                    _raw_hash_patches(blobs),
                    patch.object(
                        binder,
                        "_restore_shadow",
                        side_effect=AssertionError("raw attack constructed shadows"),
                    ),
                    self.assertRaises((binder.TDG10QA1PREF1Error, OSError)),
                ):
                    binder._raw_tree(root)

    def test_leaf_identity_race_fails_before_read_or_shadow_work(self) -> None:
        with TemporaryDirectory(prefix="qa1-pref1-race-") as directory:
            root = Path(directory)
            path = root / "leaf"
            path.write_bytes(b"race fixture\n")
            directory_fd = os.open(root, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
            try:
                with (
                    patch.object(
                        binder,
                        "_identity",
                        side_effect=[(1, 1, 1, 1, 1), (2, 2, 2, 2, 2)],
                    ),
                    patch.object(
                        binder,
                        "_restore_shadow",
                        side_effect=AssertionError("race fixture constructed shadows"),
                    ),
                    self.assertRaises(binder.TDG10QA1PREF1Error) as raised,
                ):
                    binder._read_leaf_at(
                        directory_fd, "leaf", maximum=1024, label="race leaf"
                    )
            finally:
                os.close(directory_fd)
        self.assertEqual(raised.exception.stop_id, "PREF1_LEAF_RACED")


class TDG10QA1PREF1AuthorityAndReproducerTests(unittest.TestCase):
    def test_authority_rejects_anything_other_than_twenty_blobs(self) -> None:
        def fake_git(_root: Path, *arguments: str) -> bytes:
            if arguments == ("rev-parse", "--show-toplevel"):
                return f"{ROOT}\n".encode("utf-8")
            if arguments == ("cat-file", "-t", binder.AUTHORITY_COMMIT):
                return b"commit\n"
            if arguments == (
                "show",
                "-s",
                "--format=%P",
                binder.AUTHORITY_COMMIT,
            ):
                return f"{binder.AUTHORITY_PARENT}\n".encode("ascii")
            raise AssertionError(arguments)

        with (
            patch.object(binder, "_AUTHORITY_BLOBS", binder._AUTHORITY_BLOBS[:-1]),
            patch.object(binder, "_git", side_effect=fake_git),
            self.assertRaises(binder.TDG10QA1PREF1Error) as raised,
        ):
            binder._authenticate_authority(
                ROOT, {"authority_commit": binder.AUTHORITY_COMMIT}
            )
        self.assertEqual(raised.exception.stop_id, "PREF1_AUTHORITY_BLOB_COUNT_DRIFT")

    def test_default_reproducer_missing_result_fails_without_loading_or_live_build(
        self,
    ) -> None:
        with TemporaryDirectory(prefix="qa1-pref1-missing-result-") as directory:
            output = Path(directory) / "missing.json"
            with (
                patch.object(
                    reproducer,
                    "load_binder",
                    side_effect=AssertionError("missing default path loaded binder"),
                ),
                patch.object(
                    sys,
                    "argv",
                    ["reproduce_fgc_tdg10_qa1_pref1.py", "--output", str(output)],
                ),
                self.assertRaises(SystemExit) as raised,
            ):
                reproducer.main()
            self.assertIn("compact result is absent", str(raised.exception))
            self.assertIn("never runs live binding", str(raised.exception))
            self.assertFalse(output.exists())

    def test_default_reproducer_validates_existing_compact_without_live_build(
        self,
    ) -> None:
        with TemporaryDirectory(prefix="qa1-pref1-compact-result-") as directory:
            output = Path(directory) / "compact.json"
            encoded = binder.canonical_result(compact_fixture())
            output.write_bytes(encoded)
            with (
                patch.object(reproducer, "load_binder", return_value=binder),
                patch.object(
                    binder, "COMPACT_RESULT_SHA256", sha256(encoded).hexdigest()
                ),
                patch.object(
                    binder,
                    "build_pref1_result",
                    side_effect=AssertionError("default path ran live build"),
                ),
                patch.object(
                    sys,
                    "argv",
                    ["reproduce_fgc_tdg10_qa1_pref1.py", "--output", str(output)],
                ),
            ):
                self.assertEqual(reproducer.main(), 0)


if __name__ == "__main__":
    unittest.main()
