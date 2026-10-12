"""Focused prospective authority controls for FGC-1-HLT17-SRCQ1-REC1."""

from __future__ import annotations

import ast
from contextlib import ExitStack, contextmanager
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution import hlt17_srcq1_rec1_auth1 as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_srcq1_rec1 import (  # noqa: E402
    ATTEMPT2_DERIVATION_DOCUMENT,
    FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER,
    ORIGIN_CAPTURE_SHA256,
)
from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS  # noqa: E402

MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1_auth1.py"
SCRIPT_PATH = ROOT / "scripts/reproduce_fgc_hlt17_srcq1_rec1_frz1.py"
CLI_PATH = ROOT / "scripts/qualify_fgc_hlt17_srcq1_rec1.py"
CAPS = (
    ("RK4-2049", "0x1.0000000000000p-4", "0x1.33345e70b5188p+0", "0x1.aaa90b0fb5c26p-8"),
    ("RK4-4097", "0x1.0000000000000p-5", "0x1.3333333aecf3ep+0", "0x1.aaaaaa9fefc9cp-9"),
    ("RK4-8193", "0x1.0000000000000p-6", "0x1.3333333333653p+0", "0x1.aaaaaaaaaa654p-10"),
    ("SSPRK3-4097", "0x1.0000000000000p-5", "0x1.3333333fe51c4p+0", "0x1.aaaaaa9908e70p-9"),
    ("SSPRK3-8193", "0x1.0000000000000p-6", "0x1.33333333349e5p+0", "0x1.aaaaaaaaa8b26p-10"),
    ("SSPRK3-16385", "0x1.0000000000000p-7", "0x1.3333333333346p+0", "0x1.aaaaaaaaaaa91p-11"),
)


def _git_env() -> dict[str, str]:
    env = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    env.update(
        {
            "GIT_AUTHOR_EMAIL": "srcq1-rec1@example.test",
            "GIT_AUTHOR_NAME": "srcq1-rec1-test",
            "GIT_COMMITTER_EMAIL": "srcq1-rec1@example.test",
            "GIT_COMMITTER_NAME": "srcq1-rec1-test",
            "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_CONFIG_NOSYSTEM": "1",
        }
    )
    return env


def _git(repo: Path, args: tuple[str, ...]) -> None:
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", *args],
        cwd=repo,
        check=True,
        env=_git_env(),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _git_head(repo: Path) -> str:
    return (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, env=_git_env())
        .decode("ascii")
        .strip()
    )


def _write(root: Path, relative: str, payload: bytes) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _fixture_environment() -> dict[str, str]:
    environment = {key: "synthetic" for key in authority.ENVIRONMENT_KEYS}
    environment["executable_sha256"] = "a" * 64
    environment["numpy_extension_sha256"] = "b" * 64
    return environment


def _fixture_hashes() -> dict[str, str]:
    return {
        "implementation_sha256": "c" * 64,
        "runner_sha256": "d" * 64,
        "authority_sha256": "e" * 64,
    }


def _fixture_closure() -> dict[str, object]:
    return {
        "sha256": "f" * 64,
        "file_count": 3,
        "files": (("src/recursive_horizons/x.py", "1" * 64),),
        "executable_sha256": "a" * 64,
        "numpy_extension_sha256": "b" * 64,
    }


def _cap_config():
    return {
        "member_cap": [
            {
                "member_key": key,
                "grid_spacing_hex": spacing,
                "previous_speed_upper_hex": speed,
                "cfl_maximum_hex": "0x1.0000000000000p-3",
                "requested_cap_hex": cap,
            }
            for key, spacing, speed, cap in CAPS
        ]
    }


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "reproduce_fgc_hlt17_srcq1_rec1_frz1", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_cli():
    spec = importlib.util.spec_from_file_location(
        "qualify_fgc_hlt17_srcq1_rec1_cli_auth", CLI_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _freeze_repo(
    temporary: str, *, extra: str | None = None
) -> tuple[Path, str, str, dict[str, object]]:
    repo = Path(temporary).resolve()
    subprocess.run(
        ["git", "-c", "init.defaultBranch=main", "init"],
        cwd=repo,
        check=True,
        env=_git_env(),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    _write(repo, "seed.txt", b"seed\n")
    _git(repo, ("add", "seed.txt"))
    _git(repo, ("commit", "-m", "base"))
    base = _git_head(repo)
    hashes = _fixture_hashes()
    environment = _fixture_environment()
    config = authority.expected_config(
        implementation_commit=base,
        implementation_hashes=hashes,
        source_closure_sha256="f" * 64,
        environment=environment,
    )
    raw = authority.render_config(config)
    _write(repo, authority.CONFIG_PATH, raw)
    _write(repo, authority.OWNER_PATH, b"owner\n")
    _write(repo, authority.MAKE_PATH, b"make\n")
    added = [authority.CONFIG_PATH, authority.OWNER_PATH, authority.MAKE_PATH]
    if extra is not None:
        _write(repo, extra, b"extra\n")
        added.append(extra)
    _git(repo, ("add", *added))
    _git(repo, ("commit", "-m", "freeze"))
    return repo, base, _git_head(repo), config


class SRCQ1REC1AuthorityImplementationTests(unittest.TestCase):
    def test_cap_mapping_is_exact_and_rejects_a_swapped_value(self) -> None:
        parsed = authority._cap_mapping(_cap_config())
        self.assertEqual(tuple(parsed), tuple(row[0] for row in CAPS))
        self.assertEqual(parsed["RK4-2049"], CAPS[0][3])
        self.assertEqual(
            parsed["SSPRK3-4097"],
            FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"],
        )
        altered = _cap_config()
        altered["member_cap"][0]["requested_cap_hex"] = CAPS[2][3]
        with self.assertRaisesRegex(authority.HLT17SRCQ1REC1AuthorityError, "formula"):
            authority._cap_mapping(altered)

    def test_environment_binds_executable_numpy_extension_and_kernel(self) -> None:
        observed = authority.environment_identity()
        for name in (
            "executable_sha256",
            "numpy_extension_sha256",
            "kernel_release",
            "python_version",
            "numpy_version",
            "machine",
            "executable_realpath",
            "numpy_extension_realpath",
        ):
            self.assertTrue(observed[name])
        self.assertEqual(len(observed["executable_sha256"]), 64)
        self.assertEqual(len(observed["numpy_extension_sha256"]), 64)
        self.assertEqual(tuple(observed), authority.ENVIRONMENT_KEYS)

    def test_software_same_size_race_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary).resolve() / "software.bin"
            path.write_bytes(b"AAAA")
            original_read = authority.os.read
            changed = False

            def raced_read(descriptor: int, size: int) -> bytes:
                nonlocal changed
                data = original_read(descriptor, size)
                if not changed:
                    changed = True
                    path.write_bytes(b"BBBB")
                return data

            with patch.object(authority.os, "read", side_effect=raced_read):
                with self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError,
                    "software path changed during read",
                ):
                    authority._hash_software_file(str(path))

    def test_source_closure_covers_committed_package_and_runner(self) -> None:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        record = authority.source_closure_identity(ROOT, head)
        self.assertGreater(record["file_count"], 100)
        paths = {item[0] for item in record["files"]}
        self.assertIn(authority.RUNNER_PATH, paths)
        self.assertIn(authority.HARNESS_PATH, paths)
        self.assertEqual(len(record["sha256"]), 64)

    def test_sealed_and_metadata_linked_compacts_differ_by_derivation_document(
        self,
    ) -> None:
        record = authority.authenticate_attempt2_compacts(ROOT)
        self.assertIs(record["sealed_historical_blob"]["derivation_document"], None)
        self.assertEqual(
            record["metadata_linked_blob"]["derivation_document"],
            ATTEMPT2_DERIVATION_DOCUMENT,
        )
        compact = record["compact"]
        self.assertEqual(compact["derivation_document"], ATTEMPT2_DERIVATION_DOCUMENT)
        self.assertEqual(
            record["sealed_historical_blob"]["sha256"],
            "7860491360fe058d6c5019f780103e6bf13b2c243d3c1d59ec3f74109636ab5f",
        )
        self.assertEqual(
            record["metadata_linked_blob"]["sha256"],
            "ba656a63343a87b389e388e39e142e8cfc94194c96c33326f6da0963f4504c48",
        )

    def test_raw_directory_digest_requires_exactly_one_nofollow_leaf(self) -> None:
        payload = {
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
            "source_configuration": dict(
                authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER
            ),
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
        digest = sha256(
            canonical_json_bytes({"qualification.json": sha256(raw).hexdigest()})
        ).hexdigest()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            _write(root, authority.ATTEMPT2_RAW_LEAF_RELATIVE, raw)
            with (
                patch.object(authority, "ATTEMPT2_RAW_LEAF_BYTES", len(raw)),
                patch.object(
                    authority, "ATTEMPT2_RAW_LEAF_SHA256", sha256(raw).hexdigest()
                ),
                patch.object(
                    authority, "ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256", digest
                ),
            ):
                record = authority.authenticate_attempt2_raw(root)
            self.assertEqual(record["leaf_count"], 1)
            self.assertEqual(
                record["source_configuration"],
                dict(authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER),
            )
            _write(root, f"{authority.ATTEMPT2_RAW_DIRECTORY}/extra.json", b"{}")
            with (
                patch.object(authority, "ATTEMPT2_RAW_LEAF_BYTES", len(raw)),
                patch.object(
                    authority, "ATTEMPT2_RAW_LEAF_SHA256", sha256(raw).hexdigest()
                ),
                patch.object(
                    authority, "ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256", digest
                ),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "leaf count"
                ),
            ):
                authority.authenticate_attempt2_raw(root)

    def test_raw_leaf_same_size_race_is_rejected(self) -> None:
        payload = {
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
            "source_configuration": dict(
                authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER
            ),
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
            "ascii"
        )
        digest = sha256(
            canonical_json_bytes({"qualification.json": sha256(raw).hexdigest()})
        ).hexdigest()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            leaf = root / authority.ATTEMPT2_RAW_LEAF_RELATIVE
            _write(root, authority.ATTEMPT2_RAW_LEAF_RELATIVE, raw)
            original_read = authority.os.read
            changed = False

            def raced_read(descriptor: int, size: int) -> bytes:
                nonlocal changed
                data = original_read(descriptor, size)
                if not changed:
                    changed = True
                    leaf.write_bytes(b"X" * len(raw))
                return data

            with (
                patch.object(authority, "ATTEMPT2_RAW_LEAF_BYTES", len(raw)),
                patch.object(
                    authority, "ATTEMPT2_RAW_LEAF_SHA256", sha256(raw).hexdigest()
                ),
                patch.object(
                    authority, "ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256", digest
                ),
                patch.object(authority.os, "read", side_effect=raced_read),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError,
                    "raw leaf changed during read",
                ),
            ):
                authority.authenticate_attempt2_raw(root)

    def test_output_absence_rejects_existing_staging_and_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / authority.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{authority._STAGING_PREFIX}dead").mkdir()
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "staging namespace already exists"
            ):
                authority.require_output_absent(root)
        with (
            tempfile.TemporaryDirectory() as temporary,
            tempfile.TemporaryDirectory() as foreign,
        ):
            root = Path(temporary).resolve()
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (root / authority.OUTPUT_NAMESPACE).symlink_to("missing-target")
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / authority.PRO20_EVENT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "PRO20 event namespace"
            ):
                authority.require_pro20_namespace_absent(root)

    def test_process_preflight_excludes_self_and_rejects_foreign_runners(self) -> None:
        self_pid = os.getpid()
        parent = os.getppid()
        rows = (
            (
                self_pid,
                parent,
                "python -m unittest tests/test_fgc_hlt17_srcq1_rec1_auth1",
            ),
            (parent, 1, "/sbin/launchd"),
            (4242, 1, "python scripts/run_fgc_tdg10_qa2.py --run"),
        )
        with patch.object(authority, "_process_rows", return_value=rows):
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "scientific runner"
            ):
                authority._require_process_preflight()
        safe = (
            (
                self_pid,
                parent,
                "python -m unittest tests/test_fgc_hlt17_srcq1_rec1_auth1",
            ),
            (parent, 1, "/sbin/launchd"),
        )
        with patch.object(authority, "_process_rows", return_value=safe):
            authority._require_process_preflight()
        review = safe + ((99, 1, "codex review --wait"),)
        with patch.object(authority, "_process_rows", return_value=review):
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "review"
            ):
                authority._require_process_preflight()
        pytest_row = safe + ((100, 1, "python -m pytest tests"),)
        with patch.object(authority, "_process_rows", return_value=pytest_row):
            with self.assertRaisesRegex(
                authority.HLT17SRCQ1REC1AuthorityError, "test partition"
            ):
                authority._require_process_preflight()

    @contextmanager
    def _authorize_after_git(self, environment: dict[str, str], hashes: dict[str, str]):
        compact = {
            "sealed_historical_blob": {
                "commit": "efb53b839a36a7f7f52018a44ef4eeaee42c1752",
                "sha256": "7860491360fe058d6c5019f780103e6bf13b2c243d3c1d59ec3f74109636ab5f",
                "derivation_document": None,
            },
            "metadata_linked_blob": {
                "commit": "4891feab58922e449a06585a848eb6a11df994ea",
                "sha256": "ba656a63343a87b389e388e39e142e8cfc94194c96c33326f6da0963f4504c48",
                "derivation_document": "docs/fgc-hlt17-srcq1-frz1.md",
            },
            "compact": {},
        }
        raw = {
            "raw_leaf_relative": authority.ATTEMPT2_RAW_LEAF_RELATIVE,
            "raw_leaf_sha256": "889a8d17551a1eea6ed56f43f602fefb3f61ab78637bbdd0c995905e22259d89",
            "raw_directory_payload_sha256": (
                "7b06deb68c989734f8dd669f2f1a8f4531aa90b82cc0e9ff34330db363123009"
            ),
            "raw_leaf_bytes": 45589,
            "source_configuration": dict(
                authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER
            ),
            "leaf_count": 1,
        }
        with ExitStack() as stack:
            stack.enter_context(
                patch.object(authority, "environment_identity", return_value=environment)
            )
            stack.enter_context(
                patch.object(
                    authority, "source_closure_identity", return_value=_fixture_closure()
                )
            )
            stack.enter_context(
                patch.object(
                    authority, "_committed_implementation_hashes", return_value=hashes
                )
            )
            stack.enter_context(
                patch.object(
                    authority, "_live_implementation_hashes", return_value=hashes
                )
            )
            stack.enter_context(
                patch.object(
                    authority, "authenticate_attempt2_compacts", return_value=compact
                )
            )
            stack.enter_context(
                patch.object(authority, "authenticate_attempt2_raw", return_value=raw)
            )
            stack.enter_context(patch.object(authority, "_require_process_preflight"))
            stack.enter_context(patch.object(authority, "require_output_absent"))
            stack.enter_context(
                patch.object(authority, "require_pro20_namespace_absent")
            )
            yield

    def test_authorize_accepts_direct_successor_and_returns_mapping_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head, config = _freeze_repo(temporary)
            environment = _fixture_environment()
            hashes = _fixture_hashes()
            with self._authorize_after_git(environment, hashes):
                receipt = authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        self.assertEqual(receipt["authority_commit"], head)
        self.assertEqual(receipt["implementation_commit"], base)
        self.assertEqual(
            receipt["implementation_sha256"], hashes["implementation_sha256"]
        )
        self.assertEqual(receipt["source_closure_sha256"], "f" * 64)
        self.assertEqual(receipt["environment"], environment)
        self.assertEqual(
            receipt["prospective_requested_cap_hex_by_member"],
            dict(FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER),
        )
        self.assertEqual(
            receipt["source_configuration"],
            dict(authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER),
        )
        self.assertEqual(receipt["output_namespace"], authority.OUTPUT_NAMESPACE)
        self.assertEqual(receipt["origin_capture_sha256"], ORIGIN_CAPTURE_SHA256)
        self.assertIs(receipt["boolean_flag_accepted"], False)
        self.assertIs(receipt["campaign_execution_authorized"], False)
        self.assertIs(receipt["accepted_state_advanced"], False)
        self.assertIs(receipt["endpoint_adopted_or_serialized"], False)
        self.assertIs(receipt["campaign_store_written"], False)
        self.assertEqual(config["artifact_id"], authority.ARTIFACT_ID)

    def test_altered_parent_merge_delta_dirty_and_hidden_are_rejected(self) -> None:
        environment = _fixture_environment()
        hashes = _fixture_hashes()
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "authority commit"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit="0" * 40,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, _head, _config = _freeze_repo(temporary)
            _git(repo, ("checkout", "-b", "side"))
            _write(repo, "seed.txt", b"side\n")
            _git(repo, ("commit", "-am", "side"))
            _git(repo, ("checkout", "main"))
            _git(repo, ("merge", "side", "--no-ff", "-m", "merge"))
            merge = _git_head(repo)
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "one implementation parent"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=merge,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary, extra="bonus.txt")
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "delta paths"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            _write(repo, authority.OWNER_PATH, b"dirty\n")
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "clean worktree"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            _git(repo, ("update-index", "--skip-worktree", authority.OWNER_PATH))
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "clean worktree"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            with (
                self._authorize_after_git(environment, hashes),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "exact commit"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit="HEAD",
                    implementation_sha256=hashes["implementation_sha256"],
                )

    def test_environment_source_raw_compact_output_and_process_drift_are_rejected(
        self,
    ) -> None:
        environment = _fixture_environment()
        hashes = _fixture_hashes()
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            other = dict(environment)
            other["python_version"] = "drifted"
            with (
                self._authorize_after_git(environment, hashes),
                patch.object(authority, "environment_identity", return_value=other),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "environment"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            drifted = dict(_fixture_closure())
            drifted["sha256"] = "0" * 64
            with (
                self._authorize_after_git(environment, hashes),
                patch.object(
                    authority, "source_closure_identity", return_value=drifted
                ),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "source closure"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            drifted_raw = {
                "raw_leaf_relative": authority.ATTEMPT2_RAW_LEAF_RELATIVE,
                "raw_leaf_sha256": "0" * 64,
                "raw_directory_payload_sha256": "1" * 64,
                "raw_leaf_bytes": 1,
                "source_configuration": {"RK4-2049": "2" * 64},
                "leaf_count": 1,
            }
            with (
                self._authorize_after_git(environment, hashes),
                patch.object(
                    authority, "authenticate_attempt2_raw", return_value=drifted_raw
                ),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "source-configuration"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            with (
                self._authorize_after_git(environment, hashes),
                patch.object(
                    authority,
                    "require_output_absent",
                    side_effect=authority.HLT17SRCQ1REC1AuthorityError(
                        "REC1 output namespace already exists"
                    ),
                ),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "output namespace"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head, _config = _freeze_repo(temporary)
            with (
                self._authorize_after_git(environment, hashes),
                patch.object(
                    authority,
                    "_require_process_preflight",
                    side_effect=authority.HLT17SRCQ1REC1AuthorityError(
                        "another scientific runner, test partition, or review is active"
                    ),
                ),
                self.assertRaisesRegex(
                    authority.HLT17SRCQ1REC1AuthorityError, "scientific runner"
                ),
            ):
                authority.validate_authority_delta(
                    repository_root=repo,
                    authority_commit=head,
                    implementation_sha256=hashes["implementation_sha256"],
                )

    def test_malformed_config_and_noncommit_receipts_are_rejected(self) -> None:
        with self.assertRaisesRegex(
            authority.HLT17SRCQ1REC1AuthorityError, "malformed"
        ):
            authority._parse_config(b"not = toml [")
        with self.assertRaisesRegex(
            authority.HLT17SRCQ1REC1AuthorityError, "fields differ"
        ):
            authority._parse_config(b"artifact_id = \"nope\"\n")
        with self.assertRaises(authority.HLT17SRCQ1REC1AuthorityError):
            authority.validate_authority_delta(
                repository_root=ROOT,
                authority_commit="HEAD",
                implementation_sha256="0" * 64,
            )
        with self.assertRaises(authority.HLT17SRCQ1REC1AuthorityError):
            authority.validate_authority_delta(
                repository_root=ROOT,
                authority_commit="0" * 40,
                implementation_sha256="0" * 64,
            )

    def test_emit_config_is_deterministic_and_does_not_write(self) -> None:
        environment = _fixture_environment()
        hashes = _fixture_hashes()
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, _head, _config = _freeze_repo(temporary)
            before = {
                path: (repo / path).read_bytes()
                for path in (
                    authority.CONFIG_PATH,
                    authority.OWNER_PATH,
                    authority.MAKE_PATH,
                )
            }
            with (
                patch.object(
                    authority, "environment_identity", return_value=environment
                ),
                patch.object(
                    authority, "source_closure_identity", return_value=_fixture_closure()
                ),
                patch.object(
                    authority, "_live_implementation_hashes", return_value=hashes
                ),
            ):
                first = authority.emit_config_bytes(
                    repo, implementation_commit=base
                )
                second = authority.emit_config_bytes(
                    repo, implementation_commit=base
                )
            self.assertEqual(first, second)
            self.assertTrue(first.startswith(b"artifact_id = "))
            self.assertEqual(
                (repo / authority.CONFIG_PATH).read_bytes(), before[authority.CONFIG_PATH]
            )
            parsed = authority._parse_config(first)
            self.assertEqual(parsed["implementation_commit"], base)
            self.assertEqual(parsed["origin_capture_sha256"], ORIGIN_CAPTURE_SHA256)
            with (
                patch.object(
                    authority, "environment_identity", return_value=environment
                ),
                patch.object(
                    authority, "source_closure_identity", return_value=_fixture_closure()
                ),
                patch.object(
                    authority, "_live_implementation_hashes", return_value=hashes
                ),
            ):
                _write(repo, authority.CONFIG_PATH, first)
                record = authority.check_tracked_config(repo)
            self.assertIs(record["authorized"], False)
            self.assertIs(record["executed"], False)
            self.assertIs(record["physical_source_qualification_executed"], False)

    def test_reproducer_emits_without_writing_or_authorizing(self) -> None:
        script = _load_script()
        emitted = []

        def fake_emit(root):
            emitted.append(root)
            return b'artifact_id = "FGC-1-HLT17-SRCQ1-REC1-FRZ1"\n'

        with (
            patch.object(authority, "emit_config_bytes", side_effect=fake_emit),
            patch.object(
                authority,
                "check_tracked_config",
                side_effect=AssertionError("authorize"),
            ),
            patch.object(
                authority,
                "validate_authority_delta",
                side_effect=AssertionError("authorize"),
            ),
        ):
            code = script.main(["--emit-config"])
        self.assertEqual(code, 0)
        self.assertEqual(emitted, [script.ROOT])

    def test_status_and_cli_make_no_origin_or_source_call(self) -> None:
        cli = _load_cli()
        with (
            patch.object(
                cli,
                "capture_pro20_historical_origin",
                side_effect=AssertionError("origin"),
            ),
            patch.object(
                cli,
                "read_static_factory_inputs",
                side_effect=AssertionError("factory"),
            ),
            patch.object(
                cli, "qualify_and_publish", side_effect=AssertionError("qualify")
            ),
            patch.object(
                cli, "observe_environment", side_effect=AssertionError("env")
            ),
        ):
            self.assertEqual(cli.main([]), 0)
            self.assertEqual(cli.main(["--status"]), 0)
            self.assertEqual(cli.main(["--status", "--authority-commit", "a" * 40]), 1)
            self.assertEqual(cli.main(["--authority-commit", "a" * 40]), 2)

    def test_cli_real_path_passes_authority_pinned_expectations(self) -> None:
        cli = _load_cli()
        captured: dict[str, object] = {}
        identity = {
            "implementation_id": "tdg11_c1r1_integer_exponent_direct_ring_v1",
            "mathematical_object": "exact_accumulation_reconstruction_with_debit",
            "reference_wire_artifact_id": "FGC-1-TDG11-IMP1",
            "reference_wire_evaluator_id": (
                "tdg11_imp1_exact_bernstein_with_dual_fallback_v1"
            ),
            "admission_runtime_id": "hlt17_c1r1_admission_runtime_v1",
        }
        environment = {"python_implementation": "CPython"}
        configs = {key: "a" * 64 for key in MEMBER_KEYS}

        def fake_require(**kwargs):
            return {
                "authority_commit": kwargs["authority_commit"],
                "implementation_sha256": kwargs["implementation_sha256"],
                "source_closure_sha256": "c" * 64,
                "environment": environment,
                "implementation_identity": identity,
                "source_configuration": configs,
                "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
                "output_namespace": "/tmp/srcq1-rec1-out",
                "prospective_requested_cap_hex_by_member": dict(
                    FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER
                ),
            }

        def fake_publish(*args, **kwargs):
            captured.update(kwargs)
            return {"publication": {"payload_sha256": "e" * 64}}

        with (
            patch.object(cli, "require_authority_delta", side_effect=fake_require),
            patch.object(cli, "qualify_and_publish", side_effect=fake_publish),
            patch.object(
                cli,
                "implementation_identity",
                return_value={"sha256": "b" * 64},
            ),
            patch.object(
                cli,
                "capture_pro20_historical_origin",
                side_effect=AssertionError("origin"),
            ),
        ):
            code = cli.main(
                [
                    "--authority-commit",
                    "a" * 40,
                    "--output-directory",
                    "/tmp/srcq1-rec1-out",
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(captured["expected_environment"], environment)
        self.assertEqual(captured["expected_implementation"], identity)
        self.assertEqual(captured["expected_source_configurations"], configs)
        self.assertEqual(captured["environment"], environment)
        self.assertEqual(captured["expected_origin_sha256"], ORIGIN_CAPTURE_SHA256)

    def test_neutral_git_and_no_historical_authority_imports(self) -> None:
        tree = ast.parse(MODULE_PATH.read_text("utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
                imported.extend(alias.name for alias in node.names)
        joined = " ".join(imported)
        self.assertNotIn("hlt17_srcq1_auth1", joined)
        self.assertIn("ResolveCommit", joined)
        self.assertIn("CommitParents", joined)
        self.assertIn("InspectDelta", joined)
        self.assertIn("InspectWorktree", joined)
        self.assertIn("InspectTree", joined)
        self.assertIn("ReadBlob", joined)
        source = MODULE_PATH.read_text("utf-8")
        self.assertNotIn("rev-parse", source)
        self.assertNotIn("ls-tree", source)
        self.assertNotIn("diff --name-only", source)
        self.assertNotIn("hlt17_srcq1_auth1", authority.__dict__)


if __name__ == "__main__":
    unittest.main()
