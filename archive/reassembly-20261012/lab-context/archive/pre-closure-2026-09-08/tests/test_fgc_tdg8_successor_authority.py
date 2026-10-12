"""Adversarial Git/image tests for the one-event TDG8 successor authority."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg8_successor_authority as authority


SOURCE_ROOT = Path(__file__).resolve().parents[1]
REPRODUCER_PATH = SOURCE_ROOT / "scripts/reproduce_fgc_tdg8_run1_auth1.py"


def git(root: Path, *args: str, input_raw: bytes | None = None) -> str:
    result = subprocess.run(
        ("git", "-C", str(root), *args),
        input=input_raw,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return result.stdout.decode("utf-8").strip()


def write(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)


def toml_scalar(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    raise AssertionError(f"unsupported fixture value: {value!r}")


def append_table(lines: list[str], name: str, values: dict[str, object]) -> None:
    lines.append(f"[{name}]")
    lines.extend(f"{key} = {toml_scalar(value)}" for key, value in values.items())
    lines.append("")


def config_bytes(root: Path, predecessor: str) -> bytes:
    lines = [
        f"schema_version = {authority.SCHEMA_VERSION}",
        f"artifact_id = {toml_scalar(authority.ARTIFACT_ID)}",
        f"project_version = {toml_scalar(authority.PROJECT_VERSION)}",
        f"target_protocol = {toml_scalar(authority.TARGET_PROTOCOL)}",
        f"classification = {toml_scalar(authority.CLASSIFICATION)}",
        "nonclaims = [",
        *(f"  {toml_scalar(item)}," for item in authority.NONCLAIMS),
        "]",
        "",
    ]
    append_table(
        lines,
        "predecessor",
        {
            "repair_commit": predecessor,
            "must_be_ancestor_of_authority_commit": True,
        },
    )
    for item in authority.IMMUTABLE_BINDINGS:
        lines.append("[[immutable_bindings]]")
        lines.extend(
            (
                f"name = {toml_scalar(item.name)}",
                f"artifact_id = {toml_scalar(item.artifact_id)}",
                f"path = {toml_scalar(item.path)}",
                f"sha256 = {toml_scalar(item.sha256)}",
                "",
            )
        )
    for role, relative in authority.REQUIRED_IMPLEMENTATION_INVENTORY:
        digest = sha256((root / relative).read_bytes()).hexdigest()
        lines.append("[[implementation_inventory]]")
        lines.extend(
            (
                f"role = {toml_scalar(role)}",
                f"path = {toml_scalar(relative)}",
                f"sha256 = {toml_scalar(digest)}",
                "",
            )
        )
    append_table(lines, "environment", authority.PINNED_ENVIRONMENT)
    append_table(lines, "progression", authority.PROGRESSION_CONTRACT)
    append_table(lines, "stores", authority.STORES_CONTRACT)
    append_table(lines, "scope", authority.SCOPE_CONTRACT)
    append_table(lines, "claims", authority.CLAIMS_CONTRACT)
    return ("\n".join(lines).rstrip() + "\n").encode("utf-8")


@dataclass
class RepositoryFixture:
    temporary: tempfile.TemporaryDirectory[str]
    root: Path
    predecessor: str
    config_raw: bytes
    result_raw: bytes | None
    authority_commit: str | None


class TDG8SuccessorAuthorityTests(unittest.TestCase):
    def make_repository(
        self,
        *,
        finalize: bool = True,
        missing_delta_path: str | None = None,
    ) -> RepositoryFixture:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        git(root, "init", "-q")
        git(root, "config", "user.email", "tdg8@example.invalid")
        git(root, "config", "user.name", "TDG8 test")

        # Immutable predecessor evidence and the one unchanged critical runtime
        # leaf already exist at the repaired parent.  They must not appear in
        # the prospective successor delta.
        for binding in authority.IMMUTABLE_BINDINGS:
            write(root / binding.path, (SOURCE_ROOT / binding.path).read_bytes())
        for role, relative in authority.REQUIRED_IMPLEMENTATION_INVENTORY:
            if relative not in authority.AUTHORITY_DELTA_PATHS:
                write(root / relative, f"{role}: {relative}\n".encode("utf-8"))
        write(root / "repair.txt", b"representation repair\n")
        git(root, "add", ".")
        git(root, "commit", "-qm", "repair")
        predecessor = git(root, "rev-parse", "HEAD")

        inventory = {
            relative: role
            for role, relative in authority.REQUIRED_IMPLEMENTATION_INVENTORY
        }
        for relative in authority.AUTHORITY_DELTA_PATHS:
            if relative in {authority.CONFIG_PATH, authority.RESULT_PATH}:
                continue
            if relative == missing_delta_path:
                continue
            role = inventory.get(relative, "successor-owner")
            write(root / relative, f"{role}: {relative}\n".encode("utf-8"))

        config_raw = config_bytes(root, predecessor)
        write(root / authority.CONFIG_PATH, config_raw)
        git(root, "add", ".")
        git(root, "commit", "-qm", "prospective implementation")

        if not finalize:
            return RepositoryFixture(
                temporary=temporary,
                root=root,
                predecessor=predecessor,
                config_raw=config_raw,
                result_raw=None,
                authority_commit=None,
            )

        with (
            patch.object(authority, "PREDECESSOR_REPAIR_COMMIT", predecessor),
            patch.object(
                authority,
                "observed_environment",
                return_value=dict(authority.PINNED_ENVIRONMENT),
            ),
        ):
            result = authority.build_prelaunch(config_raw, root)
        result_raw = authority.canonical_result(result)
        write(root / authority.RESULT_PATH, result_raw)
        git(root, "add", authority.RESULT_PATH)
        git(root, "commit", "-qm", "committed successor authority")
        authority_commit = git(root, "rev-parse", "HEAD")
        return RepositoryFixture(
            temporary=temporary,
            root=root,
            predecessor=predecessor,
            config_raw=config_raw,
            result_raw=result_raw,
            authority_commit=authority_commit,
        )

    def authorize(
        self,
        fixture: RepositoryFixture,
        *,
        commit: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> authority.TDG8SuccessorExecutionAuthority:
        assert fixture.result_raw is not None
        assert fixture.authority_commit is not None
        with (
            patch.object(authority, "PREDECESSOR_REPAIR_COMMIT", fixture.predecessor),
            patch.object(
                authority,
                "observed_environment",
                return_value=(
                    dict(authority.PINNED_ENVIRONMENT)
                    if environment is None
                    else environment
                ),
            ),
        ):
            return authority.authorize_execution(
                fixture.root,
                fixture.config_raw,
                fixture.result_raw,
                fixture.authority_commit if commit is None else commit,
            )

    def test_exact_committed_image_returns_typed_one_event_authority(self) -> None:
        fixture = self.make_repository()
        receipt = self.authorize(fixture)
        self.assertIsInstance(receipt, authority.TDG8SuccessorExecutionAuthority)
        self.assertEqual(receipt.campaign_id, authority.CAMPAIGN_ID)
        self.assertEqual(receipt.branch, "GR-0")
        self.assertEqual(receipt.event_count, 1)
        self.assertTrue(receipt.source_terminal_read_only)
        self.assertTrue(receipt.candidate_branches_forbidden)
        self.assertEqual(
            {item.path for item in receipt.implementation_inventory},
            {path for _role, path in authority.REQUIRED_IMPLEMENTATION_INVENTORY},
        )
        self.assertTrue(
            {
                "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
                "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py",
                "src/recursive_horizons/fgc/evolution/hlt16_member_codec.py",
            }.issubset({item.path for item in receipt.implementation_inventory})
        )

    def test_authority_requires_the_exact_committed_successor_delta(self) -> None:
        extra = self.make_repository()
        write(extra.root / "adjacent-unreviewed-source.py", b"unexpected delta\n")
        git(extra.root, "add", "adjacent-unreviewed-source.py")
        git(extra.root, "commit", "-qm", "unexpected adjacent change")
        extra.authority_commit = git(extra.root, "rev-parse", "HEAD")
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "Git delta differs"
        ):
            self.authorize(extra)

        missing = self.make_repository(missing_delta_path="Makefile")
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "Git delta differs"
        ):
            self.authorize(missing)

    def test_destination_present_blocks_prelaunch_but_not_later_authorization(
        self,
    ) -> None:
        unfinished = self.make_repository(finalize=False)
        destination = unfinished.root / authority.DESTINATION_CAMPAIGN_PATH
        destination.mkdir(parents=True)
        with (
            patch.object(
                authority, "PREDECESSOR_REPAIR_COMMIT", unfinished.predecessor
            ),
            patch.object(
                authority,
                "observed_environment",
                return_value=dict(authority.PINNED_ENVIRONMENT),
            ),
        ):
            with self.assertRaisesRegex(
                authority.TDG8SuccessorAuthorityError, "destination.*present"
            ):
                authority.build_prelaunch(unfinished.config_raw, unfinished.root)

        finished = self.make_repository()
        (finished.root / authority.DESTINATION_CAMPAIGN_PATH).mkdir(parents=True)
        receipt = self.authorize(finished)
        self.assertEqual(
            receipt.destination_campaign_path,
            authority.DESTINATION_CAMPAIGN_PATH,
        )

    def test_wrong_head_and_unrelated_ancestry_fail_closed(self) -> None:
        fixture = self.make_repository()
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "current HEAD differs"
        ):
            self.authorize(fixture, commit=fixture.predecessor)

        tree = git(fixture.root, "write-tree")
        unrelated = git(
            fixture.root,
            "commit-tree",
            tree,
            input_raw=b"unrelated root\n",
        )
        git(fixture.root, "checkout", "-q", "--detach", unrelated)
        assert fixture.result_raw is not None
        with (
            patch.object(authority, "PREDECESSOR_REPAIR_COMMIT", fixture.predecessor),
            patch.object(
                authority,
                "observed_environment",
                return_value=dict(authority.PINNED_ENVIRONMENT),
            ),
        ):
            with self.assertRaisesRegex(
                authority.TDG8SuccessorAuthorityError, "not an ancestor"
            ):
                authority.authorize_execution(
                    fixture.root,
                    fixture.config_raw,
                    fixture.result_raw,
                    unrelated,
                )

    def test_tracked_drift_is_rejected(self) -> None:
        fixture = self.make_repository()
        runner = next(
            path
            for role, path in authority.REQUIRED_IMPLEMENTATION_INVENTORY
            if role == "runner"
        )
        write(fixture.root / runner, b"tracked drift\n")
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "tracked worktree changes"
        ):
            self.authorize(fixture)

    def test_committed_historical_binding_drift_is_rejected(self) -> None:
        fixture = self.make_repository()
        changed = authority.IMMUTABLE_BINDINGS[0]
        write(fixture.root / changed.path, b"self-consistent history rewrite\n")
        git(fixture.root, "add", changed.path)
        git(fixture.root, "commit", "-qm", "rewrite bound history")
        rewritten_head = git(fixture.root, "rev-parse", "HEAD")
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "immutable binding bytes differ"
        ):
            self.authorize(fixture, commit=rewritten_head)

    def test_environment_drift_is_rejected(self) -> None:
        fixture = self.make_repository()
        changed = dict(authority.PINNED_ENVIRONMENT)
        changed["numpy_version"] = "999.0"
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "environment differs"
        ):
            self.authorize(fixture, environment=changed)

    def test_committed_claim_promotion_is_rejected(self) -> None:
        fixture = self.make_repository()
        assert fixture.result_raw is not None
        promoted = json.loads(fixture.result_raw)
        promoted["artifact_payload"]["claims"]["candidate_execution_authorized"] = True
        fixture.result_raw = authority.canonical_result(promoted)
        write(fixture.root / authority.RESULT_PATH, fixture.result_raw)
        git(fixture.root, "add", authority.RESULT_PATH)
        git(fixture.root, "commit", "-qm", "promote forbidden claim")
        fixture.authority_commit = git(fixture.root, "rev-parse", "HEAD")
        with self.assertRaisesRegex(
            authority.TDG8SuccessorAuthorityError, "compact result differs"
        ):
            self.authorize(fixture)

    def test_noncanonical_or_semantic_compact_mutation_is_rejected(self) -> None:
        fixture = self.make_repository()
        assert fixture.result_raw is not None
        with patch.object(authority, "PREDECESSOR_REPAIR_COMMIT", fixture.predecessor):
            with self.assertRaisesRegex(
                authority.TDG8SuccessorAuthorityError, "noncanonical"
            ):
                authority.validate_compact(
                    fixture.config_raw, fixture.result_raw + b"\n"
                )
            changed = json.loads(fixture.result_raw)
            changed["source_config_sha256"] = "0" * 64
            with self.assertRaisesRegex(
                authority.TDG8SuccessorAuthorityError, "compact result differs"
            ):
                authority.validate_compact(
                    fixture.config_raw, authority.canonical_result(changed)
                )

    def test_unrelated_untracked_files_and_created_destination_are_allowed(
        self,
    ) -> None:
        fixture = self.make_repository()
        write(fixture.root / "unrelated.tmp", b"not authority\n")
        (fixture.root / authority.DESTINATION_CAMPAIGN_PATH).mkdir(parents=True)
        write(
            fixture.root / authority.DESTINATION_CAMPAIGN_PATH / "created.txt",
            b"post-authority destination\n",
        )
        receipt = self.authorize(fixture)
        self.assertEqual(receipt.event, 23)

    def test_compact_reproducer_is_store_blind_and_rejects_inner_symlink(self) -> None:
        fixture = self.make_repository()
        spec = importlib.util.spec_from_file_location(
            "test_tdg8_successor_reproducer", REPRODUCER_PATH
        )
        if spec is None or spec.loader is None:
            self.fail("cannot load TDG8 successor reproducer")
        reproducer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reproducer)
        (fixture.root / authority.DESTINATION_CAMPAIGN_PATH).mkdir(parents=True)
        with (
            patch.object(authority, "PREDECESSOR_REPAIR_COMMIT", fixture.predecessor),
            patch.object(
                authority,
                "_require_path_absent",
                side_effect=AssertionError("compact verification probed a store"),
            ),
        ):
            self.assertEqual(
                reproducer._verify_compact(fixture.root), fixture.result_raw
            )

        configs = fixture.root / "configs"
        configs.rename(fixture.root / "configs-real")
        configs.symlink_to(fixture.root / "configs-real", target_is_directory=True)
        with self.assertRaises(authority.TDG8SuccessorAuthorityError):
            reproducer._read_config(fixture.root)


if __name__ == "__main__":
    unittest.main()
