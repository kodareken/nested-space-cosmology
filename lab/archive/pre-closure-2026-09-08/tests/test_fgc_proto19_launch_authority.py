"""External-Git and no-follow tests for the PRO19 first-event authority."""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import proto19_launch_authority as authority


ROLES = ("runtime", "config", "result", "documentation", "reproducer", "runner", "test")


def git(root: Path, *args: str) -> str:
    result = subprocess.run(("git", "-C", str(root), *args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return result.stdout.decode("utf-8").strip()


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def rebind_manifest(root: Path, manifest: str, path: str, old: str, new: str) -> None:
    target = root / manifest
    text = target.read_text(encoding="utf-8")
    needle = f'path = "{path}"\nsha256 = "{old}"'
    if needle not in text:
        raise AssertionError(f"missing manifest binding for {path}")
    target.write_text(text.replace(needle, f'path = "{path}"\nsha256 = "{new}"', 1), encoding="utf-8")


def sealed_mon16_result(
    *,
    source_config_sha256: str,
    predecessors: list[tuple[str, str]],
    inventory: list[tuple[str, str, str]],
) -> str:
    return json.dumps({
        "artifact_id": "FGC-1-HLT16-MON16",
        "classification": "sealed_durable_runtime_qualified_nonexecuting",
        "artifact_payload": {
            "claims": {
                "durable_runtime_implemented": True,
                "durable_runtime_qualified": True,
                "one_event_preflight_authorized": False,
                "state_advance_authorized": False,
                "GR0_calibration_completed": False,
                "candidate_execution_authorized": False,
                "physical_transition_claim_authorized": False,
            },
            "explicit_nonexecution": {
                "arrays_opened": False,
                "live_run_store_opened": False,
                "runner_invoked": False,
                "one_event_preflight_completed": False,
            },
            "coordinator_refresh_required_paths": [],
            "predecessors": [
                {"artifact_id": "fixture", "binding_status": "sealed", "path": path,
                 "sha256": digest, "observed_sha256": digest}
                for path, digest in predecessors
            ],
            "inventory": [
                {"role": role, "binding_status": "sealed", "path": path,
                 "sha256": digest, "observed_sha256": digest}
                for role, path, digest in inventory
            ],
        },
        "source_config_sha256": source_config_sha256,
    }, sort_keys=True) + "\n"


class Proto19LaunchAuthorityTests(unittest.TestCase):
    def make_repository(self) -> tuple[tempfile.TemporaryDirectory[str], Path, str, str, str]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        git(root, "init", "-q"); git(root, "config", "user.email", "test@example.invalid"); git(root, "config", "user.name", "Test")
        write(root / "freeze.txt", "freeze\n")
        git(root, "add", "freeze.txt"); git(root, "commit", "-qm", "freeze")
        freeze = git(root, "rev-parse", "HEAD")
        observed = authority.observed_environment()
        rows: list[str] = [f'schema = "{authority.MANIFEST_SCHEMA}"', "", "[environment]", 'source = "HLT16"', ""]
        for index, role in enumerate(ROLES):
            path = f"authority/{role}-{index}.txt"; payload = f"{role}\n".encode()
            write(root / path, payload.decode())
            rows.extend(("[[authority_path]]", f'role = "{role}"', f'path = "{path}"', f'sha256 = "{sha256(payload).hexdigest()}"', ""))
        predecessor_path = "authority/mon16-predecessor.txt"
        inventory_path = "authority/mon16-inventory.txt"
        write(root / predecessor_path, "predecessor\n")
        write(root / inventory_path, "inventory\n")
        predecessor_digest = sha256((root / predecessor_path).read_bytes()).hexdigest()
        inventory_digest = sha256((root / inventory_path).read_bytes()).hexdigest()
        mon16_config = "configs/fgc/fgc-1-hlt16-mon16.toml"
        config_text = "\n".join((
            'artifact_id = "FGC-1-HLT16-MON16"',
            "",
            "[[predecessors]]",
            'artifact_id = "fixture"',
            f'path = "{predecessor_path}"',
            f'sha256 = "{predecessor_digest}"',
            "",
            "[[inventory]]",
            'role = "source"',
            f'path = "{inventory_path}"',
            f'sha256 = "{inventory_digest}"',
            "",
        ))
        write(root / mon16_config, config_text)
        config_digest = sha256((root / mon16_config).read_bytes()).hexdigest()
        mon16_result = "results/fgc-1-hlt16-mon16.json"
        write(root / mon16_result, sealed_mon16_result(
            source_config_sha256=config_digest,
            predecessors=[(predecessor_path, predecessor_digest)],
            inventory=[("source", inventory_path, inventory_digest)],
        ))
        for path in (
            "results/fgc-1-pro19-frz1.json",
            "results/fgc-1-pro18-auth1.json",
        ):
            payload = b"{}\n"
            write(root / path, payload.decode())
            rows.extend(("[[authority_path]]", 'role = "result"', f'path = "{path}"', f'sha256 = "{sha256(payload).hexdigest()}"', ""))
        for role, path in (
            ("config", mon16_config),
            ("result", mon16_result),
            ("config", predecessor_path),
            ("runtime", inventory_path),
        ):
            payload = (root / path).read_bytes()
            rows.extend(("[[authority_path]]", f'role = "{role}"', f'path = "{path}"', f'sha256 = "{sha256(payload).hexdigest()}"', ""))
        rows.extend(("[environment.contract]", *(f'{key} = "{value}"' for key, value in observed.items()), ""))
        manifest = "configs/fgc/fgc-1-pro19-launch-authority.toml"
        write(root / manifest, "\n".join(rows))
        git(root, "add", "."); git(root, "commit", "-qm", "authority")
        authorization = git(root, "rev-parse", "HEAD")
        return directory, root, authorization, freeze, manifest

    def test_authorizes_only_exact_clean_descendant_image(self) -> None:
        directory, root, commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        sentinel = object()
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), \
             patch.object(authority, "construct_first_event", return_value=sentinel) as construct:
            receipt = authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)
        self.assertEqual(receipt.authorization_commit, commit)
        self.assertTrue(receipt.freeze_is_ancestor_of_authorization)
        self.assertEqual(receipt.progression_plan, sentinel)
        self.assertEqual({role for role, _path, _digest in receipt.authority_paths}, set(ROLES))
        supplied = construct.call_args.kwargs["evidence_bytes"]
        self.assertEqual(
            set(supplied),
            {
                "results/fgc-1-pro19-frz1.json",
                "results/fgc-1-pro18-auth1.json",
            },
        )
        self.assertIn(
            "results/fgc-1-hlt16-mon16.json",
            {path for _role, path, _digest in receipt.authority_paths},
        )

    def test_pending_or_claim_promoted_mon16_result_refuses_launch(self) -> None:
        directory, root, commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        result = root / "results/fgc-1-hlt16-mon16.json"
        original_sha = sha256(result.read_bytes()).hexdigest()
        pending = json.loads(result.read_text())
        pending["classification"] = "implemented_qualified_runtime_inventory_pending_integration_refresh"
        pending["artifact_payload"]["coordinator_refresh_required_paths"] = ["src/runtime.py"]
        result.write_text(json.dumps(pending, sort_keys=True) + "\n")
        git(root, "add", str(result.relative_to(root))); git(root, "commit", "-qm", "pending mon16")
        pending_commit = git(root, "rev-parse", "HEAD")
        # The old manifest's hash catches drift first; rebuild it with the
        # changed result hash to prove the semantic authority catches pending.
        manifest_data = (root / manifest).read_text()
        old_sha = original_sha
        new_sha = sha256(result.read_bytes()).hexdigest()
        (root / manifest).write_text(manifest_data.replace(old_sha, new_sha))
        git(root, "add", manifest); git(root, "commit", "-qm", "bind pending mon16")
        pending_commit = git(root, "rev-parse", "HEAD")
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), patch.object(authority, "construct_first_event", return_value=object()):
            with self.assertRaisesRegex(authority.Proto19LaunchAuthorityError, "MON16 result is not sealed"):
                authority.authorize_first_event(root, authorization_commit=pending_commit, manifest_path=manifest)

    def test_semantically_promoted_mon16_claim_refuses_launch(self) -> None:
        directory, root, _commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        result = root / "results/fgc-1-hlt16-mon16.json"
        original_sha = sha256(result.read_bytes()).hexdigest()
        promoted = json.loads(result.read_text())
        promoted["artifact_payload"]["claims"]["state_advance_authorized"] = True
        result.write_text(json.dumps(promoted, sort_keys=True) + "\n")
        old_sha = original_sha
        new_sha = sha256(result.read_bytes()).hexdigest()
        manifest_data = (root / manifest).read_text()
        (root / manifest).write_text(manifest_data.replace(old_sha, new_sha))
        git(root, "add", str(result.relative_to(root)), manifest)
        git(root, "commit", "-qm", "bind promoted mon16")
        promoted_commit = git(root, "rev-parse", "HEAD")
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), patch.object(authority, "construct_first_event", return_value=object()):
            with self.assertRaisesRegex(authority.Proto19LaunchAuthorityError, "MON16 durable-runtime claim boundary differs"):
                authority.authorize_first_event(root, authorization_commit=promoted_commit, manifest_path=manifest)

    def test_self_consistent_forged_mon16_config_and_result_cannot_hide_actual_bytes(self) -> None:
        directory, root, _commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        config_path = "configs/fgc/fgc-1-hlt16-mon16.toml"
        result_path = "results/fgc-1-hlt16-mon16.json"
        predecessor_path = "authority/mon16-predecessor.txt"
        old_config = (root / config_path).read_bytes()
        old_result = (root / result_path).read_bytes()
        forged = "0" * 64
        config = (root / config_path).read_text(encoding="utf-8")
        original_pred = sha256((root / predecessor_path).read_bytes()).hexdigest()
        (root / config_path).write_text(config.replace(original_pred, forged, 1), encoding="utf-8")
        forged_config_sha = sha256((root / config_path).read_bytes()).hexdigest()
        result = json.loads((root / result_path).read_text())
        result["source_config_sha256"] = forged_config_sha
        result["artifact_payload"]["predecessors"][0]["sha256"] = forged
        result["artifact_payload"]["predecessors"][0]["observed_sha256"] = forged
        (root / result_path).write_text(json.dumps(result, sort_keys=True) + "\n", encoding="utf-8")
        rebind_manifest(root, manifest, config_path, sha256(old_config).hexdigest(), forged_config_sha)
        rebind_manifest(root, manifest, result_path, sha256(old_result).hexdigest(), sha256((root / result_path).read_bytes()).hexdigest())
        git(root, "add", config_path, result_path, manifest); git(root, "commit", "-qm", "forge self-consistent mon16")
        commit = git(root, "rev-parse", "HEAD")
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), patch.object(authority, "construct_first_event", return_value=object()):
            with self.assertRaisesRegex(authority.Proto19LaunchAuthorityError, "actual-byte"):
                authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)

    def test_stale_mon16_source_config_hash_refuses_even_when_result_is_manifest_bound(self) -> None:
        directory, root, _commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        result_path = "results/fgc-1-hlt16-mon16.json"
        old_result = (root / result_path).read_bytes()
        stale = json.loads(old_result)
        stale["source_config_sha256"] = "0" * 64
        (root / result_path).write_text(json.dumps(stale, sort_keys=True) + "\n", encoding="utf-8")
        rebind_manifest(root, manifest, result_path, sha256(old_result).hexdigest(), sha256((root / result_path).read_bytes()).hexdigest())
        git(root, "add", result_path, manifest); git(root, "commit", "-qm", "stale mon16 config hash")
        commit = git(root, "rev-parse", "HEAD")
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), patch.object(authority, "construct_first_event", return_value=object()):
            with self.assertRaisesRegex(authority.Proto19LaunchAuthorityError, "source-config"):
                authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)

    def test_untracked_is_ignored_but_head_and_tracked_drift_stop(self) -> None:
        directory, root, commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", freeze), patch.object(authority, "construct_first_event", return_value=object()):
            write(root / "untracked.tmp", "permitted\n")
            authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)
            write(root / "authority/runtime-0.txt", "drift\n")
            with self.assertRaises(authority.Proto19LaunchAuthorityError):
                authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)
            git(root, "checkout", "-q", "--", "authority/runtime-0.txt")
            git(root, "checkout", "-q", freeze)
            with self.assertRaises(authority.Proto19LaunchAuthorityError):
                authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)

    def test_manifest_closure_environment_and_ancestry_fail_closed(self) -> None:
        directory, root, commit, freeze, manifest = self.make_repository()
        self.addCleanup(directory.cleanup)
        with self.assertRaises(authority.Proto19LaunchAuthorityError):
            authority._manifest_paths({"authority_path": []})
        # An unrelated/future-looking freeze identity cannot authorize the image.
        with patch.object(authority, "SEALED_PRO19_FREEZE_COMMIT", "f" * 40), patch.object(authority, "construct_first_event", return_value=object()):
            with self.assertRaises(authority.Proto19LaunchAuthorityError):
                authority.authorize_first_event(root, authorization_commit=commit, manifest_path=manifest)

    def test_no_follow_reader_rejects_leaf_and_inner_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); write(root / "safe.txt", "safe\n")
            os.symlink(root / "safe.txt", root / "leaf")
            os.mkdir(root / "real"); write(root / "real/value", "ok\n")
            os.symlink(root / "real", root / "inner")
            with self.assertRaises(authority.Proto19LaunchAuthorityError):
                authority._nofollow_regular_bytes(root, "leaf")
            with self.assertRaises(authority.Proto19LaunchAuthorityError):
                authority._nofollow_regular_bytes(root, "inner/value")


if __name__ == "__main__":
    unittest.main()
