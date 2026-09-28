from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from typing import Callable

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto18_auth1_inputs import build_calibration_genesis


def runner_module():
    spec = importlib.util.spec_from_file_location(
        "hlt15_runner_test", ROOT / "scripts/run_fgc_proto17_hlt15_gen1.py",
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HLT15RunnerAuthorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = runner_module()
        record = json.loads((ROOT / "results/fgc-1-pro18-pref26.json").read_text())
        cls.evidence = record["artifact_payload"]["AUTH1_input_evidence"]

    def _git(self, root: Path, *args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=root, check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    def _authority_repository(
        self,
        root: Path,
        mutate: Callable[[dict[str, object]], None] | None = None,
        mutate_pref26: Callable[[dict[str, object]], None] | None = None,
    ) -> tuple[dict[str, str], Path]:
        self._git(root, "init", "-q")
        self._git(root, "config", "user.email", "test@example.invalid")
        self._git(root, "config", "user.name", "test")
        pref26_path = root / "results/fgc-1-pro18-pref26.json"
        pref26_path.parent.mkdir(parents=True)
        pref26_record = {
            "artifact_id": "FGC-1-PRO18-PREF26",
            "artifact_payload": {
                "test_fixture": True,
                "AUTH1_input_evidence": deepcopy(self.evidence),
            },
        }
        if mutate_pref26 is not None:
            mutate_pref26(pref26_record)
        pref26_path.write_bytes(self.runner._indented(pref26_record))
        self._git(root, "add", ".")
        self._git(root, "commit", "-qm", "pref26")
        pref26_commit = self._git(root, "rev-parse", "HEAD")

        authority_config = root / "configs/fgc/fgc-1-pro18-auth1.toml"
        raw_config = root / "configs/fgc/fgc-1-pro18-frz1.toml"
        authority_config.parent.mkdir(parents=True)
        authority_config.write_text("schema_version = 1\n", encoding="utf-8")
        raw_config.write_text("schema_version = 1\n", encoding="utf-8")
        (root / "unrelated-tracked.txt").write_text(
            "authenticated image\n", encoding="utf-8",
        )
        closure: list[dict[str, str]] = []
        direct_by_path = {path: role for role, path in self.runner._DIRECT_PATHS.items()}
        for relative in sorted(self.runner._EXPECTED_CLOSURE_PATHS):
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# authenticated test source\n", encoding="utf-8")
            closure.append({
                "path": relative,
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "role": direct_by_path.get(relative, "dependency"),
            })
        role_hash = {entry["role"]: entry["sha256"] for entry in closure if entry["role"] != "dependency"}
        source_hashes = {
            "protocol_config_sha256": "1" * 64,
            "protocol_freeze_result_sha256": "2" * 64,
            "hlt13_result_sha256": "3" * 64,
            "run_plan_sha256": "4" * 64,
            "runtime_module_sha256": role_hash["runtime_module"],
            "adapter_module_sha256": role_hash["adapter_module"],
            "runner_sha256": role_hash["runner"],
        }
        environment = {"fixture": True}
        genesis = build_calibration_genesis(
            self.evidence, campaign_id="temporary-authority",
            namespace="runs/fgc-2-sf1/proto17/calibration",
            source_hashes=source_hashes, numerical_environment=environment,
        ).genesis_spec
        payload: dict[str, object] = {
            "pref26_binding": {
                "historical_commit": pref26_commit,
                "compact_result": "results/fgc-1-pro18-pref26.json",
                "compact_result_sha256": sha256(pref26_path.read_bytes()).hexdigest(),
                "AUTH1_input_evidence_sha256": self.evidence["AUTH1_input_evidence_sha256"],
            },
            "authority_inputs": {
                "authority_config": {
                    "path": "configs/fgc/fgc-1-pro18-auth1.toml",
                    "sha256": sha256(authority_config.read_bytes()).hexdigest(),
                },
                "raw_source_config": {
                    "path": "configs/fgc/fgc-1-pro18-frz1.toml",
                    "sha256": sha256(raw_config.read_bytes()).hexdigest(),
                },
            },
            "auth1_input_evidence": deepcopy(self.evidence),
            "auth1_input_evidence_sha256": self.evidence["AUTH1_input_evidence_sha256"],
            "source_hashes": source_hashes,
            "campaign_id": "temporary-authority",
            "genesis_spec": genesis,
            "genesis_spec_sha256": self.runner._digest(genesis),
            "source_closure": closure,
            "source_closure_sha256": self.runner._digest(closure),
            "numerical_environment": {
                "contract": environment,
                "canonical_sha256": self.runner._digest(environment),
            },
            "claims": {
                **{name: True for name in self.runner._TRUE_CLAIMS},
                **{name: False for name in self.runner._FALSE_CLAIMS},
            },
        }
        if mutate is not None:
            mutate(payload)
        result = {"artifact_id": "FGC-1-PRO18-AUTH1", "artifact_payload": payload}
        result_path = root / "results/fgc-1-pro18-auth1.json"
        result_path.write_bytes(self.runner._indented(result))
        self._git(root, "add", ".")
        self._git(root, "commit", "-qm", "auth1")
        commit = self._git(root, "rev-parse", "HEAD")
        launch = {
            "authorization_commit": commit,
            "authorization_result_path": "results/fgc-1-pro18-auth1.json",
            "authorization_result_sha256": sha256(result_path.read_bytes()).hexdigest(),
            "genesis_spec_sha256": payload["genesis_spec_sha256"],
        }
        self._git(root, "checkout", "--detach", "-q", commit)
        return launch, result_path

    def test_real_auth1_uses_the_owner_indented_canonical_format(self):
        raw = (ROOT / "results/fgc-1-pro18-auth1.json").read_bytes()
        parsed = self.runner._decode(raw, label="AUTH1 result", style="indented")
        self.assertEqual(parsed["artifact_id"], "FGC-1-PRO18-AUTH1")
        self.runner._validate_environment(parsed["artifact_payload"])

    def test_temporary_git_authority_uses_commit_blob_and_creates_no_namespace(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _result_path = self._authority_repository(root)
            parsed = self.runner.verify_authority_tuple(
                root, launch, check_environment=False,
            )
            self.assertEqual(parsed["payload"]["campaign_id"], "temporary-authority")
            self.assertFalse((root / "runs").exists())

    def test_worktree_substitution_after_commit_fails(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, result_path = self._authority_repository(root)
            result_path.write_bytes(result_path.read_bytes() + b"\n")
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "live authority blob differs",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=False)

    def test_descendant_or_attached_checkout_is_not_a_launch_image(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _result_path = self._authority_repository(root)
            self._git(root, "switch", "-q", "-c", "attached-test")
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "not detached",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=False)
            (root / "later.txt").write_text("later\n", encoding="utf-8")
            self._git(root, "add", "later.txt")
            self._git(root, "commit", "-qm", "later")
            later = self._git(root, "rev-parse", "HEAD")
            self._git(root, "checkout", "--detach", "-q", later)
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "not the exact AUTH1 commit",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=False)

    def test_dirty_tracked_detached_checkout_is_not_a_launch_image(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _result_path = self._authority_repository(root)
            (root / "unrelated-tracked.txt").write_text("mutated\n", encoding="utf-8")
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "not tracked-clean",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=False)

    def test_pre_reimport_checks_do_not_enumerate_untracked_future_paths(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _result_path = self._authority_repository(root)
            guarded = root / "runs/fgc-2-sf1/proto17/calibration"
            guarded.mkdir(parents=True)
            (guarded / "must-not-be-enumerated").write_text("guard\n", encoding="utf-8")
            commands: list[tuple[str, ...]] = []
            actual_run = self.runner.subprocess.run

            def recording_run(*args, **kwargs):
                command = args[0] if args else kwargs["args"]
                commands.append(tuple(str(part) for part in command))
                return actual_run(*args, **kwargs)

            with patch.object(self.runner.subprocess, "run", side_effect=recording_run):
                parsed = self.runner.verify_authority_tuple(
                    root, launch, check_environment=False,
                )
            self.assertEqual(parsed["payload"]["campaign_id"], "temporary-authority")
            self.assertFalse(any(
                "status" in command or any("untracked" in part for part in command)
                for command in commands
            ))

    def test_pref26_embedded_evidence_must_equal_auth1_evidence(self):
        def corrupt_pref26(record):
            record["artifact_payload"]["AUTH1_input_evidence"]["source_manifest_sha256"] = "0" * 64

        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _result_path = self._authority_repository(
                root, mutate_pref26=corrupt_pref26,
            )
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "embedded AUTH1 evidence",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=False)

    def test_verified_bytes_are_executed_after_live_file_changes(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = "src/recursive_horizons/fgc/evolution/verified_fixture.py"
            path = root / relative
            path.parent.mkdir(parents=True)
            path.write_text("VALUE = 'mutated'\n", encoding="utf-8")
            code = f"""
import importlib.util
from pathlib import Path
runner_path = Path({str(ROOT / 'scripts/run_fgc_proto17_hlt15_gen1.py')!r})
spec = importlib.util.spec_from_file_location('verified_runner', runner_path)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
runner._MODULE_LOAD_ORDER = ((
    'recursive_horizons.fgc.evolution.verified_fixture', {relative!r},
),)
modules = runner._install_production_modules(
    Path({str(root)!r}), {{{relative!r}: b\"VALUE = 'verified'\\n\"}},
)
print(modules['recursive_horizons.fgc.evolution.verified_fixture'].VALUE)
"""
            observed = subprocess.run(
                [sys.executable, "-c", code], check=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            self.assertEqual(observed.stdout.strip(), "verified")

    def test_claim_promotion_and_role_swap_fail(self):
        def promote_dependency_to_runner(payload):
            dependency = next(
                item for item in payload["source_closure"]
                if item["role"] == "dependency"
            )
            dependency["role"] = "runner"

        for mutation, message in (
            (
                lambda payload: payload["claims"].__setitem__("HLT15_authorized", True),
                "claims",
            ),
            (
                promote_dependency_to_runner,
                "source closure",
            ),
        ):
            with self.subTest(message=message), TemporaryDirectory() as temporary:
                root = Path(temporary)
                launch, _path = self._authority_repository(root, mutation)
                with self.assertRaises(self.runner.Proto17HLT15LaunchError):
                    self.runner.verify_authority_tuple(root, launch, check_environment=False)

    def test_live_environment_mismatch_fails_before_any_namespace(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            launch, _path = self._authority_repository(root)
            with self.assertRaisesRegex(
                self.runner.Proto17HLT15LaunchError, "live numerical environment",
            ):
                self.runner.verify_authority_tuple(root, launch, check_environment=True)
            self.assertFalse((root / "runs").exists())

    def test_authority_tuple_requires_exact_fields(self):
        with self.assertRaises(self.runner.Proto17HLT15LaunchError):
            self.runner.verify_authority_tuple(ROOT, {})


if __name__ == "__main__":
    unittest.main()
