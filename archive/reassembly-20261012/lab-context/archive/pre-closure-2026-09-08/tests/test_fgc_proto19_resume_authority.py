from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
from pathlib import Path
import tomllib
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import proto19_resume_authority as sid3
from recursive_horizons.fgc.evolution.proto19_execution_closure import (
    ExecutionClosureRecord,
    FilePin,
    ImportPin,
    InterpreterIdentity,
    NamespacePin,
)


ROOT = Path(__file__).resolve().parents[1]


def _quote(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _array(values: list[str]) -> str:
    return "[" + ", ".join(_quote(value) for value in values) + "]"


class ResumeAuthorityFixture:
    def __init__(self) -> None:
        module_paths = dict(sid3.APPROVED_IMPORT_MODULE_PATHS)
        all_paths = {
            *module_paths.values(),
            "scripts/bootstrap_fgc_pro19_sid3.py",
            *sid3.MANDATORY_STATIC_INPUT_PATHS,
        }
        self.captured = {path: (ROOT / path).read_bytes() for path in all_paths}
        pins = []
        for path in sorted(all_paths):
            kind = "input" if path in sid3.MANDATORY_STATIC_INPUT_PATHS else "python"
            pins.append(
                FilePin(path, kind, "a" * 40, sha256(self.captured[path]).hexdigest())
            )
        by_path = {item.path: item for item in pins}
        imports = tuple(
            ImportPin(module, path, by_path[path].sha256)
            for module, path in sorted(module_paths.items())
        )
        self.record = ExecutionClosureRecord(
            implementation_commit="b" * 40,
            git_object_format="sha1",
            files=tuple(pins),
            interpreter=InterpreterIdentity(
                implementation="cpython",
                implementation_version="test",
                cache_tag="cpython-test",
                hexversion=1,
                executable="/python",
                executable_realpath="/python",
                executable_sha256="c" * 64,
                prefix="/prefix",
                base_prefix="/prefix",
                exec_prefix="/prefix",
                base_exec_prefix="/prefix",
                platform="test",
                platform_release="test",
                machine="test",
                byteorder="little",
                filesystem_encoding="utf-8",
                filesystem_errors="surrogateescape",
                sysconfig_platform="test",
                soabi="test",
                flags=(("dont_write_bytecode", True), ("isolated", True)),
                base_sys_path=("/stdlib",),
            ),
            imports=imports,
            namespaces=tuple(
                NamespacePin(module, path)
                for module, path in sorted(sid3.MANDATORY_NAMESPACE_PATHS.items())
            ),
            forbidden_module_prefixes=tuple(
                sorted(sid3.MANDATORY_FORBIDDEN_MODULE_PREFIXES)
            ),
            forbidden_script_prefixes=tuple(
                sorted(sid3.MANDATORY_FORBIDDEN_SCRIPT_PREFIXES)
            ),
            authority_delta_paths=tuple(
                sorted((sid3.CONFIG_PATH, sid3.CLOSURE_PATH, sid3.RESULT_PATH))
            ),
        )
        self.closure_raw = sid3.canonical_result(self.record.to_mapping())
        self.numpy = sid3.observe_numpy_runtime_identity().mapping()
        self.config_raw = self._config()
        self.result_raw = sid3.canonical_result(
            sid3.build_sid3_result(self.config_raw, self.closure_raw)
        )

    def _config(
        self, *, safe_to_resume: bool = False, numpy_version: str | None = None
    ) -> bytes:
        imports = [item.module for item in self.record.imports]
        namespaces = [item.module for item in self.record.namespaces]
        static = [item.path for item in self.record.files if item.kind == "input"]
        numpy = dict(self.numpy)
        if numpy_version is not None:
            numpy["version"] = numpy_version
        lines = [
            "schema_version = 1",
            f"artifact_id = {_quote(sid3.ARTIFACT_ID)}",
            'project_version = "0.11.0"',
            'target_protocol = "FGC-2-SF1-PROTO18"',
            f"classification = {_quote(sid3.CLASSIFICATION)}",
            f"nonclaims = {_array(sid3.EXPECTED_NONCLAIMS)}",
            "",
            "[predecessor]",
            'artifact_id = "FGC-1-PRO19-SID2-PREF1"',
            f"sealed_commit = {_quote(sid3.SID2_COMMIT)}",
            f"config_path = {_quote(sid3.SID2_CONFIG_PATH)}",
            f"config_sha256 = {_quote(sid3.SID2_CONFIG_SHA256)}",
            f"result_path = {_quote(sid3.SID2_RESULT_PATH)}",
            f"result_sha256 = {_quote(sid3.SID2_RESULT_SHA256)}",
            "",
            "[execution_closure]",
            f"path = {_quote(sid3.CLOSURE_PATH)}",
            f"raw_sha256 = {_quote(sha256(self.closure_raw).hexdigest())}",
            f"canonical_sha256 = {_quote(self.record.canonical_sha256)}",
            f"implementation_commit = {_quote(self.record.implementation_commit)}",
            f"authority_module = {_quote(sid3.AUTHORITY_MODULE)}",
            f"runner_module = {_quote(sid3.RUNNER_MODULE)}",
            f"required_import_modules = {_array(imports)}",
            f"required_namespaces = {_array(namespaces)}",
            "required_forbidden_module_prefixes = "
            + _array(list(self.record.forbidden_module_prefixes)),
            "required_forbidden_script_prefixes = "
            + _array(list(self.record.forbidden_script_prefixes)),
            f"authority_delta_paths = {_array(list(self.record.authority_delta_paths))}",
            f"static_input_paths = {_array(static)}",
            "",
            "[progression]",
            f"original_authorization_commit = {_quote(sid3.ORIGINAL_AUTHORIZATION_COMMIT)}",
            f"original_plan_sha256 = {_quote(sid3.ORIGINAL_PLAN_SHA256)}",
            f"campaign_id = {_quote(sid3.CAMPAIGN_ID)}",
            'branch = "GR-0"',
            'amplitude = "3"',
            "event = 23",
            'start_rational = "23/16"',
            'target_rational = "3/2"',
            "candidate_branch_opened = false",
            "",
            "[recovered_generation]",
            "generation = 8",
            f"checkpoint_sha256 = {_quote(sid3.GENERATION8_CHECKPOINT_SHA256)}",
            "journal_sequence = 8",
            f"journal_sha256 = {_quote(sid3.GENERATION8_JOURNAL_SHA256)}",
            'member_key = "RK4-2049"',
            f"member_accepted_time_hex = {_quote(sid3.GENERATION8_ACCEPTED_TIME_HEX)}",
            f"member_cursor_sha256 = {_quote(sid3.GENERATION8_CURSOR_SHA256)}",
            'member_mode = "RETRY_PENDING"',
            'pending_owner = "temporal"',
            f"pending_cap_hex = {_quote(sid3.GENERATION8_PENDING_CAP_HEX)}",
            'disposition = "nonterminal"',
            "candidate_branch_opened = false",
            "physical_state_advanced_by_recovery = false",
            "",
            "[numpy_runtime]",
            *[f"{key} = {_quote(value)}" for key, value in numpy.items()],
            "",
            "[scope]",
            "raw_store_mutation = false",
            "writer_lease_acquisition = false",
            "PDE_proposal_execution = false",
            "GEN0_reimport = false",
            "trajectory_resume = false",
            "candidate_branch_opened = false",
            "shared_candidate_capable_definitions_loaded = true",
            "read_only_runtime_preflight_required = true",
            "runtime_recheck_before_mutation_required = true",
            "",
            "[claims]",
            "committed_resume_contract_frozen = true",
            "committed_resume_authority_present = true",
            "committed_resume_image_authenticated_by_artifact_alone = false",
            "exact_generation8_inputs_bound = true",
            f"safe_to_resume_trajectory = {'true' if safe_to_resume else 'false'}",
            "trajectory_resume_authorized_by_artifact_alone = false",
            "GR0_calibration_completed = false",
            "candidate_execution_authorized = false",
            "candidate_runtime_configuration_state_or_outcome_opened = false",
            "physical_result_earned = false",
            "",
        ]
        return "\n".join(lines).encode("utf-8")


class Proto19ResumeAuthorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixture = ResumeAuthorityFixture()

    def test_build_validate_and_authorize_exact_committed_image(self) -> None:
        fixture = self.fixture
        self.assertEqual(
            sid3.validate_sid3_result(
                fixture.config_raw, fixture.closure_raw, fixture.result_raw
            )["gate_status"],
            "pass",
        )
        with patch.object(
            sid3, "construct_first_event", wraps=sid3.construct_first_event
        ) as construct:
            receipt = sid3.authorize_resume_image(
                fixture.config_raw,
                fixture.closure_raw,
                fixture.result_raw,
                captured_files=fixture.captured,
                authority_commit="d" * 40,
            )
        self.assertEqual(receipt.authority_commit, "d" * 40)
        self.assertEqual(receipt.implementation_commit, "b" * 40)
        self.assertEqual(receipt.store_path, sid3.STORE_PATH)
        self.assertEqual(receipt.checkpoint_generation, 8)
        self.assertEqual(receipt.checkpoint_sha256, sid3.GENERATION8_CHECKPOINT_SHA256)
        self.assertEqual(receipt.journal_sequence, 8)
        self.assertEqual(receipt.journal_tip_sha256, sid3.GENERATION8_JOURNAL_SHA256)
        self.assertEqual(receipt.progression_plan.sha256, sid3.ORIGINAL_PLAN_SHA256)
        self.assertEqual(
            set(receipt.static_input_bytes), set(sid3.MANDATORY_STATIC_INPUT_PATHS)
        )
        supplied = construct.call_args.kwargs["evidence_bytes"]
        self.assertEqual(set(supplied), set(sid3.PROGRESSION_EVIDENCE_PATHS))
        with self.assertRaises(TypeError):
            receipt.static_input_bytes["x"] = b"x"  # type: ignore[index]

    def test_captured_byte_mutation_or_inventory_change_is_rejected(self) -> None:
        fixture = self.fixture
        path = sid3.SID2_RESULT_PATH
        changed = dict(fixture.captured)
        changed[path] = changed[path][:-1] + b" "
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "captured closure byte"):
            sid3.authorize_resume_image(
                fixture.config_raw,
                fixture.closure_raw,
                fixture.result_raw,
                captured_files=changed,
                authority_commit="d" * 40,
            )
        missing = dict(fixture.captured)
        missing.pop(path)
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "inventory"):
            sid3.authorize_resume_image(
                fixture.config_raw,
                fixture.closure_raw,
                fixture.result_raw,
                captured_files=missing,
                authority_commit="d" * 40,
            )

    def test_artifact_cannot_promote_safe_to_resume(self) -> None:
        config = self.fixture._config(safe_to_resume=True)
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "exact contract"):
            sid3.build_sid3_result(config, self.fixture.closure_raw)

    def test_numpy_runtime_mismatch_is_rejected(self) -> None:
        config = self.fixture._config(numpy_version="0.0.invalid")
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "NumPy runtime"):
            sid3.build_sid3_result(config, self.fixture.closure_raw)

    def test_package_initializers_cannot_replace_exact_synthetic_namespaces(
        self,
    ) -> None:
        fixture = self.fixture
        bad = replace(
            fixture.record,
            namespaces=(NamespacePin("scripts", "scripts"),),
        )
        closure_raw = sid3.canonical_result(bad.to_mapping())
        contract = dict(
            tomllib.loads(fixture.config_raw.decode("utf-8"))["execution_closure"]
        )
        contract.update(
            raw_sha256=sha256(closure_raw).hexdigest(),
            canonical_sha256=bad.canonical_sha256,
            required_namespaces=[item.module for item in bad.namespaces],
        )
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "namespace set"):
            sid3._execution_contract(contract, closure_raw=closure_raw, closure=bad)

    def test_matching_config_cannot_widen_A_owned_import_permission(self) -> None:
        fixture = self.fixture
        extra_path = (
            "src/recursive_horizons/fgc/evolution/proto19_progression_freeze.py"
        )
        extra_module = "recursive_horizons.fgc.evolution.proto19_progression_freeze"
        extra_bytes = (ROOT / extra_path).read_bytes()
        extra_pin = FilePin(
            extra_path,
            "python",
            "a" * 40,
            sha256(extra_bytes).hexdigest(),
        )
        bad = replace(
            fixture.record,
            files=tuple(
                sorted((*fixture.record.files, extra_pin), key=lambda item: item.path)
            ),
            imports=tuple(
                sorted(
                    (
                        *fixture.record.imports,
                        ImportPin(extra_module, extra_path, extra_pin.sha256),
                    ),
                    key=lambda item: item.module,
                )
            ),
        )
        closure_raw = sid3.canonical_result(bad.to_mapping())
        contract = dict(
            tomllib.loads(fixture.config_raw.decode("utf-8"))["execution_closure"]
        )
        contract.update(
            raw_sha256=sha256(closure_raw).hexdigest(),
            canonical_sha256=bad.canonical_sha256,
            required_import_modules=[item.module for item in bad.imports],
        )
        with self.assertRaisesRegex(
            sid3.ResumeAuthorityError, "import permission surface differs"
        ):
            sid3._execution_contract(contract, closure_raw=closure_raw, closure=bad)

    def test_noncanonical_or_duplicate_result_is_rejected(self) -> None:
        compact = self.fixture.result_raw.replace(b"\n  ", b"\n ", 1)
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "noncanonical"):
            sid3.validate_sid3_result(
                self.fixture.config_raw, self.fixture.closure_raw, compact
            )
        duplicate = self.fixture.result_raw.replace(
            b'{\n  "artifact_id"',
            b'{\n  "artifact_id": "duplicate",\n  "artifact_id"',
            1,
        )
        with self.assertRaisesRegex(sid3.ResumeAuthorityError, "malformed"):
            sid3.validate_sid3_result(
                self.fixture.config_raw, self.fixture.closure_raw, duplicate
            )


if __name__ == "__main__":
    unittest.main()
