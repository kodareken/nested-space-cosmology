"""Focused tests for the prospective TDG8 retry recovery authority."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import (
    tdg8_retry_recovery_authority as authority,
)


ROOT = Path(__file__).resolve().parents[1]


def _git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ("git", "-C", str(root), *arguments),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.stdout.strip()


class TDG8RetryRecoveryAuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config_raw = (ROOT / authority.CONFIG_PATH).read_bytes()

    def _compact_result(self) -> bytes:
        config = authority._parse_config(self.config_raw)
        inventory = [
            {"role": role, "path": path, "sha256": f"{index + 1:x}" * 64}
            for index, (role, path) in enumerate(authority.IMPLEMENTATION_INVENTORY)
        ]
        # Keep every digest exactly 64 lowercase hexadecimal characters.
        for item in inventory:
            item["sha256"] = item["sha256"][:64]
        environment = {
            "implementation": "CPython",
            "python_version": "3.14.0",
            "numpy_version": "2.3.0",
            "system": "Darwin",
            "machine": "arm64",
            "platform": "darwin",
            "byteorder": "little",
        }
        fixed = authority._fixed_live_contract()
        result = {
            "schema_version": authority.SCHEMA_VERSION,
            "artifact_id": authority.ARTIFACT_ID,
            "project_version": authority.PROJECT_VERSION,
            "target_protocol": authority.TARGET_PROTOCOL,
            "classification": authority.CLASSIFICATION,
            "gate_status": "pass",
            "source_config_sha256": authority._sha(self.config_raw),
            "artifact_payload": authority._payload(
                config, inventory, fixed, environment
            ),
        }
        return authority.canonical(result)

    def test_config_and_compact_result_close_only_the_declared_scope(self) -> None:
        config = authority._parse_config(self.config_raw)
        self.assertEqual(config["claims"], authority.CLAIMS)
        self.assertTrue(config["scope"]["candidate_branches_forbidden"])
        self.assertFalse(config["scope"]["generation_zero_reimport"])
        result = authority.validate_compact(
            self.config_raw, self._compact_result()
        )
        self.assertEqual(result["gate_status"], "pass")
        self.assertFalse(
            result["artifact_payload"]["claims"]["physical_result_earned"]
        )

    def test_diagnosis_changes_only_the_omitted_frozen_minimum_width(self) -> None:
        record = {
            "payload": {
                "evidence": {
                    "initial_time": float.fromhex("0x1.78554de5a30e0p+0"),
                    "retry_step_size": float.fromhex(
                        authority.EXPECTED_PENDING_CAP_HEX
                    ),
                }
            }
        }
        diagnosis = authority._diagnose(record, "0x1.8000000000000p+0")
        self.assertEqual(
            diagnosis["differing_fields"],
            {
                "minimum_width_hex_or_none": {
                    "buggy": None,
                    "required": authority.FROZEN_MINIMUM_WIDTH_HEX,
                }
            },
        )
        self.assertTrue(diagnosis["numerical_boundaries_equal"])
        self.assertFalse(diagnosis["threshold_or_rule_change"])

    def test_compact_mutations_cannot_promote_claims_or_change_the_anchor(self) -> None:
        original = json.loads(self._compact_result().decode("ascii"))
        mutations = []

        claim = deepcopy(original)
        claim["artifact_payload"]["claims"]["physical_result_earned"] = True
        mutations.append(claim)

        width = deepcopy(original)
        width["artifact_payload"]["diagnosis"]["differing_fields"][
            "minimum_width_hex_or_none"
        ]["required"] = "0x1.0000000000000p-29"
        mutations.append(width)

        anchor = deepcopy(original)
        anchor["artifact_payload"]["live_anchor"]["checkpoint_generation"] = 8
        mutations.append(anchor)

        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with self.assertRaises(authority.TDG8RetryRecoveryAuthorityError):
                    authority.validate_compact(
                        self.config_raw, authority.canonical(mutation)
                    )

    def test_git_image_requires_the_exact_prospective_delta(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        _git(root, "init", "-q")
        _git(root, "config", "user.email", "rcv1@example.invalid")
        _git(root, "config", "user.name", "RCV1 test")
        (root / "base.txt").write_text("base\n", encoding="utf-8")
        _git(root, "add", "base.txt")
        _git(root, "commit", "-qm", "base")
        predecessor = _git(root, "rev-parse", "HEAD")

        for name in ("one.txt", "two.txt"):
            (root / name).write_text(f"{name}\n", encoding="utf-8")
        _git(root, "add", "one.txt", "two.txt")
        _git(root, "commit", "-qm", "exact authority delta")
        exact = _git(root, "rev-parse", "HEAD")

        with (
            patch.object(authority, "ORIGINAL_EXECUTION_COMMIT", predecessor),
            patch.object(
                authority, "AUTHORITY_DELTA_PATHS", ("one.txt", "two.txt")
            ),
        ):
            authority._require_git_image(root, exact)

            (root / "extra.txt").write_text("extra\n", encoding="utf-8")
            _git(root, "add", "extra.txt")
            _git(root, "commit", "-qm", "extra path")
            extra = _git(root, "rev-parse", "HEAD")
            with self.assertRaisesRegex(
                authority.TDG8RetryRecoveryAuthorityError,
                "authority Git delta differs",
            ):
                authority._require_git_image(root, extra)

        missing = tempfile.TemporaryDirectory()
        self.addCleanup(missing.cleanup)
        missing_root = Path(missing.name)
        _git(missing_root, "init", "-q")
        _git(missing_root, "config", "user.email", "rcv1@example.invalid")
        _git(missing_root, "config", "user.name", "RCV1 test")
        (missing_root / "base.txt").write_text("base\n", encoding="utf-8")
        _git(missing_root, "add", "base.txt")
        _git(missing_root, "commit", "-qm", "base")
        missing_predecessor = _git(missing_root, "rev-parse", "HEAD")
        (missing_root / "one.txt").write_text("one\n", encoding="utf-8")
        _git(missing_root, "add", "one.txt")
        _git(missing_root, "commit", "-qm", "missing path")
        missing_commit = _git(missing_root, "rev-parse", "HEAD")
        with (
            patch.object(
                authority, "ORIGINAL_EXECUTION_COMMIT", missing_predecessor
            ),
            patch.object(
                authority, "AUTHORITY_DELTA_PATHS", ("one.txt", "two.txt")
            ),
        ):
            with self.assertRaisesRegex(
                authority.TDG8RetryRecoveryAuthorityError,
                "authority Git delta differs",
            ):
                authority._require_git_image(missing_root, missing_commit)

    @staticmethod
    def _member(*, retry_pending: bool) -> SimpleNamespace:
        return SimpleNamespace(
            descriptor_sha256=authority.ANCHOR_DESCRIPTOR_SHA256,
            cursor={
                "accepted_boundary_time": dict(authority.ANCHOR_ACCEPTED_TIME),
                "mode": "RETRY_PENDING" if retry_pending else "FRESH_READY",
            },
            pending_owner="temporal" if retry_pending else "cfl",
            pending_cap_hex=(
                authority.EXPECTED_PENDING_CAP_HEX if retry_pending else "0x1p-6"
            ),
            cfl_current=1,
            cfl_total=1,
            ledger={"cumulative_temporal_retry_count": 1 if retry_pending else 0},
        )

    @staticmethod
    def _status() -> SimpleNamespace:
        return SimpleNamespace(
            active_write=True,
            writer_state="stale_verified_lock",
            safe_to_restart=False,
        )

    def test_intermediate_sequence_eight_cut_derives_only_frozen_generation_eight(
        self,
    ) -> None:
        member = self._member(retry_pending=False)
        checkpoint = SimpleNamespace(
            authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
            plan_sha256=authority.ORIGINAL_PLAN_SHA256,
            campaign_id=authority.CAMPAIGN_ID,
            protocol=authority.TARGET_PROTOCOL,
            event=23,
            target={"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"},
            generation=authority.ANCHOR_CHECKPOINT_GENERATION,
            sha256=authority.ANCHOR_CHECKPOINT_SHA256,
            journal_sequence=6,
            journal_tip_sha256=authority.ANCHOR_CHECKPOINT_JOURNAL_TIP_SHA256,
            members={authority.MEMBER_KEY: member},
        )
        rejection = {
            "sequence": authority.ANCHOR_JOURNAL_SEQUENCE,
            "record_sha256": authority.ANCHOR_JOURNAL_SHA256,
            "kind": "tdg6_rejection",
            "payload": {"member_key": authority.MEMBER_KEY},
        }
        transition = {
            "sequence": authority.EXPECTED_JOURNAL_SEQUENCE,
            "record_sha256": authority.EXPECTED_JOURNAL_SHA256,
            "kind": "cursor_transition",
            "previous_record_sha256": authority.ANCHOR_JOURNAL_SHA256,
            "payload": {"member_key": authority.MEMBER_KEY},
        }
        recovered = SimpleNamespace(
            generation=authority.EXPECTED_CHECKPOINT_GENERATION,
            sha256=authority.EXPECTED_CHECKPOINT_SHA256,
            journal_sequence=authority.EXPECTED_JOURNAL_SEQUENCE,
            journal_tip_sha256=authority.EXPECTED_JOURNAL_SHA256,
        )
        preview = SimpleNamespace(
            records=(rejection, transition),
            predecessor_checkpoint_sha256=authority.ANCHOR_CHECKPOINT_SHA256,
            evolution_state_sha256=authority.ANCHOR_PHYSICAL_STATE_SHA256,
            physical_state_advanced=False,
            checkpoint=recovered,
        )
        snapshot = SimpleNamespace(
            checkpoint=checkpoint,
            suffix_kinds=("tdg6_rejection", "cursor_transition"),
            suffix_records=(rejection, transition),
            suffix_classification="tdg6_rejection+cursor_transition",
            staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            terminal_lock_present=False,
        )
        store = SimpleNamespace(
            authenticated_snapshot=lambda: snapshot,
            inspect_recovery=self._status,
            load_state=lambda _descriptor: SimpleNamespace(arrays={}),
        )
        with (
            patch.object(authority, "HLT16CampaignStore", return_value=store),
            patch.object(
                authority.campaign_recovery,
                "_resolve_tdg6_rejection",
                return_value=preview,
            ),
        ):
            result = authority.inspect_recovery_boundary(ROOT)
        self.assertEqual(result["state"], "generation7_transition_published")
        self.assertEqual(
            result["physical_state_sha256"],
            authority.ANCHOR_PHYSICAL_STATE_SHA256,
        )

    def test_exact_generation_eight_stale_writer_is_metadata_cleanup_only(self) -> None:
        member = self._member(retry_pending=True)
        checkpoint = SimpleNamespace(
            authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
            plan_sha256=authority.ORIGINAL_PLAN_SHA256,
            campaign_id=authority.CAMPAIGN_ID,
            protocol=authority.TARGET_PROTOCOL,
            event=23,
            target={"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"},
            generation=authority.EXPECTED_CHECKPOINT_GENERATION,
            sha256=authority.EXPECTED_CHECKPOINT_SHA256,
            parent_sha256=authority.ANCHOR_CHECKPOINT_SHA256,
            journal_sequence=authority.EXPECTED_JOURNAL_SEQUENCE,
            journal_tip_sha256=authority.EXPECTED_JOURNAL_SHA256,
            disposition="nonterminal",
            members={authority.MEMBER_KEY: member},
        )
        snapshot = SimpleNamespace(
            checkpoint=checkpoint,
            suffix_kinds=(),
            suffix_records=(),
            suffix_classification="clean",
            staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            terminal_lock_present=False,
        )
        store = SimpleNamespace(
            authenticated_snapshot=lambda: snapshot,
            inspect_recovery=self._status,
            load_state=lambda _descriptor: SimpleNamespace(
                arrays={"u": object(), "p": object(), "q": object()}
            ),
        )
        with (
            patch.object(authority, "HLT16CampaignStore", return_value=store),
            patch.object(
                authority,
                "array_content_sha256",
                return_value=authority.ANCHOR_PHYSICAL_STATE_SHA256,
            ),
        ):
            result = authority.inspect_recovery_boundary(ROOT)
        self.assertEqual(result["state"], "generation8_checkpoint_published")
        self.assertEqual(
            result["physical_state_sha256"],
            authority.ANCHOR_PHYSICAL_STATE_SHA256,
        )


if __name__ == "__main__":
    unittest.main()
