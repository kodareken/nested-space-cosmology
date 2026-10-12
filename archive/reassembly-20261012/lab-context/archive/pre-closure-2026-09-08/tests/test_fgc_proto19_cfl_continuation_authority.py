"""Adversarial tests for the one-use PRO19 CFL continuation authority."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import (
    proto19_cfl_continuation_authority as authority,
)
from recursive_horizons.fgc.evolution.proto19_launch_authority import (
    Proto19LaunchAuthorityError,
)


ROOT = Path(__file__).resolve().parents[1]


class CFLContinuationAuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config_raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        self.anchor = authority.expected_store_anchor()
        self.result_raw = authority._canonical(
            authority.build_cfl1_result(self.config_raw, self.anchor)
        )
        self.old_pro19 = b'{"artifact_id":"FGC-1-PRO19-FRZ1"}\n'
        self.old_auth1 = b'{"artifact_id":"FGC-1-PRO18-AUTH1"}\n'
        self.bytes = {
            authority.CONFIG_PATH: self.config_raw,
            authority.RESULT_PATH: self.result_raw,
            "results/fgc-1-pro19-frz1.json": self.old_pro19,
            "results/fgc-1-pro18-auth1.json": self.old_auth1,
        }
        self.paths = tuple(
            ("config" if path.startswith("configs/") else "result", path, authority._sha(raw))
            for path, raw in sorted(self.bytes.items())
        )
        self.current_plan = SimpleNamespace(branch="GR-0")
        self.original_plan = SimpleNamespace(
            sha256=authority.ORIGINAL_PLAN,
            campaign_id=authority.ORIGINAL_CAMPAIGN,
            branch="GR-0",
            amplitude="3",
        )
        self.launch = SimpleNamespace(
            progression_plan=self.current_plan,
            manifest_sha256="1" * 64,
            authority_paths=self.paths,
            environment={"python_version": "test"},
        )

    def authorize(self, *, anchor=None):
        with patch.object(
            authority, "authorize_first_event", return_value=self.launch
        ), patch.object(
            authority,
            "_nofollow_regular_bytes",
            side_effect=lambda _root, path: self.bytes[path],
        ), patch.object(
            authority, "construct_first_event", return_value=self.original_plan
        ) as construct:
            receipt = authority.authorize_cfl_continuation(
                ROOT,
                continuation_commit="2" * 40,
                launch_manifest_path="manifest.toml",
                store_anchor=self.anchor if anchor is None else anchor,
            )
        self.assertEqual(
            set(construct.call_args.kwargs["evidence_bytes"]),
            {
                "results/fgc-1-pro19-frz1.json",
                "results/fgc-1-pro18-auth1.json",
            },
        )
        self.assertEqual(
            construct.call_args.kwargs["authorization_commit"],
            authority.ORIGINAL_COMMIT,
        )
        return receipt

    def test_exact_anchor_preserves_original_plan_under_corrected_image(self) -> None:
        receipt = self.authorize()
        self.assertEqual(receipt.continuation_commit, "2" * 40)
        self.assertEqual(receipt.original_authorization_commit, authority.ORIGINAL_COMMIT)
        self.assertEqual(receipt.original_plan_sha256, authority.ORIGINAL_PLAN)
        self.assertEqual(receipt.recovery_checkpoint_generation, authority.GENERATION)
        self.assertEqual(receipt.recovery_checkpoint_sha256, authority.CHECKPOINT_SHA256)
        self.assertEqual(receipt.recovery_journal_sequence, authority.JOURNAL_SEQUENCE)
        self.assertEqual(receipt.recovery_member_key, authority.MEMBER_KEY)
        self.assertEqual(
            receipt.recovery_member_descriptor_sha256,
            authority.MEMBER_DESCRIPTOR_SHA256,
        )
        self.assertIs(receipt.progression_plan, self.original_plan)

    def test_every_store_anchor_mutation_is_rejected(self) -> None:
        mutations = {
            "authorization_commit": "3" * 40,
            "plan_sha256": "4" * 64,
            "campaign_id": "other",
            "checkpoint_generation": 5,
            "checkpoint_sha256": "5" * 64,
            "journal_sequence": 4,
            "journal_tip_sha256": "6" * 64,
            "member_descriptor_sha256": "7" * 64,
            "member_accepted_time_hex": (23 / 16).hex(),
            "member_mode": "RETRY_PENDING",
            "cfl_retry_total": 1,
            "durable_cfl_rejection_exists": True,
            "durable_terminal_exists": True,
            "suffix_classification": "journal_suffix",
            "writer_state": "foreign_live_lock",
        }
        for key, value in mutations.items():
            with self.subTest(key=key):
                changed = dict(self.anchor)
                changed[key] = value
                with self.assertRaisesRegex(
                    authority.CFLContinuationAuthorityError,
                    "store cross-binding|store anchor",
                ):
                    self.authorize(anchor=changed)

    def test_result_or_manifest_closure_drift_is_rejected(self) -> None:
        changed = dict(self.bytes)
        changed[authority.RESULT_PATH] = self.result_raw + b"\n"
        with patch.object(
            authority, "authorize_first_event", return_value=self.launch
        ), patch.object(
            authority,
            "_nofollow_regular_bytes",
            side_effect=lambda _root, path: changed[path],
        ):
            with self.assertRaisesRegex(
                authority.CFLContinuationAuthorityError,
                "authenticated image",
            ):
                authority.authorize_cfl_continuation(
                    ROOT,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )

        omitted = SimpleNamespace(
            progression_plan=self.current_plan,
            manifest_sha256="1" * 64,
            authority_paths=tuple(
                item for item in self.paths if item[1] != authority.RESULT_PATH
            ),
            environment={},
        )
        with patch.object(
            authority, "authorize_first_event", return_value=omitted
        ):
            with self.assertRaisesRegex(
                authority.CFLContinuationAuthorityError,
                "outside the closed launch manifest",
            ):
                authority.authorize_cfl_continuation(
                    ROOT,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )

    def test_current_image_and_original_plan_each_fail_closed(self) -> None:
        with patch.object(
            authority,
            "authorize_first_event",
            side_effect=Proto19LaunchAuthorityError("drift"),
        ):
            with self.assertRaisesRegex(
                authority.CFLContinuationAuthorityError,
                "corrected execution image",
            ):
                authority.authorize_cfl_continuation(
                    ROOT,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )

        bad_plan = SimpleNamespace(
            sha256="8" * 64,
            campaign_id=authority.ORIGINAL_CAMPAIGN,
            branch="GR-0",
            amplitude="3",
        )
        with patch.object(
            authority, "authorize_first_event", return_value=self.launch
        ), patch.object(
            authority,
            "_nofollow_regular_bytes",
            side_effect=lambda _root, path: self.bytes[path],
        ), patch.object(
            authority, "construct_first_event", return_value=bad_plan
        ):
            with self.assertRaisesRegex(
                authority.CFLContinuationAuthorityError,
                "original progression plan",
            ):
                authority.authorize_cfl_continuation(
                    ROOT,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )


if __name__ == "__main__":
    unittest.main()
