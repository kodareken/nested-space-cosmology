"""Focused prospective authority controls for TDG9 AC1."""

from __future__ import annotations

import inspect
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ti2 as ti2_runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ac1_authority as authority  # noqa: E402


class AC1AuthorityTests(unittest.TestCase):
    def test_config_freezes_retry3_all_channels_and_nonclaims(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        self.assertEqual(parsed["predecessor"]["commit"], authority.BASE_COMMIT)
        self.assertEqual(parsed["selection"]["retry"], 3)
        self.assertEqual(
            parsed["selection"]["channel_order"], list(authority.CHANNEL_ORDER)
        )
        self.assertNotIn(
            "historical_RK4_production_failed_channels",
            parsed["selection"],
        )
        self.assertEqual(
            parsed["selection"]["TI2_sampled_retry3_channels"],
            ["u:alpha", "u:R"],
        )
        self.assertEqual(parsed["work_budget"]["shadow_proposals"], 7)
        self.assertEqual(
            parsed["work_budget"]["maximum_stage_and_endpoint_RHS_records"], 28
        )
        self.assertFalse(parsed["work_budget"]["FAILED_OCCURRENCES_filter_applied"])
        self.assertFalse(parsed["work_budget"]["retry_4_authorized"])
        self.assertFalse(parsed["work_budget"]["LOC1_exact_localization_authorized"])
        self.assertTrue(parsed["decision"]["order_inconclusive_is_channel_failure"])
        self.assertFalse(parsed["claims"]["AC1_result_earned"])
        self.assertFalse(parsed["claims"]["production_method_earned"])

    def test_wrong_hash_retry_threshold_or_channel_order_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (authority.PREF1_RESULT_SHA256.encode(), b"0" * 64),
            (authority.TI2_RAW_TERMINAL_SHA256.encode(), b"1" * 64),
            (authority.TRANSACTION_SHA256.encode(), b"2" * 64),
            (b"retry = 3", b"retry = 4"),
            (b'minimum_observed_order = "3/2"', b'minimum_observed_order = "1"'),
            (b'"u:alpha"', b'"u:beta"'),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.AC1AuthorityError):
                    authority.parse_config(changed)

    def test_real_retry3_predecessor_fingerprint_precedes_any_compositor(self) -> None:
        with patch.object(
            ti2_runner.tdg6,
            "prepare_tdg6_gr0_compositor",
            side_effect=AssertionError("prelaunch constructed a shadow"),
        ) as compositor:
            receipts = authority.real_predecessor_fingerprints(ROOT)
        compositor.assert_not_called()
        self.assertEqual(len(receipts), 1)
        self.assertEqual(receipts[0]["retry"], 3)
        self.assertEqual(
            receipts[0]["transaction_sha256"], authority.TRANSACTION_SHA256
        )
        self.assertFalse(receipts[0]["shadow_proposal_constructed"])

    def test_compact_result_is_outcome_blind(self) -> None:
        config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
        result = authority.read_leaf(ROOT, authority.RESULT_PATH)
        parsed = authority.validate_compact(config, result)
        payload = parsed["artifact_payload"]
        self.assertFalse(payload["claims"]["AC1_result_earned"])
        self.assertFalse(payload["shadow_executed"])
        self.assertEqual(payload["shadow_proposals_constructed"], 0)
        self.assertNotIn("complete_admission_passed", payload)
        self.assertNotIn("failed_channels", payload)
        self.assertTrue(payload["replace_refs_absent"])

    def test_output_absence_is_prelaunch_only(self) -> None:
        authority.require_output_absent(ROOT)

    def test_store_identity_is_immutable(self) -> None:
        self.assertEqual(
            authority._snapshot_store(ROOT),
            (
                authority.SEALED_STORE_LEAF_COUNT,
                authority.SEALED_STORE_SNAPSHOT_SHA256,
            ),
        )

    def test_output_absence_rejects_existing_and_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / authority.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.AC1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.AC1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_delta_is_direct_successor_scoped_and_excludes_raw_qdrant(self) -> None:
        self.assertEqual(
            authority.BASE_COMMIT,
            "660f369365ab8417f4c68eded46acb9bb2b05a5b",
        )
        self.assertEqual(len(authority.DELTA_PATHS), 16)
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(any(path.startswith("runs/") for path in authority.DELTA_PATHS))
        self.assertNotIn(authority.OUTPUT_NAMESPACE, authority.DELTA_PATHS)
        self.assertNotIn("paper/recursive-horizons.md", authority.DELTA_PATHS)

    def test_authorize_does_not_require_output_absence(self) -> None:
        self.assertNotIn(
            "require_output_absent",
            inspect.getsource(authority.authorize),
        )
        self.assertIn(
            "require_output_absent",
            inspect.getsource(authority.build_prelaunch),
        )

    def test_raw_replay_receipt_schema_is_exact(self) -> None:
        receipt = authority.expected_raw_replay_receipt()
        self.assertEqual(
            set(receipt),
            {
                "accepted_time_hex",
                "attempted_width_hex",
                "descriptor_sha256",
                "historical_journal_sha256",
                "member_key",
                "retry",
                "state_sha256",
                "transaction_sha256",
            },
        )
        self.assertEqual(receipt["member_key"], "RK4-2049")
        self.assertEqual(receipt["retry"], 3)
        self.assertEqual(receipt["attempted_width_hex"], authority.WIDTH_HEX)
        self.assertEqual(receipt["accepted_time_hex"], authority.ACCEPTED_TIME_HEX)
        self.assertEqual(receipt["state_sha256"], authority.PHYSICAL_STATE_SHA256)
        self.assertEqual(
            receipt["descriptor_sha256"], authority.MEMBER_DESCRIPTOR_SHA256
        )
        self.assertEqual(receipt["transaction_sha256"], authority.TRANSACTION_SHA256)
        self.assertEqual(
            receipt["historical_journal_sha256"],
            "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59",
        )

    def test_compact_rejects_nan_and_infinity(self) -> None:
        config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token):
                with self.assertRaisesRegex(
                    authority.AC1AuthorityError,
                    "nonfinite AC1 compact JSON constant",
                ):
                    authority.validate_compact(
                        config, f'{{"gate_status": {token}}}\n'.encode()
                    )


if __name__ == "__main__":
    unittest.main()
