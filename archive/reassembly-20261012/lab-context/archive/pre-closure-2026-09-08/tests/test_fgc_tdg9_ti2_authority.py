"""Focused prospective authority controls for TDG9 TI2."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ti2 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ti2_authority as authority  # noqa: E402


class TI2AuthorityTests(unittest.TestCase):
    def test_config_freezes_recovery_encoding_and_nonclaims(self) -> None:
        parsed = authority.parse_config(
            (ROOT / authority.CONFIG_PATH).read_bytes()
        )
        self.assertEqual(
            parsed["recovery"]["TI1_authority_commit"], authority.BASE_COMMIT
        )
        self.assertEqual(
            parsed["fingerprint_encoding"]["schema"],
            "finite_builtin_binary64_hex_v1",
        )
        self.assertTrue(
            parsed["recovery"]["failure_occurs_before_prepare_shadow"]
        )
        self.assertEqual(
            parsed["semantics"]["actual_spatial_operator"],
            "inherited_RK4_2049_SBP4",
        )
        self.assertFalse(parsed["semantics"]["production_SSPRK3_comparator"])
        self.assertFalse(parsed["scope"]["PDE_state_commit_authorized"])
        self.assertFalse(parsed["claims"]["TI2_result_earned"])

    def test_wrong_recovery_hash_float_schema_or_threshold_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (authority.TI1_INVALID_STDOUT_SHA256.encode(), b"0" * 64),
            (b'first_live_float_path = "monitor.last_accepted_time"',
             b'first_live_float_path = "monitor.other"'),
            (b'schema = "finite_builtin_binary64_hex_v1"', b'schema = "decimal"'),
            (authority.TRANSACTION_SHA256_BY_RETRY[3].encode(), b"1" * 64),
            (b'minimum_observed_order = "3/2"', b'minimum_observed_order = "1"'),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.TI2AuthorityError):
                    authority.parse_config(changed)

    def test_real_three_predecessor_fingerprints_precede_any_compositor(self) -> None:
        with patch.object(
            runner.tdg6,
            "prepare_tdg6_gr0_compositor",
            side_effect=AssertionError("prelaunch constructed a shadow"),
        ) as compositor:
            receipts = authority.real_predecessor_fingerprints(ROOT)
        compositor.assert_not_called()
        self.assertEqual(
            [item["transaction_sha256"] for item in receipts],
            [authority.TRANSACTION_SHA256_BY_RETRY[index] for index in (3, 4, 5)],
        )
        self.assertTrue(
            all(not item["shadow_proposal_constructed"] for item in receipts)
        )

    def test_historical_ti1_blobs_are_authenticated_from_base_commit(self) -> None:
        authority._verify_historical_bindings(ROOT)
        proof = authority._verify_historical_source_order(ROOT)
        self.assertEqual(proof["shadow_proposals_before_serialization_failure"], 0)
        self.assertLess(proof["run_restore_line"], proof["run_prepare_shadow_line"])

    def test_old_and_new_output_absence_is_prelaunch_only(self) -> None:
        authority.ti1.require_output_absent(ROOT)
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
                authority.TI2AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.TI2AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_delta_is_direct_successor_scoped_and_excludes_raw_qdrant(self) -> None:
        self.assertEqual(authority.BASE_COMMIT, "8a32636654dbd6479fcb91a88d4557a8090041d3")
        self.assertEqual(len(authority.DELTA_PATHS), 18)
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(any(path.startswith("runs/") for path in authority.DELTA_PATHS))
        self.assertNotIn(authority.OUTPUT_NAMESPACE, authority.DELTA_PATHS)


if __name__ == "__main__":
    unittest.main()
