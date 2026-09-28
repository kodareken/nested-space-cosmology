"""Focused prospective authority controls for TDG9 UR1."""

from __future__ import annotations

import inspect
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import tdg9_ur1_authority as authority  # noqa: E402


class UR1AuthorityTests(unittest.TestCase):
    def test_config_binds_sealed_predecessors_question_and_nonclaims(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        self.assertEqual(parsed["predecessor"]["commit"], authority.BASE_COMMIT)
        self.assertEqual(
            parsed["predecessor"]["AC1_PREF1_result_sha256"],
            "173183871b2d750f98e9758cdc5f702be40b0bbb73d10e707672e0e68752935d",
        )
        self.assertEqual(
            parsed["predecessor"]["AC1_raw_manifest_sha256"],
            authority.AC1_RAW_MANIFEST_SHA256,
        )
        self.assertEqual(
            parsed["predecessor"]["AC1_raw_terminal_sha256"],
            authority.AC1_RAW_TERMINAL_SHA256,
        )
        self.assertEqual(parsed["selection"]["published_channel"], "u:R")
        self.assertEqual(parsed["selection"]["owned_row_count"], 2044)
        self.assertEqual(
            parsed["selection"]["TI2_retry3_u_R_row_sha256"],
            "b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969",
        )
        self.assertEqual(parsed["work_budget"]["shadow_proposals"], 7)
        self.assertEqual(
            parsed["work_budget"]["maximum_stage_and_endpoint_RHS_records"], 28
        )
        self.assertFalse(parsed["claims"]["UR1_result_earned"])
        self.assertFalse(parsed["claims"]["state_advance_authorized"])
        self.assertFalse(parsed["scope"]["temporal_retry_admission_authorized"])
        self.assertFalse(parsed["scope"]["fine_path_commit_authorized"])

    def test_wrong_hash_owner_formula_or_replay_receipt_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (authority.AC1_PREF1_RESULT_SHA256.encode(), b"0" * 64),
            (authority.AC1_RAW_TERMINAL_SHA256.encode(), b"1" * 64),
            (authority.TI2_RETRY3_U_R_ROW_SHA256.encode(), b"2" * 64),
            (b'retry = 3', b'retry = 4'),
            (
                b'bernstein_certification_slack_part_of_lower_clip = false',
                b'bernstein_certification_slack_part_of_lower_clip = true',
            ),
            (
                b'"raw_candidate_maximum_is_zero"',
                b'"unbounded_owner"',
            ),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.UR1AuthorityError):
                    authority.parse_config(changed)

    def test_exact_integrated_direct_successor_scope(self) -> None:
        self.assertEqual(authority.BASE_COMMIT, "22f7dc5e181af45eb43f6f19dee23fbc29f89c0c")
        self.assertEqual(len(authority.DELTA_PATHS), 16)
        self.assertEqual(
            set(authority.DELTA_PATHS),
            {
                "Makefile",
                "README.md",
                "src/recursive_horizons/fgc/evolution/tdg9_ur1_authority.py",
                "scripts/run_fgc_tdg9_ur1.py",
                "scripts/reproduce_fgc_tdg9_ur1_frz1.py",
                "scripts/check_repo.py",
                "configs/fgc/fgc-1-tdg9-ur1-frz1.toml",
                "results/fgc-1-tdg9-ur1-frz1.json",
                "results/README.md",
                "docs/claim-ledger.md",
                "docs/fgc-runtime-matrix.md",
                "docs/fgc-tdg9-ur1-frz1.md",
                "docs/research-roadmap.md",
                "tests/test_check_repo_tdg9_ur1_frz1.py",
                "tests/test_fgc_tdg9_ur1_authority.py",
                "tests/test_fgc_tdg9_ur1_runner.py",
            },
        )
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(any(path.startswith("runs/") for path in authority.DELTA_PATHS))
        self.assertTrue({"README.md", "Makefile"}.issubset(authority.DELTA_PATHS))

    def test_compact_is_canonical_and_outcome_blind(self) -> None:
        config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
        result = authority.read_leaf(ROOT, authority.RESULT_PATH)
        parsed = authority.validate_compact(config, result)
        payload = parsed["artifact_payload"]
        self.assertFalse(payload["UR1_result_earned"])
        self.assertFalse(payload["shadow_executed"])
        self.assertEqual(payload["shadow_proposals_constructed"], 0)
        self.assertNotIn("u_R", payload)
        self.assertNotIn("zero_lower_owner", payload)
        self.assertTrue(payload["UR1_output_namespace_absent_at_prelaunch"])
        self.assertEqual(
            payload["AC1_raw_identity_observed_at_prelaunch"]["leaf_count"], 2
        )

    def test_prelaunch_restores_fingerprint_without_proposal_construction(self) -> None:
        config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
        expected_raw = {
            "namespace": authority.AC1_RAW_NAMESPACE,
            "leaf_count": 2,
            "manifest_sha256": authority.AC1_RAW_MANIFEST_SHA256,
            "terminal_sha256": authority.AC1_RAW_TERMINAL_SHA256,
        }

        def fake_git(_root: Path, *arguments: str) -> bytes:
            if arguments[:2] == ("rev-parse", "--verify"):
                return (authority.BASE_COMMIT + "\n").encode()
            if arguments[:2] == ("for-each-ref", "--format=%(refname)"):
                return b""
            raise AssertionError(arguments)

        with (
            patch.object(authority, "_git", side_effect=fake_git),
            patch.object(
                authority,
                "_working",
                return_value=((), (), authority.DELTA_PATHS),
            ),
            patch.object(authority, "_verify_predecessor_bindings"),
            patch.object(authority, "_verify_live_bindings"),
            patch.object(authority, "_verify_ac1_raw", return_value=expected_raw),
            patch.object(authority, "require_output_absent") as absent,
            patch.object(
                authority,
                "real_predecessor_fingerprints",
                return_value=(authority._expected_fingerprint_receipt(),),
            ) as fingerprints,
            patch.object(authority.ac1.ti2.ti1, "_environment", return_value=authority.ENVIRONMENT),
        ):
            result = authority.build_prelaunch(
                config,
                ROOT,
                store_snapshot=(
                    authority.SEALED_STORE_LEAF_COUNT,
                    authority.SEALED_STORE_SNAPSHOT_SHA256,
                ),
            )
        absent.assert_called_once_with(ROOT)
        fingerprints.assert_called_once_with(ROOT)
        self.assertEqual(result["artifact_payload"]["shadow_proposals_constructed"], 0)
        self.assertFalse(result["artifact_payload"]["shadow_executed"])

    def test_authorize_allows_status_after_raw_but_prelaunch_requires_absence(self) -> None:
        self.assertNotIn("require_output_absent", inspect.getsource(authority.authorize))
        self.assertIn("require_output_absent", inspect.getsource(authority.build_prelaunch))

    def test_output_absence_rejects_existing_and_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / authority.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.UR1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.UR1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_compact_rejects_nonfinite_json_constants(self) -> None:
        config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token), self.assertRaisesRegex(
                authority.UR1AuthorityError,
                "nonfinite UR1 compact JSON constant",
            ):
                authority.validate_compact(config, f'{{"gate_status": {token}}}\n'.encode())


if __name__ == "__main__":
    unittest.main()
