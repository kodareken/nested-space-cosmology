"""Focused prospective authority controls for TDG10-QA1."""

from __future__ import annotations

from hashlib import sha256
import inspect
from pathlib import Path
import sys
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import tdg10_qa1_authority as authority  # noqa: E402


class QA1AuthorityTests(unittest.TestCase):
    def test_config_freezes_retry3_generation9_and_nonclaims(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        predecessor = parsed["predecessor"]
        self.assertEqual(parsed["project_version"], "0.11.0")
        self.assertEqual(parsed["target_protocol"], "FGC-2-SF1-PROTO18")
        self.assertEqual(parsed["base_commit"], authority.BASE_COMMIT)
        self.assertEqual(
            predecessor["immediate_predecessor_artifact_id"],
            "FGC-1-TDG9-UR1-PREF1",
        )
        self.assertEqual(predecessor["immediate_predecessor_commit"], authority.BASE_COMMIT)
        self.assertEqual(
            predecessor["immediate_predecessor_result_path"],
            "results/fgc-1-tdg9-ur1-pref1.json",
        )
        self.assertEqual(
            predecessor["immediate_predecessor_result_sha256"],
            "27c2aef9b033092285b868f583441842de99f61941ff57aa7e62749cc8d91300",
        )
        self.assertEqual(predecessor["generation"], 9)
        self.assertEqual(predecessor["retry"], 3)
        self.assertEqual(predecessor["prior_retry_count"], 2)
        self.assertEqual(predecessor["checkpoint_journal_sequence"], 10)
        self.assertEqual(predecessor["historical_retry3_rejection_sequence"], 11)
        self.assertEqual(
            predecessor["checkpoint_sha256"],
            "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56",
        )
        self.assertEqual(
            predecessor["historical_retry3_rejection_sha256"],
            "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a",
        )
        self.assertEqual(predecessor["cursor_mode"], "RETRY_PENDING")
        self.assertEqual(predecessor["pending_owner"], "temporal")
        self.assertEqual(parsed["formulation"]["tableau_selector"], "SSPRK3")
        self.assertEqual(
            parsed["formulation"]["actual_spatial_operator"],
            "inherited_RK4_2049_SBP4",
        )
        self.assertEqual(parsed["formulation"]["channel_count"], 18)
        self.assertTrue(parsed["formulation"]["complete_admission_is_all_of"])
        self.assertFalse(parsed["transaction"]["campaign_state_write_authorized"])
        self.assertFalse(parsed["transaction"]["PDE_state_commit_authorized"])
        self.assertFalse(parsed["transaction"]["fine_path_commit_authorized"])
        self.assertFalse(parsed["transaction"]["diagnostic_fine_endpoint_is_accepted_state"])
        self.assertTrue(
            parsed["transaction"]["diagnostic_fine_endpoint_serialized_iff_all_channels_pass"]
        )
        self.assertFalse(any(parsed["exclusions"].values()))
        self.assertEqual(int(authority.REPLAY["predecessor_generation"]), 9)
        self.assertNotEqual(int(authority.REPLAY["predecessor_generation"]), 10)
        self.assertEqual(int(authority.REPLAY["journal_sequence"]), 11)

    def test_wrong_hash_generation_retry_or_flag_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (authority.CHECKPOINT_SHA256.encode(), b"0" * 64),
            (authority.CHECKPOINT_JOURNAL_SHA256.encode(), b"1" * 64),
            (authority.HISTORICAL_RETRY3_REJECTION_SHA256.encode(), b"2" * 64),
            (authority.TRANSACTION_SHA256.encode(), b"3" * 64),
            (b"generation = 9", b"generation = 10"),
            (b"checkpoint_journal_sequence = 10", b"checkpoint_journal_sequence = 11"),
            (b"retry = 3", b"retry = 4"),
            (b'target_protocol = "FGC-2-SF1-PROTO18"', b'target_protocol = "FGC-2-SF1-PROTO19"'),
            (b'project_version = "0.11.0"', b'project_version = "0.12.0"'),
            (
                authority.UR1_PREF1_RESULT_SHA256.encode(),
                b"4" * 64,
            ),
            (b"campaign_state_write_authorized = false", b"campaign_state_write_authorized = true"),
            (b"physical_result_authorized = false", b"physical_result_authorized = true"),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.QA1AuthorityError):
                    authority.parse_config(changed)

    def test_compact_is_canonical_outcome_blind_and_fail_closed(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        expected = authority.expected_compact(config)
        raw = authority.canonical_pretty(expected)
        parsed = authority.validate_compact(config, raw)
        payload = parsed["artifact_payload"]
        self.assertEqual(parsed["gate_status"], "pass")
        self.assertFalse(payload["shadow_executed"])
        self.assertEqual(payload["shadow_proposals_constructed"], 0)
        self.assertTrue(payload["diagnostic_qualification_only"])
        self.assertFalse(payload["state_advance_authorized"])
        self.assertFalse(payload["new_protocol_id_authorized"])
        self.assertEqual(payload["target_protocol"], "FGC-2-SF1-PROTO18")
        self.assertEqual(payload["retry3_restores_generation"], 9)
        self.assertTrue(payload["retry3_never_restores_generation_10"])
        self.assertTrue(payload["checkpoint_journal_is_sequence_10"])
        self.assertTrue(payload["historical_retry3_rejection_is_sequence_11"])
        self.assertNotIn("complete_admission_passed", payload)
        self.assertNotIn("failed_channels", payload)
        self.assertNotIn("diagnostic_payload_sha256", payload)
        self.assertEqual(
            payload["implementation"],
            authority.parse_config(config)["implementation"],
        )
        self.assertEqual(payload["project_version"], "0.11.0")
        self.assertEqual(raw, authority.canonical_pretty(parsed))

    def test_compact_rejects_nonfinite_json_constants(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token), self.assertRaisesRegex(
                authority.QA1AuthorityError,
                "nonfinite QA1 compact JSON constant",
            ):
                authority.validate_compact(
                    config, f'{{"gate_status": {token}}}\n'.encode()
                )

    def test_output_absence_rejects_existing_and_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / authority.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.QA1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{authority.STAGING_PREFIX}dead").mkdir()
            with self.assertRaisesRegex(
                authority.QA1AuthorityError, "staging namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.QA1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_generation9_checkpoint_tip_is_journal_sequence_10(self) -> None:
        authority._verify_predecessor_leaves(ROOT)
        self.assertEqual(authority.PREDECESSOR_GENERATION, 9)
        self.assertEqual(authority.FORBIDDEN_RETRY3_GENERATION, 10)
        self.assertEqual(authority.CHECKPOINT_JOURNAL_SEQUENCE, 10)
        self.assertEqual(authority.HISTORICAL_RETRY3_REJECTION_SEQUENCE, 11)
        self.assertNotEqual(
            authority.CHECKPOINT_JOURNAL_SHA256,
            authority.HISTORICAL_RETRY3_REJECTION_SHA256,
        )

    def test_prelaunch_restores_fingerprint_without_proposal_construction(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()

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
            patch.object(authority, "read_leaf", return_value=b"leaf"),
            patch.object(authority, "verify_implementation_inventory", return_value=()),
            patch.object(authority, "_verify_ur1_pref1_result"),
            patch.object(authority, "_verify_predecessor_leaves"),
            patch.object(authority, "require_output_absent") as absent,
            patch.object(
                authority,
                "real_predecessor_fingerprints",
                return_value=(authority._expected_fingerprint_receipt(),),
            ) as fingerprints,
            patch.object(authority, "_environment", return_value=authority.ENVIRONMENT),
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
        self.assertEqual(result["artifact_payload"]["retry3_restores_generation"], 9)

    def test_authorize_requires_direct_successor_and_matching_bytes(self) -> None:
        self.assertNotIn(
            "require_output_absent", inspect.getsource(authority.authorize)
        )
        self.assertIn(
            "require_output_absent", inspect.getsource(authority.build_prelaunch)
        )
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        result = authority.canonical_pretty(authority.expected_compact(config))
        commit = "a" * 40

        def fake_git(_root: Path, *arguments: str) -> bytes:
            if arguments[:2] == ("rev-parse", "--verify"):
                return (commit + "\n").encode()
            if arguments[0] == "rev-list":
                return f"{commit} {authority.BASE_COMMIT}\n".encode()
            if arguments[:2] == ("for-each-ref", "--format=%(refname)"):
                return b"refs/replace/dead\n"
            raise AssertionError(arguments)

        with (
            patch.object(authority, "_git", side_effect=fake_git),
            patch.object(authority, "_environment", return_value=authority.ENVIRONMENT),
            self.assertRaisesRegex(authority.QA1AuthorityError, "replace refs"),
        ):
            authority.authorize(ROOT, config, result, commit)

    def test_implementation_block_is_required_and_bound(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        pairs = authority.verify_implementation_inventory(
            ROOT, parsed, commit=None
        )
        self.assertEqual(len(pairs), 5)
        self.assertEqual(
            [path for path, _digest in pairs],
            [path for _key, path in authority.IMPLEMENTATION_INVENTORY],
        )
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        digest = parsed["implementation"]["runner_sha256"]
        malformed = raw.replace(digest.encode(), b"not-a-canonical-sha256-digest", 1)
        self.assertNotEqual(malformed, raw)
        with self.assertRaises(authority.QA1AuthorityError):
            authority.parse_config(malformed)
        missing = raw.split(b"\n[implementation]\n")[0] + b"\n"
        with self.assertRaises(authority.QA1AuthorityError):
            authority.parse_config(missing)
        mutated = dict(parsed)
        mutated["implementation"] = dict(parsed["implementation"])
        mutated["implementation"]["runner_sha256"] = "0" * 64
        with self.assertRaisesRegex(
            authority.QA1AuthorityError, "implementation binding differs"
        ):
            authority.verify_implementation_inventory(
                ROOT, mutated, commit=None
            )

    def test_read_leaf_rejects_parent_symlink_and_leaf_race(self) -> None:
        payload = b"qa1-leaf\n"
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            real = Path(foreign) / "leaf.txt"
            real.write_bytes(payload)
            (root / "parent").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.QA1AuthorityError, "authority parent is unsafe"
            ):
                authority.read_leaf(root, "parent/leaf.txt")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "configs" / "fgc"
            target.mkdir(parents=True)
            leaf = target / "bound.toml"
            leaf.write_bytes(payload)
            leaf.unlink()
            leaf.symlink_to("/etc/passwd")
            with self.assertRaisesRegex(
                authority.QA1AuthorityError, "unsafe authority leaf"
            ):
                authority.read_leaf(root, "configs/fgc/bound.toml")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = "configs/fgc/bound.toml"
            leaf = root / relative
            leaf.parent.mkdir(parents=True)
            leaf.write_bytes(payload)
            calls = {"regular": 0}
            real_fstat = authority.os.fstat

            def racing_fstat(fd: int) -> object:
                result = real_fstat(fd)
                if authority.stat.S_ISREG(result.st_mode):
                    calls["regular"] += 1
                    if calls["regular"] == 2:
                        return SimpleNamespace(
                            st_dev=result.st_dev,
                            st_ino=result.st_ino,
                            st_size=result.st_size,
                            st_mtime_ns=result.st_mtime_ns + 1,
                            st_nlink=result.st_nlink,
                            st_mode=result.st_mode,
                        )
                return result

            with (
                patch.object(authority.os, "fstat", side_effect=racing_fstat),
                self.assertRaisesRegex(
                    authority.QA1AuthorityError, "authority leaf changed"
                ),
            ):
                authority.read_leaf(root, relative)

    def test_ur1_pref1_live_and_committed_bytes_match(self) -> None:
        authority._verify_ur1_pref1_result(ROOT)
        live = authority.read_leaf(ROOT, authority.UR1_PREF1_RESULT_PATH)
        self.assertEqual(sha256(live).hexdigest(), authority.UR1_PREF1_RESULT_SHA256)

    def test_delta_excludes_raw_store_and_qdrant(self) -> None:
        self.assertEqual(authority.BASE_COMMIT, "49c514ac283ae3ec8091a6638442da6fdaf58db0")
        self.assertEqual(len(authority.DELTA_PATHS), 20)
        self.assertEqual(
            set(authority.DELTA_PATHS),
            {
                "Makefile",
                "README.md",
                authority.CONFIG_PATH,
                "docs/claim-ledger.md",
                "docs/fgc-runtime-matrix.md",
                authority.OWNER_DOCUMENT,
                "docs/research-roadmap.md",
                authority.RESULT_PATH,
                "results/README.md",
                "scripts/check_repo.py",
                "scripts/reproduce_fgc_tdg10_qa1_frz1.py",
                authority.RUNNER_PATH,
                authority.AUTHORITY_PATH,
                authority.EXACT_ADMISSION_PATH,
                authority.EXACT_RUNTIME_PATH,
                "tests/test_check_repo_tdg10_qa1_frz1.py",
                "tests/test_fgc_tdg10_exact_complete_c_admission.py",
                "tests/test_fgc_tdg10_exact_complete_c_runtime.py",
                "tests/test_fgc_tdg10_qa1_authority.py",
                "tests/test_fgc_tdg10_qa1_runner.py",
            },
        )
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(any(path.startswith("runs/") for path in authority.DELTA_PATHS))
        self.assertNotIn(authority.OUTPUT_NAMESPACE, authority.DELTA_PATHS)

    def test_receipt_flags_are_diagnostic_only(self) -> None:
        receipt = authority.QA1Authority(authority_commit="b" * 40)
        self.assertTrue(receipt.diagnostic_qualification_only)
        self.assertFalse(receipt.state_advance_authorized)
        self.assertFalse(receipt.campaign_state_write_authorized)
        self.assertFalse(receipt.PDE_state_commit_authorized)
        self.assertFalse(receipt.fine_path_commit_authorized)
        self.assertFalse(receipt.candidate_execution_authorized)
        self.assertFalse(receipt.physical_result_authorized)
        self.assertFalse(receipt.new_protocol_id_authorized)
        self.assertEqual(receipt.target_protocol, "FGC-2-SF1-PROTO18")
        self.assertEqual(receipt.predecessor_generation, 9)
        self.assertEqual(receipt.checkpoint_journal_sequence, 10)
        self.assertEqual(receipt.historical_retry3_rejection_sequence, 11)


if __name__ == "__main__":
    unittest.main()
