"""Focused prospective authority controls for TDG10-QA2."""

from __future__ import annotations

from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import tdg10_qa2_authority as authority  # noqa: E402


class QA2AuthorityTests(unittest.TestCase):
    def test_config_freezes_retry4_and_retry5_and_nonclaims(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        predecessor = parsed["predecessor"]
        retry_4 = predecessor["retry_4"]
        retry_5 = predecessor["retry_5"]
        self.assertEqual(parsed["project_version"], "0.11.0")
        self.assertEqual(parsed["target_protocol"], "FGC-2-SF1-PROTO18")
        self.assertEqual(parsed["base_commit"], authority.BASE_COMMIT)
        self.assertEqual(
            predecessor["immediate_predecessor_artifact_id"],
            "FGC-1-TDG10-QA1-PREF1",
        )
        self.assertEqual(predecessor["immediate_predecessor_commit"], authority.BASE_COMMIT)
        self.assertEqual(
            predecessor["immediate_predecessor_result_path"],
            "results/fgc-1-tdg10-qa1-pref1.json",
        )
        self.assertEqual(
            predecessor["immediate_predecessor_result_sha256"],
            "3a107bcede479df4326aac5e194beefa1042860a8be498d7dc5d0519fdbc8603",
        )
        self.assertTrue(predecessor["immediate_predecessor_licenses_qa2_method_design"])
        self.assertFalse(predecessor["immediate_predecessor_licenses_old_member_adoption"])
        self.assertEqual(predecessor["selected_retries"], [4, 5])
        self.assertEqual(retry_4["generation"], 10)
        self.assertEqual(retry_4["retry"], 4)
        self.assertEqual(retry_4["prior_retry_count"], 3)
        self.assertEqual(retry_4["checkpoint_journal_sequence"], 12)
        self.assertEqual(retry_4["historical_rejection_sequence"], 13)
        self.assertEqual(retry_4["forbidden_predecessor_generation"], 11)
        self.assertEqual(
            retry_4["checkpoint_sha256"],
            "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187",
        )
        self.assertEqual(retry_5["generation"], 11)
        self.assertEqual(retry_5["retry"], 5)
        self.assertEqual(retry_5["prior_retry_count"], 4)
        self.assertEqual(retry_5["checkpoint_journal_sequence"], 14)
        self.assertEqual(retry_5["historical_rejection_sequence"], 15)
        self.assertEqual(retry_5["forbidden_predecessor_generation"], 12)
        self.assertEqual(
            retry_5["checkpoint_sha256"],
            "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8",
        )
        self.assertEqual(predecessor["cursor_mode"], "RETRY_PENDING")
        self.assertEqual(predecessor["pending_owner"], "temporal")
        self.assertEqual(predecessor["ledger_owner"], "TDG6TemporalLedger")
        self.assertEqual(predecessor["tracer_owner"], "NormalFlowTracers")
        self.assertEqual(predecessor["transaction_owner"], "GR0RuntimeStageTransaction")
        self.assertEqual(parsed["formulation"]["tableau_selector"], "SSPRK3")
        self.assertEqual(
            parsed["formulation"]["actual_spatial_operator"],
            "inherited_RK4_2049_SBP4",
        )
        self.assertEqual(parsed["formulation"]["channel_count"], 18)
        self.assertTrue(parsed["formulation"]["complete_admission_is_all_of"])
        self.assertTrue(
            parsed["formulation"]["both_widths_pass_licenses_production_method_design"]
        )
        self.assertTrue(
            parsed["formulation"][
                "any_width_nonpass_rejects_this_exact_remedy_on_tested_neighborhood"
            ]
        )
        self.assertFalse(parsed["transaction"]["campaign_state_write_authorized"])
        self.assertFalse(parsed["transaction"]["PDE_state_commit_authorized"])
        self.assertFalse(parsed["transaction"]["fine_path_commit_authorized"])
        self.assertFalse(parsed["transaction"]["diagnostic_fine_endpoint_serialized"])
        self.assertFalse(parsed["transaction"]["diagnostic_fine_endpoint_is_accepted_state"])
        self.assertEqual(parsed["transaction"]["width_count"], 2)
        self.assertEqual(parsed["transaction"]["maximum_stage_and_endpoint_records"], 56)
        self.assertFalse(any(parsed["exclusions"].values()))
        self.assertEqual(
            [int(item["retry"]) for item in authority.REPLAYS],
            [4, 5],
        )
        self.assertEqual(int(authority.REPLAYS[0]["predecessor_generation"]), 10)
        self.assertEqual(int(authority.REPLAYS[1]["predecessor_generation"]), 11)
        self.assertEqual(int(authority.REPLAYS[0]["journal_sequence"]), 13)
        self.assertEqual(int(authority.REPLAYS[1]["journal_sequence"]), 15)

    def test_wrong_hash_generation_retry_or_flag_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (authority.WIDTHS[0]["checkpoint_sha256"].encode(), b"0" * 64),
            (authority.WIDTHS[1]["checkpoint_sha256"].encode(), b"1" * 64),
            (authority.WIDTHS[0]["historical_rejection_sha256"].encode(), b"2" * 64),
            (authority.WIDTHS[1]["transaction_sha256"].encode(), b"3" * 64),
            (b"generation = 10", b"generation = 11"),
            (b"checkpoint_journal_sequence = 12", b"checkpoint_journal_sequence = 13"),
            (b"retry = 4", b"retry = 3"),
            (b'target_protocol = "FGC-2-SF1-PROTO18"', b'target_protocol = "FGC-2-SF1-PROTO19"'),
            (b'project_version = "0.11.0"', b'project_version = "0.12.0"'),
            (
                authority.QA1_PREF1_RESULT_SHA256.encode(),
                b"4" * 64,
            ),
            (b"campaign_state_write_authorized = false", b"campaign_state_write_authorized = true"),
            (b"physical_result_authorized = false", b"physical_result_authorized = true"),
            (b"diagnostic_fine_endpoint_serialized = false", b"diagnostic_fine_endpoint_serialized = true"),
            (b"width_robustness_passed = false", b"width_robustness_passed = true"),
            (b"old_member_adoption_authorized = false", b"old_member_adoption_authorized = true"),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.QA2AuthorityError):
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
        self.assertFalse(payload["width_robustness_passed"])
        self.assertFalse(payload["production_method_earned"])
        self.assertFalse(payload["diagnostic_fine_endpoint_serialized"])
        self.assertEqual(payload["target_protocol"], "FGC-2-SF1-PROTO18")
        self.assertEqual(payload["retry4_restores_generation"], 10)
        self.assertTrue(payload["retry4_never_restores_generation_11"])
        self.assertTrue(payload["checkpoint_journal_is_sequence_12"])
        self.assertTrue(payload["historical_retry4_rejection_is_sequence_13"])
        self.assertEqual(payload["retry5_restores_generation"], 11)
        self.assertTrue(payload["retry5_never_restores_generation_12"])
        self.assertTrue(payload["checkpoint_journal_is_sequence_14"])
        self.assertTrue(payload["historical_retry5_rejection_is_sequence_15"])
        self.assertTrue(payload["QA1_licenses_QA2_method_design_never_old_member_adoption"])
        self.assertNotIn("complete_admission_passed", payload)
        self.assertNotIn("failed_channels", payload)
        self.assertNotIn("diagnostic_payload_sha256", payload)
        self.assertEqual(
            payload["implementation"],
            authority.parse_config(config)["implementation"],
        )
        self.assertEqual(payload["project_version"], "0.11.0")
        self.assertEqual(raw, authority.canonical_pretty(parsed))
        self.assertEqual(len(payload["restored_predecessor_fingerprints"]), 2)

    def test_compact_rejects_nonfinite_json_constants(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token), self.assertRaisesRegex(
                authority.QA2AuthorityError,
                "nonfinite QA2 compact JSON constant",
            ):
                authority.validate_compact(
                    config, f'{{"gate_status": {token}}}\n'.encode()
                )

    def test_output_absence_rejects_existing_and_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / authority.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.QA2AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = (root / authority.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{authority.STAGING_PREFIX}dead").mkdir()
            with self.assertRaisesRegex(
                authority.QA2AuthorityError, "staging namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.QA2AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_generation_tips_match_historical_rejections(self) -> None:
        synthetic_store = "synthetic-store"
        synthetic_widths: list[dict[str, object]] = []
        leaves: dict[str, bytes] = {}
        for ordinal, original in enumerate(authority.WIDTHS, start=1):
            width = dict(original)
            checkpoint_sha256 = f"{ordinal:x}" * 64
            journal_sha256 = f"{ordinal + 2:x}" * 64
            rejection_sha256 = f"{ordinal + 4:x}" * 64
            checkpoint = json.dumps(
                {
                    "generation": width["generation"],
                    "checkpoint_sha256": checkpoint_sha256,
                    "journal_sequence": width["checkpoint_journal_sequence"],
                    "journal_tip_sha256": journal_sha256,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            journal = json.dumps(
                {
                    "kind": "tdg6_accept",
                    "record_sha256": journal_sha256,
                    "sequence": width["checkpoint_journal_sequence"],
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            rejection = json.dumps(
                {
                    "kind": "tdg6_rejection",
                    "record_sha256": rejection_sha256,
                    "sequence": width["historical_rejection_sequence"],
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            width.update(
                checkpoint_sha256=checkpoint_sha256,
                checkpoint_raw_sha256=sha256(checkpoint).hexdigest(),
                checkpoint_journal_sha256=journal_sha256,
                checkpoint_journal_raw_sha256=sha256(journal).hexdigest(),
                historical_rejection_sha256=rejection_sha256,
                historical_rejection_raw_sha256=sha256(rejection).hexdigest(),
            )
            generation = int(width["generation"])
            checkpoint_path = (
                f"{synthetic_store}/checkpoints/"
                f"{generation:020d}-{checkpoint_sha256}.json"
            )
            journal_path = (
                f"{synthetic_store}/journal/"
                f"{int(width['checkpoint_journal_sequence']):020d}-"
                f"{journal_sha256}.journal"
            )
            rejection_path = (
                f"{synthetic_store}/journal/"
                f"{int(width['historical_rejection_sequence']):020d}-"
                f"{rejection_sha256}.journal"
            )
            leaves.update(
                {
                    checkpoint_path: checkpoint,
                    journal_path: journal,
                    rejection_path: rejection,
                }
            )
            synthetic_widths.append(width)

        def synthetic_leaf(_root: Path, relative: str) -> bytes:
            return leaves[relative]

        with (
            patch.object(authority, "STORE_PATH", synthetic_store),
            patch.object(authority, "WIDTHS", tuple(synthetic_widths)),
            patch.object(authority, "read_leaf", side_effect=synthetic_leaf) as reader,
        ):
            authority._verify_predecessor_leaves(Path("/synthetic-unavailable"))
        self.assertEqual(reader.call_count, 12)
        self.assertEqual(authority.WIDTHS[0]["generation"], 10)
        self.assertEqual(authority.WIDTHS[0]["checkpoint_journal_sequence"], 12)
        self.assertEqual(authority.WIDTHS[0]["historical_rejection_sequence"], 13)
        self.assertEqual(authority.WIDTHS[1]["generation"], 11)
        self.assertEqual(authority.WIDTHS[1]["checkpoint_journal_sequence"], 14)
        self.assertEqual(authority.WIDTHS[1]["historical_rejection_sequence"], 15)
        self.assertEqual(authority.FORBIDDEN_RETRY3_GENERATION, 9)
        self.assertNotEqual(
            authority.WIDTHS[0]["checkpoint_journal_sha256"],
            authority.WIDTHS[0]["historical_rejection_sha256"],
        )
        self.assertNotEqual(
            authority.WIDTHS[1]["checkpoint_journal_sha256"],
            authority.WIDTHS[1]["historical_rejection_sha256"],
        )

    def test_prelaunch_restores_fingerprints_without_proposal_construction(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        expected = [
            authority._expected_fingerprint_receipt(width) for width in authority.WIDTHS
        ]

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
            patch.object(authority, "_verify_qa1_pref1_result"),
            patch.object(authority, "_verify_predecessor_leaves"),
            patch.object(authority, "require_output_absent") as absent,
            patch.object(
                authority,
                "real_predecessor_fingerprints",
                return_value=tuple(expected),
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
        self.assertFalse(result["artifact_payload"]["width_robustness_passed"])
        self.assertEqual(result["artifact_payload"]["retry4_restores_generation"], 10)
        self.assertEqual(result["artifact_payload"]["retry5_restores_generation"], 11)

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
            self.assertRaisesRegex(authority.QA2AuthorityError, "replace refs"),
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
        with self.assertRaises(authority.QA2AuthorityError):
            authority.parse_config(malformed)
        missing = raw.split(b"\n[implementation]\n")[0] + b"\n"
        with self.assertRaises(authority.QA2AuthorityError):
            authority.parse_config(missing)
        mutated = dict(parsed)
        mutated["implementation"] = dict(parsed["implementation"])
        mutated["implementation"]["runner_sha256"] = "0" * 64
        with self.assertRaisesRegex(
            authority.QA2AuthorityError, "implementation binding differs"
        ):
            authority.verify_implementation_inventory(
                ROOT, mutated, commit=None
            )

    def test_read_leaf_rejects_parent_symlink_and_leaf_race(self) -> None:
        payload = b"qa2-leaf\n"
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            real = Path(foreign) / "leaf.txt"
            real.write_bytes(payload)
            (root / "parent").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.QA2AuthorityError, "authority parent is unsafe"
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
                authority.QA2AuthorityError, "unsafe authority leaf"
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
                patch.object(authority.qa1.os, "fstat", side_effect=racing_fstat),
                self.assertRaisesRegex(
                    authority.QA2AuthorityError, "authority leaf changed"
                ),
            ):
                authority.read_leaf(root, relative)

    def test_qa1_pref1_live_and_committed_bytes_match(self) -> None:
        authority._verify_qa1_pref1_result(ROOT)
        live = authority.read_leaf(ROOT, authority.QA1_PREF1_RESULT_PATH)
        self.assertEqual(sha256(live).hexdigest(), authority.QA1_PREF1_RESULT_SHA256)

    def test_delta_excludes_raw_store_and_qdrant(self) -> None:
        self.assertEqual(authority.BASE_COMMIT, "a7915925a23ceaff6a491b8d08c572a70d242985")
        self.assertEqual(len(authority.DELTA_PATHS), 16)
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
                "scripts/reproduce_fgc_tdg10_qa2_frz1.py",
                authority.RUNNER_PATH,
                authority.AUTHORITY_PATH,
                "tests/test_check_repo_tdg10_qa2_frz1.py",
                "tests/test_fgc_tdg10_qa2_authority.py",
                "tests/test_fgc_tdg10_qa2_runner.py",
            },
        )
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(any(path.startswith("runs/") for path in authority.DELTA_PATHS))
        self.assertNotIn(authority.OUTPUT_NAMESPACE, authority.DELTA_PATHS)
        self.assertNotIn("docs/fgc-tdg10-qa1-frz1.md", authority.DELTA_PATHS)
        self.assertNotIn("docs/fgc-tdg10-qa1-pref1.md", authority.DELTA_PATHS)

    def test_receipt_flags_are_diagnostic_only(self) -> None:
        receipt = authority.QA2Authority(authority_commit="b" * 40)
        self.assertTrue(receipt.diagnostic_qualification_only)
        self.assertFalse(receipt.state_advance_authorized)
        self.assertFalse(receipt.campaign_state_write_authorized)
        self.assertFalse(receipt.PDE_state_commit_authorized)
        self.assertFalse(receipt.fine_path_commit_authorized)
        self.assertFalse(receipt.diagnostic_fine_endpoint_serialized)
        self.assertFalse(receipt.candidate_execution_authorized)
        self.assertFalse(receipt.physical_result_authorized)
        self.assertFalse(receipt.new_protocol_id_authorized)
        self.assertFalse(receipt.production_method_earned)
        self.assertFalse(receipt.width_robustness_passed)
        self.assertEqual(receipt.target_protocol, "FGC-2-SF1-PROTO18")
        self.assertEqual(receipt.selected_retries, (4, 5))


if __name__ == "__main__":
    unittest.main()
