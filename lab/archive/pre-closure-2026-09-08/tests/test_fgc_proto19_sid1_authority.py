from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution import proto19_sid1_authority as sid1  # noqa: E402


class SID1AuthorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        try:
            cls.anchor = sid1.derive_live_anchor(REPOSITORY)
        except sid1.SID1AuthorityError as exc:
            raise unittest.SkipTest(f"SID1 live anchor unavailable: {exc}")
        cls.config_raw = (REPOSITORY / sid1.CONFIG_PATH).read_bytes()

    def test_live_anchor_matches_exact_static_contract(self) -> None:
        self.assertEqual(self.anchor, sid1.expected_anchor())

    def test_result_is_deterministic(self) -> None:
        first = sid1.build_sid1_result(self.config_raw, self.anchor)
        second = sid1.build_sid1_result(self.config_raw, self.anchor)
        self.assertEqual(first, second)

    def test_anchor_rejects_each_hash_identity_mutation(self) -> None:
        for field in ("checkpoint_sha256", "checkpoint_journal_tip_sha256", "tdg6_suffix_sha256", "descriptor_sha256", "evolution_state_sha256"):
            mutated = copy.deepcopy(self.anchor)
            mutated[field] = "0" * 64
            with self.subTest(field=field):
                with self.assertRaises(sid1.SID1AuthorityError):
                    sid1.build_sid1_result(self.config_raw, mutated)

    def test_recursive_decoder_rejects_missing_nested_envelope(self) -> None:
        root = REPOSITORY / sid1.STORE_PATH
        raw = sid1._read(root, f"journal/{7:020d}-{sid1.TDG6_SUFFIX_SHA256}.journal", "test suffix")
        evidence = sid1._json(raw, "test suffix")["payload"]["evidence"]
        mutated = copy.deepcopy(evidence)
        del mutated["channel_admissions"][0]["outer_difference"]["envelope"]
        with self.assertRaises(sid1.SID1AuthorityError):
            sid1.decode_tdg6_evidence(mutated)

    def test_recursive_decoder_rejects_channel_order_mutation(self) -> None:
        root = REPOSITORY / sid1.STORE_PATH
        raw = sid1._read(root, f"journal/{7:020d}-{sid1.TDG6_SUFFIX_SHA256}.journal", "test suffix")
        evidence = sid1._json(raw, "test suffix")["payload"]["evidence"]
        mutated = copy.deepcopy(evidence)
        mutated["channel_admissions"][0], mutated["channel_admissions"][1] = mutated["channel_admissions"][1], mutated["channel_admissions"][0]
        with self.assertRaises(sid1.SID1AuthorityError):
            sid1.decode_tdg6_evidence(mutated)

    def test_live_anchor_rejects_suffix_parent_chain_mutation(self) -> None:
        root = REPOSITORY / sid1.STORE_PATH
        suffix_relative = f"journal/{7:020d}-{sid1.TDG6_SUFFIX_SHA256}.journal"
        original = sid1._read(root, suffix_relative, "test suffix")
        altered = json.loads(original.decode("ascii"))
        altered["previous_record_sha256"] = "0" * 64
        body = dict(altered)
        body.pop("record_sha256")
        altered["record_sha256"] = sid1._sha(sid1.canonical(body))
        mutated = sid1.canonical(altered)
        original_reader = sid1._read

        def reader(path_root, relative, label):
            if relative == suffix_relative:
                return mutated
            return original_reader(path_root, relative, label)

        with patch.object(sid1, "_read", side_effect=reader):
            with self.assertRaises(sid1.SID1AuthorityError):
                sid1.derive_live_anchor(REPOSITORY)

    def test_live_anchor_rejects_any_later_store_inventory(self) -> None:
        original = sid1._safe_inventory

        def inventory(root, relative, label):
            values = original(root, relative, label)
            if relative == "journal":
                return values + ("00000000000000000008-" + "0" * 64 + ".journal",)
            return values

        with patch.object(sid1, "_safe_inventory", side_effect=inventory):
            with self.assertRaisesRegex(sid1.SID1AuthorityError, "store moved"):
                sid1.derive_live_anchor(REPOSITORY)

    def test_live_anchor_is_read_only_for_its_explicit_leaves(self) -> None:
        root = REPOSITORY / sid1.STORE_PATH
        leaves = (
            f"checkpoints/{7:020d}-{sid1.CHECKPOINT_SHA256}.json",
            f"journal/{6:020d}-{sid1.CFL_RECORD_SHA256}.journal",
            f"journal/{7:020d}-{sid1.TDG6_SUFFIX_SHA256}.journal",
            f"states/{sid1.DESCRIPTOR_SHA256}.json",
            f"payloads/{sid1.SEMANTIC_PAYLOAD_SHA256}.npz",
        )
        before = {leaf: sid1._sha(sid1._read(root, leaf, "before")) for leaf in leaves}
        self.assertEqual(sid1.derive_live_anchor(REPOSITORY), sid1.expected_anchor())
        after = {leaf: sid1._sha(sid1._read(root, leaf, "after")) for leaf in leaves}
        self.assertEqual(after, before)

    def test_config_rejects_scope_promotion(self) -> None:
        import tomllib
        config = tomllib.loads(self.config_raw.decode("utf-8"))
        config["scope"]["raw_store_mutation"] = True
        with self.assertRaises(sid1.SID1AuthorityError):
            sid1._validate_config(config)

    def _authorization_fixture(self):
        result_raw = sid1._canonical_result(
            sid1.build_sid1_result(self.config_raw, self.anchor)
        )
        old_pro19 = b'{"artifact_id":"FGC-1-PRO19-FRZ1"}\n'
        old_auth1 = b'{"artifact_id":"FGC-1-PRO18-AUTH1"}\n'
        values = {
            sid1.CONFIG_PATH: self.config_raw,
            sid1.RESULT_PATH: result_raw,
            "results/fgc-1-pro19-frz1.json": old_pro19,
            "results/fgc-1-pro18-auth1.json": old_auth1,
        }
        paths = tuple(
            (
                "config" if path.startswith("configs/") else "result",
                path,
                sid1._sha(raw),
            )
            for path, raw in sorted(values.items())
        )
        launch = SimpleNamespace(
            progression_plan=SimpleNamespace(branch="GR-0"),
            manifest_sha256="1" * 64,
            authority_paths=paths,
            environment={"python_version": "test"},
        )
        original = SimpleNamespace(
            sha256=sid1.PLAN_SHA256,
            campaign_id=sid1.CAMPAIGN_ID,
            branch="GR-0",
            amplitude="3",
        )
        return values, launch, original

    def test_committed_image_authorizes_only_the_exact_sid1_edge(self) -> None:
        values, launch, original = self._authorization_fixture()
        with patch.object(
            sid1, "authorize_first_event", return_value=launch
        ), patch.object(
            sid1,
            "_nofollow_regular_bytes",
            side_effect=lambda _root, path: values[path],
        ), patch.object(
            sid1, "construct_first_event", return_value=original
        ):
            receipt = sid1.authorize_sid1_recovery(
                REPOSITORY,
                continuation_commit="2" * 40,
                launch_manifest_path="manifest.toml",
                store_anchor=self.anchor,
            )
        self.assertEqual(receipt.recovery_checkpoint_generation, 7)
        self.assertEqual(receipt.recovery_checkpoint_sha256, sid1.CHECKPOINT_SHA256)
        self.assertEqual(receipt.suffix_sequence, 7)
        self.assertEqual(receipt.suffix_sha256, sid1.TDG6_SUFFIX_SHA256)
        self.assertEqual(receipt.descriptor_sha256, sid1.DESCRIPTOR_SHA256)
        self.assertEqual(
            receipt.evolution_state_sha256, sid1.EVOLUTION_STATE_SHA256
        )
        self.assertIs(receipt.progression_plan, original)

    def test_committed_image_rejects_missing_sid1_binding_or_wrong_original_plan(self) -> None:
        values, launch, original = self._authorization_fixture()
        omitted = SimpleNamespace(
            **{
                **vars(launch),
                "authority_paths": tuple(
                    item
                    for item in launch.authority_paths
                    if item[1] != sid1.RESULT_PATH
                ),
            }
        )
        with patch.object(
            sid1, "authorize_first_event", return_value=omitted
        ):
            with self.assertRaisesRegex(
                sid1.SID1AuthorityError, "outside the launch manifest"
            ):
                sid1.authorize_sid1_recovery(
                    REPOSITORY,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )

        wrong = SimpleNamespace(**{**vars(original), "sha256": "0" * 64})
        with patch.object(
            sid1, "authorize_first_event", return_value=launch
        ), patch.object(
            sid1,
            "_nofollow_regular_bytes",
            side_effect=lambda _root, path: values[path],
        ), patch.object(
            sid1, "construct_first_event", return_value=wrong
        ):
            with self.assertRaisesRegex(
                sid1.SID1AuthorityError, "original progression plan"
            ):
                sid1.authorize_sid1_recovery(
                    REPOSITORY,
                    continuation_commit="2" * 40,
                    launch_manifest_path="manifest.toml",
                    store_anchor=self.anchor,
                )


if __name__ == "__main__":
    unittest.main()
