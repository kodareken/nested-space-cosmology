"""Focused and adversarial controls for the independent PRO20 PREF1 certificate."""

from __future__ import annotations

import ast
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import (  # noqa: E402
    pro20_ev1_pref1_certificate as certificate,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17_MEMBER_KEYS,
)


SCRIPT = ROOT / "scripts/reproduce_fgc_pro20_ev1_pref1.py"
CERT = "recursive_horizons.fgc.evolution.pro20_ev1_pref1_certificate"
FORBIDDEN_OWNERS = (
    "pro20_ev1_runtime",
    "pro20_ev1_protocol",
    "pro20_ev1_store",
    "pro20_ev1_authority",
)


def _load_script():
    spec = importlib.util.spec_from_file_location("pref1_reproducer", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_bundles() -> dict[str, object]:
    return {key: object() for key in HLT17_MEMBER_KEYS}


class PREF1CertificateTests(unittest.TestCase):
    def test_compact_check_is_absent_without_tracked_files(self) -> None:
        with TemporaryDirectory() as tmp:
            with self.assertRaises(certificate.PRO20EV1PREF1CertificateError) as raised:
                certificate.verify_compact(Path(tmp))
        self.assertIn("absent", str(raised.exception))

    def test_validate_compact_accepts_pinned_certificate(self) -> None:
        config = certificate.emit_config_bytes()
        result = certificate.canonical_result(certificate.expected_compact())
        parsed = certificate.validate_compact_result(config, result)
        self.assertEqual(
            parsed["artifact_payload"]["classification"],
            certificate.CLASSIFICATION,
        )
        self.assertEqual(
            sha256_hex := __import__("hashlib").sha256(result).hexdigest(),
            certificate.COMPACT_RESULT_SHA256,
        )
        del sha256_hex
        with patch.object(certificate, "COMPACT_RESULT_SHA256", None):
            with self.assertRaisesRegex(
                certificate.PRO20EV1PREF1CertificateError, "not canonical"
            ):
                certificate.validate_compact_result(config, b" " + result)

    def test_check_mode_never_accesses_git_raw_store_source_or_planck(self) -> None:
        config = certificate.emit_config_bytes()
        result = certificate.canonical_result(certificate.expected_compact())
        with TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            (root / "configs/fgc").mkdir(parents=True)
            (root / "results").mkdir()
            (root / certificate.CONFIG_PATH).write_bytes(config)
            (root / certificate.RESULT_PATH).write_bytes(result)
            with (
                patch.object(certificate, "git_read") as git_read,
                patch.object(certificate, "inventory_historical_runs") as runs,
                patch.object(certificate, "snapshot_store", create=True) as snap,
                patch.object(certificate, "bind_terminal_store", create=True) as store,
                patch.object(certificate, "capture_pro20_historical_origin", create=True),
                patch.object(certificate, "build_pro20_runtime_origin", create=True),
                patch.object(certificate, "authenticate_planck") as planck,
                patch.object(certificate, "rebuild_independent_seed_bundles") as rebuild,
            ):
                certificate.verify_compact(root)
            git_read.assert_not_called()
            runs.assert_not_called()
            snap.assert_not_called()
            store.assert_not_called()
            planck.assert_not_called()
            rebuild.assert_not_called()
            self.assertFalse((root / "runs").exists())

    def test_default_script_is_check_blind_and_does_not_reconstruct(self) -> None:
        script = _load_script()
        with (
            patch.object(certificate, "bind_live") as live,
            patch.object(certificate, "git_read") as git_read,
            patch.object(certificate, "inventory_historical_runs") as runs,
        ):
            code = script.main([])
        self.assertEqual(code, 0)
        live.assert_not_called()
        git_read.assert_not_called()
        runs.assert_not_called()

    def test_explicit_check_flag_stays_compact_only(self) -> None:
        script = _load_script()
        with patch.object(certificate, "bind_live") as live:
            code = script.main(["--check"])
        self.assertEqual(code, 0)
        live.assert_not_called()

    def test_altered_compact_authority_is_rejected(self) -> None:
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["authority"]["authority_commit"] = "0" * 40
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.validate_compact_result(
                certificate.emit_config_bytes(),
                certificate.canonical_result(compact),
            )

    def test_scientific_nonpass_or_common_event_class_is_rejected(self) -> None:
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["classification"] = (
            "independently_bound_scientific_nonpass"
        )
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.validate_compact_result(
                certificate.emit_config_bytes(),
                certificate.canonical_result(compact),
            )
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["classification"] = (
            "independently_bound_six_member_first_event_complete"
        )
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.validate_compact_result(
                certificate.emit_config_bytes(),
                certificate.canonical_result(compact),
            )

    def test_promoted_physics_or_wave2_claims_are_rejected(self) -> None:
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["claims"]["physics_claimed"] = True
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.validate_compact_result(
                certificate.emit_config_bytes(),
                certificate.canonical_result(compact),
            )
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["pro21_wave2_authorized"] = True
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.validate_compact_result(
                certificate.emit_config_bytes(),
                certificate.canonical_result(compact),
            )

    def test_classification_requires_every_independent_predicate(self) -> None:
        facts = {
            "authority_bound": True,
            "freeze_bound": True,
            "root_bound": True,
            "store_bound": True,
            "terminal_bound": True,
            "historical_runs_bound": True,
            "planck_bound": True,
            "environment_bound": True,
            "source_bound": True,
            "seed_independently_rebuilt": True,
            "resource_exhausted_terminal": True,
            "no_common_event": True,
            "partial_accepted_state": True,
            "runner_label_sufficient": False,
        }
        self.assertEqual(
            certificate.classification_if_established(facts),
            certificate.CLASSIFICATION,
        )
        for name in ("resource_exhausted_terminal", "no_common_event", "seed_independently_rebuilt"):
            broken = dict(facts)
            broken[name] = False
            with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
                certificate.classification_if_established(broken)
        promoted = dict(facts)
        promoted["runner_label_sufficient"] = True
        with self.assertRaises(certificate.PRO20EV1PREF1CertificateError):
            certificate.classification_if_established(promoted)

    def test_partial_accepted_state_is_explicit_in_compact(self) -> None:
        payload = certificate.expected_compact_payload()
        self.assertIs(payload["progression"]["partial_accepted_state"], True)
        self.assertIs(payload["progression"]["accepted_state_advanced"], True)
        self.assertIs(payload["progression"]["common_event_completed"], False)
        self.assertEqual(
            payload["progression"]["accepted_generation_in_canonical_member_order"],
            [14, 28, 39, 29, 36, 0],
        )
        self.assertEqual(payload["progression"]["member_keys"], list(HLT17_MEMBER_KEYS))
        self.assertEqual(payload["progression"]["terminal_reason"], "max RSS exceeded")
        self.assertEqual(
            payload["progression"]["ssprk3_16385_accepted_time_hex"],
            "0x1.7000000000000p+0",
        )
        self.assertEqual(
            payload["seed"]["identity_stream_sha256"],
            certificate.SEED_IDENTITY_STREAM_SHA256,
        )
        self.assertIs(payload["production_campaign_store_present"], True)
        self.assertIs(payload["production_store_written_by_event"], True)
        self.assertEqual(payload["accepted_fine_states_persisted"], 146)
        self.assertIs(payload["binder_mutated_store"], False)
        self.assertIs(payload["binder_wrote_campaign_store"], False)
        self.assertIs(payload["binder_serialized_or_adopted_endpoint"], False)
        self.assertIs(payload["pro21_wave2_authorized"], False)
        self.assertIs(payload["claims"]["calibration_authorized"], False)
        self.assertIs(payload["claims"]["wave2_authorized"], False)
        self.assertIs(payload["claims"]["common_event_completed"], False)
        for key in certificate.FORBIDDEN_UNQUALIFIED_KEYS:
            self.assertNotIn(key, payload)
            self.assertNotIn(key, payload["claims"])
            self.assertNotIn(key, certificate.expected_config())

    def test_unqualified_store_write_fields_are_rejected(self) -> None:
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["campaign_store_written"] = False
        with patch.object(certificate, "COMPACT_RESULT_SHA256", None):
            with self.assertRaisesRegex(
                certificate.PRO20EV1PREF1CertificateError, "unqualified store field"
            ):
                certificate.validate_compact_result(
                    certificate.emit_config_bytes(),
                    certificate.canonical_result(compact),
                )
        compact = deepcopy(certificate.expected_compact())
        compact["artifact_payload"]["claims"]["endpoint_adopted_or_serialized"] = False
        with patch.object(certificate, "COMPACT_RESULT_SHA256", None):
            with self.assertRaisesRegex(
                certificate.PRO20EV1PREF1CertificateError, "unqualified store field"
            ):
                certificate.validate_compact_result(
                    certificate.emit_config_bytes(),
                    certificate.canonical_result(compact),
                )

    def test_import_graph_excludes_production_decision_owners(self) -> None:
        path = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/pro20_ev1_pref1_certificate.py"
        )
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        for forbidden in FORBIDDEN_OWNERS:
            self.assertFalse(
                any(forbidden in module for module in imported),
                f"certificate imports forbidden production owner {forbidden}",
            )
        self.assertTrue(any("pro20_ev1_pref1_binder" in module for module in imported))

    def test_synthetic_live_path_is_no_write_and_does_not_import_runners(self) -> None:
        env = {key: "x" for key in certificate.ENVIRONMENT_KEYS}
        env["executable_sha256"] = "a" * 64
        env["numpy_extension_sha256"] = "b" * 64
        terminal = object()
        bundles = _synthetic_bundles()
        imported = set(sys.modules)
        with TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            before = list(root.rglob("*"))
            with (
                patch.object(
                    certificate,
                    "environment_sha256",
                    return_value=certificate.ENVIRONMENT_SHA256,
                ),
                patch.object(certificate, "authenticate_authority") as auth,
                patch.object(certificate, "bind_terminal_store") as store,
                patch.object(certificate, "_compare_terminal") as compare,
                patch.object(certificate, "bind_independent_seed") as seed,
                patch.object(certificate, "inventory_historical_runs") as runs,
                patch.object(certificate, "authenticate_planck") as planck,
                patch.object(certificate, "rebuild_independent_seed_bundles") as rebuild,
                patch.object(certificate, "store_identity") as identity,
            ):
                compare.return_value = {"outcome": {}, "progression": {}}
                seed.return_value = {
                    "identity_stream_sha256": certificate.SEED_IDENTITY_STREAM_SHA256,
                    "independently_reconstructed_origin_matches_raw_seed": True,
                }
                runs.return_value = {
                    "file_count": certificate.HISTORICAL_RUNS_FILE_COUNT,
                    "byte_count": certificate.HISTORICAL_RUNS_BYTE_COUNT,
                    "digest": certificate.HISTORICAL_RUNS_DIGEST,
                    "excluded_namespaces": list(
                        certificate.EXCLUDED_HISTORICAL_NAMESPACES
                    ),
                }
                planck.return_value = {}
                identity.return_value = {"lexical_sha256": "c" * 64}
                result = certificate.bind_live(
                    root,
                    environment=env,
                    skip_git=True,
                    skip_store=True,
                    skip_source=True,
                    terminal=terminal,
                    seed_bundles=bundles,
                    historical_runs=runs.return_value,
                    planck={},
                )
            after = list(root.rglob("*"))
        self.assertEqual(before, after)
        self.assertFalse(result["written"])
        self.assertEqual(result["classification"], certificate.CLASSIFICATION)
        auth.assert_not_called()
        store.assert_not_called()
        rebuild.assert_not_called()
        leaked = [
            name
            for name in sys.modules
            if name not in imported
            and any(needle in name for needle in FORBIDDEN_OWNERS)
        ]
        self.assertEqual(leaked, [])

    def test_live_store_identity_is_unchanged_when_namespace_present(self) -> None:
        live = ROOT
        namespace = live.joinpath(*certificate.PRODUCTION_NAMESPACE.split("/"))
        if not namespace.is_dir():
            self.skipTest("closed production namespace is not present")
        before = certificate.store_identity(live)
        env = {key: "x" for key in certificate.ENVIRONMENT_KEYS}
        env["executable_sha256"] = "a" * 64
        env["numpy_extension_sha256"] = "b" * 64
        terminal = object()
        bundles = _synthetic_bundles()
        with (
            patch.object(
                certificate,
                "environment_sha256",
                return_value=certificate.ENVIRONMENT_SHA256,
            ),
            patch.object(certificate, "_compare_terminal") as compare,
            patch.object(certificate, "bind_independent_seed") as seed,
            patch.object(
                certificate,
                "inventory_historical_runs",
                return_value={
                    "file_count": certificate.HISTORICAL_RUNS_FILE_COUNT,
                    "byte_count": certificate.HISTORICAL_RUNS_BYTE_COUNT,
                    "digest": certificate.HISTORICAL_RUNS_DIGEST,
                    "excluded_namespaces": list(
                        certificate.EXCLUDED_HISTORICAL_NAMESPACES
                    ),
                },
            ),
            patch.object(certificate, "authenticate_planck", return_value={}),
        ):
            compare.return_value = {"outcome": {}, "progression": {}}
            seed.return_value = {
                "identity_stream_sha256": certificate.SEED_IDENTITY_STREAM_SHA256,
                "independently_reconstructed_origin_matches_raw_seed": True,
            }
            certificate.bind_live(
                live,
                environment=env,
                skip_git=True,
                skip_store=True,
                skip_source=True,
                terminal=terminal,
                seed_bundles=bundles,
                historical_runs={
                    "file_count": certificate.HISTORICAL_RUNS_FILE_COUNT,
                    "byte_count": certificate.HISTORICAL_RUNS_BYTE_COUNT,
                    "digest": certificate.HISTORICAL_RUNS_DIGEST,
                },
                planck={},
            )
        after = certificate.store_identity(live)
        self.assertEqual(before, after)
        self.assertEqual(before["lexical_sha256"], certificate.STORE_LEXICAL_SHA256)
        self.assertEqual(before["leaf_count"], certificate.STORE_LEAF_COUNT)
        self.assertEqual(before["byte_count"], certificate.STORE_BYTE_COUNT)


if __name__ == "__main__":
    unittest.main()
