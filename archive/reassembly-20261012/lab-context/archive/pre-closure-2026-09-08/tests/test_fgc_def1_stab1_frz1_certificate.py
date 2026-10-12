from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
import inspect
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc import def1_stab1_frz1_certificate as cert  # noqa: E402
from recursive_horizons.fgc.def1_stab1 import DEF1_BOOLEAN_NAMES  # noqa: E402
from recursive_horizons.fgc.def1_stab1_freeze_contract import (  # noqa: E402
    build_def1_stab1_freeze_contract,
)


CONFIG = ROOT / cert.CONFIG_PATH
RESULT = ROOT / cert.RESULT_PATH
REPRODUCER = ROOT / "scripts/reproduce_fgc_def1_stab1_frz1.py"
FORBIDDEN_OWNER_MODULES = (
    "recursive_horizons.fgc.def1_stab1",
    "recursive_horizons.fgc.def1_geometry_error",
    "recursive_horizons.fgc.def1_stab1_providers",
    "recursive_horizons.fgc.def1_stab1_qualification",
    "recursive_horizons.fgc.def1_stab1_freeze_contract",
)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _load_result() -> dict:
    return json.loads(RESULT.read_text(encoding="ascii"))


class Def1Stab1Frz1CertificateTests(unittest.TestCase):
    def test_tracked_bytes_match_live_reconstruction(self) -> None:
        config_raw, result_raw = cert.compose_canonical_artifacts(ROOT)
        self.assertEqual(config_raw, CONFIG.read_bytes())
        self.assertEqual(result_raw, RESULT.read_bytes())
        self.assertEqual(_sha(config_raw), cert.CONFIG_SHA256)
        self.assertEqual(_sha(result_raw), cert.RESULT_SHA256)

    def test_compact_check_accepts_the_tracked_bundle(self) -> None:
        result = cert.verify_compact(ROOT)
        self.assertEqual(result["artifact_id"], cert.ARTIFACT_ID)
        self.assertEqual(result["classification"], cert.CLASSIFICATION)
        self.assertEqual(result["base_commit"], cert.BASE_COMMIT)
        self.assertTrue(result["claims"]["conversion_instrument_contract_complete"])
        self.assertFalse(result["def1_error_map_passed"])
        self.assertFalse(result["def1_booleans_evaluated"])
        self.assertTrue(result["licenses_separate_independent_pref1_binder"])
        self.assertFalse(result["pref1_implemented"])
        self.assertEqual(result["pref1_artifact_id"], cert.PREF1_ARTIFACT_ID)

    def test_true_and_false_claims_are_complete(self) -> None:
        result = _load_result()
        claims = result["claims"]
        self.assertEqual(set(claims), set(cert.TRUE_CLAIMS) | set(cert.FALSE_CLAIMS))
        for name in cert.TRUE_CLAIMS:
            self.assertIs(claims[name], True, name)
        for name in cert.FALSE_CLAIMS:
            self.assertIs(claims[name], False, name)
        self.assertEqual(tuple(result["def1_boolean_names"]), DEF1_BOOLEAN_NAMES)
        self.assertEqual(set(result["def1_booleans"]), set(DEF1_BOOLEAN_NAMES))
        self.assertTrue(all(value is None for value in result["def1_booleans"].values()))

    def test_controls_cover_every_required_instrument(self) -> None:
        controls = _load_result()["controls"]
        self.assertEqual(
            set(controls),
            {
                "conversion_identities",
                "provider_route_coverage",
                "q_error_assembly",
                "misner_sharp",
                "activation",
                "margins",
                "base_to_adm",
                "imp1_to_q",
            },
        )
        coverage = controls["provider_route_coverage"]
        self.assertEqual(coverage["route_count"], 21)
        self.assertEqual(coverage["component_count"], 14)
        self.assertTrue(coverage["positive_controls_passed"])
        self.assertTrue(coverage["injected_failures_refused"])
        q_error = controls["q_error_assembly"]
        self.assertEqual(q_error["mixed_total"], {"numerator": 37, "denominator": 30})
        self.assertTrue(q_error["measured_q_argument_refused"])
        self.assertFalse(q_error["used_measured_q"])
        self.assertTrue(controls["misner_sharp"]["typed_veto_control_passed"])
        self.assertTrue(controls["activation"]["control_dominance_passed"])
        self.assertTrue(controls["margins"]["trappedness_strict_pass"])
        self.assertTrue(controls["margins"]["complete_q_strict_pass"])
        self.assertTrue(controls["base_to_adm"]["roundtrip_holds"])
        self.assertEqual(controls["base_to_adm"]["geometry_input_count"], 26)
        self.assertEqual(controls["imp1_to_q"]["channel_count"], 18)
        self.assertFalse(controls["imp1_to_q"]["def1_error_map_passed"])
        contract = build_def1_stab1_freeze_contract()
        self.assertEqual(
            _load_result()["freeze_contract_payload_sha256"],
            contract.payload_sha256,
        )

    def test_ordinary_check_is_reconstruction_git_and_trajectory_blind(self) -> None:
        verify_source = inspect.getsource(cert.verify_compact)
        validate_source = inspect.getsource(cert.validate_compact_result)
        reproducer = ast.parse(REPRODUCER.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(reproducer):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.col_offset == 0:
                imported.append(node.module or "")
        self.assertNotIn("compose_canonical_artifacts", verify_source)
        self.assertNotIn("compose_canonical_artifacts", validate_source)
        self.assertNotIn("git_read", verify_source)
        self.assertNotIn("runs/", verify_source)
        self.assertNotIn("campaign", verify_source.lower())
        self.assertNotIn("ls-files", verify_source)
        self.assertFalse(any("compose_canonical_artifacts" in item for item in imported))
        with patch.object(
            cert,
            "compose_canonical_artifacts",
            side_effect=AssertionError("live reconstruction"),
        ):
            cert.verify_compact(ROOT)

    def test_check_subprocess_never_loads_owner_or_campaign_modules(self) -> None:
        script = """
import sys
from pathlib import Path
root = Path(%r)
sys.path.insert(0, str(root / "src"))
from recursive_horizons.fgc.def1_stab1_frz1_certificate import verify_compact
verify_compact(root)
forbidden = {
    "recursive_horizons.fgc.def1_stab1",
    "recursive_horizons.fgc.def1_geometry_error",
    "recursive_horizons.fgc.def1_stab1_providers",
    "recursive_horizons.fgc.def1_stab1_qualification",
    "recursive_horizons.fgc.def1_stab1_freeze_contract",
    "recursive_horizons.fgc.evolution.pro20_ev1_runtime",
    "recursive_horizons.fgc.evolution.hlt17_campaign_store",
}
loaded = forbidden.intersection(sys.modules)
assert not loaded, sorted(loaded)
assert not any(name.startswith("git") for name in sys.modules)
""" % str(ROOT)
        completed = subprocess.run(
            [sys.executable, "-I", "-B", "-c", script],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    def test_reproducer_check_does_not_call_compose(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-I", "-B", str(REPRODUCER), "--check"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["mode"], "check")
        self.assertFalse(payload["def1_error_map_passed"])

    def _mutated(self, config: bytes, result: bytes) -> None:
        with (
            patch.object(cert, "CONFIG_SHA256", _sha(config)),
            patch.object(cert, "RESULT_SHA256", _sha(result)),
        ):
            cert.validate_compact_result(config, result)

    def test_altered_canonical_bytes_fail_closed_on_hash(self) -> None:
        config = CONFIG.read_bytes()
        result = RESULT.read_bytes()
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "config SHA-256"):
            cert.validate_compact_result(config + b"\n", result)
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "compact SHA-256"):
            cert.validate_compact_result(config, result.replace(b"false", b"true", 1))

    def test_owner_hash_and_route_coverage_mutations_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = json.loads(RESULT.read_bytes())
        poisoned = config.replace(
            b"7f51be0fd4b16ad1aac61596e9a5fa4ec99c29a1bf403e82aeb205f287b4e0d6",
            b"0" * 64,
            1,
        )
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "config SHA-256"):
            cert.validate_compact_result(poisoned, RESULT.read_bytes())
        result["controls"]["provider_route_coverage"]["positive_controls_passed"] = False
        mutated = cert.canonical_result(result)
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "positive route coverage"):
            self._mutated(config, mutated)

    def test_q_assembly_mass_flux_margin_and_promotion_mutations_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        cases = (
            (("controls", "q_error_assembly", "used_measured_q"), True, "used_measured_q"),
            (
                ("controls", "misner_sharp", "typed_veto_control_passed"),
                False,
                "mass-flux veto",
            ),
            (
                ("controls", "margins", "trappedness_strict_pass"),
                False,
                "trappedness pass",
            ),
            (("claims", "def1_error_map_passed"), True, "def1_error_map_passed"),
            (("def1_error_map_passed",), True, "def1_error_map_passed"),
            (("def1_booleans_evaluated",), True, "def1_booleans_evaluated"),
            (("holdout_authorized",), True, "holdout_authorized"),
            (("physical_claimed",), True, "physical_claimed"),
        )
        for path, value, needle in cases:
            mutated = deepcopy(result)
            cursor = mutated
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            with self.subTest(path=path):
                with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, needle):
                    self._mutated(config, cert.canonical_result(mutated))

    def test_partial_outputs_and_evaluated_booleans_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        result["controls"].pop("imp1_to_q")
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "control inventory"):
            self._mutated(config, cert.canonical_result(result))
        result = _load_result()
        result["def1_booleans"]["resolved_activation"] = True
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "DEF1 boolean was evaluated"):
            self._mutated(config, cert.canonical_result(result))
        result = _load_result()
        result["controls"].pop("q_error_assembly")
        result["controls"].pop("misner_sharp")
        with self.assertRaisesRegex(cert.Def1Stab1Frz1Error, "control inventory"):
            self._mutated(config, cert.canonical_result(result))

    def test_compose_import_graph_stays_inside_def1_owners(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/def1_stab1_frz1_certificate.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        top_imports: list[str] = []
        for node in tree.body:
            if isinstance(node, ast.Import):
                top_imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                top_imports.append(node.module or "")
        self.assertTrue(any("evidence_io" in item for item in top_imports))
        for forbidden in FORBIDDEN_OWNER_MODULES:
            self.assertFalse(any(forbidden in item for item in top_imports))
        self.assertFalse(any("campaign" in item.lower() for item in top_imports))
        self.assertFalse(any("pro20_ev1_runtime" in item for item in top_imports))


if __name__ == "__main__":
    unittest.main()
