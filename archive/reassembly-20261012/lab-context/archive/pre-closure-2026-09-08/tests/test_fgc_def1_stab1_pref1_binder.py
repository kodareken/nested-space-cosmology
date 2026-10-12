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

from recursive_horizons.fgc import def1_stab1_pref1_binder as pref1  # noqa: E402
from recursive_horizons.fgc.def1_stab1 import DEF1_BOOLEAN_NAMES  # noqa: E402


CONFIG = ROOT / pref1.CONFIG_PATH
RESULT = ROOT / pref1.RESULT_PATH
REPRODUCER = ROOT / "scripts/reproduce_fgc_def1_stab1_pref1.py"
BINDER = ROOT / "src/recursive_horizons/fgc/def1_stab1_pref1_binder.py"
FORBIDDEN_FRZ1_MODULES = (
    "recursive_horizons.fgc.def1_stab1_frz1_certificate",
    "recursive_horizons.fgc.def1_stab1_freeze_contract",
)
FORBIDDEN_OWNER_MODULES = (
    "recursive_horizons.fgc.def1_stab1",
    "recursive_horizons.fgc.def1_geometry_error",
    "recursive_horizons.fgc.def1_stab1_providers",
    "recursive_horizons.fgc.def1_stab1_qualification",
    *FORBIDDEN_FRZ1_MODULES,
)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _load_result() -> dict:
    return json.loads(RESULT.read_text(encoding="ascii"))


class Def1Stab1Pref1BinderTests(unittest.TestCase):
    def test_tracked_bytes_match_live_reconstruction(self) -> None:
        config_raw, result_raw = pref1.compose_canonical_artifacts(ROOT)
        self.assertEqual(config_raw, CONFIG.read_bytes())
        self.assertEqual(result_raw, RESULT.read_bytes())
        self.assertEqual(_sha(config_raw), pref1.CONFIG_SHA256)
        self.assertEqual(_sha(result_raw), pref1.RESULT_SHA256)

    def test_compact_check_accepts_the_tracked_bundle(self) -> None:
        result = pref1.verify_compact(ROOT)
        self.assertEqual(result["artifact_id"], pref1.ARTIFACT_ID)
        self.assertEqual(result["classification"], pref1.CLASSIFICATION)
        self.assertEqual(result["freeze_commit"], pref1.FREEZE_COMMIT)
        self.assertEqual(result["freeze_parent"], pref1.FREEZE_PARENT)
        self.assertIs(result["def1_error_map_passed"], True)
        self.assertIs(result["map_readiness_only"], True)
        self.assertIs(result["def1_booleans_evaluated"], False)
        self.assertIs(result["candidate_or_control_trajectory_read"], False)
        self.assertIs(result["actual_error_bounds_evaluated"], False)
        self.assertIs(result["holdout_authorized"], False)
        self.assertIs(result["rob1_passed"], False)
        self.assertIs(result["physical_claimed"], False)
        self.assertIs(result["mechanism_claimed"], False)

    def test_true_and_false_claims_are_complete(self) -> None:
        result = _load_result()
        claims = result["claims"]
        self.assertEqual(set(claims), set(pref1.TRUE_CLAIMS) | set(pref1.FALSE_CLAIMS))
        for name in pref1.TRUE_CLAIMS:
            self.assertIs(claims[name], True, name)
        for name in pref1.FALSE_CLAIMS:
            self.assertIs(claims[name], False, name)
        self.assertEqual(tuple(result["def1_boolean_names"]), DEF1_BOOLEAN_NAMES)
        self.assertEqual(set(result["def1_booleans"]), set(DEF1_BOOLEAN_NAMES))
        self.assertTrue(all(value is None for value in result["def1_booleans"].values()))

    def test_map_readiness_is_not_a_trajectory_result(self) -> None:
        result = _load_result()
        self.assertIs(result["def1_error_map_passed"], True)
        self.assertIs(result["map_readiness_only"], True)
        self.assertIs(result["actual_error_bounds_evaluated"], False)
        self.assertIs(result["candidate_or_control_trajectory_read"], False)
        self.assertIs(result["controls"]["imp1_to_q"]["conditional_premise"], True)
        self.assertIs(result["controls"]["imp1_to_q"]["def1_error_map_passed"], False)
        self.assertIn("map_readiness_only", " ".join(result["nonclaims"]))
        self.assertIn("not a trajectory result", " ".join(result["nonclaims"]))

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
                "typed_refusal",
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
        self.assertTrue(q_error["missing_premises_refused"])
        self.assertFalse(q_error["used_measured_q"])
        self.assertTrue(controls["misner_sharp"]["typed_veto_control_passed"])
        self.assertTrue(controls["activation"]["control_dominance_passed"])
        self.assertTrue(controls["margins"]["trappedness_strict_pass"])
        self.assertTrue(controls["margins"]["trappedness_boundary_fails"])
        self.assertTrue(controls["margins"]["complete_q_strict_pass"])
        self.assertTrue(controls["margins"]["complete_q_boundary_fails"])
        self.assertTrue(controls["base_to_adm"]["roundtrip_holds"])
        self.assertEqual(controls["base_to_adm"]["geometry_input_count"], 26)
        self.assertEqual(controls["imp1_to_q"]["channel_count"], 18)
        self.assertTrue(controls["typed_refusal"]["c_only_gauge_refused"])

    def test_remaining_work_taxonomy_is_independently_reconstructed(self) -> None:
        taxonomy = _load_result()["remaining_work_taxonomy"]
        self.assertEqual(
            taxonomy["reconstructed_from_providers"],
            list(pref1.PINNED_REMAINING_GATE_WORK),
        )
        self.assertEqual(len(taxonomy["reconstructed_from_providers"]), 9)
        self.assertIn(
            "independent qualification of conversion/error-map contract",
            " ".join(taxonomy["closed_by_this_binder"]),
        )
        self.assertIn("COL1", " ".join(taxonomy["later_application_not_map_readiness"]))
        self.assertIn(
            "typed refusal",
            " ".join(taxonomy["operational_premises"]),
        )
        self.assertIn(
            "universal nonlinear PDE theorem",
            " ".join(taxonomy["preserved_nonclaims"]),
        )

    def test_ordinary_check_is_reconstruction_git_and_trajectory_blind(self) -> None:
        verify_source = inspect.getsource(pref1.verify_compact)
        validate_source = inspect.getsource(pref1.validate_compact_result)
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
        self.assertNotIn("git_read", validate_source)
        self.assertNotIn("runs/", verify_source)
        self.assertNotIn("campaign", verify_source.lower())
        self.assertNotIn("ls-files", verify_source)
        self.assertFalse(any("compose_canonical_artifacts" in item for item in imported))
        with (
            patch.object(
                pref1,
                "compose_canonical_artifacts",
                side_effect=AssertionError("live reconstruction"),
            ),
            patch.object(
                pref1,
                "git_read",
                side_effect=AssertionError("git"),
            ),
        ):
            pref1.verify_compact(ROOT)

    def test_check_subprocess_never_loads_owner_frz1_or_campaign_modules(self) -> None:
        script = """
import sys
from pathlib import Path
root = Path(%r)
sys.path.insert(0, str(root / "src"))
from recursive_horizons.fgc.def1_stab1_pref1_binder import verify_compact
verify_compact(root)
forbidden = {
    "recursive_horizons.fgc.def1_stab1",
    "recursive_horizons.fgc.def1_geometry_error",
    "recursive_horizons.fgc.def1_stab1_providers",
    "recursive_horizons.fgc.def1_stab1_qualification",
    "recursive_horizons.fgc.def1_stab1_frz1_certificate",
    "recursive_horizons.fgc.def1_stab1_freeze_contract",
    "recursive_horizons.fgc.evolution.pro20_ev1_runtime",
    "recursive_horizons.fgc.evolution.hlt17_campaign_store",
}
loaded = forbidden.intersection(sys.modules)
assert not loaded, sorted(loaded)
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
        self.assertTrue(payload["def1_error_map_passed"])
        self.assertTrue(payload["map_readiness_only"])
        self.assertFalse(payload["holdout_authorized"])

    def test_binder_source_never_imports_frz1_certificate_or_contract(self) -> None:
        tree = ast.parse(BINDER.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        for forbidden in FORBIDDEN_FRZ1_MODULES:
            self.assertFalse(any(forbidden in item for item in imported), forbidden)
        self.assertFalse(any("campaign" in item.lower() for item in imported))
        self.assertFalse(any("pro20_ev1_runtime" in item for item in imported))

    def test_binder_module_has_no_write_or_publish_calls(self) -> None:
        source = BINDER.read_text(encoding="utf-8")
        self.assertNotIn("publish_exclusive_file", source)
        self.assertNotIn("write_text", source)
        self.assertNotIn("write_bytes", source)
        self.assertNotIn("open(", inspect.getsource(pref1.compose_canonical_artifacts))
        self.assertNotIn("open(", inspect.getsource(pref1.verify_compact))

    def _mutated(self, config: bytes, result: bytes) -> None:
        with (
            patch.object(pref1, "CONFIG_SHA256", _sha(config)),
            patch.object(pref1, "RESULT_SHA256", _sha(result)),
        ):
            pref1.validate_compact_result(config, result)

    def test_altered_canonical_bytes_fail_closed_on_hash(self) -> None:
        config = CONFIG.read_bytes()
        result = RESULT.read_bytes()
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "config SHA-256"):
            pref1.validate_compact_result(config + b"\n", result)
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "compact SHA-256"):
            pref1.validate_compact_result(config, result.replace(b"true", b"fals", 1))

    def test_owner_hash_and_route_coverage_mutations_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = json.loads(RESULT.read_bytes())
        poisoned = config.replace(
            b"7f51be0fd4b16ad1aac61596e9a5fa4ec99c29a1bf403e82aeb205f287b4e0d6",
            b"0" * 64,
            1,
        )
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "config SHA-256"):
            pref1.validate_compact_result(poisoned, RESULT.read_bytes())
        result["controls"]["provider_route_coverage"]["positive_controls_passed"] = False
        mutated = pref1.canonical_result(result)
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "positive route coverage"):
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
            (("claims", "map_readiness_only"), False, "map_readiness_only"),
            (("def1_booleans_evaluated",), True, "def1_booleans_evaluated"),
            (("holdout_authorized",), True, "holdout_authorized"),
            (("physical_claimed",), True, "physical_claimed"),
            (("actual_error_bounds_evaluated",), True, "actual_error_bounds_evaluated"),
            (("rob1_passed",), True, "rob1_passed"),
            (("mechanism_claimed",), True, "mechanism_claimed"),
        )
        for path, value, needle in cases:
            mutated = deepcopy(result)
            cursor = mutated
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            with self.subTest(path=path):
                with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, needle):
                    self._mutated(config, pref1.canonical_result(mutated))

    def test_fake_map_readiness_without_map_readiness_only_fails_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        result["map_readiness_only"] = False
        result["claims"]["map_readiness_only"] = False
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "map_readiness_only"):
            self._mutated(config, pref1.canonical_result(result))

    def test_partial_outputs_and_evaluated_booleans_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        result["controls"].pop("imp1_to_q")
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "control inventory"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["def1_booleans"]["resolved_activation"] = True
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "DEF1 boolean was evaluated"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["controls"].pop("q_error_assembly")
        result["controls"].pop("misner_sharp")
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "control inventory"):
            self._mutated(config, pref1.canonical_result(result))

    def test_altered_freeze_commit_parent_and_delta_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        result["freeze_commit"] = "0" * 40
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "freeze commit"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["freeze_parent"] = "0" * 40
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "freeze parent"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["freeze_evidence"]["delta"] = result["freeze_evidence"]["delta"][:-1]
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "freeze evidence delta"):
            self._mutated(config, pref1.canonical_result(result))

    def test_conditional_premise_and_measured_q_semantics_fail_closed(self) -> None:
        config = CONFIG.read_bytes()
        result = _load_result()
        result["controls"]["imp1_to_q"]["conditional_premise"] = False
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "IMP1 conditional premise"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["claims"]["universal_pde_theorem"] = True
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "universal_pde_theorem"):
            self._mutated(config, pref1.canonical_result(result))
        result = _load_result()
        result["controls"]["typed_refusal"]["c_only_gauge_refused"] = False
        with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "C-only refusal"):
            self._mutated(config, pref1.canonical_result(result))

    def test_live_bind_rejects_altered_freeze_parent(self) -> None:
        real_git = pref1._git

        def poisoned(root: Path, operation: object) -> object:
            if type(operation).__name__ == "CommitParents":
                return ("0" * 40,)
            return real_git(root, operation)

        with patch.object(pref1, "_git", side_effect=poisoned):
            with self.assertRaisesRegex(pref1.Def1Stab1Pref1Error, "freeze parent"):
                pref1.compose_canonical_artifacts(ROOT)

    def test_top_level_imports_exclude_def1_owners_and_frz1_decision(self) -> None:
        tree = ast.parse(BINDER.read_text(encoding="utf-8"))
        top_imports: list[str] = []
        for node in tree.body:
            if isinstance(node, ast.Import):
                top_imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                top_imports.append(node.module or "")
        self.assertTrue(any("evidence_io" in item for item in top_imports))
        for forbidden in FORBIDDEN_OWNER_MODULES:
            self.assertFalse(any(forbidden in item for item in top_imports), forbidden)


if __name__ == "__main__":
    unittest.main()
