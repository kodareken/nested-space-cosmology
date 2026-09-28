from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from recursive_horizons.fgc import sgb1_ctl1_sol1_frz1_certificate as cert

ROOT = Path(__file__).resolve().parents[1]


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _script():
    script_path = ROOT / "scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py"
    spec = importlib.util.spec_from_file_location("sol1_frz1_script", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SOL1FRZ1CertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = (ROOT / cert.CONFIG_PATH).read_bytes()
        self.result = (ROOT / cert.RESULT_PATH).read_bytes()

    def test_tracked_bundle_is_canonical_and_pinned(self) -> None:
        parsed = cert.verify_compact(ROOT)
        self.assertEqual(_sha(self.config), cert.CONFIG_SHA256)
        self.assertEqual(_sha(self.result), cert.RESULT_SHA256)
        self.assertEqual(parsed["classification"], cert.CLASSIFICATION)

    def test_tracked_bytes_match_direct_reconstruction(self) -> None:
        config, result = cert.compose_canonical_artifacts(ROOT)
        self.assertEqual(config, self.config)
        self.assertEqual(result, self.result)

    def test_nominal_branch_and_control_hashes_are_exact(self) -> None:
        parsed = cert.verify_compact(ROOT)
        nominal = parsed["nominal"]
        control = parsed["control"]
        self.assertEqual(nominal["classification"], "interval_inconclusive")
        self.assertEqual(
            nominal["obstruction"], "prefix_endpoint_not_inside_declared_tube"
        )
        self.assertEqual(nominal["unique_tubes"], 0)
        self.assertFalse(nominal["tiles_remaining_support"])
        self.assertEqual(nominal["contract_sha256"], cert.NOMINAL_CONTRACT_SHA256)
        self.assertEqual(control["contract_sha256"], cert.CONTROL_CONTRACT_SHA256)

    def test_claim_inventory_is_complete_and_nonpromoting(self) -> None:
        claims = cert.verify_compact(ROOT)["claims"]
        self.assertEqual(claims, cert.expected_claims())
        for name in cert.TRUE_CLAIMS:
            self.assertIs(claims[name], True)
        for name in cert.FALSE_CLAIMS:
            self.assertIs(claims[name], False)

    def test_hash_mutations_fail_closed(self) -> None:
        with self.assertRaises(cert.SGBLSOL1FRZ1Error):
            cert.validate_compact_result(self.config + b" ", self.result)
        with self.assertRaises(cert.SGBLSOL1FRZ1Error):
            cert.validate_compact_result(self.config, self.result[:-2] + b"x\n")

    def test_reclassification_and_claim_promotion_fail_closed(self) -> None:
        parsed = json.loads(self.result)
        for operation in ("classification", "obstruction", "unique", "claim"):
            changed = deepcopy(parsed)
            if operation == "classification":
                changed["nominal"]["classification"] = "unique_local_affine_orbit"
            elif operation == "obstruction":
                changed["nominal"]["obstruction"] = "jacobian_diagonal_contains_zero"
            elif operation == "unique":
                changed["nominal"]["unique_tubes"] = 1
            else:
                changed["claims"]["SGBL_branch_owned_and_healthy"] = True
            raw = cert.canonical_result(changed)
            with patch.object(cert, "RESULT_SHA256", _sha(raw)):
                with self.assertRaises(cert.SGBLSOL1FRZ1Error):
                    cert.validate_compact_result(self.config, raw)

    def test_policy_and_contract_hash_mutations_fail_closed(self) -> None:
        parsed = json.loads(self.result)
        for key, value in (
            ("tube_radius", "1/16"),
            ("continuation_steps", 17),
            ("nominal_contract", "0" * 64),
            ("control_contract", "1" * 64),
        ):
            changed = deepcopy(parsed)
            if key == "nominal_contract":
                changed["nominal"]["contract_sha256"] = value
            elif key == "control_contract":
                changed["control"]["contract_sha256"] = value
            else:
                changed["policy"][key] = value
            raw = cert.canonical_result(changed)
            with patch.object(cert, "RESULT_SHA256", _sha(raw)):
                with self.assertRaises(cert.SGBLSOL1FRZ1Error):
                    cert.validate_compact_result(self.config, raw)

    def test_compact_check_never_reconstructs_or_reads_other_paths(self) -> None:
        with patch.object(
            cert, "compose_canonical_artifacts", side_effect=AssertionError("reconstruct")
        ):
            cert.verify_compact(ROOT)

    def test_compact_bundle_passes_without_source_or_runs(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            for relative, raw in ((cert.CONFIG_PATH, self.config), (cert.RESULT_PATH, self.result)):
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(raw)
            cert.verify_compact(root)
            self.assertFalse((root / "src").exists())
            self.assertFalse((root / "runs").exists())

    def test_default_script_does_not_call_compose(self) -> None:
        module = _script()
        with patch.object(cert, "compose_canonical_artifacts") as compose:
            self.assertEqual(module.main([]), 0)
        compose.assert_not_called()

    def test_certificate_top_level_imports_are_compact_safe(self) -> None:
        source_path = ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_sol1_frz1_certificate.py"
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        imports = []
        for node in tree.body:
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        self.assertFalse(any("sgb1_ctl1_sol1" in item for item in imports))
        self.assertFalse(any("campaign" in item or "run_fgc" in item for item in imports))


if __name__ == "__main__":
    unittest.main()
