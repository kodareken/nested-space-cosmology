from __future__ import annotations

import ast
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.def1_stab1 import (  # noqa: E402
    DEF1_BOOLEAN_NAMES,
    ERROR_BUDGET_COMPONENTS,
)
from recursive_horizons.fgc.def1_stab1_freeze_contract import (  # noqa: E402
    ARTIFACT_ID,
    OWNER_PATHS,
    PREF1_ARTIFACT_ID,
    build_def1_stab1_freeze_contract,
)
from recursive_horizons.fgc.def1_stab1_qualification import (  # noqa: E402
    Def1Stab1QualificationError,
    IMP1_18_CHANNEL_ORDER,
    PROVIDER_ROUTES,
)
from recursive_horizons.fgc.def1_geometry_error import INPUT_NAMES  # noqa: E402


class Def1Stab1FreezeContractTests(unittest.TestCase):
    def test_contract_is_complete_but_every_gate_remains_closed(self) -> None:
        contract = build_def1_stab1_freeze_contract()
        payload = dict(contract.payload)
        self.assertEqual(payload["artifact_id"], ARTIFACT_ID)
        self.assertEqual(payload["pref1_artifact_id"], PREF1_ARTIFACT_ID)
        self.assertEqual(payload["geometry_input_order"], list(INPUT_NAMES))
        self.assertEqual(payload["imp1_channel_order"], list(IMP1_18_CHANNEL_ORDER))
        self.assertEqual(
            payload["error_budget_component_order"], list(ERROR_BUDGET_COMPONENTS)
        )
        self.assertEqual(payload["def1_boolean_names"], list(DEF1_BOOLEAN_NAMES))
        self.assertEqual(payload["coverage_control_count"], len(PROVIDER_ROUTES))
        self.assertEqual(payload["coverage_component_count"], 14)
        self.assertTrue(contract.conversion_instrument_contract_complete)
        for name in (
            "def1_error_map_passed",
            "candidate_or_control_trajectory_read",
            "global_pde_error_certified",
            "def1_booleans_evaluated",
            "rob1_passed",
            "holdout_authorized",
            "physical_claimed",
        ):
            self.assertFalse(getattr(contract, name))
        self.assertEqual(
            contract.payload_sha256,
            sha256(canonical_json_bytes(payload)).hexdigest(),
        )

    def test_owner_hashes_bind_unique_live_files(self) -> None:
        payload = dict(build_def1_stab1_freeze_contract().payload)
        hashes = payload["owner_hashes"]
        self.assertEqual(set(hashes), set(OWNER_PATHS))
        for relative, digest in hashes.items():
            self.assertEqual(digest, sha256((ROOT / relative).read_bytes()).hexdigest())

    def test_replace_cannot_promote_or_change_cached_payload(self) -> None:
        contract = build_def1_stab1_freeze_contract()
        for name in (
            "def1_error_map_passed",
            "candidate_or_control_trajectory_read",
            "global_pde_error_certified",
            "def1_booleans_evaluated",
            "rob1_passed",
            "holdout_authorized",
            "physical_claimed",
        ):
            with self.subTest(name=name):
                with self.assertRaises(Def1Stab1QualificationError):
                    replace(contract, **{name: True})
        changed = dict(contract.payload)
        changed["geometry_input_order"] = list(reversed(INPUT_NAMES))
        with self.assertRaisesRegex(
            Def1Stab1QualificationError, "payload differs"
        ):
            replace(contract, payload=changed)
        with self.assertRaisesRegex(
            Def1Stab1QualificationError, "payload hash differs"
        ):
            replace(contract, payload_sha256="0" * 64)

    def test_import_graph_has_no_trajectory_or_execution_owner(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/def1_stab1_freeze_contract.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        for forbidden in (
            "col1",
            "rob1",
            "run_fgc",
            "campaign_store",
            "pro20_ev1_runtime",
        ):
            self.assertFalse(any(forbidden in item.lower() for item in imports))


if __name__ == "__main__":
    unittest.main()
