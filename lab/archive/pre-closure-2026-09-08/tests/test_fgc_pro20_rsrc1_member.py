"""One-member RSRC1 physical-factory scope controls."""

from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import pro20_rsrc1_member as member


ROOT = Path(__file__).resolve().parents[1]


class RSRC1MemberFactoryTests(unittest.TestCase):
    def test_invalid_member_refuses_before_static_decode(self) -> None:
        with patch.object(member, "_decode_inputs", side_effect=AssertionError()), self.assertRaises(
            member.Pro20RSRC1MemberError
        ):
            member.build_rsrc1_static_shell("RK4-16385", static_input_bytes={})

    def test_static_factory_constructs_exactly_one_selected_shell(self) -> None:
        fake = SimpleNamespace(key="RK4-2049")
        cal9 = {
            "physical_inputs": {},
            "numerics": {},
            "universal_thresholds": {},
            "primary_method": {},
            "comparator_method": {},
        }
        with (
            patch.object(member, "_decode_inputs", return_value=(cal9, {}, {}, {})),
            patch.object(
                member,
                "_validated_sources",
                return_value=(cal9, {("RK4", 2049): {}}, {}),
            ),
            patch.object(member, "_mapping", side_effect=lambda value, _label: value),
            patch.object(member, "_build_member", return_value=fake) as build,
            patch.object(member, "_template_identity", return_value={}),
            patch.object(
                member,
                "_digest",
                return_value=member._PERSISTED_TEMPLATE_SHA256["RK4-2049"],
            ),
        ):
            observed = member.build_rsrc1_static_shell(
                "RK4-2049", static_input_bytes={}
            )
        self.assertIs(observed, fake)
        build.assert_called_once()
        self.assertEqual(build.call_args.kwargs["point_count"], 2049)

    def test_new_campaign_identity_does_not_reuse_closed_store_identity(self) -> None:
        from recursive_horizons.fgc.evolution.pro20_source_factory import (
            CAMPAIGN_ID as CLOSED_CAMPAIGN_ID,
        )

        self.assertNotEqual(member.CAMPAIGN_ID, CLOSED_CAMPAIGN_ID)
        self.assertEqual(member.PROTOCOL, "FGC-2-SF1-PROTO19")

    def test_source_contains_no_six_member_construction_loop_or_store(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        joined = " ".join(imports)
        self.assertNotIn("pro20_ev1_store", joined)
        self.assertNotIn("campaign_store", joined)
        source = path.read_text(encoding="utf-8")
        self.assertNotIn("for key in MEMBER_KEYS", source)
        self.assertNotIn("publish_", source)
        self.assertNotIn("write_bytes", source)


if __name__ == "__main__":
    unittest.main()
