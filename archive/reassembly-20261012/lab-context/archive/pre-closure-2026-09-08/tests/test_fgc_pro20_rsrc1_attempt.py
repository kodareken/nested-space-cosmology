"""Store-blind RSRC1 one-attempt parity and failure controls."""

from __future__ import annotations

import ast
from pathlib import Path
import unittest

from recursive_horizons.fgc.evolution import pro20_rsrc1_attempt as isolated
from recursive_horizons.fgc.evolution.pro20_ev1_runtime import (
    execute_scheduled_attempt as historical_execute,
)
from tests.test_fgc_pro20_ev1_protocol import ControllableRHS
from tests.test_fgc_pro20_ev1_runtime import _checkpoint


ROOT = Path(__file__).resolve().parents[1]


def _summary(result) -> tuple[object, ...]:
    bundle = result.successor_bundle
    return (
        result.kind,
        result.member_key,
        result.accepted_state_advanced,
        result.executable_retry_cursor,
        result.physical_state_preserved,
        result.structural_class,
        result.terminal_evidence,
        bundle.descriptor,
        bundle.payload,
        bundle.cursor_bytes,
    )


class RSRC1AttemptParityTests(unittest.TestCase):
    def test_fresh_accepted_and_retry_relations_match_closed_runtime(self) -> None:
        for mode in ("accepted", "source", "cfl", "impulse"):
            with self.subTest(mode=mode):
                historical = _checkpoint()
                current = _checkpoint()
                if mode != "accepted":
                    for checkpoint in (historical, current):
                        rhs = checkpoint.member.operator
                        self.assertIsInstance(rhs, ControllableRHS)
                        rhs.impulse_from = checkpoint.member.time
                        rhs.mode = mode
                expected = historical_execute(historical)
                observed = isolated.execute_scheduled_attempt(current)
                self.assertEqual(_summary(observed), _summary(expected))

    def test_cap_and_pending_plan_match_closed_runtime(self) -> None:
        from recursive_horizons.fgc.evolution import pro20_ev1_runtime as historical

        left = _checkpoint()
        right = _checkpoint()
        self.assertEqual(
            isolated.fresh_requested_cap(left).hex(),
            historical.fresh_requested_cap(right).hex(),
        )
        self.assertIsNone(isolated.pending_plan_from_cursor(left.cursor))

    def test_resource_and_premise_stops_return_unchanged_bytes(self) -> None:
        predecessor = _checkpoint().generation_bundle()
        for key, evidence in (
            ("resource_stop", {"reason": "isolated peak"}),
            ("premise_stop", {"reason": "identity drift"}),
        ):
            result = isolated.classify_attempt(
                member_key="RK4-2049",
                predecessor_bundle=predecessor,
                successor_bundle=predecessor,
                accepted_state_advanced=False,
                **{key: evidence},
            )
            with self.subTest(key=key):
                self.assertEqual(result.successor_bundle.descriptor, predecessor.descriptor)
                self.assertEqual(result.successor_bundle.payload, predecessor.payload)
                self.assertEqual(result.successor_bundle.cursor_bytes, predecessor.cursor_bytes)
                self.assertFalse(result.accepted_state_advanced)
                self.assertFalse(result.executable_retry_cursor)

    def test_store_and_closed_runner_are_absent_from_the_module(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        joined = " ".join(imports)
        self.assertNotIn("pro20_ev1_runtime", joined)
        self.assertNotIn("pro20_ev1_store", joined)
        self.assertNotIn("campaign_store", joined)
        source = path.read_text(encoding="utf-8")
        self.assertNotIn("publish_", source)
        self.assertNotIn("write_bytes", source)


if __name__ == "__main__":
    unittest.main()
