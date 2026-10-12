"""Focused no-advance tests for the PROTO19 GR-0 reconstruction bridge."""
from __future__ import annotations

import ast
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.proto19_progression_inputs import (
    AMPLITUDE, BRANCH, RESTART_TIME, TARGET_TIME,
)
from tests.fgc_gen0_fixture import reconstruct_sealed_gen0_members


class Proto19ProgressionInputTests(unittest.TestCase):
    def test_public_scope_is_hardcoded_gr0_event_23_to_24(self) -> None:
        self.assertEqual((BRANCH, AMPLITUDE, RESTART_TIME, TARGET_TIME), ("GR-0", "3", 23.0 / 16.0, 1.5))

    def test_source_forbids_candidate_branch_imports_and_mutation_apis(self) -> None:
        source = Path(__file__).parents[1] / "src/recursive_horizons/fgc/evolution/proto19_progression_inputs.py"
        tree = ast.parse(source.read_text("utf-8"))
        names = {node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)}
        self.assertFalse(any("fgcqr" in name.lower() or "sgbl" in name.lower() for name in names))
        text = source.read_text("utf-8")
        for forbidden in ("holdout", "materialize_generation_zero", "advance_to", "commit_tdg", "candidate_execution"):
            self.assertNotIn(forbidden, text)

    def test_actual_sealed_corpus_reconstructs_six_members_without_advance(self) -> None:
        result = reconstruct_sealed_gen0_members(ROOT)
        self.assertEqual(tuple(result.members), ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"))
        for member in result.members.values():
            self.assertEqual(member.time, RESTART_TIME)
            self.assertEqual(member.temporal_ledger.last_accepted_time, RESTART_TIME)
            self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 0)
            self.assertEqual(member.temporal_ledger.cumulative_temporal_retry_count, 0)
            self.assertEqual(member.temporal_ledger.accumulated_debit_vector, (0.0,) * 18)

    def test_legacy_raw_restore_is_never_called(self) -> None:
        from scripts import run_fgc_gr0_calibration_v13 as proto13
        with patch.object(proto13, "_restore_frozen_members", side_effect=AssertionError("legacy raw open")):
            result = reconstruct_sealed_gen0_members(ROOT)
        self.assertEqual(len(result.members), 6)


if __name__ == "__main__":
    unittest.main()
