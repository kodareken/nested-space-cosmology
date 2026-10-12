from pathlib import Path
import sys
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto16_transition_theorem import (
    RECEIPT_FIELDS,
    audit,
    receipt_fixture,
)


class TransitionTheoremTests(unittest.TestCase):
    def test_static_completeness_audit(self):
        with (ROOT / "configs/fgc/fgc-2-sf1-protocol-v16.toml").open("rb") as handle:
            value = audit(tomllib.load(handle))
        self.assertTrue(value["receipt"]["receipt_acyclic"])
        self.assertFalse(value["successor_schema"]["exact_derivation_frozen"])
        self.assertFalse(value["genesis_schema"]["exact_derivation_frozen"])
        self.assertTrue(
            value["conclusion"]["PROTO17_successor_freeze_design_authorized"]
        )

    def test_receipt_excludes_successors(self):
        record = receipt_fixture()
        fields = tuple(key for key in record["payload"] if key != "sequence")
        self.assertEqual(fields, RECEIPT_FIELDS)
        self.assertFalse(
            any("new_" in key or "successor" in key for key in record["payload"])
        )
