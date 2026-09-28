"""Focused attacks for the prospective TDG8 replay-evidence freeze."""
from __future__ import annotations

from pathlib import Path
import unittest
from unittest.mock import patch

from scripts import reproduce_fgc_tdg8_frz1 as freeze


ROOT = Path(__file__).resolve().parents[1]


class TDG8ReplayFreezeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.raw = (ROOT / freeze.CONFIG).read_bytes()

    def test_exact_freeze_reproduces_without_run_namespace_access(self) -> None:
        result = freeze.build_result(self.raw)
        self.assertEqual(result["artifact_id"], freeze.ARTIFACT)
        self.assertFalse(result["artifact_payload"]["claims"]["repair_applied"])
        self.assertFalse(
            result["artifact_payload"]["claims"]["successor_campaign_authorized"]
        )

    def test_predecessor_and_conditional_repair_mutations_fail(self) -> None:
        mutations = (
            (b'terminal_generation = 9', b'terminal_generation = 8'),
            (
                b'replacement_comparison = "p15._json_safe(seen[0]) != dict(evidence)"',
                b'replacement_comparison = "seen[0] != dict(evidence)"',
            ),
            (b'preserve_array_order = true', b'preserve_array_order = false'),
            (b'event_count = 1', b'event_count = 2'),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                with self.assertRaises(freeze.TDG8FreezeError):
                    freeze.build_result(self.raw.replace(old, new))

    def test_immutable_source_drift_fails_closed(self) -> None:
        real = freeze._git_show

        def changed(repository: Path, commit: str, path: str) -> bytes:
            raw = real(repository, commit, path)
            if path.endswith("hlt16_progression_attempt.py"):
                return raw + b"\n"
            return raw

        with patch.object(freeze, "_git_show", side_effect=changed):
            with self.assertRaises(freeze.TDG8FreezeError):
                freeze.build_result(self.raw)


if __name__ == "__main__":
    unittest.main()
