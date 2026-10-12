"""Compact repository-checker coverage for sealed HLT16/PRO19 prelaunch."""
from __future__ import annotations

import unittest

from scripts import check_repo


class HLT16PrelaunchRepositoryCheckTests(unittest.TestCase):
    def test_mon16_compact_certificate_is_sealed(self) -> None:
        self.assertEqual(check_repo._check_hlt16_mon16_result(), [])

    def test_pro19_prelaunch_manifest_is_sealed(self) -> None:
        self.assertEqual(check_repo._check_pro19_prelaunch_result(), [])


if __name__ == "__main__":
    unittest.main()
