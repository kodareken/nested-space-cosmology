"""Focused prospective authority controls for TDG9 LOC2."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

from recursive_horizons.fgc.evolution import tdg9_loc2_authority as authority  # noqa: E402


class LOC2AuthorityTests(unittest.TestCase):
    def test_config_audit_budget_and_nonclaims_are_exact(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        self.assertEqual(parsed["coefficient_depth_audit"], authority.COEFFICIENT_DEPTH_AUDIT)
        self.assertEqual(parsed["coefficient_depth_audit"]["maxima"]["proof_bits"], 200)  # type: ignore[index]
        self.assertEqual(parsed["implementation"]["global_proof_bit_ceiling"], 9216)  # type: ignore[index]
        self.assertEqual(parsed["work_budget"]["component_cubic_count"], 367_920)  # type: ignore[index]
        self.assertFalse(parsed["scope"]["candidate_branches_authorized"])  # type: ignore[index]

    def test_absence_rejects_existing_and_broken_symlink_parents(self) -> None:
        for broken in (False, True):
            with self.subTest(broken=broken), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "runs").mkdir()
                destination = root / "missing" if broken else root / "foreign"
                if not broken:
                    destination.mkdir()
                os.symlink(destination, root / "runs" / "fgc-2-sf1")
                with self.assertRaisesRegex(authority.LOC2AuthorityError, "output parent is unsafe"):
                    authority._require_absent(root)  # type: ignore[attr-defined]

    def test_compact_result_rejects_claim_mutation(self) -> None:
        raw = (ROOT / authority.RESULT_PATH).read_bytes()
        value = authority.validate_compact((ROOT / authority.CONFIG_PATH).read_bytes(), raw)
        value["artifact_payload"]["claims"]["physical_result_earned"] = True  # type: ignore[index]
        with self.assertRaises(authority.LOC2AuthorityError):
            authority.validate_compact(
                (ROOT / authority.CONFIG_PATH).read_bytes(),
                authority.canonical_pretty(value),
            )

    def test_committed_binding_mutation_fails_closed(self) -> None:
        with (
            patch.object(authority, "_git", return_value=b"mutated"),
            self.assertRaisesRegex(
                authority.LOC2AuthorityError, "committed LOC2 binding differs"
            ),
        ):
            authority._verify_committed_bindings(ROOT, "f" * 40)  # type: ignore[attr-defined]

    def test_authorization_absence_gate_rejects_existing_namespace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.LOC2AuthorityError, "output namespace already exists"
            ):
                authority._require_absent(root)  # type: ignore[attr-defined]


if __name__ == "__main__":
    unittest.main()
