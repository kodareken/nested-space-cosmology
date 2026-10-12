from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.reference_branch_principal_identity import (
    abstract_branch_block_identity,
    reference_branch_principal_identity_certificate,
)


class ReferenceBranchPrincipalIdentityTests(unittest.TestCase):
    def test_universal_structural_composition(self) -> None:
        certificate = reference_branch_principal_identity_certificate()
        self.assertTrue(certificate["live_definition_route"]["ref1_full_rows_call_unredefined_residuals"])
        self.assertTrue(certificate["live_definition_route"]["red1_second_jet_block_differentiates_residuals"])
        self.assertTrue(certificate["gauge_extension"]["universal_formal_certificate"]["connection_variation_route_equals_projector_route"])
        self.assertTrue(certificate["gauge_extension"]["scalar_rows_unmodified_by_definition"])
        self.assertTrue(certificate["branch_reduction"]["abstract_block_identity_exact"])
        self.assertFalse(certificate["theorem_scope"]["characteristics_or_hyperbolicity_proven"])

    def test_sign_and_orientation_mutations_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "differ"):
            abstract_branch_block_identity(mutate_coefficient_sign=True)
        with self.assertRaisesRegex(ValueError, "differ"):
            abstract_branch_block_identity(mutate_lower_sign=True)


if __name__ == "__main__":
    unittest.main()
