from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.reference_principal_identity import (  # noqa: E402
    reference_principal_identity_certificate,
)


class ReferencePrincipalIdentityTests(unittest.TestCase):
    def test_universal_formal_index_identity(self) -> None:
        certificate = reference_principal_identity_certificate()
        self.assertTrue(
            certificate["connection_variation_route_equals_projector_route"]
        )
        self.assertTrue(
            certificate["all_sixteen_contravariant_tensor_differences_zero"]
        )
        self.assertTrue(
            certificate["reference_connection_absent_from_principal_variation"]
        )
        self.assertFalse(
            certificate["uses_fixtures_sampling_intervals_or_floating_point"]
        )
        self.assertEqual(len(certificate["canonical_identity_sha256"]), 64)
        self.assertTrue(
            all(
                count > 0
                for row in certificate["connection_route_term_counts"]
                for count in row
            )
        )

    def test_wrong_connection_sign_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "differs"):
            reference_principal_identity_certificate(
                connection_difference_sign=Q(1)
            )

    def test_nonexact_mutation_input_is_rejected(self) -> None:
        with self.assertRaisesRegex(TypeError, "Fraction"):
            reference_principal_identity_certificate(  # type: ignore[arg-type]
                connection_difference_sign=-1
            )


if __name__ == "__main__":
    unittest.main()
