from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.modified_harmonic_auxiliary_identity import (
    auxiliary_sector_identity_certificate,
    hat_projector_contraction_identity,
    noether_extension_composition,
    radial_hat_omitted_row_identity,
    tilde_pure_gauge_contraction_identity,
    tilde_right_kernel_composition,
)


class AuxiliaryIdentityTests(unittest.TestCase):
    def test_conditional_hat_completion_identity(self) -> None:
        certificate = auxiliary_sector_identity_certificate()
        self.assertTrue(certificate["hat_row_completion"]["all_six_radial_rows_zero"])
        self.assertTrue(certificate["hat_row_completion"]["raised_covector_chart_components_nonzero"])
        self.assertFalse(certificate["nonclaims"]["uniform_eigenframe_proven"])
        self.assertTrue(certificate["noether_extension_composition"]["action_diffeomorphism_invariance_is_named_premise_not_machine_derived"])
        self.assertTrue(certificate["canonical_chi_principal_identity"]["principal_dependency_ledger"]["metric_and_phi_rows_have_no_chi_second_jet_dependence"])
        self.assertTrue(certificate["tilde_pure_gauge_right_kernel"]["full_symbol_right_kernels_at_tilde_null"])

    def test_chart_and_hat_null_mutations_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonzero raised"):
            radial_hat_omitted_row_identity(xi_t_up=Q(0))
        with self.assertRaisesRegex(ValueError, "hat-null"):
            radial_hat_omitted_row_identity(hat_null=Q(1))
        with self.assertRaisesRegex(ValueError, "selected"):
            radial_hat_omitted_row_identity(selected_metric_tr=Q(1))
        with self.assertRaisesRegex(ValueError, "contraction differs"):
            hat_projector_contraction_identity(mutate_trace_sign=1)
        with self.assertRaisesRegex(ValueError, "diffeomorphism"):
            noether_extension_composition(action_diffeomorphism_invariance=False)
        with self.assertRaisesRegex(ValueError, "tilde pure-gauge contraction"):
            tilde_pure_gauge_contraction_identity(mutate_connection_sign=1)
        with self.assertRaisesRegex(ValueError, "diffeomorphism"):
            tilde_right_kernel_composition(action_diffeomorphism_invariance=False)


if __name__ == "__main__":
    unittest.main()
