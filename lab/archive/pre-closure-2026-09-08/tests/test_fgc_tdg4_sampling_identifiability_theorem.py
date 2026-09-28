from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg4_sampling_identifiability_theorem import (
    TDG4_THEOREM_INSUFFICIENT_PREMISES,
    TDG4_THEOREM_SUFFICIENT_ROUTES,
    TDG4_THEOREM_TOP_BINS,
    factorization_control,
    factorization_lemma_certificate,
    sampling_identifiability_theorem_certificate,
    smooth_kernel_witness_theorem,
    uniform_alias_crosscheck,
)


class FGCTDG4SamplingIdentifiabilityTheoremTests(unittest.TestCase):
    def test_factorization_lemma_controls_both_directions(self) -> None:
        certificate = factorization_lemma_certificate()
        controls = certificate["exact_matrix_controls"]
        self.assertTrue(certificate["all_factorization_controls_pass"])
        self.assertTrue(controls["injective"]["sample_map_is_injective"])
        self.assertTrue(
            controls["factorable_noninjective"]["target_factors_through_samples"]
        )
        self.assertFalse(
            controls["nonidentifiable"]["target_annihilates_sampling_kernel"]
        )
        self.assertTrue(
            controls["nonidentifiable"][
                "unrestricted_kernel_amplitude_makes_target_unbounded"
            ]
        )

    def test_smooth_witness_handles_nonuniform_nodes_and_every_bin(self) -> None:
        nodes = (0, Fraction(1, 10), Fraction(2, 5), Fraction(9, 10), 1)
        for frequency in TDG4_THEOREM_TOP_BINS:
            witness = smooth_kernel_witness_theorem(
                nodes, target_bin=frequency
            )
            self.assertEqual(witness["selected_gap"], ["2/5", "9/10"])
            self.assertEqual(witness["support_interval"], ["17/30", "11/15"])
            self.assertTrue(
                witness["zero_on_a_neighborhood_of_every_sample"]
            )
            self.assertTrue(
                witness[
                    "target_coefficient_nonzero_for_every_nonzero_amplitude"
                ]
            )
            self.assertTrue(
                witness[
                    "both_declared_powers_unbounded_as_abs_A_tends_to_infinity"
                ]
            )

    def test_uniform_alias_is_exact_for_each_declared_bin(self) -> None:
        for frequency in TDG4_THEOREM_TOP_BINS:
            alias = uniform_alias_crosscheck(
                sample_count=64, target_bin=frequency
            )
            self.assertTrue(alias["all_sample_values_exactly_zero"])
            self.assertEqual(alias["carrier_frequency"], 63)
            self.assertEqual(alias["partner_bin"], 126 - frequency)
            self.assertEqual(alias["target_complex_coefficient"], "1/4")
            self.assertEqual(alias["partner_complex_coefficient"], "-1/4")
            self.assertEqual(alias["positive_bin_field_power"], "1/8")

    def test_any_finite_densification_retains_a_kernel_witness(self) -> None:
        for sample_count in (2, 3, 64, 129):
            nodes = tuple(
                Fraction(index, sample_count - 1)
                for index in range(sample_count)
            )
            certificate = sampling_identifiability_theorem_certificate(
                nodes, target_bins=(28,)
            )
            self.assertTrue(
                certificate["sampling_identifiability_theorem_completed"]
            )
            self.assertFalse(
                certificate["finite_samples_alone_bound_field_power"]
            )
            self.assertFalse(
                certificate["finite_samples_alone_bound_derivative_power"]
            )

    def test_complete_certificate_preserves_scope_and_assumptions(self) -> None:
        certificate = sampling_identifiability_theorem_certificate()
        self.assertTrue(
            certificate["sampling_identifiability_theorem_completed"]
        )
        self.assertEqual(certificate["sample_count"], 64)
        self.assertEqual(
            tuple(certificate["target_bins"]), TDG4_THEOREM_TOP_BINS
        )
        self.assertEqual(
            tuple(certificate["sufficient_assumption_routes"]),
            TDG4_THEOREM_SUFFICIENT_ROUTES,
        )
        self.assertEqual(
            tuple(certificate["insufficient_premises"]),
            TDG4_THEOREM_INSUFFICIENT_PREMISES,
        )
        self.assertFalse(
            certificate[
                "actual_unsampled_power_of_the_PROTO13_histories_inferred"
            ]
        )
        self.assertFalse(certificate["PDE_solution_manifold_restriction_proved"])
        self.assertFalse(certificate["actual_terminal_histories_consumed"])
        self.assertFalse(certificate["campaign_checkpoint_loaded"])
        self.assertFalse(certificate["state_advanced"])

    def test_invalid_inputs_fail_closed(self) -> None:
        invalid_calls = (
            lambda: smooth_kernel_witness_theorem((0, 1), target_bin=0),
            lambda: smooth_kernel_witness_theorem((0, 1, 1), target_bin=28),
            lambda: smooth_kernel_witness_theorem(
                (0, Fraction(3, 2)), target_bin=28
            ),
            lambda: uniform_alias_crosscheck(sample_count=2, target_bin=1),
            lambda: uniform_alias_crosscheck(sample_count=64, target_bin=63),
            lambda: factorization_control(((1, 0),), ((1, 0, 0),)),
            lambda: factorization_control(((1, 0), (1,)), ((1, 0),)),
        )
        for call in invalid_calls:
            with self.assertRaises((TypeError, ValueError)):
                call()


if __name__ == "__main__":
    unittest.main()
