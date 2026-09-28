from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg4_sampling_identifiability import (
    TDG4_INSUFFICIENT_PREMISES,
    TDG4_SAMPLE_COUNT,
    TDG4_SUFFICIENT_ASSUMPTION_ROUTES,
    TDG4_TOP_BINS,
    finite_dimensional_sampling_assessment,
    generic_smooth_kernel_witness,
    lipschitz_coefficient_error_bound,
    sampling_identifiability_preflight,
    uniform_alias_witness,
)


class TDG4SamplingIdentifiabilityTests(unittest.TestCase):
    def test_smooth_between_sample_witness_is_universal_for_every_target_bin(self) -> None:
        nodes = tuple(Fraction(index, 63) for index in range(64))
        for target_bin in TDG4_TOP_BINS:
            witness = generic_smooth_kernel_witness(
                nodes, target_bin=target_bin
            )
            self.assertTrue(witness["all_sample_values_are_exactly_zero"])
            self.assertTrue(witness["all_derivatives_at_sample_nodes_are_zero"])
            self.assertTrue(
                witness["target_coefficient_nonzero_for_nonzero_amplitude"]
            )
            self.assertTrue(witness["sample_fibre_target_power_is_unbounded"])
            self.assertEqual(witness["support_interval"], ["1/189", "2/189"])

    def test_exact_uniform_alias_hits_each_declared_top_bin(self) -> None:
        for target_bin in TDG4_TOP_BINS:
            witness = uniform_alias_witness(
                sample_count=TDG4_SAMPLE_COUNT,
                target_bin=target_bin,
            )
            self.assertEqual(witness["carrier_frequency"], 63)
            self.assertEqual(witness["partner_bin"], 126 - target_bin)
            self.assertEqual(witness["target_complex_coefficient"], "1/4")
            self.assertEqual(
                witness["unit_amplitude_positive_bin_field_power"], "1/8"
            )
            self.assertTrue(witness["all_sample_values_are_exactly_zero"])
            self.assertTrue(
                witness["target_power_unbounded_under_amplitude_scaling"]
            )

    def test_exact_row_space_criterion_separates_factorable_and_kernel_target(self) -> None:
        nonidentifiable = finite_dimensional_sampling_assessment(
            ((1, 0, 0), (0, 1, 0)),
            ((0, 0, 1),),
        )
        self.assertFalse(nonidentifiable["target_annihilates_sampling_kernel"])
        self.assertTrue(
            nonidentifiable["unrestricted_sample_fibre_target_is_unbounded"]
        )

        factorable = finite_dimensional_sampling_assessment(
            ((1, 0, 0), (0, 1, 0)),
            ((2, 3, 0),),
        )
        self.assertTrue(factorable["target_annihilates_sampling_kernel"])
        self.assertTrue(factorable["target_factors_through_samples"])
        self.assertFalse(factorable["sample_map_is_injective"])

        injective = finite_dimensional_sampling_assessment(
            ((1, 0), (0, 1)),
            ((2, 3),),
        )
        self.assertTrue(injective["sample_map_is_injective"])
        self.assertTrue(injective["target_factors_through_samples"])

    def test_quantitative_derivative_bound_yields_exact_coefficient_enclosure(self) -> None:
        evidence = lipschitz_coefficient_error_bound(
            (Fraction(0), Fraction(1, 3), Fraction(1)),
            derivative_bound=Fraction(6),
        )
        self.assertEqual(evidence["sum_gap_squares"], "5/9")
        self.assertEqual(evidence["maximum_gap"], "2/3")
        self.assertEqual(evidence["coefficient_error_bound"], "5/3")
        self.assertEqual(evidence["coarser_max_gap_bound"], "2")
        self.assertTrue(
            evidence["derivative_bound_is_an_external_quantitative_premise"]
        )
        self.assertTrue(evidence["samples_alone_do_not_supply_the_derivative_bound"])

    def test_finite_densification_without_regular_control_remains_nonidentifying(self) -> None:
        for count in (64, 128, 1024):
            nodes = tuple(Fraction(index, count - 1) for index in range(count))
            witness = generic_smooth_kernel_witness(nodes, target_bin=29)
            self.assertEqual(witness["sample_count"], count)
            self.assertTrue(
                witness["finite_densification_alone_remains_nonidentifying"]
            )

    def test_preflight_keeps_every_physical_and_replacement_gate_closed(self) -> None:
        evidence = sampling_identifiability_preflight()
        self.assertTrue(evidence["all_exact_controls_pass"])
        self.assertEqual(evidence["sample_count"], 64)
        self.assertEqual(tuple(evidence["top_bins"]), TDG4_TOP_BINS)
        self.assertEqual(
            tuple(evidence["sufficient_assumption_routes"]),
            TDG4_SUFFICIENT_ASSUMPTION_ROUTES,
        )
        self.assertEqual(
            tuple(evidence["insufficient_premises"]), TDG4_INSUFFICIENT_PREMISES
        )
        self.assertFalse(evidence["actual_terminal_histories_consumed"])
        self.assertFalse(evidence["replacement_temporal_admission_defined"])
        self.assertFalse(evidence["PROTO14_frozen"])

    def test_invalid_nodes_bins_matrices_and_bounds_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            generic_smooth_kernel_witness((0, Fraction(1, 2), Fraction(1, 2), 1), target_bin=29)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            generic_smooth_kernel_witness((0.0, 1.0), target_bin=29)
        with self.assertRaisesRegex(ValueError, "positive"):
            generic_smooth_kernel_witness((0, 1), target_bin=0)
        with self.assertRaisesRegex(ValueError, "strictly below"):
            uniform_alias_witness(sample_count=64, target_bin=63)
        with self.assertRaisesRegex(ValueError, "domain dimension"):
            finite_dimensional_sampling_assessment(((1, 0),), ((1, 0, 0),))
        with self.assertRaisesRegex(ValueError, "spanning"):
            lipschitz_coefficient_error_bound(
                (Fraction(1, 10), Fraction(1)), derivative_bound=1
            )
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            lipschitz_coefficient_error_bound((0, 1), derivative_bound=-1)


if __name__ == "__main__":
    unittest.main()
