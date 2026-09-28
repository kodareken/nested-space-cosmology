from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.reproduce_fgc_hyp1_reduction import load_config  # noqa: E402
from scripts.reproduce_fgc_hyp1_mhg_implicit import load_config as load_imp1_config  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_propagation import (  # noqa: E402
    GaugeConstraintJet2, VectorJet2, covariant_vector_derivatives,
    gauge_propagation_principal_blocks, reference_gauge_propagation_operator,
    modified_harmonic_propagation_certificate, required_prop1_nonclaims,
)
import recursive_horizons.fgc.modified_harmonic_propagation as propagation_module  # noqa: E402
from recursive_horizons.fgc.reference_connection import flat_spherical_annulus_reference  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_reference import physical_connection_data  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import state_from_generalized_adm_pg_fixture  # noqa: E402


class PropagationOperatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures = load_config().configuration["fixtures"]
        cls.reference = flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2))

    def probe(self):
        return GaugeConstraintJet2(VectorJet2(Q(2, 3), Q(-1, 5), Q(3, 7), Q(1, 9), Q(-1, 8), Q(2, 11)), VectorJet2(Q(-1, 4), Q(2, 9), Q(-3, 10), Q(1, 12), Q(2, 13), Q(-1, 14)))

    def test_strict_input_and_spherical_components(self):
        with self.assertRaises(ValueError):
            GaugeConstraintJet2.from_mapping({"C^t": {"value": 0}})
        with self.assertRaises(TypeError):
            VectorJet2(0.5)  # type: ignore[arg-type]
        state = state_from_generalized_adm_pg_fixture(self.fixtures[0])
        derivative = covariant_vector_derivatives(GaugeConstraintJet2(VectorJet2(0), VectorJet2(0)), physical_connection_data(state))
        self.assertTrue(all(value == 0 for row in derivative["covariant_first"] for value in row))

    def test_flat_zero_and_nonzero_controls(self):
        state = state_from_generalized_adm_pg_fixture(load_imp1_config().fixture)
        zero = reference_gauge_propagation_operator(state, GaugeConstraintJet2(VectorJet2(0), VectorJet2(0)), reference=self.reference, coordinate_radius=Q(4), hat_normal_factor=Q(9))
        self.assertTrue(zero["two_routes_agree_exactly"])
        self.assertEqual(zero["coordinate_divergence"], (Q(0),) * 4)
        nonzero = reference_gauge_propagation_operator(state, self.probe(), reference=self.reference, coordinate_radius=Q(4), hat_normal_factor=Q(9))
        self.assertTrue(nonzero["two_routes_agree_exactly"])
        self.assertNotEqual(nonzero["coordinate_divergence"], (Q(0),) * 4)

    def test_activated_lower_order_paths_and_nonclaims(self):
        state = state_from_generalized_adm_pg_fixture(self.fixtures[3])
        data = reference_gauge_propagation_operator(state, self.probe(), reference=self.reference, coordinate_radius=Q(9, 2), hat_normal_factor=Q(9))
        self.assertTrue(data["two_routes_agree_exactly"])
        self.assertNotEqual(data["expanded_gradient_F_term"], (Q(0),) * 4)
        self.assertNotEqual(data["expanded_projector_derivative_term"], (Q(0),) * 4)
        self.assertNotEqual(data["expanded_curvature_commutator_term"], (Q(0),) * 4)
        self.assertNotEqual(data["expanded_connection_wave_contribution"], (Q(0),) * 4)
        self.assertTrue(data["gauge_input_is_independent"])
        self.assertFalse(data["direct_metric_derived_residual_divergence_evaluated"])
        self.assertTrue(all(v is False for v in required_prop1_nonclaims().values()))

    def test_all_spherical_principal_blocks_regress_to_mhg1(self):
        state = state_from_generalized_adm_pg_fixture(self.fixtures[3])
        blocks = gauge_propagation_principal_blocks(
            state, reference=self.reference, coordinate_radius=Q(9, 2),
            hat_normal_factor=Q(9),
        )
        self.assertTrue(blocks["spherical_t_r_columns_match_mhg1_exactly"])
        self.assertTrue(blocks["angular_independent_C_columns_not_represented"])
        self.assertEqual(
            set(blocks["C_value_first_second_partial_blocks"]),
            {"value", "dt", "dr", "dtt", "dtr", "drr"},
        )
        self.assertEqual(blocks["coefficient_block_shape"], (4, 2))
        self.assertEqual(blocks["represented_C_input_order"], ("C^t", "C^r"))
        self.assertTrue(
            all(
                len(block) == 4 and all(len(row) == 2 for row in block)
                for block in blocks["C_value_first_second_partial_blocks"].values()
            )
        )

    def test_route_mutation_is_detected_by_independent_derivative_assemblies(self):
        state = state_from_generalized_adm_pg_fixture(self.fixtures[3])
        original = propagation_module.covariant_vector_derivatives

        def corrupt_expanded_covariant_second(gauge, physical):
            data = dict(original(gauge, physical))
            second = [
                [list(component_row) for component_row in direction_rows]
                for direction_rows in data["covariant_second"]
            ]
            second[0][0][0] += Q(1)
            data["covariant_second"] = tuple(
                tuple(tuple(row) for row in direction_rows)
                for direction_rows in second
            )
            return data

        with patch.object(
            propagation_module,
            "covariant_vector_derivatives",
            side_effect=corrupt_expanded_covariant_second,
        ):
            mutated = reference_gauge_propagation_operator(
                state,
                self.probe(),
                reference=self.reference,
                coordinate_radius=Q(9, 2),
                hat_normal_factor=Q(9),
            )
        self.assertFalse(mutated["two_routes_agree_exactly"])
        self.assertNotEqual(
            mutated["coordinate_divergence"], mutated["expanded_divergence"]
        )

    def test_complete_prop1_certificate_is_exact_and_fail_closed(self):
        reduction = load_config()
        implicit = load_imp1_config()
        activated = next(
            fixture
            for fixture in reduction.configuration["fixtures"]
            if fixture["fixture_id"] == "FGCQR_activated_generic"
        )
        certificate = modified_harmonic_propagation_certificate(
            {
                "flat_fixture": implicit.fixture,
                "activated_fixture": activated,
                "flat_probe": GaugeConstraintJet2(
                    VectorJet2(Q(1, 3), Q(2, 5), Q(-3, 7), Q(5, 11), Q(-7, 13), Q(11, 17)),
                    VectorJet2(Q(-2, 3), Q(3, 5), Q(4, 7), Q(-6, 11), Q(8, 13), Q(-10, 17)),
                ),
                "activated_probe": GaugeConstraintJet2(
                    VectorJet2(Q(2, 7), Q(-3, 8), Q(5, 9), Q(-7, 10), Q(11, 12), Q(-13, 14)),
                    VectorJet2(Q(-3, 7), Q(4, 9), Q(-5, 11), Q(8, 13), Q(-9, 14), Q(15, 16)),
                ),
                "reference": self.reference,
                "flat_coordinate_radius": Q(4),
                "activated_coordinate_radius": Q(9, 2),
                "tilde_normal_factor": Q(4),
                "hat_normal_factor": Q(9),
            }
        )
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            certificate["conditional_lower_order_reference_gauge_operator_derived"]
        )
        self.assertTrue(all(certificate["verified_exact_checks"].values()))
        self.assertEqual(
            certificate["formulation"]["output_order"],
            ("nu=t", "nu=r", "nu=theta", "nu=phi"),
        )
        self.assertTrue(all(value is False for value in certificate["nonclaims"].values()))

        bad = {
            "flat_fixture": implicit.fixture,
            "activated_fixture": activated,
        }
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            modified_harmonic_propagation_certificate(bad)


if __name__ == "__main__":
    unittest.main()
