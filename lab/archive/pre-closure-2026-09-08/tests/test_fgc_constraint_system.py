from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config  # noqa: E402
from recursive_horizons.fgc.constraint_system import (  # noqa: E402
    acceleration_independence_control,
    conditional_constraint_closure_statement,
    constraint_acceleration_structure_certificate,
    constraint_monitor_contract,
    gauge_constraint_characteristics,
    gauge_constraint_monitor_bundle,
    physical_constraint_term_certificate,
    physical_to_normal_gauge_map,
    reduction_constraint_monitor_bundle,
    relative_cancellation_monitor,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_reduction_subsidiary import (  # noqa: E402
    ReductionDifferentialFieldJet,
    ReductionDifferentialState,
)
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_gauge_constraint,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import Jet2, SphericalState  # noqa: E402


def _nonflat_state() -> SphericalState:
    return SphericalState(
        h_tt=Jet2(-Q(3, 2), Q(1, 13), -Q(1, 17), Q(1, 19), Q(1, 23), -Q(1, 29)),
        h_tr=Jet2(Q(1, 5), -Q(1, 31), Q(1, 37), Q(1, 41), -Q(1, 43), Q(1, 47)),
        h_rr=Jet2(Q(4, 3), Q(1, 53), -Q(1, 59), Q(1, 61), Q(1, 67), -Q(1, 71)),
        areal_radius=Jet2(Q(5, 2), -Q(1, 73), Q(1, 79), Q(1, 83), -Q(1, 89), Q(1, 97)),
        phi=Jet2(Q(1, 7), Q(1, 101), -Q(1, 103), Q(1, 107), Q(1, 109), -Q(1, 113)),
        chi=Jet2(-Q(1, 11), -Q(1, 127), Q(1, 131), Q(1, 137), -Q(1, 139), Q(1, 149)),
        planck_mass=Q(1),
        beta=Q(1, 3),
        mu=Q(1),
        g4=Q(1),
        eta=Q(1, 5),
        branch="FGC-QR",
    )


class ConstraintSystemTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        fixture = load_config().flat_fixture
        cls.datum = activated_compatible_state(fixture)
        cls.compatible = cls.datum["state"]
        cls.compatible_gauge = modified_harmonic_gauge_constraint(
            cls.compatible,
            reference=cls.datum["reference"],
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
        )
        cls.nonflat = _nonflat_state()
        cls.nonflat_gauge = modified_harmonic_gauge_constraint(
            cls.nonflat,
            reference=flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2)),
            coordinate_radius=Q(5, 2),
            tilde_normal_factor=Q(4),
        )

    def test_structural_and_exact_acceleration_independence(self) -> None:
        structure = constraint_acceleration_structure_certificate()
        self.assertTrue(structure["coordinate_time_accelerations_absent_structurally"])
        self.assertEqual(structure["Gauss_Bonnet_sector"]["Hamiltonian_support_count"], 36)
        self.assertTrue(
            structure["Gauss_Bonnet_sector"][
                "Hamiltonian_curvature_indices_all_spatial"
            ]
        )
        control = acceleration_independence_control([self.compatible, self.nonflat])
        self.assertEqual(control["state_count"], 2)
        self.assertTrue(control["all_exact_controls_passed"])
        self.assertTrue(
            all(
                direction["both_exactly_zero"]
                for record in control["records"]
                for direction in record["directions"]
            )
        )

    def test_physical_terms_recompose_and_monitor_nonzero_constraint(self) -> None:
        compatible = physical_constraint_term_certificate(self.compatible)
        self.assertTrue(compatible["Hamiltonian_exact_recomposition"])
        self.assertTrue(compatible["momentum_exact_recomposition"])
        self.assertEqual(compatible["Hamiltonian"], 0)
        nonflat = physical_constraint_term_certificate(self.nonflat)
        self.assertTrue(nonflat["Hamiltonian_exact_recomposition"])
        self.assertTrue(nonflat["momentum_exact_recomposition"])
        self.assertTrue(
            nonflat["Hamiltonian"] != 0 or nonflat["momentum"] != 0
        )
        self.assertLessEqual(
            nonflat["Hamiltonian_monitor"]["relative_cancellation"], Q(1)
        )

    def test_general_physical_to_normal_gauge_map(self) -> None:
        for state, gauge in (
            (self.compatible, self.compatible_gauge),
            (self.nonflat, self.nonflat_gauge),
        ):
            result = physical_to_normal_gauge_map(
                state, gauge, hat_normal_factor=Q(9)
            )
            self.assertEqual(result["computed_matrix"], result["analytic_matrix"])
            self.assertTrue(result["determinant_strictly_negative"])
            self.assertTrue(
                result["normal_derivative_equivalence"][
                    "H_and_M_zero_iff_normal_nabla_C_zero_on_REF1_shell"
                ]
            )

        changed = replace(
            self.nonflat,
            h_tr=replace(self.nonflat.h_tr, value=Q(2, 5)),
        )
        gauge = modified_harmonic_gauge_constraint(
            changed,
            reference=flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2)),
            coordinate_radius=Q(5, 2),
            tilde_normal_factor=Q(4),
        )
        result = physical_to_normal_gauge_map(changed, gauge, hat_normal_factor=Q(9))
        self.assertNotEqual(result["analytic_matrix"][1][0], 0)

    def test_gauge_and_reduction_characteristics_and_monitors(self) -> None:
        characteristics = gauge_constraint_characteristics(
            self.compatible, hat_normal_factor=Q(9)
        )
        self.assertEqual(
            tuple(root["exact"] for root in characteristics["root_isolations"]),
            (-Q(1, 3), Q(1, 3)),
        )
        self.assertEqual(characteristics["outer_incoming_component_count_at_control"], 2)
        self.assertEqual(characteristics["inner_incoming_component_count_at_control"], 2)
        self.assertFalse(characteristics["constraint_preserving_boundary_map_to_main_fields_derived"])

        gauge_monitors = gauge_constraint_monitor_bundle(self.nonflat_gauge)
        self.assertTrue(gauge_monitors["all_components_recomposed_exactly"])

        jet = ReductionDifferentialFieldJet(
            u_t=Q(2), p=Q(2), q=Q(3), u_r=Q(3),
            q_t=Q(5), p_r=Q(5), u_tr=Q(5),
        )
        reduction = ReductionDifferentialState(
            **{field: jet for field in ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi")}
        )
        monitors = reduction_constraint_monitor_bundle(reduction)
        self.assertEqual(monitors["subsidiary_radial_principal_speed"], 0)
        self.assertTrue(
            all(
                record["subsidiary_identity_d_t_C_plus_d_r_D_minus_K"]["raw_value"] == 0
                for record in monitors["records"]
            )
        )

    def test_monitor_and_closure_nonclaim_boundaries(self) -> None:
        zero = relative_cancellation_monitor((Q(0), Q(0)), expected_value=Q(0))
        self.assertTrue(zero["zero_scale_convention_applied"])
        cancellation = relative_cancellation_monitor((Q(3), -Q(2)))
        self.assertEqual(cancellation["relative_cancellation"], Q(1, 5))
        with self.assertRaisesRegex(ValueError, "reconstruct"):
            relative_cancellation_monitor((Q(1), Q(2)), expected_value=Q(4))

        contract = constraint_monitor_contract()
        self.assertTrue(contract["smallness_on_one_run_is_not_a_propagation_proof"])
        statement = conditional_constraint_closure_statement()
        self.assertTrue(
            statement["complete_boundary_free_spherical_constraint_system_closed_conditionally"]
        )
        for open_item in (
            "smooth_solution_existence_supplied",
            "regular_center_supplied",
            "compatible_initial_hypersurface_supplied",
            "constraint_preserving_boundary_map_supplied",
            "initial_boundary_value_problem_supplied",
            "evolution_authorized",
        ):
            self.assertFalse(statement[open_item])


if __name__ == "__main__":
    unittest.main()
