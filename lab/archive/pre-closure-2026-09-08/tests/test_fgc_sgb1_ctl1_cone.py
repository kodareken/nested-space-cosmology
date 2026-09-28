from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    esf_background,
    principal_coefficient_tensors,
)
from recursive_horizons.fgc.exact_interval import Interval, interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_cone import (  # noqa: E402
    ALGEBRAIC_CURVATURE_GENERATOR_COUNT,
    DEFAULT_ARC_COUNT,
    FORBIDDEN_HEALTH_IMPORTS,
    HESSIAN_GENERATOR_COUNT,
    KULKARNI_NOMIZU_CANDIDATE_COUNT,
    ORTHONORMAL_MINKOWSKI_CHART,
    REGRESSION_DIRECTIONS,
    SGBLConeLimits,
    SGBLConeStop,
    SGBLPrincipalBackgroundBox,
    SYMMETRIC_INDEX_PAIRS,
    _contains_float_matrix,
    _interval_inverse,
    sgbl_angular_regression_imaginary_parts,
    sgbl_enclose_principal_coefficients,
    sgbl_enclose_principal_cone,
    sgbl_esf_exact_reference,
    sgbl_esf_principal_box,
    sgbl_nonflat_source_principal_box,
    sgbl_principal_box_from_state,
    sgbl_pure_gauge_identity_defect,
    sgbl_structural_generator_basis,
    sgbl_structural_identities,
)
from recursive_horizons.fgc.sgb1_ctl1_principal import (  # noqa: E402
    sgbl_covariant_principal_background,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    SGBLSourceInputs,
    sgbl_source_solve,
    sgbl_source_state,
)


def _fixture_a() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(2, Q(1, 3), Q(-1, 2), Q(1, 5), Q(-1, 7)),
        shift=(Q(1, 3), Q(-2, 5), Q(1, 4), Q(-1, 6), Q(2, 9)),
        radial_metric=(Q(3, 2), Q(1, 8), Q(-1, 3), Q(1, 7), Q(-1, 4)),
        areal_radius=(4, Q(-1, 5), Q(3, 2), Q(1, 9), Q(-2, 5)),
        phi=(Q(1, 2), Q(2, 3), Q(-4, 5), Q(1, 4), Q(-1, 8)),
        chi=(Q(-3, 4), Q(1, 2), Q(3, 5), Q(-2, 7), Q(1, 6)),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


class SGBLConeInputTests(unittest.TestCase):
    def test_chart_sign_and_hyp2_sample_count_are_refused(self) -> None:
        with self.assertRaises(SGBLConeStop) as chart:
            SGBLPrincipalBackgroundBox(
                planck_mass=2, alpha_gb=Q(-1, 4), chart="metric_dtt"
            )
        self.assertEqual(chart.exception.reason, "cone_chart_error")
        with self.assertRaises(SGBLConeStop) as sign:
            SGBLPrincipalBackgroundBox(
                planck_mass=interval(-1, 1), alpha_gb=Q(-1, 4)
            )
        self.assertEqual(sign.exception.reason, "sign_chart_error")
        with self.assertRaises(SGBLConeStop) as copied:
            SGBLConeLimits(arc_count=512)
        self.assertEqual(copied.exception.reason, "resource_limit")
        with self.assertRaises(ValueError):
            SGBLPrincipalBackgroundBox(planck_mass=2, alpha_gb=0)
        control = sgbl_esf_principal_box(
            planck_mass=2, alpha_gb=0, allow_zero_coupling_control=True
        )
        self.assertEqual(control.gb_coupling_prime, Interval.singleton(0))
        self.assertEqual(control.effective_planck_prime, Interval.singleton(0))
        self.assertEqual(control.chart, ORTHONORMAL_MINKOWSKI_CHART)
        self.assertNotEqual(DEFAULT_ARC_COUNT, 512)

    def test_width_cap_is_a_resource_stop_not_a_fitted_floor(self) -> None:
        with self.assertRaises(SGBLConeStop) as stopped:
            sgbl_esf_principal_box(
                planck_mass=2,
                alpha_gb=Q(-1, 4),
                riemann_half_width=2,
                limits=SGBLConeLimits(max_parameter_half_width=1),
            )
        self.assertEqual(stopped.exception.reason, "resource_limit")


class SGBLConeGeneratorBasisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.basis = sgbl_structural_generator_basis()

    def test_complete_exact_basis_has_rank_ten_and_twenty(self) -> None:
        basis = self.basis
        self.assertEqual(len(SYMMETRIC_INDEX_PAIRS), 10)
        self.assertEqual(HESSIAN_GENERATOR_COUNT, 10)
        self.assertEqual(ALGEBRAIC_CURVATURE_GENERATOR_COUNT, 20)
        self.assertEqual(KULKARNI_NOMIZU_CANDIDATE_COUNT, 55)
        self.assertEqual(basis["hessian_generator_count"], 10)
        self.assertEqual(basis["hessian_generator_rank"], 10)
        self.assertEqual(basis["hessian_labels"], SYMMETRIC_INDEX_PAIRS)
        self.assertEqual(len(basis["hessian_generators"]), 10)
        self.assertTrue(basis["hessian_symmetric"])
        self.assertEqual(basis["algebraic_curvature_generator_count"], 20)
        self.assertEqual(basis["algebraic_curvature_generator_rank"], 20)
        self.assertEqual(len(basis["algebraic_curvature_generators"]), 20)
        self.assertEqual(len(basis["algebraic_curvature_labels"]), 20)
        self.assertEqual(basis["kulkarni_nomizu_candidate_count"], 55)
        self.assertEqual(basis["kulkarni_nomizu_span_rank"], 20)
        self.assertTrue(basis["pair_antisymmetric"])
        self.assertTrue(basis["pair_exchange"])
        self.assertTrue(basis["first_bianchi"])
        self.assertTrue(basis["candidates_all_algebraic"])
        self.assertTrue(basis["representative_r0101_in_span"])
        self.assertIs(basis["representative_generator_sufficient"], False)
        self.assertTrue(basis["complete_exact_basis"])
        self.assertEqual(basis["hessian_labels"][0], (0, 0))

    def test_mutation_or_omission_of_the_basis_fails_closed(self) -> None:
        with self.assertRaises(SGBLConeStop) as omitted_hess:
            sgbl_structural_generator_basis(omit_hessian_index=0)
        self.assertEqual(omitted_hess.exception.reason, "broken_gauge")
        self.assertEqual(omitted_hess.exception.payload["hessian_generator_count"], 9)
        with self.assertRaises(SGBLConeStop) as omitted_curv:
            sgbl_structural_generator_basis(omit_curvature_index=3)
        self.assertEqual(omitted_curv.exception.reason, "broken_gauge")
        self.assertEqual(
            omitted_curv.exception.payload["algebraic_curvature_generator_count"], 19
        )
        with self.assertRaises(SGBLConeStop) as bianchi:
            sgbl_structural_generator_basis(mutate_curvature_bianchi=True)
        self.assertEqual(bianchi.exception.reason, "broken_gauge")
        self.assertFalse(bianchi.exception.payload["first_bianchi"])
        with self.assertRaises(SGBLConeStop) as hessian:
            sgbl_structural_generator_basis(mutate_hessian_symmetry=True)
        self.assertEqual(hessian.exception.reason, "broken_gauge")
        self.assertFalse(hessian.exception.payload["hessian_symmetric"])

    def test_formula_identities_hold_on_every_generator(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        payload = sgbl_structural_identities(box)
        self.assertTrue(payload["structural_hessian_symmetric"])
        self.assertTrue(payload["structural_pure_gauge_identity"])
        self.assertTrue(payload["checked_esf_generator"])
        self.assertEqual(payload["checked_hessian_generators"], 10)
        self.assertEqual(payload["checked_curvature_generators"], 20)
        self.assertEqual(payload["hessian_generator_rank"], 10)
        self.assertEqual(payload["algebraic_curvature_generator_rank"], 20)
        self.assertTrue(payload["pair_antisymmetric"])
        self.assertTrue(payload["pair_exchange"])
        self.assertTrue(payload["first_bianchi"])
        self.assertIs(payload["representative_generator_sufficient"], False)
        self.assertTrue(box.curvature_free)
        self.assertEqual(payload["checked_curvature_generators"], 20)


class SGBLConeExactReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        cls.record = sgbl_enclose_principal_cone(cls.box)

    def test_esf_singleton_proves_the_continuum_cone_and_keeps_aggregate_health_closed(self) -> None:
        record = self.record
        self.assertEqual(record.classification, "continuum_cone_proved")
        self.assertTrue(record.real_complete_basis)
        self.assertTrue(record.positive_symmetrizer)
        self.assertTrue(record.strongly_hyperbolic)
        self.assertTrue(record.continuum_direction_bound_independent_of_n)
        self.assertGreater(record.kinetic_singular_value_lower, 0)
        self.assertEqual(record.physical_energy_coercivity_lower, 2)
        self.assertEqual(record.companion_deformation_upper, 0)
        self.assertTrue(
            all(margin > 0 for margin in record.cluster_resolvent_margins.values())
        )
        self.assertEqual(record.exact_eigenframe_rank, 24)
        self.assertTrue(record.structural_hessian_symmetric)
        self.assertTrue(record.structural_pure_gauge_identity)
        self.assertIs(record.representative_generator_sufficient, False)
        self.assertEqual(
            sgbl_structural_identities(self.box)["checked_hessian_generators"], 10
        )
        self.assertEqual(
            sgbl_structural_identities(self.box)["checked_curvature_generators"], 20
        )
        self.assertEqual(record.gauge_identity_defect, Interval.singleton(0))
        self.assertEqual(record.action_hessian_symmetry_defect, Interval.singleton(0))
        self.assertEqual(record.box.effective_planck_prime, Interval.singleton(0))
        self.assertEqual(record.action["beta"], 0)
        self.assertEqual(record.action["eta"], 0)
        self.assertEqual(record.action["F_prime"], 0)
        self.assertEqual(
            record.esf_reference["classification"],
            "exact_unnormalized_radial_eigenframe_and_Bauer_Fike_resolvent",
        )
        self.assertFalse(record.esf_reference["proof_uses_floating_svd"])
        self.assertTrue(record.esf_reference["arc_samples_are_regression_only"])
        self.assertTrue(record.esf_reference["ulp_allowance_not_used"])
        self.assertTrue(record.esf_reference["hyp2_sample_count_not_used"])
        self.assertEqual(record.esf_reference["eigenframe_rank"], 24)
        for cluster in record.esf_reference["clusters"].values():
            self.assertEqual(cluster["exact_nullity"], 4)
            self.assertEqual(cluster["eigenvalue_count"], 4)
        self.assertEqual(
            record.esf_reference["physical_energy"]["physical_plus"][
                "sylvester_leading_minors"
            ],
            (Q(4), Q(16), Q(32), Q(64)),
        )
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(
            record.quantitative_all_covector_weak_coupling_health_envelope_passed,
            False,
        )
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(record, representative_generator_sufficient=True)
        with self.assertRaises(ValueError):
            replace(record, real_complete_basis=False)

    def test_coefficient_enclosure_is_outward_and_grounded_in_sgbl(self) -> None:
        enclosure = self.record.coefficient_enclosure
        numpy_tensors = principal_coefficient_tensors(esf_background(2.0))
        self.assertEqual(enclosure.time_time[0][0], Interval.singleton(-18))
        self.assertEqual(enclosure.ungauged_time_time[10][10], Interval.singleton(-1))
        self.assertEqual(enclosure.ungauged_time_time[11][11], Interval.singleton(-1))
        mid = np.array(
            [[float(entry.midpoint()) for entry in row] for row in enclosure.time_time]
        )
        np.testing.assert_allclose(mid, numpy_tensors.time_time, rtol=0.0, atol=1.0e-10)
        hess = self.box.hessian_gb_lower
        alpha = self.box.alpha_gb
        for first in range(4):
            for second in range(4):
                self.assertEqual(
                    hess[first][second],
                    alpha * self.box.hessian_phi_lower[first][second],
                )
        self.assertTrue(self.box.curvature_free)
        self.assertTrue(self.box.is_singleton)

    def test_angular_regression_is_weaker_than_the_continuum_bound(self) -> None:
        imag = sgbl_angular_regression_imaginary_parts(self.box)
        self.assertEqual(tuple(imag), REGRESSION_DIRECTIONS)
        for direction, value in imag.items():
            self.assertLess(value, 1.0e-10)
            self.assertIn(direction, REGRESSION_DIRECTIONS)
        self.assertEqual(
            self.record.angular_regression_directions, REGRESSION_DIRECTIONS
        )

    def test_exact_reference_has_rank_24_frame_and_positive_sylvester_energy(self) -> None:
        reference = sgbl_esf_exact_reference(planck_mass=2)
        self.assertEqual(reference["eigenframe_rank"], 24)
        self.assertFalse(reference["proof_uses_floating_svd"])
        self.assertEqual(reference["physical_energy_coercivity_lower"], 2)
        self.assertEqual(reference["physical_energy_operator_norm_upper"], 4)
        self.assertGreater(reference["condition_number_2_upper"], 1)
        nullities = [
            cluster["exact_nullity"]
            for cluster in reference["clusters"].values()
        ]
        self.assertEqual(nullities, [4, 4, 4, 4, 4, 4])
        self.assertTrue(
            all(cluster["certified_lower_bound"] > 0 for cluster in reference["clusters"].values())
        )


class SGBLConeInconclusiveAndStopTests(unittest.TestCase):
    def test_wide_planck_box_is_inconclusive_not_a_kinetic_failure(self) -> None:
        record = sgbl_enclose_principal_cone(
            sgbl_esf_principal_box(
                planck_mass=interval(1, 3),
                alpha_gb=Q(-1, 4),
            )
        )
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertEqual(record.inconclusive_reason, "neumann_rho_not_below_one")
        self.assertFalse(record.real_complete_basis)
        self.assertFalse(record.positive_symmetrizer)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        with self.assertRaisesRegex(ValueError, "inconclusive cone"):
            replace(record, real_complete_basis=True, positive_symmetrizer=True)

    def test_resource_limit_stops_before_a_fitted_margin(self) -> None:
        box = sgbl_esf_principal_box(
            planck_mass=2,
            alpha_gb=Q(-1, 4),
            limits=SGBLConeLimits(max_residual_evaluations=1),
        )
        with self.assertRaises(SGBLConeStop) as stopped:
            sgbl_enclose_principal_coefficients(box)
        self.assertEqual(stopped.exception.reason, "resource_limit")

    def test_lost_kinetic_is_a_singular_center_not_wrapping(self) -> None:
        zero = tuple(tuple(Interval.singleton(0) for _ in range(12)) for _ in range(12))
        with self.assertRaises(SGBLConeStop) as stopped:
            _interval_inverse(zero)
        self.assertEqual(stopped.exception.reason, "lost_kinetic")

    def test_broken_gauge_mutation_is_strictly_positive(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        true_defect = sgbl_pure_gauge_identity_defect(box)
        mutated = sgbl_pure_gauge_identity_defect(box, mutate_sign=-1)
        self.assertEqual(true_defect, Interval.singleton(0))
        self.assertGreater(mutated.abs_upper(), 0)
        self.assertFalse(mutated.contains_zero() and mutated.abs_upper() == 0)

    def test_nonflat_source_fixture_contains_numpy_adapter_and_is_exactly_obstructed(self) -> None:
        box = sgbl_nonflat_source_principal_box()
        self.assertFalse(box.curvature_free)
        self.assertTrue(box.is_singleton)
        enclosure = sgbl_enclose_principal_coefficients(box)
        state = sgbl_source_state(_fixture_a(), sgbl_source_solve(_fixture_a()))
        numpy_tensors = principal_coefficient_tensors(
            sgbl_covariant_principal_background(state)
        )
        self.assertTrue(
            _contains_float_matrix(enclosure.time_time, numpy_tensors.time_time)
        )
        self.assertTrue(
            _contains_float_matrix(
                enclosure.time_space[0], numpy_tensors.time_space[0]
            )
        )
        record = sgbl_enclose_principal_cone(box)
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertEqual(
            record.inconclusive_reason, "resolvent_margin_not_strictly_positive"
        )
        self.assertGreater(record.companion_deformation_upper, 1)
        self.assertTrue(record.structural_hessian_symmetric)
        self.assertTrue(record.structural_pure_gauge_identity)
        self.assertIs(record.representative_generator_sufficient, False)
        identities = sgbl_structural_identities(box)
        self.assertEqual(identities["structural_pure_gauge_identity"], True)
        self.assertEqual(identities["checked_hessian_generators"], 10)
        self.assertEqual(identities["checked_curvature_generators"], 20)
        self.assertIs(identities["representative_generator_sufficient"], False)
        self.assertEqual(record.exact_eigenframe_rank, 24)
        self.assertFalse(record.real_complete_basis)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)

    def test_small_nonflat_product_box_is_not_a_family_health_pass(self) -> None:
        box = sgbl_nonflat_source_principal_box(riemann_half_width=Q(1, 1 << 20))
        self.assertFalse(box.is_singleton)
        self.assertFalse(box.curvature_free)
        try:
            record = sgbl_enclose_principal_cone(box)
        except SGBLConeStop as stopped:
            self.assertEqual(stopped.reason, "resource_limit")
            return
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertFalse(record.real_complete_basis)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)

    def test_adapter_box_keeps_the_linear_branch_derivative_contract(self) -> None:
        state = sgbl_source_state(_fixture_a(), sgbl_source_solve(_fixture_a()))
        background = sgbl_covariant_principal_background(state)
        box = sgbl_principal_box_from_state(state)
        self.assertEqual(box.effective_planck_prime, Interval.singleton(0))
        self.assertEqual(box.chart, ORTHONORMAL_MINKOWSKI_CHART)
        self.assertAlmostEqual(
            float(box.gb_coupling_prime.midpoint()),
            background.gb_coupling_prime,
            places=12,
        )
        self.assertAlmostEqual(
            float(box.effective_planck_squared.midpoint()),
            background.effective_planck_squared,
            places=12,
        )
        reconstructed = np.array(
            [
                [
                    float((box.alpha_gb * box.hessian_phi_lower[a][b]).midpoint())
                    for b in range(4)
                ]
                for a in range(4)
            ]
        )
        np.testing.assert_allclose(
            reconstructed, background.hessian_gb_lower, rtol=0.0, atol=1.0e-12
        )
        self.assertFalse(hasattr(box, "SGBL_branch_owned_and_healthy"))

    def test_module_does_not_import_hyp2_pass_bits_or_qr_builder(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_cone as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertTrue(set(FORBIDDEN_HEALTH_IMPORTS).isdisjoint(imported))
        self.assertNotIn("WeakCouplingThresholds", imported)
        self.assertNotIn("weak_coupling_health_certificate", imported)
        self.assertNotIn("esf_reference_cluster_certificate", imported)
        self.assertNotIn("background_from_spherical_state", imported)
        source = Path(owner.__file__).read_text()
        self.assertNotIn("sample_count=512", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)


if __name__ == "__main__":
    unittest.main()
