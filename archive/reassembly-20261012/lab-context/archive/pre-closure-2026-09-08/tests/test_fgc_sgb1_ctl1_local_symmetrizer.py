from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_cone import (  # noqa: E402
    ORTHONORMAL_MINKOWSKI_CHART,
    SGBLConeLimits,
    SGBLConeStop,
    SGBLPrincipalBackgroundBox,
    sgbl_esf_principal_box,
    sgbl_nonflat_source_principal_box,
    sgbl_pure_gauge_identity_defect,
    sgbl_structural_generator_basis,
    sgbl_structural_identities,
)
from recursive_horizons.fgc.sgb1_ctl1_local_symmetrizer import (  # noqa: E402
    ALL_DIRECTION_INCOMPLETE,
    FAMILY_BOX_API,
    FORBIDDEN_HEALTH_IMPORTS,
    INCOMPLETE_CLASSIFICATION,
    INCONCLUSIVE_CLASSIFICATION,
    LOCAL_DIRECTION,
    LOCAL_THEOREM,
    NONFLAT_OBSTRUCTION,
    PROVED_CLASSIFICATION,
    SGBLLocalSymmetrizerLimits,
    SGBLLocalSymmetrizerRecord,
    SGBLLocalSymmetrizerStop,
    sgbl_assert_distinct_cluster_centers,
    sgbl_certify_companion_eigenframe,
    sgbl_coefficientwise_physical_scalar_identities,
    sgbl_injected_coalesced_speeds,
    sgbl_injected_complex_companion,
    sgbl_injected_defective_companion,
    sgbl_krawczyk_simple_eigenpair,
    sgbl_local_physical_energy,
    sgbl_local_symmetrizer_certificate,
    sgbl_local_symmetrizer_family_contract,
    sgbl_local_symmetrizer_from_family_box,
    sgbl_seed_floating_spectra,
)


class SGBLLocalSymmetrizerInputTests(unittest.TestCase):
    def test_chart_sign_and_copied_hyp2_sample_count_are_refused(self) -> None:
        with self.assertRaises(SGBLConeStop) as chart:
            SGBLPrincipalBackgroundBox(
                planck_mass=2, alpha_gb=Q(-1, 4), chart="metric_dtt"
            )
        self.assertEqual(chart.exception.reason, "cone_chart_error")
        with self.assertRaises(SGBLConeStop) as sign:
            SGBLPrincipalBackgroundBox(planck_mass=interval(-1, 1), alpha_gb=Q(-1, 4))
        self.assertEqual(sign.exception.reason, "sign_chart_error")
        with self.assertRaises(SGBLConeStop) as copied:
            SGBLConeLimits(arc_count=512)
        self.assertEqual(copied.exception.reason, "resource_limit")

    def test_local_limits_are_declared_caps_not_fitted_gaps(self) -> None:
        with self.assertRaises(ValueError):
            SGBLLocalSymmetrizerLimits(max_seed_candidates=0)
        with self.assertRaises(ValueError):
            SGBLLocalSymmetrizerLimits(krawczyk_displacement=0)
        limits = SGBLLocalSymmetrizerLimits(max_seed_denominator=8)
        self.assertEqual(limits.max_seed_denominator, 8)
        self.assertGreater(limits.krawczyk_displacement, 0)


class SGBLLocalSymmetrizerReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        cls.record = sgbl_local_symmetrizer_certificate(cls.box)

    def test_reference_proves_the_point_local_24_mode_cone(self) -> None:
        record = self.record
        self.assertEqual(record.classification, PROVED_CLASSIFICATION)
        self.assertEqual(record.theorem, LOCAL_THEOREM)
        self.assertEqual(record.certified_direction, LOCAL_DIRECTION)
        self.assertTrue(record.real_complete_basis)
        self.assertTrue(record.positive_symmetrizer)
        self.assertTrue(record.strongly_hyperbolic)
        self.assertEqual(record.exact_eigenframe_rank, 24)
        self.assertEqual(record.certified_mode_count, 24)
        self.assertEqual(len(record.clusters), 6)
        for cluster in record.clusters.values():
            self.assertEqual(cluster["exact_nullity"], 4)
            self.assertEqual(cluster["eigenvalue_count"], 4)
            self.assertEqual(cluster["method"], "exact_rational_kernel")
            self.assertTrue(cluster["seed_only_floating_spectra"])
        self.assertIsNotNone(record.riesz_projectors)
        for projector in record.riesz_projectors.values():
            self.assertTrue(projector["idempotent"])
            self.assertEqual(projector["trace"], 4)
            self.assertTrue(projector["commutes_with_companion"])
            self.assertTrue(projector["pure_eigenprojector"])
            self.assertEqual(projector["method"], "exact_V_I_J_V_inverse")
        self.assertGreater(record.kinetic_singular_value_lower, 0)
        self.assertEqual(record.physical_energy_coercivity_lower, 2)
        self.assertTrue(record.structural_hessian_symmetric)
        self.assertTrue(record.structural_pure_gauge_identity)
        self.assertIs(record.representative_generator_sufficient, False)
        self.assertEqual(
            sgbl_structural_identities(self.box)["checked_hessian_generators"], 10
        )
        self.assertEqual(
            sgbl_structural_identities(self.box)["checked_curvature_generators"], 20
        )
        self.assertEqual(record.action["beta"], 0)
        self.assertEqual(record.action["eta"], 0)
        self.assertEqual(record.action["F_prime"], 0)
        self.assertEqual(record.box.chart, ORTHONORMAL_MINKOWSKI_CHART)
        self.assertTrue(self.box.curvature_free)
        self.assertTrue(self.box.is_singleton)

    def test_floating_spectra_seed_but_do_not_prove(self) -> None:
        self.assertFalse(self.record.proof_uses_floating_eig)
        self.assertFalse(self.record.proof_uses_sampled_directions)
        self.assertFalse(self.record.proof_uses_fitted_gap)
        self.assertFalse(self.record.proof_uses_esf_perturbation)
        from recursive_horizons.fgc.sgb1_ctl1_local_symmetrizer import (
            _radial_companion_and_action,
        )

        companion, _a, _b, _inv, _time = _radial_companion_and_action(self.box)
        seed = sgbl_seed_floating_spectra(companion)
        self.assertFalse(seed["proof_uses_floating_eig"])
        self.assertEqual(seed["role"], "candidate_seed_only")
        self.assertLess(seed["max_imaginary_part"], 1.0e-10)

    def test_aggregate_health_stays_false_and_cannot_be_attached(self) -> None:
        record = self.record
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
        with self.assertRaises(ValueError):
            replace(record, proof_uses_floating_eig=True)
        with self.assertRaises(ValueError):
            replace(record, family_covering_proved=True)

    def test_all_direction_is_typed_incomplete_without_a_24_mode_construction(
        self,
    ) -> None:
        self.assertEqual(self.record.all_direction_status, ALL_DIRECTION_INCOMPLETE)
        identities = sgbl_coefficientwise_physical_scalar_identities(self.box)
        self.assertTrue(identities["holds"])
        self.assertTrue(identities["determining_quadratic_grid"])
        self.assertTrue(identities["sampled_directions_are_not_the_proof"])
        self.assertFalse(identities["complete_24_mode_all_direction"])
        self.assertTrue(self.record.coefficientwise_physical_scalar_identities)


class SGBLLocalSymmetrizerNonflatTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.box = sgbl_nonflat_source_principal_box()
        cls.record = sgbl_local_symmetrizer_certificate(cls.box)

    def test_exact_nonflat_fixture_a_is_one_typed_obstruction(self) -> None:
        record = self.record
        self.assertFalse(self.box.curvature_free)
        self.assertTrue(self.box.is_singleton)
        self.assertEqual(record.classification, INCOMPLETE_CLASSIFICATION)
        self.assertEqual(record.obstruction, NONFLAT_OBSTRUCTION)
        self.assertFalse(record.real_complete_basis)
        self.assertFalse(record.positive_symmetrizer)
        self.assertLess(record.certified_mode_count, 24)
        self.assertGreaterEqual(record.certified_mode_count, 16)
        self.assertEqual(record.clusters["tilde_plus"]["exact_nullity"], 4)
        self.assertEqual(record.clusters["hat_minus"]["exact_nullity"], 4)
        self.assertEqual(record.clusters["physical_plus"]["exact_nullity"], 1)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertEqual(record.all_direction_status, ALL_DIRECTION_INCOMPLETE)
        self.assertTrue(record.structural_hessian_symmetric)
        self.assertTrue(record.structural_pure_gauge_identity)

    def test_small_interval_box_is_inconclusive_not_health(self) -> None:
        box = sgbl_esf_principal_box(
            planck_mass=2,
            alpha_gb=Q(-1, 4),
            riemann_half_width=Q(1, 1 << 20),
        )
        self.assertFalse(box.is_singleton)
        record = sgbl_local_symmetrizer_certificate(box)
        self.assertEqual(record.classification, INCONCLUSIVE_CLASSIFICATION)
        self.assertEqual(record.inconclusive_reason, "interval_companion_not_singleton")
        self.assertFalse(record.real_complete_basis)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)


class SGBLLocalSymmetrizerFamilyAPITests(unittest.TestCase):
    def test_family_box_entrypoint_is_exposed_and_does_not_cover(self) -> None:
        contract = sgbl_local_symmetrizer_family_contract()
        self.assertEqual(contract["accepted_slot"], "SGBLPrincipalBackgroundBox")
        self.assertEqual(
            contract["later_owner"],
            "uniform_family_cell_product_box_local_symmetrizer",
        )
        self.assertIs(contract["family_covering_proved"], False)
        self.assertIs(contract["aggregate_health"], False)
        self.assertEqual(contract["entrypoint"], FAMILY_BOX_API["entrypoint"])
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        record = sgbl_local_symmetrizer_from_family_box(box)
        self.assertTrue(record.family_box_entrypoint)
        self.assertFalse(record.family_covering_proved)
        self.assertEqual(record.later_family_owner, contract["later_owner"])
        self.assertEqual(record.classification, PROVED_CLASSIFICATION)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)


class SGBLLocalSymmetrizerAttackTests(unittest.TestCase):
    def test_defective_jordan_block_is_geometric_multiplicity_one(self) -> None:
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_certify_companion_eigenframe(sgbl_injected_defective_companion())
        self.assertEqual(stopped.exception.reason, "defective_cluster")
        self.assertLess(stopped.exception.payload["geometric_multiplicity"], 2)

    def test_complex_block_has_no_real_eigenvalue(self) -> None:
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_certify_companion_eigenframe(sgbl_injected_complex_companion())
        self.assertEqual(stopped.exception.reason, "complex_spectrum")
        self.assertLess(stopped.exception.payload["discriminant"], 0)

    def test_coalesced_centres_are_refused(self) -> None:
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_assert_distinct_cluster_centers(sgbl_injected_coalesced_speeds())
        self.assertEqual(stopped.exception.reason, "cluster_coalescence")
        distinct = sgbl_assert_distinct_cluster_centers(
            (Q(-1), Q(-1, 2), Q(-1, 3), Q(1, 3), Q(1, 2), Q(1))
        )
        self.assertEqual(len(distinct), 6)

    def test_eigenframe_omission_fails_closed(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_local_symmetrizer_certificate(box, omit_eigenframe_columns=4)
        self.assertEqual(stopped.exception.reason, "eigenframe_incomplete")

    def test_symmetry_mutation_is_broken_gauge(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        true_defect = sgbl_pure_gauge_identity_defect(box)
        mutated = sgbl_pure_gauge_identity_defect(box, mutate_sign=-1)
        self.assertEqual(true_defect.abs_upper(), 0)
        self.assertGreater(mutated.abs_upper(), 0)
        with self.assertRaises(SGBLConeStop) as omitted:
            sgbl_structural_generator_basis(omit_hessian_index=0)
        self.assertEqual(omitted.exception.reason, "broken_gauge")

    def test_coercivity_sign_flip_fails_closed(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_local_physical_energy(box, mutate_sign=-1)
        self.assertEqual(stopped.exception.reason, "coercivity_failed")
        with self.assertRaises(SGBLLocalSymmetrizerStop) as certificate:
            sgbl_local_symmetrizer_certificate(box, mutate_energy_sign=-1)
        self.assertEqual(certificate.exception.reason, "coercivity_failed")

    def test_resource_limit_stops_before_a_fitted_margin(self) -> None:
        box = sgbl_esf_principal_box(
            planck_mass=2,
            alpha_gb=Q(-1, 4),
            limits=SGBLConeLimits(max_residual_evaluations=1),
        )
        with self.assertRaises(SGBLLocalSymmetrizerStop) as stopped:
            sgbl_local_symmetrizer_certificate(box)
        self.assertEqual(stopped.exception.reason, "resource_limit")
        with self.assertRaises(SGBLLocalSymmetrizerStop) as candidates:
            sgbl_certify_companion_eigenframe(
                tuple(
                    tuple(Q(row + column + 1) for column in range(24))
                    for row in range(24)
                ),
                limits=SGBLLocalSymmetrizerLimits(max_seed_candidates=1),
            )
        self.assertEqual(candidates.exception.reason, "resource_limit")

    def test_tampered_digest_is_refused(self) -> None:
        box = sgbl_esf_principal_box(planck_mass=2, alpha_gb=Q(-1, 4))
        record = sgbl_local_symmetrizer_certificate(box)
        with self.assertRaises(ValueError):
            SGBLLocalSymmetrizerRecord(
                box=record.box,
                branch=record.branch,
                action=dict(record.action),
                coefficient_enclosure=record.coefficient_enclosure,
                theorem=record.theorem,
                certified_direction=record.certified_direction,
                clusters=dict(record.clusters),
                riesz_projectors=dict(record.riesz_projectors),
                kinetic_rho_infinity=record.kinetic_rho_infinity,
                kinetic_singular_value_lower=record.kinetic_singular_value_lower,
                physical_energy_coercivity_lower=record.physical_energy_coercivity_lower,
                exact_eigenframe_rank=record.exact_eigenframe_rank,
                certified_mode_count=record.certified_mode_count,
                real_complete_basis=record.real_complete_basis,
                positive_symmetrizer=record.positive_symmetrizer,
                all_direction_status=record.all_direction_status,
                coefficientwise_physical_scalar_identities=(
                    record.coefficientwise_physical_scalar_identities
                ),
                classification=record.classification,
                inconclusive_reason=record.inconclusive_reason,
                obstruction=record.obstruction,
                family_covering_proved=False,
                family_box_entrypoint=False,
                sha256="0" * 64,
                structural_hessian_symmetric=True,
                structural_pure_gauge_identity=True,
                proof_uses_floating_eig=False,
                proof_uses_sampled_directions=False,
                proof_uses_fitted_gap=False,
                proof_uses_esf_perturbation=False,
            )


class SGBLLocalSymmetrizerKrawczykTests(unittest.TestCase):
    def test_simple_2x2_eigenpair_is_krawczyk_certified(self) -> None:
        matrix = ((Q(0), Q(1)), (Q(-2), Q(3)))
        payload = sgbl_krawczyk_simple_eigenpair(
            matrix,
            lambda_seed=Q(1),
            vector_seed=(Q(1), Q(1)),
        )
        self.assertIsNotNone(payload)
        self.assertTrue(payload["krawczyk_image_strictly_inside_displacement_box"])
        self.assertLess(payload["rho_infinity_upper_bound"], 1)
        self.assertEqual(payload["method"], "rational_interval_parametric_Krawczyk")
        self.assertFalse(payload.get("proof_uses_floating_eig", False))


class SGBLLocalSymmetrizerIsolationTests(unittest.TestCase):
    def test_module_does_not_import_hyp2_pass_bits_or_the_esf_perturbation_owner(
        self,
    ) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_local_symmetrizer as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertTrue(set(FORBIDDEN_HEALTH_IMPORTS).isdisjoint(imported))
        self.assertNotIn("sgbl_enclose_principal_cone", imported)
        self.assertNotIn("esf_reference_cluster_certificate", imported)
        self.assertNotIn("WeakCouplingThresholds", imported)
        source = Path(owner.__file__).read_text()
        self.assertNotIn("sample_count=512", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertNotIn("proof_uses_fitted_gap=True", source)


if __name__ == "__main__":
    unittest.main()
