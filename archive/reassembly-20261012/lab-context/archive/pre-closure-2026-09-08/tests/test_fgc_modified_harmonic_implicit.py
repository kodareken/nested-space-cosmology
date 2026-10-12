from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.reproduce_fgc_hyp1_mhg_implicit import load_config  # noqa: E402
from scripts.reproduce_fgc_hyp1_reduction import (  # noqa: E402
    load_config as load_reduction_config,
)
from recursive_horizons.fgc.exact_tangent import (  # noqa: E402
    FirstTangent,
    primal_and_tangent,
)
from recursive_horizons.fgc.modified_harmonic_implicit import (  # noqa: E402
    EXPECTED_FGCQR_PARAMETERS,
    exact_full_residual_acceleration_jacobian,
    modified_harmonic_implicit_certificate,
    required_imp1_nonclaims,
)
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_full_residuals,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    BASE_FIELD_ORDER,
    state_from_generalized_adm_pg_fixture,
)


class FGCModifiedHarmonicImplicitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_config()
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=cls.config.radial_domain_minimum
        )
        cls.adapter = {
            "fixture": cls.config.fixture,
            "reference": cls.reference,
            "coordinate_radius": cls.config.coordinate_radius,
            "tilde_normal_factor": cls.config.tilde_normal_factor,
            "hat_normal_factor": cls.config.hat_normal_factor,
        }
        cls.certificate = modified_harmonic_implicit_certificate(cls.adapter)

    def test_exact_first_tangent_arithmetic_and_validation(self) -> None:
        x = FirstTangent.seed(Q(3, 2), Q(5, 7))
        value = (x**3 + 2 * x - 1) / (x + 1)
        primal = (Q(3, 2) ** 3 + 2 * Q(3, 2) - 1) / (Q(3, 2) + 1)
        derivative = (
            ((3 * Q(3, 2) ** 2 + 2) * (Q(3, 2) + 1))
            - (Q(3, 2) ** 3 + 2 * Q(3, 2) - 1)
        ) / (Q(3, 2) + 1) ** 2 * Q(5, 7)
        self.assertEqual(primal_and_tangent(value), (primal, derivative))
        with self.assertRaisesRegex(TypeError, "exact rational"):
            FirstTangent(0.5)  # type: ignore[arg-type]
        with self.assertRaisesRegex(ZeroDivisionError, "zero primal"):
            FirstTangent.seed(Q(0)).reciprocal()

    def test_flat_fgcqr_root_and_exact_ift_certificate(self) -> None:
        certificate = self.certificate
        solution = certificate["reference_solution"]
        self.assertEqual(solution["model_id"], "FGC-QR")
        self.assertEqual(
            solution["action_parameters"], EXPECTED_FGCQR_PARAMETERS
        )
        self.assertEqual(solution["effective_planck_coefficient"], Q(4))
        self.assertTrue(all(value == 0 for value in solution["full_residual_vector"]))
        self.assertTrue(
            all(value == 0 for value in solution["gauge_extension_residual_vector"])
        )
        self.assertTrue(all(value == 0 for value in solution["constraint_up"]))
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            certificate[
                "local_smooth_implicit_coordinate_time_acceleration_branch_proven"
            ]
        )
        self.assertEqual(
            certificate["implicit_branch_neighborhood"],
            "exists_but_not_quantified",
        )

    def test_full_residual_jacobian_is_exact_mhg1_kinetic_block(self) -> None:
        acceleration = self.certificate["acceleration_map"]
        expected = (
            (Q(36), Q(0), Q(9), Q(9), Q(0), Q(0)),
            (Q(0), Q(72), Q(0), Q(0), Q(0), Q(0)),
            (Q(4), Q(0), Q(1), Q(-1), Q(0), Q(0)),
            (Q(64), Q(0), Q(-16), Q(0), Q(0), Q(0)),
            (Q(0), Q(0), Q(0), Q(0), Q(-1), Q(0)),
            (Q(0), Q(0), Q(0), Q(0), Q(0), Q(-1)),
        )
        self.assertEqual(acceleration["jacobian"], expected)
        self.assertEqual(
            acceleration["jacobian"],
            acceleration["mhg1_coordinate_time_kinetic_block"],
        )
        self.assertEqual(acceleration["jacobian_determinant"], Q(-165888))
        self.assertEqual(acceleration["jacobian_rank"], 6)
        self.assertTrue(acceleration["left_inverse_identity_exact"])
        self.assertTrue(acceleration["right_inverse_identity_exact"])
        self.assertTrue(acceleration["seeded_primals_reproduce_base_residual"])

    def test_forward_tangent_propagates_at_nonflat_activated_fixture(self) -> None:
        fixture = load_reduction_config().configuration["fixtures"][3]
        state = state_from_generalized_adm_pg_fixture(fixture)
        tangent = exact_full_residual_acceleration_jacobian(
            state,
            reference=self.reference,
            coordinate_radius=Q(9, 2),
            tilde_normal_factor=Q(4),
            hat_normal_factor=Q(9),
        )
        self.assertTrue(tangent["seeded_primals_reproduce_base_residual"])
        self.assertEqual(len(tangent["jacobian"]), 6)
        self.assertTrue(all(len(row) == 6 for row in tangent["jacobian"]))

        # For fixed lower jets the complete REF1 residual is at most quadratic
        # in each second-jet slot. An ordinary-rational symmetric slope is
        # therefore an exact derivative oracle, independent of FirstTangent.
        step = Q(1)
        for column, field in enumerate(BASE_FIELD_ORDER):
            jet = getattr(state, field)
            plus_state = replace(
                state,
                **{field: replace(jet, dtt=jet.dtt + step)},
            )
            minus_state = replace(
                state,
                **{field: replace(jet, dtt=jet.dtt - step)},
            )
            plus = modified_harmonic_full_residuals(
                plus_state,
                reference=self.reference,
                coordinate_radius=Q(9, 2),
                tilde_normal_factor=Q(4),
                hat_normal_factor=Q(9),
            )["full_residual_vector"]
            minus = modified_harmonic_full_residuals(
                minus_state,
                reference=self.reference,
                coordinate_radius=Q(9, 2),
                tilde_normal_factor=Q(4),
                hat_normal_factor=Q(9),
            )["full_residual_vector"]
            exact_central_column = tuple(
                (plus_row - minus_row) / (2 * step)
                for plus_row, minus_row in zip(plus, minus, strict=True)
            )
            tangent_column = tuple(
                tangent["jacobian"][row][column]
                for row in range(len(tangent["jacobian"]))
            )
            self.assertEqual(tangent_column, exact_central_column)

    def test_root_and_configuration_mutations_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            modified_harmonic_implicit_certificate({**self.adapter, "extra": True})

        changed_parameter = deepcopy(self.config.fixture)
        changed_parameter["action_parameters"]["ricci_coupling"] = Q(-1, 5)
        with self.assertRaisesRegex(ValueError, "flat FGC-QR root failed"):
            modified_harmonic_implicit_certificate(
                {**self.adapter, "fixture": changed_parameter}
            )

        changed_radius = deepcopy(self.config.fixture)
        changed_radius["state"]["areal_radius"]["dr"] = Q(2)
        with self.assertRaisesRegex(ValueError, "flat FGC-QR root failed"):
            modified_harmonic_implicit_certificate(
                {**self.adapter, "fixture": changed_radius}
            )

        changed_acceleration = deepcopy(self.config.fixture)
        changed_acceleration["state"]["phi"]["dtt"] = Q(1)
        with self.assertRaisesRegex(ValueError, "flat FGC-QR root failed"):
            modified_harmonic_implicit_certificate(
                {**self.adapter, "fixture": changed_acceleration}
            )

    def test_nonclaims_remain_false(self) -> None:
        self.assertTrue(all(value is False for value in required_imp1_nonclaims().values()))
        self.assertEqual(self.certificate["nonclaims"], required_imp1_nonclaims())


if __name__ == "__main__":
    unittest.main()
