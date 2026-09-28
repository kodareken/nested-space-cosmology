from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from decimal import Decimal
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_extension_residual,
    modified_harmonic_full_residuals,
    modified_harmonic_gauge_constraint,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.sgb1_ctl1_constraints import (  # noqa: E402
    SGBLConstraintSolveStop,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceCoefficients,
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    SOURCE_EQUATION_ORDER,
    TILDE_NORMAL_FACTOR,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_solve,
    sgbl_source_state,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    residuals,
)


Q = Fraction
ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
EXPECTED_DET_A = Q(6519803496633, 6553600000)
EXPECTED_DET_B = Q(57291897169490646528, 7593631744384765625)
EXPECTED_DET_A_ZERO_COUPLING = Q(2187, 2)
HOLDOUT = (Q(2), Q(-1), Q(1, 2), Q(3), Q(-2), Q(1))
REFERENCE = flat_spherical_annulus_reference(
    radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
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


def _fixture_b() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=3,
        alpha=(Q(5, 2), Q(-2, 7), Q(3, 4), Q(-1, 8), Q(2, 11)),
        shift=(Q(-2, 5), Q(1, 6), Q(-3, 7), Q(2, 9), Q(-1, 5)),
        radial_metric=(Q(4, 3), Q(-1, 4), Q(2, 5), Q(-1, 9), Q(3, 8)),
        areal_radius=(Q(7, 2), Q(3, 8), Q(-5, 4), Q(2, 7), Q(1, 10)),
        phi=(Q(-2, 3), Q(-5, 8), Q(1, 2), Q(-3, 11), Q(2, 13)),
        chi=(Q(4, 5), Q(-1, 3), Q(-2, 9), Q(3, 10), Q(-1, 12)),
        planck_mass=Q(3, 2),
        scalar_mass=Q(2, 3),
        quartic_coupling=Q(5, 7),
        alpha_gb=Q(2, 5),
    )


def _flat_inputs(*, alpha_gb: Fraction | int = Q(-1, 4)) -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(Q(5, 2), 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=alpha_gb,
    )


def _schwarzschild_inputs() -> SGBLSourceInputs:
    # Static areal-gauge Schwarzschild, r_s=3, R=r=4.
    return SGBLSourceInputs(
        coordinate_radius=4,
        alpha=(Q(1, 2), 0, Q(3, 16), 0, Q(-21, 128)),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(2, 0, Q(-3, 4), 0, Q(39, 32)),
        areal_radius=(4, 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


def _quadratic_nodes() -> tuple[tuple[int, ...], ...]:
    nodes = [ZERO6]
    for index in range(6):
        plus = [0] * 6
        plus[index] = 1
        minus = [0] * 6
        minus[index] = -1
        nodes.append(tuple(plus))
        nodes.append(tuple(minus))
    for i, j in combinations(range(6), 2):
        pair = [0] * 6
        pair[i] = 1
        pair[j] = 1
        nodes.append(tuple(pair))
    assert len(nodes) == 28
    return tuple(nodes)


def _unredefined_residual(
    point: SGBLSourceInputs,
    accelerations: tuple[Fraction | int, ...],
) -> tuple[Fraction, ...]:
    result = modified_harmonic_full_residuals(
        sgbl_source_state(point, accelerations),
        reference=REFERENCE,
        coordinate_radius=point.coordinate_radius,
        tilde_normal_factor=TILDE_NORMAL_FACTOR,
        hat_normal_factor=HAT_NORMAL_FACTOR,
    )
    return tuple(result["unredefined_residual_vector"])


def _assembled_warped_residual(
    point: SGBLSourceInputs,
    accelerations: tuple[Fraction | int, ...],
) -> tuple[Fraction, ...]:
    state = sgbl_source_state(point, accelerations)
    unredefined = residuals(state, use_warped=True)
    gauge = modified_harmonic_gauge_constraint(
        state,
        reference=REFERENCE,
        coordinate_radius=point.coordinate_radius,
        tilde_normal_factor=TILDE_NORMAL_FACTOR,
    )
    extension = modified_harmonic_extension_residual(
        state, gauge=gauge, hat_normal_factor=HAT_NORMAL_FACTOR,
    )
    metric = unredefined["metric"]
    covariant = extension["covariant"]
    return (
        metric[0][0] + covariant[0][0],
        metric[0][1] + covariant[0][1],
        metric[1][1] + covariant[1][1],
        metric[2][2] + covariant[2][2],
        unredefined["phi"],
        unredefined["chi"],
    )


def _zero_chi_momenta(point: SGBLSourceInputs) -> SGBLSourceInputs:
    chi = point.chi
    return replace(point, chi=Jet2(chi.value, dt=0, dr=chi.dr, dtr=0, drr=chi.drr))


def _zero_phi_momenta(point: SGBLSourceInputs) -> SGBLSourceInputs:
    phi = point.phi
    return replace(point, phi=Jet2(phi.value, dt=0, dr=phi.dr, dtr=0, drr=phi.drr))


class SGBLSourceInputTests(unittest.TestCase):
    def test_inputs_normalize_builtin_integers_and_are_frozen(self) -> None:
        point = _flat_inputs()
        self.assertEqual(type(point.coordinate_radius), Fraction)
        self.assertEqual(type(point.alpha_gb), Fraction)
        self.assertTrue(all(type(getattr(point.alpha, name)) is Fraction for name in (
            "value", "dt", "dr", "dtt", "dtr", "drr",
        )))
        self.assertEqual(point.alpha.dtt, Q(0))
        with self.assertRaises(FrozenInstanceError):
            point.coordinate_radius = Q(3)
        self.assertEqual(replace(point, alpha_gb=0).alpha_gb, Q(0))

    def test_jet2_with_zero_dtt_is_accepted(self) -> None:
        point = _flat_inputs()
        jet = Jet2(1, dt=Q(1, 2), dr=Q(-1, 3), dtt=0, dtr=Q(1, 4), drr=Q(-1, 5))
        rebuilt = replace(point, alpha=jet)
        self.assertEqual(rebuilt.alpha, jet)
        self.assertEqual(rebuilt.alpha.dtt, Q(0))

    def test_nonzero_dtt_is_refused_rather_than_dropped(self) -> None:
        point = _flat_inputs()
        with self.assertRaisesRegex(ValueError, "dtt"):
            replace(point, phi=Jet2(0, dtt=1))
        with self.assertRaisesRegex(ValueError, "dtt"):
            replace(point, alpha=(1, 0, 0, 0, 0, 1))

    def test_every_scalar_input_refuses_nonexact_aliases(self) -> None:
        point = _flat_inputs()
        scalar_fields = (
            "coordinate_radius", "planck_mass", "scalar_mass",
            "quartic_coupling", "alpha_gb",
        )
        for field, invalid in product(
            scalar_fields, (True, False, 1.0, float("nan"), float("inf"), "1", Decimal(1)),
        ):
            with self.subTest(field=field, invalid=invalid):
                with self.assertRaisesRegex(TypeError, field):
                    replace(point, **{field: invalid})

    def test_every_jet_component_refuses_nonexact_aliases(self) -> None:
        point = _flat_inputs()
        for field in ("alpha", "shift", "radial_metric", "areal_radius", "phi", "chi"):
            for index, invalid in product(
                range(5), (True, False, 1.0, "1", Decimal(1)),
            ):
                values = [1, 0, 0, 0, 0]
                values[index] = invalid
                with self.subTest(field=field, index=index, invalid=invalid):
                    with self.assertRaisesRegex(TypeError, field):
                        replace(point, **{field: tuple(values)})

    def test_positive_lapse_lambda_radius_and_action_parameters(self) -> None:
        point = _flat_inputs()
        for name, invalid in product(
            ("coordinate_radius", "planck_mass", "scalar_mass", "quartic_coupling"),
            (0, -1, -Q(1, 5)),
        ):
            with self.subTest(name=name, invalid=invalid):
                with self.assertRaisesRegex(ValueError, name):
                    replace(point, **{name: invalid})
        with self.assertRaisesRegex(ValueError, "alpha"):
            replace(point, alpha=(0, 0, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, "radial_metric"):
            replace(point, radial_metric=(0, 0, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, "areal_radius"):
            replace(point, areal_radius=(0, 0, 1, 0, 0))

    def test_point_requires_the_declared_input_type(self) -> None:
        for invalid in (None, {}, 1, Q(1)):
            with self.assertRaisesRegex(TypeError, "SGBLSourceInputs"):
                sgbl_source_coefficients(invalid)

    def test_coefficient_data_are_exact_and_frozen(self) -> None:
        identity = tuple(
            tuple(Q(int(row == column)) for column in range(6))
            for row in range(6)
        )
        coefficients = SGBLSourceCoefficients(constant=(1, 0, 0, 0, 0, 0), jacobian=identity)
        self.assertTrue(all(type(entry) is Fraction for entry in coefficients.constant))
        self.assertTrue(
            all(type(entry) is Fraction for row in coefficients.jacobian for entry in row)
        )
        self.assertEqual(coefficients.jacobian_determinant, Q(1))
        self.assertEqual(coefficients.acceleration_order, SOURCE_ACCELERATION_ORDER)
        self.assertEqual(coefficients.equation_order, SOURCE_EQUATION_ORDER)
        with self.assertRaises(FrozenInstanceError):
            coefficients.constant = ZERO6
        for invalid in (True, False, 1.0, "1"):
            with self.assertRaisesRegex(TypeError, "constant"):
                replace(coefficients, constant=(invalid, 0, 0, 0, 0, 0))

    def test_public_records_cannot_forge_identities_or_solved_claims(self) -> None:
        identity = tuple(
            tuple(Q(int(row == column)) for column in range(6))
            for row in range(6)
        )
        coefficients = SGBLSourceCoefficients(constant=ZERO6, jacobian=identity)
        for forged in (
            dict(quadratic_vanished=True),
            dict(identity_justified=True),
            dict(solved_residual_is_zero=True),
            dict(jacobian_determinant=Q(1)),
        ):
            with self.subTest(forged=forged):
                with self.assertRaises(TypeError):
                    SGBLSourceCoefficients(constant=ZERO6, jacobian=identity, **forged)
                with self.assertRaises(TypeError):
                    replace(coefficients, **forged)
        singular = replace(coefficients, jacobian=tuple((0,) * 6 for _ in range(6)))
        self.assertEqual(singular.jacobian_determinant, Q(0))
        self.assertEqual(singular.evaluate(HOLDOUT), ZERO6)

    def test_evaluation_requires_exact_accelerations_even_when_rows_vanish(self) -> None:
        coefficients = SGBLSourceCoefficients(
            constant=ZERO6,
            jacobian=tuple((0,) * 6 for _ in range(6)),
        )
        self.assertEqual(coefficients.evaluate(HOLDOUT), ZERO6)
        for invalid in (True, False, 1.0, "1"):
            with self.assertRaisesRegex(TypeError, "accelerations"):
                coefficients.evaluate((invalid, 0, 0, 0, 0, 0))


class SGBLSourceExactTests(unittest.TestCase):
    def test_actual_off_shell_sgbl_point_reaches_the_typed_singular_source_stop(self) -> None:
        # The metric is flat, but phi_r is nonzero: this is deliberately
        # off shell, not a physical solution or a geometric singularity.
        point = SGBLSourceInputs(
            coordinate_radius=Q(5, 2), alpha=(1, 0, 0, 0, 0),
            shift=(0, 0, 0, 0, 0), radial_metric=(1, 0, 0, 0, 0),
            areal_radius=(Q(5, 2), 0, 1, 0, 0),
            phi=(0, 0, -5, 0, 0), chi=(0, 0, 0, 0, 0),
            planck_mass=2, scalar_mass=1, quartic_coupling=1,
            alpha_gb=-Q(1, 4),
        )
        coefficients = sgbl_source_coefficients(point)
        self.assertEqual(coefficients.jacobian_determinant, 0)
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            sgbl_source_solve(point)
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")

    @classmethod
    def setUpClass(cls) -> None:
        cls.point_a = _fixture_a()
        cls.point_b = _fixture_b()
        cls.coefficients_a = sgbl_source_coefficients(cls.point_a)
        cls.coefficients_b = sgbl_source_coefficients(cls.point_b)

    def test_seven_evaluations_match_the_affine_identity_on_both_fixtures(self) -> None:
        for name, point, coefficients, expected_det in (
            ("A", self.point_a, self.coefficients_a, EXPECTED_DET_A),
            ("B", self.point_b, self.coefficients_b, EXPECTED_DET_B),
        ):
            with self.subTest(fixture=name):
                origin = sgbl_source_residual(point, ZERO6)
                self.assertEqual(coefficients.constant, origin)
                for column in range(6):
                    seed = list(ZERO6)
                    seed[column] = 1
                    basis = sgbl_source_residual(point, seed)
                    self.assertEqual(
                        tuple(coefficients.jacobian[row][column] for row in range(6)),
                        tuple(entry - origin[row] for row, entry in enumerate(basis)),
                    )
                self.assertEqual(coefficients.jacobian_determinant, expected_det)
                self.assertEqual(len(coefficients.jacobian), 6)
                self.assertTrue(all(len(row) == 6 for row in coefficients.jacobian))
                self.assertNotEqual(coefficients.jacobian_determinant, Q(0))

    def test_complete_reconstruction_at_arbitrary_and_holdout_accelerations(self) -> None:
        samples = (HOLDOUT, (Q(1, 2), -3, Q(4, 5), -Q(1, 7), 2, Q(-5, 3)))
        for name, point, coefficients in (
            ("A", self.point_a, self.coefficients_a),
            ("B", self.point_b, self.coefficients_b),
        ):
            for accelerations in samples:
                with self.subTest(fixture=name, accelerations=accelerations):
                    self.assertEqual(
                        coefficients.evaluate(accelerations),
                        sgbl_source_residual(point, accelerations),
                    )

    def test_solved_accelerations_zero_the_actual_full_residual(self) -> None:
        for name, point, coefficients in (
            ("A", self.point_a, self.coefficients_a),
            ("B", self.point_b, self.coefficients_b),
        ):
            with self.subTest(fixture=name):
                algebraic = coefficients.solve()
                self.assertTrue(all(type(entry) is Fraction for entry in algebraic))
                self.assertEqual(coefficients.evaluate(algebraic), ZERO6)
                direct = sgbl_source_residual(point, algebraic)
                self.assertEqual(direct, ZERO6)
                self.assertEqual(sgbl_source_solve(point), algebraic)

    def test_unredefined_six_row_determinant_vanishes(self) -> None:
        for name, point in (("A", self.point_a), ("B", self.point_b)):
            with self.subTest(fixture=name):
                origin = _unredefined_residual(point, ZERO6)
                columns = []
                for index in range(6):
                    seed = list(ZERO6)
                    seed[index] = 1
                    evaluated = _unredefined_residual(point, seed)
                    columns.append(tuple(
                        entry - origin[row] for row, entry in enumerate(evaluated)
                    ))
                jacobian = tuple(
                    tuple(columns[column][row] for column in range(6))
                    for row in range(6)
                )
                unredefined = SGBLSourceCoefficients(constant=origin, jacobian=jacobian)
                self.assertEqual(unredefined.jacobian_determinant, Q(0))
                with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
                    unredefined.solve()
                self.assertEqual(stopped.exception.reason, "singular_source_jacobian")
                self.assertIs(stopped.exception.coefficients, unredefined)

    def test_zero_coupling_limit_on_A_is_not_the_canonical_gr0_collapse_control(self) -> None:
        point = replace(self.point_a, alpha_gb=0)
        self.assertNotEqual(point.phi.value, Q(0))
        coefficients = sgbl_source_coefficients(point)
        self.assertEqual(coefficients.jacobian_determinant, EXPECTED_DET_A_ZERO_COUPLING)
        self.assertEqual(sgbl_source_state(point).branch, "GR-0")
        self.assertEqual(sgbl_source_residual(point, coefficients.solve()), ZERO6)

    def test_warped_assembly_matches_the_mhg_residual_instrument(self) -> None:
        for name, point, accelerations in (
            ("A0", self.point_a, ZERO6),
            ("A_holdout", self.point_a, HOLDOUT),
            ("B0", self.point_b, ZERO6),
        ):
            with self.subTest(name=name):
                self.assertEqual(
                    _assembled_warped_residual(point, accelerations),
                    sgbl_source_residual(point, accelerations),
                )

    def test_scalar_and_chi_momenta_are_retained(self) -> None:
        point = self.point_a
        full = self.coefficients_a
        without_chi = sgbl_source_coefficients(_zero_chi_momenta(point))
        without_phi = sgbl_source_coefficients(_zero_phi_momenta(point))
        self.assertNotEqual(point.chi.dt, Q(0))
        self.assertNotEqual(point.chi.dtr, Q(0))
        self.assertNotEqual(point.phi.dt, Q(0))
        self.assertNotEqual(point.phi.dtr, Q(0))
        self.assertNotEqual(full.constant, without_chi.constant)
        self.assertNotEqual(full.constant, without_phi.constant)
        # Chi stress has no accelerations, so J need not change with Pi_chi.
        # Phi momenta enter the GB/Hessian mixing and do change J.
        self.assertNotEqual(full.jacobian, without_phi.jacobian)
        self.assertEqual(full.jacobian[5][5], without_chi.jacobian[5][5])
        self.assertNotEqual(full.jacobian[4][4], Q(0))
        self.assertNotEqual(full.jacobian[5][5], Q(0))

    def test_source_jacobian_is_not_the_annular_constraint_jacobian(self) -> None:
        self.assertEqual(SGBLSourceJacobianSolveStop.reason, "singular_source_jacobian")
        self.assertNotEqual(
            SGBLSourceJacobianSolveStop.reason,
            SGBLConstraintSolveStop.reason,
        )
        self.assertEqual(len(self.coefficients_a.jacobian), 6)
        self.assertEqual(SOURCE_ACCELERATION_ORDER[2], "lambda_tt")
        self.assertEqual(SOURCE_ACCELERATION_ORDER[3], "R_tt")

    def test_singular_source_jacobian_refuses_even_a_compatible_zero_row(self) -> None:
        coefficients = SGBLSourceCoefficients(
            constant=ZERO6,
            jacobian=tuple((0,) * 6 for _ in range(6)),
        )
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            coefficients.solve()
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")
        self.assertIs(stopped.exception.coefficients, coefficients)
        self.assertEqual(coefficients.evaluate((1, 2, 3, 4, 5, 6)), ZERO6)

    def test_arbitrarily_small_nonzero_determinants_have_no_invented_floor(self) -> None:
        for n in (1, 10 ** 40):
            jacobian = tuple(
                tuple(Q(1, n) if row == column else Q(0) for column in range(6))
                for row in range(6)
            )
            coefficients = SGBLSourceCoefficients(
                constant=(1, 0, 0, 0, 0, 0),
                jacobian=jacobian,
            )
            self.assertEqual(coefficients.jacobian_determinant, Q(1, n) ** 6)
            accelerations = coefficients.solve()
            self.assertEqual(accelerations, (-Q(n), 0, 0, 0, 0, 0))
            self.assertEqual(coefficients.evaluate(accelerations), ZERO6)


class SGBLSourceQuadraticSupportTests(unittest.TestCase):
    def test_all_twenty_one_quadratic_coefficients_vanish_on_both_fixtures(self) -> None:
        nodes = _quadratic_nodes()
        plus = nodes[1:13:2]
        minus = nodes[2:14:2]
        pairs = nodes[13:]
        pair_indices = tuple(combinations(range(6), 2))
        self.assertEqual(len(plus), 6)
        self.assertEqual(len(minus), 6)
        self.assertEqual(len(pairs), 15)
        for name, point, coefficients in (
            ("A", _fixture_a(), sgbl_source_coefficients(_fixture_a())),
            ("B", _fixture_b(), sgbl_source_coefficients(_fixture_b())),
        ):
            with self.subTest(fixture=name):
                values = tuple(sgbl_source_residual(point, node) for node in nodes)
                origin = values[0]
                self.assertEqual(origin, coefficients.constant)
                quadratic = []
                for index in range(6):
                    diagonal = tuple(
                        (plus_value + minus_value) / 2 - origin[row]
                        for row, (plus_value, minus_value) in enumerate(
                            zip(values[1 + 2 * index], values[2 + 2 * index], strict=True)
                        )
                    )
                    quadratic.append(diagonal)
                    linear = tuple(
                        (plus_value - minus_value) / 2
                        for plus_value, minus_value in zip(
                            values[1 + 2 * index], values[2 + 2 * index], strict=True
                        )
                    )
                    self.assertEqual(
                        linear,
                        tuple(coefficients.jacobian[row][index] for row in range(6)),
                    )
                for (i, j), pair_residual in zip(pair_indices, values[13:], strict=True):
                    off_diagonal = tuple(
                        pair_residual[row]
                        - values[1 + 2 * i][row]
                        - values[1 + 2 * j][row]
                        + origin[row]
                        for row in range(6)
                    )
                    quadratic.append(off_diagonal)
                self.assertEqual(len(quadratic), 21)
                for row in range(6):
                    self.assertEqual(
                        tuple(entry[row] for entry in quadratic),
                        (Q(0),) * 21,
                    )
                self.assertEqual(
                    coefficients.evaluate(HOLDOUT),
                    sgbl_source_residual(point, HOLDOUT),
                )


class SGBLSourceControlTests(unittest.TestCase):
    def test_flat_vacuum_is_an_exact_zero_source_for_nonzero_coupling(self) -> None:
        point = _flat_inputs(alpha_gb=Q(-1, 4))
        self.assertEqual(sgbl_source_residual(point, ZERO6), ZERO6)
        coefficients = sgbl_source_coefficients(point)
        self.assertEqual(coefficients.solve(), ZERO6)
        self.assertEqual(sgbl_source_solve(point), ZERO6)
        unredefined = residuals(sgbl_source_state(point), use_warped=True)
        self.assertEqual(unredefined["GB"], Q(0))
        self.assertEqual(unredefined["phi"], Q(0))

    def test_flat_zero_coupling_control_remains_exact(self) -> None:
        point = _flat_inputs(alpha_gb=0)
        self.assertEqual(sgbl_source_state(point).branch, "GR-0")
        self.assertEqual(sgbl_source_residual(point, ZERO6), ZERO6)
        self.assertEqual(sgbl_source_solve(point), ZERO6)

    def test_schwarzschild_with_phi_zero_is_not_a_full_sgbl_solution(self) -> None:
        point = _schwarzschild_inputs()
        state = sgbl_source_state(point)
        unredefined = residuals(state, use_warped=True)
        self.assertEqual(point.phi.value, Q(0))
        self.assertNotEqual(point.alpha_gb, Q(0))
        self.assertEqual(unredefined["R"], Q(0))
        self.assertEqual(unredefined["GB"], Q(27, 1024))
        self.assertEqual(unredefined["phi"], point.alpha_gb * unredefined["GB"])
        self.assertEqual(unredefined["chi"], Q(0))
        self.assertTrue(all(
            unredefined["metric"][a][b] == 0 for a in range(4) for b in range(4)
        ))
        full = sgbl_source_residual(point, ZERO6)
        self.assertEqual(full[4], point.alpha_gb * unredefined["GB"])
        self.assertNotEqual(full, ZERO6)


if __name__ == "__main__":
    unittest.main()
