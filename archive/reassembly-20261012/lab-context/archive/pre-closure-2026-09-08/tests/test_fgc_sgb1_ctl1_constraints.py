from __future__ import annotations

from dataclasses import FrozenInstanceError, fields, replace
from decimal import Decimal
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.chi_principal_identity import FormalPolynomial  # noqa: E402
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    gauge_compatible_metric_time_derivatives,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    physical_constraint_projections,
)
from recursive_horizons.fgc.sgb1_ctl1_constraints import (  # noqa: E402
    SGBLAnnularInputs,
    SGBLConstraintCoefficients,
    SGBLConstraintSolveStop,
    sgbl_constraint_coefficients,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    SphericalState,
    _p,
    _ricci,
    _scalar,
)


Q = Fraction


class _Polynomial:
    """Test-only operator adapter for the repository's exact formal algebra."""

    def __init__(self, value: int | Fraction | FormalPolynomial | _Polynomial = 0):
        self.formal = (
            value.formal if isinstance(value, _Polynomial)
            else value if isinstance(value, FormalPolynomial)
            else FormalPolynomial.constant(value)
        )

    @classmethod
    def atom(cls, name: str) -> _Polynomial:
        return cls(FormalPolynomial.atom(name))

    def __add__(self, other):
        return _Polynomial(self.formal + _Polynomial(other).formal)

    __radd__ = __add__

    def __neg__(self):
        return _Polynomial(-self.formal)

    def __sub__(self, other):
        return self + (-_Polynomial(other))

    def __rsub__(self, other):
        return _Polynomial(other) - self

    def __mul__(self, other):
        return _Polynomial(self.formal * _Polynomial(other).formal)

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction):
        return self * (Q(1) / Q(scalar))

    def __pow__(self, exponent: int):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("test polynomial powers must be nonnegative integers")
        return prod((self for _ in range(exponent)), start=_Polynomial(1))

    def __eq__(self, other):
        if not isinstance(other, (_Polynomial, int, Fraction)):
            return NotImplemented
        return self.formal == _Polynomial(other).formal

    def __repr__(self):
        return self.formal.canonical()

    def substitute(self, values: dict):
        return sum((
            coefficient * prod(
                (values.get(atom, _Polynomial.atom(atom)) for atom in monomial),
                start=_Polynomial(1),
            )
            for monomial, coefficient in self.formal.terms
        ), start=_Polynomial())

    def coefficient(self, lambda_degree: int, k_degree: int):
        return _Polynomial(FormalPolynomial(tuple(
            (tuple(atom for atom in monomial if atom not in ("L_r", "k_r")), value)
            for monomial, value in self.formal.terms
            if (monomial.count("L_r"), monomial.count("k_r")) == (lambda_degree, k_degree)
        )))


def _vacuum_inputs() -> SGBLAnnularInputs:
    return SGBLAnnularInputs(
        radius=1, radial_metric=1, angular_extrinsic_curvature=0,
        phi=0, phi_r=0, phi_rr=0, phi_pi=0, phi_pi_r=0,
        chi_r=0, chi_pi=0, planck_mass=1,
        scalar_mass=1, quartic_coupling=1, alpha_gb=1,
    )


def _nontrivial_inputs() -> tuple[SGBLAnnularInputs, ...]:
    return (
        SGBLAnnularInputs(
            radius=Q(5, 2), radial_metric=Q(4, 3),
            angular_extrinsic_curvature=Q(1, 17),
            phi=Q(1, 29), phi_r=-Q(1, 31), phi_rr=Q(1, 37),
            phi_pi=Q(1, 41), phi_pi_r=-Q(1, 43),
            chi_r=Q(1, 53), chi_pi=-Q(1, 61),
            planck_mass=2, scalar_mass=3, quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        ),
        SGBLAnnularInputs(
            radius=Q(11, 3), radial_metric=Q(7, 6),
            angular_extrinsic_curvature=-Q(2, 19),
            phi=-Q(2, 37), phi_r=Q(1, 41), phi_rr=-Q(1, 43),
            phi_pi=-Q(1, 47), phi_pi_r=Q(1, 53),
            chi_r=-Q(2, 61), chi_pi=Q(1, 71),
            planck_mass=Q(3, 2), scalar_mass=Q(2, 3),
            quartic_coupling=Q(5, 7), alpha_gb=Q(2, 5),
        ),
    )


def _red1_pair(
    point: SGBLAnnularInputs,
    Lr: Fraction | int,
    kr: Fraction | int,
    *,
    free_second_jet: Fraction = Q(0),
) -> tuple[Fraction, Fraction]:
    """Independent full tensor route with the actual ID1 slice/gauge jets."""
    r, L, k = point.radius, point.radial_metric, point.angular_extrinsic_curvature
    h_tt_dt, h_tr_dt = gauge_compatible_metric_time_derivatives(r, L, Lr)
    Lrr = free_second_jet
    # The radial derivative of ID1's gauge-compatible h_tr,t along the slice.
    h_tr_dtr = L * Lr / r - (L**2 - 1) / (2 * r**2) + Lrr / (4 * L) - Lr**2 / (4 * L**2)
    state = SphericalState(
        h_tt=Jet2(-1, dt=h_tt_dt, dtt=free_second_jet),
        h_tr=Jet2(0, dt=h_tr_dt, dtt=2 * free_second_jet, dtr=h_tr_dtr),
        h_rr=Jet2(
            L**2, dt=4 * L**2 * k, dr=2 * L * Lr,
            dtt=3 * free_second_jet,
            dtr=8 * L * Lr * k + 4 * L**2 * kr,
            drr=2 * Lr**2 + 2 * L * Lrr,
        ),
        areal_radius=Jet2(r, dt=-r * k, dr=1, dtt=5 * free_second_jet, dtr=-k - r * kr),
        phi=Jet2(
            point.phi, dt=point.phi_pi, dr=point.phi_r,
            dtt=7 * free_second_jet, dtr=point.phi_pi_r, drr=point.phi_rr,
        ),
        chi=Jet2(
            -Q(1, 47), dt=point.chi_pi, dr=point.chi_r,
            dtt=11 * free_second_jet, dtr=point.chi_pi / 7, drr=point.chi_r / 11,
        ),
        planck_mass=point.planck_mass, mu=point.scalar_mass,
        g4=point.quartic_coupling, alpha=point.alpha_gb,
        beta=0, eta=0, branch="SGB-L" if point.alpha_gb else "GR-0",
    )
    result = physical_constraint_projections(state)
    return result["H"], result["M"]


class SGBLConstraintInputTests(unittest.TestCase):
    def test_inputs_normalize_builtin_integers_and_are_frozen(self) -> None:
        point = _vacuum_inputs()
        self.assertTrue(all(type(getattr(point, field.name)) is Fraction for field in fields(point)))
        with self.assertRaises(FrozenInstanceError):
            point.radius = Q(2)
        self.assertEqual(replace(point, alpha_gb=0).alpha_gb, Q(0))

    def test_every_input_refuses_nonexact_aliases(self) -> None:
        point = _vacuum_inputs()
        for field, invalid in product(fields(point), (True, False, 1.0, float("nan"), float("inf"), "1", Decimal(1))):
            with self.subTest(field=field.name, invalid=invalid):
                with self.assertRaisesRegex(TypeError, field.name):
                    replace(point, **{field.name: invalid})

    def test_annular_and_action_parameters_are_strictly_positive(self) -> None:
        point = _vacuum_inputs()
        for name, invalid in product(
            ("radius", "radial_metric", "planck_mass", "scalar_mass", "quartic_coupling"),
            (0, -1, -Q(1, 5)),
        ):
            with self.subTest(name=name, invalid=invalid):
                with self.assertRaisesRegex(ValueError, name):
                    replace(point, **{name: invalid})

    def test_point_requires_the_declared_input_type(self) -> None:
        for invalid in (None, {}, 1, Q(1)):
            with self.assertRaisesRegex(TypeError, "SGBLAnnularInputs"):
                sgbl_constraint_coefficients(invalid)

    def test_coefficient_data_are_exact_and_frozen(self) -> None:
        coefficients = SGBLConstraintCoefficients(
            hamiltonian_constant=1, hamiltonian_lambda_r=2,
            momentum_constant=3, momentum_k_r=4,
        )
        self.assertTrue(all(type(getattr(coefficients, field.name)) is Fraction for field in fields(coefficients)))
        self.assertEqual(coefficients.jacobian_diagonal, (Q(2), Q(4)))
        self.assertEqual(coefficients.jacobian_determinant, Q(8))
        with self.assertRaises(FrozenInstanceError):
            coefficients.momentum_k_r = Q(0)
        for field, invalid in product(fields(coefficients), (True, False, 1.0, "1")):
            with self.subTest(field=field.name, invalid=invalid):
                with self.assertRaisesRegex(TypeError, field.name):
                    replace(coefficients, **{field.name: invalid})

    def test_evaluation_requires_exact_derivatives_even_when_rows_vanish(self) -> None:
        coefficients = SGBLConstraintCoefficients(
            hamiltonian_constant=0, hamiltonian_lambda_r=0,
            momentum_constant=0, momentum_k_r=0,
        )
        valid = dict(radial_metric_derivative=1, angular_extrinsic_curvature_derivative=Q(2, 3))
        self.assertEqual(coefficients.evaluate(**valid), (Q(0), Q(0)))
        for name, invalid in product(valid, (True, False, 1.0, "1")):
            with self.subTest(name=name, invalid=invalid):
                with self.assertRaisesRegex(TypeError, name):
                    coefficients.evaluate(**(valid | {name: invalid}))


class SGBLConstraintExactTests(unittest.TestCase):
    def test_coefficients_and_complete_pair_match_unredefined_red1(self) -> None:
        for point in _nontrivial_inputs():
            with self.subTest(alpha_gb=point.alpha_gb):
                coefficients = sgbl_constraint_coefficients(point)
                origin = _red1_pair(point, 0, 0)
                lambda_basis = _red1_pair(point, 1, 0)
                k_basis = _red1_pair(point, 0, 1)
                self.assertEqual(coefficients.hamiltonian_constant, origin[0])
                self.assertEqual(coefficients.momentum_constant, origin[1])
                self.assertEqual(coefficients.hamiltonian_lambda_r, lambda_basis[0] - origin[0])
                self.assertEqual(coefficients.momentum_k_r, k_basis[1] - origin[1])
                self.assertEqual(lambda_basis[1], origin[1])
                self.assertEqual(k_basis[0], origin[0])
                Lr, kr = -Q(1, 19), Q(1, 23)
                self.assertEqual(
                    coefficients.evaluate(radial_metric_derivative=Lr, angular_extrinsic_curvature_derivative=kr),
                    _red1_pair(point, Lr, kr),
                )

    def test_solved_derivatives_zero_the_independent_tensor_constraints(self) -> None:
        for point in _nontrivial_inputs():
            coefficients = sgbl_constraint_coefficients(point)
            Lr, kr = coefficients.solve()
            self.assertIs(type(Lr), Fraction)
            self.assertIs(type(kr), Fraction)
            self.assertEqual(
                coefficients.evaluate(radial_metric_derivative=Lr, angular_extrinsic_curvature_derivative=kr),
                (Q(0), Q(0)),
            )
            self.assertEqual(_red1_pair(point, Lr, kr), (Q(0), Q(0)))

    def test_other_normal_and_radial_second_jets_cancel_from_constraints(self) -> None:
        point = _nontrivial_inputs()[0]
        Lr, kr = Q(3, 29), -Q(1, 31)
        expected = sgbl_constraint_coefficients(point).evaluate(
            radial_metric_derivative=Lr, angular_extrinsic_curvature_derivative=kr,
        )
        self.assertEqual(_red1_pair(point, Lr, kr, free_second_jet=Q(17, 13)), expected)

    def test_gr_limit_has_exact_adm_coefficients_with_both_scalar_stresses(self) -> None:
        point = replace(_nontrivial_inputs()[1], alpha_gb=0)
        r, L, k = point.radius, point.radial_metric, point.angular_extrinsic_curvature
        rho = (
            (point.phi_pi**2 + point.chi_pi**2) / 2
            + (point.phi_r**2 + point.chi_r**2) / (2 * L**2)
            + point.scalar_mass**2 * point.phi**2 / 2
            + point.quartic_coupling * point.phi**4 / 4
        )
        expected = SGBLConstraintCoefficients(
            hamiltonian_constant=point.planck_mass**2 * ((1 - 1 / L**2) / r**2 - 3 * k**2) - rho,
            hamiltonian_lambda_r=2 * point.planck_mass**2 / (L**3 * r),
            momentum_constant=6 * point.planck_mass**2 * k / r - point.phi_pi * point.phi_r - point.chi_pi * point.chi_r,
            momentum_k_r=2 * point.planck_mass**2,
        )
        self.assertEqual(sgbl_constraint_coefficients(point), expected)
        self.assertEqual(
            _red1_pair(point, Q(3, 5), -Q(2, 7)),
            expected.evaluate(radial_metric_derivative=Q(3, 5), angular_extrinsic_curvature_derivative=-Q(2, 7)),
        )

    def test_flat_vacuum_controls_keep_the_planck_mass_square(self) -> None:
        for alpha in (0, -Q(1, 4), Q(2, 3)):
            point = replace(_vacuum_inputs(), radius=Q(5, 2), planck_mass=2, alpha_gb=alpha)
            coefficients = sgbl_constraint_coefficients(point)
            self.assertEqual(coefficients.jacobian_diagonal, (Q(16, 5), Q(8)))
            self.assertEqual(coefficients.solve(), (Q(0), Q(0)))
            self.assertEqual(_red1_pair(point, 0, 0), (Q(0), Q(0)))

    def test_independent_chi_stress_is_not_dropped(self) -> None:
        point = _nontrivial_inputs()[0]
        with_chi = sgbl_constraint_coefficients(point)
        without_chi = sgbl_constraint_coefficients(replace(point, chi_r=0, chi_pi=0))
        self.assertEqual(with_chi.jacobian_diagonal, without_chi.jacobian_diagonal)
        self.assertEqual(
            with_chi.hamiltonian_constant - without_chi.hamiltonian_constant,
            -(point.chi_pi**2 + point.chi_r**2 / point.radial_metric**2) / 2,
        )
        self.assertEqual(with_chi.momentum_constant - without_chi.momentum_constant, -point.chi_pi * point.chi_r)

    def test_declared_singular_witness_requires_planck_mass_two(self) -> None:
        point = replace(_vacuum_inputs(), radius=Q(5, 2), planck_mass=2, alpha_gb=-Q(1, 4), phi_r=-5)
        coefficients = sgbl_constraint_coefficients(point)
        self.assertEqual(coefficients.jacobian_diagonal, (Q(0), Q(0)))
        self.assertEqual(coefficients.jacobian_determinant, Q(0))
        with self.assertRaises(SGBLConstraintSolveStop) as stopped:
            coefficients.solve()
        self.assertEqual(stopped.exception.reason, "singular_constraint_jacobian")
        self.assertEqual(stopped.exception.zero_diagonals, ("H_L", "M_k"))
        self.assertIs(stopped.exception.coefficients, coefficients)
        self.assertEqual(_red1_pair(point, Q(3, 7), Q(2, 9)), (-Q(25, 2), Q(0)))
        mass_one = sgbl_constraint_coefficients(replace(point, planck_mass=1))
        self.assertEqual(mass_one.jacobian_diagonal, (-Q(12, 5), -Q(6)))
        self.assertEqual(mass_one.solve(), (-Q(125, 24), Q(0)))

    def test_each_zero_diagonal_refuses_even_compatible_underdetermined_rows(self) -> None:
        controls = (
            (dict(phi_r=Q(1, 4)), (Q(0), -Q(2)), (Q(159, 32), -Q(2)), ("H_L",)),
            (dict(phi_r=Q(1, 4), phi_rr=Q(159, 256)), (Q(0), -Q(2)), (Q(0), -Q(2)), ("H_L",)),
            (dict(phi_r=Q(1, 8)), (Q(1), Q(0)), (Q(127, 128), Q(2)), ("M_k",)),
            (dict(phi_r=Q(1, 8), phi_pi_r=Q(1, 4)), (Q(1), Q(0)), (Q(127, 128), Q(0)), ("M_k",)),
        )
        for changes, diagonal, constants, zero_names in controls:
            with self.subTest(changes=changes):
                point = replace(_vacuum_inputs(), angular_extrinsic_curvature=1, **changes)
                coefficients = sgbl_constraint_coefficients(point)
                self.assertEqual(coefficients.jacobian_diagonal, diagonal)
                self.assertEqual((coefficients.hamiltonian_constant, coefficients.momentum_constant), constants)
                with self.assertRaises(SGBLConstraintSolveStop) as stopped:
                    coefficients.solve()
                self.assertEqual(stopped.exception.zero_diagonals, zero_names)

    def test_arbitrarily_small_nonzero_diagonals_have_no_invented_floor(self) -> None:
        for N in (1, 10**40):
            point = replace(_vacuum_inputs(), phi_pi=1, phi_r=Q(2 * N - 1, 16 * N))
            coefficients = sgbl_constraint_coefficients(point)
            self.assertEqual(coefficients.jacobian_diagonal, (Q(1, N), Q(1, N)))
            Lr, kr = coefficients.solve()
            self.assertEqual(Lr, Q(N, 2) + Q((2 * N - 1)**2, 512 * N))
            self.assertEqual(kr, Q(2 * N - 1, 16))
            self.assertEqual(
                coefficients.evaluate(radial_metric_derivative=Lr, angular_extrinsic_curvature_derivative=kr),
                (Q(0), Q(0)),
            )


class SGBLConstraintSymbolicTests(unittest.TestCase):
    def test_formal_algebra_adapter_retains_exact_symbolic_coefficients(self) -> None:
        x, y = _Polynomial.atom("L_r"), _Polynomial.atom("k_r")
        polynomial = (x + y)**2 - x**2 - 2 * x * y
        self.assertEqual(polynomial, y**2)
        self.assertEqual(polynomial.substitute({"k_r": Q(2, 3)}), Q(4, 9))
        self.assertEqual((x / 3 + 2 * y + Q(4, 7)).coefficient(1, 0), Q(1, 3))
        self.assertEqual((x / 3 + 2 * y + Q(4, 7)).coefficient(0, 0), Q(4, 7))
        self.assertEqual((x + 1).substitute({"L_r": 2 * y}), 2 * y + 1)

    def test_double_dual_and_complete_polynomial_support_are_symbolic_identities(self) -> None:
        # No production coefficient formula is used to assemble this tensor.
        A, B, C, D, E, X, Y, Z, W = map(_Polynomial.atom, "A B C D E X Y Z W".split())
        riemann = [[[[_Polynomial() for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
        for a, b, c, d, value in (
            (0, 1, 0, 1, D), (0, 2, 0, 2, E), (0, 3, 0, 3, E),
            (1, 2, 1, 2, A), (1, 3, 1, 3, A), (2, 3, 2, 3, B),
            (0, 2, 1, 2, C), (0, 3, 1, 3, C),
        ):
            for i, j, sign1 in ((a, b, 1), (b, a, -1)):
                for m, n, sign2 in ((c, d, 1), (d, c, -1)):
                    riemann[i][j][m][n] = sign1 * sign2 * value
                    riemann[m][n][i][j] = sign1 * sign2 * value
        metric = [[Q(0) if a != b else Q(-1 if a == 0 else 1) for b in range(4)] for a in range(4)]
        ricci = _ricci(riemann, metric)
        scalar = _scalar(ricci, metric)
        double_dual = _p(riemann, ricci, scalar, metric)
        einstein = [[ricci[a][b] - metric[a][b] * scalar / 2 for b in range(4)] for a in range(4)]
        for indices, expected in (
            ((0, 1, 0, 1), -B), ((0, 2, 0, 2), -A), ((0, 3, 0, 3), -A),
            ((0, 1, 1, 0), B), ((0, 2, 1, 2), -C), ((0, 3, 1, 3), -C),
        ):
            a, b, c, d = indices
            self.assertEqual(double_dual[a][b][c][d], expected)
        self.assertEqual(einstein[0][0], 2 * A + B)
        self.assertEqual(einstein[0][1], 2 * C)
        hessian = [[W, Z, 0, 0], [Z, X, 0, 0], [0, 0, Y, 0], [0, 0, 0, Y]]

        def contraction(a: int, b: int):
            return sum(
                double_dual[a][c][b][d] * metric[c][e] * metric[d][f] * hessian[e][f]
                for c, d, e, f in product(range(4), repeat=4)
            )

        gb00, gb01 = contraction(0, 0), contraction(0, 1)
        self.assertEqual(gb00, -B * X - 2 * A * Y)
        self.assertEqual(gb01, -B * Z - 2 * C * Y)
        # Independent inverse-radius/metric atoms cover all r>0,L>0 without
        # symbolic rational-function division or a sample-based degree bound.
        ir, iL, k, Lr, kr, phi, pr, prr, pi, pir, xr, xpi, mpl, mu, g4, alpha = map(
            _Polynomial.atom,
            "inverse_r inverse_L k L_r k_r phi phi_r phi_rr Pi_phi Pi_phi_r chi_r Pi_chi Mpl mu g4 alpha_gb".split(),
        )
        components = {
            "A": Lr * iL**3 * ir - 2 * k**2,
            "B": k**2 + (1 - iL**2) * ir**2,
            "X": prr * iL**2 - Lr * pr * iL**3 - 2 * k * pi,
            "Y": k * pi + pr * iL**2 * ir,
        }
        rho = (pi**2 + xpi**2) / 2 + (pr**2 + xr**2) * iL**2 / 2 + mu**2 * phi**2 / 2 + g4 * phi**4 / 4
        H = (mpl**2 * einstein[0][0] - rho + 8 * alpha * gb00).substitute(components)
        # The already proved G01/GB01 are homogeneous of degree one in C,Z.
        # Thus M=L*E_orth01 is obtained by substituting L*C and L*Z; the
        # scalar coordinate stress is Pi_phi*phi_r + Pi_chi*chi_r.
        coordinate_components = components | {"C": kr + 3 * k * ir, "Z": pir - 2 * k * pr}
        M = (mpl**2 * einstein[0][1] + 8 * alpha * gb01).substitute(coordinate_components) - pi * pr - xpi * xr
        for expression, support in ((H, {(0, 0), (1, 0)}), (M, {(0, 0), (0, 1)})):
            self.assertEqual({
                (monomial.count("L_r"), monomial.count("k_r"))
                for monomial, _ in expression.formal.terms
            }, support)
            self.assertFalse(any(
                atom in ("D", "E", "W")
                for monomial, _ in expression.formal.terms for atom in monomial
            ))
        symbolic_coefficients = (
            H.coefficient(0, 0), H.coefficient(1, 0),
            M.coefficient(0, 0), M.coefficient(0, 1),
        )
        self.assertEqual(
            symbolic_coefficients[1],
            2 * mpl**2 * iL**3 * ir + 8 * alpha * iL**3 * (
                k**2 * pr - 2 * k * pi * ir + pr * ir**2 - 3 * pr * iL**2 * ir**2
            ),
        )
        self.assertEqual(symbolic_coefficients[3], 2 * mpl**2 - 16 * alpha * k * pi - 16 * alpha * pr * iL**2 * ir)
        for point in _nontrivial_inputs():
            substitution = dict(zip(
                "inverse_r inverse_L k phi phi_r phi_rr Pi_phi Pi_phi_r chi_r Pi_chi Mpl mu g4 alpha_gb".split(),
                (1 / point.radius, 1 / point.radial_metric, point.angular_extrinsic_curvature,
                 point.phi, point.phi_r, point.phi_rr, point.phi_pi, point.phi_pi_r,
                 point.chi_r, point.chi_pi, point.planck_mass, point.scalar_mass,
                 point.quartic_coupling, point.alpha_gb),
                strict=True,
            ))
            coefficients = sgbl_constraint_coefficients(point)
            for symbolic, field in zip(symbolic_coefficients, fields(coefficients), strict=True):
                self.assertEqual(symbolic.substitute(substitution), getattr(coefficients, field.name))


if __name__ == "__main__":
    unittest.main()
