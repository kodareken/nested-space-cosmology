"""Exact annular SGB-L constraints on the ID1 maximal polar--areal slice.

The slice has unit lapse, zero shift, R=r, K^r_r=-2k and K^theta_theta=k;
in particular R_t=-r*k, not zero.  For f(phi)=alpha_gb*phi and beta=eta=0,
RED1 gives H=E_tt and the *coordinate* projection M=E_tr (not E_orth01).
The component proof and remaining obligations live in docs/fgc-sgb1-ctl1.md.

This is finite pointwise algebra only.  It does not construct an initial-data
family, continue through r=0, integrate a radial ODE, or certify branch health.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from fractions import Fraction


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLAnnularInputs:
    """Exact slice data excluding the two derivatives to be solved.

    Every value accepts only Fraction or built-in int and is stored as a
    Fraction.  Radius, radial_metric, planck_mass, scalar_mass (mu), and
    quartic_coupling (g4) are positive, matching RED1's parameter convention.
    planck_mass is the mass, not its square.  alpha_gb=0 is allowed as an
    algebraic GR-limit control, not as a production SGB-L branch.
    """

    radius: Fraction | int
    radial_metric: Fraction | int
    angular_extrinsic_curvature: Fraction | int
    phi: Fraction | int
    phi_r: Fraction | int
    phi_rr: Fraction | int
    phi_pi: Fraction | int
    phi_pi_r: Fraction | int
    chi_r: Fraction | int
    chi_pi: Fraction | int
    planck_mass: Fraction | int
    scalar_mass: Fraction | int
    quartic_coupling: Fraction | int
    alpha_gb: Fraction | int

    def __post_init__(self) -> None:
        for field in fields(self):
            object.__setattr__(self, field.name, _fraction(field.name, getattr(self, field.name)))
        for name in ("radius", "radial_metric", "planck_mass", "scalar_mass", "quartic_coupling"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLConstraintCoefficients:
    """H0, H_L, M0, M_k with H=H0+H_L*L_r and M=M0+M_k*k_r."""

    hamiltonian_constant: Fraction | int
    hamiltonian_lambda_r: Fraction | int
    momentum_constant: Fraction | int
    momentum_k_r: Fraction | int

    def __post_init__(self) -> None:
        for field in fields(self):
            object.__setattr__(self, field.name, _fraction(field.name, getattr(self, field.name)))

    @property
    def jacobian_diagonal(self) -> tuple[Fraction, Fraction]:
        return self.hamiltonian_lambda_r, self.momentum_k_r

    @property
    def jacobian_determinant(self) -> Fraction:
        return self.hamiltonian_lambda_r * self.momentum_k_r

    def evaluate(
        self,
        *,
        radial_metric_derivative: Fraction | int,
        angular_extrinsic_curvature_derivative: Fraction | int,
    ) -> tuple[Fraction, Fraction]:
        """Return the complete (H, M) at arbitrary exact (L_r, k_r)."""
        Lr = _fraction("radial_metric_derivative", radial_metric_derivative)
        kr = _fraction("angular_extrinsic_curvature_derivative", angular_extrinsic_curvature_derivative)
        return (
            self.hamiltonian_constant + self.hamiltonian_lambda_r * Lr,
            self.momentum_constant + self.momentum_k_r * kr,
        )

    def solve(self) -> tuple[Fraction, Fraction]:
        """Return the exact local (L_r, k_r) root, or refuse a zero diagonal.

        Nonzero is an exact algebraic condition, not a quantitative Jacobian
        floor, a neighborhood/existence theorem, or a kinetic-health bound.
        Even a compatible underdetermined zero row is refused.
        """
        if self.hamiltonian_lambda_r == 0 or self.momentum_k_r == 0:
            raise SGBLConstraintSolveStop(self)
        return (
            -self.hamiltonian_constant / self.hamiltonian_lambda_r,
            -self.momentum_constant / self.momentum_k_r,
        )


class SGBLConstraintSolveStop(ArithmeticError):
    """An exactly singular annular constraint Jacobian; no state is advanced."""

    reason = "singular_constraint_jacobian"

    def __init__(self, coefficients: SGBLConstraintCoefficients) -> None:
        self.coefficients = coefficients
        self.zero_diagonals = tuple(
            name
            for name, value in zip(("H_L", "M_k"), coefficients.jacobian_diagonal, strict=True)
            if value == 0
        )
        super().__init__(f"singular constraint Jacobian: {', '.join(self.zero_diagonals)} = 0")


def sgbl_constraint_coefficients(point: SGBLAnnularInputs) -> SGBLConstraintCoefficients:
    """Specialize the complete RED1 normal projections without fitting samples.

    In an orthonormal frame, P0101=-B, P0202=P0303=-A,
    P0110=B, P0212=P0313=-C.  Raising the mixed Hessian index and retaining
    RED1's +8*alpha_gb*P*Hess(phi) residual term gives

        H = Mpl^2*(2*A+B) - rho - 8*alpha_gb*(B*X+2*A*Y)
        M = 2*Mpl^2*L*C - Pi_phi*phi_r - Pi_chi*chi_r
            - 8*alpha_gb*L*(B*Z+2*C*Y).

    A=A0+A_L*L_r and X=X0+X_L*L_r never multiply one another; only
    C=C0+C_k*k_r contains k_r.  This proves the affine, diagonal support.
    """
    if type(point) is not SGBLAnnularInputs:
        raise TypeError("point must be SGBLAnnularInputs")
    r = point.radius
    L = point.radial_metric
    k = point.angular_extrinsic_curvature
    pr = point.phi_r
    pi = point.phi_pi
    mpl2 = point.planck_mass**2
    alpha = point.alpha_gb

    a0 = -2 * k**2
    a_l = 1 / (L**3 * r)
    b = k**2 + (1 - 1 / L**2) / r**2
    c0 = 3 * k / (L * r)
    c_k = 1 / L
    x0 = point.phi_rr / L**2 - 2 * k * pi
    x_l = -pr / L**3
    y = k * pi + pr / (L**2 * r)
    z = (point.phi_pi_r - 2 * k * pr) / L
    rho = (
        (pi**2 + point.chi_pi**2) / 2
        + (pr**2 + point.chi_r**2) / (2 * L**2)
        + point.scalar_mass**2 * point.phi**2 / 2
        + point.quartic_coupling * point.phi**4 / 4
    )
    return SGBLConstraintCoefficients(
        hamiltonian_constant=mpl2 * (2 * a0 + b) - rho - 8 * alpha * (b * x0 + 2 * a0 * y),
        hamiltonian_lambda_r=2 * mpl2 * a_l - 8 * alpha * (b * x_l + 2 * a_l * y),
        momentum_constant=(
            2 * mpl2 * L * c0 - pi * pr - point.chi_pi * point.chi_r
            - 8 * alpha * L * (b * z + 2 * c0 * y)
        ),
        momentum_k_r=2 * mpl2 * L * c_k - 16 * alpha * L * c_k * y,
    )
