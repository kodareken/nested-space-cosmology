"""Prospective trajectory-free ``D=1-C`` Nagumo barrier obstruction.

This owner does not edit the Picard, Taylor, family, or initial-health
instruments.  It binds the exact failure of a whole-domain invariant
barrier ``{C<1}`` for the affine SGB-L constraint vector field on the
nominal ``A_chi=3`` family.

The barrier is ``D=1-C``.  After the affine H/M substitution the radial
derivative is evaluated on the surface ``D=0``, parameterized by both
signs of ``k=±1/(r lambda)``.  One outward witness is enough to reject a
trajectory-free Nagumo proof.  That rejection is not a trapped orbit and
does not reject SGB-L.

Aggregate health, ``FRZ1``, ``PREF1``, execution, and holdout remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from json import dumps
from types import MappingProxyType
from typing import Any, Mapping

from .exact_interval import Interval, interval
from .sgb1_ctl1_constraints import SGBLAnnularInputs, sgbl_constraint_coefficients
from .sgb1_ctl1_initial_health import (
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    sgbl_interval_affine_coefficients,
    sgbl_misner_sharp_compactness_interval,
    sgbl_misner_sharp_compactness_radial_derivative,
)


Q = Fraction
INSTRUMENT_ID = "FGC-1-SGB1-CTL1-TRAP-BARRIER"
CLASSIFICATION = "barrier_nonpass_not_trajectory_trapping"
NOMINAL_CHI_AMPLITUDE = Q(3)
SMALL_AMPLITUDE_CONTROL = Q(1, 8)
WITNESS_RADIUS = Q(209, 20)
WITNESS_LAMBDA = Q(1)
WITNESS_K = Q(20, 209)
OUTWARD_K_SIGN = 1
BARRIER_K_SIGNS = (1, -1)
WITNESS_D_R_STRICT_UPPER = -Q(3, 8)
NEIGHBORHOOD_D_R_STRICT_UPPER = -Q(1, 20)
NEIGHBORHOOD_RADIUS_OPEN_LEFT = Q(52, 5)
NEIGHBORHOOD_RADIUS_OPEN_RIGHT = Q(21, 2)
DECLARED_SCALAR_MASS_MAJORANT = Q(6)
AFFINE_DENOMINATOR = "H_L*M_k"
SUBSTITUTED_D_R_FORMULA = "D_r=-2*r*k^2+2*r^2*k*M0/M_k+2*H0/(H_L*lambda^3)"
CHAIN_RULE_C_R_FORMULA = "C_r=2*r*k^2+2*r^2*k*k_r+2*lambda_r/lambda^3"
BARRIER_FORMULA = "D=1-C=lambda^{-2}-r^2*k^2"
VACUUM_D_R_IDENTITY = "D_r=1/r"
NULL_EXPANSION_OUTGOING_FORMULA = "theta_+=2*(-k+1/(r*lambda))"
NULL_EXPANSION_INGOING_FORMULA = "theta_-=2*(-k-1/(r*lambda))"
MISSING_THEOREM = "whole_domain_Nagumo_invariance_of_D_equals_1_minus_C_on_the_nominal_A_chi_3_barrier"
FALSE_CLAIMS = (
    "SGBL_branch_owned_and_healthy",
    "FRZ1",
    "PREF1",
    "holdout_authorized",
    "execution_authorized",
    "whole_domain_nagumo_invariant",
    "trajectory_trapping_claimed",
    "actual_orbit_traps",
    "sgbl_model_rejected",
    "copied_gr0_or_fgcqr_health_evidence",
)
BARRIER_STOP_REASONS = frozenset(
    {
        "formula_mutation",
        "sign_mutation",
        "denominator_mutation",
        "witness_mutation",
        "family_mutation",
        "claim_mutation",
        "tampered_contract",
        "jacobian_diagonal_contains_zero",
        "domain_error",
    }
)
ZERO = Interval.singleton(0)
ONE = Interval.singleton(1)
TWO = Interval.singleton(2)


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    raise TypeError("value must be an Interval or exact rational")


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _interval_text(value: Interval) -> tuple[str, str]:
    return (_fraction_text(value.lower), _fraction_text(value.upper))


def _k_sign(value: object) -> int:
    if type(value) is not int or value not in BARRIER_K_SIGNS:
        raise SGBLTrapBarrierStop(
            "sign_mutation",
            "barrier k sign must be the exact integer +1 or -1",
            {"k_sign": value},
        )
    return value


class SGBLTrapBarrierStop(ValueError):
    """Typed contract stop of the prospective Nagumo barrier instrument."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in BARRIER_STOP_REASONS:
            raise ValueError("unknown trap-barrier stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


def sgbl_misner_sharp_barrier_value(
    radius: Interval | Fraction | int,
    radial_metric: Interval | Fraction | int,
    angular_extrinsic_curvature: Interval | Fraction | int,
) -> Interval:
    """Exact barrier ``D=1-C=lambda^{-2}-r^2 k^2``."""

    compactness = sgbl_misner_sharp_compactness_interval(
        _as_interval(radius),
        _as_interval(radial_metric),
        _as_interval(angular_extrinsic_curvature),
    )
    return ONE - compactness


def sgbl_exact_barrier_value(
    radius: Fraction | int,
    radial_metric: Fraction | int,
    angular_extrinsic_curvature: Fraction | int,
) -> Fraction:
    """Exact rational ``D=lambda^{-2}-r^2 k^2``."""

    radius = _fraction("radius", radius)
    radial_metric = _fraction("radial_metric", radial_metric)
    k = _fraction("angular_extrinsic_curvature", angular_extrinsic_curvature)
    if radius <= 0 or radial_metric <= 0:
        raise SGBLTrapBarrierStop(
            "domain_error",
            "barrier value requires positive r and lambda",
        )
    return 1 / radial_metric**2 - radius**2 * k**2


def sgbl_parameterize_barrier_k(
    radius: Interval | Fraction | int,
    radial_metric: Interval | Fraction | int,
    k_sign: int,
) -> Interval:
    """Parameterize ``D=0`` by ``k=sigma/(r lambda)`` for ``sigma=±1``."""

    sign = _k_sign(k_sign)
    radius_box = _as_interval(radius)
    lambda_box = _as_interval(radial_metric)
    if not radius_box.strictly_positive() or not lambda_box.strictly_positive():
        raise SGBLTrapBarrierStop(
            "domain_error",
            "barrier parameterization requires positive r and lambda",
        )
    return Interval.singleton(sign) / (radius_box * lambda_box)


def sgbl_exact_barrier_k(
    radius: Fraction | int,
    radial_metric: Fraction | int,
    k_sign: int,
) -> Fraction:
    """Exact rational ``k=sigma/(r lambda)`` on ``D=0``."""

    sign = _k_sign(k_sign)
    radius = _fraction("radius", radius)
    radial_metric = _fraction("radial_metric", radial_metric)
    if radius <= 0 or radial_metric <= 0:
        raise SGBLTrapBarrierStop(
            "domain_error",
            "exact barrier k requires positive r and lambda",
        )
    return Fraction(sign) / (radius * radial_metric)


def sgbl_polar_areal_null_expansions(
    radius: Interval | Fraction | int,
    radial_metric: Interval | Fraction | int,
    angular_extrinsic_curvature: Interval | Fraction | int,
) -> tuple[Interval, Interval]:
    """Unit-lapse, zero-shift, ``R=r`` future null expansions of the 2-sphere.

    ``theta_+=2(-k+1/(r lambda))`` and ``theta_-=2(-k-1/(r lambda))``.
    """

    radius_box = _as_interval(radius)
    lambda_box = _as_interval(radial_metric)
    k_box = _as_interval(angular_extrinsic_curvature)
    if not radius_box.strictly_positive() or not lambda_box.strictly_positive():
        raise SGBLTrapBarrierStop(
            "domain_error",
            "null expansions require positive r and lambda",
        )
    inverse = (radius_box * lambda_box).reciprocal()
    return TWO * (-k_box + inverse), TWO * (-k_box - inverse)


def sgbl_null_expansion_barrier_identity(
    radius: Interval | Fraction | int,
    radial_metric: Interval | Fraction | int,
    angular_extrinsic_curvature: Interval | Fraction | int,
) -> dict[str, Interval]:
    """Exact identity ``theta_+ theta_- = -(4/r^2) D``.

    Vanishing ``D`` is therefore equivalent to a vanishing null expansion.
    """

    radius_box = _as_interval(radius)
    outgoing, ingoing = sgbl_polar_areal_null_expansions(
        radius_box,
        radial_metric,
        angular_extrinsic_curvature,
    )
    barrier = sgbl_misner_sharp_barrier_value(
        radius_box, radial_metric, angular_extrinsic_curvature
    )
    product = outgoing * ingoing
    formula = Interval.singleton(-4) * barrier / (radius_box**2)
    return {
        "theta_plus": outgoing,
        "theta_minus": ingoing,
        "barrier": barrier,
        "expansion_product": product,
        "formula": formula,
        "residual": product - formula,
    }


def sgbl_affine_substituted_d_r(
    radius: Interval | Fraction | int,
    radial_metric: Interval | Fraction | int,
    angular_extrinsic_curvature: Interval | Fraction | int,
    spec: SGBLExactInitialSlice,
    *,
    formula: str = SUBSTITUTED_D_R_FORMULA,
    denominator: str = AFFINE_DENOMINATOR,
) -> dict[str, Interval]:
    """``D_r`` after ``lambda_r=-H0/H_L`` and ``k_r=-M0/M_k``.

    The substituted formula is identical to ``D_r=-C_r`` with the chain-rule
    compactness derivative.  A mutated formula or dropped ``H_L M_k``
    denominator is refused.
    """

    if formula != SUBSTITUTED_D_R_FORMULA:
        raise SGBLTrapBarrierStop(
            "formula_mutation",
            "Nagumo D_r must use the affine H/M substitution formula",
            {"formula": formula},
        )
    if denominator != AFFINE_DENOMINATOR:
        raise SGBLTrapBarrierStop(
            "denominator_mutation",
            "affine substitution denominator must be the product H_L M_k",
            {"denominator": denominator},
        )
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    radius_box = _as_interval(radius)
    lambda_box = _as_interval(radial_metric)
    k_box = _as_interval(angular_extrinsic_curvature)
    try:
        coefficients = sgbl_interval_affine_coefficients(
            radius=radius_box,
            radial_metric=lambda_box,
            angular_extrinsic_curvature=k_box,
            spec=spec,
        )
    except (SGBLInitialHealthStop, ZeroDivisionError) as exc:
        raise SGBLTrapBarrierStop(
            "domain_error",
            "affine H/M coefficients could not be enclosed on the barrier box",
            {"error": str(exc)},
        ) from exc
    hamiltonian_constant = coefficients["hamiltonian_constant"]
    hamiltonian_lambda_r = coefficients["hamiltonian_lambda_r"]
    momentum_constant = coefficients["momentum_constant"]
    momentum_k_r = coefficients["momentum_k_r"]
    if hamiltonian_lambda_r.contains_zero() or momentum_k_r.contains_zero():
        raise SGBLTrapBarrierStop(
            "jacobian_diagonal_contains_zero",
            "H_L and M_k must exclude zero on the barrier evaluation",
            {
                "hamiltonian_lambda_r": _interval_text(hamiltonian_lambda_r),
                "momentum_k_r": _interval_text(momentum_k_r),
            },
        )
    lambda_r = -hamiltonian_constant / hamiltonian_lambda_r
    k_r = -momentum_constant / momentum_k_r
    compactness_r = sgbl_misner_sharp_compactness_radial_derivative(
        radius_box, lambda_box, k_box, lambda_r, k_r
    )
    chain_rule = -compactness_r
    substituted = (
        -TWO * radius_box * (k_box**2)
        + TWO * (radius_box**2) * k_box * (momentum_constant / momentum_k_r)
        + TWO * hamiltonian_constant / (hamiltonian_lambda_r * (lambda_box**3))
    )
    if chain_rule != substituted:
        raise SGBLTrapBarrierStop(
            "formula_mutation",
            "substituted D_r does not match -C_r after affine H/M substitution",
            {
                "chain_rule": _interval_text(chain_rule),
                "substituted": _interval_text(substituted),
            },
        )
    compactness = sgbl_misner_sharp_compactness_interval(radius_box, lambda_box, k_box)
    barrier = ONE - compactness
    return {
        "compactness": compactness,
        "barrier": barrier,
        "hamiltonian_constant": hamiltonian_constant,
        "hamiltonian_lambda_r": hamiltonian_lambda_r,
        "momentum_constant": momentum_constant,
        "momentum_k_r": momentum_k_r,
        "lambda_r": lambda_r,
        "k_r": k_r,
        "compactness_r": compactness_r,
        "d_r": substituted,
    }


def sgbl_vacuum_substituted_d_r(
    radius: Fraction | int,
    radial_metric: Fraction | int,
    angular_extrinsic_curvature: Fraction | int,
) -> Fraction:
    """Exact vacuum ``D_r`` from the annular kernel with vanishing scalars."""

    radius = _fraction("radius", radius)
    radial_metric = _fraction("radial_metric", radial_metric)
    k = _fraction("angular_extrinsic_curvature", angular_extrinsic_curvature)
    coefficients = sgbl_constraint_coefficients(
        SGBLAnnularInputs(
            radius=radius,
            radial_metric=radial_metric,
            angular_extrinsic_curvature=k,
            phi=0,
            phi_r=0,
            phi_rr=0,
            phi_pi=0,
            phi_pi_r=0,
            chi_r=0,
            chi_pi=0,
            planck_mass=2,
            scalar_mass=3,
            quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        )
    )
    if coefficients.hamiltonian_lambda_r == 0 or coefficients.momentum_k_r == 0:
        raise SGBLTrapBarrierStop(
            "jacobian_diagonal_contains_zero",
            "vacuum H_L and M_k must exclude zero",
        )
    return (
        -2 * radius * k**2
        + 2
        * radius**2
        * k
        * (coefficients.momentum_constant / coefficients.momentum_k_r)
        + 2
        * coefficients.hamiltonian_constant
        / (coefficients.hamiltonian_lambda_r * radial_metric**3)
    )


def sgbl_vacuum_barrier_d_r_identity(
    radius: Fraction | int,
    radial_metric: Fraction | int,
    k_sign: int,
) -> dict[str, Fraction]:
    """On vacuum ``D=0``, the substituted derivative is exactly ``1/r``."""

    radius = _fraction("radius", radius)
    radial_metric = _fraction("radial_metric", radial_metric)
    k = sgbl_exact_barrier_k(radius, radial_metric, k_sign)
    barrier = sgbl_exact_barrier_value(radius, radial_metric, k)
    if barrier != 0:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "vacuum identity requires the D=0 parameterization",
            {"barrier": _fraction_text(barrier)},
        )
    derivative = sgbl_vacuum_substituted_d_r(radius, radial_metric, k)
    identity = 1 / radius
    if derivative != identity:
        raise SGBLTrapBarrierStop(
            "formula_mutation",
            "vacuum substituted D_r is not exactly 1/r on D=0",
            {
                "d_r": _fraction_text(derivative),
                "identity": _fraction_text(identity),
            },
        )
    return {
        "radius": radius,
        "radial_metric": radial_metric,
        "k": k,
        "k_sign": Fraction(k_sign),
        "barrier": barrier,
        "d_r": derivative,
        "identity": identity,
    }


def _require_nominal_family(spec: SGBLExactInitialSlice) -> SGBLExactInitialSlice:
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    nominal = SGBLExactInitialSlice()
    if spec != nominal:
        raise SGBLTrapBarrierStop(
            "family_mutation",
            "obstruction bind requires the frozen nominal A_chi=3 family",
            {
                "chi_amplitude": _fraction_text(spec.chi_amplitude),
                "expected": _fraction_text(NOMINAL_CHI_AMPLITUDE),
            },
        )
    return spec


def _require_false_claims(claims: Mapping[str, Any] | None) -> dict[str, bool]:
    flags = {name: False for name in FALSE_CLAIMS}
    if claims is None:
        return flags
    if not isinstance(claims, Mapping):
        raise TypeError("claims must be a mapping")
    for name, value in claims.items():
        if name not in flags:
            raise SGBLTrapBarrierStop(
                "claim_mutation",
                "unknown promoting claim is not part of the barrier contract",
                {"claim": name},
            )
        if value is not False:
            raise SGBLTrapBarrierStop(
                "claim_mutation",
                "barrier contract cannot promote aggregate or trapping claims",
                {"claim": name, "value": value},
            )
        flags[name] = False
    return flags


def _require_outward_witness(
    *,
    radius: object,
    radial_metric: object,
    k_sign: object,
    angular_extrinsic_curvature: object | None,
) -> tuple[Fraction, Fraction, Fraction, int]:
    radius = _fraction("radius", radius)
    radial_metric = _fraction("radial_metric", radial_metric)
    sign = _k_sign(k_sign)
    if sign != OUTWARD_K_SIGN:
        raise SGBLTrapBarrierStop(
            "sign_mutation",
            "outward Nagumo witness must use the positive barrier k sign",
            {"k_sign": sign},
        )
    parameterized = sgbl_exact_barrier_k(radius, radial_metric, sign)
    if angular_extrinsic_curvature is None:
        k = parameterized
    else:
        k = _fraction("angular_extrinsic_curvature", angular_extrinsic_curvature)
        if k != parameterized:
            raise SGBLTrapBarrierStop(
                "sign_mutation",
                "supplied k does not match the D=0 parameterization for this sign",
                {
                    "k": _fraction_text(k),
                    "parameterized": _fraction_text(parameterized),
                },
            )
    if radius != WITNESS_RADIUS or radial_metric != WITNESS_LAMBDA or k != WITNESS_K:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "obstruction bind requires the frozen outward witness "
            "r=209/20, lambda=1, k=20/209",
            {
                "radius": _fraction_text(radius),
                "radial_metric": _fraction_text(radial_metric),
                "k": _fraction_text(k),
            },
        )
    if sgbl_exact_barrier_value(radius, radial_metric, k) != 0:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "outward witness must lie on D=0",
        )
    return radius, radial_metric, k, sign


def _scalar_mass_potential_majorant(spec: SGBLExactInitialSlice) -> Fraction:
    """``V<=mu^2 A_phi^2/2+g4 A_phi^4/4`` because the compact bump is in ``[0,1]``."""

    phi_amp = spec.phi_amplitude
    return spec.scalar_mass**2 * phi_amp**2 / 2 + spec.quartic_coupling * phi_amp**4 / 4


def _evaluate_controls(
    spec: SGBLExactInitialSlice,
    *,
    radius: Fraction,
    radial_metric: Fraction,
    k: Fraction,
    k_sign: int,
) -> dict[str, Any]:
    opposite_k = sgbl_exact_barrier_k(radius, radial_metric, -k_sign)
    opposite = sgbl_affine_substituted_d_r(radius, radial_metric, opposite_k, spec)
    neighborhood_radius = interval(
        NEIGHBORHOOD_RADIUS_OPEN_LEFT,
        NEIGHBORHOOD_RADIUS_OPEN_RIGHT,
    )
    neighborhood_k = sgbl_parameterize_barrier_k(
        neighborhood_radius, radial_metric, k_sign
    )
    neighborhood = sgbl_affine_substituted_d_r(
        neighborhood_radius, radial_metric, neighborhood_k, spec
    )
    if not neighborhood["d_r"].upper < NEIGHBORHOOD_D_R_STRICT_UPPER:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "declared open neighborhood does not prove D_r < -1/20",
            {"neighborhood_d_r": _interval_text(neighborhood["d_r"])},
        )
    small_spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
    small = sgbl_affine_substituted_d_r(radius, radial_metric, k, small_spec)
    if not small["d_r"].strictly_positive():
        raise SGBLTrapBarrierStop(
            "family_mutation",
            "named A_chi=1/8 control must remain inward at the witness",
            {"small_amplitude_d_r": _interval_text(small["d_r"])},
        )
    dropped = sgbl_affine_substituted_d_r(
        radius,
        radial_metric,
        k,
        SGBLExactInitialSlice(phi_amplitude=0),
    )
    majorant_spec = SGBLExactInitialSlice(scalar_mass=DECLARED_SCALAR_MASS_MAJORANT)
    majorant = sgbl_affine_substituted_d_r(radius, radial_metric, k, majorant_spec)
    vacuum = {
        sign: sgbl_vacuum_barrier_d_r_identity(radius, radial_metric, sign)
        for sign in BARRIER_K_SIGNS
    }
    expansions = sgbl_polar_areal_null_expansions(radius, radial_metric, k)
    identity = sgbl_null_expansion_barrier_identity(radius, radial_metric, k)
    if identity["residual"] != ZERO:
        raise SGBLTrapBarrierStop(
            "formula_mutation",
            "null-expansion product identity missed the origin at the witness",
        )
    return {
        "opposite_k": opposite_k,
        "opposite": opposite,
        "neighborhood_radius": neighborhood_radius,
        "neighborhood_k": neighborhood_k,
        "neighborhood": neighborhood,
        "small": small,
        "dropped": dropped,
        "majorant": majorant,
        "vacuum": vacuum,
        "theta_plus": expansions[0],
        "theta_minus": expansions[1],
        "null_identity": identity,
        "mass_term_majorant": _scalar_mass_potential_majorant(spec),
    }


def _contract_payload(
    spec: SGBLExactInitialSlice,
    evidence: Mapping[str, Any],
    *,
    flags: Mapping[str, bool],
) -> dict[str, Any]:
    witness = evidence["witness"]
    controls = evidence["controls"]
    vacuum = controls["vacuum"][OUTWARD_K_SIGN]
    return {
        "INSTRUMENT_ID": INSTRUMENT_ID,
        "classification": CLASSIFICATION,
        "formula": SUBSTITUTED_D_R_FORMULA,
        "barrier": BARRIER_FORMULA,
        "denominator": AFFINE_DENOMINATOR,
        "missing_theorem": MISSING_THEOREM,
        "chi_amplitude": _fraction_text(spec.chi_amplitude),
        "witness_radius": _fraction_text(WITNESS_RADIUS),
        "witness_lambda": _fraction_text(WITNESS_LAMBDA),
        "witness_k": _fraction_text(WITNESS_K),
        "witness_k_sign": OUTWARD_K_SIGN,
        "barrier_value": list(_interval_text(witness["barrier"])),
        "compactness": list(_interval_text(witness["compactness"])),
        "hamiltonian_lambda_r": list(_interval_text(witness["hamiltonian_lambda_r"])),
        "momentum_k_r": list(_interval_text(witness["momentum_k_r"])),
        "jacobian_diagonals_exclude_zero": True,
        "d_r": list(_interval_text(witness["d_r"])),
        "d_r_strictly_below": _fraction_text(WITNESS_D_R_STRICT_UPPER),
        "k_signs": list(BARRIER_K_SIGNS),
        "opposite_sign_k": _fraction_text(controls["opposite_k"]),
        "opposite_sign_d_r": list(_interval_text(controls["opposite"]["d_r"])),
        "neighborhood_open_radius": [
            _fraction_text(NEIGHBORHOOD_RADIUS_OPEN_LEFT),
            _fraction_text(NEIGHBORHOOD_RADIUS_OPEN_RIGHT),
        ],
        "neighborhood_d_r": list(_interval_text(controls["neighborhood"]["d_r"])),
        "neighborhood_d_r_strictly_below": _fraction_text(
            NEIGHBORHOOD_D_R_STRICT_UPPER
        ),
        "vacuum_d_r": _fraction_text(vacuum["d_r"]),
        "vacuum_d_r_identity": VACUUM_D_R_IDENTITY,
        "small_amplitude": _fraction_text(SMALL_AMPLITUDE_CONTROL),
        "small_amplitude_d_r": list(_interval_text(controls["small"]["d_r"])),
        "small_amplitude_inward": True,
        "scalar_mass_majorant": _fraction_text(DECLARED_SCALAR_MASS_MAJORANT),
        "scalar_mass_majorant_d_r": list(_interval_text(controls["majorant"]["d_r"])),
        "phi_dropped_d_r": list(_interval_text(controls["dropped"]["d_r"])),
        "scalar_mass_majorant_closes": False,
        "null_expansion_outgoing": list(_interval_text(controls["theta_plus"])),
        "null_expansion_ingoing": list(_interval_text(controls["theta_minus"])),
        "null_expansion_barrier_equivalent": True,
        "sampled_nodes_are_not_the_certificate": True,
        **{name: bool(flags[name]) for name in FALSE_CLAIMS},
    }


def sgbl_trap_barrier_contract_sha256(payload: Mapping[str, Any]) -> str:
    """SHA-256 of the compact nonpromoting barrier payload."""

    raw = dumps(dict(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256(raw.encode("ascii")).hexdigest()


def _compute_barrier_evidence(
    spec: SGBLExactInitialSlice,
    *,
    radius: Fraction,
    radial_metric: Fraction,
    k: Fraction,
    k_sign: int,
) -> dict[str, Any]:
    witness = sgbl_affine_substituted_d_r(radius, radial_metric, k, spec)
    if witness["barrier"] != ZERO:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "outward witness enclosure must be the singleton D=0",
            {"barrier": _interval_text(witness["barrier"])},
        )
    if witness["compactness"] != ONE:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "outward witness enclosure must be the singleton C=1",
            {"compactness": _interval_text(witness["compactness"])},
        )
    if not witness["d_r"].upper < WITNESS_D_R_STRICT_UPPER:
        raise SGBLTrapBarrierStop(
            "witness_mutation",
            "outward witness must prove D_r < -3/8",
            {"d_r": _interval_text(witness["d_r"])},
        )
    controls = _evaluate_controls(
        spec,
        radius=radius,
        radial_metric=radial_metric,
        k=k,
        k_sign=k_sign,
    )
    if controls["majorant"]["d_r"].lower >= 0 or controls["dropped"]["d_r"].lower >= 0:
        raise SGBLTrapBarrierStop(
            "claim_mutation",
            "a scalar-mass majorant or dropped-phi evaluation cannot close the barrier",
        )
    return {"witness": witness, "controls": controls}


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLTrapBarrierRecord:
    """Exact outward Nagumo nonpass on the frozen witness and neighborhood.

    Local ``A_chi=1/8`` inwardness and vacuum ``D_r=1/r`` are controls.  They
    do not promote health, freeze, preference, holdout, or an SGB-L rejection.
    """

    slice: SGBLExactInitialSlice
    classification: str
    witness_radius: Fraction
    witness_lambda: Fraction
    witness_k: Fraction
    witness_k_sign: int
    compactness: Interval
    barrier_value: Interval
    hamiltonian_lambda_r: Interval
    momentum_k_r: Interval
    d_r: Interval
    opposite_sign_k: Fraction
    opposite_sign_d_r: Interval
    neighborhood_radius: Interval
    neighborhood_d_r: Interval
    vacuum_d_r: Fraction
    small_amplitude_d_r: Interval
    scalar_mass_majorant_d_r: Interval
    phi_dropped_d_r: Interval
    theta_plus: Interval
    theta_minus: Interval
    scalar_mass_majorant_closes: bool
    contract_payload: Mapping[str, Any]
    contract_sha256: str

    def __post_init__(self) -> None:
        spec = _require_nominal_family(self.slice)
        if self.classification != CLASSIFICATION:
            raise SGBLTrapBarrierStop(
                "claim_mutation",
                "classification must remain barrier_nonpass_not_trajectory_trapping",
            )
        if self.witness_k_sign != OUTWARD_K_SIGN:
            raise SGBLTrapBarrierStop(
                "sign_mutation",
                "record k sign must be the frozen outward sign",
            )
        if (
            self.witness_radius != WITNESS_RADIUS
            or self.witness_lambda != WITNESS_LAMBDA
            or self.witness_k != WITNESS_K
        ):
            raise SGBLTrapBarrierStop(
                "witness_mutation",
                "record witness must remain r=209/20, lambda=1, k=20/209",
            )
        if self.scalar_mass_majorant_closes is not False:
            raise SGBLTrapBarrierStop(
                "claim_mutation",
                "scalar-mass majorant cannot be recorded as closing the barrier",
            )
        evidence = _compute_barrier_evidence(
            spec,
            radius=self.witness_radius,
            radial_metric=self.witness_lambda,
            k=self.witness_k,
            k_sign=self.witness_k_sign,
        )
        witness = evidence["witness"]
        controls = evidence["controls"]
        expected = {
            "compactness": witness["compactness"],
            "barrier_value": witness["barrier"],
            "hamiltonian_lambda_r": witness["hamiltonian_lambda_r"],
            "momentum_k_r": witness["momentum_k_r"],
            "d_r": witness["d_r"],
            "opposite_sign_k": controls["opposite_k"],
            "opposite_sign_d_r": controls["opposite"]["d_r"],
            "neighborhood_radius": controls["neighborhood_radius"],
            "neighborhood_d_r": controls["neighborhood"]["d_r"],
            "vacuum_d_r": controls["vacuum"][OUTWARD_K_SIGN]["d_r"],
            "small_amplitude_d_r": controls["small"]["d_r"],
            "scalar_mass_majorant_d_r": controls["majorant"]["d_r"],
            "phi_dropped_d_r": controls["dropped"]["d_r"],
            "theta_plus": controls["theta_plus"],
            "theta_minus": controls["theta_minus"],
        }
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise SGBLTrapBarrierStop(
                    "tampered_contract",
                    f"record field {name} does not match the affine barrier evaluation",
                )
        flags = {name: False for name in FALSE_CLAIMS}
        payload = _contract_payload(spec, evidence, flags=flags)
        if dict(self.contract_payload) != payload:
            raise SGBLTrapBarrierStop(
                "tampered_contract",
                "contract payload does not match the nonpromoting barrier record",
            )
        digest = sgbl_trap_barrier_contract_sha256(payload)
        if self.contract_sha256 != digest:
            raise SGBLTrapBarrierStop(
                "tampered_contract",
                "contract hash does not match the nonpromoting payload",
            )
        if any(dict(self.contract_payload)[name] is not False for name in FALSE_CLAIMS):
            raise SGBLTrapBarrierStop(
                "claim_mutation",
                "nonpromoting contract payload cannot set aggregate or trapping claims true",
            )

    @property
    def jacobian_diagonals_exclude_zero(self) -> bool:
        return (
            not self.hamiltonian_lambda_r.contains_zero()
            and not self.momentum_k_r.contains_zero()
        )

    @property
    def outward_witness(self) -> bool:
        return self.d_r.upper < WITNESS_D_R_STRICT_UPPER

    @property
    def neighborhood_outward(self) -> bool:
        return self.neighborhood_d_r.upper < NEIGHBORHOOD_D_R_STRICT_UPPER

    @property
    def small_amplitude_inward(self) -> bool:
        return self.small_amplitude_d_r.strictly_positive()

    @property
    def null_expansion_barrier_equivalent(self) -> bool:
        return self.theta_plus == ZERO

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def whole_domain_nagumo_invariant(self) -> bool:
        return False

    @property
    def trajectory_trapping_claimed(self) -> bool:
        return False

    @property
    def actual_orbit_traps(self) -> bool:
        return False

    @property
    def sgbl_model_rejected(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False

    @property
    def holdout_authorized(self) -> bool:
        return False

    @property
    def execution_authorized(self) -> bool:
        return False

    @property
    def copied_gr0_or_fgcqr_health_evidence(self) -> bool:
        return False


def sgbl_bind_nagumo_barrier(
    spec: SGBLExactInitialSlice | None = None,
    *,
    radius: object = WITNESS_RADIUS,
    radial_metric: object = WITNESS_LAMBDA,
    k_sign: object = OUTWARD_K_SIGN,
    angular_extrinsic_curvature: object | None = None,
    formula: object = SUBSTITUTED_D_R_FORMULA,
    denominator: object = AFFINE_DENOMINATOR,
    claims: Mapping[str, Any] | None = None,
) -> SGBLTrapBarrierRecord:
    """Bind the frozen outward Nagumo nonpass, or refuse a mutation."""

    if spec is None:
        spec = SGBLExactInitialSlice()
    spec = _require_nominal_family(spec)
    if formula != SUBSTITUTED_D_R_FORMULA:
        raise SGBLTrapBarrierStop(
            "formula_mutation",
            "Nagumo D_r must use the affine H/M substitution formula",
            {"formula": formula},
        )
    if denominator != AFFINE_DENOMINATOR:
        raise SGBLTrapBarrierStop(
            "denominator_mutation",
            "affine substitution denominator must be the product H_L M_k",
            {"denominator": denominator},
        )
    flags = _require_false_claims(claims)
    radius, radial_metric, k, sign = _require_outward_witness(
        radius=radius,
        radial_metric=radial_metric,
        k_sign=k_sign,
        angular_extrinsic_curvature=angular_extrinsic_curvature,
    )
    evidence = _compute_barrier_evidence(
        spec,
        radius=radius,
        radial_metric=radial_metric,
        k=k,
        k_sign=sign,
    )
    witness = evidence["witness"]
    controls = evidence["controls"]
    payload = _contract_payload(spec, evidence, flags=flags)
    return SGBLTrapBarrierRecord(
        slice=spec,
        classification=CLASSIFICATION,
        witness_radius=radius,
        witness_lambda=radial_metric,
        witness_k=k,
        witness_k_sign=sign,
        compactness=witness["compactness"],
        barrier_value=witness["barrier"],
        hamiltonian_lambda_r=witness["hamiltonian_lambda_r"],
        momentum_k_r=witness["momentum_k_r"],
        d_r=witness["d_r"],
        opposite_sign_k=controls["opposite_k"],
        opposite_sign_d_r=controls["opposite"]["d_r"],
        neighborhood_radius=controls["neighborhood_radius"],
        neighborhood_d_r=controls["neighborhood"]["d_r"],
        vacuum_d_r=controls["vacuum"][OUTWARD_K_SIGN]["d_r"],
        small_amplitude_d_r=controls["small"]["d_r"],
        scalar_mass_majorant_d_r=controls["majorant"]["d_r"],
        phi_dropped_d_r=controls["dropped"]["d_r"],
        theta_plus=controls["theta_plus"],
        theta_minus=controls["theta_minus"],
        scalar_mass_majorant_closes=False,
        contract_payload=MappingProxyType(payload),
        contract_sha256=sgbl_trap_barrier_contract_sha256(payload),
    )


def sgbl_trap_barrier_aggregate_flags(
    record: SGBLTrapBarrierRecord,
) -> dict[str, Any]:
    """Aggregate flags stay false.  This is not a health, FRZ, or PREF gate."""

    if type(record) is not SGBLTrapBarrierRecord:
        raise TypeError("record must be SGBLTrapBarrierRecord")
    return {
        "classification": record.classification,
        "outward_witness": record.outward_witness,
        "neighborhood_outward": record.neighborhood_outward,
        "small_amplitude_inward": record.small_amplitude_inward,
        "scalar_mass_majorant_closes": record.scalar_mass_majorant_closes,
        "null_expansion_barrier_equivalent": (record.null_expansion_barrier_equivalent),
        "sampled_nodes_are_not_the_certificate": True,
        "contract_sha256": record.contract_sha256,
        **{name: False for name in FALSE_CLAIMS},
    }


__all__ = [
    "AFFINE_DENOMINATOR",
    "BARRIER_FORMULA",
    "BARRIER_K_SIGNS",
    "CLASSIFICATION",
    "DECLARED_SCALAR_MASS_MAJORANT",
    "INSTRUMENT_ID",
    "MISSING_THEOREM",
    "NEIGHBORHOOD_D_R_STRICT_UPPER",
    "NEIGHBORHOOD_RADIUS_OPEN_LEFT",
    "NEIGHBORHOOD_RADIUS_OPEN_RIGHT",
    "NOMINAL_CHI_AMPLITUDE",
    "OUTWARD_K_SIGN",
    "SMALL_AMPLITUDE_CONTROL",
    "SUBSTITUTED_D_R_FORMULA",
    "VACUUM_D_R_IDENTITY",
    "WITNESS_D_R_STRICT_UPPER",
    "WITNESS_K",
    "WITNESS_LAMBDA",
    "WITNESS_RADIUS",
    "SGBLTrapBarrierRecord",
    "SGBLTrapBarrierStop",
    "sgbl_affine_substituted_d_r",
    "sgbl_bind_nagumo_barrier",
    "sgbl_exact_barrier_k",
    "sgbl_exact_barrier_value",
    "sgbl_misner_sharp_barrier_value",
    "sgbl_null_expansion_barrier_identity",
    "sgbl_parameterize_barrier_k",
    "sgbl_polar_areal_null_expansions",
    "sgbl_trap_barrier_aggregate_flags",
    "sgbl_trap_barrier_contract_sha256",
    "sgbl_vacuum_barrier_d_r_identity",
    "sgbl_vacuum_substituted_d_r",
]
