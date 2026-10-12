"""Branch-owned SGB-L weak-field, topological, Schwarzschild, and established controls.

These are exact or outward-rational algebraic controls of the linear branch
``f(phi)=alpha_gb*phi``, ``beta=eta=0``.  They do not borrow FGC-QR health,
do not open holdout, and do not claim a full backreacted SGB-L solution.

The established control is the known regular decoupling-limit linear-GB
scalar on a fixed Schwarzschild metric.  Its unredefined scalar residual
reduces to minus the production potential; equivalently
``Box phi + alpha_gb GB = 0`` by outward rational algebra.  Metric residuals
remain those of the probe scalar stress.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Any, Mapping, Sequence

from .sgb1_ctl1_source import (
    SGBLSourceInputs,
    sgbl_source_residual,
    sgbl_source_state,
)
from .spherical_reduction import residuals


Q = Fraction
_ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
SCHWARZSCHILD_RADIUS = Q(3)
AREAL_RADIUS = Q(4)
LINEAR_GB_COUPLING = -Q(1, 4)
PLANCK_MASS = Q(2)
SCALAR_MASS = Q(3)
QUARTIC_COUPLING = Q(1, 2)


def _matrix_is_zero(matrix: Sequence[Sequence[object]]) -> bool:
    return all(entry == 0 for row in matrix for entry in row)


def _schwarzschild_metric_jets() -> dict[str, tuple[Fraction, Fraction, Fraction, Fraction, Fraction]]:
    """Static areal-gauge Schwarzschild jets at ``r_s=3``, ``R=r=4``."""

    return {
        "alpha": (Q(1, 2), Q(0), Q(3, 16), Q(0), Q(-21, 128)),
        "shift": (Q(0), Q(0), Q(0), Q(0), Q(0)),
        "radial_metric": (Q(2), Q(0), Q(-3, 4), Q(0), Q(39, 32)),
        "areal_radius": (AREAL_RADIUS, Q(0), Q(1), Q(0), Q(0)),
    }


def _source_point(
    *,
    phi: tuple[object, ...],
    chi: tuple[object, ...] = (0, 0, 0, 0, 0),
    alpha_gb: Fraction = LINEAR_GB_COUPLING,
    scalar_mass: Fraction = SCALAR_MASS,
    quartic_coupling: Fraction = QUARTIC_COUPLING,
    metric: str = "schwarzschild",
) -> SGBLSourceInputs:
    if metric == "minkowski":
        jets = {
            "alpha": (1, 0, 0, 0, 0),
            "shift": (0, 0, 0, 0, 0),
            "radial_metric": (1, 0, 0, 0, 0),
            "areal_radius": (Q(5, 2), 0, 1, 0, 0),
        }
        radius = Q(5, 2)
    elif metric == "schwarzschild":
        jets = _schwarzschild_metric_jets()
        radius = AREAL_RADIUS
    else:
        raise ValueError("unknown control metric")
    return SGBLSourceInputs(
        coordinate_radius=radius,
        alpha=jets["alpha"],
        shift=jets["shift"],
        radial_metric=jets["radial_metric"],
        areal_radius=jets["areal_radius"],
        phi=phi,
        chi=chi,
        planck_mass=PLANCK_MASS,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
        alpha_gb=alpha_gb,
    )


def schwarzschild_gauss_bonnet_exact(
    schwarzschild_radius: Fraction | int,
    areal_radius: Fraction | int,
) -> Fraction:
    """Exact Ricci-flat identity ``GB=12 r_s^2 / r^6``."""

    radius_s = Fraction(schwarzschild_radius)
    radius = Fraction(areal_radius)
    if radius_s <= 0 or radius <= 0:
        raise ValueError("Schwarzschild Gauss-Bonnet control requires positive radii")
    return 12 * radius_s * radius_s / radius**6


def decoupling_limit_linear_gb_scalar_jets(
    *,
    schwarzschild_radius: Fraction | int = SCHWARZSCHILD_RADIUS,
    areal_radius: Fraction | int = AREAL_RADIUS,
    alpha_gb: Fraction | int = LINEAR_GB_COUPLING,
) -> dict[str, Fraction]:
    """Regular decaying linear-GB scalar on Schwarzschild, potential omitted.

    With ``gamma=alpha_gb`` and integration constants chosen so that ``phi'``
    is finite at the horizon and ``phi->0`` at infinity,

        phi(r) = 4 gamma (1/(r_s r) + 1/(2 r^2) + r_s/(3 r^3)),
        phi'(r) = -4 gamma (1/(r_s r^2) + 1/r^3 + r_s/r^4),
        phi''(r) = 4 gamma (2/(r_s r^3) + 3/r^4 + 4 r_s/r^5).

    Then ``(1/r^2) (r^2 (1-r_s/r) phi')' = -12 gamma r_s^2 / r^6 = -gamma GB``.
    """

    radius_s = Fraction(schwarzschild_radius)
    radius = Fraction(areal_radius)
    gamma = Fraction(alpha_gb)
    value = 4 * gamma * (
        1 / (radius_s * radius) + 1 / (2 * radius**2) + radius_s / (3 * radius**3)
    )
    first = -4 * gamma * (
        1 / (radius_s * radius**2) + 1 / radius**3 + radius_s / radius**4
    )
    second = 4 * gamma * (
        2 / (radius_s * radius**3) + 3 / radius**4 + 4 * radius_s / radius**5
    )
    box = -12 * gamma * radius_s * radius_s / radius**6
    gauss_bonnet = schwarzschild_gauss_bonnet_exact(radius_s, radius)
    return {
        "value": value,
        "dr": first,
        "drr": second,
        "box_phi": box,
        "gauss_bonnet": gauss_bonnet,
        "box_phi_plus_alpha_gb_GB": box + gamma * gauss_bonnet,
        "schwarzschild_radius": radius_s,
        "areal_radius": radius,
        "alpha_gb": gamma,
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLControlRecord:
    """One branch-owned algebraic control.  Aggregate health remains false."""

    name: str
    classification: str
    local_identity_holds: bool
    payload: Mapping[str, Any]
    label: str
    not_a_full_backreacted_solution: bool

    def __post_init__(self) -> None:
        if not self.local_identity_holds:
            raise ValueError("control records are fail-closed: the identity must hold")
        object.__setattr__(self, "payload", MappingProxyType(dict(self.payload)))

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def holdout_authorized(self) -> bool:
        return False

    @property
    def copied_fgcqr_health_evidence(self) -> bool:
        return False


def sgbl_weak_field_control() -> SGBLControlRecord:
    """Minkowski: GB vanishes and the unredefined scalar is ``Box phi - V'(phi)``."""

    constant = _source_point(metric="minkowski", phi=(Q(1, 16), 0, 0, 0, 0))
    gradient = _source_point(metric="minkowski", phi=(Q(1, 16), 0, Q(1, 32), 0, 0))
    constant_state = sgbl_source_state(constant)
    gradient_state = sgbl_source_state(gradient)
    constant_res = residuals(constant_state, use_warped=True)
    gradient_res = residuals(gradient_state, use_warped=True)
    radius = Q(5, 2)
    box_gradient = 2 / radius * Q(1, 32)
    potential_constant = SCALAR_MASS**2 * Q(1, 16) + QUARTIC_COUPLING * (Q(1, 16) ** 3)
    potential_gradient = potential_constant
    if constant_res["GB"] != 0 or gradient_res["GB"] != 0:
        raise ValueError("weak-field Minkowski control lost vanishing Gauss-Bonnet")
    if constant_res["phi"] != -potential_constant:
        raise ValueError("constant weak-field scalar residual is not -V'")
    if gradient_res["phi"] != box_gradient - potential_gradient:
        raise ValueError("radial weak-field scalar residual is not Box phi - V'")
    if constant_res["chi"] != 0 or gradient_res["chi"] != 0:
        raise ValueError("weak-field chi residual must vanish")
    return SGBLControlRecord(
        name="weak_field",
        classification="minkowski_linear_scalar_and_vanishing_gb",
        local_identity_holds=True,
        label="weak-field control; GB is quadratic in curvature and vanishes on Minkowski",
        not_a_full_backreacted_solution=True,
        payload={
            "GB_constant": constant_res["GB"],
            "GB_gradient": gradient_res["GB"],
            "phi_residual_constant": constant_res["phi"],
            "phi_residual_gradient": gradient_res["phi"],
            "box_phi_gradient": box_gradient,
            "potential_prime": potential_constant,
            "metric_residual_is_minus_scalar_stress": True,
        },
    )


def sgbl_constant_coupling_topological_control() -> SGBLControlRecord:
    """Constant ``phi``: ``Hess f=0``, so the metric GB stress vanishes."""

    point = _source_point(metric="schwarzschild", phi=(Q(1, 16), 0, 0, 0, 0))
    state = sgbl_source_state(point)
    unredefined = residuals(state, use_warped=True)
    gauss_bonnet = schwarzschild_gauss_bonnet_exact(SCHWARZSCHILD_RADIUS, AREAL_RADIUS)
    potential_prime = SCALAR_MASS**2 * Q(1, 16) + QUARTIC_COUPLING * (Q(1, 16) ** 3)
    expected_phi = -potential_prime + LINEAR_GB_COUPLING * gauss_bonnet
    if unredefined["GB"] != gauss_bonnet:
        raise ValueError("topological control lost the exact Schwarzschild GB identity")
    if unredefined["R"] != 0:
        raise ValueError("topological control requires Ricci-flat Schwarzschild")
    if not _matrix_is_zero(unredefined["gb_residual_term"]):
        raise ValueError("constant coupling must make the metric GB stress vanish")
    if unredefined["phi"] != expected_phi:
        raise ValueError("constant-phi scalar residual is not -V' + alpha_gb GB")
    return SGBLControlRecord(
        name="constant_coupling_topological",
        classification="constant_phi_vanishing_metric_gb_stress",
        local_identity_holds=True,
        label=(
            "constant-coupling/topological control: f'=alpha_gb is constant and "
            "Hess phi=0, so 4d GB is topological in the metric equations"
        ),
        not_a_full_backreacted_solution=True,
        payload={
            "GB": unredefined["GB"],
            "Ricci_scalar": unredefined["R"],
            "gb_residual_term_zero": True,
            "phi_residual": unredefined["phi"],
            "expected_phi_residual": expected_phi,
            "chi_residual": unredefined["chi"],
        },
    )


def sgbl_schwarzschild_vacuum_control() -> SGBLControlRecord:
    """``phi=0`` Schwarzschild is Ricci-flat with unredefined scalar ``alpha_gb GB``."""

    point = _source_point(metric="schwarzschild", phi=(0, 0, 0, 0, 0))
    state = sgbl_source_state(point)
    unredefined = residuals(state, use_warped=True)
    gauss_bonnet = schwarzschild_gauss_bonnet_exact(SCHWARZSCHILD_RADIUS, AREAL_RADIUS)
    full = sgbl_source_residual(point, _ZERO6)
    if unredefined["R"] != 0 or unredefined["GB"] != gauss_bonnet:
        raise ValueError("Schwarzschild control lost Ricci-flat GB identity")
    if not _matrix_is_zero(unredefined["metric"]):
        raise ValueError("phi=0 Schwarzschild metric residual must vanish")
    if unredefined["phi"] != LINEAR_GB_COUPLING * gauss_bonnet:
        raise ValueError("phi=0 Schwarzschild scalar residual must be alpha_gb GB")
    if full[4] != unredefined["phi"] or full == _ZERO6:
        raise ValueError("complete MHG residual must retain the unredefined scalar source")
    return SGBLControlRecord(
        name="schwarzschild",
        classification="ricci_flat_schwarzschild_not_a_full_sgbl_solution",
        local_identity_holds=True,
        label="Schwarzschild control: metric vacuum with unredefined scalar alpha_gb*GB",
        not_a_full_backreacted_solution=True,
        payload={
            "GB": unredefined["GB"],
            "Ricci_scalar": unredefined["R"],
            "phi_residual": unredefined["phi"],
            "chi_residual": unredefined["chi"],
            "metric_residual_zero": True,
            "complete_mhg_scalar_row": full[4],
        },
    )


def sgbl_established_decoupling_control() -> SGBLControlRecord:
    """Regular decoupling-limit linear-GB scalar on fixed Schwarzschild.

    This is not a backreacted solution.  Production ``mu`` and ``g4`` remain
    in the unredefined residual as ``-V'(phi)``; the geometric identity
    ``Box phi + alpha_gb GB = 0`` is exact.
    """

    jets = decoupling_limit_linear_gb_scalar_jets()
    if jets["box_phi_plus_alpha_gb_GB"] != 0:
        raise ValueError("outward rational Box+alpha_gb GB identity failed")
    point = _source_point(
        metric="schwarzschild",
        phi=(jets["value"], 0, jets["dr"], 0, jets["drr"]),
    )
    state = sgbl_source_state(point)
    unredefined = residuals(state, use_warped=True)
    potential_prime = (
        SCALAR_MASS**2 * jets["value"] + QUARTIC_COUPLING * jets["value"] ** 3
    )
    geometric = unredefined["phi"] + potential_prime
    if unredefined["GB"] != jets["gauss_bonnet"]:
        raise ValueError("established control lost the exact Schwarzschild GB identity")
    if geometric != 0:
        raise ValueError(
            "unredefined scalar residual is not -V'; Box phi + alpha_gb GB failed"
        )
    if unredefined["chi"] != 0:
        raise ValueError("established control chi residual must vanish")
    if _matrix_is_zero(unredefined["metric"]):
        raise ValueError("probe scalar stress must source a nonzero metric residual")
    return SGBLControlRecord(
        name="established_sgb_decoupling",
        classification="decoupling_limit_linear_gb_scalar_on_schwarzschild",
        local_identity_holds=True,
        label=(
            "decoupling/established control: regular linear-GB scalar on "
            "Schwarzschild; not a full backreacted solution"
        ),
        not_a_full_backreacted_solution=True,
        payload={
            "phi": jets["value"],
            "phi_r": jets["dr"],
            "phi_rr": jets["drr"],
            "box_phi": jets["box_phi"],
            "GB": unredefined["GB"],
            "box_phi_plus_alpha_gb_GB": jets["box_phi_plus_alpha_gb_GB"],
            "unredefined_phi_residual": unredefined["phi"],
            "potential_prime": potential_prime,
            "unredefined_phi_plus_potential_prime": geometric,
            "metric_residual_zero": False,
            "chi_residual": unredefined["chi"],
            "backreacted": False,
        },
    )


def sgbl_branch_controls() -> dict[str, Any]:
    """Run every required control.  Aggregate health stays false."""

    records = (
        sgbl_weak_field_control(),
        sgbl_constant_coupling_topological_control(),
        sgbl_schwarzschild_vacuum_control(),
        sgbl_established_decoupling_control(),
    )
    return {
        "controls": {record.name: record for record in records},
        "order": tuple(record.name for record in records),
        "all_local_identities_hold": all(record.local_identity_holds for record in records),
        "SGBL_branch_owned_and_healthy": False,
        "holdout_authorized": False,
        "FRZ1": False,
        "PREF1": False,
        "established_control_is_decoupling_not_backreacted": True,
        "copied_fgcqr_health_evidence": False,
    }


__all__ = [
    "SGBLControlRecord",
    "decoupling_limit_linear_gb_scalar_jets",
    "schwarzschild_gauss_bonnet_exact",
    "sgbl_branch_controls",
    "sgbl_constant_coupling_topological_control",
    "sgbl_established_decoupling_control",
    "sgbl_schwarzschild_vacuum_control",
    "sgbl_weak_field_control",
]
