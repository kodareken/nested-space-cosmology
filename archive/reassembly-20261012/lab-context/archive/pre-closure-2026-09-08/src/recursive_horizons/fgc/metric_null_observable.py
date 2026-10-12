"""Exact local metric-null and affine Raychaudhuri diagnostics.

The observable acts only on the physical spherical metric.  It checks a
supplied affine radial null generator, evaluates both null expansions, and
verifies the complete twist-free spherical Raychaudhuri identity at one
two-jet.  Analytic and synthetic controls are useful before an evolution
exists, but they are not COL1 or DEF1 evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .modified_harmonic_reference import physical_connection_data
from .spherical_reduction import Jet2, SphericalState, direct_4d_curvature


Q = Fraction
N = 4


@dataclass(frozen=True, slots=True)
class SphericalADMGeometryJet:
    """Exact two-jet of ``alpha``, shift ``v``, ``lambda``, and areal radius."""

    lapse: Jet2
    shift: Jet2
    radial_scale: Jet2
    areal_radius: Jet2

    def __post_init__(self) -> None:
        for name in ("lapse", "shift", "radial_scale", "areal_radius"):
            if not isinstance(getattr(self, name), Jet2):
                raise TypeError(f"{name} must be a Jet2")
        if self.lapse.value <= 0 or self.radial_scale.value <= 0:
            raise ValueError("lapse and radial scale must be strictly positive")
        if self.areal_radius.value <= 0:
            raise ValueError("areal radius must be strictly positive")

    def physical_state(self) -> SphericalState:
        """Construct the exact physical metric without introducing FGC sources."""

        lam2 = self.radial_scale * self.radial_scale
        zero = Jet2.constant(0)
        return SphericalState(
            h_tt=-(self.lapse * self.lapse) + lam2 * self.shift * self.shift,
            h_tr=lam2 * self.shift,
            h_rr=lam2,
            areal_radius=self.areal_radius,
            phi=zero,
            chi=zero,
            branch="GR-0",
            planck_mass=Q(1),
            beta=Q(0),
            alpha=Q(0),
            eta=Q(0),
        )


def _derivative_jet(jet: Jet2, direction: str) -> Jet2:
    if direction == "t":
        return Jet2(jet.dt, dt=jet.dtt, dr=jet.dtr)
    if direction == "r":
        return Jet2(jet.dr, dt=jet.dtr, dr=jet.drr)
    raise ValueError("direction must be t or r")


def _metric_contraction(metric, left, right) -> Fraction:
    return sum(
        metric[a][b] * left[a] * right[b]
        for a in range(N)
        for b in range(N)
    )


def _ricci_from_covariant_riemann(riemann, inverse):
    return tuple(
        tuple(
            sum(
                inverse[a][c] * riemann[a][b][c][d]
                for a in range(N)
                for c in range(N)
            )
            for d in range(N)
        )
        for b in range(N)
    )


def _classification(theta_plus: Fraction, theta_minus: Fraction) -> str:
    if theta_plus == 0 or theta_minus == 0:
        return "marginal"
    if theta_plus < 0 and theta_minus < 0:
        return "trapped"
    if theta_plus > 0 and theta_minus > 0:
        return "anti_trapped"
    return "normal"


def metric_null_raychaudhuri_point_certificate(
    geometry: SphericalADMGeometryJet,
    *,
    branch: str,
    affine_scale: Jet2,
) -> dict[str, Any]:
    """Verify an exact affine physical-null Raychaudhuri point diagnostic."""

    if not isinstance(geometry, SphericalADMGeometryJet):
        raise TypeError("geometry must be a SphericalADMGeometryJet")
    if branch not in {"outgoing", "ingoing"}:
        raise ValueError("branch must be outgoing or ingoing")
    if not isinstance(affine_scale, Jet2) or affine_scale.value <= 0:
        raise ValueError("affine scale must be a positive Jet2")
    sign = Q(1) if branch == "outgoing" else Q(-1)
    state = geometry.physical_state()
    connection = physical_connection_data(state)
    metric = connection.metric

    inverse_lapse = Q(1) / geometry.lapse
    inverse_radial = Q(1) / geometry.radial_scale
    ell_t = inverse_lapse
    ell_r = -geometry.shift * inverse_lapse + sign * inverse_radial
    k_t = affine_scale * ell_t
    k_r = affine_scale * ell_r
    zero = Jet2.constant(0)
    k_jets = (k_t, k_r, zero, zero)
    k = tuple(component.value for component in k_jets)

    ell_plus = (
        ell_t.value,
        (-geometry.shift * inverse_lapse + inverse_radial).value,
        Q(0),
        Q(0),
    )
    ell_minus = (
        ell_t.value,
        (-geometry.shift * inverse_lapse - inverse_radial).value,
        Q(0),
        Q(0),
    )
    plus_null = _metric_contraction(metric, ell_plus, ell_plus)
    minus_null = _metric_contraction(metric, ell_minus, ell_minus)
    cross = _metric_contraction(metric, ell_plus, ell_minus)
    null_residual = _metric_contraction(metric, k, k)
    if plus_null != 0 or minus_null != 0 or cross != -2 or null_residual != 0:
        raise ValueError("physical radial null-frame normalization failed")

    acceleration = []
    for upper in range(N):
        derivative = k[0] * k_jets[upper].dt + k[1] * k_jets[upper].dr
        connection_term = sum(
            connection.christoffel[upper][lower_one][lower_two]
            * k[lower_one]
            * k[lower_two]
            for lower_one in range(N)
            for lower_two in range(N)
        )
        acceleration.append(derivative + connection_term)
    if any(value != 0 for value in acceleration):
        raise ValueError("supplied radial null generator is not affine")

    radius = geometry.areal_radius
    radius_t = _derivative_jet(radius, "t")
    radius_r = _derivative_jet(radius, "r")
    plus_directional_radius = ell_t * radius_t + (
        -geometry.shift * inverse_lapse + inverse_radial
    ) * radius_r
    minus_directional_radius = ell_t * radius_t + (
        -geometry.shift * inverse_lapse - inverse_radial
    ) * radius_r
    theta_plus = (Q(2) * plus_directional_radius / radius).value
    theta_minus = (Q(2) * minus_directional_radius / radius).value

    directional_radius = k_t * radius_t + k_r * radius_r
    theta_jet = Q(2) * directional_radius / radius
    theta = theta_jet.value
    lhs = k[0] * theta_jet.dt + k[1] * theta_jet.dr

    riemann, inverse_metric, _ = direct_4d_curvature(state)
    ricci = _ricci_from_covariant_riemann(riemann, inverse_metric)
    ricci_null = sum(
        ricci[a][b] * k[a] * k[b] for a in range(N) for b in range(N)
    )
    expansion_term = -theta * theta / 2
    shear_squared = Q(0)
    twist_squared = Q(0)
    ricci_term = -ricci_null
    rhs = expansion_term - shear_squared + twist_squared + ricci_term
    residual = lhs - rhs
    if residual != 0:
        raise ValueError("exact affine Raychaudhuri identity failed")

    return {
        "classification": "exact_local_physical_metric_affine_radial_null_Raychaudhuri_diagnostic_not_COL1_or_DEF1",
        "physical_metric_only": True,
        "modified_harmonic_auxiliary_cones_used": False,
        "null_frame": {
            "ell_plus": ell_plus,
            "ell_minus": ell_minus,
            "ell_plus_null_residual": plus_null,
            "ell_minus_null_residual": minus_null,
            "cross_normalization": cross,
            "theta_plus": theta_plus,
            "theta_minus": theta_minus,
            "round_sphere_classification": _classification(theta_plus, theta_minus),
        },
        "affine_generator": {
            "branch": branch,
            "scale_value": affine_scale.value,
            "k": k,
            "null_residual": null_residual,
            "affine_residual": tuple(acceleration),
            "affine_exact": True,
        },
        "raychaudhuri": {
            "theta": theta,
            "dtheta_dlambda_from_trajectory_jet": lhs,
            "minus_half_theta_squared": expansion_term,
            "minus_shear_squared": -shear_squared,
            "twist_squared": twist_squared,
            "minus_R_ab_k_a_k_b": ricci_term,
            "complete_rhs": rhs,
            "lhs_minus_rhs": residual,
            "identity_exact": True,
            "spherical_radial_shear_specialization": True,
            "radial_null_hypersurface_twist_free_specialization": True,
            "locally_defocusing_at_point": rhs > 0,
        },
        "nonclaims": {
            "HYP1_complete": False,
            "constraint_preserving_ACT1_IBVP_proven": False,
            "COL1_evolution_solution_consumed": False,
            "finite_trapped_interval_derived": False,
            "resolution_independent_positive_Raychaudhuri_margin_derived": False,
            "metric_null_affine_defocusing_DEF1_derived": False,
            "ROB1_open_parameter_neighborhood_derived": False,
            "non_spherical_shear_robustness_derived": False,
            "finite_invariant_transition_surface_derived": False,
            "trapped_to_anti_trapped_transition_derived": False,
            "singularity_resolution_derived": False,
            "physical_wall_derived": False,
        },
    }
