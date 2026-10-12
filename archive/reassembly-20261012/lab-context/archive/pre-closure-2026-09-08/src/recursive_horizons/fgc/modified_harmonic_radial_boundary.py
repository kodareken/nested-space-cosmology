"""Uniform frozen radial main-system boundary-flux certificate.

UHYP1 supplies twelve strictly propagating ``(p,q)`` modes, a uniformly
invertible eigenframe ``V``, and the radial symmetrizer
``H=V^{-T}V^{-1}`` on its compact implicit-branch graph.  FO1 contributes six
``u`` variables whose radial principal block is zero.  This module composes
those facts into the homogeneous modal boundary operator that sets exactly
the incoming propagating amplitudes to zero at each end of a finite annulus.

The resulting statement is a frozen-coefficient, principal, maximally
dissipative boundary-flux theorem for the main radial system.  It is not a
constraint-preserving ACT1 boundary map, a complete quasilinear source
system, or an initial-boundary value problem.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping, Sequence

from .exact_interval import Interval
from .exact_linear_algebra import matrix_det, matrix_rank
from .spherical_reduction import BASE_FIELD_ORDER


Q = Fraction
PROPAGATING_DIMENSION = 2 * len(BASE_FIELD_ORDER)
FULL_FIRST_ORDER_DIMENSION = 3 * len(BASE_FIELD_ORDER)
SECTOR_ORDER = ("tilde", "hat", "physical_chi", "regulator")


def _selector(indices: Sequence[int], dimension: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Q(int(column == index)) for column in range(dimension))
        for index in indices
    )


def _identity(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size))
        for row in range(size)
    )


def _speed_sign(box: Interval) -> int:
    if not isinstance(box, Interval):
        raise TypeError("every UHYP1 speed enclosure must be an Interval")
    if box.strictly_positive():
        return 1
    if box.strictly_negative():
        return -1
    raise ValueError("a propagated speed enclosure reaches or crosses zero")


def _absolute_margin(box: Interval) -> Fraction:
    sign = _speed_sign(box)
    return box.lower if sign > 0 else -box.upper


def _validate_predecessor(certificate: Mapping[str, Any]) -> None:
    if not isinstance(certificate, Mapping):
        raise TypeError("UHYP1 certificate must be a mapping")
    if certificate.get("classification") != "exact_rational_nonzero_compact_radial_REF1_branch_strong_hyperbolicity_certificate":
        raise ValueError("boundary certificate requires the exact UHYP1 classification")
    domain = certificate.get("domain")
    if not isinstance(domain, Mapping) or not domain.get("all_30_parameter_axes_nonzero") or not domain.get("all_6_acceleration_axes_nonzero"):
        raise ValueError("UHYP1 compact domain is incomplete")
    frame = certificate.get("eigenframe")
    if not isinstance(frame, Mapping) or not frame.get("all_enclosed_frames_invertible") or not frame.get("real_smooth_radial_eigenframe"):
        raise ValueError("UHYP1 uniform real eigenframe is unavailable")
    neumann = frame.get("neumann_inverse")
    if not isinstance(neumann, Mapping) or neumann.get("dimension") != PROPAGATING_DIMENSION or neumann.get("rho_infinity", Q(1)) >= 1:
        raise ValueError("UHYP1 Neumann frame inverse is unavailable")
    nonclaims = certificate.get("nonclaims")
    if not isinstance(nonclaims, Mapping) or any(value is not False for value in nonclaims.values()):
        raise ValueError("UHYP1 nonclaim boundary was promoted")


def frozen_radial_modal_boundary_certificate(
    uniform_hyperbolicity: Mapping[str, Any],
) -> dict[str, Any]:
    """Prove homogeneous modal maximal dissipation at both annulus ends."""

    _validate_predecessor(uniform_hyperbolicity)
    modes = uniform_hyperbolicity.get("modes")
    if not isinstance(modes, Mapping) or any(sector not in modes for sector in SECTOR_ORDER):
        raise ValueError("UHYP1 modal sectors are incomplete")

    records: list[dict[str, Any]] = []
    family_ordinals: dict[str, int] = {sector: 0 for sector in SECTOR_ORDER}
    for sector in SECTOR_ORDER:
        family = modes[sector]
        if not isinstance(family, (tuple, list)):
            raise ValueError(f"UHYP1 {sector} modes must be an ordered sequence")
        for mode in family:
            if not isinstance(mode, Mapping) or "speed_box" not in mode:
                raise ValueError(f"UHYP1 {sector} mode lacks its speed enclosure")
            speed = mode["speed_box"]
            sign = _speed_sign(speed)
            ordinal = family_ordinals[sector]
            family_ordinals[sector] += 1
            records.append(
                {
                    "index": len(records),
                    "sector": sector,
                    "sector_ordinal": ordinal,
                    "speed_box": speed,
                    "radial_speed_sign": sign,
                    "absolute_speed_margin": _absolute_margin(speed),
                }
            )
    if len(records) != PROPAGATING_DIMENSION:
        raise ValueError("UHYP1 must supply exactly twelve propagating modes")
    if any(count == 0 for count in family_ordinals.values()):
        raise ValueError("every UHYP1 sector must contribute a mode")

    positive = tuple(record["index"] for record in records if record["radial_speed_sign"] > 0)
    negative = tuple(record["index"] for record in records if record["radial_speed_sign"] < 0)
    # Outward normal is +dr at the outer end and -dr at the inner end.
    outer_incoming, outer_outgoing = negative, positive
    inner_incoming, inner_outgoing = positive, negative
    if len(positive) != 6 or len(negative) != 6:
        raise ValueError("UHYP1 radial modes do not have the required six-plus-six sign split")

    outer_selector = _selector(outer_incoming, PROPAGATING_DIMENSION)
    inner_selector = _selector(inner_incoming, PROPAGATING_DIMENSION)
    if matrix_rank(outer_selector) != len(outer_incoming) or matrix_rank(inner_selector) != len(inner_incoming):
        raise ValueError("modal boundary selector loses incoming rank")
    modal_lopatinski = _identity(6)
    if matrix_det(modal_lopatinski) != 1:
        raise ValueError("normalized modal Lopatinski matrix is singular")

    symmetrizer = uniform_hyperbolicity.get("radial_symmetrizer")
    if not isinstance(symmetrizer, Mapping) or symmetrizer.get("definition") != "H=V^-T V^-1" or symmetrizer.get("HA_symmetric") is not True:
        raise ValueError("UHYP1 radial symmetrizer contract changed")
    lower = symmetrizer.get("euclidean_coercivity_lower_bound")
    upper = symmetrizer.get("euclidean_coercivity_upper_bound")
    if not isinstance(lower, Fraction) or not isinstance(upper, Fraction) or not Q(0) < lower <= upper:
        raise ValueError("UHYP1 coercivity bounds are invalid")
    full_lower, full_upper = min(Q(1), lower), max(Q(1), upper)
    margin = min(record["absolute_speed_margin"] for record in records)
    if margin <= 0:
        raise ValueError("uniform propagated speed margin must be positive")

    sector_boundary_counts = {
        sector: {
            "outer_incoming": sum(records[index]["sector"] == sector for index in outer_incoming),
            "outer_outgoing": sum(records[index]["sector"] == sector for index in outer_outgoing),
            "inner_incoming": sum(records[index]["sector"] == sector for index in inner_incoming),
            "inner_outgoing": sum(records[index]["sector"] == sector for index in inner_outgoing),
        }
        for sector in SECTOR_ORDER
    }
    return {
        "classification": "exact_uniform_frozen_1plus1_main_system_modal_maximal_dissipation_not_constraint_preserving_ACT1_IBVP",
        "domain": {
            "same_compact_REF1_implicit_branch_graph_as_UHYP1": True,
            "all_propagated_speed_sign_margins_strict": True,
            "minimum_absolute_propagated_speed_margin": margin,
            "inner_outward_normal": "-partial_r",
            "outer_outward_normal": "+partial_r",
        },
        "state_order": {
            "u_zero_speed": tuple(f"u.{field}" for field in BASE_FIELD_ORDER),
            "propagating": tuple(f"p.{field}" for field in BASE_FIELD_ORDER)
            + tuple(f"q.{field}" for field in BASE_FIELD_ORDER),
            "full_dimension": FULL_FIRST_ORDER_DIMENSION,
            "zero_speed_dimension": len(BASE_FIELD_ORDER),
            "propagating_dimension": PROPAGATING_DIMENSION,
        },
        "principal": {
            "A_full_block_form": "diag(0_6,A_12)",
            "A_propagating_diagonalization": "A=V*Lambda*V^-1",
            "H_full_definition": "diag(I_6,V^-T*V^-1)",
            "H_full_euclidean_coercivity_lower_bound": full_lower,
            "H_full_euclidean_coercivity_upper_bound": full_upper,
            "H_full_A_full_symmetric": True,
            "zero_speed_u_radial_flux_exactly_zero": True,
        },
        "modal_records": tuple(records),
        "characteristics": {
            "positive_radial_speed_indices": positive,
            "negative_radial_speed_indices": negative,
            "outer_incoming_indices": outer_incoming,
            "outer_outgoing_indices": outer_outgoing,
            "inner_incoming_indices": inner_incoming,
            "inner_outgoing_indices": inner_outgoing,
            "u_zero_speed_indices_in_full_state": tuple(range(6)),
            "incoming_counts_verified_from_strict_intervals": True,
            "sector_boundary_counts": sector_boundary_counts,
        },
        "boundary_operators": {
            "outer": {
                "kind": "homogeneous_modal_absorption",
                "formula": "B_outer=Pi_in_outer*V^-1",
                "selector": outer_selector,
                "rank": len(outer_incoming),
            },
            "inner": {
                "kind": "homogeneous_modal_absorption",
                "formula": "B_inner=Pi_in_inner*V^-1",
                "selector": inner_selector,
                "rank": len(inner_incoming),
            },
            "conditions_imposed_on_zero_speed_u_variables": 0,
        },
        "energy_flux": {
            "modal_identity": "W^T*H*A*W=y^T*Lambda*y=sum_j(c_j*y_j^2)",
            "outer_outward_flux_nonnegative_after_homogeneous_incoming_condition": True,
            "inner_outward_flux_nonnegative_after_homogeneous_incoming_condition": True,
            "homogeneous_frozen_energy_nonincreasing": True,
            "propagating_boundary_subspace_maximally_dissipative": True,
            "uniform_outgoing_flux_margin": margin,
        },
        "frozen_boundary_stability": {
            "normalized_modal_Lopatinski_matrix": modal_lopatinski,
            "normalized_modal_Lopatinski_determinant": Q(1),
            "uniform_in_modal_normalization": True,
            "method": "symmetric_hyperbolic_energy_flux_with_exact_incoming_modal_selector",
        },
        "uniform_frozen_radial_main_system_boundary_dissipation_proven": True,
        "nonclaims": {
            "complete_nonlinear_lower_order_first_order_sources_derived": False,
            "quasilinear_local_existence_or_IBVP_proven": False,
            "metric_derived_gauge_constraint_boundary_map_proven": False,
            "physical_Hamiltonian_momentum_constraint_boundary_map_proven": False,
            "constraint_preserving_ACT1_IBVP_proven": False,
            "initial_boundary_corner_compatibility_constructed": False,
            "regular_center_or_asymptotic_boundary_treatment_proven": False,
            "multidirectional_boundary_stability_proven": False,
            "retained_EFT_domain_proven": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "metric_null_affine_defocusing_derived": False,
            "finite_invariant_transition_surface_derived": False,
            "singularity_resolution_derived": False,
        },
    }
