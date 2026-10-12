"""UHYP1: a compact exact radial strong-hyperbolicity certificate.

This module is deliberately a small, fail-closed consumer of the preceding
REF1 interval-principal, MODE1, and auxiliary-identity gates.  Its domain is
the graph of the already enclosed implicit acceleration branch, not a PDE
solution space.  In particular, it makes no constraint, IBVP, EFT, or
multidirectional assertion.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence

from .exact_interval import Interval, interval
from .exact_interval_krawczyk import parametric_krawczyk_inclusion
from .exact_interval_linear_algebra import center_preconditioned_neumann_inverse
from .exact_interval_polynomial import quadratic_root_bracket
from .exact_linear_algebra import matrix_inverse
from .modified_harmonic_auxiliary_identity import auxiliary_sector_identity_certificate
from .modified_harmonic_first_order import exact_full_residual_first_order_argument_jacobian
from .modified_harmonic_interval_principal import compact_ref1_interval_principal_certificate
from .modified_harmonic_modes import auxiliary_polynomial_mode_bases, exact_comp1_acceleration_root
from .modified_harmonic_physical_modes import physical_mode_certificate, quadratic_mode_chart
from .spherical_reduction import BASE_FIELD_ORDER, SphericalState


Q = Fraction
N = len(BASE_FIELD_ORDER)


def _i(value: Fraction | int) -> Interval:
    return interval(Q(value))


def _box(value: Fraction, width: Fraction) -> Interval:
    return interval(value - width, value + width)


def _require_positive(name: str, value: Fraction) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be strictly positive")


def _center_blocks(state: SphericalState, *, reference: Any, radius: Fraction) -> dict[str, tuple[tuple[Fraction, ...], ...]]:
    return exact_full_residual_first_order_argument_jacobian(
        state, reference=reference, coordinate_radius=radius,
        tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
    )["blocks"]


def _p_box(blocks: dict[str, Sequence[Sequence[Interval]]], c: Interval) -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(
        blocks["q_r"][r][j] - c * blocks["p_r"][r][j] + c**2 * blocks["p_t"][r][j]
        for j in range(N)
    ) for r in range(N))


def _p_point(blocks: dict[str, Sequence[Sequence[Fraction]]], c: Fraction) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(blocks["q_r"][r][j] - c * blocks["p_r"][r][j] + c*c * blocks["p_t"][r][j]
                       for j in range(N)) for r in range(N))


def _dp_box(blocks: dict[str, Sequence[Sequence[Interval]]], c: Interval) -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(-blocks["p_r"][r][j] + 2*c*blocks["p_t"][r][j]
                       for j in range(N)) for r in range(N))


def _residue_at(residue: tuple[Fraction, Fraction], c: Fraction) -> Fraction:
    return residue[0] + residue[1] * c


def _interval_inverse_2(metric: tuple[Interval, Interval, Interval]) -> tuple[Interval, Interval, Interval]:
    htt, htr, hrr = metric
    determinant = htt*hrr - htr*htr
    if determinant.contains_zero():
        raise ValueError("metric-value box loses inverse")
    return hrr/determinant, -htr/determinant, htt/determinant


def _aux_inverse(gtt: Interval, gtr: Interval, grr: Interval, factor: Fraction) -> tuple[Interval, Interval, Interval]:
    # n^a n^b=-g^{a0}g^{b0}/g^{00}; all divisions are fail-closed.
    if gtt.contains_zero():
        raise ValueError("physical inverse gtt contains zero")
    ntt, ntr, nrr = -gtt*gtt/gtt, -gtt*gtr/gtt, -gtr*gtr/gtt
    return gtt-(factor-1)*ntt, gtr-(factor-1)*ntr, grr-(factor-1)*nrr


def _quadratic_from_inverse(inv: tuple[Interval, Interval, Interval]) -> tuple[Interval, Interval, Interval]:
    gtt, gtr, grr = inv
    return grr, -2*gtr, gtt


def _root_pair(coefficients: tuple[Interval, Interval, Interval], center_roots: tuple[Fraction, Fraction], width: Fraction) -> tuple[dict[str, Any], dict[str, Any]]:
    answer = tuple(quadratic_root_bracket(coefficients, _box(root, width)) for root in center_roots)
    if answer[0]["root_interval"].upper >= answer[1]["root_interval"].lower:
        raise ValueError("declared cone-root brackets overlap")
    return answer  # type: ignore[return-value]


def _disjoint(name: str, left: Interval, right: Interval) -> None:
    if not (left.upper < right.lower or right.upper < left.lower):
        raise ValueError(f"{name} brackets are not disjoint")


def _lift(c: Interval, second: Sequence[Interval]) -> tuple[Interval, ...]:
    return tuple(-c*x for x in second) + tuple(second)


def _point_lift(c: Fraction, second: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(-c*x for x in second) + tuple(second)


def _hat_mode(
    blocks_box: dict[str, Sequence[Sequence[Interval]]], blocks0: dict[str, Sequence[Sequence[Fraction]]],
    cbracket: Interval, c0: Fraction, free: int, scale: Fraction, displacement_width: Fraction,
) -> dict[str, Any]:
    rows, pivots = (1, 3, 4), (0, 1, 2)
    p0 = _p_point(blocks0, c0)
    pivot0 = tuple(tuple(p0[r][j] for j in pivots) for r in rows)
    inv = matrix_inverse(pivot0)
    free_vector = [Q(0)] * N
    free_vector[free] = scale
    x0 = tuple(-sum(inv[a][b] * p0[rows[b]][free] * scale for b in range(3)) for a in range(3))
    p = _p_box(blocks_box, cbracket)
    residual = tuple(sum((p[r][j] * free_vector[j] for j in range(N)), _i(0)) +
                     sum((p[r][pivots[j]] * x0[j] for j in range(3)), _i(0)) for r in rows)
    jac = tuple(tuple(p[r][col] for col in pivots) for r in rows)
    displacement = tuple(_box(Q(0), displacement_width) for _ in range(3))
    kraw = parametric_krawczyk_inclusion(center_inverse=inv, residual_at_center=residual,
                                         jacobian_box=jac, displacement_box=displacement)
    image = kraw["krawczyk_displacement_image_box"]
    second: list[Interval] = [_i(0) for _ in range(N)]
    for col, value, displacement_value in zip(pivots, x0, image, strict=True):
        second[col] = _i(value) + displacement_value
    second[free] = _i(scale)
    return {"free_column": free, "scale": scale, "center_speed": c0, "speed_box": cbracket,
            "center_second_order_vector": tuple(
                x0[pivots.index(i)] if i in pivots else free_vector[i] for i in range(N)),
            "second_order_vector_box": tuple(second), "linear_pivot_krawczyk": kraw,
            "first_order_vector_box": _lift(cbracket, second),
            "full_six_row_reason": "selected_rows_linear_Krawczyk_plus_scalar_rows_plus_conditional_AUX1_hat_Bianchi_completion"}


def _regulator_mode(
    blocks_box: dict[str, Sequence[Sequence[Interval]]], blocks0: dict[str, Sequence[Sequence[Fraction]]],
    seed: Sequence[Fraction], c0: Fraction, displacement_width: Fraction,
) -> dict[str, Any]:
    # x=(c,h_tt,h_tr,h_rr,R); phi=1, chi=0.  Rows 0..4 are direct full REF1 rows.
    if len(seed) != N or seed[4] != 1 or seed[5] != 0:
        raise ValueError("regulator chart must use phi=1 and chi=0")
    p0 = _p_point(blocks0, c0)
    dp0 = tuple(tuple(-blocks0["p_r"][r][j] + 2*c0*blocks0["p_t"][r][j] for j in range(N)) for r in range(N))
    jac0 = tuple(tuple(
        sum((dp0[r][j]*seed[j] for j in range(N)), Q(0)) if column == 0 else p0[r][column-1]
        for column in range(5)) for r in range(5))
    inv = matrix_inverse(jac0)
    # At nonroot c0, the Newton correction is deliberately part of the variable box.
    cbox = _box(c0, displacement_width)
    p = _p_box(blocks_box, cbox)
    dp = _dp_box(blocks_box, cbox)
    vbox = [_box(seed[j], displacement_width) for j in range(4)] + [_i(1), _i(0)]
    # Krawczyk's centre residual is F(x0,p), hence c and v are held fixed;
    # cbox/vbox occur only in D_xF(X,p), not a second time in this term.
    p_center = _p_box(blocks_box, _i(c0))
    residual = tuple(sum((p_center[r][j]*_i(seed[j]) for j in range(N)), _i(0)) for r in range(5))
    jac: list[tuple[Interval, ...]] = []
    for r in range(5):
        jac.append(tuple(
            sum((dp[r][j]*vbox[j] for j in range(N)), _i(0)) if column == 0 else p[r][column-1]
            for column in range(5)))
    displacement = tuple(_box(Q(0), displacement_width) for _ in range(5))
    kraw = parametric_krawczyk_inclusion(center_inverse=inv, residual_at_center=residual,
                                         jacobian_box=tuple(jac), displacement_box=displacement)
    image = kraw["krawczyk_displacement_image_box"]
    speed = _i(c0) + image[0]
    second = tuple(_i(seed[j]) + image[j+1] for j in range(4)) + (_i(1), _i(0))
    return {"center_seed_speed": c0, "speed_box": speed, "center_second_order_vector": tuple(seed),
            "second_order_vector_box": second, "nonlinear_full_five_row_krawczyk": kraw,
            "first_order_vector_box": _lift(speed, second),
            "sixth_row_reason": "universal_chi_principal_identity"}


def uniform_radial_hyperbolicity_certificate(
    unresolved_comp1_state: SphericalState, *, reference: Any, coordinate_radius: Fraction | int = Q(4),
    parameter_half_width: Fraction = Q(1, 2**300), acceleration_half_width: Fraction = Q(1, 2**290),
    root_bracket_half_width: Fraction = Q(1, 2**200), mode_displacement_width: Fraction = Q(1, 2**140),
    regulator_displacement_width: Fraction = Q(1, 2**68),
) -> dict[str, Any]:
    """Certify a nonzero compact radial eigenframe on the solved COMP1 graph.

    The theorem is conditional exactly on the input interval AD/Krawczyk
    enclosures and on the named action-diffeomorphism premise in AUX1.
    """
    radius = Q(coordinate_radius)
    for name, value in (("parameter_half_width", parameter_half_width), ("acceleration_half_width", acceleration_half_width),
                        ("root_bracket_half_width", root_bracket_half_width), ("mode_displacement_width", mode_displacement_width),
                        ("regulator_displacement_width", regulator_displacement_width)):
        _require_positive(name, value)
    root = exact_comp1_acceleration_root(unresolved_comp1_state, reference=reference, coordinate_radius=radius,
                                         tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
    center = root["solved_state"]
    principal = compact_ref1_interval_principal_certificate(center, reference=reference, coordinate_radius=radius,
        tilde_normal_factor=Q(4), hat_normal_factor=Q(9), parameter_half_width=parameter_half_width,
        acceleration_half_width=acceleration_half_width)
    blocks_box = principal["box"]["blocks"]
    blocks0 = _center_blocks(center, reference=reference, radius=radius)
    # Value-only metric enclosure; all thirty remaining parameter directions remain in the REF1 AD box.
    metric = (_box(center.h_tt.value, parameter_half_width), _box(center.h_tr.value, parameter_half_width),
              _box(center.h_rr.value, parameter_half_width))
    physical = _quadratic_from_inverse(_interval_inverse_2(metric))
    tilde = _quadratic_from_inverse(_aux_inverse(*_interval_inverse_2(metric), Q(4)))
    hat = _quadratic_from_inverse(_aux_inverse(*_interval_inverse_2(metric), Q(9)))
    cones = {"physical": _root_pair(physical, (Q(-1), Q(1)), root_bracket_half_width),
             "tilde": _root_pair(tilde, (Q(-1,2), Q(1,2)), root_bracket_half_width),
             "hat": _root_pair(hat, (Q(-1,3), Q(1,3)), root_bracket_half_width)}
    for family_a in cones:
        for family_b in cones:
            if family_a < family_b:
                for a, b in zip(cones[family_a], cones[family_b], strict=True):
                    _disjoint(f"{family_a}/{family_b}", a["root_interval"], b["root_interval"])
    auxiliary = auxiliary_sector_identity_certificate()
    atlas = auxiliary_polynomial_mode_bases(center, tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
    physical_atlas = physical_mode_certificate(center)
    # Tilde generators are exact polynomials; evaluate their universal kernels on each root box.
    tilde_modes = []
    for record in cones["tilde"]:
        c = record["root_interval"]
        tilde_modes.extend({"speed_box": c, "second_order_vector_box": tuple(_i(-2)*c if k == 0 else _i(1) if k == 1 else _i(0) for k in range(N)),
                            "first_order_vector_box": None, "source": "universal_tilde_pure_gauge"} if which == 0 else
                           {"speed_box": c, "second_order_vector_box": tuple(_i(0) if k == 0 else -c if k == 1 else _i(2) if k == 2 else _i(0) for k in range(N)),
                            "first_order_vector_box": None, "source": "universal_tilde_pure_gauge"}
                           for which in range(2))
    for mode in tilde_modes:
        mode["first_order_vector_box"] = _lift(mode["speed_box"], mode["second_order_vector_box"])
    hat_modes = []
    for c0, record in zip((Q(-1,3), Q(1,3)), cones["hat"], strict=True):
        hat_modes.append(_hat_mode(blocks_box, blocks0, record["root_interval"], c0, 3, Q(1), mode_displacement_width))
        hat_modes.append(_hat_mode(blocks_box, blocks0, record["root_interval"], c0, 4, Q(1,2**22), mode_displacement_width))
    physical_modes = []
    for record in cones["physical"]:
        second = tuple(_i(1) if j == 5 else _i(0) for j in range(N))
        physical_modes.append({"speed_box": record["root_interval"], "second_order_vector_box": second,
                               "first_order_vector_box": _lift(record["root_interval"], second),
                               "source": "universal_chi_identity", "full_six_row_reason": "universal_chi_principal_identity"})
    full_symbol = __import__("recursive_horizons.fgc.modified_harmonic", fromlist=["modified_harmonic_symbol"]).modified_harmonic_symbol(center, tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
    chart = quadratic_mode_chart(full_symbol, physical_atlas["regulator_factor"], required_chart=(0, 4))
    regulator_modes = []
    for c0 in (Q(-1), Q(1)):
        seed = tuple(_residue_at(value, c0) for value in chart["second_order_vector"])
        regulator_modes.append(_regulator_mode(blocks_box, blocks0, seed, c0, regulator_displacement_width))
    for mode, record in zip(regulator_modes, cones["physical"], strict=True):
        _disjoint("regulator/physical", mode["speed_box"], record["root_interval"])
        if not ((mode["speed_box"].strictly_positive() and record["root_interval"].strictly_positive()) or
                (mode["speed_box"].strictly_negative() and record["root_interval"].strictly_negative())):
            raise ValueError("regulator speed enclosure does not retain its declared sign")
    # Assemble the actual 12 interval columns in paired-sector ordering.
    columns = [m["first_order_vector_box"] for m in tilde_modes + hat_modes + physical_modes + regulator_modes]
    if len(columns) != 12:
        raise ValueError("UHYP1 requires exactly twelve propagated first-order modes")
    vbox = tuple(tuple(columns[col][row] for col in range(12)) for row in range(12))
    # Exact point preconditioner is obtained from the declared central representatives.
    point_columns = []
    for mode in tilde_modes + hat_modes + physical_modes + regulator_modes:
        c = mode.get("center_speed", mode.get("center_seed_speed"))
        if c is None:
            c = mode["speed_box"].midpoint()
        second = mode["center_second_order_vector"] if "center_second_order_vector" in mode else tuple(x.midpoint() for x in mode["second_order_vector_box"])
        point_columns.append(_point_lift(c, second))
    v0 = tuple(tuple(point_columns[col][row] for col in range(12)) for row in range(12))
    frame = center_preconditioned_neumann_inverse(vbox, matrix_inverse(v0))
    vnorm = max(sum(entry.abs_upper() for entry in row) for row in vbox)
    vinvnorm = frame["center_inverse_infinity_norm"] / (1-frame["rho_infinity"])
    coercivity_lower = Q(1, 12) / (vnorm*vnorm)
    coercivity_upper = Q(12) * vinvnorm*vinvnorm
    return {
        "classification": "exact_rational_nonzero_compact_radial_REF1_branch_strong_hyperbolicity_certificate",
        "domain": {"parameter_half_width": parameter_half_width, "acceleration_half_width": acceleration_half_width,
                   "all_30_parameter_axes_nonzero": True, "all_6_acceleration_axes_nonzero": True,
                   "branch": "unique_REF1_acceleration_root_graph_from_compact_Krawczyk"},
        "center_root": {k:v for k,v in root.items() if k != "solved_state"},
        "interval_principal": principal,
        "cone_quadratics": {"physical": physical, "tilde": tilde, "hat": hat, "root_brackets": cones},
        "universal_auxiliary_identities": auxiliary,
        "point_atlases": {"auxiliary": atlas, "physical": physical_atlas},
        "modes": {"tilde": tilde_modes, "hat": hat_modes, "physical_chi": physical_modes, "regulator": regulator_modes,
                  "correlated_eigenmode_routes": {
                    "tilde": "cone-root bracket plus universal conditional pure-gauge right identity; never interval-residual containment",
                    "hat": "three selected-row linear Krawczyk solves plus scalar rows and conditional AUX1 Bianchi row completion; never interval-residual containment",
                    "physical_chi": "metric-null root bracket plus universal chi principal identity",
                    "regulator": "five-row normalized nonlinear parametric Krawczyk inclusion plus universal chi principal identity"}},
        "eigenframe": {"interval_columns": columns, "exact_center_frame": v0, "neumann_inverse": frame,
                       "all_enclosed_frames_invertible": True, "real_smooth_radial_eigenframe": True},
        "radial_symmetrizer": {"definition": "H=V^-T V^-1", "HA_symmetric": True,
                               "euclidean_coercivity_lower_bound": coercivity_lower,
                               "euclidean_coercivity_upper_bound": coercivity_upper},
        "nonclaims": {"multidirectional_strong_hyperbolicity_proven": False, "constraint_propagation_proven": False,
                      "IBVP_proven": False, "retained_EFT_domain_proven": False, "evolution_authorized": False,
                      "collapse_solution_derived": False, "metric_null_affine_defocusing_derived": False,
                      "singularity_resolution_derived": False},
    }
