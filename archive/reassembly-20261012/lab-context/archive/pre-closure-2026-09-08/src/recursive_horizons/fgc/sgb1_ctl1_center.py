"""SGB-L initial-centre profile and series controls.

The empty Minkowski buffer of the matched family is proved here as an
*initial* regular centre: even primitives, elementary flatness
``A=lambda`` through two time derivatives, and exact vanishing Laurent
residuals.  That is not a qualification of the evolving centre after
matter arrives.

Pointwise RED1 evaluation requires ``r>0``.  Substituting a tiny radius,
or borrowing the source algebra's frozen annulus minimum ``1/2``, is not
a centre chart.  Neutral ``LaurentSeries`` and ``SeriesJet2`` instruments
are reused; the FGC-QR-only profile validator and its certificate are
not.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping, Sequence

from .regular_center import (
    CERTIFIED_MAX_POWER,
    CERTIFIED_MIN_POWER,
    LaurentSeries,
    RegularSeriesState,
    SeriesJet2,
)
from .spherical_reduction import (
    Jet2,
    residuals,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
PRIMITIVE_ORDER = ("alpha", "lambda", "A", "V", "phi", "chi")
TIME_LAYER_ORDER = ("value", "dt", "dtt")
REGULAR_EQUATION_ORDER = (
    "metric_tt",
    "metric_tr_over_r",
    "metric_rr",
    "metric_angular_over_r_squared",
    "scalar_phi",
    "scalar_chi",
)
FIRST_GRID_RADII = (Q(1, 16), Q(1, 32), Q(1, 64))
EVEN_JET_MINIMUM_ERROR_RATIO = Q(15, 4)
FORMAL_SERIES_IDENTITY_KEYS = (
    "profile_id",
    "regular_equation_order",
    "center_limits",
    "regular_equation_series",
    "curvature_center",
)
DECLARED_ACTION = {
    "branch": "SGB-L",
    "planck_mass": Q(2),
    "beta": Q(0),
    "mu": Q(3),
    "g4": Q(1, 2),
    "alpha": -Q(1, 4),
    "eta": Q(0),
}


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _coefficient_layer(name: str, value: object) -> dict[int, Fraction]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a power-to-rational mapping")
    output: dict[int, Fraction] = {}
    for power, coefficient in value.items():
        if isinstance(power, bool) or not isinstance(power, int) or power < 0:
            raise ValueError(f"{name} powers must be nonnegative integers")
        if power > CERTIFIED_MAX_POWER + 4:
            raise ValueError(f"{name} power exceeds the declared profile window")
        exact = _fraction(f"{name}[{power}]", coefficient)
        if exact:
            output[power] = exact
    return output


def validate_sgbl_initial_center_profile(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Validate even parity, elementary flatness, and the declared SGB-L action."""

    if not isinstance(profile, Mapping):
        raise TypeError("SGB-L initial-centre profile must be a mapping")
    expected = {"profile_id", "fields", "action"}
    if set(profile) != expected:
        raise ValueError("SGB-L initial-centre profile keys differ")
    profile_id = profile["profile_id"]
    if not isinstance(profile_id, str) or not profile_id:
        raise ValueError("SGB-L initial-centre profile requires a nonempty identifier")
    fields = profile["fields"]
    if not isinstance(fields, Mapping) or set(fields) != set(PRIMITIVE_ORDER):
        raise ValueError("SGB-L initial-centre primitive field set differs")
    parsed: dict[str, dict[str, dict[int, Fraction]]] = {}
    for field in PRIMITIVE_ORDER:
        layers = fields[field]
        if not isinstance(layers, Mapping) or set(layers) != set(TIME_LAYER_ORDER):
            raise ValueError(f"{field} must define value, dt, and dtt layers")
        parsed[field] = {}
        for layer in TIME_LAYER_ORDER:
            coefficients = _coefficient_layer(f"fields.{field}.{layer}", layers[layer])
            if any(power % 2 for power in coefficients):
                raise ValueError(f"{field}.{layer} violates even centre parity")
            parsed[field][layer] = coefficients
    for field in ("alpha", "lambda", "A"):
        if parsed[field]["value"].get(0, Q(0)) <= 0:
            raise ValueError(f"{field}(0) must be strictly positive")
    for layer in TIME_LAYER_ORDER:
        if parsed["A"][layer].get(0, Q(0)) != parsed["lambda"][layer].get(0, Q(0)):
            raise ValueError(f"elementary flatness A=lambda fails in {layer} layer")
    action = profile["action"]
    required_action = {"branch", "planck_mass", "beta", "mu", "g4", "alpha", "eta"}
    if not isinstance(action, Mapping) or set(action) != required_action:
        raise ValueError("SGB-L initial-centre action keys differ")
    exact_action = {
        key: _fraction(f"action.{key}", action[key])
        for key in ("planck_mass", "beta", "mu", "g4", "alpha", "eta")
    }
    branch = action["branch"]
    if branch != "SGB-L" or exact_action != {key: DECLARED_ACTION[key] for key in exact_action}:
        raise ValueError("initial-centre control must use the declared linear SGB-L action")
    return {
        "profile_id": profile_id,
        "fields": parsed,
        "action": {"branch": branch, **exact_action},
        "parity_valid": True,
        "elementary_flatness_value_dt_dtt_valid": True,
        "evolving_center_after_matter_arrives": False,
    }


def sgbl_empty_minkowski_center_profile() -> dict[str, Any]:
    """Exact empty initial-centre buffer: Minkowski with vanishing scalars."""

    def unit_layers() -> dict[str, dict[int, Fraction]]:
        return {"value": {0: Q(1)}, "dt": {0: Q(0)}, "dtt": {0: Q(0)}}

    def zero_layers() -> dict[str, dict[int, Fraction]]:
        return {layer: {0: Q(0)} for layer in TIME_LAYER_ORDER}

    return {
        "profile_id": "sgbl_empty_minkowski_initial_center",
        "fields": {
            "alpha": unit_layers(),
            "lambda": unit_layers(),
            "A": unit_layers(),
            "V": zero_layers(),
            "phi": zero_layers(),
            "chi": zero_layers(),
        },
        "action": dict(DECLARED_ACTION),
    }


def sgbl_nontrivial_regular_center_profile() -> dict[str, Any]:
    """Even Minkowski-background scalars used only as a regularity control."""

    return {
        "profile_id": "sgbl_nontrivial_regular_initial_center",
        "fields": {
            "alpha": {"value": {0: Q(1)}, "dt": {}, "dtt": {}},
            "lambda": {"value": {0: Q(1)}, "dt": {}, "dtt": {}},
            "A": {"value": {0: Q(1)}, "dt": {}, "dtt": {}},
            "V": {"value": {}, "dt": {}, "dtt": {}},
            "phi": {
                "value": {0: Q(1, 20), 2: Q(1, 80), 4: Q(-1, 2000)},
                "dt": {},
                "dtt": {},
            },
            "chi": {
                "value": {0: Q(-1, 25), 2: Q(1, 90), 4: Q(1, 2500)},
                "dt": {},
                "dtt": {},
            },
        },
        "action": dict(DECLARED_ACTION),
    }


def _series_layer(parsed: Mapping[str, Any], field: str, layer: str) -> LaurentSeries:
    return LaurentSeries.from_mapping(parsed["fields"][field][layer])


def _nonzero_coefficient_map(layer: Mapping[int, Fraction]) -> dict[int, Fraction]:
    return {power: coefficient for power, coefficient in layer.items() if coefficient}


def _field_is_constant(parsed: Mapping[str, Any], field: str, value: Fraction) -> bool:
    layers = parsed["fields"][field]
    expected = {} if value == 0 else {0: value}
    if _nonzero_coefficient_map(layers["value"]) != expected:
        return False
    return not any(_nonzero_coefficient_map(layers[layer]) for layer in ("dt", "dtt"))


def _field_is_identically_zero(parsed: Mapping[str, Any], field: str) -> bool:
    return _field_is_constant(parsed, field, Q(0))


def _laurent_is_zero(series: Any) -> bool:
    if isinstance(series, LaurentSeries):
        return series.is_zero()
    return series == 0


def _profile_geometry_is_empty_minkowski(parsed: Mapping[str, Any]) -> bool:
    """Emptiness is a coefficient fact, not a profile_id prefix."""

    return (
        _field_is_constant(parsed, "alpha", Q(1))
        and _field_is_constant(parsed, "lambda", Q(1))
        and _field_is_constant(parsed, "A", Q(1))
        and _field_is_identically_zero(parsed, "V")
        and _field_is_identically_zero(parsed, "phi")
        and _field_is_identically_zero(parsed, "chi")
    )


def sgbl_initial_center_series_state(
    profile: Mapping[str, Any],
) -> tuple[RegularSeriesState, dict[str, SeriesJet2]]:
    """Build the regular SGB-L series state from even primitives."""

    parsed = validate_sgbl_initial_center_profile(profile)
    primitive: dict[str, SeriesJet2] = {}
    for field in PRIMITIVE_ORDER:
        primitive[field] = SeriesJet2.from_time_layers(
            _series_layer(parsed, field, "value"),
            _series_layer(parsed, field, "dt"),
            _series_layer(parsed, field, "dtt"),
        )
    radius = SeriesJet2.from_time_layers(
        LaurentSeries.monomial(1, 1), LaurentSeries.zero(), LaurentSeries.zero()
    )
    shift = radius * primitive["V"]
    areal_radius = radius * primitive["A"]
    lam2 = primitive["lambda"] * primitive["lambda"]
    state = RegularSeriesState(
        h_tt=-(primitive["alpha"] * primitive["alpha"]) + lam2 * shift * shift,
        h_tr=lam2 * shift,
        h_rr=lam2,
        areal_radius=areal_radius,
        phi=primitive["phi"],
        chi=primitive["chi"],
        planck_mass=parsed["action"]["planck_mass"],
        beta=parsed["action"]["beta"],
        mu=parsed["action"]["mu"],
        g4=parsed["action"]["g4"],
        alpha=parsed["action"]["alpha"],
        eta=parsed["action"]["eta"],
        branch="SGB-L",
    )
    return state, primitive


def _regularized_residual_vector(unredefined: Mapping[str, Any], radius: Any) -> tuple[Any, ...]:
    metric = unredefined["metric"]
    return (
        metric[0][0],
        metric[0][1] / radius,
        metric[1][1],
        metric[2][2] / radius**2,
        unredefined["phi"],
        unredefined["chi"],
    )


def sgbl_initial_center_series(
    profile: Mapping[str, Any],
) -> dict[str, Any]:
    """Exact Laurent residual series of an SGB-L initial-centre profile.

    The empty Minkowski buffer is the family centre.  A nontrivial regular
    profile is a series/first-grid instrument, not evolving-centre evidence.
    """

    parsed = validate_sgbl_initial_center_profile(profile)
    state, _primitive = sgbl_initial_center_series_state(profile)
    unredefined = residuals(state)  # type: ignore[arg-type]
    radius = LaurentSeries.monomial(1, 1)
    regular_vector = _regularized_residual_vector(unredefined, radius)
    negative = tuple(
        series.has_negative_power(minimum=CERTIFIED_MIN_POWER)
        if hasattr(series, "has_negative_power")
        else False
        for series in regular_vector
    )
    if any(negative):
        raise ValueError("SGB-L initial-centre series has a certified negative power")
    ricci = unredefined["R"]
    gb = unredefined["GB"]
    center_limits = tuple(series.coefficient(0) for series in regular_vector)
    curvature_center = {
        "Ricci_scalar": ricci.coefficient(0),
        "Gauss_Bonnet": gb.coefficient(0),
    }
    residuals_empty = all(_laurent_is_zero(series) for series in regular_vector) and all(
        value == 0 for value in curvature_center.values()
    )
    return {
        "profile_id": parsed["profile_id"],
        "elementary_flatness_value_dt_dtt_valid": True,
        "parity_valid": True,
        "regular_equation_order": REGULAR_EQUATION_ORDER,
        "center_limits": center_limits,
        "regular_equation_series": {
            name: series.coefficient_map()
            for name, series in zip(REGULAR_EQUATION_ORDER, regular_vector)
        },
        "curvature_center": curvature_center,
        "certified_negative_powers_absent": True,
        "initial_empty_center": (
            _profile_geometry_is_empty_minkowski(parsed) and residuals_empty
        ),
        "evolving_center_after_matter_arrives": False,
        "not_an_fgc_qr_center_certificate": True,
        "nonclaims": {
            "evolving_center_after_matter_arrives_qualified": False,
            "SGBL_branch_owned_and_healthy": False,
            "source_runtime": False,
            "interval_invertibility": False,
        },
    }


def _six_exact_limits(name: str, value: object) -> tuple[Fraction, ...]:
    if not isinstance(value, (tuple, list)) or len(value) != len(REGULAR_EQUATION_ORDER):
        raise ValueError(
            f"{name} must be the complete {len(REGULAR_EQUATION_ORDER)}-equation exact vector"
        )
    return tuple(_fraction(f"{name}[{index}]", item) for index, item in enumerate(value))


def _validate_formal_series_cache(
    formal_series: Mapping[str, Any],
    computed: Mapping[str, Any],
) -> None:
    """Reject incomplete, foreign, or altered cached series; never trust a name-only cache."""

    if not isinstance(formal_series, Mapping):
        raise TypeError("formal centre series must be a mapping")
    missing = [key for key in FORMAL_SERIES_IDENTITY_KEYS if key not in formal_series]
    if missing:
        raise ValueError(
            "formal centre series is incomplete; missing "
            + ", ".join(missing)
        )
    if formal_series.get("profile_id") != computed["profile_id"]:
        raise ValueError("formal centre series belongs to a different profile")
    if tuple(formal_series.get("regular_equation_order") or ()) != REGULAR_EQUATION_ORDER:
        raise ValueError("formal centre series must declare the complete six-equation order")
    cached_limits = _six_exact_limits("center_limits", formal_series["center_limits"])
    if cached_limits != computed["center_limits"]:
        raise ValueError("formal centre series center_limits differ from the profile series")
    cached_equations = formal_series["regular_equation_series"]
    if not isinstance(cached_equations, Mapping) or set(cached_equations) != set(REGULAR_EQUATION_ORDER):
        raise ValueError("formal centre series must include every regular equation")
    if dict(cached_equations) != dict(computed["regular_equation_series"]):
        raise ValueError("formal centre series equation coefficients differ from the profile series")
    cached_curvature = formal_series["curvature_center"]
    if not isinstance(cached_curvature, Mapping) or dict(cached_curvature) != dict(computed["curvature_center"]):
        raise ValueError("formal centre series curvature_center differs from the profile series")


def sgbl_elementary_flatness_defect(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Leading conical and odd-scalar poles, without passing a failed profile."""

    fields = profile.get("fields") if isinstance(profile, Mapping) else None
    if not isinstance(fields, Mapping):
        raise TypeError("defect ledger requires profile fields")

    def coefficient(field: str, power: int) -> Fraction:
        layers = fields.get(field)
        if not isinstance(layers, Mapping):
            raise ValueError(f"defect ledger lacks {field}")
        values = layers.get("value")
        if not isinstance(values, Mapping):
            raise ValueError(f"defect ledger lacks {field}.value")
        return _fraction(f"{field}.value[{power}]", values.get(power, Q(0)))

    A0 = coefficient("A", 0)
    lambda0 = coefficient("lambda", 0)
    if A0 == 0 or lambda0 == 0:
        raise ValueError("defect ledger requires nonzero A0 and lambda0")
    phi1 = coefficient("phi", 1)
    chi1 = coefficient("chi", 1)
    return {
        "Ricci_scalar_r_minus_2_from_conical_defect": 2 * (1 - A0**2 / lambda0**2) / A0**2,
        "box_phi_r_minus_1_from_odd_scalar_term": 2 * phi1 / lambda0**2,
        "box_chi_r_minus_1_from_odd_scalar_term": 2 * chi1 / lambda0**2,
        "elementary_flatness_holds": A0 == lambda0,
        "scalar_odd_defects_absent": phi1 == 0 and chi1 == 0,
    }


def _polynomial_derivative_value(
    coefficients: Mapping[int, Fraction],
    radius: Fraction,
    derivative_order: int,
) -> Fraction:
    if derivative_order not in {0, 1, 2}:
        raise ValueError("only radial derivatives through second order are supported")
    total = Q(0)
    for power, coefficient in coefficients.items():
        if power < derivative_order:
            continue
        factor = Q(1)
        for offset in range(derivative_order):
            factor *= power - offset
        total += coefficient * factor * radius ** (power - derivative_order)
    return total


def sgbl_pointwise_center_state(
    profile: Mapping[str, Any],
    *,
    coordinate_radius: int | Fraction,
) -> Any:
    """RED1 state at a positive radius.  ``r=0`` is refused, not replaced by epsilon."""

    parsed = validate_sgbl_initial_center_profile(profile)
    radius = _fraction("coordinate_radius", coordinate_radius)
    if radius <= 0:
        raise ValueError(
            "pointwise RED1 evaluation requires r>0; the initial centre is a Laurent chart"
        )

    def jet(field: str, *, multiply_by_r: bool = False) -> Jet2:
        layers = parsed["fields"][field]
        shifted = {
            layer: (
                {power + 1: value for power, value in layers[layer].items()}
                if multiply_by_r
                else dict(layers[layer])
            )
            for layer in TIME_LAYER_ORDER
        }
        return Jet2(
            _polynomial_derivative_value(shifted["value"], radius, 0),
            _polynomial_derivative_value(shifted["dt"], radius, 0),
            _polynomial_derivative_value(shifted["value"], radius, 1),
            _polynomial_derivative_value(shifted["dtt"], radius, 0),
            _polynomial_derivative_value(shifted["dt"], radius, 1),
            _polynomial_derivative_value(shifted["value"], radius, 2),
        )

    jets = {
        "alpha": jet("alpha"),
        "shift": jet("V", multiply_by_r=True),
        "lambda": jet("lambda"),
        "areal_radius": jet("A", multiply_by_r=True),
        "phi": jet("phi"),
        "chi": jet("chi"),
    }
    action = parsed["action"]
    return state_from_generalized_adm_pg_fixture(
        {
            "model_id": "SGB-L",
            "action_parameters": {
                "planck_mass": action["planck_mass"],
                "ricci_coupling": action["beta"],
                "scalar_mass": action["mu"],
                "quartic_coupling": action["g4"],
                "linear_gb_coupling": action["alpha"],
                "quadratic_gb_coupling": action["eta"],
            },
            "state": {
                name: {
                    component: getattr(value, component)
                    for component in ("value", "dt", "dr", "dtt", "dtr", "drr")
                }
                for name, value in jets.items()
            },
        }
    )


def sgbl_pointwise_regular_equations(
    profile: Mapping[str, Any],
    *,
    coordinate_radius: int | Fraction,
) -> tuple[Fraction, ...]:
    """Regularized RED1 residual at ``r>0``.  This is not the source annulus."""

    radius = _fraction("coordinate_radius", coordinate_radius)
    state = sgbl_pointwise_center_state(profile, coordinate_radius=radius)
    unredefined = residuals(state)
    return tuple(_regularized_residual_vector(unredefined, radius))


def sgbl_initial_center_first_grid(
    profile: Mapping[str, Any],
    *,
    radii: Sequence[int | Fraction] = FIRST_GRID_RADII,
    formal_series: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """First positive-grid approach to the Laurent centre limits.

    Radii must be positive and halve exactly.  The source algebra's frozen
    annulus minimum ``1/2`` is not a licence to evaluate ``r=0`` by a small
    substitute.  This is not a PDE or time-method convergence claim.
    """

    exact_radii = tuple(_fraction(f"radii[{index}]", value) for index, value in enumerate(radii))
    if len(exact_radii) < 3 or any(value <= 0 for value in exact_radii):
        raise ValueError("centre first-grid requires at least three positive radii")
    if any(exact_radii[index + 1] * 2 != exact_radii[index] for index in range(len(exact_radii) - 1)):
        raise ValueError("centre first-grid radii must halve exactly")
    computed = sgbl_initial_center_series(profile)
    if formal_series is not None:
        _validate_formal_series_cache(formal_series, computed)
    limits = computed["center_limits"]
    if len(limits) != len(REGULAR_EQUATION_ORDER):
        raise ValueError("computed centre series must provide the complete six-equation vector")
    values = tuple(
        sgbl_pointwise_regular_equations(profile, coordinate_radius=radius)
        for radius in exact_radii
    )
    records: list[dict[str, Any]] = []
    for component, (name, limit) in enumerate(zip(REGULAR_EQUATION_ORDER, limits)):
        errors = tuple(abs(row[component] - limit) for row in values)
        ratios: list[Fraction | None] = []
        for coarse, fine in zip(errors, errors[1:]):
            ratios.append(None if fine == 0 else coarse / fine)
        passed = all(
            (fine == 0 and coarse == 0)
            or (fine > 0 and coarse / fine >= EVEN_JET_MINIMUM_ERROR_RATIO)
            for coarse, fine in zip(errors, errors[1:])
        )
        records.append(
            {
                "component": name,
                "exact_center_limit": limit,
                "values": tuple(row[component] for row in values),
                "absolute_errors": errors,
                "coarse_to_fine_error_ratios": tuple(ratios),
                "at_least_second_order_control_passed": passed,
            }
        )
    if not all(record["at_least_second_order_control_passed"] for record in records):
        raise ValueError("SGB-L initial-centre first-grid convergence control failed")
    return {
        "profile_id": computed["profile_id"],
        "radii": exact_radii,
        "halving_sequence": True,
        "minimum_required_error_ratio": EVEN_JET_MINIMUM_ERROR_RATIO,
        "records": tuple(records),
        "all_regular_equation_controls_passed": True,
        "uses_source_annulus_minimum_one_half": False,
        "not_a_PDE_or_time_method_convergence_claim": True,
        "evolving_center_after_matter_arrives": False,
    }
