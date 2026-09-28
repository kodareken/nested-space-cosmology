#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-RED1 spherical-reduction preflight record."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "src"))

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-reduction.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-reduction.json"
SPHERICAL_REDUCTION_SOURCE = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py"
)
FGC_PUBLIC_API_SOURCE = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py"
)
HYP1_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-reduction.md"

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config


@dataclass(frozen=True, slots=True)
class HYP1ReductionConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    action_config: Path
    action_config_relative: str
    variation_config: Path
    variation_config_relative: str
    configuration: dict[str, Any]
    source_sha256: str


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _exact_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{name} keys differ: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )


def _repository_path(name: str, value: Any) -> tuple[Path, str]:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty repository-relative path")
    relative = Path(value)
    if relative.is_absolute():
        raise ValueError(f"{name} must be repository relative")
    resolved = (REPOSITORY / relative).resolve()
    try:
        canonical = resolved.relative_to(REPOSITORY)
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if canonical != relative:
        raise ValueError(f"{name} must be canonical and traversal free")
    if not resolved.is_file():
        raise ValueError(f"{name} does not exist")
    return resolved, canonical.as_posix()


def _positive_integer(name: str, value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _fraction_string(name: str, value: Any, *, nonzero: bool = False) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical exact rational string")
    try:
        fraction = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical exact rational string") from exc
    canonical = (
        str(fraction.numerator)
        if fraction.denominator == 1
        else f"{fraction.numerator}/{fraction.denominator}"
    )
    if value != canonical:
        raise ValueError(f"{name} must use canonical exact rational encoding")
    if nonzero and fraction == 0:
        raise ValueError(f"{name} must be nonzero")
    return value


def _bounded_fraction_string(
    name: str, value: Any, maximum: int, *, nonzero: bool = False
) -> str:
    encoded = _fraction_string(name, value, nonzero=nonzero)
    fraction = Fraction(encoded)
    if abs(fraction.numerator) > maximum or fraction.denominator > maximum:
        raise ValueError(f"{name} exceeds the frozen rational bound")
    return encoded


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1ReductionConfig:
    """Load the frozen HYP1-RED1 fixture with no implicit defaults."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "action_config",
            "variation_config",
            "formulation",
            "ordering",
            "fixture_contract",
            "fixtures",
        },
    )
    if raw["schema_version"] != 1:
        raise ValueError("schema_version must equal 1")
    if raw["artifact_id"] != "FGC-1-HYP1-RED1":
        raise ValueError("artifact_id must equal FGC-1-HYP1-RED1")
    if raw["project_version"] != "0.11.0":
        raise ValueError("project_version must equal 0.11.0")
    if raw["metric_signature"] != "-+++":
        raise ValueError("metric_signature must equal -+++")
    if raw["riemann_convention"] != "plus_partial_mu_gamma_nu":
        raise ValueError("riemann_convention must equal plus_partial_mu_gamma_nu")
    action_path, action_relative = _repository_path(
        "action_config", raw["action_config"]
    )
    if action_relative != "configs/fgc/fgc-1-action-gate.toml":
        raise ValueError("action_config must name the frozen ACT1 configuration")
    variation_path, variation_relative = _repository_path(
        "variation_config", raw["variation_config"]
    )
    if variation_relative != "configs/fgc/fgc-1-metric-variation.toml":
        raise ValueError("variation_config must name the frozen VAR1 configuration")
    action_config = load_action_config(action_path)
    variation_config = load_variation_config(variation_path)
    if action_config.project_version != raw["project_version"]:
        raise ValueError("ACT1 and HYP1-RED1 project versions differ")
    if action_config.metric_signature != raw["metric_signature"]:
        raise ValueError("ACT1 and HYP1-RED1 metric signatures differ")
    if variation_config.project_version != raw["project_version"]:
        raise ValueError("VAR1 and HYP1-RED1 project versions differ")
    if variation_config.metric_signature != raw["metric_signature"]:
        raise ValueError("VAR1 and HYP1-RED1 metric signatures differ")
    if variation_config.riemann_convention != raw["riemann_convention"]:
        raise ValueError("VAR1 and HYP1-RED1 Riemann conventions differ")
    if variation_config.action_config_relative != action_relative:
        raise ValueError("VAR1 must reference the same frozen ACT1 configuration")
    action_models = {model.model_id.value: model for model in action_config.models}

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "gauge_id",
            "primary_formulation",
            "crosscheck_formulation",
            "time_slices",
            "radial_gauge",
            "physical_metric",
            "equation_modification",
        },
    )
    expected_formulation = {
        "gauge_id": "generalized_radial_adm_dynamic_lambda_unfixed_areal_radius_v1",
        "primary_formulation": "generalized_radial_adm_dynamic_lambda_unfixed_areal_radius",
        "crosscheck_formulation": "general_2plus2_warped_product",
        "time_slices": "spacelike_horizon_penetrating_compatible_when_regular",
        "radial_gauge": "unfixed_in_red1",
        "physical_metric": "g_ab_matter_coupled_metric",
        "equation_modification": "none_unredefined_covariant_equations",
    }
    if dict(formulation) != expected_formulation:
        raise ValueError("formulation is not the frozen HYP1-RED1 formulation")

    ordering = _mapping("ordering", raw["ordering"])
    _exact_keys(
        "ordering",
        ordering,
        {
            "field_order",
            "independent_equation_projection_order",
            "candidate_constraint_projection_order",
            "principal_variable_order",
        },
    )
    expected_field_order = ["alpha", "shift", "lambda", "areal_radius", "phi", "chi"]
    expected_independent_projection_order = [
        "metric_tt",
        "metric_tr",
        "metric_rr",
        "metric_theta_theta",
        "scalar_phi",
        "scalar_chi",
    ]
    expected_constraint_order = ["metric_tt", "metric_tr"]
    if ordering["field_order"] != expected_field_order:
        raise ValueError("ordering.field_order is not frozen")
    if (
        ordering["independent_equation_projection_order"]
        != expected_independent_projection_order
    ):
        raise ValueError("ordering.independent_equation_projection_order is not frozen")
    if ordering["candidate_constraint_projection_order"] != expected_constraint_order:
        raise ValueError("ordering.candidate_constraint_projection_order is not frozen")
    if ordering["principal_variable_order"] != expected_field_order:
        raise ValueError("ordering.principal_variable_order is not frozen")

    contract = _mapping("fixture_contract", raw["fixture_contract"])
    _exact_keys(
        "fixture_contract",
        contract,
        {
            "rational_encoding",
            "coordinate_dimension",
            "fixture_count",
            "maximum_abs_numerator_or_denominator",
            "regular_areal_radius_required",
        },
    )
    if (
        contract["rational_encoding"]
        != "canonical_fraction_string_numerator_over_denominator"
    ):
        raise ValueError("fixture_contract.rational_encoding is not frozen")
    if contract["coordinate_dimension"] != 2:
        raise ValueError("fixture_contract.coordinate_dimension must equal 2")
    fixture_count = _positive_integer(
        "fixture_contract.fixture_count", contract["fixture_count"]
    )
    maximum = _positive_integer(
        "fixture_contract.maximum_abs_numerator_or_denominator",
        contract["maximum_abs_numerator_or_denominator"],
    )
    if contract["regular_areal_radius_required"] is not True:
        raise ValueError("fixture_contract.regular_areal_radius_required must be true")

    fixtures = raw["fixtures"]
    if not isinstance(fixtures, list) or len(fixtures) != fixture_count:
        raise ValueError("fixtures must be the frozen non-empty fixture list")
    expected_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    expected_models = ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"]
    state_keys = {
        "alpha",
        "shift",
        "lambda",
        "areal_radius",
        "phi",
        "chi",
        "alpha_t",
        "alpha_r",
        "alpha_tt",
        "alpha_tr",
        "alpha_rr",
        "shift_t",
        "shift_r",
        "shift_tt",
        "shift_tr",
        "shift_rr",
        "lambda_t",
        "lambda_r",
        "lambda_tt",
        "lambda_tr",
        "lambda_rr",
        "areal_radius_t",
        "areal_radius_r",
        "areal_radius_tt",
        "areal_radius_tr",
        "areal_radius_rr",
        "phi_t",
        "phi_r",
        "phi_tt",
        "phi_tr",
        "phi_rr",
        "chi_t",
        "chi_r",
        "chi_tt",
        "chi_tr",
        "chi_rr",
    }
    normalized_fixtures: list[dict[str, Any]] = []
    for index, fixture in enumerate(fixtures):
        fixture_mapping = _mapping(f"fixtures[{index}]", fixture)
        _exact_keys(
            f"fixtures[{index}]",
            fixture_mapping,
            {"fixture_id", "model_id", "purpose", "action_parameters", "state"},
        )
        if fixture_mapping["fixture_id"] != expected_ids[index]:
            raise ValueError("fixture identifiers are not frozen")
        if fixture_mapping["model_id"] != expected_models[index]:
            raise ValueError("fixture model order is not frozen")
        if (
            not isinstance(fixture_mapping["purpose"], str)
            or not fixture_mapping["purpose"]
        ):
            raise ValueError("fixture purpose must be a non-empty string")
        action_parameters = _mapping(
            f"fixtures[{index}].action_parameters", fixture_mapping["action_parameters"]
        )
        action_parameter_keys = {
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
            "ricci_coupling",
            "linear_gb_coupling",
            "quadratic_gb_coupling",
        }
        _exact_keys(
            f"fixtures[{index}].action_parameters",
            action_parameters,
            action_parameter_keys,
        )
        normalized_parameters = {
            key: _bounded_fraction_string(
                f"fixtures[{index}].action_parameters.{key}",
                action_parameters[key],
                maximum,
            )
            for key in sorted(action_parameter_keys)
        }
        action_model = action_models[fixture_mapping["model_id"]]
        action_parameter_names = {
            "planck_mass": "planck_mass",
            "scalar_mass": "scalar_mass",
            "quartic_coupling": "quartic_coupling",
            "ricci_coupling": "ricci_coupling",
            "linear_gb_coupling": "linear_gb_coupling",
            "quadratic_gb_coupling": "quadratic_gb_coupling",
        }
        expected_parameters = {
            key: Fraction(str(getattr(action_model, source_name)))
            for key, source_name in action_parameter_names.items()
        }
        if {
            key: Fraction(value) for key, value in normalized_parameters.items()
        } != expected_parameters:
            raise ValueError(
                "fixture action parameters must match frozen ACT1 model values"
            )
        state = _mapping(f"fixtures[{index}].state", fixture_mapping["state"])
        _exact_keys(f"fixtures[{index}].state", state, state_keys)
        normalized_state = {
            key: _bounded_fraction_string(
                f"fixtures[{index}].state.{key}", state[key], maximum
            )
            for key in sorted(state_keys)
        }
        if (
            Fraction(normalized_state["alpha"]) <= 0
            or Fraction(normalized_state["lambda"]) <= 0
        ):
            raise ValueError("fixture lapse and radial metric factor must be positive")
        if Fraction(normalized_state["areal_radius"]) <= 0:
            raise ValueError("fixture areal radius must be positive and regular")
        normalized_fixtures.append(
            {
                "fixture_id": fixture_mapping["fixture_id"],
                "model_id": fixture_mapping["model_id"],
                "purpose": fixture_mapping["purpose"],
                "action_parameters": normalized_parameters,
                "state": normalized_state,
            }
        )

    configuration = {
        "formulation": expected_formulation,
        "ordering": {
            "field_order": expected_field_order,
            "independent_equation_projection_order": expected_independent_projection_order,
            "candidate_constraint_projection_order": expected_constraint_order,
            "principal_variable_order": expected_field_order,
        },
        "fixtures": normalized_fixtures,
    }
    return HYP1ReductionConfig(
        schema_version=1,
        artifact_id="FGC-1-HYP1-RED1",
        project_version="0.11.0",
        metric_signature="-+++",
        riemann_convention="plus_partial_mu_gamma_nu",
        action_config=action_path,
        action_config_relative=action_relative,
        variation_config=variation_path,
        variation_config_relative=variation_relative,
        configuration=configuration,
        source_sha256=sha256(source).hexdigest(),
    )


def _source_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_json_value(value: Any) -> Any:
    """Convert exact core values to stable JSON without admitting floats."""

    if isinstance(value, Fraction):
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if isinstance(value, Mapping):
        return {str(key): _canonical_json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_canonical_json_value(item) for item in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise ValueError(
        f"HYP1 certificate contains unsupported JSON value {type(value).__name__}"
    )


def _required_nonclaims() -> dict[str, bool]:
    return {
        "first_order_reduction_constructed": False,
        "constraint_propagation_proven": False,
        "kinetic_matrix_invertible_on_retained_domain": False,
        "all_characteristics_real_on_retained_domain": False,
        "complete_characteristic_basis_or_symmetrizer_proven": False,
        "nonlinear_strong_hyperbolicity_proven": False,
        "full_coupled_quadratic_action_diagonalized": False,
        "nonlinear_ghost_freedom_proven": False,
        "evolution_authorized": False,
        "finite_collapse_solution_constructed": False,
        "metric_null_defocusing_derived": False,
        "singularity_resolution_proven": False,
        "child_spacetime_constructed": False,
        "shear_robustness_proven": False,
        "global_flux_entropy_closure_derived": False,
        "dark_matter_derived": False,
        "dark_energy_derived": False,
        "variable_speed_of_light_derived": False,
        "particle_spectrum_derived": False,
        "theory_of_everything_derived": False,
        "observational_confirmation_obtained": False,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, object]:
    """Build a source-bound HYP1-RED1 preflight record from the core API."""

    from recursive_horizons.fgc.spherical_reduction import (
        spherical_reduction_certificate,
    )

    config_path = config_path.resolve()
    config = load_config(config_path)
    if not SPHERICAL_REDUCTION_SOURCE.is_file() or not FGC_PUBLIC_API_SOURCE.is_file():
        raise ValueError(
            "HYP1 spherical-reduction implementation source does not exist"
        )
    if not HYP1_DOCUMENT.is_file():
        raise ValueError("HYP1 spherical-reduction owner document does not exist")
    # The core intentionally owns only the executable fixtures and identifiers.
    # Formulation and ordering are frozen and source-bound by this shell.
    certificate = spherical_reduction_certificate(
        {
            "schema_version": config.schema_version,
            "artifact_id": config.artifact_id,
            "fixtures": config.configuration["fixtures"],
        }
    )
    if not isinstance(certificate, Mapping):
        raise ValueError("HYP1 spherical-reduction certificate must be a mapping")
    aggregate = certificate.get("aggregate")
    if not isinstance(aggregate, Mapping):
        raise ValueError("HYP1 spherical-reduction aggregate is missing")
    checks = {
        "all_direct_warped_regular_jet_residuals_zero": aggregate.get(
            "all_direct_warped_exact"
        )
        is True,
        "all_core_nonclaims_false": aggregate.get("all_nonclaims_false") is True,
        "frozen_fixture_count_retained": aggregate.get("fixture_count")
        == len(config.configuration["fixtures"]),
        "all_fixture_principal_matrices_six_by_eighteen": aggregate.get(
            "all_principal_matrices_six_by_eighteen"
        )
        is True,
        "all_chi_principal_factors_metric_null": aggregate.get(
            "all_chi_principal_factors_metric_null"
        )
        is True,
        "activated_mixing_control_present_and_exercised": aggregate.get(
            "activated_mixing_control_present_and_exercised"
        )
        is True,
        "pg_schwarzschild_vacuum_metric_residual_zero": aggregate.get(
            "schwarzschild_control_present_and_exact"
        )
        is True,
        "flat_flrw_analytic_curvature_control_exact": aggregate.get(
            "flat_flrw_curvature_control_exact"
        )
        is True,
    }
    if not all(checks.values()):
        raise ValueError("HYP1 spherical-reduction exact checks did not pass")
    certificate = dict(certificate)
    certificate["verified_exact_checks"] = checks
    certificate = _canonical_json_value(certificate)
    nonclaims = _required_nonclaims()
    try:
        source_config = config_path.relative_to(REPOSITORY).as_posix()
    except ValueError:
        source_config = str(config_path)
    return {
        "schema_version": config.schema_version,
        "project_version": config.project_version,
        "artifact_id": config.artifact_id,
        "artifact": "fgc1_spherical_reduction_and_principal_part_preflight_certificate",
        "classification": "exact_spherical_reduction_and_principal_part_preflight_not_hyperbolicity_or_evolution",
        "development_status": "pre_publication_gate",
        "source_config": source_config,
        "source_configs_sha256": {
            source_config: config.source_sha256,
            config.action_config_relative: _source_digest(config.action_config),
            config.variation_config_relative: _source_digest(config.variation_config),
        },
        "implementation_sources_sha256": {
            "src/recursive_horizons/fgc/spherical_reduction.py": _source_digest(
                SPHERICAL_REDUCTION_SOURCE
            ),
            "src/recursive_horizons/fgc/__init__.py": _source_digest(
                FGC_PUBLIC_API_SOURCE
            ),
            "scripts/reproduce_fgc_hyp1_reduction.py": _source_digest(
                Path(__file__).resolve()
            ),
        },
        "derivation_document": "docs/fgc-hyp1-reduction.md",
        "derivation_document_sha256": _source_digest(HYP1_DOCUMENT),
        "formulation": config.configuration["formulation"],
        "ordering": config.configuration["ordering"],
        "fixtures": config.configuration["fixtures"],
        "reduction_certificate": dict(certificate),
        "gate_status": {
            "spherical_covariant_equations_evaluated_at_local_jets": True,
            "spherical_reduction_crosschecked_at_regular_jets": True,
            "independent_and_candidate_constraint_projection_order_frozen": True,
            "radial_constraint_projections_identified": False,
            "uneliminated_second_order_principal_part_extracted": True,
            "full_fgc1_health_gate_passed": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": nonclaims,
        "scope": "Reduces the unredefined ACT1/VAR1 covariant equations in one frozen generalized-radial ADM spherical formulation with a dynamical areal radius and extracts an exact pointwise second-order principal part on bounded rational fixtures. The chart is horizon-penetrating-compatible only where alpha, lambda, shift, areal radius, and the coordinate map remain regular; RED1 does not prove that condition. It does not construct a first-order system, prove constraint propagation, establish kinetic invertibility or real complete characteristics on a retained domain, authorize evolution, or derive collapse, defocusing, singularity resolution, a child spacetime, dark sectors, variable c, particle spectrum, or a theory of everything.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(record(arguments.config), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
