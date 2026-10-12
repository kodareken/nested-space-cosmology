#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-CON1-COMP1 compatible-local-datum certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, is_dataclass
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

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-con1-comp1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-con1-comp1.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-con1-comp1.md"

from scripts.reproduce_fgc_hyp1_dom1_qift1 import (
    IMPLEMENTATION_FILES as QIFT1_IMPLEMENTATION_FILES,
    load_canonical_result as load_qift1_result,
    load_config as load_qift1_config,
)
from recursive_horizons.fgc.exact_interval import Interval
from recursive_horizons.fgc.modified_harmonic_constraints import (
    ACCELERATION_ORDER as CON1_ACCELERATION_ORDER,
    PHYSICAL_PROJECTION_ORDER as CON1_PROJECTION_ORDER,
    compatible_constraint_certificate,
    required_con1_comp1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    SPHERICAL_FLAT_REFERENCE_ID,
    flat_spherical_annulus_reference,
)


Q = Fraction
IMPLEMENTATION_FILES = tuple(
    dict.fromkeys(
        QIFT1_IMPLEMENTATION_FILES
        + (
            REPOSITORY
            / "src"
            / "recursive_horizons"
            / "fgc"
            / "modified_harmonic_constraints.py",
            Path(__file__).resolve(),
        )
    )
)


@dataclass(frozen=True, slots=True)
class HYP1CON1COMP1Config:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    qift1_config_path: Path
    qift1_config_relative: str
    qift1_result_path: Path
    qift1_result_relative: str
    qift1_first_order_path: Path
    qift1_first_order_relative: str
    qift1_flat_fixture: dict[str, Any]
    qift1_unique_acceleration_root_passed: bool
    predecessor_source_paths: dict[str, Path]
    predecessor_source_relatives: dict[str, str]
    radial_domain_minimum: Fraction
    coordinate_radius: Fraction
    phi_value: Fraction
    phi_radial_derivative: Fraction
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
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


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a canonical rational string")
    if isinstance(value, int):
        result = Q(value)
    elif isinstance(value, str) and value:
        try:
            result = Q(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be a canonical rational string") from exc
        if value != _fraction_text(result):
            raise ValueError(f"{name} must use canonical rational encoding")
    else:
        raise ValueError(f"{name} must be an exact rational value")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _repository_path(name: str, value: Any) -> tuple[Path, str]:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a repository-relative path")
    relative = Path(value)
    if relative.is_absolute():
        raise ValueError(f"{name} must be repository relative")
    resolved = (REPOSITORY / relative).resolve()
    try:
        canonical = resolved.relative_to(REPOSITORY)
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if canonical != relative or not resolved.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return resolved, canonical.as_posix()


def _string_list(name: str, value: Any, expected: tuple[str, ...]) -> None:
    if not isinstance(value, list) or value != list(expected):
        raise ValueError(f"{name} is not the frozen ordered string list")


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1CON1COMP1Config:
    """Load the strict CON1-COMP1 schema and reconcile the QIFT1 predecessor."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _exact_keys(
        "root",
        raw,
        {
            "schema_version", "artifact_id", "project_version", "metric_signature",
            "riemann_convention", "quantified_domain_config", "quantified_domain_result",
            "reference", "activated_local_parameter_datum", "formulation",
            "polynomial_certificate", "qift1_contract", "proof_contract", "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-CON1-COMP1",
        "project_version": "0.11.0",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }
    for name, expected in expected_scalars.items():
        if raw[name] != expected:
            raise ValueError(f"{name} must equal {expected}")

    qift1_config_path, qift1_config_relative = _repository_path(
        "quantified_domain_config", raw["quantified_domain_config"]
    )
    if qift1_config_relative != "configs/fgc/fgc-1-hyp1-dom1-qift1.toml":
        raise ValueError("quantified_domain_config must name the frozen QIFT1 predecessor")
    qift1 = load_qift1_config(qift1_config_path)
    qift1_result_path, qift1_result_relative = _repository_path(
        "quantified_domain_result", raw["quantified_domain_result"]
    )
    if qift1_result_relative != "results/fgc-1-hyp1-dom1-qift1.json":
        raise ValueError("quantified_domain_result must name the frozen QIFT1 result")
    qift1_result = load_qift1_result(qift1_result_path)
    if (
        qift1_result.get("artifact_id") != "FGC-1-HYP1-DOM1-QIFT1"
        or qift1_result.get("project_version") != raw["project_version"]
        or qift1.project_version != raw["project_version"]
        or qift1.metric_signature != raw["metric_signature"]
        or qift1.riemann_convention != raw["riemann_convention"]
    ):
        raise ValueError("CON1-COMP1 and QIFT1 provenance or conventions differ")

    reference = _mapping("reference", raw["reference"])
    _exact_keys(reference_name := "reference", reference, {"reference_id", "radial_domain_minimum", "coordinate_radius", "center_included"})
    if reference["reference_id"] != SPHERICAL_FLAT_REFERENCE_ID or reference["center_included"] is not False:
        raise ValueError("CON1-COMP1 reference annulus is not frozen")
    radial_minimum = _fraction("reference.radial_domain_minimum", reference["radial_domain_minimum"], positive=True)
    coordinate_radius = _fraction("reference.coordinate_radius", reference["coordinate_radius"], positive=True)
    if radial_minimum != qift1.radial_domain_minimum or coordinate_radius != qift1.coordinate_radius:
        raise ValueError("CON1-COMP1 and QIFT1 reference data differ")

    datum = _mapping("activated_local_parameter_datum", raw["activated_local_parameter_datum"])
    _exact_keys("activated_local_parameter_datum", datum, {"fixture_id", "model_id", "phi_value", "phi_radial_derivative", "derived_second_derivatives", "datum_status"})
    if datum["fixture_id"] != "FGCQR_qift1_activated_local_parameter" or datum["model_id"] != "FGC-QR":
        raise ValueError("CON1-COMP1 activated datum identity is not frozen")
    if datum["derived_second_derivatives"] != ["h_tt.drr", "areal_radius.drr"]:
        raise ValueError("CON1-COMP1 derived second-derivative order is not frozen")
    if datum["datum_status"] != "exact_local_parameter_datum_not_an_initial_slice_or_solution_family":
        raise ValueError("CON1-COMP1 datum status is not frozen")
    phi_value = _fraction("activated_local_parameter_datum.phi_value", datum["phi_value"], positive=True)
    if phi_value != Q(1, 131072):
        raise ValueError("CON1-COMP1 phi datum is not frozen")
    phi_radial_derivative = _fraction(
        "activated_local_parameter_datum.phi_radial_derivative",
        datum["phi_radial_derivative"],
        positive=True,
    )
    if phi_radial_derivative != Q(1, 131072):
        raise ValueError("CON1-COMP1 phi radial derivative is not frozen")

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys("formulation", formulation, {"physical_equations", "gauge_constraint", "projection_formula", "projection_order", "acceleration_definition", "acceleration_order", "normal_gauge_extension_rows", "normal_gauge_derivative_order", "tilde_normal_factor", "hat_normal_factor"})
    expected_formulation = {
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_constraint": "C^a=-tilde_g^(bc)(Gamma^a_bc-bar_Gamma^a_bc)",
        "projection_formula": "H=E_tt-2sE_tr+s^2E_rr; M=E_tr-sE_rr; s=h_tr/h_rr",
        "acceleration_definition": "a=p_t_in_frozen_base_field_order",
        "normal_gauge_extension_rows": ["metric_tt_mhg", "metric_tr_mhg"],
        "normal_gauge_derivative_order": ["nabla_t_C^t", "nabla_t_C^r"],
    }
    for name, expected in expected_formulation.items():
        if formulation[name] != expected:
            raise ValueError(f"formulation.{name} is not frozen")
    _string_list("formulation.projection_order", formulation["projection_order"], CON1_PROJECTION_ORDER)
    _string_list("formulation.acceleration_order", formulation["acceleration_order"], CON1_ACCELERATION_ORDER)
    tilde = _fraction("formulation.tilde_normal_factor", formulation["tilde_normal_factor"], positive=True)
    hat = _fraction("formulation.hat_normal_factor", formulation["hat_normal_factor"], positive=True)
    if tilde != qift1.tilde_normal_factor or hat != qift1.hat_normal_factor or not 1 < tilde < hat:
        raise ValueError("CON1-COMP1 auxiliary factors differ from QIFT1")

    polynomial = _mapping("polynomial_certificate", raw["polynomial_certificate"])
    _exact_keys("polynomial_certificate", polynomial, {"variable", "maximum_degree", "evaluation_nodes", "require_unredefined_H_w_w_acceleration_independent", "require_unredefined_M_w_r_acceleration_independent", "require_unredefined_projection_rows_zero"})
    if (
        polynomial["variable"] != "a=p_t" or polynomial["maximum_degree"] != 2
        or polynomial["evaluation_nodes"] != ["zero", "plus_basis", "minus_basis", "pair_basis"]
        or any(polynomial[key] is not True for key in polynomial if key.startswith("require_"))
    ):
        raise ValueError("CON1-COMP1 polynomial certificate is not frozen")

    qift_contract = _mapping("qift1_contract", raw["qift1_contract"])
    _exact_keys("qift1_contract", qift_contract, {"parameter_dimension", "parameter_half_width", "acceleration_dimension", "acceleration_half_width", "require_datum_strictly_inside_parameter_box", "require_unique_acceleration_root"})
    if (
        qift_contract["parameter_dimension"] != 30
        or _fraction("qift1_contract.parameter_half_width", qift_contract["parameter_half_width"], positive=True) != qift1.parameter_half_width
        or qift_contract["acceleration_dimension"] != 6
        or _fraction("qift1_contract.acceleration_half_width", qift_contract["acceleration_half_width"], positive=True) != qift1.acceleration_half_width
        or qift_contract["require_datum_strictly_inside_parameter_box"] is not True
        or qift_contract["require_unique_acceleration_root"] is not True
    ):
        raise ValueError("CON1-COMP1 QIFT1 contract differs from its predecessor")
    qift_status = qift1_result.get("gate_status", {})
    if qift_status.get("quantified_full_dimensional_local_implicit_branch_box_passed") is not True:
        raise ValueError("QIFT1 unique-branch gate is not passed")

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "require_metric_defined_C_zero_exact",
        "require_all_non_normal_nabla_C_zero_exact",
        "require_fo1_reduction_rows_zero_exact",
        "require_normal_shift_zero_exact",
        "require_unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "require_ref1_scalar_equations_unmodified_exact",
        "require_exact_normal_gauge_to_extension_map_invertible",
        "require_full_mhg_root_implies_normal_nabla_C_zero",
        "require_gauge_extension_zero_at_full_mhg_root",
        "require_unredefined_full_equations_at_local_root",
        "require_no_floating_point_or_sampled_proof",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    if any(value is not True for value in proof.values()):
        raise ValueError("every CON1-COMP1 proof contract entry must be true")
    open_gates = _mapping("open_gates", raw["open_gates"])
    _exact_keys("open_gates", open_gates, set(required_con1_comp1_nonclaims()))
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every CON1-COMP1 open gate must remain false")

    return HYP1CON1COMP1Config(
        schema_version=raw["schema_version"], artifact_id=raw["artifact_id"], project_version=raw["project_version"],
        metric_signature=raw["metric_signature"], riemann_convention=raw["riemann_convention"],
        qift1_config_path=qift1_config_path, qift1_config_relative=qift1_config_relative,
        qift1_result_path=qift1_result_path, qift1_result_relative=qift1_result_relative,
        qift1_first_order_path=qift1.first_order_path,
        qift1_first_order_relative=qift1.first_order_relative,
        qift1_flat_fixture=qift1.flat_fixture,
        qift1_unique_acceleration_root_passed=qift_status.get(
            "quantified_full_dimensional_local_implicit_branch_box_passed"
        ) is True,
        predecessor_source_paths=qift1.predecessor_source_paths, predecessor_source_relatives=qift1.predecessor_source_relatives,
        radial_domain_minimum=radial_minimum, coordinate_radius=coordinate_radius, phi_value=phi_value,
        phi_radial_derivative=phi_radial_derivative,
        tilde_normal_factor=tilde, hat_normal_factor=hat, configuration=dict(raw), source_sha256=sha256(source).hexdigest(),
    )


def _serialize_exact(value: Any) -> Any:
    if isinstance(value, Interval):
        return {"lower": _fraction_text(value.lower), "upper": _fraction_text(value.upper)}
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _serialize_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serialize_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialize_exact(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise ValueError(f"certificate contains unsupported value {type(value).__name__}")


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _object_sha256(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def _reject_duplicate_json_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _canonical_json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        payload = json.loads(source, object_pairs_hook=_reject_duplicate_json_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("CON1-COMP1 result must be valid unique-key JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError("CON1-COMP1 result must be a JSON object")
    if source != _canonical_json_text(payload):
        raise ValueError("CON1-COMP1 result must use canonical sorted indented JSON encoding")
    return payload


def _theorem_composition(premises: Mapping[str, bool]) -> dict[str, Any]:
    """Close only the explicitly listed local composition implications."""

    required = {
        "qift1_unique_acceleration_root_available",
        "datum_strictly_inside_qift1_parameter_box",
        "metric_defined_C_zero_exact",
        "all_non_normal_nabla_C_components_zero_exact",
        "fo1_reduction_rows_zero_exact",
        "unredefined_H_w_w_and_M_w_r_acceleration_independent_and_zero",
        "normal_shift_zero_exact",
        "unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "ref1_scalar_equations_unmodified_exact",
        "normal_gauge_to_extension_map_invertible_exact",
        "qift_root_is_external_predecessor_not_solved_here",
    }
    if set(premises) != required or any(type(value) is not bool for value in premises.values()):
        raise ValueError("CON1-COMP1 theorem premises must be the frozen Boolean ledger")
    result = dict(premises)
    normal_root_premises = (
        "qift1_unique_acceleration_root_available",
        "datum_strictly_inside_qift1_parameter_box",
        "unredefined_H_w_w_and_M_w_r_acceleration_independent_and_zero",
        "metric_defined_C_zero_exact",
        "all_non_normal_nabla_C_components_zero_exact",
        "normal_shift_zero_exact",
        "unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "ref1_scalar_equations_unmodified_exact",
        "normal_gauge_to_extension_map_invertible_exact",
    )
    result["full_mhg_root_implies_normal_nabla_C_zero"] = all(
        premises[name] for name in normal_root_premises
    )
    result["gauge_extension_zero_at_full_mhg_root"] = all(
        premises[name] for name in normal_root_premises
    )
    result["unredefined_full_equations_at_local_root"] = all(
        premises[name] for name in normal_root_premises
    )
    return {
        "premises": dict(premises),
        "all_machine_checked_premises_pass": all(premises.values()),
        **{key: result[key] for key in result if key not in premises},
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    certificate = compatible_constraint_certificate(
        {
            "qift1_record": load_qift1_result(config.qift1_result_path),
            "qift1_configuration": {"flat_fixture": config.qift1_flat_fixture},
            "reference_radius": config.coordinate_radius,
            "activated_phi": config.phi_value,
            "activated_phi_radial_derivative": config.phi_radial_derivative,
            "parameter_half_width": _fraction(
                "qift1_contract.parameter_half_width",
                config.configuration["qift1_contract"]["parameter_half_width"],
                positive=True,
            ),
            "acceleration_half_width": _fraction(
                "qift1_contract.acceleration_half_width",
                config.configuration["qift1_contract"]["acceleration_half_width"],
                positive=True,
            ),
            "formula": config.configuration["formulation"]["projection_formula"],
            "projection_order": CON1_PROJECTION_ORDER,
            "polynomial_degree": config.configuration["polynomial_certificate"]["maximum_degree"],
            "evaluation_nodes": tuple(config.configuration["polynomial_certificate"]["evaluation_nodes"]),
        }
    )
    if not certificate.get("all_declared_exact_checks_pass"):
        raise ValueError("CON1-COMP1 exact compatibility certificate failed")
    serialized = _serialize_exact(certificate)
    predecessor_configs = tuple(config.predecessor_source_paths.values()) + (
        config.qift1_first_order_path,
        config.qift1_config_path,
        config_path.resolve(),
    )
    predecessor_results = (
        REPOSITORY / "results" / "fgc-1-action-gate.json",
        REPOSITORY / "results" / "fgc-1-metric-variation.json",
        REPOSITORY / "results" / "fgc-1-hyp1-reduction.json",
        REPOSITORY / "results" / "fgc-1-hyp1-symbol.json",
        REPOSITORY / "results" / "fgc-1-hyp1-modified-harmonic.json",
        REPOSITORY / "results" / "fgc-1-hyp1-mhg-reference.json",
        REPOSITORY / "results" / "fgc-1-hyp1-mhg-implicit.json",
        REPOSITORY / "results" / "fgc-1-hyp1-mhg-propagation.json",
        REPOSITORY / "results" / "fgc-1-hyp1-fo1-rc1.json",
        config.qift1_result_path,
    )
    if any(not path.is_file() for path in predecessor_results):
        raise ValueError("CON1-COMP1 predecessor result is absent")
    source_configs = {
        "action": config.predecessor_source_relatives["action_config"],
        "variation": config.predecessor_source_relatives["variation_config"],
        "reduction": config.predecessor_source_relatives["reduction_config"],
        "symbol": config.predecessor_source_relatives["symbol_config"],
        "modified_harmonic": config.predecessor_source_relatives["modified_harmonic_config"],
        "reference": config.predecessor_source_relatives["reference_config"],
        "implicit": config.predecessor_source_relatives["implicit_config"],
        "propagation": config.predecessor_source_relatives["propagation_config"],
        "first_order": config.qift1_first_order_relative,
        "quantified_domain": config.qift1_config_relative,
        "compatible_data": _relative(config_path),
    }
    local = certificate["local_exact_checks"]
    premises = {
        "qift1_unique_acceleration_root_available": config.qift1_unique_acceleration_root_passed,
        "datum_strictly_inside_qift1_parameter_box": local[
            "activated_parameter_point_strictly_inside_qift1_z_box"
        ],
        "metric_defined_C_zero_exact": local["metric_defined_C_zero"],
        "all_non_normal_nabla_C_components_zero_exact": local[
            "metric_defined_all_non_normal_nabla_C_components_zero"
        ],
        "fo1_reduction_rows_zero_exact": local["reduction_rows_zero"],
        "unredefined_H_w_w_and_M_w_r_acceleration_independent_and_zero": local[
            "all_56_unredefined_constraint_polynomial_coefficients_zero"
        ],
        "normal_shift_zero_exact": local["normal_shift_zero"],
        "unredefined_Ett_and_Etr_zero_for_every_acceleration": local[
            "unredefined_Ett_and_Etr_zero_for_every_acceleration"
        ],
        "ref1_scalar_equations_unmodified_exact": local[
            "ref1_scalar_equations_unmodified"
        ],
        "normal_gauge_to_extension_map_invertible_exact": local[
            "normal_gauge_extension_map_nonsingular"
        ],
        "qift_root_is_external_predecessor_not_solved_here": local[
            "qift_root_is_external_predecessor_not_solved_here"
        ],
    }
    theorem_composition = _theorem_composition(premises)
    if not theorem_composition["all_machine_checked_premises_pass"] or not all(
        value for key, value in theorem_composition.items() if key != "premises"
    ):
        raise ValueError("CON1-COMP1 theorem composition checks failed")
    return {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact activated local compatible-constraint datum gate",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "exact_local_compatible_constraint_datum_complete_propagation_initial_data_hyperbolicity_eft_and_evolution_gates_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": source_configs,
        "source_config_sha256": {_relative(path): _file_sha256(path) for path in predecessor_configs},
        "source_results_sha256": {_relative(path): _file_sha256(path) for path in predecessor_results},
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {_relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES},
        "reference": _serialize_exact(config.configuration["reference"]),
        "activated_local_parameter_datum": _serialize_exact(config.configuration["activated_local_parameter_datum"]),
        "formulation": _serialize_exact(config.configuration["formulation"]),
        "polynomial_certificate_contract": _serialize_exact(config.configuration["polynomial_certificate"]),
        "qift1_contract": _serialize_exact(config.configuration["qift1_contract"]),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "compatible_constraint_certificate": serialized,
        "compatible_constraint_certificate_sha256": _object_sha256(serialized),
        "theorem_composition": theorem_composition,
        "gate_status": {
            "exact_activated_local_compatible_constraint_datum_passed": True,
            "metric_derived_gauge_constraint_propagation_passed": False,
            "complete_reduction_constraint_system_propagation_passed": False,
            "physical_initial_constraints_solved": False,
            "uniform_radial_strong_hyperbolicity_box_passed": False,
            "retained_eft_validity_envelope_passed": False,
            "full_fgc1_hyp1_health_gate_passed": False,
            "evolution_authorized": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": serialized["nonclaims"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = record(arguments.config)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(_canonical_json_text(payload), encoding="utf-8")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
