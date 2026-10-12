#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-FO1-RC1 exact first-order-lift certificate."""

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

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-fo1-rc1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-fo1-rc1.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-fo1-rc1.md"
IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_linear_algebra.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_tangent.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_symbol.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "reference_connection.py",
    REPOSITORY
    / "src"
    / "recursive_horizons"
    / "fgc"
    / "modified_harmonic_reference.py",
    REPOSITORY
    / "src"
    / "recursive_horizons"
    / "fgc"
    / "modified_harmonic_implicit.py",
    REPOSITORY
    / "src"
    / "recursive_horizons"
    / "fgc"
    / "modified_harmonic_propagation.py",
    REPOSITORY
    / "src"
    / "recursive_horizons"
    / "fgc"
    / "modified_harmonic_first_order.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py",
    REPOSITORY / "scripts" / "reproduce_fgc_action.py",
    REPOSITORY / "scripts" / "reproduce_fgc_metric_variation.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_reduction.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_symbol.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_modified_harmonic.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_mhg_reference.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_mhg_implicit.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_mhg_propagation.py",
    Path(__file__).resolve(),
)

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config
from scripts.reproduce_fgc_hyp1_reduction import load_config as load_reduction_config
from scripts.reproduce_fgc_hyp1_symbol import load_config as load_symbol_config
from scripts.reproduce_fgc_hyp1_modified_harmonic import (
    load_config as load_modified_harmonic_config,
)
from scripts.reproduce_fgc_hyp1_mhg_reference import (
    load_config as load_reference_config,
)
from scripts.reproduce_fgc_hyp1_mhg_implicit import (
    load_config as load_implicit_config,
)
from scripts.reproduce_fgc_hyp1_mhg_propagation import (
    load_config as load_propagation_config,
)
from recursive_horizons.fgc.modified_harmonic_first_order import (
    FO1_EQUATION_ORDER,
    FO1_PHYSICAL_ARGUMENT_ORDER,
    FO1_STATE_ORDER,
    modified_harmonic_first_order_certificate,
    required_fo1_rc1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    SPHERICAL_FLAT_REFERENCE_ID,
    flat_spherical_annulus_reference,
)


Q = Fraction
FIXTURE_IDS = ("FGCQR_flat_reference_vacuum", "FGCQR_activated_generic")


@dataclass(frozen=True, slots=True)
class HYP1FO1RC1Config:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    source_paths: dict[str, Path]
    source_relatives: dict[str, str]
    radial_domain_minimum: Fraction
    coordinate_radii: dict[str, Fraction]
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    flat_fixture: dict[str, Any]
    activated_fixture: dict[str, Any]
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
        canonical = (
            str(result.numerator)
            if result.denominator == 1
            else f"{result.numerator}/{result.denominator}"
        )
        if value != canonical:
            raise ValueError(f"{name} must use canonical rational encoding")
    else:
        raise ValueError(f"{name} must be an exact rational value")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _positive_integer(name: str, value: Any) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return Q(value)


def _string_list(name: str, value: Any, expected: tuple[str, ...]) -> None:
    if not isinstance(value, list) or value != list(expected):
        raise ValueError(f"{name} is not the frozen ordered string list")


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
    if canonical != relative:
        raise ValueError(f"{name} must be canonical and traversal free")
    if not resolved.is_file():
        raise ValueError(f"{name} does not exist")
    return resolved, canonical.as_posix()


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1FO1RC1Config:
    """Load FO1-RC1 and reconcile its exact inputs through PROP1."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    path_keys = {
        "action_config",
        "variation_config",
        "reduction_config",
        "symbol_config",
        "modified_harmonic_config",
        "reference_config",
        "implicit_config",
        "propagation_config",
    }
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            *path_keys,
            "reference",
            "coordinate_radii",
            "fixtures",
            "formulation",
            "proof_contract",
            "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-FO1-RC1",
        "project_version": "0.11.0",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }
    for name, expected in expected_scalars.items():
        if raw[name] != expected:
            raise ValueError(f"{name} must equal {expected}")

    expected_paths = {
        "action_config": "configs/fgc/fgc-1-action-gate.toml",
        "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction_config": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol_config": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic_config": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference_config": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
        "implicit_config": "configs/fgc/fgc-1-hyp1-mhg-implicit.toml",
        "propagation_config": "configs/fgc/fgc-1-hyp1-mhg-propagation.toml",
    }
    source_paths: dict[str, Path] = {}
    source_relatives: dict[str, str] = {}
    for name, expected in expected_paths.items():
        resolved, relative = _repository_path(name, raw[name])
        if relative != expected:
            raise ValueError(f"{name} must name the frozen predecessor")
        source_paths[name] = resolved
        source_relatives[name] = relative

    predecessors = (
        ("ACT1", load_action_config(source_paths["action_config"])),
        ("VAR1", load_variation_config(source_paths["variation_config"])),
        ("RED1", load_reduction_config(source_paths["reduction_config"])),
        ("SYM1", load_symbol_config(source_paths["symbol_config"])),
        (
            "MHG1",
            load_modified_harmonic_config(
                source_paths["modified_harmonic_config"]
            ),
        ),
        ("REF1", load_reference_config(source_paths["reference_config"])),
        ("IMP1", load_implicit_config(source_paths["implicit_config"])),
        ("PROP1", load_propagation_config(source_paths["propagation_config"])),
    )
    for name, predecessor in predecessors:
        if predecessor.project_version != raw["project_version"]:
            raise ValueError(f"{name} and FO1-RC1 project versions differ")
        if predecessor.metric_signature != raw["metric_signature"]:
            raise ValueError(f"{name} and FO1-RC1 metric signatures differ")
    for name, predecessor in predecessors[1:]:
        if predecessor.riemann_convention != raw["riemann_convention"]:
            raise ValueError(f"{name} and FO1-RC1 Riemann conventions differ")

    reduction = predecessors[2][1]
    reference_predecessor = predecessors[5][1]
    implicit_predecessor = predecessors[6][1]
    propagation_predecessor = predecessors[7][1]

    reference = _mapping("reference", raw["reference"])
    _exact_keys(
        "reference",
        reference,
        {"reference_id", "radial_domain_minimum", "center_included"},
    )
    if reference["reference_id"] != SPHERICAL_FLAT_REFERENCE_ID:
        raise ValueError("reference.reference_id is not frozen")
    if reference["center_included"] is not False:
        raise ValueError("FO1-RC1 must not claim that its annulus includes the center")
    radial_minimum = _fraction(
        "reference.radial_domain_minimum",
        reference["radial_domain_minimum"],
        positive=True,
    )
    if not (
        radial_minimum == reference_predecessor.radial_domain_minimum
        == implicit_predecessor.radial_domain_minimum
        == propagation_predecessor.radial_domain_minimum
    ):
        raise ValueError("FO1-RC1 reference annulus differs from predecessors")

    radii = _mapping("coordinate_radii", raw["coordinate_radii"])
    _exact_keys("coordinate_radii", radii, set(FIXTURE_IDS))
    coordinate_radii = {
        fixture_id: _fraction(
            f"coordinate_radii.{fixture_id}", radii[fixture_id], positive=True
        )
        for fixture_id in FIXTURE_IDS
    }
    if any(radius <= radial_minimum for radius in coordinate_radii.values()):
        raise ValueError("every FO1-RC1 coordinate radius must lie inside the annulus")
    if (
        coordinate_radii[FIXTURE_IDS[0]] != implicit_predecessor.coordinate_radius
        or coordinate_radii != propagation_predecessor.coordinate_radii
        or coordinate_radii[FIXTURE_IDS[1]]
        != reference_predecessor.coordinate_radii[FIXTURE_IDS[1]]
    ):
        raise ValueError("FO1-RC1 coordinate radii differ from predecessors")

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "physical_equations",
            "gauge_equations",
            "first_order_kind",
            "state_definition",
            "state_order",
            "equation_order",
            "physical_argument_order",
            "radial_reduction_constraint",
            "radial_reduction_constraint_identity",
            "differentiation_method",
            "flat_branch_theorem",
            "flat_linearization",
            "tilde_normal_factor",
            "hat_normal_factor",
            "nonlinear_acceleration_solver_implemented",
        },
    )
    frozen_formulation = {
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_equations": "REF1_full_reference_connection_modified_harmonic",
        "first_order_kind": "exact_local_first_order_differential_algebraic_lift",
        "state_definition": "U=(u,p=partial_t_u,q=partial_r_u)",
        "radial_reduction_constraint": "C_r=q-partial_r_u",
        "radial_reduction_constraint_identity": "partial_t_C_r=partial_t_q-partial_r_p=0_on_declared_kinematic_rows",
        "differentiation_method": "exact_first_tangent_forward_ad",
        "flat_branch_theorem": "finite_dimensional_real_implicit_function_theorem",
        "flat_linearization": "D_A_at_z0=-J_p_t_inverse*D_z_R_REF1_at_z0",
    }
    for name, expected in frozen_formulation.items():
        if formulation[name] != expected:
            raise ValueError(f"formulation.{name} is not frozen")
    _string_list("formulation.state_order", formulation["state_order"], FO1_STATE_ORDER)
    _string_list(
        "formulation.equation_order", formulation["equation_order"], FO1_EQUATION_ORDER
    )
    _string_list(
        "formulation.physical_argument_order",
        formulation["physical_argument_order"],
        FO1_PHYSICAL_ARGUMENT_ORDER,
    )
    tilde_factor = _positive_integer(
        "formulation.tilde_normal_factor", formulation["tilde_normal_factor"]
    )
    hat_factor = _positive_integer(
        "formulation.hat_normal_factor", formulation["hat_normal_factor"]
    )
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("formulation requires 1 < tilde factor < hat factor")
    for predecessor in (
        reference_predecessor,
        implicit_predecessor,
        propagation_predecessor,
    ):
        if (
            tilde_factor != predecessor.tilde_normal_factor
            or hat_factor != predecessor.hat_normal_factor
        ):
            raise ValueError("FO1-RC1 auxiliary factors differ from predecessors")
    if formulation["nonlinear_acceleration_solver_implemented"] is not False:
        raise ValueError("FO1-RC1 cannot promote a nonlinear acceleration solver")

    fixtures = _mapping("fixtures", raw["fixtures"])
    _exact_keys(
        "fixtures",
        fixtures,
        {
            "flat_background_fixture_id",
            "activated_background_fixture_id",
            "activated_background_solution_status",
        },
    )
    if fixtures["flat_background_fixture_id"] != FIXTURE_IDS[0]:
        raise ValueError("fixtures.flat_background_fixture_id is not frozen")
    if fixtures["activated_background_fixture_id"] != FIXTURE_IDS[1]:
        raise ValueError("fixtures.activated_background_fixture_id is not frozen")
    if (
        fixtures["activated_background_solution_status"]
        != "off_shell_control_not_a_solution"
    ):
        raise ValueError("FO1-RC1 activated fixture must remain explicitly off shell")
    flat_fixture = implicit_predecessor.fixture
    activated_fixture = next(
        (
            fixture
            for fixture in reduction.configuration["fixtures"]
            if fixture["fixture_id"] == FIXTURE_IDS[1]
        ),
        None,
    )
    if activated_fixture is None:
        raise ValueError("RED1 does not contain the frozen FO1-RC1 activated fixture")

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "require_exact_eighteen_variable_and_equation_orders",
        "require_exact_first_order_two_jet_roundtrip_on_both_controls",
        "require_exact_complete_ref1_residual_reproduction_on_both_controls",
        "require_exact_kinematic_rows_and_radial_reduction_constraint_identity",
        "require_exact_all_argument_forward_tangent_primal_regression",
        "require_exact_weighted_tangent_superposition",
        "require_exact_flat_imp1_acceleration_block_regression",
        "require_exact_flat_acceleration_determinant_rank_and_inverse",
        "require_exact_complete_flat_implicit_branch_derivative",
        "require_exact_flat_first_order_principal_regression_to_mhg1",
        "require_exact_activated_off_shell_nontrivial_control",
        "require_no_nonlinear_acceleration_solve",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    if any(value is not True for value in proof.values()):
        raise ValueError("every FO1-RC1 proof contract entry must be true")

    open_gates = _mapping("open_gates", raw["open_gates"])
    _exact_keys("open_gates", open_gates, set(required_fo1_rc1_nonclaims()))
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every FO1-RC1 open gate must remain false")

    return HYP1FO1RC1Config(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        source_paths=source_paths,
        source_relatives=source_relatives,
        radial_domain_minimum=radial_minimum,
        coordinate_radii=coordinate_radii,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
        flat_fixture=flat_fixture,
        activated_fixture=activated_fixture,
        configuration=dict(raw),
        source_sha256=sha256(source).hexdigest(),
    )


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _serialize_exact(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
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
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return sha256(canonical.encode("utf-8")).hexdigest()


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=config.radial_domain_minimum
    )
    certificate = modified_harmonic_first_order_certificate(
        {
            "flat_fixture": config.flat_fixture,
            "activated_fixture": config.activated_fixture,
            "reference": reference,
            "flat_coordinate_radius": config.coordinate_radii[FIXTURE_IDS[0]],
            "activated_coordinate_radius": config.coordinate_radii[FIXTURE_IDS[1]],
            "tilde_normal_factor": config.tilde_normal_factor,
            "hat_normal_factor": config.hat_normal_factor,
        }
    )
    if not (
        certificate["all_declared_exact_checks_pass"]
        and certificate["exact_local_first_order_dae_lift_derived"]
        and certificate["complete_flat_root_linearized_first_order_map_derived"]
        and certificate["kinematic_radial_reduction_constraint_identity_derived"]
    ):
        raise ValueError("FO1-RC1 exact certificate failed a required proof check")
    serialized = _serialize_exact(certificate)
    source_paths = tuple(config.source_paths.values()) + (config_path.resolve(),)
    return {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact local first-order DAE lift, flat-root linearization, and kinematic radial reduction-constraint identity gate",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "local_first_order_dae_lift_flat_root_linearization_and_kinematic_radial_constraint_identity_complete_nonlinear_constraint_domain_and_evolution_gates_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": {
            "action": config.source_relatives["action_config"],
            "variation": config.source_relatives["variation_config"],
            "reduction": config.source_relatives["reduction_config"],
            "symbol": config.source_relatives["symbol_config"],
            "modified_harmonic": config.source_relatives[
                "modified_harmonic_config"
            ],
            "reference": config.source_relatives["reference_config"],
            "implicit": config.source_relatives["implicit_config"],
            "propagation": config.source_relatives["propagation_config"],
            "first_order": _relative(config_path),
        },
        "source_config_sha256": {
            _relative(path): _file_sha256(path) for path in source_paths
        },
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {
            _relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES
        },
        "reference": _serialize_exact(config.configuration["reference"]),
        "coordinate_radii": _serialize_exact(
            config.configuration["coordinate_radii"]
        ),
        "fixtures": _serialize_exact(config.configuration["fixtures"]),
        "formulation": _serialize_exact(config.configuration["formulation"]),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "input_fixture_ids": [
            config.flat_fixture["fixture_id"],
            config.activated_fixture["fixture_id"],
        ],
        "input_model_ids": [
            config.flat_fixture["model_id"],
            config.activated_fixture["model_id"],
        ],
        "first_order_certificate": serialized,
        "first_order_certificate_sha256": _object_sha256(serialized),
        "gate_status": {
            "exact_local_first_order_dae_lift_passed": True,
            "complete_flat_root_linearized_first_order_map_passed": True,
            "kinematic_radial_reduction_constraint_identity_passed": True,
            "explicit_nonlinear_acceleration_map_passed": False,
            "complete_metric_gauge_physical_constraint_propagation_passed": False,
            "uniform_open_domain_hyperbolicity_gate_passed": False,
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
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
