#!/usr/bin/env python3
"""Regenerate FGC-1-HYP1-BND1-MD1 from exact predecessor records."""

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
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp1-bnd1-md1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp1-bnd1-md1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp1-bnd1-md1.md"

from scripts.reproduce_fgc_hyp1_con2_mprop1 import (  # noqa: E402
    load_canonical_result as load_con2_result,
    load_config as load_con2_config,
)
from scripts.reproduce_fgc_hyp1_dom3_uhyp1 import (  # noqa: E402
    load_canonical_result as load_uhyp1_result,
    load_config as load_uhyp1_config,
)
from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config as load_fo1_config  # noqa: E402
from recursive_horizons.fgc.exact_interval import Interval, interval  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_radial_boundary import (  # noqa: E402
    frozen_radial_modal_boundary_certificate,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HYP1-BND1-MD1"
PROJECT_VERSION = "0.11.0"
NONCLAIM_NAMES = (
    "complete_nonlinear_lower_order_first_order_sources_derived",
    "quasilinear_local_existence_or_IBVP_proven",
    "metric_derived_gauge_constraint_boundary_map_proven",
    "physical_Hamiltonian_momentum_constraint_boundary_map_proven",
    "constraint_preserving_ACT1_IBVP_proven",
    "initial_boundary_corner_compatibility_constructed",
    "regular_center_or_asymptotic_boundary_treatment_proven",
    "multidirectional_boundary_stability_proven",
    "retained_EFT_domain_proven",
    "evolution_authorized",
    "collapse_solution_derived",
    "metric_null_affine_defocusing_derived",
    "finite_invariant_transition_surface_derived",
    "singularity_resolution_derived",
)
NONCLAIMS = {name: False for name in NONCLAIM_NAMES}
PREDECESSOR_PATHS = {
    "uniform_domain_config": "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml",
    "uniform_domain_result": "results/fgc-1-hyp1-dom3-uhyp1.json",
    "metric_subsidiary_config": "configs/fgc/fgc-1-hyp1-con2-mprop1.toml",
    "metric_subsidiary_result": "results/fgc-1-hyp1-con2-mprop1.json",
    "first_order_config": "configs/fgc/fgc-1-hyp1-fo1-rc1.toml",
    "first_order_result": "results/fgc-1-hyp1-fo1-rc1.json",
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/exact_interval.py",
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/modified_harmonic_radial_boundary.py",
        "src/recursive_horizons/fgc/modified_harmonic_reduction_subsidiary.py",
        "scripts/reproduce_fgc_hyp1_bnd1_md1.py",
    )
)


@dataclass(frozen=True, slots=True)
class Config:
    raw: dict[str, Any]
    paths: dict[str, Path]
    source_sha256: str


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(
            f"{name} keys differ: missing={sorted(expected - set(value))}, "
            f"extra={sorted(set(value) - expected)}"
        )


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if not isinstance(value, str) or value != expected:
        raise ValueError(f"{name} must name {expected}")
    candidate = (REPOSITORY / Path(value)).resolve()
    try:
        relative = candidate.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if relative != value or not candidate.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return candidate


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if value != _fraction_text(answer):
        raise ValueError(f"{name} must be a canonical rational string")
    return answer


def _interval(name: str, value: Any) -> Interval:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an exact interval record")
    _keys(name, value, {"lower", "upper"})
    lower = _fraction(f"{name}.lower", value["lower"])
    upper = _fraction(f"{name}.upper", value["upper"])
    if lower > upper:
        raise ValueError(f"{name} has reversed endpoints")
    return interval(lower, upper)


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported exact certificate type {type(value).__name__}")


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _canonical_json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _load_canonical_predecessor_result(path: Path, name: str) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} predecessor must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical_json_text(value):
        raise ValueError(f"{name} predecessor must use canonical sorted JSON")
    return value


def _verify_direct_config_binding(
    name: str,
    result: Mapping[str, Any],
    config_path: Path,
) -> None:
    ledger = result.get("source_config_sha256")
    relative = _rel(config_path)
    if not isinstance(ledger, Mapping) or ledger.get(relative) != _sha(config_path):
        raise ValueError(f"{name} result is stale against its direct source config")


def _validate_result_header(
    name: str,
    result: Mapping[str, Any],
    artifact_id: str,
    passed_gate: str,
) -> None:
    if result.get("artifact_id") != artifact_id or result.get("project_version") != PROJECT_VERSION:
        raise ValueError(f"{name} result identity or version differs")
    gates = result.get("gate_status")
    if not isinstance(gates, Mapping) or gates.get(passed_gate) is not True:
        raise ValueError(f"{name} required predecessor gate is not passed")
    nonclaims = result.get("nonclaims")
    if not isinstance(nonclaims, Mapping) or any(value is not False for value in nonclaims.values()):
        raise ValueError(f"{name} predecessor nonclaim was promoted")


def load_config(path: Path = DEFAULT_CONFIG) -> Config:
    source = path.resolve().read_bytes()
    raw = tomllib.loads(source.decode("utf-8"))
    _keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            *PREDECESSOR_PATHS,
            "boundary",
            "proof_contract",
            "open_gates",
        },
    )
    identity = {
        key: raw[key]
        for key in (
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
        )
    }
    if identity != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("BND1-MD1 identity, version, or conventions differ")
    paths = {
        name: _path(name, raw[name], relative)
        for name, relative in PREDECESSOR_PATHS.items()
    }
    # Parse each predecessor configuration through its owning strict parser.
    load_uhyp1_config(paths["uniform_domain_config"])
    load_con2_config(paths["metric_subsidiary_config"])
    load_fo1_config(paths["first_order_config"])
    _keys(
        "boundary",
        raw["boundary"],
        {
            "spatial_domain",
            "inner_outward_normal",
            "outer_outward_normal",
            "operator",
            "expected_inner_incoming_propagating_modes",
            "expected_outer_incoming_propagating_modes",
            "expected_zero_speed_u_variables",
            "conditions_on_zero_speed_u_variables",
        },
    )
    if raw["boundary"] != {
        "spatial_domain": "finite_spherical_annulus",
        "inner_outward_normal": "-partial_r",
        "outer_outward_normal": "+partial_r",
        "operator": "homogeneous_incoming_modal_absorption",
        "expected_inner_incoming_propagating_modes": 6,
        "expected_outer_incoming_propagating_modes": 6,
        "expected_zero_speed_u_variables": 6,
        "conditions_on_zero_speed_u_variables": 0,
    }:
        raise ValueError("BND1-MD1 boundary orientation/operator contract differs")
    _keys(
        "proof_contract",
        raw["proof_contract"],
        {
            "same_compact_REF1_branch_graph_as_UHYP1",
            "strict_interval_speed_sign_classification",
            "exact_incoming_selector_rank",
            "exact_modal_energy_flux_identity",
            "normalized_modal_Lopatinski_identity",
            "uniform_propagating_symmetrizer_coercivity",
            "complete_kinematic_reduction_subsidiary_consumed",
            "zero_incoming_kinematic_reduction_constraint_fields",
            "frozen_1plus1_principal_main_system_only",
            "exact_rational_arithmetic_only",
        },
    )
    if any(value is not True for value in raw["proof_contract"].values()):
        raise ValueError("every BND1-MD1 proof-contract gate must remain true")
    if raw["open_gates"] != NONCLAIMS:
        raise ValueError("BND1-MD1 open gate contract differs")
    return Config(raw=dict(raw), paths=paths, source_sha256=sha256(source).hexdigest())


def _rehydrate_uhyp1(result: Mapping[str, Any]) -> dict[str, Any]:
    certificate = result.get("uniform_radial_certificate")
    if not isinstance(certificate, Mapping):
        raise ValueError("UHYP1 result lacks its uniform radial certificate")
    modes = certificate.get("modes")
    if not isinstance(modes, Mapping):
        raise ValueError("UHYP1 result lacks ordered modal records")
    live_modes: dict[str, list[dict[str, Any]]] = {}
    for sector in ("tilde", "hat", "physical_chi", "regulator"):
        family = modes.get(sector)
        if not isinstance(family, list):
            raise ValueError(f"UHYP1 serialized {sector} family is not ordered")
        live_modes[sector] = [
            {
                "speed_box": _interval(
                    f"uniform_radial_certificate.modes.{sector}[{index}].speed_box",
                    mode.get("speed_box") if isinstance(mode, Mapping) else None,
                )
            }
            for index, mode in enumerate(family)
        ]
    domain = certificate.get("domain")
    frame = certificate.get("eigenframe")
    symmetrizer = certificate.get("radial_symmetrizer")
    nonclaims = certificate.get("nonclaims")
    if not all(isinstance(item, Mapping) for item in (domain, frame, symmetrizer, nonclaims)):
        raise ValueError("UHYP1 serialized domain/frame/symmetrizer/scope is incomplete")
    neumann = frame.get("neumann_inverse")
    if not isinstance(neumann, Mapping):
        raise ValueError("UHYP1 serialized Neumann inverse is absent")
    return {
        "classification": certificate.get("classification"),
        "domain": dict(domain),
        "modes": live_modes,
        "eigenframe": {
            "all_enclosed_frames_invertible": frame.get("all_enclosed_frames_invertible"),
            "real_smooth_radial_eigenframe": frame.get("real_smooth_radial_eigenframe"),
            "neumann_inverse": {
                "dimension": neumann.get("dimension"),
                "rho_infinity": _fraction(
                    "uniform_radial_certificate.eigenframe.neumann_inverse.rho_infinity",
                    neumann.get("rho_infinity"),
                ),
            },
        },
        "radial_symmetrizer": {
            "definition": symmetrizer.get("definition"),
            "HA_symmetric": symmetrizer.get("HA_symmetric"),
            "euclidean_coercivity_lower_bound": _fraction(
                "uniform_radial_certificate.radial_symmetrizer.euclidean_coercivity_lower_bound",
                symmetrizer.get("euclidean_coercivity_lower_bound"),
            ),
            "euclidean_coercivity_upper_bound": _fraction(
                "uniform_radial_certificate.radial_symmetrizer.euclidean_coercivity_upper_bound",
                symmetrizer.get("euclidean_coercivity_upper_bound"),
            ),
        },
        "nonclaims": dict(nonclaims),
    }


def _kinematic_reduction_compatibility(result: Mapping[str, Any]) -> dict[str, Any]:
    controls = result.get("controls")
    if not isinstance(controls, Mapping) or not controls:
        raise ValueError("CON2-MPROP1 controls are absent")
    summaries: dict[str, Any] = {}
    field_order: list[str] | None = None
    for name, control in controls.items():
        reduction = control.get("reduction_subsidiary") if isinstance(control, Mapping) else None
        if not isinstance(reduction, Mapping):
            raise ValueError(f"CON2-MPROP1 {name} reduction subsidiary is absent")
        order = reduction.get("field_order")
        speeds = reduction.get("subsidiary_characteristic_speeds")
        speeds_are_exact_zero = (
            isinstance(speeds, list)
            and len(speeds) == 6
            and all(
                (value == 0 if isinstance(value, Fraction) else _fraction("reduction speed", value) == 0)
                for value in speeds
            )
        )
        if not isinstance(order, list) or len(order) != 6 or not speeds_are_exact_zero:
            raise ValueError(f"CON2-MPROP1 {name} reduction characteristic contract differs")
        if field_order is None:
            field_order = order
        elif field_order != order:
            raise ValueError("CON2-MPROP1 controls use different field orders")
        required = (
            reduction.get("complete_kinematic_reduction_subsidiary_system_derived") is True
            and reduction.get("off_shell_identity_exact") is True
            and reduction.get("on_shell_radial_constraint_is_time_constant") is True
            and reduction.get("incoming_reduction_constraint_fields_at_inner_boundary") == []
            and reduction.get("incoming_reduction_constraint_fields_at_outer_boundary") == []
        )
        if not required:
            raise ValueError(f"CON2-MPROP1 {name} reduction compatibility gate failed")
        summaries[str(name)] = {
            "off_shell_identity": "partial_t_C+partial_r_D-K=0",
            "on_D_K_shell": "partial_t_C=0",
            "characteristic_speeds": tuple(Q(0) for _ in range(6)),
            "incoming_fields_inner": [],
            "incoming_fields_outer": [],
        }
    return {
        "classification": "complete_1plus1_kinematic_reduction_zero_speed_compatibility_not_full_ACT1_constraints",
        "field_order": field_order,
        "controls": summaries,
        "independent_reduction_constraint_boundary_conditions_required": 0,
        "compatible_with_no_boundary_conditions_on_zero_speed_u_variables": True,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    uhyp1 = load_uhyp1_result(config.paths["uniform_domain_result"])
    con2 = load_con2_result(config.paths["metric_subsidiary_result"])
    fo1 = _load_canonical_predecessor_result(
        config.paths["first_order_result"], "FO1-RC1"
    )
    _validate_result_header(
        "UHYP1",
        uhyp1,
        "FGC-1-HYP1-DOM3-UHYP1",
        "compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",
    )
    _validate_result_header(
        "CON2-MPROP1",
        con2,
        "FGC-1-HYP1-CON2-MPROP1",
        "complete_kinematic_1plus1_reduction_subsidiary",
    )
    _validate_result_header(
        "FO1-RC1",
        fo1,
        "FGC-1-HYP1-FO1-RC1",
        "exact_local_first_order_dae_lift_passed",
    )
    _verify_direct_config_binding("UHYP1", uhyp1, config.paths["uniform_domain_config"])
    _verify_direct_config_binding("CON2-MPROP1", con2, config.paths["metric_subsidiary_config"])
    _verify_direct_config_binding("FO1-RC1", fo1, config.paths["first_order_config"])
    boundary = frozen_radial_modal_boundary_certificate(_rehydrate_uhyp1(uhyp1))
    if boundary.get("nonclaims") != NONCLAIMS:
        raise ValueError("boundary core and artifact nonclaim contracts differ")
    kinematic = _kinematic_reduction_compatibility(con2)
    characteristics = boundary["characteristics"]
    operators = boundary["boundary_operators"]
    declared = config.raw["boundary"]
    if (
        len(characteristics["inner_incoming_indices"])
        != declared["expected_inner_incoming_propagating_modes"]
        or len(characteristics["outer_incoming_indices"])
        != declared["expected_outer_incoming_propagating_modes"]
        or boundary["state_order"]["zero_speed_dimension"]
        != declared["expected_zero_speed_u_variables"]
        or operators["conditions_imposed_on_zero_speed_u_variables"]
        != declared["conditions_on_zero_speed_u_variables"]
    ):
        raise ValueError("derived boundary dimensions differ from the frozen contract")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": "exact_uniform_frozen_1plus1_radial_main_system_maximal_dissipation_plus_kinematic_reduction_compatibility_not_constraint_preserving_ACT1_IBVP",
        "project_version": PROJECT_VERSION,
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): config.source_sha256},
        "predecessor_sha256": {
            _rel(path): _sha(path) for path in config.paths.values()
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "frozen_configuration": _serial(config.raw),
        "uniform_frozen_boundary_certificate": _serial(boundary),
        "kinematic_reduction_compatibility": _serial(kinematic),
        "gate_status": {
            "uniform_frozen_radial_main_system_boundary_dissipation_passed": True,
            "kinematic_reduction_zero_speed_compatibility_passed": True,
            "constraint_preserving_ACT1_IBVP_proven": False,
            "evolution_authorized": False,
        },
        "nonclaims": NONCLAIMS,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("BND1-MD1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical_json_text(value):
        raise ValueError("BND1-MD1 result must use canonical sorted JSON")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(
        _canonical_json_text(record(arguments.config)), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
