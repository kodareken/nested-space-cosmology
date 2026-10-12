"""Prospective compact freeze for the candidate-blind DEF1-STAB1 instrument.

Ordinary ``verify_compact`` reads only the tracked config and result.  It does
not reconstruct from live owners, inspect Git, open ``runs/``, read a
trajectory, or execute a runner.  ``compose_canonical_artifacts`` is explicit
initial construction from the current DEF1 owner modules.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    EvidenceIOError,
    UnsafePathError,
    canonical_json_bytes,
    read_regular_file,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-DEF1-STAB1-FRZ1"
PARTIAL_INSTRUMENT_ID = "FGC-1-DEF1-STAB1"
PREF1_ARTIFACT_ID = "FGC-1-DEF1-STAB1-PREF1"
CLASSIFICATION = (
    "candidate_blind_conversion_error_map_instrument_freeze_no_trajectory"
)
PROJECT_VERSION = "0.11.0"
BASE_COMMIT = "c4feb941e408a7c72913b49071801b7531183c1a"
CONFIG_PATH = "configs/fgc/fgc-1-def1-stab1-frz1.toml"
RESULT_PATH = "results/fgc-1-def1-stab1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-def1-stab1-frz1.md"
CONFIG_SHA256 = "c88676a78cc64d94f9dc7afadb557c1bdf4af8cc09e2fde259d4f3210bc5059f"
RESULT_SHA256 = "c2e5c8c2177d7f0bcf414f448131c0abadcab28a5b5c6b2351dd216dddb5c126"
_SHA_LENGTH = 64
_HEX = frozenset("0123456789abcdef")

TRUE_CLAIMS = (
    "conversion_instrument_contract_complete",
    "directly_exercised_conversion_identities",
    "directly_exercised_provider_route_coverage",
    "directly_exercised_q_error_assembly",
    "directly_exercised_misner_sharp_enclosure",
    "directly_exercised_activation_controls",
    "directly_exercised_trappedness_and_complete_q_margins",
    "directly_exercised_base_to_adm_roundtrip",
    "directly_exercised_imp1_to_q_debit",
    "licenses_separate_independent_pref1_binder",
)
FALSE_CLAIMS = (
    "def1_error_map_passed",
    "def1_booleans_evaluated",
    "candidate_or_control_trajectory_read",
    "global_pde_error_certified",
    "used_measured_q",
    "richardson_treated_as_global_pde_error",
    "imp1_admission_debit_treated_as_global_pde_error",
    "rob1_passed",
    "holdout_authorized",
    "physical_claimed",
    "mechanism_claimed",
    "pref1_implemented",
    "universal_pde_theorem",
    "fitted_threshold",
    "measured_q_dependent_error_shrinkage",
)
NONCLAIMS = (
    "FRZ1 proves only that the candidate-blind conversion/error-map "
    "instrument contract is complete and frozen.",
    "FRZ1 contains no candidate or control trajectory values.",
    "FRZ1 does not establish a global nonlinear PDE error theorem.",
    "Richardson and IMP1 debits remain conditional premises.",
    "FRZ1 does not evaluate the nine DEF1 booleans or ROB1.",
    "FRZ1 does not set DEF1_error_map_passed or authorize holdout.",
    "FRZ1 licenses only a separately implemented independent "
    "FGC-1-DEF1-STAB1-PREF1 binder.",
)


class Def1Stab1Frz1Error(ValueError):
    """The compact DEF1-STAB1-FRZ1 freeze could not be established."""


def _fail(message: str) -> None:
    raise Def1Stab1Frz1Error(message)


def _sha(raw: bytes) -> str:
    if type(raw) is not bytes:
        _fail("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_digest(value: object, label: str) -> str:
    if (
        type(value) is not str
        or len(value) != _SHA_LENGTH
        or any(character not in _HEX for character in value)
    ):
        _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _pinned_digest(value: str, label: str) -> str:
    if value.startswith("PLACEHOLDER_"):
        _fail(f"{label} has not been bound")
    return _require_digest(value, label)


def _mapping(value: object, label: str) -> dict[str, Any]:
    if type(value) is not dict:
        _fail(f"{label} must be an object")
    return value


def _boolean(value: object, expected: bool, label: str) -> bool:
    if type(value) is not bool or value is not expected:
        _fail(f"{label} differs")
    return value


def _same(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        _fail(f"{label} differs")


def _root(repository: Path) -> Path:
    path = Path(repository)
    if not path.is_absolute():
        _fail("repository root is not canonical")
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise Def1Stab1Frz1Error("repository root cannot be resolved") from error
    if resolved != path:
        _fail("repository root traverses a symlink")
    return path


def _read_tracked(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise Def1Stab1Frz1Error(f"tracked FRZ1 leaf is absent: {relative}") from error


def expected_claims() -> dict[str, bool]:
    claims = {name: True for name in TRUE_CLAIMS}
    claims.update({name: False for name in FALSE_CLAIMS})
    return claims


def _toml_scalar(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if type(value) is str:
        return json.dumps(value, ensure_ascii=True)
    if type(value) is list:
        return json.dumps(value, ensure_ascii=True, separators=(",", ":"))
    _fail("unsupported FRZ1 config scalar")
    raise AssertionError("unreachable")


def render_config(config: Mapping[str, object]) -> bytes:
    mapping = dict(config)
    tables = ("owner_hashes", "conversion_identity_hashes", "claims")
    lines = []
    for key, value in mapping.items():
        if key in tables:
            continue
        lines.append(f"{key} = {_toml_scalar(value)}")
    for table in tables:
        lines.append(f"\n[{table}]")
        nested = _mapping(mapping.get(table), table)
        for key, value in nested.items():
            rendered_key = json.dumps(key, ensure_ascii=True)
            lines.append(f"{rendered_key} = {_toml_scalar(value)}")
    return ("\n".join(lines) + "\n").encode("ascii")


def canonical_result(value: object) -> bytes:
    return canonical_json_bytes(value) + b"\n"


def _unevaluated_booleans(names: tuple[str, ...]) -> dict[str, None]:
    return {name: None for name in names}


def compose_canonical_artifacts(repository: Path | None = None) -> tuple[bytes, bytes]:
    """Reconstruct the compact freeze from live DEF1 owner modules.

    This path is initial construction only.  Ordinary ``--check`` must not
    call it.
    """

    del repository
    from fractions import Fraction as Q

    from recursive_horizons.fgc.def1_geometry_error import INPUT_NAMES
    from recursive_horizons.fgc.def1_stab1 import (
        DEF1_BOOLEAN_NAMES,
        ERROR_BUDGET_COMPONENTS,
        GR0_PHI_POLICY_CANONICAL_ZERO,
        PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
        QComponentPremise,
        SOURCE_SUPPLIED_CERTIFIED,
        STOPPED_MASS_FLUX_INCONSISTENCY,
        PREMISE_STATUS_PROVEN,
        ActivationAssessment,
        MatchedControlSample,
        assemble_q_error_budget,
        assess_activation,
        assess_mass_flux_ledger,
        complete_q_margin_passed,
        inject_component_input_error,
        inverse_base_metric_from_adm,
        misner_sharp_mass,
        trappedness_margin_passed,
        zero_error_premises,
    )
    from recursive_horizons.fgc.def1_stab1_freeze_contract import (
        OWNER_PATHS,
        build_def1_stab1_freeze_contract,
    )
    from recursive_horizons.fgc.def1_stab1_qualification import (
        CONVERSION_BASE_TO_ADM,
        CONVERSION_COVERAGE_MATRIX,
        CONVERSION_FO1_ADM_GEOMETRY,
        CONVERSION_IMP1_TO_Q,
        CONVERSION_RED1_MHG2,
        IMP1_18_CHANNEL_ORDER,
        PROVIDER_ROUTES,
        convert_base_to_adm_geometry_slots,
        convert_imp1_channels_to_q,
        conversion_identity_digest,
        execute_qualification_coverage_matrix,
        fo1_adm_geometry_slot_map,
        invert_base_metric_two_jets,
        red1_mhg2_equation_identity_map,
    )
    from recursive_horizons.fgc.exact_interval import Interval
    from recursive_horizons.fgc.spherical_reduction import Jet2

    contract = build_def1_stab1_freeze_contract()
    if not contract.conversion_instrument_contract_complete:
        _fail("freeze contract is incomplete")
    for name in FALSE_CLAIMS:
        if name in contract.payload and contract.payload[name] is not False:
            _fail(f"freeze contract promoted {name}")

    equation_map = red1_mhg2_equation_identity_map()
    geometry_map = fo1_adm_geometry_slot_map()
    if equation_map.claims_col1_values or geometry_map.claims_col1_values:
        _fail("conversion identities claimed COL1 values")
    conversion_hashes = {
        CONVERSION_RED1_MHG2: conversion_identity_digest(CONVERSION_RED1_MHG2),
        CONVERSION_FO1_ADM_GEOMETRY: conversion_identity_digest(
            CONVERSION_FO1_ADM_GEOMETRY
        ),
        CONVERSION_BASE_TO_ADM: conversion_identity_digest(CONVERSION_BASE_TO_ADM),
        CONVERSION_IMP1_TO_Q: conversion_identity_digest(CONVERSION_IMP1_TO_Q),
        CONVERSION_COVERAGE_MATRIX: conversion_identity_digest(
            CONVERSION_COVERAGE_MATRIX
        ),
    }
    _same(
        conversion_hashes,
        dict(contract.payload["conversion_identity_hashes"]),
        "conversion identity hashes",
    )

    coverage = execute_qualification_coverage_matrix()
    if len(coverage.routes) != len(PROVIDER_ROUTES):
        _fail("provider-route coverage is incomplete")
    if any(
        not item.positive_control_passed or not item.injected_failure_refused
        for item in coverage.routes
    ):
        _fail("a provider route lacked positive or injected-failure coverage")
    if coverage.trajectory_values_evaluated or coverage.def1_booleans_evaluated:
        _fail("coverage evaluated trajectory values or DEF1 booleans")
    if coverage.def1_error_map_passed or coverage.used_measured_q:
        _fail("coverage promoted the error-map gate or used measured Q")

    zero_assembled = assemble_q_error_budget(zero_error_premises())
    if zero_assembled.exact_total != 0 or zero_assembled.used_measured_q:
        _fail("zero Q-error control used a measured signal")
    try:
        assemble_q_error_budget(zero_error_premises(), measured_q=Q(1, 1000))
    except Exception as error:
        measured_q_refused = "measured Q" in str(error)
    else:
        measured_q_refused = False
    if not measured_q_refused:
        _fail("Q-error assembly accepted a measured-Q argument")

    interpolation = QComponentPremise(
        component="interpolation",
        sensitivity=Interval(-Q(2), Q(5, 2)),
        input_error=Interval(Q(0), Q(1, 3)),
        status=PREMISE_STATUS_PROVEN,
        source=SOURCE_SUPPLIED_CERTIFIED,
    )
    spatial = QComponentPremise(
        component="spatial_temporal",
        sensitivity=Interval.singleton(2),
        input_error=Interval.singleton(Q(1, 5)),
        status=PREMISE_STATUS_PROVEN,
        source=SOURCE_SUPPLIED_CERTIFIED,
    )
    mixed_premises = []
    for premise in zero_error_premises():
        if premise.component == "interpolation":
            mixed_premises.append(interpolation)
        elif premise.component == "spatial_temporal":
            mixed_premises.append(spatial)
        else:
            mixed_premises.append(premise)
    mixed_assembled = assemble_q_error_budget(tuple(mixed_premises))
    expected_mixed = Q(5, 2) * Q(1, 3) + Q(2) * Q(1, 5)
    if mixed_assembled.exact_total != expected_mixed:
        _fail("componentwise Q-error product sum differs")
    if mixed_assembled.used_measured_q or mixed_assembled.global_pde_error_certified:
        _fail("mixed Q-error assembly promoted a PDE or measured-Q claim")
    injected = inject_component_input_error(tuple(mixed_premises), "interpolation", Q(1, 7))
    injected_assembled = assemble_q_error_budget(injected)
    if injected_assembled.exact_total != expected_mixed + Q(5, 2) * Q(1, 7):
        _fail("injected Q-error increment depended on a measured signal")

    inverse = inverse_base_metric_from_adm(1, 0, 1)
    gradient = (Q(0), Q(1))
    mixed = ((0, 0), (0, 0))
    minkowski_mass = misner_sharp_mass(4, inverse, gradient)
    if minkowski_mass != 0:
        _fail("Minkowski Misner-Sharp mass control failed")
    valid_ledger = assess_mass_flux_ledger(
        radius=4,
        inverse_metric=inverse,
        radius_derivatives=gradient,
        coupling_F=1,
        mixed_equation_residual=mixed,
        mixed_einstein=mixed,
        delta_mass=Q(1, 2),
        integrated_flux=Q(1, 2),
        residual_enclosure=0,
    )
    if not valid_ledger.ledger_valid:
        _fail("exact conservation enclosure control failed")
    vetoed = assess_mass_flux_ledger(
        radius=4,
        inverse_metric=inverse,
        radius_derivatives=gradient,
        coupling_F=1,
        mixed_equation_residual=mixed,
        delta_mass=Q(1, 2),
        integrated_flux=0,
        residual_enclosure=Q(1, 4),
    )
    if (
        vetoed.ledger_valid
        or vetoed.protocol_token != PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN
        or vetoed.veto_class != STOPPED_MASS_FLUX_INCONSISTENCY
    ):
        _fail("typed mass-flux veto control failed")

    activation = assess_activation(
        max_abs_phi=Q(5, 2),
        activation_error=Q(1, 2),
        s_ref=Q(1),
        gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
        matched_controls=(
            MatchedControlSample("GR-0", Q(0), Q(0)),
            MatchedControlSample("SGB-L", Q(1, 4), Q(0)),
        ),
    )
    if not isinstance(activation, ActivationAssessment):
        _fail("activation control did not return the frozen DTO")
    if not activation.threshold_passed or not activation.control_dominance_passed:
        _fail("activation/control comparison control failed")
    dominated = assess_activation(
        max_abs_phi=Q(5, 2),
        activation_error=Q(1, 2),
        s_ref=Q(1),
        gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
        matched_controls=(
            MatchedControlSample("GR-0", Q(0), Q(0)),
            MatchedControlSample("SGB-L", Q(2), Q(0)),
        ),
    )
    if dominated.control_dominance_passed:
        _fail("overlapping control upper bound was treated as dominance")

    error = Q(1, 8)
    trapped_pass = trappedness_margin_passed(
        -4 * error - Q(1, 2**40),
        -4 * error - Q(1, 2**40),
        error,
        error,
    )
    trapped_fail = trappedness_margin_passed(-4 * error, -4 * error, error, error)
    q_error = Q(1, 16)
    complete_pass = complete_q_margin_passed(4 * q_error + Q(1, 2**40), q_error)
    complete_fail = complete_q_margin_passed(4 * q_error, q_error)
    if not trapped_pass or trapped_fail or not complete_pass or complete_fail:
        _fail("trappedness or complete-Q margin helper control failed")

    alpha = Jet2(Q(2), Q(1, 3), -Q(1, 5), Q(1, 7), -Q(1, 11), Q(1, 13))
    shift = Jet2(Q(1, 4), -Q(1, 6), Q(1, 8), -Q(1, 10), Q(1, 12), -Q(1, 14))
    radial = Jet2(Q(3, 2), Q(1, 9), -Q(1, 15), Q(1, 17), Q(1, 19), -Q(1, 21))
    radius = Jet2(Q(4), -Q(1, 7), Q(5, 4), Q(1, 23), -Q(1, 25), Q(1, 27))
    radial_squared = radial * radial
    base_inverse = invert_base_metric_two_jets(
        h_tt=-(alpha * alpha) + radial_squared * shift * shift,
        h_tr=radial_squared * shift,
        h_rr=radial_squared,
        areal_radius=radius,
        lapse_root=alpha.value,
        radial_scale_root=radial.value,
    )
    if not base_inverse.roundtrip_holds:
        _fail("BASE-to-ADM roundtrip control failed")
    packed = convert_base_to_adm_geometry_slots(
        base_inverse,
        tangent=(("k.t", Q(5, 2)), ("k.r", -Q(1, 4))),
    )
    if not packed.complete or tuple(name for name, _ in packed.supplied) != INPUT_NAMES:
        _fail("BASE-to-ADM control did not produce all 26 geometry inputs")

    channel_debits = tuple(
        (name, Q(index + 1, 1024)) for index, name in enumerate(IMP1_18_CHANNEL_ORDER)
    )
    lipschitz = tuple((name, Q(1, 32)) for name in IMP1_18_CHANNEL_ORDER)
    imp1 = convert_imp1_channels_to_q(
        channel_debits=channel_debits,
        lipschitz=lipschitz,
        context="DEF1-STAB1-FRZ1 independent IMP1-to-Q control",
        unit="declared complete-Q control units",
    )
    if (
        imp1.used_measured_q
        or imp1.global_pde_error_certified
        or imp1.def1_error_map_passed
        or imp1.imp1_admission_debit_treated_as_global_pde_error
    ):
        _fail("IMP1-to-Q control promoted a PDE or error-map claim")
    _same(
        {
            "numerator": imp1.additive_q.numerator,
            "denominator": imp1.additive_q.denominator,
        },
        dict(contract.payload["imp1_control_additive_q"]),
        "IMP1 additive Q",
    )

    claims = expected_claims()
    owner_hashes = dict(contract.payload["owner_hashes"])
    if set(owner_hashes) != set(OWNER_PATHS):
        _fail("owner-hash inventory differs")
    config = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "base_commit": BASE_COMMIT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": PARTIAL_INSTRUMENT_ID,
        "pref1_artifact_id": PREF1_ARTIFACT_ID,
        "freeze_contract_payload_sha256": contract.payload_sha256,
        "conversion_instrument_contract_complete": True,
        "def1_error_map_passed": False,
        "def1_booleans_evaluated": False,
        "candidate_or_control_trajectory_read": False,
        "global_pde_error_certified": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
        "pref1_implemented": False,
        "def1_boolean_names": list(DEF1_BOOLEAN_NAMES),
        "error_budget_component_order": list(ERROR_BUDGET_COMPONENTS),
        "owner_hashes": owner_hashes,
        "conversion_identity_hashes": conversion_hashes,
        "claims": claims,
    }
    config_raw = render_config(config)
    result = {
        "artifact_id": ARTIFACT_ID,
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "base_commit": BASE_COMMIT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": PARTIAL_INSTRUMENT_ID,
        "pref1_artifact_id": PREF1_ARTIFACT_ID,
        "config_sha256": _sha(config_raw),
        "freeze_contract_payload_sha256": contract.payload_sha256,
        "controls": {
            "conversion_identities": {
                "red1_mhg2_pair_count": len(equation_map.pairs),
                "red1_mhg2_conversion_identity": equation_map.conversion_identity,
                "fo1_adm_conversion_identity": geometry_map.conversion_identity,
                "hashes": conversion_hashes,
            },
            "provider_route_coverage": {
                "route_count": len(coverage.routes),
                "component_count": len(coverage.components),
                "positive_controls_passed": all(
                    item.positive_control_passed for item in coverage.routes
                ),
                "injected_failures_refused": all(
                    item.injected_failure_refused for item in coverage.routes
                ),
                "routes": [
                    [item.component, item.route_id, item.conversion]
                    for item in coverage.routes
                ],
            },
            "q_error_assembly": {
                "zero_total": {
                    "numerator": 0,
                    "denominator": 1,
                },
                "mixed_total": {
                    "numerator": mixed_assembled.exact_total.numerator,
                    "denominator": mixed_assembled.exact_total.denominator,
                },
                "injected_total": {
                    "numerator": injected_assembled.exact_total.numerator,
                    "denominator": injected_assembled.exact_total.denominator,
                },
                "measured_q_argument_refused": True,
                "used_measured_q": False,
                "global_pde_error_certified": False,
            },
            "misner_sharp": {
                "minkowski_mass": {
                    "numerator": 0,
                    "denominator": 1,
                },
                "ledger_valid_control_passed": True,
                "typed_veto_control_passed": True,
                "protocol_token": PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
                "veto_class": STOPPED_MASS_FLUX_INCONSISTENCY,
            },
            "activation": {
                "threshold_passed": True,
                "control_dominance_passed": True,
                "overlapping_control_fails_dominance": True,
                "canonical_gr0_phi_zero": True,
            },
            "margins": {
                "trappedness_strict_pass": True,
                "trappedness_boundary_fails": True,
                "complete_q_strict_pass": True,
                "complete_q_boundary_fails": True,
            },
            "base_to_adm": {
                "roundtrip_holds": True,
                "geometry_inputs_complete": True,
                "geometry_input_count": len(INPUT_NAMES),
            },
            "imp1_to_q": {
                "channel_count": len(IMP1_18_CHANNEL_ORDER),
                "additive_q": {
                    "numerator": imp1.additive_q.numerator,
                    "denominator": imp1.additive_q.denominator,
                },
                "conditional_premise": True,
                "def1_error_map_passed": False,
                "used_measured_q": False,
            },
        },
        "claims": claims,
        "nonclaims": list(NONCLAIMS),
        "def1_boolean_names": list(DEF1_BOOLEAN_NAMES),
        "def1_booleans_evaluated": False,
        "def1_booleans": _unevaluated_booleans(DEF1_BOOLEAN_NAMES),
        "def1_error_map_passed": False,
        "candidate_or_control_trajectory_read": False,
        "global_pde_error_certified": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
        "mechanism_claimed": False,
        "pref1_implemented": False,
        "licenses_separate_independent_pref1_binder": True,
    }
    return config_raw, canonical_result(result)


def expected_config_from_bytes(raw: bytes) -> dict[str, Any]:
    try:
        parsed = tomllib.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Def1Stab1Frz1Error("FRZ1 config is malformed") from error
    return _mapping(parsed, "FRZ1 config")


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    if type(config_raw) is not bytes or type(result_raw) is not bytes:
        _fail("compact bytes differ")
    config_digest = _pinned_digest(CONFIG_SHA256, "config SHA-256")
    result_digest = _pinned_digest(RESULT_SHA256, "result SHA-256")
    if _sha(config_raw) != config_digest:
        _fail("FRZ1 config SHA-256 differs")
    if _sha(result_raw) != result_digest:
        _fail("FRZ1 compact SHA-256 differs")
    config = expected_config_from_bytes(config_raw)
    if render_config(config) != config_raw:
        _fail("FRZ1 config bytes are not canonical")
    try:
        result = json.loads(result_raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise Def1Stab1Frz1Error("FRZ1 compact is not JSON") from error
    result = _mapping(result, "FRZ1 compact")
    try:
        if canonical_result(result) != result_raw:
            _fail("FRZ1 compact bytes are not canonical")
    except CanonicalJSONError as error:
        raise Def1Stab1Frz1Error("FRZ1 compact bytes are not canonical") from error
    _same(config.get("artifact_id"), ARTIFACT_ID, "config artifact")
    _same(result.get("artifact_id"), ARTIFACT_ID, "result artifact")
    _same(config.get("classification"), CLASSIFICATION, "config classification")
    _same(result.get("classification"), CLASSIFICATION, "result classification")
    _same(config.get("base_commit"), BASE_COMMIT, "config base commit")
    _same(result.get("base_commit"), BASE_COMMIT, "result base commit")
    _same(
        config.get("pref1_artifact_id"),
        PREF1_ARTIFACT_ID,
        "config PREF1 license",
    )
    _same(
        result.get("pref1_artifact_id"),
        PREF1_ARTIFACT_ID,
        "result PREF1 license",
    )
    _same(result.get("config_sha256"), config_digest, "embedded config hash")
    claims = _mapping(result.get("claims"), "result claims")
    config_claims = _mapping(config.get("claims"), "config claims")
    for name in TRUE_CLAIMS:
        _boolean(claims.get(name), True, name)
        _boolean(config_claims.get(name), True, f"config.{name}")
        if name in result:
            _boolean(result.get(name), True, f"result.{name}")
    for name in FALSE_CLAIMS:
        _boolean(claims.get(name), False, name)
        _boolean(config_claims.get(name), False, f"config.{name}")
        if name in result:
            _boolean(result.get(name), False, f"result.{name}")
        if name in config:
            _boolean(config.get(name), False, f"config.{name}")
    _same(claims, expected_claims(), "result claims")
    _same(config_claims, claims, "config claims")
    booleans = _mapping(result.get("def1_booleans"), "DEF1 booleans")
    names = result.get("def1_boolean_names")
    if type(names) is not list or any(type(item) is not str for item in names):
        _fail("DEF1 boolean names differ")
    if set(booleans) != set(names):
        _fail("unevaluated DEF1 boolean inventory differs")
    if any(value is not None for value in booleans.values()):
        _fail("a DEF1 boolean was evaluated")
    _boolean(result.get("def1_booleans_evaluated"), False, "def1_booleans_evaluated")
    controls = _mapping(result.get("controls"), "controls")
    required_controls = {
        "conversion_identities",
        "provider_route_coverage",
        "q_error_assembly",
        "misner_sharp",
        "activation",
        "margins",
        "base_to_adm",
        "imp1_to_q",
    }
    if set(controls) != required_controls:
        _fail("control inventory is partial or broadened")
    coverage = _mapping(controls.get("provider_route_coverage"), "route coverage")
    _boolean(coverage.get("positive_controls_passed"), True, "positive route coverage")
    _boolean(
        coverage.get("injected_failures_refused"),
        True,
        "injected-failure route coverage",
    )
    q_assembly = _mapping(controls.get("q_error_assembly"), "Q assembly")
    _boolean(q_assembly.get("measured_q_argument_refused"), True, "measured-Q refusal")
    _boolean(q_assembly.get("used_measured_q"), False, "used_measured_q")
    mass_flux = _mapping(controls.get("misner_sharp"), "Misner-Sharp")
    _boolean(mass_flux.get("typed_veto_control_passed"), True, "mass-flux veto")
    activation = _mapping(controls.get("activation"), "activation")
    _boolean(activation.get("threshold_passed"), True, "activation threshold")
    _boolean(activation.get("control_dominance_passed"), True, "control dominance")
    margins = _mapping(controls.get("margins"), "margins")
    _boolean(margins.get("trappedness_strict_pass"), True, "trappedness pass")
    _boolean(margins.get("complete_q_strict_pass"), True, "complete-Q pass")
    base_to_adm = _mapping(controls.get("base_to_adm"), "BASE-to-ADM")
    _boolean(base_to_adm.get("roundtrip_holds"), True, "BASE-to-ADM roundtrip")
    imp1 = _mapping(controls.get("imp1_to_q"), "IMP1-to-Q")
    _boolean(imp1.get("def1_error_map_passed"), False, "IMP1 error-map flag")
    owner_hashes = _mapping(config.get("owner_hashes"), "owner hashes")
    if len(owner_hashes) != 4:
        _fail("owner-hash inventory is partial")
    for digest in owner_hashes.values():
        _require_digest(digest, "owner hash")
    return result


def verify_compact(repository: Path) -> dict[str, Any]:
    """Validate tracked compact bytes only. No repository, raw, source, or trajectory I/O."""

    root = _root(Path(repository))
    return validate_compact_result(
        _read_tracked(root, CONFIG_PATH),
        _read_tracked(root, RESULT_PATH),
    )


def emit_config_bytes(repository: Path | None = None) -> bytes:
    config_raw, _result_raw = compose_canonical_artifacts(repository)
    return config_raw


def emit_result_bytes(repository: Path | None = None) -> bytes:
    _config_raw, result_raw = compose_canonical_artifacts(repository)
    return result_raw


__all__ = [
    "ARTIFACT_ID",
    "BASE_COMMIT",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "CONFIG_SHA256",
    "Def1Stab1Frz1Error",
    "FALSE_CLAIMS",
    "NONCLAIMS",
    "OWNER_DOCUMENT",
    "PREF1_ARTIFACT_ID",
    "RESULT_PATH",
    "RESULT_SHA256",
    "TRUE_CLAIMS",
    "compose_canonical_artifacts",
    "emit_config_bytes",
    "emit_result_bytes",
    "expected_claims",
    "validate_compact_result",
    "verify_compact",
]
