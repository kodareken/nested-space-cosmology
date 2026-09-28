"""Independent compact binder for candidate-blind DEF1 error-map readiness.

Ordinary ``verify_compact`` reads only the tracked PREF1 config and result.
It does not reconstruct from live owners, inspect Git, open ``runs/``, read a
trajectory, or execute a runner.  ``compose_canonical_artifacts`` is explicit
initial construction: it authenticates the sealed FRZ1 commit and independently
rebuilds the conversion/error-map contract from low-level DEF1 owners.

This module does not import FRZ1 certificate or freeze-contract decision code
and never writes state.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    ReadBlob,
    ResolveCommit,
    UnsafePathError,
    canonical_json_bytes,
    git_read,
    read_regular_file,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-DEF1-STAB1-PREF1"
PARTIAL_INSTRUMENT_ID = "FGC-1-DEF1-STAB1"
FREEZE_ARTIFACT_ID = "FGC-1-DEF1-STAB1-FRZ1"
CLASSIFICATION = (
    "independently_bound_candidate_blind_error_map_readiness_only_no_trajectory"
)
PROJECT_VERSION = "0.11.0"
FREEZE_COMMIT = "32438805126266999a089396f2025e94f8e6bad4"
FREEZE_PARENT = "c4feb941e408a7c72913b49071801b7531183c1a"
CONFIG_PATH = "configs/fgc/fgc-1-def1-stab1-pref1.toml"
RESULT_PATH = "results/fgc-1-def1-stab1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-def1-stab1-pref1.md"
FREEZE_CONFIG_PATH = "configs/fgc/fgc-1-def1-stab1-frz1.toml"
FREEZE_RESULT_PATH = "results/fgc-1-def1-stab1-frz1.json"
FREEZE_CONFIG_SHA256 = (
    "c88676a78cc64d94f9dc7afadb557c1bdf4af8cc09e2fde259d4f3210bc5059f"
)
FREEZE_RESULT_SHA256 = (
    "c2e5c8c2177d7f0bcf414f448131c0abadcab28a5b5c6b2351dd216dddb5c126"
)
CONFIG_SHA256 = "524eddde1757bc871b413e54c962a6f0e4591b990202080648fceb460d742fe3"
RESULT_SHA256 = "b1a536a660294c54fa9dcd746005140ee415309724a59dbfd2470f25f32d20bf"
_SHA_LENGTH = 64
_HEX = frozenset("0123456789abcdef")

OWNER_PATHS = (
    "src/recursive_horizons/fgc/def1_stab1.py",
    "src/recursive_horizons/fgc/def1_geometry_error.py",
    "src/recursive_horizons/fgc/def1_stab1_providers.py",
    "src/recursive_horizons/fgc/def1_stab1_qualification.py",
)

FRZ1_DELTA: tuple[tuple[str, str], ...] = (
    ("M", "CHANGELOG.md"),
    ("M", "PLAN.md"),
    ("M", "configs/fgc/artifact-catalog.json"),
    ("A", "configs/fgc/fgc-1-def1-stab1-frz1.toml"),
    ("M", "docs/active-code-map.md"),
    ("M", "docs/claim-ledger.md"),
    ("A", "docs/fgc-def1-stab1-frz1.md"),
    ("M", "docs/fgc-def1-stab1.md"),
    ("M", "docs/fgc-runtime-matrix.md"),
    ("M", "docs/research-roadmap.md"),
    ("M", "mk/current-foundation.mk"),
    ("M", "results/README.md"),
    ("A", "results/fgc-1-def1-stab1-frz1.json"),
    ("M", "scripts/build_artifact_catalog.py"),
    ("M", "scripts/repo_checks/catalog.py"),
    ("M", "scripts/repo_checks/core.py"),
    ("A", "scripts/reproduce_fgc_def1_stab1_frz1.py"),
    ("A", "src/recursive_horizons/fgc/def1_stab1_frz1_certificate.py"),
    ("A", "tests/test_check_repo_def1_stab1_frz1.py"),
    ("A", "tests/test_fgc_def1_stab1_frz1_certificate.py"),
    ("M", "tests/test_fgc_wave0_instrument_routing.py"),
    ("M", "tests/test_phase_minus1_artifact_catalog.py"),
    ("M", "tests/test_phase_minus1_make_routing.py"),
)

PINNED_REMAINING_GATE_WORK: tuple[str, ...] = (
    "spatial_temporal Richardson and IMP1 remain declared conditional local "
    "estimators under frozen hypotheses, not independently qualified "
    "continuum enclosures",
    "gauge C-to-jet inverse is not derived internally; C-only evidence refuses",
    "affine Gronwall factor is a supplied finite-interval stability premise",
    "interpolation second-order remainder requires enclosed second derivatives "
    "covering the claimed slots",
    "extraction continuous-Q remainder requires a derivative bound on "
    "unsampled affine gaps; positive samples are not that bound",
    "boundary numerical error is not implied by physical causality",
    "conservation-to-Q conversion is not a global flux theorem; mass-flux "
    "veto remains independent",
    "COL1 trajectory binding and DEF1-PREF1 classification are later; this "
    "map gate is not routed after holdout",
    "DEF1_error_map_passed remains false until independent qualification of "
    "these conversions; a universal nonlinear PDE theorem is not a prerequisite",
)

OPERATIONAL_PREMISES: tuple[str, ...] = (
    "conditional Richardson and IMP1 local estimators under frozen hypotheses",
    "caller-supplied C-to-jet inverse or complete gauge extension residual; "
    "C-only evidence refuses",
    "caller-supplied inverse-J times six-residual; absent inverse-J refuses",
    "supplied affine Gronwall finite-interval stability premise",
    "interpolation second-order remainder requires enclosed second derivatives",
    "extraction continuous-Q remainder requires a derivative bound on "
    "unsampled affine gaps",
    "boundary numerical error is not implied by physical causality; typed "
    "refusal on absent guards",
    "conservation-to-Q conversion is not a global flux theorem; mass-flux "
    "veto remains independent",
    "typed refusal on missing, extra, duplicate, negative, or nonfinite premises",
)

LATER_APPLICATION: tuple[str, ...] = (
    "COL1 trajectory binding and DEF1-PREF1 nine-boolean classification",
    "actual trajectory error radii and measured Q",
)

CLOSED_BY_THIS_BINDER: tuple[str, ...] = (
    "independent qualification of conversion/error-map contract for "
    "pre-holdout map readiness",
)

PRESERVED_NONCLAIMS: tuple[str, ...] = (
    "a universal nonlinear PDE theorem is not a prerequisite and is not claimed",
    "no actual trajectory error radius is claimed",
    "map_readiness_only; DEF1_error_map_passed is not a trajectory result",
)

TRUE_CLAIMS = (
    "def1_error_map_passed",
    "map_readiness_only",
    "independently_reconstructed_conversion_identities",
    "independently_reconstructed_provider_route_coverage",
    "independently_reconstructed_q_error_assembly",
    "independently_reconstructed_misner_sharp_enclosure",
    "independently_reconstructed_activation_controls",
    "independently_reconstructed_trappedness_and_complete_q_margins",
    "independently_reconstructed_base_to_adm_roundtrip",
    "independently_reconstructed_imp1_to_q_debit",
    "independently_reconstructed_owner_hashes",
    "independently_reconstructed_remaining_work_taxonomy",
    "freeze_commit_authenticated",
    "freeze_config_result_hashes_authenticated",
    "freeze_bytes_challenged_not_trusted",
    "operational_premises_make_map_ready",
    "typed_refusal_on_absent_inputs",
)
FALSE_CLAIMS = (
    "def1_booleans_evaluated",
    "candidate_or_control_trajectory_read",
    "actual_error_bounds_evaluated",
    "global_pde_error_certified",
    "used_measured_q",
    "richardson_treated_as_global_pde_error",
    "imp1_admission_debit_treated_as_global_pde_error",
    "rob1_passed",
    "holdout_authorized",
    "physical_claimed",
    "mechanism_claimed",
    "universal_pde_theorem",
    "fitted_threshold",
    "measured_q_dependent_error_shrinkage",
    "gr0_case_eligible",
    "sgbl_branch_owned_and_healthy",
)
NONCLAIMS = (
    "PREF1 binds candidate-blind pre-holdout error-map readiness only.",
    "DEF1_error_map_passed=true means map_readiness_only; it is not a "
    "trajectory result.",
    "All nine trajectory DEF1 booleans remain null and unevaluated.",
    "No actual trajectory error radius or measured Q is consumed.",
    "Conditional Richardson/IMP1, caller-supplied C-to-jet/complete residuals, "
    "inverse-J/Gronwall/derivative premises, and typed refusal on absent inputs "
    "make the map operational before holdout.",
    "A universal nonlinear PDE theorem is not claimed and is not a prerequisite.",
    "Holdout, ROB1, mechanism, and physics remain closed.",
    "GR0_case_eligible and SGBL_branch_owned_and_healthy remain separate gates.",
    "PRO20-PREF1 remains the critical frontier; RSRC1 remains the critical-path "
    "next design.",
)

FORBIDDEN_FRZ1_MODULES = (
    "recursive_horizons.fgc.def1_stab1_frz1_certificate",
    "recursive_horizons.fgc.def1_stab1_freeze_contract",
)


class Def1Stab1Pref1Error(ValueError):
    """The independent DEF1-STAB1-PREF1 binder could not be established."""


def _fail(message: str) -> None:
    raise Def1Stab1Pref1Error(message)


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
        raise Def1Stab1Pref1Error("repository root cannot be resolved") from error
    if resolved != path:
        _fail("repository root traverses a symlink")
    return path


def _read_tracked(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise Def1Stab1Pref1Error(f"tracked PREF1 leaf is absent: {relative}") from error


def expected_claims() -> dict[str, bool]:
    claims = {name: True for name in TRUE_CLAIMS}
    claims.update({name: False for name in FALSE_CLAIMS})
    return claims


def remaining_work_taxonomy() -> dict[str, list[str]]:
    return {
        "reconstructed_from_providers": list(PINNED_REMAINING_GATE_WORK),
        "operational_premises": list(OPERATIONAL_PREMISES),
        "later_application_not_map_readiness": list(LATER_APPLICATION),
        "closed_by_this_binder": list(CLOSED_BY_THIS_BINDER),
        "preserved_nonclaims": list(PRESERVED_NONCLAIMS),
    }


def freeze_delta_pairs() -> list[list[str]]:
    return [[status, path] for status, path in FRZ1_DELTA]


def _toml_scalar(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if type(value) is str:
        return json.dumps(value, ensure_ascii=True)
    if type(value) is list:
        return json.dumps(value, ensure_ascii=True, separators=(",", ":"))
    _fail("unsupported PREF1 config scalar")
    raise AssertionError("unreachable")


def render_config(config: Mapping[str, object]) -> bytes:
    mapping = dict(config)
    tables = (
        "owner_hashes",
        "conversion_identity_hashes",
        "remaining_work_taxonomy",
        "freeze_evidence",
        "claims",
    )
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


def _git(root: Path, operation: object) -> object:
    try:
        return git_read(root, operation)  # type: ignore[arg-type]
    except EvidenceIOError as error:
        raise Def1Stab1Pref1Error("immutable freeze query failed") from error


def _authenticate_freeze(root: Path) -> tuple[bytes, bytes]:
    resolved = _git(root, ResolveCommit(FREEZE_COMMIT))
    _same(resolved, FREEZE_COMMIT, "freeze commit identity")
    parents = _git(root, CommitParents(FREEZE_COMMIT))
    if parents != (FREEZE_PARENT,):
        _fail("freeze parent differs")
    delta = _git(root, InspectDelta(commit=FREEZE_COMMIT, against=FREEZE_PARENT))
    if type(delta) is not tuple:
        _fail("freeze delta is not a tuple")
    observed = tuple(sorted((item.status, item.path) for item in delta))
    expected = tuple(sorted(FRZ1_DELTA))
    if observed != expected:
        _fail("freeze authority delta differs")
    config_blob = _git(root, ReadBlob(FREEZE_COMMIT, FREEZE_CONFIG_PATH))
    result_blob = _git(root, ReadBlob(FREEZE_COMMIT, FREEZE_RESULT_PATH))
    if type(config_blob) is not bytes or type(result_blob) is not bytes:
        _fail("freeze compact objects are not bytes")
    if _sha(config_blob) != FREEZE_CONFIG_SHA256:
        _fail("freeze config SHA-256 differs")
    if _sha(result_blob) != FREEZE_RESULT_SHA256:
        _fail("freeze result SHA-256 differs")
    try:
        live_config = read_regular_file(root, FREEZE_CONFIG_PATH)
        live_result = read_regular_file(root, FREEZE_RESULT_PATH)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise Def1Stab1Pref1Error("tracked freeze leaf is absent") from error
    if live_config != config_blob:
        _fail("live freeze config differs from Git object")
    if live_result != result_blob:
        _fail("live freeze result differs from Git object")
    return config_blob, result_blob


def _parse_freeze_config(raw: bytes) -> dict[str, Any]:
    try:
        parsed = tomllib.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Def1Stab1Pref1Error("freeze config is malformed") from error
    return _mapping(parsed, "freeze config")


def _parse_freeze_result(raw: bytes) -> dict[str, Any]:
    try:
        parsed = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise Def1Stab1Pref1Error("freeze compact is not JSON") from error
    return _mapping(parsed, "freeze compact")


def _challenge_freeze_bytes(
    freeze_config: Mapping[str, Any],
    freeze_result: Mapping[str, Any],
    *,
    owner_hashes: Mapping[str, str],
    conversion_hashes: Mapping[str, str],
    routes: list[list[str]],
    mixed_total: Mapping[str, int],
    injected_total: Mapping[str, int],
    imp1_additive: Mapping[str, int],
) -> None:
    _same(freeze_config.get("artifact_id"), FREEZE_ARTIFACT_ID, "freeze config artifact")
    _same(freeze_result.get("artifact_id"), FREEZE_ARTIFACT_ID, "freeze result artifact")
    _same(freeze_config.get("base_commit"), FREEZE_PARENT, "freeze config base commit")
    _same(freeze_result.get("base_commit"), FREEZE_PARENT, "freeze result base commit")
    _same(
        freeze_config.get("pref1_artifact_id"),
        ARTIFACT_ID,
        "freeze config PREF1 license",
    )
    _same(
        freeze_result.get("pref1_artifact_id"),
        ARTIFACT_ID,
        "freeze result PREF1 license",
    )
    _boolean(freeze_config.get("def1_error_map_passed"), False, "freeze config map gate")
    _boolean(freeze_result.get("def1_error_map_passed"), False, "freeze result map gate")
    _boolean(
        freeze_config.get("pref1_implemented"),
        False,
        "freeze config pref1_implemented",
    )
    _boolean(
        freeze_result.get("pref1_implemented"),
        False,
        "freeze result pref1_implemented",
    )
    freeze_owner_hashes = _mapping(freeze_config.get("owner_hashes"), "freeze owner hashes")
    _same(freeze_owner_hashes, dict(owner_hashes), "freeze owner hashes")
    freeze_conversion = _mapping(
        freeze_config.get("conversion_identity_hashes"),
        "freeze conversion hashes",
    )
    _same(freeze_conversion, dict(conversion_hashes), "freeze conversion hashes")
    controls = _mapping(freeze_result.get("controls"), "freeze controls")
    coverage = _mapping(controls.get("provider_route_coverage"), "freeze routes")
    _same(coverage.get("routes"), routes, "freeze route inventory")
    _same(coverage.get("route_count"), 21, "freeze route count")
    q_assembly = _mapping(controls.get("q_error_assembly"), "freeze Q assembly")
    _same(q_assembly.get("mixed_total"), dict(mixed_total), "freeze mixed Q total")
    _same(
        q_assembly.get("injected_total"),
        dict(injected_total),
        "freeze injected Q total",
    )
    imp1 = _mapping(controls.get("imp1_to_q"), "freeze IMP1")
    _same(imp1.get("additive_q"), dict(imp1_additive), "freeze IMP1 additive Q")
    freeze_claims = _mapping(freeze_result.get("claims"), "freeze claims")
    _boolean(
        freeze_claims.get("conversion_instrument_contract_complete"),
        True,
        "freeze conversion-complete label",
    )
    _boolean(
        freeze_claims.get("licenses_separate_independent_pref1_binder"),
        True,
        "freeze PREF1 license label",
    )
    booleans = _mapping(freeze_result.get("def1_booleans"), "freeze DEF1 booleans")
    if any(value is not None for value in booleans.values()):
        _fail("freeze evaluated a DEF1 boolean")


def compose_canonical_artifacts(repository: Path | None = None) -> tuple[bytes, bytes]:
    """Reconstruct PREF1 from Git-authenticated FRZ1 bytes and live DEF1 owners.

    This path is initial construction only.  Ordinary ``--check`` must not
    call it.  It never reads ``runs/`` or trajectory data and never writes.
    """

    if repository is None:
        _fail("live PREF1 construction requires a repository root")
    root = _root(Path(repository))
    freeze_config_raw, freeze_result_raw = _authenticate_freeze(root)
    freeze_config = _parse_freeze_config(freeze_config_raw)
    freeze_result = _parse_freeze_result(freeze_result_raw)

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
        Def1Stab1Error,
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
    from recursive_horizons.fgc.def1_stab1_providers import (
        REMAINING_GATE_WORK,
        Def1Stab1ProviderError,
        gauge_constraint_provider,
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
        owner_file_sha256,
        red1_mhg2_equation_identity_map,
    )
    from recursive_horizons.fgc.exact_interval import Interval
    from recursive_horizons.fgc.spherical_reduction import Jet2

    if tuple(REMAINING_GATE_WORK) != PINNED_REMAINING_GATE_WORK:
        _fail("remaining-work taxonomy drifted from the independent pin")

    equation_map = red1_mhg2_equation_identity_map()
    geometry_map = fo1_adm_geometry_slot_map()
    if equation_map.claims_col1_values or geometry_map.claims_col1_values:
        _fail("conversion identities claimed COL1 values")
    if len(equation_map.pairs) != 6:
        _fail("RED1/MHG2 pair count differs")
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

    coverage = execute_qualification_coverage_matrix()
    if len(coverage.routes) != 21 or len(coverage.routes) != len(PROVIDER_ROUTES):
        _fail("provider-route coverage is incomplete")
    if any(
        not item.positive_control_passed or not item.injected_failure_refused
        for item in coverage.routes
    ):
        _fail("a provider route lacked positive or injected-failure coverage")
    if coverage.trajectory_values_evaluated or coverage.def1_booleans_evaluated:
        _fail("coverage evaluated trajectory values or DEF1 booleans")
    if coverage.used_measured_q or coverage.global_pde_error_certified:
        _fail("coverage promoted a PDE or measured-Q claim")

    zero_assembled = assemble_q_error_budget(zero_error_premises())
    if zero_assembled.exact_total != 0 or zero_assembled.used_measured_q:
        _fail("zero Q-error control used a measured signal")
    try:
        assemble_q_error_budget(zero_error_premises(), measured_q=Q(1, 1000))
    except Def1Stab1Error as error:
        measured_q_refused = "measured Q" in str(error)
    else:
        measured_q_refused = False
    if not measured_q_refused:
        _fail("Q-error assembly accepted a measured-Q argument")
    try:
        assemble_q_error_budget(zero_error_premises()[:-1])
    except Def1Stab1Error as error:
        missing_refused = "missing=" in str(error)
    else:
        missing_refused = False
    if not missing_refused:
        _fail("Q-error assembly accepted a partial premise inventory")
    try:
        gauge_constraint_provider(
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            context="PREF1 independent C-only refusal",
            unit="declared complete-Q control units",
            provenance="PREF1 C-only refusal control",
            constraint_C=(Q(1, 2), Q(1, 3)),
        )
    except Def1Stab1ProviderError as error:
        c_only_refused = "C-only" in str(error)
    else:
        c_only_refused = False
    if not c_only_refused:
        _fail("gauge provider accepted C-only evidence")

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
        context="DEF1-STAB1-PREF1 independent IMP1-to-Q control",
        unit="declared complete-Q control units",
    )
    if (
        imp1.used_measured_q
        or imp1.global_pde_error_certified
        or imp1.def1_error_map_passed
        or imp1.imp1_admission_debit_treated_as_global_pde_error
    ):
        _fail("IMP1-to-Q control promoted a PDE or error-map claim")

    owner_hashes = {path: owner_file_sha256(path) for path in OWNER_PATHS}
    if set(owner_hashes) != set(OWNER_PATHS) or len(owner_hashes) != 4:
        _fail("owner-hash inventory differs")

    independent_routes = [
        [item.component, item.route_id, item.conversion] for item in coverage.routes
    ]
    mixed_total = {
        "numerator": mixed_assembled.exact_total.numerator,
        "denominator": mixed_assembled.exact_total.denominator,
    }
    injected_total = {
        "numerator": injected_assembled.exact_total.numerator,
        "denominator": injected_assembled.exact_total.denominator,
    }
    imp1_additive = {
        "numerator": imp1.additive_q.numerator,
        "denominator": imp1.additive_q.denominator,
    }
    _challenge_freeze_bytes(
        freeze_config,
        freeze_result,
        owner_hashes=owner_hashes,
        conversion_hashes=conversion_hashes,
        routes=independent_routes,
        mixed_total=mixed_total,
        injected_total=injected_total,
        imp1_additive=imp1_additive,
    )

    claims = expected_claims()
    taxonomy = remaining_work_taxonomy()
    freeze_evidence = {
        "artifact_id": FREEZE_ARTIFACT_ID,
        "commit": FREEZE_COMMIT,
        "parent": FREEZE_PARENT,
        "config_sha256": FREEZE_CONFIG_SHA256,
        "result_sha256": FREEZE_RESULT_SHA256,
        "delta": freeze_delta_pairs(),
        "parsed_def1_error_map_passed": False,
        "challenged_against_independent_reconstruction": True,
    }
    config = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "freeze_commit": FREEZE_COMMIT,
        "freeze_parent": FREEZE_PARENT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": FREEZE_ARTIFACT_ID,
        "partial_instrument_id": PARTIAL_INSTRUMENT_ID,
        "def1_error_map_passed": True,
        "map_readiness_only": True,
        "def1_booleans_evaluated": False,
        "candidate_or_control_trajectory_read": False,
        "actual_error_bounds_evaluated": False,
        "global_pde_error_certified": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
        "mechanism_claimed": False,
        "def1_boolean_names": list(DEF1_BOOLEAN_NAMES),
        "error_budget_component_order": list(ERROR_BUDGET_COMPONENTS),
        "owner_hashes": owner_hashes,
        "conversion_identity_hashes": conversion_hashes,
        "remaining_work_taxonomy": taxonomy,
        "freeze_evidence": freeze_evidence,
        "claims": claims,
    }
    config_raw = render_config(config)
    result = {
        "artifact_id": ARTIFACT_ID,
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "freeze_commit": FREEZE_COMMIT,
        "freeze_parent": FREEZE_PARENT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": FREEZE_ARTIFACT_ID,
        "partial_instrument_id": PARTIAL_INSTRUMENT_ID,
        "config_sha256": _sha(config_raw),
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_result_sha256": FREEZE_RESULT_SHA256,
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
                "routes": independent_routes,
            },
            "q_error_assembly": {
                "zero_total": {"numerator": 0, "denominator": 1},
                "mixed_total": mixed_total,
                "injected_total": injected_total,
                "measured_q_argument_refused": True,
                "missing_premises_refused": True,
                "used_measured_q": False,
                "global_pde_error_certified": False,
            },
            "misner_sharp": {
                "minkowski_mass": {"numerator": 0, "denominator": 1},
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
                "additive_q": imp1_additive,
                "conditional_premise": True,
                "def1_error_map_passed": False,
                "used_measured_q": False,
            },
            "typed_refusal": {
                "measured_q_argument_refused": True,
                "missing_premises_refused": True,
                "c_only_gauge_refused": True,
            },
        },
        "remaining_work_taxonomy": taxonomy,
        "freeze_evidence": freeze_evidence,
        "claims": claims,
        "nonclaims": list(NONCLAIMS),
        "def1_boolean_names": list(DEF1_BOOLEAN_NAMES),
        "def1_booleans_evaluated": False,
        "def1_booleans": _unevaluated_booleans(DEF1_BOOLEAN_NAMES),
        "def1_error_map_passed": True,
        "map_readiness_only": True,
        "candidate_or_control_trajectory_read": False,
        "actual_error_bounds_evaluated": False,
        "global_pde_error_certified": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
        "mechanism_claimed": False,
        "gr0_case_eligible": False,
        "sgbl_branch_owned_and_healthy": False,
    }
    return config_raw, canonical_result(result)


def expected_config_from_bytes(raw: bytes) -> dict[str, Any]:
    try:
        parsed = tomllib.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Def1Stab1Pref1Error("PREF1 config is malformed") from error
    return _mapping(parsed, "PREF1 config")


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    if type(config_raw) is not bytes or type(result_raw) is not bytes:
        _fail("compact bytes differ")
    config_digest = _pinned_digest(CONFIG_SHA256, "config SHA-256")
    result_digest = _pinned_digest(RESULT_SHA256, "compact SHA-256")
    if _sha(config_raw) != config_digest:
        _fail("PREF1 config SHA-256 differs")
    if _sha(result_raw) != result_digest:
        _fail("PREF1 compact SHA-256 differs")
    config = expected_config_from_bytes(config_raw)
    if render_config(config) != config_raw:
        _fail("PREF1 config bytes are not canonical")
    try:
        result = json.loads(result_raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise Def1Stab1Pref1Error("PREF1 compact is not JSON") from error
    result = _mapping(result, "PREF1 compact")
    try:
        if canonical_result(result) != result_raw:
            _fail("PREF1 compact bytes are not canonical")
    except CanonicalJSONError as error:
        raise Def1Stab1Pref1Error("PREF1 compact bytes are not canonical") from error
    _same(config.get("artifact_id"), ARTIFACT_ID, "config artifact")
    _same(result.get("artifact_id"), ARTIFACT_ID, "result artifact")
    _same(config.get("classification"), CLASSIFICATION, "config classification")
    _same(result.get("classification"), CLASSIFICATION, "result classification")
    _same(config.get("freeze_commit"), FREEZE_COMMIT, "config freeze commit")
    _same(result.get("freeze_commit"), FREEZE_COMMIT, "result freeze commit")
    _same(config.get("freeze_parent"), FREEZE_PARENT, "config freeze parent")
    _same(result.get("freeze_parent"), FREEZE_PARENT, "result freeze parent")
    _same(
        config.get("predecessor_artifact_id"),
        FREEZE_ARTIFACT_ID,
        "config predecessor",
    )
    _same(
        result.get("predecessor_artifact_id"),
        FREEZE_ARTIFACT_ID,
        "result predecessor",
    )
    _same(result.get("config_sha256"), config_digest, "embedded config hash")
    _same(
        result.get("freeze_config_sha256"),
        FREEZE_CONFIG_SHA256,
        "embedded freeze config hash",
    )
    _same(
        result.get("freeze_result_sha256"),
        FREEZE_RESULT_SHA256,
        "embedded freeze result hash",
    )
    claims = _mapping(result.get("claims"), "result claims")
    config_claims = _mapping(config.get("claims"), "config claims")
    for name in TRUE_CLAIMS:
        _boolean(claims.get(name), True, name)
        _boolean(config_claims.get(name), True, f"config.{name}")
        if name in result:
            _boolean(result.get(name), True, f"result.{name}")
        if name in config:
            _boolean(config.get(name), True, f"config.{name}")
    for name in FALSE_CLAIMS:
        _boolean(claims.get(name), False, name)
        _boolean(config_claims.get(name), False, f"config.{name}")
        if name in result:
            _boolean(result.get(name), False, f"result.{name}")
        if name in config:
            _boolean(config.get(name), False, f"config.{name}")
    _same(claims, expected_claims(), "result claims")
    _same(config_claims, claims, "config claims")
    _boolean(result.get("def1_error_map_passed"), True, "def1_error_map_passed")
    _boolean(result.get("map_readiness_only"), True, "map_readiness_only")
    _boolean(config.get("def1_error_map_passed"), True, "config.def1_error_map_passed")
    _boolean(config.get("map_readiness_only"), True, "config.map_readiness_only")
    booleans = _mapping(result.get("def1_booleans"), "DEF1 booleans")
    names = result.get("def1_boolean_names")
    if type(names) is not list or any(type(item) is not str for item in names):
        _fail("DEF1 boolean names differ")
    if set(booleans) != set(names):
        _fail("unevaluated DEF1 boolean inventory differs")
    if any(value is not None for value in booleans.values()):
        _fail("a DEF1 boolean was evaluated")
    _boolean(result.get("def1_booleans_evaluated"), False, "def1_booleans_evaluated")
    _boolean(
        result.get("actual_error_bounds_evaluated"),
        False,
        "actual_error_bounds_evaluated",
    )
    _boolean(result.get("holdout_authorized"), False, "holdout_authorized")
    _boolean(result.get("physical_claimed"), False, "physical_claimed")
    _boolean(result.get("mechanism_claimed"), False, "mechanism_claimed")
    _boolean(result.get("rob1_passed"), False, "rob1_passed")
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
        "typed_refusal",
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
    _same(coverage.get("route_count"), 21, "route count")
    q_assembly = _mapping(controls.get("q_error_assembly"), "Q assembly")
    _boolean(q_assembly.get("measured_q_argument_refused"), True, "measured-Q refusal")
    _boolean(q_assembly.get("used_measured_q"), False, "used_measured_q")
    _boolean(q_assembly.get("missing_premises_refused"), True, "missing-premise refusal")
    mass_flux = _mapping(controls.get("misner_sharp"), "Misner-Sharp")
    _boolean(mass_flux.get("typed_veto_control_passed"), True, "mass-flux veto")
    activation = _mapping(controls.get("activation"), "activation")
    _boolean(activation.get("threshold_passed"), True, "activation threshold")
    _boolean(activation.get("control_dominance_passed"), True, "control dominance")
    margins = _mapping(controls.get("margins"), "margins")
    _boolean(margins.get("trappedness_strict_pass"), True, "trappedness pass")
    _boolean(margins.get("complete_q_strict_pass"), True, "complete-Q pass")
    _boolean(margins.get("trappedness_boundary_fails"), True, "trappedness boundary")
    _boolean(margins.get("complete_q_boundary_fails"), True, "complete-Q boundary")
    base_to_adm = _mapping(controls.get("base_to_adm"), "BASE-to-ADM")
    _boolean(base_to_adm.get("roundtrip_holds"), True, "BASE-to-ADM roundtrip")
    imp1 = _mapping(controls.get("imp1_to_q"), "IMP1-to-Q")
    _boolean(imp1.get("def1_error_map_passed"), False, "IMP1 error-map flag")
    _boolean(imp1.get("conditional_premise"), True, "IMP1 conditional premise")
    typed = _mapping(controls.get("typed_refusal"), "typed refusal")
    _boolean(typed.get("c_only_gauge_refused"), True, "C-only refusal")
    taxonomy = _mapping(result.get("remaining_work_taxonomy"), "remaining-work taxonomy")
    _same(taxonomy, remaining_work_taxonomy(), "remaining-work taxonomy")
    config_taxonomy = _mapping(
        config.get("remaining_work_taxonomy"),
        "config remaining-work taxonomy",
    )
    _same(config_taxonomy, taxonomy, "config remaining-work taxonomy")
    freeze_evidence = _mapping(result.get("freeze_evidence"), "freeze evidence")
    _same(freeze_evidence.get("commit"), FREEZE_COMMIT, "freeze evidence commit")
    _same(freeze_evidence.get("parent"), FREEZE_PARENT, "freeze evidence parent")
    _same(freeze_evidence.get("delta"), freeze_delta_pairs(), "freeze evidence delta")
    _boolean(
        freeze_evidence.get("challenged_against_independent_reconstruction"),
        True,
        "freeze bytes challenged",
    )
    _boolean(
        freeze_evidence.get("parsed_def1_error_map_passed"),
        False,
        "parsed freeze map gate",
    )
    owner_hashes = _mapping(config.get("owner_hashes"), "owner hashes")
    if len(owner_hashes) != 4:
        _fail("owner-hash inventory is partial")
    for digest in owner_hashes.values():
        _require_digest(digest, "owner hash")
    return result


def verify_compact(repository: Path) -> dict[str, Any]:
    """Validate tracked compact bytes only. No Git, owner, raw, or trajectory I/O."""

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
    "CLASSIFICATION",
    "CONFIG_PATH",
    "CONFIG_SHA256",
    "Def1Stab1Pref1Error",
    "FALSE_CLAIMS",
    "FORBIDDEN_FRZ1_MODULES",
    "FREEZE_COMMIT",
    "FREEZE_CONFIG_SHA256",
    "FREEZE_PARENT",
    "FREEZE_RESULT_SHA256",
    "FRZ1_DELTA",
    "NONCLAIMS",
    "OWNER_DOCUMENT",
    "PINNED_REMAINING_GATE_WORK",
    "RESULT_PATH",
    "RESULT_SHA256",
    "TRUE_CLAIMS",
    "compose_canonical_artifacts",
    "emit_config_bytes",
    "emit_result_bytes",
    "expected_claims",
    "remaining_work_taxonomy",
    "validate_compact_result",
    "verify_compact",
]
