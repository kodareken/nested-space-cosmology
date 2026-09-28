"""Compact freeze for the prospective SOL1 orbit-local continuation.

Ordinary :func:`verify_compact` reads only the tracked config and result.
Explicit construction imports the committed SOL1 owner, re-executes the
frozen nominal and control paths, and writes nothing itself.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from recursive_horizons.evidence_io import (
    EvidenceIOError,
    UnsafePathError,
    canonical_json_bytes,
    read_regular_file,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1-FRZ1"
SOL1_ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1"
PREF1_ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1-PREF1"
CLASSIFICATION = "prospective_sol1_frozen_nominal_interval_inconclusive_no_health"
PROJECT_VERSION = "0.11.0"
BASE_COMMIT = "2be811b"
CONFIG_PATH = "configs/fgc/fgc-1-sgb1-ctl1-sol1-frz1.toml"
RESULT_PATH = "results/fgc-1-sgb1-ctl1-sol1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-sgb1-ctl1-sol1-frz1.md"
CONFIG_SHA256 = "db2cc1f2dc1bb27b6a39df8c49c9e44bcd90aeeb7029c70f6923c6e7d65596b5"
RESULT_SHA256 = "51923aeabc2c4e6e6f04db6f35812efa9f7675debe9303c2695e4ac4b5f5fa4a"
NOMINAL_CONTRACT_SHA256 = (
    "5fc44f7023eb82632309c6c15120f2e62e5e6483911e18589579576b36654f4b"
)
CONTROL_CONTRACT_SHA256 = (
    "661f06f9e6d6fa91d64f1767cf8b8b810c6ce17a4e9bd7ccf8a293b72ed946b9"
)
OWNER_PATHS = (
    "src/recursive_horizons/fgc/sgb1_ctl1_sol1.py",
    "docs/fgc-sgb1-ctl1-sol1.md",
    "tests/test_fgc_sgb1_ctl1_sol1.py",
)
EXPECTED_POLICY = {
    "nominal_chi_amplitude": "3",
    "prefix_radius": "10749/1024",
    "support_maximum": "14",
    "continuation_steps": 16,
    "step_width": "3587/16384",
    "tube_radius": "1/8",
    "max_bisection_depth": 0,
    "max_rational_bit_length": 16384,
}
TRUE_CLAIMS = (
    "sol1_instrument_contract_frozen",
    "nominal_policy_directly_executed",
    "nominal_interval_inconclusive_preserved",
    "policy_prefix_and_predecessors_bound",
    "small_amplitude_control_directly_executed",
    "licenses_separate_independent_pref1_binder",
)
FALSE_CLAIMS = (
    "policy_retuned",
    "nominal_family_replaced",
    "on_orbit_compactness_nonpass",
    "on_orbit_jacobian_nonpass",
    "scientific_nonpass",
    "SGBL_branch_owned_and_healthy",
    "sol1_pref1_implemented",
    "holdout_authorized",
    "execution_authorized",
    "trajectory_trapping_claimed",
    "actual_orbit_traps",
    "sgbl_model_rejected",
    "spacetime_ibvp_theorem",
    "mechanism_claimed",
    "physical_claimed",
)
NONCLAIMS = (
    "FRZ1 freezes the SOL1 instrument and its scoped nominal outcome only.",
    "The nominal A_chi=3 endpoint does not fit the prospectively declared tube.",
    "That outcome is interval wrapping, not trapping or a Jacobian nonpass.",
    "The tube policy and physical family are not retuned after the outcome.",
    "A_chi=1/8 is a control only and never replaces the nominal family.",
    "Branch health, holdout, execution, model rejection, and physics remain closed.",
    "FRZ1 licenses only a separately implemented independent SOL1-PREF1 binder.",
)


class SGBLSOL1FRZ1Error(ValueError):
    """The SOL1 compact freeze could not be established."""


def _fail(message: str) -> None:
    raise SGBLSOL1FRZ1Error(message)


def _sha(raw: bytes) -> str:
    if type(raw) is not bytes:
        _fail("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _pinned(value: str, label: str) -> str:
    if value.startswith("PLACEHOLDER_"):
        _fail(f"{label} is not pinned")
    if len(value) != 64 or any(ch not in "0123456789abcdef" for ch in value):
        _fail(f"{label} is not a lowercase SHA-256 digest")
    return value


def _root(repository: Path) -> Path:
    root = Path(repository)
    if not root.is_absolute():
        _fail("repository root is not canonical")
    resolved = root.resolve(strict=True)
    if resolved != root:
        _fail("repository root traverses a symlink")
    return root


def _read(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise SGBLSOL1FRZ1Error(f"tracked SOL1-FRZ1 leaf is absent: {relative}") from error


def expected_claims() -> dict[str, bool]:
    claims = {name: True for name in TRUE_CLAIMS}
    claims.update({name: False for name in FALSE_CLAIMS})
    return claims


def _scalar(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if type(value) is str:
        return json.dumps(value, ensure_ascii=True)
    if type(value) is list:
        return json.dumps(value, ensure_ascii=True, separators=(",", ":"))
    _fail("unsupported config scalar")
    raise AssertionError("unreachable")


def render_config(config: Mapping[str, object]) -> bytes:
    mapping = dict(config)
    lines: list[str] = []
    for key, value in mapping.items():
        if key in {"owner_hashes", "claims"}:
            continue
        lines.append(f"{key} = {_scalar(value)}")
    for table in ("owner_hashes", "claims"):
        lines.append(f"\n[{table}]")
        nested = mapping.get(table)
        if type(nested) is not dict:
            _fail(f"{table} is not an object")
        for key, value in nested.items():
            lines.append(f"{json.dumps(key)} = {_scalar(value)}")
    return ("\n".join(lines) + "\n").encode("ascii")


def canonical_result(value: object) -> bytes:
    return canonical_json_bytes(value) + b"\n"


def _summary(record: object) -> dict[str, object]:
    return {
        "classification": record.classification,
        "obstruction": record.obstruction,
        "contract_sha256": record.contract_sha256,
        "tube_count": len(record.tubes),
        "unique_tubes": sum(1 for tube in record.tubes if tube.unique_local_affine_orbit),
        "coverage_left": str(record.coverage_left),
        "coverage_right": str(record.coverage_right),
        "tiles_remaining_support": record.tiles_remaining_support,
        "unique_affine_orbit_on_remaining_support": (
            record.unique_affine_orbit_on_remaining_support
        ),
    }


def compose_canonical_artifacts(repository: Path) -> tuple[bytes, bytes]:
    """Directly reconstruct the prospective compact freeze. No Git or raw data."""

    root = _root(repository)
    from recursive_horizons.fgc.sgb1_ctl1_sol1 import (
        DECLARED_CONTINUATION_STEPS,
        DECLARED_MAX_BISECTION_DEPTH,
        DECLARED_MAX_RATIONAL_BIT_LENGTH,
        DECLARED_STEP_WIDTH,
        DECLARED_SUPPORT_MAXIMUM,
        DECLARED_TUBE_RADIUS,
        NOMINAL_CHI_AMPLITUDE,
        PREFIX_RADIUS,
        SMALL_AMPLITUDE_CONTROL,
        sgbl_bind_orbit_local_picard_lindelof,
        sgbl_orbit_local_small_amplitude_control,
    )

    nominal = sgbl_bind_orbit_local_picard_lindelof()
    control = sgbl_orbit_local_small_amplitude_control()
    if nominal.contract_sha256 != NOMINAL_CONTRACT_SHA256:
        _fail("nominal SOL1 contract hash differs")
    if control.contract_sha256 != CONTROL_CONTRACT_SHA256:
        _fail("small-amplitude control hash differs")
    if (
        nominal.classification != "interval_inconclusive"
        or nominal.obstruction != "prefix_endpoint_not_inside_declared_tube"
        or any(tube.unique_local_affine_orbit for tube in nominal.tubes)
    ):
        _fail("nominal SOL1 branch was reclassified")
    if nominal.slice.chi_amplitude != NOMINAL_CHI_AMPLITUDE:
        _fail("nominal SOL1 family differs")
    if control.slice.chi_amplitude != SMALL_AMPLITUDE_CONTROL:
        _fail("small-amplitude control differs")

    owner_hashes = {relative: _sha(_read(root, relative)) for relative in OWNER_PATHS}
    claims = expected_claims()
    observed_policy = {
        "nominal_chi_amplitude": str(NOMINAL_CHI_AMPLITUDE),
        "prefix_radius": str(PREFIX_RADIUS),
        "support_maximum": str(DECLARED_SUPPORT_MAXIMUM),
        "continuation_steps": DECLARED_CONTINUATION_STEPS,
        "step_width": str(DECLARED_STEP_WIDTH),
        "tube_radius": str(DECLARED_TUBE_RADIUS),
        "max_bisection_depth": DECLARED_MAX_BISECTION_DEPTH,
        "max_rational_bit_length": DECLARED_MAX_RATIONAL_BIT_LENGTH,
    }
    if observed_policy != EXPECTED_POLICY:
        _fail("live SOL1 policy differs from the compact freeze")
    policy = dict(EXPECTED_POLICY)
    config = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "base_commit": BASE_COMMIT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": SOL1_ARTIFACT_ID,
        "pref1_artifact_id": PREF1_ARTIFACT_ID,
        "nominal_contract_sha256": NOMINAL_CONTRACT_SHA256,
        "control_contract_sha256": CONTROL_CONTRACT_SHA256,
        "policy": [f"{key}={value}" for key, value in policy.items()],
        "owner_hashes": owner_hashes,
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
        "predecessor_artifact_id": SOL1_ARTIFACT_ID,
        "pref1_artifact_id": PREF1_ARTIFACT_ID,
        "config_sha256": _sha(config_raw),
        "policy": policy,
        "nominal": _summary(nominal),
        "control": _summary(control),
        "owner_hashes": owner_hashes,
        "claims": claims,
        "nonclaims": list(NONCLAIMS),
    }
    return config_raw, canonical_result(result)


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        parsed = tomllib.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise SGBLSOL1FRZ1Error("SOL1-FRZ1 config is malformed") from error
    if type(parsed) is not dict:
        _fail("SOL1-FRZ1 config is not an object")
    return parsed


def _parse_result(raw: bytes) -> dict[str, Any]:
    try:
        parsed = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SGBLSOL1FRZ1Error("SOL1-FRZ1 result is malformed") from error
    if type(parsed) is not dict or canonical_result(parsed) != raw:
        _fail("SOL1-FRZ1 result is not canonical")
    return parsed


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    if _sha(config_raw) != _pinned(CONFIG_SHA256, "config SHA-256"):
        _fail("SOL1-FRZ1 config hash differs")
    if _sha(result_raw) != _pinned(RESULT_SHA256, "result SHA-256"):
        _fail("SOL1-FRZ1 result hash differs")
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    if result.get("config_sha256") != CONFIG_SHA256:
        _fail("SOL1-FRZ1 result does not bind the config")
    for mapping, label in ((config, "config"), (result, "result")):
        if mapping.get("artifact_id") != ARTIFACT_ID:
            _fail(f"{label} artifact differs")
        if mapping.get("classification") != CLASSIFICATION:
            _fail(f"{label} classification differs")
        if mapping.get("base_commit") != BASE_COMMIT:
            _fail(f"{label} base commit differs")
        if mapping.get("predecessor_artifact_id") != SOL1_ARTIFACT_ID:
            _fail(f"{label} predecessor differs")
        if mapping.get("pref1_artifact_id") != PREF1_ARTIFACT_ID:
            _fail(f"{label} PREF1 routing differs")
        if mapping.get("claims") != expected_claims():
            _fail(f"{label} claims differ")
    if config.get("nominal_contract_sha256") != NOMINAL_CONTRACT_SHA256:
        _fail("config nominal contract hash differs")
    if config.get("control_contract_sha256") != CONTROL_CONTRACT_SHA256:
        _fail("config control contract hash differs")
    expected_policy_lines = [f"{key}={value}" for key, value in EXPECTED_POLICY.items()]
    if config.get("policy") != expected_policy_lines or result.get("policy") != EXPECTED_POLICY:
        _fail("SOL1-FRZ1 policy differs")
    if config.get("owner_hashes") != result.get("owner_hashes"):
        _fail("SOL1-FRZ1 owner hashes differ between config and result")
    if set(config.get("owner_hashes", ())) != set(OWNER_PATHS):
        _fail("SOL1-FRZ1 owner inventory differs")
    nominal = result.get("nominal")
    control = result.get("control")
    if type(nominal) is not dict or type(control) is not dict:
        _fail("SOL1-FRZ1 result omits nominal or control evidence")
    if (
        nominal.get("classification") != "interval_inconclusive"
        or nominal.get("obstruction") != "prefix_endpoint_not_inside_declared_tube"
        or nominal.get("contract_sha256") != NOMINAL_CONTRACT_SHA256
        or nominal.get("unique_tubes") != 0
        or nominal.get("tiles_remaining_support") is not False
    ):
        _fail("nominal SOL1 evidence differs")
    if control.get("contract_sha256") != CONTROL_CONTRACT_SHA256:
        _fail("SOL1 control evidence differs")
    if result.get("nonclaims") != list(NONCLAIMS):
        _fail("SOL1-FRZ1 nonclaims differ")
    return result


def verify_compact(repository: Path) -> dict[str, Any]:
    """Validate tracked bytes only; no Git, source reconstruction, or raw data."""

    root = _root(repository)
    return validate_compact_result(_read(root, CONFIG_PATH), _read(root, RESULT_PATH))


__all__ = [
    "ARTIFACT_ID",
    "BASE_COMMIT",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "CONFIG_SHA256",
    "CONTROL_CONTRACT_SHA256",
    "FALSE_CLAIMS",
    "EXPECTED_POLICY",
    "NOMINAL_CONTRACT_SHA256",
    "OWNER_DOCUMENT",
    "OWNER_PATHS",
    "PREF1_ARTIFACT_ID",
    "RESULT_PATH",
    "RESULT_SHA256",
    "SGBLSOL1FRZ1Error",
    "TRUE_CLAIMS",
    "compose_canonical_artifacts",
    "expected_claims",
    "validate_compact_result",
    "verify_compact",
]
