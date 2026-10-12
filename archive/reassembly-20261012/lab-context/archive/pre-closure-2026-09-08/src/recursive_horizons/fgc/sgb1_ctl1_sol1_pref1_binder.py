"""Independent compact binder for the frozen SOL1 scoped outcome."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from recursive_horizons.evidence_io import (
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    ReadBlob,
    UnsafePathError,
    canonical_json_bytes,
    git_read,
    read_regular_file,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1-PREF1"
FREEZE_ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1-FRZ1"
SOL1_ARTIFACT_ID = "FGC-1-SGB1-CTL1-SOL1"
CLASSIFICATION = "independently_bound_sol1_nominal_interval_inconclusive_no_health"
PROJECT_VERSION = "0.11.0"
FREEZE_COMMIT = "a3fdbf5176d7edb5ab90ebe49c136780985e0f3e"
FREEZE_PARENT = "2be811bd7868904e00e2acddac75186cde169408"
FREEZE_CONFIG_PATH = "configs/fgc/fgc-1-sgb1-ctl1-sol1-frz1.toml"
FREEZE_RESULT_PATH = "results/fgc-1-sgb1-ctl1-sol1-frz1.json"
FREEZE_CONFIG_SHA256 = "db2cc1f2dc1bb27b6a39df8c49c9e44bcd90aeeb7029c70f6923c6e7d65596b5"
FREEZE_RESULT_SHA256 = "51923aeabc2c4e6e6f04db6f35812efa9f7675debe9303c2695e4ac4b5f5fa4a"
CONFIG_PATH = "configs/fgc/fgc-1-sgb1-ctl1-sol1-pref1.toml"
RESULT_PATH = "results/fgc-1-sgb1-ctl1-sol1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-sgb1-ctl1-sol1-pref1.md"
CONFIG_SHA256 = "2d820babe3b3141ca0cc29eb20c608bdbf9db074e5f1c7380157ca695846aa7a"
RESULT_SHA256 = "d02452f15f802bf68465a3122617cb34b21c02c913657f94c0316ecb5c071446"
NOMINAL_CONTRACT_SHA256 = (
    "5fc44f7023eb82632309c6c15120f2e62e5e6483911e18589579576b36654f4b"
)
CONTROL_CONTRACT_SHA256 = (
    "661f06f9e6d6fa91d64f1767cf8b8b810c6ce17a4e9bd7ccf8a293b72ed946b9"
)
PINNED_OWNER_PATHS = (
    "src/recursive_horizons/fgc/sgb1_ctl1_sol1.py",
    "docs/fgc-sgb1-ctl1-sol1.md",
    "tests/test_fgc_sgb1_ctl1_sol1.py",
)
OWNER_PATHS = PINNED_OWNER_PATHS
FREEZE_DELTA = (
    ("M", "CHANGELOG.md"),
    ("M", "PLAN.md"),
    ("M", "configs/fgc/artifact-catalog.json"),
    ("A", FREEZE_CONFIG_PATH),
    ("M", "docs/active-code-map.md"),
    ("M", "docs/claim-ledger.md"),
    ("M", "docs/fgc-runtime-matrix.md"),
    ("A", "docs/fgc-sgb1-ctl1-sol1-frz1.md"),
    ("M", "docs/fgc-sgb1-ctl1.md"),
    ("M", "docs/research-roadmap.md"),
    ("M", "mk/current-foundation.mk"),
    ("M", "results/README.md"),
    ("A", FREEZE_RESULT_PATH),
    ("M", "scripts/build_artifact_catalog.py"),
    ("M", "scripts/repo_checks/catalog.py"),
    ("M", "scripts/repo_checks/core.py"),
    ("A", "scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py"),
    ("A", "src/recursive_horizons/fgc/sgb1_ctl1_sol1_frz1_certificate.py"),
    ("A", "tests/test_check_repo_sgb1_ctl1_sol1_frz1.py"),
    ("A", "tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py"),
    ("M", "tests/test_fgc_wave0_instrument_routing.py"),
    ("M", "tests/test_phase_minus1_artifact_catalog.py"),
    ("M", "tests/test_phase_minus1_make_routing.py"),
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
    "freeze_commit_authenticated",
    "freeze_delta_authenticated",
    "freeze_bytes_challenged_not_trusted",
    "independently_reconstructed_policy_and_prefix",
    "independently_reconstructed_predecessor_hashes",
    "independently_reconstructed_nominal_outcome",
    "independently_reconstructed_small_amplitude_control",
    "nominal_interval_inconclusive_bound",
)
FALSE_CLAIMS = (
    "policy_retuned",
    "nominal_family_replaced",
    "on_orbit_compactness_nonpass",
    "on_orbit_jacobian_nonpass",
    "scientific_nonpass",
    "SGBL_branch_owned_and_healthy",
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
    "PREF1 independently binds only the frozen SOL1 scoped outcome.",
    "The nominal endpoint is wider than the frozen tube; no policy is retuned.",
    "The result is interval wrapping, not compactness or Jacobian nonpass.",
    "A_chi=1/8 remains a control and does not replace the nominal family.",
    "A later matched-data redesign must be separately frozen.",
    "Branch health, holdout, execution, trapping, model rejection, and physics remain closed.",
)


class SGBLSOL1PREF1Error(ValueError):
    """The independent SOL1 binder could not be established."""


def _fail(message: str) -> None:
    raise SGBLSOL1PREF1Error(message)


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
        raise SGBLSOL1PREF1Error(f"tracked SOL1-PREF1 leaf is absent: {relative}") from error


def _git(root: Path, operation: object) -> object:
    try:
        return git_read(root, operation)  # type: ignore[arg-type]
    except EvidenceIOError as error:
        raise SGBLSOL1PREF1Error("freeze Git authority cannot be inspected") from error


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


def _authenticate_freeze(root: Path) -> tuple[bytes, bytes, dict[str, str]]:
    parents = _git(root, CommitParents(commit=FREEZE_COMMIT))
    if tuple(parents) != (FREEZE_PARENT,):
        _fail("SOL1-FRZ1 parent differs")
    delta = _git(root, InspectDelta(commit=FREEZE_COMMIT, against=FREEZE_PARENT))
    observed = tuple(sorted((item.status, item.path) for item in delta))
    if observed != tuple(sorted(FREEZE_DELTA)):
        _fail("SOL1-FRZ1 exact delta differs")
    config = _git(root, ReadBlob(commit=FREEZE_COMMIT, path=FREEZE_CONFIG_PATH))
    result = _git(root, ReadBlob(commit=FREEZE_COMMIT, path=FREEZE_RESULT_PATH))
    if type(config) is not bytes or type(result) is not bytes:
        _fail("SOL1-FRZ1 blobs are not bytes")
    if _sha(config) != FREEZE_CONFIG_SHA256 or _sha(result) != FREEZE_RESULT_SHA256:
        _fail("SOL1-FRZ1 compact hashes differ")
    owner_hashes: dict[str, str] = {}
    for relative in OWNER_PATHS:
        blob = _git(root, ReadBlob(commit=FREEZE_PARENT, path=relative))
        if type(blob) is not bytes:
            _fail("SOL1 owner blob is not bytes")
        owner_hashes[relative] = _sha(blob)
    return config, result, owner_hashes


def _challenge_freeze(config_raw: bytes, result_raw: bytes) -> None:
    try:
        config = tomllib.loads(config_raw.decode("ascii"))
        result = json.loads(result_raw.decode("ascii"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError, json.JSONDecodeError) as error:
        raise SGBLSOL1PREF1Error("freeze bytes are malformed") from error
    if config.get("artifact_id") != FREEZE_ARTIFACT_ID:
        _fail("freeze config artifact differs")
    if result.get("artifact_id") != FREEZE_ARTIFACT_ID:
        _fail("freeze result artifact differs")
    nominal = result.get("nominal")
    if type(nominal) is not dict:
        _fail("freeze nominal evidence is absent")
    if (
        nominal.get("classification") != "interval_inconclusive"
        or nominal.get("obstruction") != "prefix_endpoint_not_inside_declared_tube"
        or nominal.get("contract_sha256") != NOMINAL_CONTRACT_SHA256
        or nominal.get("unique_tubes") != 0
    ):
        _fail("freeze nominal evidence was promoted or reclassified")
    if result.get("control", {}).get("contract_sha256") != CONTROL_CONTRACT_SHA256:
        _fail("freeze control evidence differs")


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
    }


def compose_canonical_artifacts(repository: Path) -> tuple[bytes, bytes]:
    """Authenticate FRZ1 and independently reconstruct SOL1. Writes nothing."""

    root = _root(repository)
    if OWNER_PATHS != PINNED_OWNER_PATHS:
        _fail("SOL1 owner inventory differs")
    freeze_config, freeze_result, owner_hashes = _authenticate_freeze(root)
    _challenge_freeze(freeze_config, freeze_result)
    from recursive_horizons.fgc.sgb1_ctl1_sol1 import (
        DECLARED_CONTINUATION_STEPS,
        DECLARED_MAX_BISECTION_DEPTH,
        DECLARED_MAX_RATIONAL_BIT_LENGTH,
        DECLARED_STEP_WIDTH,
        DECLARED_SUPPORT_MAXIMUM,
        DECLARED_TUBE_RADIUS,
        NOMINAL_CHI_AMPLITUDE,
        PREFIX_RADIUS,
        sgbl_bind_orbit_local_picard_lindelof,
        sgbl_orbit_local_small_amplitude_control,
    )

    nominal = sgbl_bind_orbit_local_picard_lindelof()
    control = sgbl_orbit_local_small_amplitude_control()
    policy = {
        "nominal_chi_amplitude": str(NOMINAL_CHI_AMPLITUDE),
        "prefix_radius": str(PREFIX_RADIUS),
        "support_maximum": str(DECLARED_SUPPORT_MAXIMUM),
        "continuation_steps": DECLARED_CONTINUATION_STEPS,
        "step_width": str(DECLARED_STEP_WIDTH),
        "tube_radius": str(DECLARED_TUBE_RADIUS),
        "max_bisection_depth": DECLARED_MAX_BISECTION_DEPTH,
        "max_rational_bit_length": DECLARED_MAX_RATIONAL_BIT_LENGTH,
    }
    if policy != EXPECTED_POLICY:
        _fail("independent SOL1 policy differs")
    if (
        nominal.contract_sha256 != NOMINAL_CONTRACT_SHA256
        or nominal.classification != "interval_inconclusive"
        or nominal.obstruction != "prefix_endpoint_not_inside_declared_tube"
        or any(tube.unique_local_affine_orbit for tube in nominal.tubes)
    ):
        _fail("independent nominal reconstruction differs")
    if control.contract_sha256 != CONTROL_CONTRACT_SHA256:
        _fail("independent control reconstruction differs")
    claims = expected_claims()
    config = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "freeze_commit": FREEZE_COMMIT,
        "freeze_parent": FREEZE_PARENT,
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": FREEZE_ARTIFACT_ID,
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
        "freeze_commit": FREEZE_COMMIT,
        "freeze_parent": FREEZE_PARENT,
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_result_sha256": FREEZE_RESULT_SHA256,
        "freeze_delta": [list(item) for item in FREEZE_DELTA],
        "owner_document": OWNER_DOCUMENT,
        "predecessor_artifact_id": FREEZE_ARTIFACT_ID,
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
        raise SGBLSOL1PREF1Error("SOL1-PREF1 config is malformed") from error
    return parsed


def _parse_result(raw: bytes) -> dict[str, Any]:
    try:
        parsed = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SGBLSOL1PREF1Error("SOL1-PREF1 result is malformed") from error
    if type(parsed) is not dict or canonical_result(parsed) != raw:
        _fail("SOL1-PREF1 result is not canonical")
    return parsed


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    if _sha(config_raw) != _pinned(CONFIG_SHA256, "config SHA-256"):
        _fail("SOL1-PREF1 config hash differs")
    if _sha(result_raw) != _pinned(RESULT_SHA256, "result SHA-256"):
        _fail("SOL1-PREF1 result hash differs")
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    if result.get("config_sha256") != CONFIG_SHA256:
        _fail("SOL1-PREF1 result does not bind config")
    for mapping, label in ((config, "config"), (result, "result")):
        if mapping.get("artifact_id") != ARTIFACT_ID:
            _fail(f"{label} artifact differs")
        if mapping.get("classification") != CLASSIFICATION:
            _fail(f"{label} classification differs")
        if mapping.get("freeze_commit") != FREEZE_COMMIT:
            _fail(f"{label} freeze commit differs")
        if mapping.get("freeze_parent") != FREEZE_PARENT:
            _fail(f"{label} freeze parent differs")
        if mapping.get("predecessor_artifact_id") != FREEZE_ARTIFACT_ID:
            _fail(f"{label} predecessor differs")
        if mapping.get("claims") != expected_claims():
            _fail(f"{label} claims differ")
    if config.get("policy") != [f"{k}={v}" for k, v in EXPECTED_POLICY.items()]:
        _fail("config policy differs")
    if result.get("policy") != EXPECTED_POLICY:
        _fail("result policy differs")
    if config.get("owner_hashes") != result.get("owner_hashes"):
        _fail("owner hashes differ")
    if set(config.get("owner_hashes", ())) != set(PINNED_OWNER_PATHS):
        _fail("owner inventory differs")
    nominal = result.get("nominal")
    control = result.get("control")
    if type(nominal) is not dict or type(control) is not dict:
        _fail("SOL1-PREF1 omits nominal or control evidence")
    if (
        nominal.get("classification") != "interval_inconclusive"
        or nominal.get("obstruction") != "prefix_endpoint_not_inside_declared_tube"
        or nominal.get("contract_sha256") != NOMINAL_CONTRACT_SHA256
        or nominal.get("unique_tubes") != 0
    ):
        _fail("independent nominal evidence differs")
    if control.get("contract_sha256") != CONTROL_CONTRACT_SHA256:
        _fail("independent control evidence differs")
    if result.get("freeze_config_sha256") != FREEZE_CONFIG_SHA256:
        _fail("freeze config hash differs")
    if result.get("freeze_result_sha256") != FREEZE_RESULT_SHA256:
        _fail("freeze result hash differs")
    if result.get("freeze_delta") != [list(item) for item in FREEZE_DELTA]:
        _fail("freeze delta differs")
    if result.get("nonclaims") != list(NONCLAIMS):
        _fail("SOL1-PREF1 nonclaims differ")
    return result


def verify_compact(repository: Path) -> dict[str, Any]:
    """Validate compact bytes only; no Git, source, or reconstruction."""

    root = _root(repository)
    return validate_compact_result(_read(root, CONFIG_PATH), _read(root, RESULT_PATH))


__all__ = [
    "ARTIFACT_ID",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "CONFIG_SHA256",
    "CONTROL_CONTRACT_SHA256",
    "FALSE_CLAIMS",
    "FREEZE_COMMIT",
    "FREEZE_CONFIG_SHA256",
    "FREEZE_DELTA",
    "FREEZE_PARENT",
    "FREEZE_RESULT_SHA256",
    "NOMINAL_CONTRACT_SHA256",
    "OWNER_PATHS",
    "PINNED_OWNER_PATHS",
    "RESULT_PATH",
    "RESULT_SHA256",
    "SGBLSOL1PREF1Error",
    "TRUE_CLAIMS",
    "compose_canonical_artifacts",
    "expected_claims",
    "validate_compact_result",
    "verify_compact",
]
