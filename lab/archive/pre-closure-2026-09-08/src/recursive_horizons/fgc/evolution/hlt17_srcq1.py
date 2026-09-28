"""Prospective HLT17 physical-source/origin/resource qualification core.

This is implementation A: a deterministic no-write harness. It does not
publish a campaign store, adopt an endpoint, or advance accepted state.
Ordinary construction and status inspection do not call the physical source.
The one real six-member measurement is a later coordinator-owned freeze-B
invocation of the explicit CLI, after exact-delta authority is bound.

Do not import historical versioned GR-0 calibration runners. The preserved
static factory is the only GR-0 construction path.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
from hashlib import sha256
import importlib
import math
from pathlib import Path
import platform
import re
import resource
import sys
import time
from types import MappingProxyType
from typing import Any, Callable, Mapping

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    PostpublicationUncertainty,
    PrepublicationError,
    UnsafePathError,
    canonical_json_bytes,
    load_canonical_json,
    publish_exclusive_directory,
    read_regular_file,
)

from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RuntimeResourceStop,
    hlt17_implementation_identity,
    require_hlt17_implementation_identity,
    require_hlt17_prepared,
)
from .hlt17_imp1_cursor import (
    HLT17_CFL_RETRY_CAP,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SOURCE_RETRY_CAP,
    encode_tdg7_plan,
    production_member,
)
from .numerical_engine import EvolutionRHS, EvolutionState, array_content_sha256
from .proto5_runtime import CFLRetryRequired
from .pro20_origin import (
    MEMBER_KEYS,
    ORIGIN_TIME_HEX,
    Pro20OriginCapture,
    revalidate_pro20_origin_capture,
)
from .pro20_source_factory import (
    CAMPAIGN_ID,
    EVENT_TARGET,
    PROTOCOL,
    build_pro20_runtime_origin,
)
from .proto19_gr0_static_factory import STATIC_INPUT_PATHS
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS
from .tdg11_imp1_ledger import inherited_tdg6_ledger


SCHEMA = "FGC-1-HLT17-SRCQ1-qualification-v1"
ARTIFACT_ID = "FGC-1-HLT17-SRCQ1"
IMPLEMENTATION_RELATIVE = "src/recursive_horizons/fgc/evolution/hlt17_srcq1.py"
RESULT_LEAF = "qualification.json"
AUTHORITY_VALIDATOR_MODULE = "recursive_horizons.fgc.evolution.hlt17_srcq1_auth1"
AUTHORITY_VALIDATOR_NAME = "validate_authority_delta"
AUTHORITY_SEAM = (
    "authority/delta validation is the coordinator freeze-B seam "
    f"{AUTHORITY_VALIDATOR_MODULE}.{AUTHORITY_VALIDATOR_NAME}; "
    "implementation A does not accept a boolean skip/force flag"
)
EXPECTED_SPATIAL_ORDER = MappingProxyType({"RK4": 4, "SSPRK3": 2})
STATIC_INPUT_SHA256 = MappingProxyType(
    {
        "configs/fgc/fgc-1-cal9-run1.toml": (
            "8e6fafff638ffd49549efab8b6f969a01964c915116db9ce3698953201172176"
        ),
        "configs/fgc/fgc-1-rsp2-run1.toml": (
            "8213341f6cc5dd1d044a1ee9ae68b9ce7c7445ef9b933bcf8bcc944d30e4abd1"
        ),
        "results/fgc-1-hlt10-mon10.json": (
            "630d0d1843130180cb49dd237fe32e83e0336e7e44e7c684e4baf406b1883fe5"
        ),
        "results/fgc-1-rsp2-frz1.json": (
            "a6c147d7b369387e70318ee42eb56c97332702a227fd93d7ad20a1cfe2c1aa22"
        ),
    }
)
SOURCE_DIAGNOSTIC_KEYS = (
    "source_residual_infinity",
    "source_raw_gate_passed",
    "source_refinement_iterations",
    "source_residual_decreased_monotonically",
    "kinetic_condition_infinity",
    "reduction_constraint_infinity",
    "coordinate_speed_upper",
    "minimum_lapse",
    "minimum_radial_metric",
    "minimum_areal_radius_away_from_center",
)
FORBIDDEN_ENDPOINT_KEYS = frozenset(
    {
        "corrected_endpoint",
        "shadow_endpoint",
        "fine_endpoint",
        "outer_endpoint",
        "medium_endpoint",
        "fine_endpoint_state",
        "shadow_endpoint_state",
        "corrected_endpoint_state",
        "committed_state",
        "adopted_endpoint",
        "endpoint_arrays",
        "physical_endpoint",
    }
)
_COMMIT = re.compile(r"\A[0-9a-f]{40}(?:[0-9a-f]{24})?\Z")
_SHA = re.compile(r"\A[0-9a-f]{64}\Z")
GiB = 1024 * 1024 * 1024


class HLT17SRCQ1Error(ValueError):
    """Typed SRCQ1 stop; not a campaign classification."""

    outcome = "invalid"
    published = False


class HLT17SRCQ1CapError(HLT17SRCQ1Error):
    """Prospective initial cap inputs or authority mapping differ."""

    outcome = "invalid"


class HLT17SRCQ1SourceError(HLT17SRCQ1Error):
    outcome = "source"


class HLT17SRCQ1ResourceError(HLT17SRCQ1Error):
    outcome = "resource"


class HLT17SRCQ1IdentityError(HLT17SRCQ1Error):
    outcome = "invalid"


class HLT17SRCQ1AuthoritySeamError(HLT17SRCQ1Error):
    outcome = "invalid"


class HLT17SRCQ1PublicationError(HLT17SRCQ1Error):
    outcome = "invalid"

    def __init__(self, message: str, *, published: bool = False) -> None:
        super().__init__(message)
        self.published = published
        if published:
            self.outcome = "inconclusive"


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise HLT17SRCQ1Error("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if type(value) is not str or _SHA.fullmatch(value) is None:
        raise HLT17SRCQ1Error(f"{label} must be a lowercase SHA-256 digest")
    return value


def _require_commit(value: object, label: str) -> str:
    if type(value) is not str or _COMMIT.fullmatch(value) is None:
        raise HLT17SRCQ1AuthoritySeamError(
            f"{label} must be an exact lowercase 40- or 64-digit commit"
        )
    return value


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or type(value) is bool:
        raise HLT17SRCQ1Error(f"{label} must be a mapping")
    return dict(value)


def recommended_resource_ceilings() -> "ResourceCeilings":
    """Conservative resource guards from preserved controls, not fitted thresholds.

    Final numbers remain coordinator freeze-B values. Retry/enclosure ceilings
    stay the already frozen 32 and 8N/16N production tables.
    """

    return ResourceCeilings(
        max_rss_bytes=4 * GiB,
        max_wall_seconds_per_member=MappingProxyType(
            {key: 900.0 for key in MEMBER_KEYS}
        ),
        max_total_wall_seconds=3600.0,
        coordinator_freeze_required=True,
    )


@dataclass(frozen=True, slots=True)
class ResourceCeilings:
    """Process resource guards. These are not scientific tolerances."""

    max_rss_bytes: int
    max_wall_seconds_per_member: Mapping[str, float]
    max_total_wall_seconds: float
    coordinator_freeze_required: bool = True

    def __post_init__(self) -> None:
        if type(self.max_rss_bytes) is not int or self.max_rss_bytes <= 0:
            raise HLT17SRCQ1ResourceError("max_rss_bytes must be a positive integer")
        if type(self.max_total_wall_seconds) is not float or not (
            math.isfinite(self.max_total_wall_seconds) and self.max_total_wall_seconds > 0.0
        ):
            raise HLT17SRCQ1ResourceError("max_total_wall_seconds must be a positive finite float")
        if type(self.coordinator_freeze_required) is not bool:
            raise TypeError("coordinator_freeze_required must be a built-in bool")
        walls = _mapping(self.max_wall_seconds_per_member, "max_wall_seconds_per_member")
        if tuple(walls) != MEMBER_KEYS:
            raise HLT17SRCQ1ResourceError("per-member wall ceilings are not in six-member order")
        checked: dict[str, float] = {}
        for key in MEMBER_KEYS:
            value = walls[key]
            if type(value) is not float or not math.isfinite(value) or value <= 0.0:
                raise HLT17SRCQ1ResourceError(f"{key} wall ceiling must be a positive finite float")
            checked[key] = value
        object.__setattr__(self, "max_wall_seconds_per_member", MappingProxyType(checked))

    def as_mapping(self) -> dict[str, object]:
        return {
            "max_rss_bytes": self.max_rss_bytes,
            "max_wall_seconds_per_member": {
                key: self.max_wall_seconds_per_member[key] for key in MEMBER_KEYS
            },
            "max_total_wall_seconds": self.max_total_wall_seconds,
            "coordinator_freeze_required": self.coordinator_freeze_required,
            "scientific_threshold_fitted": False,
            "source_retry_cap": HLT17_SOURCE_RETRY_CAP,
            "cfl_retry_cap": HLT17_CFL_RETRY_CAP,
            "production_enclosure_limits": [
                {
                    "member_key": item.member_key,
                    "point_count": item.point_count,
                    "owned_rows": item.owned_rows,
                    "fallback_d01": item.fallback_d01,
                    "fallback_d12": item.fallback_d12,
                }
                for item in HLT17_PRODUCTION_MEMBERS
            ],
        }


def observe_environment() -> dict[str, str]:
    """Observe the live arithmetic image. Status/import never call this."""

    import numpy as np

    configuration = np.show_config(mode="dicts")
    if type(configuration) is not dict:
        raise HLT17SRCQ1IdentityError("NumPy build configuration is unavailable")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if type(dependencies) is dict else {}
    if type(blas) is not dict:
        raise HLT17SRCQ1IdentityError("NumPy BLAS configuration is unavailable")
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": str(np.__version__),
        "blas_name": str(blas.get("name") or ""),
        "blas_version": str(blas.get("version") or "unknown"),
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "executable": sys.executable,
    }
    for key, value in observed.items():
        if type(value) is not str or not value:
            raise HLT17SRCQ1IdentityError(f"environment {key} is empty")
    return observed


def current_rss_bytes() -> int:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = int(usage.ru_maxrss)
    if rss < 0:
        raise HLT17SRCQ1ResourceError("RSS observation is negative")
    if sys.platform == "darwin":
        return rss
    return rss * 1024


def implementation_identity(repository_root: Path | None = None) -> dict[str, object]:
    identity = hlt17_implementation_identity()
    record: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "path": IMPLEMENTATION_RELATIVE,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "c1r1_actual": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "c1r1_is_not_imp1_execution": True,
        "campaign_execution_authorized": False,
    }
    if identity["implementation_id"] == identity["reference_wire_artifact_id"]:
        raise HLT17SRCQ1IdentityError("C1R1 actual identity collapsed into the IMP1 wire")
    if repository_root is not None:
        raw = read_regular_file(Path(repository_root), IMPLEMENTATION_RELATIVE)
        record["sha256"] = _sha256_hex(raw)
    return record


def plan_status() -> dict[str, object]:
    """Default CLI/status record. Performs no I/O, source call, or write."""

    identity = hlt17_implementation_identity()
    ceilings = recommended_resource_ceilings()
    return {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "role": "implementation_A_no_write_qualification_core",
        "physical_six_member_qualification_executed": False,
        "campaign_execution_authorized": False,
        "store_publication_authorized": False,
        "endpoint_adoption_authorized": False,
        "authority_delta_validation": AUTHORITY_SEAM,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "event_target_hex": EVENT_TARGET.hex(),
        "origin_time_hex": ORIGIN_TIME_HEX,
        "prospective_requested_cap_formula": (
            "cfl_maximum*grid_spacing/inherited_previous_speed_upper"
        ),
        "member_keys": list(MEMBER_KEYS),
        "static_input_paths": list(STATIC_INPUT_PATHS),
        "static_input_sha256": dict(STATIC_INPUT_SHA256),
        "c1r1_actual": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "resource_ceilings": ceilings.as_mapping(),
        "one_real_run": {
            "requires": ["exact --authority-commit", "absent --output-directory"],
            "forbids": [
                "boolean skip/force authority flag",
                "existing namespace",
                "cap retune after a result",
                "informal retry",
                "endpoint serialization",
                "campaign store publication",
            ],
            "authority_module": AUTHORITY_VALIDATOR_MODULE,
            "authority_symbol": AUTHORITY_VALIDATOR_NAME,
        },
    }


def load_coordinator_authority_delta_validator() -> Callable[..., Mapping[str, object]]:
    """Coordinator freeze-B seam. Missing module is a typed stop, not a skip."""

    try:
        module = importlib.import_module(AUTHORITY_VALIDATOR_MODULE)
    except ImportError as error:
        raise HLT17SRCQ1AuthoritySeamError(AUTHORITY_SEAM) from error
    validator = getattr(module, AUTHORITY_VALIDATOR_NAME, None)
    if not callable(validator):
        raise HLT17SRCQ1AuthoritySeamError(AUTHORITY_SEAM)
    return validator


def require_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
    validator: Callable[..., object] | None = None,
) -> dict[str, object]:
    commit = _require_commit(authority_commit, "authority_commit")
    digest = _require_sha(implementation_sha256, "implementation_sha256")
    loaded = validator
    if loaded is None:
        loaded = load_coordinator_authority_delta_validator()
    if not callable(loaded):
        raise HLT17SRCQ1AuthoritySeamError(AUTHORITY_SEAM)
    receipt = loaded(
        repository_root=Path(repository_root),
        authority_commit=commit,
        implementation_sha256=digest,
    )
    if type(receipt) is bool:
        raise HLT17SRCQ1AuthoritySeamError(
            "authority/delta validation must not be a boolean flag"
        )
    mapping = _mapping(receipt, "authority delta receipt")
    mapping.setdefault("authority_commit", commit)
    mapping.setdefault("implementation_sha256", digest)
    mapping["boolean_flag_accepted"] = False
    return mapping


def read_static_factory_inputs(repository_root: Path) -> dict[str, bytes]:
    """Safe no-follow reads of the four frozen static factory inputs."""

    root = Path(repository_root)
    payloads: dict[str, bytes] = {}
    for relative in STATIC_INPUT_PATHS:
        try:
            raw = read_regular_file(root, relative)
        except UnsafePathError as error:
            raise HLT17SRCQ1Error(f"static factory input is absent or unsafe: {relative}") from error
        digest = _sha256_hex(raw)
        expected = STATIC_INPUT_SHA256[relative]
        if digest != expected:
            raise HLT17SRCQ1IdentityError(f"static factory input identity differs: {relative}")
        payloads[relative] = raw
    if tuple(payloads) != STATIC_INPUT_PATHS:
        raise HLT17SRCQ1IdentityError("static factory input order differs")
    return payloads


def static_input_identities(payloads: Mapping[str, bytes]) -> dict[str, str]:
    mapping = _mapping(payloads, "static inputs")
    if tuple(mapping) != STATIC_INPUT_PATHS:
        raise HLT17SRCQ1IdentityError("static factory input inventory differs")
    identities: dict[str, str] = {}
    for relative in STATIC_INPUT_PATHS:
        raw = mapping[relative]
        if type(raw) is not bytes:
            raise HLT17SRCQ1Error(f"static factory input is not bytes: {relative}")
        digest = _sha256_hex(raw)
        if digest != STATIC_INPUT_SHA256[relative]:
            raise HLT17SRCQ1IdentityError(f"static factory input identity differs: {relative}")
        identities[relative] = digest
    return identities


def _require_requested_cap_mapping(value: object) -> dict[str, float]:
    mapping = _mapping(value, "prospective requested caps")
    if tuple(mapping) != MEMBER_KEYS:
        raise HLT17SRCQ1CapError("prospective requested-cap order/inventory differs")
    result: dict[str, float] = {}
    remaining = EVENT_TARGET - float.fromhex(ORIGIN_TIME_HEX)
    for key in MEMBER_KEYS:
        text = mapping[key]
        if type(text) is not str or text.lower() != text:
            raise HLT17SRCQ1CapError(f"{key} requested cap must be canonical binary64 hex")
        try:
            cap = float.fromhex(text)
        except ValueError as error:
            raise HLT17SRCQ1CapError(f"{key} requested cap is malformed") from error
        if not math.isfinite(cap) or cap <= 0.0 or cap > remaining or cap.hex() != text:
            raise HLT17SRCQ1CapError(f"{key} requested cap is outside its frozen domain")
        result[key] = cap
    return result


def derive_prospective_requested_caps(
    origin: object,
    runtime: object,
    *,
    expected_cap_hex_by_member: Mapping[str, str] | None = None,
) -> tuple[dict[str, float], dict[str, dict[str, str]]]:
    """Derive initial caps from authenticated causal speed and frozen CFL/grid.

    The HLT15 FRESH_READY cursor has no executed plan. This is a new
    prospective control, not a fabricated historical plan or a remaining-
    interval fallback. TDG7 separately chooses the executed lattice width.
    """

    members = getattr(runtime, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        raise HLT17SRCQ1CapError("runtime origin member order differs for cap derivation")
    captured_members = getattr(origin, "members", None)
    if not isinstance(captured_members, (tuple, list)):
        raise HLT17SRCQ1CapError("captured origin members differ for cap derivation")
    captured = {item.member_key: item for item in captured_members}
    if tuple(captured) != MEMBER_KEYS:
        raise HLT17SRCQ1CapError("captured origin member order differs for cap derivation")
    expected = (
        None
        if expected_cap_hex_by_member is None
        else _require_requested_cap_mapping(expected_cap_hex_by_member)
    )
    caps: dict[str, float] = {}
    facts: dict[str, dict[str, str]] = {}
    for key in MEMBER_KEYS:
        try:
            causal = load_canonical_json(captured[key].causal_bytes)
        except (CanonicalJSONError, TypeError, ValueError) as error:
            raise HLT17SRCQ1CapError(f"{key} retained causal bytes differ") from error
        causal = _mapping(causal, f"{key} retained causal state")
        speed = causal.get("previous_speed_upper")
        if isinstance(speed, bool) or not isinstance(speed, (int, float)):
            raise HLT17SRCQ1CapError(f"{key} inherited causal speed differs")
        previous_speed = float(speed)
        view = members[key]
        transaction = view.member.transaction
        spacing = float(transaction.grid_spacing)
        cfl = float(transaction.cfl_maximum)
        if not all(math.isfinite(item) and item > 0.0 for item in (previous_speed, spacing, cfl)):
            raise HLT17SRCQ1CapError(f"{key} cap inputs are outside their positive finite domain")
        cap = cfl * spacing / previous_speed
        if not math.isfinite(cap) or cap <= 0.0:
            raise HLT17SRCQ1CapError(f"{key} derived requested cap differs")
        if expected is not None and cap.hex() != expected[key].hex():
            raise HLT17SRCQ1CapError(f"{key} authority requested cap differs from retained inputs")
        caps[key] = cap
        facts[key] = {
            "formula": "cfl_maximum*grid_spacing/inherited_previous_speed_upper",
            "cfl_maximum_hex": cfl.hex(),
            "grid_spacing_hex": spacing.hex(),
            "inherited_previous_speed_upper_hex": previous_speed.hex(),
            "requested_cap_hex": cap.hex(),
        }
    return caps, facts


def refuse_endpoint_payload(value: object, *, path: str = "result") -> None:
    """Refuse serialized corrected/shadow endpoints. Flags named *_committable stay."""

    try:
        import numpy as np
    except ImportError:
        np = None  # type: ignore[assignment]
    if np is not None and isinstance(value, np.ndarray):
        raise HLT17SRCQ1Error(f"{path} serializes an array endpoint")
    if type(value) is bytes or type(value) is bytearray:
        raise HLT17SRCQ1Error(f"{path} serializes raw endpoint bytes")
    if isinstance(value, Mapping):
        for key, item in value.items():
            if type(key) is not str:
                raise HLT17SRCQ1Error(f"{path} has a non-text key")
            if key in FORBIDDEN_ENDPOINT_KEYS:
                raise HLT17SRCQ1Error(f"{path}.{key} serializes a forbidden endpoint")
            if key.endswith("_endpoint") and not key.endswith("_committable"):
                raise HLT17SRCQ1Error(f"{path}.{key} serializes a forbidden endpoint")
            refuse_endpoint_payload(item, path=f"{path}.{key}")
        return
    if type(value) is list:
        for index, item in enumerate(value):
            refuse_endpoint_payload(item, path=f"{path}[{index}]")


def _identity_hash(value: object, label: str) -> str:
    if hasattr(value, "as_mapping") and callable(value.as_mapping):
        payload = value.as_mapping()
    elif is_dataclass(value) and not isinstance(value, type):
        payload = asdict(value)
    elif isinstance(value, Mapping):
        payload = dict(value)
    else:
        raise HLT17SRCQ1Error(f"{label} cannot be fingerprinted")
    return _sha256_hex(canonical_json_bytes(payload))


def _tracer_fingerprint(tracers: object) -> str:
    positions = getattr(tracers, "positions", None)
    proper_times = getattr(tracers, "proper_times", None)
    if positions is None or proper_times is None:
        raise HLT17SRCQ1Error("tracer fingerprint requires positions and proper times")
    return array_content_sha256(positions, proper_times)


def member_fingerprints(member: object, binding: object) -> dict[str, object]:
    state = member.state
    if type(state) is not EvolutionState:
        raise HLT17SRCQ1Error("member state type differs")
    transaction = member.transaction
    configuration = getattr(binding, "configuration_sha256", None)
    inherited = inherited_tdg6_ledger(member.temporal_ledger)
    return {
        "physical_state_sha256": array_content_sha256(state.u, state.p, state.q),
        "time_hex": float(member.time).hex(),
        "step_index": int(member.step_index),
        "transaction_serial": int(member.transaction_serial),
        "source_retry_count": int(member.source_retry_count),
        "CFL_retry_count": int(member.CFL_retry_count),
        "source_configuration_sha256": _require_sha(
            configuration, "source configuration"
        ),
        "monitor_sha256": _identity_hash(transaction.state, "monitor"),
        "causal_sha256": _identity_hash(transaction.causal_state, "causal"),
        "tracer_sha256": _tracer_fingerprint(member.tracers),
        "tdg6_inherited_time_hex": float(inherited.last_accepted_time).hex(),
    }


def check_method_owned_identity(member_key: str, member: object) -> dict[str, object]:
    spec = production_member(member_key)
    identity = member.identity
    method_label = identity.method_label
    spatial_order = int(identity.spatial_order)
    expected_order = EXPECTED_SPATIAL_ORDER[method_label]
    if identity.point_count != spec.point_count:
        raise HLT17SRCQ1IdentityError(f"{member_key} point_count differs")
    if identity.integrator_id != spec.method:
        raise HLT17SRCQ1IdentityError(f"{member_key} integrator differs")
    if spatial_order != expected_order:
        raise HLT17SRCQ1IdentityError(f"{member_key} spatial operator differs")
    if method_label == "SSPRK3" and spatial_order == 4:
        raise HLT17SRCQ1IdentityError("SSPRK3-on-SBP4 hybrid is forbidden")
    return {
        "member_key": member_key,
        "method_label": method_label,
        "integrator_id": identity.integrator_id,
        "point_count": identity.point_count,
        "owned_rows": spec.owned_rows,
        "spatial_order": spatial_order,
        "fallback_d01": spec.fallback_d01,
        "fallback_d12": spec.fallback_d12,
    }


def check_inherited_counters(member_key: str, member: object, captured: object) -> None:
    expected = {
        "step_index": int(captured.step_index),
        "transaction_serial": int(captured.transaction_serial),
        "source_retry_count": int(captured.source_retry_count),
        "CFL_retry_count": int(captured.cfl_retry_count),
    }
    actual = {
        "step_index": int(member.step_index),
        "transaction_serial": int(member.transaction_serial),
        "source_retry_count": int(member.source_retry_count),
        "CFL_retry_count": int(member.CFL_retry_count),
    }
    if actual != expected:
        raise HLT17SRCQ1IdentityError(f"{member_key} inherited counters differ")
    inherited = inherited_tdg6_ledger(member.temporal_ledger)
    if float(inherited.last_accepted_time).hex() != ORIGIN_TIME_HEX:
        raise HLT17SRCQ1IdentityError(f"{member_key} inherited TDG6 time differs")
    accepted = captured.accepted_time
    if isinstance(accepted, Mapping) and accepted.get("binary64_hex") != ORIGIN_TIME_HEX:
        raise HLT17SRCQ1IdentityError(f"{member_key} captured origin time differs")


def evaluate_bound_source_once(
    *,
    binding: object,
    time: float,
    state: EvolutionState,
    transaction: object,
    tracers: object,
) -> dict[str, object]:
    """One origin-time RHS evaluation. Input/monitor/causal/tracer must not change."""

    if hasattr(binding, "validate"):
        binding.validate()
    if type(state) is not EvolutionState:
        raise HLT17SRCQ1SourceError("source evaluation state type differs")
    before_state = array_content_sha256(state.u, state.p, state.q)
    before_monitor = _identity_hash(transaction.state, "monitor")
    before_causal = _identity_hash(transaction.causal_state, "causal")
    before_tracer = _tracer_fingerprint(tracers)
    if hasattr(binding, "configuration_sha256"):
        before_config = binding.configuration_sha256
    else:
        before_config = None
    rhs = binding.rhs(time, state)
    if type(rhs) is not EvolutionRHS:
        raise HLT17SRCQ1SourceError("bound source did not return EvolutionRHS")
    if array_content_sha256(state.u, state.p, state.q) != before_state:
        raise HLT17SRCQ1SourceError("bound source mutated its input state")
    if _identity_hash(transaction.state, "monitor") != before_monitor:
        raise HLT17SRCQ1SourceError("bound source mutated monitor state")
    if _identity_hash(transaction.causal_state, "causal") != before_causal:
        raise HLT17SRCQ1SourceError("bound source mutated causal state")
    if _tracer_fingerprint(tracers) != before_tracer:
        raise HLT17SRCQ1SourceError("bound source mutated tracer state")
    if before_config is not None and binding.configuration_sha256 != before_config:
        raise HLT17SRCQ1SourceError("source configuration changed during evaluation")
    if hasattr(binding, "validate"):
        binding.validate()
    diagnostics = dict(rhs.diagnostics)
    missing = [name for name in SOURCE_DIAGNOSTIC_KEYS if name not in diagnostics]
    if missing:
        raise HLT17SRCQ1SourceError(
            "source diagnostics are incomplete: " + ",".join(missing)
        )
    raw_gate = diagnostics["source_raw_gate_passed"]
    if type(raw_gate) is not bool or raw_gate is not True:
        raise HLT17SRCQ1SourceError("source raw gate failed")
    monotonic = diagnostics["source_residual_decreased_monotonically"]
    if type(monotonic) is not bool or monotonic is not True:
        raise HLT17SRCQ1SourceError("source residual monotonicity failed")
    for name in (
        "source_residual_infinity",
        "kinetic_condition_infinity",
        "reduction_constraint_infinity",
        "coordinate_speed_upper",
        "minimum_lapse",
        "minimum_radial_metric",
        "minimum_areal_radius_away_from_center",
    ):
        value = diagnostics[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise HLT17SRCQ1SourceError(f"{name} is not a finite scalar")
        number = float(value)
        if not math.isfinite(number):
            raise HLT17SRCQ1SourceError(f"{name} is not finite")
        diagnostics[name] = number
    if diagnostics["minimum_lapse"] <= 0.0:
        raise HLT17SRCQ1SourceError("nonpositive lapse")
    if diagnostics["minimum_radial_metric"] <= 0.0:
        raise HLT17SRCQ1SourceError("nonpositive radial metric")
    if diagnostics["minimum_areal_radius_away_from_center"] <= 0.0:
        raise HLT17SRCQ1SourceError("nonpositive areal radius away from center")
    encoded = {
        key: (
            value.hex()
            if type(value) is float
            else value
            if type(value) in (bool, int, str)
            else str(value)
        )
        for key, value in diagnostics.items()
        if type(key) is str
    }
    return {
        "diagnostics": encoded,
        "source_health": "recorded",
        "raw_health": "recorded",
        "kinetic_health": "recorded",
        "causal_sha256": before_causal,
        "constraint_health": "recorded",
        "input_state_sha256": before_state,
        "monitor_sha256": before_monitor,
        "tracer_sha256": before_tracer,
        "mutated_input": False,
        "mutated_monitor": False,
        "mutated_causal": False,
        "mutated_tracer": False,
    }


def eighteen_channel_decisions(assessment: Mapping[str, object]) -> list[dict[str, object]]:
    channels = assessment.get("channels")
    if type(channels) is not list or len(channels) != len(TDG6_COMPLETE_STATE_CHANNELS):
        raise HLT17SRCQ1Error("C1R1 assessment is not eighteen complete-state channels")
    records: list[dict[str, object]] = []
    for expected, item in zip(TDG6_COMPLETE_STATE_CHANNELS, channels, strict=True):
        channel = _mapping(item, "channel assessment")
        name = channel.get("channel")
        if name != expected:
            raise HLT17SRCQ1Error("eighteen-channel order differs")
        corrected = _mapping(channel.get("corrected"), "corrected path")
        decision = _mapping(corrected.get("decision"), "channel decision")
        records.append(
            {
                "channel": expected,
                "admission_passed": channel.get("admission_passed"),
                "classification": decision.get("classification"),
                "public_fine_debit": channel.get("public_fine_debit"),
                "gate_debit": channel.get("gate_debit"),
                "extra_enclosure_debit": channel.get("extra_enclosure_debit"),
            }
        )
    return records


def c1r1_shadow_record(prepared: object) -> dict[str, object]:
    bound = require_hlt17_prepared(prepared)
    identity = hlt17_implementation_identity()
    if bound.implementation_id != identity["implementation_id"]:
        raise HLT17SRCQ1IdentityError("prepared implementation is not the fixed C1R1 object")
    if bound.reference_wire_evaluator_id != identity["reference_wire_evaluator_id"]:
        raise HLT17SRCQ1IdentityError("prepared IMP1 reference wire differs")
    assessment = bound.assessment.as_mapping()
    record = {
        "prepared": True,
        "implementation_identity": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "c1r1_is_not_imp1_execution": True,
        "plan": encode_tdg7_plan(bound.plan),
        "eighteen_channel_decisions": eighteen_channel_decisions(assessment),
        "admission_passed": bound.assessment.admission_passed,
        "public_fine_debits": list(assessment["public_fine_debits"]),
        "receipt_sha256": bound.receipt_sha256,
        "assessment_sha256": bound.assessment.assessment_sha256,
        "family_sha256": bound.assessment.family_sha256,
        "corrected_endpoint_committable": False,
        "outer_or_medium_committable": False,
        "endpoint_adopted": False,
        "campaign_execution_authorized": False,
    }
    refuse_endpoint_payload(record)
    return record


def classify_prepare_outcome(result: object, prepared_record: Mapping[str, object] | None) -> str:
    disposition = getattr(result, "disposition", None)
    if disposition in {"source_retry_required", "source_retry_exhausted"}:
        return "source"
    if disposition in {"cfl_retry_required", "cfl_retry_exhausted"}:
        return "cfl"
    if disposition in {"temporal_retry_required", "temporal_retry_exhausted"}:
        return "temporal"
    if disposition == "invalid_premise":
        return "invalid"
    if disposition == "accepted_fine":
        raise HLT17SRCQ1Error("qualification adopted an endpoint")
    if disposition != "prepared_fresh":
        return "inconclusive"
    if not prepared_record:
        return "inconclusive"
    decisions = prepared_record.get("eighteen_channel_decisions")
    if type(decisions) is not list:
        return "inconclusive"
    classifications = [item.get("classification") for item in decisions]
    if any(item == "order_inconclusive" for item in classifications):
        return "inconclusive"
    if prepared_record.get("admission_passed") is True:
        return "prepared"
    return "temporal"


def _check_resources(
    *,
    member_key: str,
    ceilings: ResourceCeilings,
    wall_seconds: float,
    total_wall_seconds: float,
    rss: int,
) -> None:
    if rss > ceilings.max_rss_bytes:
        raise HLT17SRCQ1ResourceError(f"{member_key} exceeded RSS ceiling")
    limit = ceilings.max_wall_seconds_per_member[member_key]
    if wall_seconds > limit:
        raise HLT17SRCQ1ResourceError(f"{member_key} exceeded per-member wall ceiling")
    if total_wall_seconds > ceilings.max_total_wall_seconds:
        raise HLT17SRCQ1ResourceError("qualification exceeded total wall ceiling")


def _captured_member(origin: Pro20OriginCapture, member_key: str):
    matches = [item for item in origin.members if item.member_key == member_key]
    if len(matches) != 1:
        raise HLT17SRCQ1IdentityError(f"captured origin member differs: {member_key}")
    return matches[0]


def _require_expected_identity(
    actual: Mapping[str, object],
    expected: Mapping[str, object] | None,
    label: str,
) -> None:
    if expected is None:
        return
    wanted = _mapping(expected, label)
    for key, value in wanted.items():
        if actual.get(key) != value:
            raise HLT17SRCQ1IdentityError(f"foreign {label}: {key}")


def qualify_physical_source_origin(
    origin: object,
    *,
    static_input_bytes: Mapping[str, bytes],
    source_closure_sha256: str,
    environment: Mapping[str, str],
    resource_ceilings: ResourceCeilings | None = None,
    runtime_origin: object | None = None,
    runtime_origin_builder: Callable[..., object] | None = None,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
    shadow_preparer: Callable[..., object] | None = None,
    clock: Callable[[], float] | None = None,
    rss_bytes: Callable[[], int] | None = None,
    expected_origin_sha256: str | None = None,
    expected_environment: Mapping[str, str] | None = None,
    expected_implementation: Mapping[str, str] | None = None,
    expected_source_configurations: Mapping[str, str] | None = None,
    prospective_requested_cap_hex_by_member: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Qualify six restored members one at a time. No write, adopt, or retry."""

    if isinstance(origin, Pro20OriginCapture):
        revalidate_pro20_origin_capture(origin)
    elif not hasattr(origin, "capture_sha256") or not hasattr(origin, "members"):
        raise HLT17SRCQ1IdentityError("origin capture type differs")
    closure = _require_sha(source_closure_sha256, "source closure")
    env = _mapping(environment, "environment")
    _require_expected_identity(env, expected_environment, "environment")
    if expected_origin_sha256 is not None:
        if origin.capture_sha256 != _require_sha(expected_origin_sha256, "expected origin"):
            raise HLT17SRCQ1IdentityError("foreign origin identity")
    implementation = hlt17_implementation_identity()
    if expected_implementation is not None:
        require_hlt17_implementation_identity(
            expected_implementation, label="expected implementation"
        )
        if dict(expected_implementation) != implementation:
            raise HLT17SRCQ1IdentityError("foreign implementation identity")
    static_identities = static_input_identities(static_input_bytes)
    ceilings = recommended_resource_ceilings() if resource_ceilings is None else resource_ceilings
    if type(ceilings) is not ResourceCeilings:
        raise HLT17SRCQ1ResourceError("resource ceilings type differs")
    builder = runtime_origin_builder or build_pro20_runtime_origin
    now = clock or time.perf_counter
    rss = rss_bytes or current_rss_bytes
    started = now()
    initial_rss = int(rss())
    runtime = runtime_origin
    if runtime is None:
        try:
            runtime = builder(
                origin,
                static_input_bytes=static_input_bytes,
                source_closure_sha256=closure,
            )
        except Exception as error:
            raise HLT17SRCQ1SourceError(
                "physical static-factory/origin construction failed"
            ) from error
    members = getattr(runtime, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        raise HLT17SRCQ1IdentityError("runtime origin member order differs")
    if getattr(runtime, "captured_origin_sha256", None) != origin.capture_sha256:
        raise HLT17SRCQ1IdentityError("runtime origin capture identity differs")
    if getattr(runtime, "source_closure_sha256", None) != closure:
        raise HLT17SRCQ1IdentityError("runtime origin source-closure identity differs")
    requested_caps, requested_cap_facts = derive_prospective_requested_caps(
        origin,
        runtime,
        expected_cap_hex_by_member=prospective_requested_cap_hex_by_member,
    )
    evaluate = source_evaluator or (
        lambda view: evaluate_bound_source_once(
            binding=view.source_binding,
            time=float(view.member.time),
            state=view.member.state,
            transaction=view.member.transaction,
            tracers=view.member.tracers,
        )
    )
    prepare = shadow_preparer
    expected_configs = (
        None
        if expected_source_configurations is None
        else _mapping(expected_source_configurations, "expected source configurations")
    )
    member_records: list[dict[str, object]] = []
    halted: str | None = None
    for key in MEMBER_KEYS:
        view = members[key]
        if getattr(view, "member_key", None) != key:
            raise HLT17SRCQ1IdentityError(f"{key} runtime member identity differs")
        captured = _captured_member(origin, key)
        member = view.member
        binding = view.source_binding
        owned = check_method_owned_identity(key, member)
        check_inherited_counters(key, member, captured)
        cap = requested_caps[key]
        configuration = _require_sha(
            getattr(binding, "configuration_sha256", None),
            f"{key} source configuration",
        )
        if expected_configs is not None:
            expected_config = expected_configs.get(key)
            if expected_config != configuration:
                raise HLT17SRCQ1IdentityError(f"foreign source configuration: {key}")
        if configuration != getattr(view, "source_configuration_sha256", configuration):
            raise HLT17SRCQ1IdentityError(f"{key} source configuration identity differs")
        before = member_fingerprints(member, binding)
        member_started = now()
        rss_before = int(rss())
        try:
            _check_resources(
                member_key=key,
                ceilings=ceilings,
                wall_seconds=0.0,
                total_wall_seconds=now() - started,
                rss=rss_before,
            )
            source_record = dict(evaluate(view))
            checkpoint = view.checkpoint
            if hasattr(checkpoint, "agree"):
                checkpoint.agree()
            if prepare is not None:
                prepared_result = prepare(view, cap)
            else:
                prepared_result = checkpoint.prepare(requested_cap=cap)
            if getattr(prepared_result, "accepted_state_advanced", False) is True:
                raise HLT17SRCQ1Error("C1R1 shadow prepare advanced accepted state")
            prepared = getattr(prepared_result, "prepared", None)
            c1r1_record = None
            if prepared is not None:
                c1r1_record = c1r1_shadow_record(prepared)
            outcome = classify_prepare_outcome(prepared_result, c1r1_record)
            if outcome == "source" or outcome == "cfl":
                halted = outcome
            elif outcome == "invalid":
                halted = "invalid"
            elif outcome in {"temporal", "inconclusive"}:
                halted = outcome
            wall = now() - member_started
            rss_after = int(rss())
            _check_resources(
                member_key=key,
                ceilings=ceilings,
                wall_seconds=wall,
                total_wall_seconds=now() - started,
                rss=rss_after,
            )
            after = member_fingerprints(member, binding)
            if after != before:
                raise HLT17SRCQ1SourceError(
                    f"{key} fingerprints changed after no-adopt prepare"
                )
            member_records.append(
                {
                    "member_key": key,
                    "outcome": outcome,
                    "requested_cap_hex": cap.hex(),
                    "owned": owned,
                    "source_evaluation": source_record,
                    "c1r1": c1r1_record,
                    "prepare_disposition": getattr(prepared_result, "disposition", None),
                    "fingerprints_before": before,
                    "fingerprints_after": after,
                    "wall_seconds": wall,
                    "rss_bytes_before": rss_before,
                    "rss_bytes_after": rss_after,
                    "accepted_state_advanced": False,
                    "endpoint_adopted": False,
                }
            )
            if halted is not None:
                break
        except HLT17SRCQ1Error as error:
            wall = now() - member_started
            rss_after = int(rss())
            member_records.append(
                {
                    "member_key": key,
                    "outcome": error.outcome,
                    "reason": str(error),
                    "requested_cap_hex": cap.hex(),
                    "owned": owned,
                    "fingerprints_before": before,
                    "wall_seconds": wall,
                    "rss_bytes_before": rss_before,
                    "rss_bytes_after": rss_after,
                    "accepted_state_advanced": False,
                    "endpoint_adopted": False,
                }
            )
            halted = error.outcome
            break
        except (C1R1RefinementPathStop, C1R1RuntimeResourceStop, C1R1CoordinateLatticeStop) as error:
            if isinstance(error, C1R1RuntimeResourceStop):
                outcome = "resource"
            elif isinstance(error, C1R1CoordinateLatticeStop):
                outcome = "invalid"
            elif isinstance(error, C1R1RefinementPathStop) and isinstance(
                error.cause, CFLRetryRequired
            ):
                outcome = "cfl"
            else:
                outcome = "source"
            wall = now() - member_started
            member_records.append(
                {
                    "member_key": key,
                    "outcome": outcome,
                    "reason": str(error),
                    "requested_cap_hex": cap.hex(),
                    "owned": owned,
                    "fingerprints_before": before,
                    "wall_seconds": wall,
                    "rss_bytes_before": rss_before,
                    "rss_bytes_after": int(rss()),
                    "accepted_state_advanced": False,
                    "endpoint_adopted": False,
                }
            )
            halted = outcome
            break
    result = {
        "schema": SCHEMA,
        "artifact_id": ARTIFACT_ID,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "implementation": {
            **implementation_identity(),
            **implementation,
        },
        "environment": env,
        "input": {
            "static_input_sha256": static_identities,
            "source_closure_sha256": closure,
        },
        "origin": {
            "capture_sha256": origin.capture_sha256,
            "generation1_bridge_content_id": origin.generation1_bridge_content_id,
            "origin_time_hex": ORIGIN_TIME_HEX,
        },
        "source_configuration": {
            key: members[key].source_configuration_sha256 for key in MEMBER_KEYS
            if key in members and hasattr(members[key], "source_configuration_sha256")
        },
        "prospective_initial_requested_caps": requested_cap_facts,
        "resource_ceilings": ceilings.as_mapping(),
        "resource_facts": {
            "total_wall_seconds": now() - started,
            "rss_bytes_before_factory": initial_rss,
            "rss_bytes": int(rss()),
        },
        "members": member_records,
        "member_order": list(MEMBER_KEYS),
        "halted_outcome": halted,
        "completed_member_count": len(member_records),
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
        "store_published": False,
        "campaign_execution_authorized": False,
        "c1r1_actual_versus_imp1_reference_wire": {
            "c1r1_implementation_id": implementation["implementation_id"],
            "imp1_reference_wire_artifact_id": implementation["reference_wire_artifact_id"],
            "distinct": True,
        },
    }
    refuse_endpoint_payload(result)
    canonical_json_bytes(result)
    return result


def publish_qualification_result(
    output_directory: Path,
    result: Mapping[str, object],
    *,
    fault_hook: Callable[[str], None] | None = None,
) -> dict[str, object]:
    refuse_endpoint_payload(result)
    payload = canonical_json_bytes(dict(result))
    destination = Path(output_directory)
    if destination.exists():
        raise HLT17SRCQ1PublicationError("output namespace already exists")
    parent = destination.parent
    name = destination.name
    if not parent.exists() or not parent.is_dir():
        raise HLT17SRCQ1PublicationError("output parent directory is absent")
    try:
        receipt = publish_exclusive_directory(
            parent,
            name,
            {RESULT_LEAF: payload},
            _fault_hook=fault_hook,
        )
    except PostpublicationUncertainty as error:
        raise HLT17SRCQ1PublicationError(
            "qualification publication is uncertain", published=True
        ) from error
    except PrepublicationError as error:
        raise HLT17SRCQ1PublicationError("qualification publication failed before visibility") from error
    return {
        "relative_path": receipt.relative_path,
        "kind": receipt.kind,
        "payload_sha256": receipt.payload_sha256,
        "published": True,
    }


def qualify_and_publish(
    origin: object,
    *,
    repository_root: Path,
    output_directory: Path,
    authority_commit: str,
    static_input_bytes: Mapping[str, bytes],
    source_closure_sha256: str,
    environment: Mapping[str, str],
    authority_delta_validator: Callable[..., object] | None = None,
    resource_ceilings: ResourceCeilings | None = None,
    runtime_origin: object | None = None,
    runtime_origin_builder: Callable[..., object] | None = None,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
    shadow_preparer: Callable[..., object] | None = None,
    clock: Callable[[], float] | None = None,
    rss_bytes: Callable[[], int] | None = None,
    expected_origin_sha256: str | None = None,
    expected_environment: Mapping[str, str] | None = None,
    expected_implementation: Mapping[str, str] | None = None,
    expected_source_configurations: Mapping[str, str] | None = None,
    fault_hook: Callable[[str], None] | None = None,
) -> dict[str, object]:
    root = Path(repository_root)
    implementation = implementation_identity(root)
    digest = _require_sha(implementation.get("sha256"), "implementation sha256")
    authority = require_authority_delta(
        repository_root=root,
        authority_commit=authority_commit,
        implementation_sha256=digest,
        validator=authority_delta_validator,
    )
    requested_caps = authority.get("prospective_requested_cap_hex_by_member")
    if not isinstance(requested_caps, Mapping):
        raise HLT17SRCQ1AuthoritySeamError(
            "authority receipt omits prospective requested-cap mapping"
        )
    _require_requested_cap_mapping(requested_caps)
    expected_output = authority.get("output_namespace")
    if type(expected_output) is not str or not expected_output:
        raise HLT17SRCQ1AuthoritySeamError(
            "authority receipt omits the qualification output namespace"
        )
    output = Path(output_directory)
    expected_path = Path(expected_output)
    expected_path = expected_path if expected_path.is_absolute() else root / expected_path
    if not output.is_absolute() or output != expected_path:
        raise HLT17SRCQ1AuthoritySeamError(
            "qualification output directory differs from authority receipt"
        )
    try:
        result = qualify_physical_source_origin(
            origin,
            static_input_bytes=static_input_bytes,
            source_closure_sha256=source_closure_sha256,
            environment=environment,
            resource_ceilings=resource_ceilings,
            runtime_origin=runtime_origin,
            runtime_origin_builder=runtime_origin_builder,
            source_evaluator=source_evaluator,
            shadow_preparer=shadow_preparer,
            clock=clock,
            rss_bytes=rss_bytes,
            expected_origin_sha256=expected_origin_sha256,
            expected_environment=expected_environment,
            expected_implementation=expected_implementation,
            expected_source_configurations=expected_source_configurations,
            prospective_requested_cap_hex_by_member=requested_caps,
        )
    except HLT17SRCQ1Error as error:
        result = {
            "schema": SCHEMA,
            "artifact_id": ARTIFACT_ID,
            "protocol": PROTOCOL,
            "campaign_id": CAMPAIGN_ID,
            "halted_outcome": error.outcome,
            "reason": str(error),
            "implementation": dict(implementation),
            "environment": dict(environment),
            "input": {"source_closure_sha256": source_closure_sha256},
            "origin": {"capture_sha256": getattr(origin, "capture_sha256", None)},
            "resource_ceilings": (
                recommended_resource_ceilings()
                if resource_ceilings is None
                else resource_ceilings
            ).as_mapping(),
            "completed_member_count": 0,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
            "physical_source_qualification_attempted": True,
        }
    result["authority_commit"] = _require_commit(authority_commit, "authority_commit")
    result["authority_delta"] = authority
    result["implementation"]["sha256"] = digest
    publication = publish_qualification_result(
        Path(output_directory), result, fault_hook=fault_hook
    )
    result["publication"] = publication
    return result


__all__ = [
    "ARTIFACT_ID",
    "AUTHORITY_SEAM",
    "AUTHORITY_VALIDATOR_MODULE",
    "AUTHORITY_VALIDATOR_NAME",
    "HLT17SRCQ1AuthoritySeamError",
    "HLT17SRCQ1CapError",
    "HLT17SRCQ1Error",
    "HLT17SRCQ1IdentityError",
    "HLT17SRCQ1PublicationError",
    "HLT17SRCQ1ResourceError",
    "HLT17SRCQ1SourceError",
    "RESOURCE_CEILINGS_ARE_COORDINATOR_FREEZE_VALUES",
    "ResourceCeilings",
    "SCHEMA",
    "STATIC_INPUT_SHA256",
    "c1r1_shadow_record",
    "check_inherited_counters",
    "check_method_owned_identity",
    "current_rss_bytes",
    "evaluate_bound_source_once",
    "derive_prospective_requested_caps",
    "implementation_identity",
    "load_coordinator_authority_delta_validator",
    "member_fingerprints",
    "observe_environment",
    "plan_status",
    "publish_qualification_result",
    "qualify_and_publish",
    "qualify_physical_source_origin",
    "read_static_factory_inputs",
    "recommended_resource_ceilings",
    "refuse_endpoint_payload",
    "require_authority_delta",
    "static_input_identities",
]


RESOURCE_CEILINGS_ARE_COORDINATOR_FREEZE_VALUES = True
