"""Store-blind per-member seed handoff for prospective PRO20 RSRC1.

Six fresh children may each reconstruct one HLT15-origin member, but no store
is published until a future parent has independently reauthenticated a complete
six-member cohort.  Resource stops at this layer are pre-seed outcomes, not
campaign terminals.  This module opens no production namespace and writes no
state.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import time
from types import MappingProxyType
from typing import Callable, Mapping

from recursive_horizons.evidence_io import canonical_json_bytes, read_regular_file

from .hlt17_member_checkpoint import HLT17GenerationBundle
from .hlt17_member_codec import HLT17MemberEncoding
from .pro20_origin import (
    MEMBER_KEYS,
    ORIGIN_TIME_HEX,
    capture_pro20_historical_origin,
)
from .pro20_rsrc1_isolation import (
    CHILD_MAX_PEAK_RSS_BYTES,
    ChildMetrics,
    CompletedChild,
    MAX_CHILD_WALL_SECONDS,
    MAX_HANDOFF_BYTES,
    MAX_STDERR_BYTES,
    PRO20RSRC1ForensicUncertainty,
    decode_bundle,
    encode_bundle,
    parse_canonical_handoff,
)
from .pro20_rsrc1_member import CAMPAIGN_ID, EVENT_TARGET, build_rsrc1_member
from .proto19_gr0_static_factory import STATIC_INPUT_PATHS


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
SEED_REQUEST_SCHEMA = "FGC-1-PRO20-EV1-RSRC1-seed-request-v1"
SEED_RESPONSE_SCHEMA = "FGC-1-PRO20-EV1-RSRC1-seed-response-v1"


class PRO20RSRC1SeedError(ValueError):
    """The per-member seed request, response, or cohort differs."""


def _sha(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise PRO20RSRC1SeedError("seed payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if type(value) is not str or len(value) != 64 or value.lower() != value:
        raise PRO20RSRC1SeedError(f"{label} is not lowercase SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise PRO20RSRC1SeedError(f"{label} is not lowercase SHA-256") from error
    return value


def _member_key(value: object) -> str:
    if type(value) is not str or value not in MEMBER_KEYS:
        raise PRO20RSRC1SeedError("seed member is outside the live cohort")
    return value


def _exact(value: object, fields: frozenset[str], label: str) -> dict[str, object]:
    if type(value) is not dict or set(value) != fields:
        raise PRO20RSRC1SeedError(f"{label} fields differ")
    return dict(value)


@dataclass(frozen=True, slots=True)
class IsolatedSeedRequest:
    member_key: str
    authority_receipt_sha256: str
    origin_capture_sha256: str
    generation1_bridge_content_id: str
    source_closure_sha256: str
    source_configuration_sha256: str

    def __post_init__(self) -> None:
        _member_key(self.member_key)
        for value, label in (
            (self.authority_receipt_sha256, "authority receipt"),
            (self.origin_capture_sha256, "origin capture"),
            (self.generation1_bridge_content_id, "generation-one bridge"),
            (self.source_closure_sha256, "source closure"),
            (self.source_configuration_sha256, "source configuration"),
        ):
            _require_sha(value, label)

    @property
    def mapping_without_id(self) -> dict[str, object]:
        return {
            "authority_receipt_sha256": self.authority_receipt_sha256,
            "generation1_bridge_content_id": self.generation1_bridge_content_id,
            "member_key": self.member_key,
            "origin_capture_sha256": self.origin_capture_sha256,
            "schema": SEED_REQUEST_SCHEMA,
            "source_closure_sha256": self.source_closure_sha256,
            "source_configuration_sha256": self.source_configuration_sha256,
        }

    @property
    def request_sha256(self) -> str:
        return _sha(canonical_json_bytes(self.mapping_without_id))

    @property
    def mapping(self) -> dict[str, object]:
        return {**self.mapping_without_id, "request_sha256": self.request_sha256}


_SEED_REQUEST_FIELDS = frozenset(
    {
        "authority_receipt_sha256",
        "generation1_bridge_content_id",
        "member_key",
        "origin_capture_sha256",
        "request_sha256",
        "schema",
        "source_closure_sha256",
        "source_configuration_sha256",
    }
)


def encode_seed_request(request: IsolatedSeedRequest) -> bytes:
    request.__post_init__()
    return canonical_json_bytes(request.mapping)


def decode_seed_request(raw: bytes) -> IsolatedSeedRequest:
    mapping = _exact(
        parse_canonical_handoff(raw, "seed request"),
        _SEED_REQUEST_FIELDS,
        "seed request",
    )
    if mapping["schema"] != SEED_REQUEST_SCHEMA:
        raise PRO20RSRC1SeedError("seed request schema differs")
    request = IsolatedSeedRequest(
        member_key=_member_key(mapping["member_key"]),
        authority_receipt_sha256=_require_sha(
            mapping["authority_receipt_sha256"], "authority receipt"
        ),
        origin_capture_sha256=_require_sha(
            mapping["origin_capture_sha256"], "origin capture"
        ),
        generation1_bridge_content_id=_require_sha(
            mapping["generation1_bridge_content_id"], "generation-one bridge"
        ),
        source_closure_sha256=_require_sha(
            mapping["source_closure_sha256"], "source closure"
        ),
        source_configuration_sha256=_require_sha(
            mapping["source_configuration_sha256"], "source configuration"
        ),
    )
    if mapping["request_sha256"] != request.request_sha256:
        raise PRO20RSRC1SeedError("seed request digest differs")
    return request


@dataclass(frozen=True, slots=True)
class SeedOutcome:
    member_key: str
    request_sha256: str
    bundle: HLT17GenerationBundle | None
    construction_sha256: str | None
    resource_evidence: Mapping[str, object] | None
    metrics: ChildMetrics

    def __post_init__(self) -> None:
        _member_key(self.member_key)
        _require_sha(self.request_sha256, "seed request")
        if not isinstance(self.metrics, ChildMetrics):
            raise TypeError("seed metrics must be ChildMetrics")
        successful = self.bundle is not None
        if successful:
            self.bundle.__post_init__()
            _require_sha(self.construction_sha256, "seed construction")
            if self.resource_evidence is not None:
                raise PRO20RSRC1SeedError("successful seed carries resource evidence")
        else:
            if self.construction_sha256 is not None:
                raise PRO20RSRC1SeedError("resource seed carries construction identity")
            if type(self.resource_evidence) is not dict:
                raise PRO20RSRC1SeedError("resource seed omitted evidence")
            if set(self.resource_evidence) != {"limit_bytes", "reason", "scope"}:
                raise PRO20RSRC1SeedError("resource seed evidence fields differ")
            if self.resource_evidence["limit_bytes"] != CHILD_MAX_PEAK_RSS_BYTES:
                raise PRO20RSRC1SeedError("resource seed limit differs")
        object.__setattr__(
            self,
            "resource_evidence",
            None if self.resource_evidence is None else dict(self.resource_evidence),
        )

    @property
    def successful(self) -> bool:
        return self.bundle is not None


_SEED_RESPONSE_FIELDS = frozenset(
    {
        "bundle",
        "child_peak_rss_bytes",
        "child_wall_seconds_hex",
        "construction_sha256",
        "member_key",
        "outcome",
        "request_sha256",
        "resource_evidence",
        "schema",
    }
)


def encode_seed_response(
    request: IsolatedSeedRequest,
    *,
    bundle: HLT17GenerationBundle,
    construction_sha256: str,
    peak_rss_bytes: int,
    wall_seconds: float,
) -> bytes:
    request.__post_init__()
    bundle.__post_init__()
    construction = _require_sha(construction_sha256, "seed construction")
    metrics = ChildMetrics(int(peak_rss_bytes), float(wall_seconds))
    over = metrics.peak_rss_bytes > CHILD_MAX_PEAK_RSS_BYTES
    response = {
        "bundle": None if over else encode_bundle(bundle),
        "child_peak_rss_bytes": metrics.peak_rss_bytes,
        "child_wall_seconds_hex": metrics.wall_seconds.hex(),
        "construction_sha256": None if over else construction,
        "member_key": request.member_key,
        "outcome": "resource_stop" if over else "seed_bundle",
        "request_sha256": request.request_sha256,
        "resource_evidence": (
            {
                "limit_bytes": CHILD_MAX_PEAK_RSS_BYTES,
                "reason": "isolated seed child peak RSS exceeded",
                "scope": "one_seed_member",
            }
            if over
            else None
        ),
        "schema": SEED_RESPONSE_SCHEMA,
    }
    return canonical_json_bytes(response)


def encode_seed_memory_stop(
    request: IsolatedSeedRequest,
    *,
    peak_rss_bytes: int,
) -> bytes:
    request.__post_init__()
    metrics = ChildMetrics(int(peak_rss_bytes), 0.0)
    return canonical_json_bytes(
        {
            "bundle": None,
            "child_peak_rss_bytes": metrics.peak_rss_bytes,
            "child_wall_seconds_hex": 0.0.hex(),
            "construction_sha256": None,
            "member_key": request.member_key,
            "outcome": "resource_stop",
            "request_sha256": request.request_sha256,
            "resource_evidence": {
                "limit_bytes": CHILD_MAX_PEAK_RSS_BYTES,
                "reason": "isolated seed child memory allocation failed",
                "scope": "one_seed_member",
            },
            "schema": SEED_RESPONSE_SCHEMA,
        }
    )


def decode_seed_response(
    request: IsolatedSeedRequest,
    raw: bytes,
) -> SeedOutcome:
    request.__post_init__()
    mapping = _exact(
        parse_canonical_handoff(raw, "seed response"),
        _SEED_RESPONSE_FIELDS,
        "seed response",
    )
    if mapping["schema"] != SEED_RESPONSE_SCHEMA:
        raise PRO20RSRC1SeedError("seed response schema differs")
    if mapping["request_sha256"] != request.request_sha256:
        raise PRO20RSRC1SeedError("seed response request digest differs")
    if mapping["member_key"] != request.member_key:
        raise PRO20RSRC1SeedError("seed response member differs")
    if type(mapping["child_peak_rss_bytes"]) is not int:
        raise PRO20RSRC1SeedError("seed response peak RSS differs")
    if type(mapping["child_wall_seconds_hex"]) is not str:
        raise PRO20RSRC1SeedError("seed response wall differs")
    try:
        wall = float.fromhex(mapping["child_wall_seconds_hex"])
    except ValueError as error:
        raise PRO20RSRC1SeedError("seed response wall differs") from error
    metrics = ChildMetrics(mapping["child_peak_rss_bytes"], wall)
    if mapping["outcome"] == "seed_bundle":
        if metrics.peak_rss_bytes > CHILD_MAX_PEAK_RSS_BYTES:
            raise PRO20RSRC1SeedError("over-ceiling seed response retained a bundle")
        if mapping["resource_evidence"] is not None:
            raise PRO20RSRC1SeedError("successful seed response carries resource evidence")
        outcome = SeedOutcome(
            member_key=request.member_key,
            request_sha256=request.request_sha256,
            bundle=decode_bundle(mapping["bundle"]),
            construction_sha256=_require_sha(
                mapping["construction_sha256"], "seed construction"
            ),
            resource_evidence=None,
            metrics=metrics,
        )
    elif mapping["outcome"] == "resource_stop":
        if mapping["bundle"] is not None or mapping["construction_sha256"] is not None:
            raise PRO20RSRC1SeedError("resource seed response retained an endpoint")
        outcome = SeedOutcome(
            member_key=request.member_key,
            request_sha256=request.request_sha256,
            bundle=None,
            construction_sha256=None,
            resource_evidence=mapping["resource_evidence"],
            metrics=metrics,
        )
    else:
        raise PRO20RSRC1SeedError("seed response outcome differs")
    return outcome


def reduce_completed_seed_child(
    request: IsolatedSeedRequest,
    completed: CompletedChild,
) -> SeedOutcome:
    """Accept only one silent, bounded, zero-exit canonical seed response."""

    request.__post_init__()
    completed.__post_init__()
    if completed.returncode != 0:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated seed child exited without a valid response"
        )
    if completed.stderr:
        if len(completed.stderr) > MAX_STDERR_BYTES:
            raise PRO20RSRC1ForensicUncertainty(
                "isolated seed child stderr was unbounded"
            )
        raise PRO20RSRC1ForensicUncertainty("isolated seed child emitted stderr")
    if len(completed.stdout) > MAX_HANDOFF_BYTES:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated seed child response was unbounded"
        )
    if completed.wall_seconds > MAX_CHILD_WALL_SECONDS:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated seed child exceeded its wall allowance"
        )
    try:
        return decode_seed_response(request, completed.stdout)
    except (PRO20RSRC1SeedError, TypeError, ValueError) as error:
        raise PRO20RSRC1ForensicUncertainty(
            "isolated seed child response failed parent reauthentication"
        ) from error


def _load_static_inputs(root: Path) -> dict[str, bytes]:
    return {relative: read_regular_file(root, relative) for relative in STATIC_INPUT_PATHS}


def execute_seed_child_request(
    repository_root: Path,
    request_raw: bytes,
    *,
    origin_loader: Callable[[Path], object] = capture_pro20_historical_origin,
    static_loader: Callable[[Path], Mapping[str, bytes]] = _load_static_inputs,
    member_builder: Callable[..., object] = build_rsrc1_member,
    peak_rss: Callable[[], int],
    now: Callable[[], float] = time.monotonic,
) -> bytes:
    """Construct one physical seed in memory and return bytes; write nothing."""

    root = Path(repository_root)
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise PRO20RSRC1SeedError("seed child repository root is not canonical")
    request = decode_seed_request(request_raw)
    started = float(now())
    origin = origin_loader(root)
    record = member_builder(
        origin,
        member_key=request.member_key,
        static_input_bytes=static_loader(root),
        source_closure_sha256=request.source_closure_sha256,
    )
    for observed, expected, label in (
        (getattr(record, "captured_origin_sha256", None), request.origin_capture_sha256, "origin"),
        (
            getattr(record, "generation1_bridge_content_id", None),
            request.generation1_bridge_content_id,
            "generation-one bridge",
        ),
        (getattr(record, "source_closure_sha256", None), request.source_closure_sha256, "source closure"),
        (
            getattr(record, "source_configuration_sha256", None),
            request.source_configuration_sha256,
            "source configuration",
        ),
    ):
        if observed != expected:
            raise PRO20RSRC1SeedError(f"seed child {label} identity differs")
    bundle = getattr(record, "bundle", None)
    construction = getattr(record, "construction_sha256", None)
    if not isinstance(bundle, HLT17GenerationBundle):
        raise PRO20RSRC1SeedError("seed child omitted its generation bundle")
    return encode_seed_response(
        request,
        bundle=bundle,
        construction_sha256=construction,
        peak_rss_bytes=peak_rss(),
        wall_seconds=float(now()) - started,
    )


@dataclass(frozen=True, slots=True)
class SeedCohort:
    bundles: Mapping[str, HLT17GenerationBundle]
    construction_sha256_by_member: Mapping[str, str]
    cohort_sha256: str

    def __post_init__(self) -> None:
        bundles = MappingProxyType(dict(self.bundles))
        constructions = MappingProxyType(dict(self.construction_sha256_by_member))
        object.__setattr__(self, "bundles", bundles)
        object.__setattr__(self, "construction_sha256_by_member", constructions)
        if tuple(bundles) != MEMBER_KEYS or tuple(constructions) != MEMBER_KEYS:
            raise PRO20RSRC1SeedError("seed cohort order differs")
        _require_sha(self.cohort_sha256, "seed cohort")


def bind_seed_cohort(outcomes: Mapping[str, SeedOutcome]) -> SeedCohort:
    """Require six fresh successful children before any future store publish."""

    if not isinstance(outcomes, Mapping) or tuple(outcomes) != MEMBER_KEYS:
        raise PRO20RSRC1SeedError("seed outcome cohort order differs")
    bundles: dict[str, HLT17GenerationBundle] = {}
    constructions: dict[str, str] = {}
    digest_rows: dict[str, object] = {}
    for key in MEMBER_KEYS:
        outcome = outcomes[key]
        if not isinstance(outcome, SeedOutcome) or outcome.member_key != key:
            raise PRO20RSRC1SeedError("seed outcome member differs")
        outcome.__post_init__()
        if not outcome.successful or outcome.bundle is None:
            raise PRO20RSRC1SeedError("seed cohort contains a resource stop")
        encoding = HLT17MemberEncoding(
            outcome.bundle.descriptor,
            outcome.bundle.payload,
        )
        identity = encoding.metadata["runtime_identity"]
        if (
            identity.get("campaign_id") != CAMPAIGN_ID
            or identity.get("member_key") != key
            or identity.get("accepted_generation") != 0
            or encoding.metadata.get("accepted_time_hex") != ORIGIN_TIME_HEX
            or outcome.bundle.cursor_mapping.get("event_target_hex")
            != EVENT_TARGET.hex()
        ):
            raise PRO20RSRC1SeedError("seed bundle identity differs")
        construction = _require_sha(outcome.construction_sha256, "seed construction")
        bundles[key] = outcome.bundle
        constructions[key] = construction
        digest_rows[key] = {
            "construction_sha256": construction,
            "cursor_sha256": outcome.bundle.cursor_sha256,
            "descriptor_sha256": outcome.bundle.descriptor_sha256,
            "payload_sha256": outcome.bundle.payload_sha256,
            "physical_state_sha256": outcome.bundle.physical_state_sha256,
        }
    cohort_sha256 = _sha(
        canonical_json_bytes(
            {
                "artifact_id": ARTIFACT_ID,
                "campaign_id": CAMPAIGN_ID,
                "members": digest_rows,
                "schema": "FGC-1-PRO20-EV1-RSRC1-seed-cohort-v1",
            }
        )
    )
    return SeedCohort(bundles, constructions, cohort_sha256)


__all__ = [
    "ARTIFACT_ID",
    "IsolatedSeedRequest",
    "PRO20RSRC1SeedError",
    "SEED_REQUEST_SCHEMA",
    "SEED_RESPONSE_SCHEMA",
    "SeedCohort",
    "SeedOutcome",
    "bind_seed_cohort",
    "decode_seed_request",
    "decode_seed_response",
    "encode_seed_memory_stop",
    "encode_seed_request",
    "encode_seed_response",
    "execute_seed_child_request",
    "reduce_completed_seed_child",
]
