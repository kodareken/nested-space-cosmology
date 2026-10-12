"""Parent-only prospective RSRC1 first-event coordinator.

The parent owns authority, the fresh durable store, scheduling and publication.
It never constructs a live PDE member.  Each seed and scheduled attempt is a
fresh child-process handoff that the parent reauthenticates before publication.
No authority module or live Make target exists in this implementation slice.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Callable, Mapping

from recursive_horizons.evidence_io import read_regular_file

from .pro20_origin import MEMBER_KEYS
from .pro20_rsrc1_attempt import classify_attempt
from .pro20_rsrc1_isolation import (
    CompletedChild,
    IsolatedAttemptRequest,
    PRO20RSRC1ForensicUncertainty,
    check_parent_resource_premises,
    encode_attempt_request,
    launch_child_process,
    reduce_completed_child,
)
from .pro20_rsrc1_member import CAMPAIGN_ID
from .pro20_rsrc1_protocol import (
    PRODUCTION_NAMESPACE,
    PRO20RSRC1AuthorityReceipt,
    PRO20RSRC1ValidatedGeneration,
    build_authority_receipt,
    next_required_member_key,
)
from .pro20_rsrc1_seed import (
    IsolatedSeedRequest,
    SeedOutcome,
    bind_seed_cohort,
    encode_seed_request,
    reduce_completed_seed_child,
)
from .pro20_rsrc1_store import (
    PRODUCTION_CONTAINER,
    PRO20RSRC1CampaignStore,
    PRO20RSRC1CampaignStoreError,
    PRO20RSRC1PostpublicationUncertainty,
    PRO20RSRC1TerminalStore,
    PRO20RSRC1UnclosedSession,
    PRO20RSRC1WriterCapability,
)


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
FUTURE_FREEZE_ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1-FRZ1"
SCHEMA = "FGC-1-PRO20-EV1-RSRC1-parent-runtime-v1"
AUTHORITY_MODULE = "recursive_horizons.fgc.evolution.pro20_rsrc1_authority"
AUTHORITY_FUNCTION = "validate_authority_delta"
CHILD_SCRIPT = "scripts/run_fgc_pro20_rsrc1_child.py"
OWNER_TOKEN = "pro20-rsrc1-first-event-one-shot"
MAX_TOTAL_WALL_SECONDS = 86400.0
MAX_NAMESPACE_BYTES = 17179869184
MIN_FREE_DISK_BYTES = 34359738368
MAX_PUBLISHED_ATTEMPTS = 4096
PLANNING_ACCEPTED_STEP_ESTIMATE = 224
IMPLEMENTATION_PATHS = (
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_authority.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_protocol.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_store.py",
    CHILD_SCRIPT,
    "scripts/run_fgc_pro20_rsrc1.py",
)
_VM_PAGE = re.compile(r"page size of (?P<size>[0-9]+) bytes")
_VM_ROW = re.compile(r"^(?P<name>Pages [A-Za-z ]+):\s+(?P<count>[0-9]+)\.\s*$")


class PRO20RSRC1RuntimeError(ValueError):
    """Typed prospective coordinator stop; never a scientific classification."""

    phase = "runtime"
    outcome = "invalid"
    published = False
    abandon_writer = False
    physical_source_constructed = False


class PRO20RSRC1AuthorityError(PRO20RSRC1RuntimeError):
    phase = "authority"


class PRO20RSRC1PremiseError(PRO20RSRC1RuntimeError):
    phase = "premise"


class PRO20RSRC1ResourceError(PRO20RSRC1RuntimeError):
    phase = "resource"
    outcome = "resource"


class PRO20RSRC1SeedResourceError(PRO20RSRC1ResourceError):
    """A pre-store seed resource stop; no namespace or terminal exists."""


class PRO20RSRC1ForensicError(PRO20RSRC1RuntimeError):
    phase = "forensic"
    outcome = "inconclusive"
    abandon_writer = True


class PRO20RSRC1PostpublicationError(PRO20RSRC1ForensicError):
    phase = "postpublication"
    published = True


def _fail(message: str, cls: type[PRO20RSRC1RuntimeError]) -> None:
    raise cls(message)


def _canonical_root(root: Path) -> Path:
    if not isinstance(root, Path) or not root.is_absolute():
        _fail("repository root is not an absolute pathlib path", PRO20RSRC1PremiseError)
    try:
        resolved = root.resolve(strict=True)
    except OSError as error:
        raise PRO20RSRC1PremiseError("repository root cannot be resolved") from error
    if resolved != root:
        _fail("repository root traverses a symlink", PRO20RSRC1PremiseError)
    return root


def _sha(raw: bytes) -> str:
    if type(raw) is not bytes:
        _fail("implementation payload is not immutable bytes", PRO20RSRC1PremiseError)
    return sha256(raw).hexdigest()


def implementation_identity(repository_root: Path) -> dict[str, object]:
    root = _canonical_root(repository_root)
    hashes = {
        relative: _sha(read_regular_file(root, relative))
        for relative in IMPLEMENTATION_PATHS
    }
    digest = _sha(
        "".join(f"{path}\0{hashes[path]}\n" for path in IMPLEMENTATION_PATHS).encode(
            "ascii"
        )
    )
    return {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "paths": list(IMPLEMENTATION_PATHS),
        "sha256_by_path": hashes,
        "implementation_sha256": digest,
        "campaign_execution_authorized": False,
    }


def current_parent_rss_bytes() -> int:
    try:
        completed = subprocess.run(
            ("ps", "-o", "rss=", "-p", str(os.getpid())),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=5.0,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PRO20RSRC1ResourceError("parent current RSS cannot be observed") from error
    text = completed.stdout.strip()
    if completed.returncode != 0 or not text.isdecimal():
        _fail("parent current RSS cannot be observed", PRO20RSRC1ResourceError)
    return int(text) * 1024


def available_memory_bytes() -> int:
    if sys.platform != "darwin":
        _fail("available-memory observer is not the frozen Darwin route", PRO20RSRC1ResourceError)
    try:
        completed = subprocess.run(
            ("vm_stat",),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=5.0,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PRO20RSRC1ResourceError("available memory cannot be observed") from error
    lines = completed.stdout.splitlines()
    if completed.returncode != 0 or not lines:
        _fail("available memory cannot be observed", PRO20RSRC1ResourceError)
    page_match = _VM_PAGE.search(lines[0])
    if page_match is None:
        _fail("vm_stat page size is malformed", PRO20RSRC1ResourceError)
    page_size = int(page_match.group("size"))
    counts: dict[str, int] = {}
    for line in lines[1:]:
        matched = _VM_ROW.match(line)
        if matched is not None:
            counts[matched.group("name")] = int(matched.group("count"))
    names = (
        "Pages free",
        "Pages inactive",
        "Pages speculative",
        "Pages purgeable",
    )
    if any(name not in counts for name in names):
        _fail("vm_stat available-page inventory differs", PRO20RSRC1ResourceError)
    return page_size * sum(counts[name] for name in names)


def free_disk_bytes(path: Path) -> int:
    try:
        info = os.statvfs(os.fspath(path))
    except OSError as error:
        raise PRO20RSRC1ResourceError("free disk cannot be observed") from error
    answer = int(info.f_bavail) * int(info.f_frsize)
    if answer < 0:
        _fail("free disk observation is negative", PRO20RSRC1ResourceError)
    return answer


def namespace_byte_count(root: Path) -> int:
    target = root.joinpath(*PRODUCTION_CONTAINER.split("/"))
    try:
        target.lstat()
    except FileNotFoundError:
        return 0
    except OSError as error:
        raise PRO20RSRC1ResourceError("namespace cannot be observed") from error
    if not target.is_dir() or target.is_symlink():
        _fail("namespace is not a safe directory", PRO20RSRC1ResourceError)
    total = 0
    for path in (target, *target.rglob("*")):
        metadata = path.lstat()
        if path.is_symlink():
            _fail("namespace contains a symlink", PRO20RSRC1ResourceError)
        if path.is_file():
            if metadata.st_nlink != 1:
                _fail("namespace contains a hard link", PRO20RSRC1ResourceError)
            total += int(metadata.st_size)
        elif not path.is_dir():
            _fail("namespace contains a nonregular node", PRO20RSRC1ResourceError)
    return total


def _require_namespace_absent(root: Path) -> None:
    target = root.joinpath(*PRODUCTION_CONTAINER.split("/"))
    try:
        target.lstat()
    except FileNotFoundError:
        parent = target.parent
        if not parent.exists():
            return
        names = os.listdir(parent)
        if any(name.startswith(f".{target.name}.") for name in names):
            _fail("RSRC1 staging namespace already exists", PRO20RSRC1PremiseError)
        return
    except OSError as error:
        raise PRO20RSRC1PremiseError("RSRC1 namespace cannot be inspected") from error
    _fail("RSRC1 namespace already exists; one-shot is consumed", PRO20RSRC1PremiseError)


def _resource_check(
    *,
    root: Path,
    started: float,
    published_attempts: int,
    now: Callable[[], float],
    parent_rss: Callable[[], int],
    available_memory: Callable[[], int],
    disk_free: Callable[[Path], int],
    namespace_bytes: Callable[[Path], int],
) -> None:
    if float(now()) - float(started) > MAX_TOTAL_WALL_SECONDS:
        _fail("max total wall exceeded", PRO20RSRC1ResourceError)
    try:
        check_parent_resource_premises(
            parent_current_rss_bytes=parent_rss(),
            available_memory_bytes=available_memory(),
        )
    except ValueError as error:
        raise PRO20RSRC1ResourceError(str(error)) from error
    if disk_free(root) < MIN_FREE_DISK_BYTES:
        _fail("free disk is below the RSRC1 minimum", PRO20RSRC1ResourceError)
    if namespace_bytes(root) > MAX_NAMESPACE_BYTES:
        _fail("RSRC1 namespace byte ceiling exceeded", PRO20RSRC1ResourceError)
    if published_attempts > MAX_PUBLISHED_ATTEMPTS:
        _fail("published-attempt ceiling exceeded", PRO20RSRC1ResourceError)


def _load_authority_validator() -> Callable[..., Mapping[str, object]]:
    from importlib import import_module

    try:
        module = import_module(AUTHORITY_MODULE)
        validator = getattr(module, AUTHORITY_FUNCTION)
    except (ImportError, AttributeError) as error:
        raise PRO20RSRC1AuthorityError(
            "RSRC1 exact-delta authority is not implemented"
        ) from error
    if not callable(validator):
        _fail("RSRC1 authority validator is not callable", PRO20RSRC1AuthorityError)
    return validator


def _authority(
    *,
    root: Path,
    authority_commit: str | None,
    implementation_sha256: str,
    validator: Callable[..., Mapping[str, object]] | None,
) -> dict[str, object]:
    if authority_commit is None:
        _fail("RSRC1 run refuses absent authority", PRO20RSRC1AuthorityError)
    loaded = validator or _load_authority_validator()
    try:
        value = loaded(
            repository_root=root,
            authority_commit=authority_commit,
            implementation_sha256=implementation_sha256,
        )
    except PRO20RSRC1RuntimeError:
        raise
    except Exception as error:
        raise PRO20RSRC1AuthorityError("RSRC1 authority validation failed") from error
    if type(value) is bool or not isinstance(value, Mapping):
        _fail("RSRC1 authority receipt is not a mapping", PRO20RSRC1AuthorityError)
    mapping = dict(value)
    identities = mapping.get("receipt_identities")
    source_config = mapping.get("source_configuration")
    historical = mapping.get("historical_source_tree")
    if not isinstance(identities, Mapping):
        _fail("RSRC1 authority identities are absent", PRO20RSRC1AuthorityError)
    if not isinstance(source_config, Mapping) or tuple(source_config) != MEMBER_KEYS:
        _fail("RSRC1 source-configuration inventory differs", PRO20RSRC1AuthorityError)
    if not isinstance(historical, Mapping):
        _fail("RSRC1 historical source identity is absent", PRO20RSRC1AuthorityError)
    if mapping.get("campaign_id") != CAMPAIGN_ID:
        _fail("RSRC1 campaign identity differs", PRO20RSRC1AuthorityError)
    if mapping.get("output_namespace") != PRODUCTION_NAMESPACE:
        _fail("RSRC1 authority namespace differs", PRO20RSRC1AuthorityError)
    if identities.get("implementation_sha256") != implementation_sha256:
        _fail("RSRC1 authority implementation identity differs", PRO20RSRC1AuthorityError)
    return mapping


def require_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
    validator: Callable[..., Mapping[str, object]] | None = None,
) -> dict[str, object]:
    """Public read-only authority seam used by the guarded status command."""

    return _authority(
        root=_canonical_root(Path(repository_root)),
        authority_commit=authority_commit,
        implementation_sha256=implementation_sha256,
        validator=validator,
    )


def _receipt(authority: Mapping[str, object]) -> PRO20RSRC1AuthorityReceipt:
    identities = authority["receipt_identities"]
    return build_authority_receipt(
        campaign_id=CAMPAIGN_ID,
        authority_sha256=str(identities["authority_sha256"]),
        implementation_sha256=str(identities["implementation_sha256"]),
        config_sha256=str(identities["config_sha256"]),
        source_sha256=str(identities["source_sha256"]),
        origin_sha256=str(identities["origin_sha256"]),
        environment_sha256=str(identities["environment_sha256"]),
    )


def _seed_request(
    member_key: str,
    *,
    receipt: PRO20RSRC1AuthorityReceipt,
    authority: Mapping[str, object],
) -> IsolatedSeedRequest:
    historical = authority["historical_source_tree"]
    return IsolatedSeedRequest(
        member_key=member_key,
        authority_receipt_sha256=receipt.sha256,
        origin_capture_sha256=str(receipt.mapping["origin_sha256"]),
        generation1_bridge_content_id=str(historical["generation1_bridge_content_id"]),
        source_closure_sha256=str(receipt.mapping["source_sha256"]),
        source_configuration_sha256=str(authority["source_configuration"][member_key]),
    )


def _attempt_request(
    member_key: str,
    predecessor,
    *,
    receipt: PRO20RSRC1AuthorityReceipt,
    authority: Mapping[str, object],
) -> IsolatedAttemptRequest:
    historical = authority["historical_source_tree"]
    return IsolatedAttemptRequest(
        member_key=member_key,
        predecessor_bundle=predecessor,
        authority_receipt_sha256=receipt.sha256,
        origin_capture_sha256=str(receipt.mapping["origin_sha256"]),
        generation1_bridge_content_id=str(historical["generation1_bridge_content_id"]),
        source_closure_sha256=str(receipt.mapping["source_sha256"]),
        source_configuration_sha256=str(authority["source_configuration"][member_key]),
    )


def _child_command(root: Path) -> tuple[str, ...]:
    script = root / CHILD_SCRIPT
    if not script.is_file():
        _fail("RSRC1 child script is absent", PRO20RSRC1PremiseError)
    return (sys.executable, "-I", "-B", os.fspath(script))


def _launch_seed(root: Path, request: IsolatedSeedRequest) -> CompletedChild:
    return launch_child_process(
        _child_command(root),
        repository_root=root,
        request_raw=encode_seed_request(request),
    )


def _launch_attempt(root: Path, request: IsolatedAttemptRequest) -> CompletedChild:
    return launch_child_process(
        _child_command(root),
        repository_root=root,
        request_raw=encode_attempt_request(request),
    )


@dataclass
class RuntimeHooks:
    """Injectable authority, child and resource seams for bounded tests."""

    authority_validator: Callable[..., Mapping[str, object]] | None = None
    seed_child: Callable[[Path, IsolatedSeedRequest], CompletedChild] = _launch_seed
    attempt_child: Callable[[Path, IsolatedAttemptRequest], CompletedChild] = _launch_attempt
    store_cls: type[PRO20RSRC1CampaignStore] = PRO20RSRC1CampaignStore
    now: Callable[[], float] = time.monotonic
    parent_rss: Callable[[], int] = current_parent_rss_bytes
    available_memory: Callable[[], int] = available_memory_bytes
    disk_free: Callable[[Path], int] = free_disk_bytes
    namespace_bytes: Callable[[Path], int] = namespace_byte_count
    implementation_identity_fn: Callable[[Path], Mapping[str, object]] = (
        implementation_identity
    )
    seed_requests: list[str] = field(default_factory=list)
    attempt_requests: list[str] = field(default_factory=list)


def _publish_resource_terminal(
    store: PRO20RSRC1CampaignStore,
    capability: PRO20RSRC1WriterCapability,
    last: PRO20RSRC1ValidatedGeneration,
    reason: str,
) -> PRO20RSRC1ValidatedGeneration:
    member_key = next_required_member_key(last)
    classification = classify_attempt(
        member_key=member_key,
        predecessor_bundle=last.bundles[member_key],
        successor_bundle=last.bundles[member_key],
        accepted_state_advanced=False,
        resource_stop={"reason": reason, "scope": "rsrc1_parent"},
    )
    return store.publish_attempt(
        capability,
        kind=classification.kind,
        member_key=member_key,
        successor_bundle=classification.successor_bundle,
        terminal_evidence=classification.terminal_evidence,
    )


def run_first_event(
    repository_root: Path,
    *,
    authority_commit: str | None = None,
    hooks: RuntimeHooks | None = None,
    host: str = "local",
) -> dict[str, object]:
    """Prospective one-shot parent. Default authority is intentionally absent."""

    hooks = RuntimeHooks() if hooks is None else hooks
    root = _canonical_root(repository_root)
    started = hooks.now()
    store: PRO20RSRC1CampaignStore | None = None
    capability: PRO20RSRC1WriterCapability | None = None
    published = False
    try:
        _require_namespace_absent(root)
        _resource_check(
            root=root,
            started=started,
            published_attempts=0,
            now=hooks.now,
            parent_rss=hooks.parent_rss,
            available_memory=hooks.available_memory,
            disk_free=hooks.disk_free,
            namespace_bytes=hooks.namespace_bytes,
        )
        identity = hooks.implementation_identity_fn(root)
        authority = _authority(
            root=root,
            authority_commit=authority_commit,
            implementation_sha256=str(identity["implementation_sha256"]),
            validator=hooks.authority_validator,
        )
        receipt = _receipt(authority)
        seed_outcomes: dict[str, SeedOutcome] = {}
        for key in MEMBER_KEYS:
            _resource_check(
                root=root,
                started=started,
                published_attempts=0,
                now=hooks.now,
                parent_rss=hooks.parent_rss,
                available_memory=hooks.available_memory,
                disk_free=hooks.disk_free,
                namespace_bytes=hooks.namespace_bytes,
            )
            request = _seed_request(key, receipt=receipt, authority=authority)
            hooks.seed_requests.append(key)
            try:
                completed = hooks.seed_child(root, request)
                outcome = reduce_completed_seed_child(request, completed)
            except PRO20RSRC1ForensicUncertainty as error:
                raise PRO20RSRC1ForensicError(str(error)) from error
            if not outcome.successful:
                reason = str(outcome.resource_evidence.get("reason"))
                raise PRO20RSRC1SeedResourceError(reason)
            seed_outcomes[key] = outcome
        cohort = bind_seed_cohort(seed_outcomes)
        seed_cohort_sha256 = cohort.cohort_sha256
        store = hooks.store_cls.publish_seed_store(
            root,
            receipt=receipt,
            bundles=cohort.bundles,
        )
        del seed_outcomes, cohort
        published = True
        capability = store.acquire_writer(owner_token=OWNER_TOKEN, host=host)
        last = store.authenticate_latest()
        attempts: list[dict[str, object]] = []
        while True:
            if last.terminal or last.first_event_complete:
                break
            try:
                _resource_check(
                    root=root,
                    started=started,
                    published_attempts=len(attempts),
                    now=hooks.now,
                    parent_rss=hooks.parent_rss,
                    available_memory=hooks.available_memory,
                    disk_free=hooks.disk_free,
                    namespace_bytes=hooks.namespace_bytes,
                )
            except PRO20RSRC1ResourceError as error:
                last = _publish_resource_terminal(store, capability, last, str(error))
                attempts.append(
                    {
                        "kind": last.kind,
                        "member_key": last.changed_member_key,
                        "store_generation": last.store_generation,
                        "terminal": True,
                    }
                )
                break
            if len(attempts) >= MAX_PUBLISHED_ATTEMPTS:
                last = _publish_resource_terminal(
                    store, capability, last, "published-attempt ceiling"
                )
                break
            member_key = next_required_member_key(last)
            request = _attempt_request(
                member_key,
                last.bundles[member_key],
                receipt=receipt,
                authority=authority,
            )
            hooks.attempt_requests.append(member_key)
            try:
                completed = hooks.attempt_child(root, request)
                classification, metrics = reduce_completed_child(request, completed)
            except PRO20RSRC1ForensicUncertainty as error:
                raise PRO20RSRC1ForensicError(str(error)) from error
            last = store.publish_attempt(
                capability,
                kind=classification.kind,
                member_key=classification.member_key,
                successor_bundle=classification.successor_bundle,
                terminal_evidence=classification.terminal_evidence,
            )
            attempts.append(
                {
                    "kind": classification.kind,
                    "member_key": classification.member_key,
                    "accepted_state_advanced": classification.accepted_state_advanced,
                    "child_peak_rss_bytes": metrics.peak_rss_bytes,
                    "child_wall_seconds_hex": metrics.wall_seconds.hex(),
                    "store_generation": last.store_generation,
                    "terminal": last.terminal,
                    "first_event_complete": last.first_event_complete,
                }
            )
            del completed, classification, metrics, request
        store.close_writer(capability)
        capability = None
        return {
            "artifact_id": ARTIFACT_ID,
            "schema": SCHEMA,
            "campaign_id": CAMPAIGN_ID,
            "namespace": PRODUCTION_NAMESPACE,
            "seed_cohort_sha256": seed_cohort_sha256,
            "seed_children": list(MEMBER_KEYS),
            "attempts": attempts,
            "published_attempts": len(attempts),
            "terminal": last.terminal,
            "first_event_complete": last.first_event_complete,
            "disposition": last.disposition,
            "checkpoint_sha256": last.checkpoint_sha256,
            "journal_sha256": last.journal_sha256,
            "unclosed_session": False,
            "calibration_eligible": False,
            "campaign_execution_authorized": False,
            "physics_claimed": False,
            "wall_seconds": hooks.now() - started,
        }
    except PRO20RSRC1PostpublicationUncertainty as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise PRO20RSRC1PostpublicationError(str(error)) from error
    except PRO20RSRC1ForensicError:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise
    except (PRO20RSRC1TerminalStore, PRO20RSRC1UnclosedSession) as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise PRO20RSRC1ForensicError(str(error)) from error
    except PRO20RSRC1CampaignStoreError as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        wrapped = PRO20RSRC1PostpublicationError if published else PRO20RSRC1PremiseError
        raise wrapped(str(error)) from error
    except PRO20RSRC1RuntimeError:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise
    except Exception as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise PRO20RSRC1ForensicError(str(error)) from error


def plan_status() -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "classification": "prospective_parent_store_runtime_no_authority",
        "campaign_id": CAMPAIGN_ID,
        "namespace": PRODUCTION_NAMESPACE,
        "future_freeze_artifact_id": FUTURE_FREEZE_ARTIFACT_ID,
        "seed_child_per_member": True,
        "attempt_child_per_schedule": True,
        "parent_is_sole_publisher": True,
        "authority_implemented": False,
        "live_target_present": False,
        "campaign_execution_authorized": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "calibration_eligible": False,
        "mechanism_claimed": False,
        "physics_claimed": False,
    }


__all__ = [
    "ARTIFACT_ID",
    "AUTHORITY_FUNCTION",
    "AUTHORITY_MODULE",
    "CHILD_SCRIPT",
    "FUTURE_FREEZE_ARTIFACT_ID",
    "IMPLEMENTATION_PATHS",
    "MAX_NAMESPACE_BYTES",
    "MAX_PUBLISHED_ATTEMPTS",
    "MAX_TOTAL_WALL_SECONDS",
    "MIN_FREE_DISK_BYTES",
    "PRO20RSRC1AuthorityError",
    "PRO20RSRC1ForensicError",
    "PRO20RSRC1PostpublicationError",
    "PRO20RSRC1PremiseError",
    "PRO20RSRC1ResourceError",
    "PRO20RSRC1RuntimeError",
    "PRO20RSRC1SeedResourceError",
    "RuntimeHooks",
    "available_memory_bytes",
    "current_parent_rss_bytes",
    "free_disk_bytes",
    "implementation_identity",
    "namespace_byte_count",
    "plan_status",
    "require_authority_delta",
    "run_first_event",
]
