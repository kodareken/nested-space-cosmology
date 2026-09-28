#!/usr/bin/env python3
"""Run or inspect the single fresh TDG8 GR-0 successor event.

The command owns exactly one event, ``23 -> 24``, in the fixed
``runs/fgc-2-sf1/proto19/calibration`` namespace.  It can bootstrap that
namespace only from the authenticated generation-zero source selected by the
TDG8 successor authority.  Once generation one exists, it uses persisted
HLT16 state plus fresh GR-0 static shells and never reopens the historical raw
source.  No candidate branch, calibration conclusion, or physical result is
available through this command.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
# Select the tracked ``scripts`` owner from the repository first, while
# deliberately removing the exact src entry so that this reused owner installs
# its own guarded [src, repository] order before it imports runtime modules.
# This keeps ignored untracked paths from becoming import authority.
sys.path[:] = [entry for entry in sys.path if entry != _SOURCE_ROOT]
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_pro19_event1 as event1  # noqa: E402

# A cached owner has already imported its dependencies, so normalize the
# remaining imports explicitly as well.  The two exact authority roots win
# over the working directory, PYTHONPATH aliases, and site packages.
sys.path[:] = [entry for entry in sys.path if entry not in {_SOURCE_ROOT, str(ROOT)}]
sys.path[:0] = [_SOURCE_ROOT, str(ROOT)]

from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_successor_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_successor_runtime as tdg8  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import (  # noqa: E402
    MEMBER_KEYS,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)


RUNNER_ID = "FGC-1-TDG8-SUCCESSOR-EVENT-RUNNER"
SCHEMA = "FGC-1-TDG8-successor-event-runner-v1"
EXPECTED_DESTINATION_RELATIVE = Path("runs/fgc-2-sf1/proto19/calibration")
EXPECTED_SOURCE_RELATIVE = Path("runs/fgc-2-sf1/proto17/calibration")
EXPECTED_PROTOCOL = "FGC-2-SF1-PROTO18"
EXPECTED_BRANCH = "GR-0"
EXPECTED_AMPLITUDE = "3"
EXPECTED_EVENT = 23
EXPECTED_START = ("23/16", (23 / 16).hex())
EXPECTED_TARGET = ("3/2", (3 / 2).hex())
EXPECTED_SUCCESSOR_TARGET = ("25/16", (25 / 16).hex())
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class TDG8SuccessorRunnerError(RuntimeError):
    """The fixed successor event cannot be run without ambiguity."""


def _canonical(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _git(root: Path, *arguments: str) -> subprocess.CompletedProcess[bytes]:
    """Run one read-only Git query without ambient Git authority."""

    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    try:
        return subprocess.run(
            (
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                "-C",
                str(root),
                *arguments,
            ),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            env=environment,
        )
    except OSError as exc:
        raise TDG8SuccessorRunnerError("Git authority query could not start") from exc


def _fixed_root(root: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        raise TDG8SuccessorRunnerError("repository root is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise TDG8SuccessorRunnerError("repository root is unsafe")
    canonical = ROOT.resolve()
    if supplied.resolve() != canonical:
        raise TDG8SuccessorRunnerError("successor runner requires its fixed repository root")
    return canonical


def _require_git_authority(root: Path, authorization_commit: str) -> str:
    """Require the supplied commit to be live HEAD with no tracked drift.

    Untracked files are intentionally outside this check.  They do not alter
    the committed execution image, and the fixed campaign/store validators own
    every namespace leaf they consume.
    """

    if not isinstance(authorization_commit, str) or not _COMMIT.fullmatch(
        authorization_commit
    ):
        raise TDG8SuccessorRunnerError("authority commit is not a full lowercase Git id")
    head_result = _git(root, "rev-parse", "--verify", "HEAD")
    try:
        head = head_result.stdout.decode("ascii").strip()
    except UnicodeDecodeError as exc:
        raise TDG8SuccessorRunnerError("repository HEAD is not ASCII") from exc
    if head_result.returncode != 0 or not _COMMIT.fullmatch(head):
        raise TDG8SuccessorRunnerError("repository HEAD cannot be resolved")
    if head != authorization_commit:
        raise TDG8SuccessorRunnerError("authority commit differs from live HEAD")
    tracked = _git(
        root,
        "diff",
        "--quiet",
        "--no-ext-diff",
        "--ignore-submodules=all",
        "HEAD",
        "--",
    )
    if tracked.returncode not in {0, 1}:
        raise TDG8SuccessorRunnerError("tracked-tree cleanliness cannot be determined")
    if tracked.returncode == 1:
        raise TDG8SuccessorRunnerError("tracked tree differs from authority commit")
    return head


def _require_execution_authority(
    receipt: object,
    *,
    requested_commit: str,
    config_raw: bytes,
    result_raw: bytes,
    config: Mapping[str, object],
) -> authority.TDG8SuccessorExecutionAuthority:
    """Validate every runner-consumed field of the typed authority receipt."""

    if not isinstance(receipt, authority.TDG8SuccessorExecutionAuthority):
        raise TDG8SuccessorRunnerError("TDG8 execution authority type differs")
    try:
        configured_inventory = config["implementation_inventory"]
        if not isinstance(configured_inventory, list):
            raise TypeError
        expected_inventory = tuple(
            authority.ImplementationBinding(
                role=str(item["role"]),
                path=str(item["path"]),
                sha256=str(item["sha256"]),
            )
            for item in configured_inventory
        )
    except (KeyError, TypeError) as exc:
        raise TDG8SuccessorRunnerError(
            "TDG8 execution implementation binding differs"
        ) from exc
    observed = (
        receipt.authority_commit,
        receipt.predecessor_repair_commit,
        receipt.config_sha256,
        receipt.result_sha256,
        receipt.campaign_id,
        receipt.target_protocol,
        receipt.branch,
        receipt.amplitude,
        receipt.event,
        receipt.start_rational,
        receipt.target_rational,
        receipt.event_count,
        receipt.source_campaign_path,
        receipt.destination_campaign_path,
        receipt.source_terminal_read_only,
        receipt.candidate_branches_forbidden,
        receipt.immutable_bindings,
        receipt.implementation_inventory,
        dict(receipt.environment),
    )
    expected = (
        requested_commit,
        authority.PREDECESSOR_REPAIR_COMMIT,
        sha256(config_raw).hexdigest(),
        sha256(result_raw).hexdigest(),
        tdg8.CAMPAIGN_ID,
        EXPECTED_PROTOCOL,
        EXPECTED_BRANCH,
        EXPECTED_AMPLITUDE,
        EXPECTED_EVENT,
        EXPECTED_START[0],
        EXPECTED_TARGET[0],
        1,
        EXPECTED_SOURCE_RELATIVE.as_posix(),
        EXPECTED_DESTINATION_RELATIVE.as_posix(),
        True,
        True,
        authority.IMMUTABLE_BINDINGS,
        expected_inventory,
        authority.PINNED_ENVIRONMENT,
    )
    if observed != expected:
        raise TDG8SuccessorRunnerError("TDG8 execution authority scope differs")
    return receipt


def _execution_authority(
    root: Path,
    authorization_commit: str,
) -> authority.TDG8SuccessorExecutionAuthority:
    """Consume the compact committed authority without inspecting a store."""

    try:
        config_raw = authority._read_repository_leaf(
            root,
            authority.CONFIG_PATH,
            label="TDG8 successor execution config",
            maximum_bytes=1024 * 1024,
        )
        result_raw = authority._read_repository_leaf(
            root,
            authority.RESULT_PATH,
            label="TDG8 successor execution result",
            maximum_bytes=4 * 1024 * 1024,
        )
        receipt = authority.authorize_execution(
            root,
            config_raw,
            result_raw,
            authorization_commit,
        )
        config = authority._parse_config(config_raw)
    except authority.TDG8SuccessorAuthorityError as exc:
        raise TDG8SuccessorRunnerError("TDG8 execution authority failed") from exc
    return _require_execution_authority(
        receipt,
        requested_commit=authorization_commit,
        config_raw=config_raw,
        result_raw=result_raw,
        config=config,
    )


def _fixed_destination(root: Path) -> Path:
    relative = Path(tdg8.DESTINATION_RELATIVE)
    if (
        relative != EXPECTED_DESTINATION_RELATIVE
        or relative.is_absolute()
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise TDG8SuccessorRunnerError("TDG8 destination authority differs")
    destination = root.joinpath(*relative.parts)
    # ``strict=False`` still resolves every existing parent.  This prevents a
    # parent symlink from redirecting the fixed campaign namespace.
    if destination.resolve(strict=False) != destination:
        raise TDG8SuccessorRunnerError("TDG8 destination is redirected")
    return destination


def _destination_state(destination: Path) -> str:
    try:
        metadata = destination.lstat()
    except FileNotFoundError:
        return "absent"
    except OSError as exc:
        raise TDG8SuccessorRunnerError("TDG8 destination cannot be inspected") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise TDG8SuccessorRunnerError("TDG8 destination is not a direct directory")
    return "directory"


def _time_pair(value: object) -> tuple[str, str]:
    try:
        return str(value.rational), str(value.binary64_hex)
    except AttributeError as exc:
        raise TDG8SuccessorRunnerError("successor plan time identity differs") from exc


def _require_plan_scope(plan: object, authorization_commit: str) -> Any:
    try:
        members = tuple(item.member_key for item in plan.members)
        observed = (
            plan.authorization_commit,
            plan.protocol_artifact_id,
            plan.campaign_id,
            plan.branch,
            plan.amplitude,
            plan.common_event_index,
            _time_pair(plan.start_time),
            _time_pair(plan.target_time),
            members,
        )
    except (AttributeError, TypeError) as exc:
        raise TDG8SuccessorRunnerError("successor progression plan is malformed") from exc
    expected = (
        authorization_commit,
        EXPECTED_PROTOCOL,
        tdg8.CAMPAIGN_ID,
        EXPECTED_BRANCH,
        EXPECTED_AMPLITUDE,
        EXPECTED_EVENT,
        EXPECTED_START,
        EXPECTED_TARGET,
        MEMBER_KEYS,
    )
    if observed != expected:
        raise TDG8SuccessorRunnerError("successor plan crossed its GR-0 event scope")
    digest = getattr(plan, "sha256", None)
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise TDG8SuccessorRunnerError("successor plan digest differs")
    return plan


def _require_source_scope(source: object) -> Mapping[str, object]:
    try:
        members = source.members
        identity = (source.branch, source.amplitude, tuple(members))
    except (AttributeError, TypeError) as exc:
        raise TDG8SuccessorRunnerError("authenticated generation-zero source differs") from exc
    if (
        not isinstance(members, Mapping)
        or identity != (EXPECTED_BRANCH, EXPECTED_AMPLITUDE, MEMBER_KEYS)
    ):
        raise TDG8SuccessorRunnerError("authenticated source crossed its GR-0 scope")
    return members


def _require_checkpoint_scope(checkpoint: object, plan: object) -> None:
    try:
        event = checkpoint.event
        target = checkpoint.target
        disposition = checkpoint.disposition
        generation = checkpoint.generation
        member_keys = tuple(checkpoint.members)
        identity = (
            checkpoint.authorization_commit,
            checkpoint.plan_sha256,
            checkpoint.campaign_id,
            checkpoint.protocol,
        )
    except (AttributeError, TypeError) as exc:
        raise TDG8SuccessorRunnerError("successor checkpoint is malformed") from exc
    if (
        not isinstance(target, Mapping)
        or isinstance(generation, bool)
        or not isinstance(generation, int)
        or generation < 1
        or identity != (
            plan.authorization_commit,
            plan.sha256,
            tdg8.CAMPAIGN_ID,
            EXPECTED_PROTOCOL,
        )
        or member_keys != MEMBER_KEYS
        or event not in {23, 24}
    ):
        raise TDG8SuccessorRunnerError("successor checkpoint crossed its event scope")
    if (
        event == 23
        and disposition not in {"nonterminal", "scientific_terminal", "invalid_terminal"}
    ) or (event == 24 and disposition != "event_complete"):
        raise TDG8SuccessorRunnerError("successor checkpoint disposition differs")
    expected_target = EXPECTED_TARGET if event == 23 else EXPECTED_SUCCESSOR_TARGET
    if (target.get("rational"), target.get("binary64_hex")) != expected_target:
        raise TDG8SuccessorRunnerError("successor checkpoint target differs")


def _require_terminal_result(result: object, plan: object) -> dict[str, Any]:
    if not isinstance(result, Mapping):
        raise TDG8SuccessorRunnerError("successor result is not object-valued")
    answer = dict(result)
    event = answer.get("event")
    state = answer.get("state")
    target = answer.get("target")
    if (
        answer.get("authorization_commit") != plan.authorization_commit
        or answer.get("plan_sha256") != plan.sha256
        or answer.get("campaign_id") != tdg8.CAMPAIGN_ID
        or answer.get("candidate_branch_opened") is not False
        or answer.get("calibration_result_earned") is not False
        or answer.get("physical_result_earned") is not False
    ):
        raise TDG8SuccessorRunnerError("successor result crossed its claim scope")
    if not isinstance(target, Mapping):
        raise TDG8SuccessorRunnerError("successor result target differs")
    if event == 24:
        if state != "event_one_complete" or (
            target.get("rational"),
            target.get("binary64_hex"),
        ) != EXPECTED_SUCCESSOR_TARGET:
            raise TDG8SuccessorRunnerError("event 24 result is not the single completed edge")
    elif event == 23:
        if state not in {"scientific_terminal", "invalid_terminal"} or (
            target.get("rational"),
            target.get("binary64_hex"),
        ) != EXPECTED_TARGET:
            raise TDG8SuccessorRunnerError("event 23 result is neither complete nor terminal")
    else:
        raise TDG8SuccessorRunnerError("successor runner advanced beyond one event")
    return {
        **answer,
        "schema": SCHEMA,
        "runner_id": RUNNER_ID,
        "destination": EXPECTED_DESTINATION_RELATIVE.as_posix(),
        "branch": EXPECTED_BRANCH,
        "amplitude": EXPECTED_AMPLITUDE,
        "event_limit": 1,
    }


def run_successor_event(
    root: Path,
    *,
    authorization_commit: str,
    attempt_runner: Callable[..., Any] | None = None,
) -> dict[str, Any]:
    """Create or resume exactly the fresh TDG8 event-23 successor."""

    repository = _fixed_root(root)
    execution_authority = _execution_authority(repository, authorization_commit)
    authority_commit = execution_authority.authority_commit
    _require_git_authority(repository, authority_commit)
    destination = _fixed_destination(repository)
    state = _destination_state(destination)
    plan = _require_plan_scope(
        tdg8.construct_successor_plan(repository, authority_commit),
        authority_commit,
    )

    source: object | None = None
    if state == "absent":
        source = tdg8.authenticate_and_reconstruct(repository)
        _require_source_scope(source)
        # Recheck the committed image immediately before the first mutation.
        _require_git_authority(repository, authority_commit)
        tdg8.install_generation_zero_prefix(
            repository,
            destination,
            source,
            plan.sha256,
        )

    store = HLT16CampaignStore(destination)
    bootstrap_prefix = event1._recoverable_bootstrap_prefix(store)
    if state == "absent" and not bootstrap_prefix:
        raise TDG8SuccessorRunnerError(
            "fresh generation-zero installation did not retain its sole bootstrap prefix"
        )
    if bootstrap_prefix:
        if source is None:
            source = tdg8.authenticate_and_reconstruct(repository)
            _require_source_scope(source)
        _require_git_authority(repository, authority_commit)
        checkpoint = runtime.initialize_or_recover_generation_one(plan, store, source)
        templates = _require_source_scope(source)
    else:
        # Once a generation-one checkpoint exists, no historical source array
        # is reopened.  HLT16 owns the persisted state and the factory supplies
        # only the fixed GR-0 types/operators into which it is restored.
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
        try:
            templates = build_static_gr0_shells(repository)
        except (OSError, TypeError, ValueError) as exc:
            raise TDG8SuccessorRunnerError("persisted GR-0 shell construction failed") from exc
        if tuple(templates) != MEMBER_KEYS:
            raise TDG8SuccessorRunnerError("persisted GR-0 shell order differs")

    _require_checkpoint_scope(checkpoint, plan)
    _require_git_authority(repository, authority_commit)
    result = event1._advance_authenticated_event(
        plan,
        store,
        checkpoint,
        templates,
        attempt_runner=attempt_runner,
    )
    return _require_terminal_result(result, plan)


def collect_status(root: Path, *, authorization_commit: str) -> dict[str, Any]:
    """Return one read-only view of the fixed successor namespace."""

    repository = _fixed_root(root)
    execution_authority = _execution_authority(repository, authorization_commit)
    authority_commit = execution_authority.authority_commit
    _require_git_authority(repository, authority_commit)
    destination = _fixed_destination(repository)
    plan = _require_plan_scope(
        tdg8.construct_successor_plan(repository, authority_commit),
        authority_commit,
    )
    if _destination_state(destination) == "absent":
        return {
            "schema": SCHEMA,
            "runner_id": RUNNER_ID,
            "state": "destination_absent",
            "destination": EXPECTED_DESTINATION_RELATIVE.as_posix(),
            "authorization_commit": authority_commit,
            "plan_sha256": plan.sha256,
            "campaign_id": tdg8.CAMPAIGN_ID,
            "safe_to_start": True,
            "safe_to_restart": False,
            "candidate_branch_opened": False,
        }

    store = HLT16CampaignStore(destination)
    if event1._recoverable_bootstrap_prefix(store):
        return {
            "schema": SCHEMA,
            "runner_id": RUNNER_ID,
            "state": "authenticated_bootstrap_prefix",
            "destination": EXPECTED_DESTINATION_RELATIVE.as_posix(),
            "authorization_commit": authority_commit,
            "plan_sha256": plan.sha256,
            "campaign_id": tdg8.CAMPAIGN_ID,
            "safe_to_start": True,
            "safe_to_restart": False,
            "candidate_branch_opened": False,
        }

    snapshot = store.authenticated_snapshot()
    recovery = store.inspect_recovery()
    checkpoint = snapshot.checkpoint
    _require_checkpoint_scope(checkpoint, plan)
    if (
        recovery.authorization_commit != authority_commit
        or recovery.plan_sha256 != plan.sha256
        or recovery.checkpoint_sha256 != checkpoint.sha256
    ):
        raise TDG8SuccessorRunnerError("status recovery authority differs")
    recovery_mapping = asdict(recovery) if is_dataclass(recovery) else dict(vars(recovery))
    return {
        "schema": SCHEMA,
        "runner_id": RUNNER_ID,
        "state": recovery.state,
        "destination": EXPECTED_DESTINATION_RELATIVE.as_posix(),
        "authorization_commit": authority_commit,
        "plan_sha256": plan.sha256,
        "campaign_id": tdg8.CAMPAIGN_ID,
        "checkpoint": {
            "generation": checkpoint.generation,
            "sha256": checkpoint.sha256,
            "event": checkpoint.event,
            "disposition": checkpoint.disposition,
        },
        "recovery": recovery_mapping,
        "safe_to_start": False,
        "safe_to_restart": recovery.safe_to_restart,
        "candidate_branch_opened": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--status", action="store_true")
    args = parser.parse_args()
    try:
        result = (
            run_successor_event(ROOT, authorization_commit=args.authority_commit)
            if args.run
            else collect_status(ROOT, authorization_commit=args.authority_commit)
        )
    except (
        TDG8SuccessorRunnerError,
        event1.Proto19Event1RunnerError,
        HLT16CampaignStoreError,
        OSError,
        ValueError,
    ) as exc:
        print(
            _canonical(
                {
                    "schema": SCHEMA,
                    "runner_id": RUNNER_ID,
                    "state": "invalid_implementation_or_nonconverged_run",
                    "detail": str(exc),
                    "safe_to_restart": False,
                    "candidate_branch_opened": False,
                    "physical_result_earned": False,
                }
            )
        )
        return 2
    print(_canonical(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
