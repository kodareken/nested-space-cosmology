"""Fresh PROTO18 successor bootstrap from the sealed PREF28 terminal.

The invalid generation-nine terminal remains immutable.  This module first
authenticates that terminal, selects its original HLT16 generation-one
ancestor, and restores the six accepted ``t = 23/16`` members into static GR-0
shells.  A successor namespace is seeded separately with the byte-identical
PREF27 generation-zero prefix; the fresh campaign identity begins only at the
HLT16 generation-one bridge.

There is deliberately no proposal, PDE-advance, event-commit, or terminal
repair entry point here.
"""
from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
from typing import Any, Iterable, Mapping

from . import proto17_hlt15_runtime as hlt15
from . import proto19_pref28_binder as pref28
from .hlt16_campaign_runtime import (
    HLT16CampaignRuntimeError,
    restore_member_with_overlay,
)
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .hlt16_state_store import HLT16StateStoreError, read_nofollow
from .proto17_pure_construction import (
    MEMBER_KEYS,
    Proto17ConstructionError,
    build_genesis,
    canonical,
    digest,
)
from .proto18_auth1_inputs import (
    Proto18Auth1InputError,
    source_identities_from_pref26_result,
)
from .proto18_pref27_binder import Pref27Evidence
from .proto19_gr0_static_factory import (
    Proto19GR0StaticFactoryError,
    build_static_gr0_shells,
)
from .proto19_progression_contract import ProgressionPlan, construct_first_event
from .proto19_progression_inputs import (
    ReconstructedGR0MemberSet,
    _check_member,
)


CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
SOURCE_RELATIVE = Path("runs/fgc-2-sf1/proto17/calibration")
DESTINATION_RELATIVE = Path("runs/fgc-2-sf1/proto19/calibration")

GENERATION_ONE_CHECKPOINT_SHA256 = (
    "23f4a1ac64af5e1463862bf67ca240980921611a4a77d31050b0f92c81d44895"
)
GENERATION_ZERO_CHECKPOINT_SHA256 = (
    "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
)
GENERATION_ZERO_RECEIPT_SHA256 = (
    "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf"
)
PREF27_TREE_SHA256 = (
    "95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585"
)

_PREF27_RESULT = "results/fgc-1-pro18-pref27.json"
_PREF27_RESULT_SHA256 = (
    "102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83"
)
_AUTH1_RESULT = "results/fgc-1-pro18-auth1.json"
_AUTH1_RESULT_SHA256 = (
    "dd967568fbae0bd19e063af94ea8e2d53480c98afec5c78bd86601c20faf4efa"
)
_PREF26_RESULT = "results/fgc-1-pro18-pref26.json"
_PREF26_RESULT_SHA256 = (
    "7f56db1555d84d08fdf055ebbe450b6ac8eafc5acff94290e42b5d4d00d47395"
)
_GEN0_CHECKPOINT_RELATIVE = (
    SOURCE_RELATIVE
    / "checkpoints"
    / f"{0:020d}-{GENERATION_ZERO_CHECKPOINT_SHA256}.json"
)
_GEN0_RECEIPT_RELATIVE = SOURCE_RELATIVE / "receipts/generation-zero.json"


class TDG8SuccessorRuntimeError(ValueError):
    """The sealed predecessor cannot seed one fresh PROTO18 successor."""


def _repository_root(root: Path) -> Path:
    repository = Path(os.path.abspath(os.fspath(root)))
    try:
        info = repository.lstat()
    except OSError as exc:
        raise TDG8SuccessorRuntimeError("repository root is absent") from exc
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise TDG8SuccessorRuntimeError("repository root is unsafe")
    return repository


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        value: dict[str, object] = {}
        for key, item in items:
            if key in value:
                raise ValueError(key)
            value[key] = item
        return value

    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=reject_duplicates
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8SuccessorRuntimeError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8SuccessorRuntimeError(f"{label} is not object-valued")
    return value


def _read_exact_json(
    root: Path, relative: str | Path, expected_sha256: str, label: str
) -> tuple[dict[str, Any], bytes]:
    relative_text = Path(relative).as_posix()
    try:
        raw = read_nofollow(root, relative_text, label)
    except HLT16StateStoreError as exc:
        raise TDG8SuccessorRuntimeError(f"{label} cannot be read safely") from exc
    if sha256(raw).hexdigest() != expected_sha256:
        raise TDG8SuccessorRuntimeError(f"{label} bytes differ")
    return _decode_json(raw, label), raw


def _pref27_evidence(
    root: Path,
) -> tuple[Pref27Evidence, Mapping[str, Mapping[str, Any]]]:
    """Rebuild PREF27's exact eight-leaf view from tracked compact evidence."""
    pref27_record, _ = _read_exact_json(
        root, _PREF27_RESULT, _PREF27_RESULT_SHA256, "PREF27 result"
    )
    auth1_record, _ = _read_exact_json(
        root, _AUTH1_RESULT, _AUTH1_RESULT_SHA256, "AUTH1 result"
    )
    pref26_record, _ = _read_exact_json(
        root, _PREF26_RESULT, _PREF26_RESULT_SHA256, "PREF26 result"
    )
    payload = pref27_record.get("artifact_payload")
    auth_payload = auth1_record.get("artifact_payload")
    pref26_payload = pref26_record.get("artifact_payload")
    if (
        pref27_record.get("artifact_id") != "FGC-1-PRO18-PREF27"
        or auth1_record.get("artifact_id") != "FGC-1-PRO18-AUTH1"
        or pref26_record.get("artifact_id") != "FGC-1-PRO18-PREF26"
        or not isinstance(payload, Mapping)
        or not isinstance(auth_payload, Mapping)
        or not isinstance(pref26_payload, Mapping)
    ):
        raise TDG8SuccessorRuntimeError("PREF27 authority identity differs")

    binding = payload.get("authority_binding")
    genesis_spec = auth_payload.get("genesis_spec")
    auth_evidence = auth_payload.get("auth1_input_evidence")
    pref26_binding = auth_payload.get("pref26_binding")
    if (
        not isinstance(binding, Mapping)
        or not isinstance(genesis_spec, Mapping)
        or not isinstance(auth_evidence, Mapping)
        or not isinstance(pref26_binding, Mapping)
        or binding.get("authorization_result_path") != _AUTH1_RESULT
        or binding.get("authorization_result_sha256") != _AUTH1_RESULT_SHA256
        or auth_payload.get("genesis_spec_sha256") != binding.get("genesis_spec_sha256")
        or digest(genesis_spec) != binding.get("genesis_spec_sha256")
        or pref26_binding.get("compact_result") != _PREF26_RESULT
        or pref26_binding.get("compact_result_sha256") != _PREF26_RESULT_SHA256
        or pref26_payload.get("AUTH1_input_evidence") != auth_evidence
    ):
        raise TDG8SuccessorRuntimeError("tracked PREF27/AUTH1/PREF26 binding differs")

    try:
        derived = build_genesis(genesis_spec)
    except (Proto17ConstructionError, TypeError, ValueError) as exc:
        raise TDG8SuccessorRuntimeError("AUTH1 GenesisSpec does not derive") from exc
    checkpoint = derived["checkpoint"]
    if (
        checkpoint["checkpoint_sha256"] != GENERATION_ZERO_CHECKPOINT_SHA256
        or derived["genesis_spec"]["namespace"] != SOURCE_RELATIVE.as_posix()
    ):
        raise TDG8SuccessorRuntimeError("AUTH1 generation-zero identity differs")

    rows = payload.get("store_leaves")
    if not isinstance(rows, list) or len(rows) != 8:
        raise TDG8SuccessorRuntimeError("PREF27 leaf inventory differs")
    leaves: list[dict[str, str]] = []
    for row in rows:
        if (
            not isinstance(row, Mapping)
            or set(row) != {"path", "raw_sha256", "content_sha256"}
            or any(not isinstance(row[name], str) for name in row)
        ):
            raise TDG8SuccessorRuntimeError("PREF27 leaf record differs")
        leaves.append({name: str(row[name]) for name in row})
    ordered = tuple(sorted(leaves, key=lambda row: row["path"]))
    if tuple(leaves) != ordered or digest(list(ordered)) != payload.get("store_tree_sha256"):
        raise TDG8SuccessorRuntimeError("PREF27 leaf tree digest differs")
    if payload.get("store_tree_sha256") != PREF27_TREE_SHA256:
        raise TDG8SuccessorRuntimeError("PREF27 leaf tree identity differs")

    state_paths = {
        (SOURCE_RELATIVE / "states" / f"{address}.json").as_posix()
        for address in derived["store_plan"]["state_objects"]
    }
    expected_paths = state_paths | {
        _GEN0_CHECKPOINT_RELATIVE.as_posix(),
        _GEN0_RECEIPT_RELATIVE.as_posix(),
    }
    if {row["path"] for row in ordered} != expected_paths:
        raise TDG8SuccessorRuntimeError("PREF27 leaf paths differ")

    raw_by_path: dict[str, bytes] = {}
    for row in ordered:
        try:
            raw = read_nofollow(root, row["path"], f"PREF27 leaf {row['path']}")
        except HLT16StateStoreError as exc:
            raise TDG8SuccessorRuntimeError("PREF27 leaf cannot be read safely") from exc
        if sha256(raw).hexdigest() != row["raw_sha256"]:
            raise TDG8SuccessorRuntimeError(f"PREF27 leaf bytes differ: {row['path']}")
        if row["path"] in state_paths and row["content_sha256"] != row["raw_sha256"]:
            raise TDG8SuccessorRuntimeError("PREF27 state address differs")
        raw_by_path[row["path"]] = raw

    checkpoint_raw = raw_by_path[_GEN0_CHECKPOINT_RELATIVE.as_posix()]
    receipt_raw = raw_by_path[_GEN0_RECEIPT_RELATIVE.as_posix()]
    stored_checkpoint = _decode_json(checkpoint_raw, "generation-zero checkpoint")
    receipt = _decode_json(receipt_raw, "generation-zero receipt")
    if canonical(stored_checkpoint) != checkpoint_raw or canonical(receipt) != receipt_raw:
        raise TDG8SuccessorRuntimeError("generation-zero base is noncanonical")
    if stored_checkpoint != checkpoint:
        raise TDG8SuccessorRuntimeError("generation-zero checkpoint differs from AUTH1")
    if (
        payload.get("checkpoint_sha256") != GENERATION_ZERO_CHECKPOINT_SHA256
        or payload.get("receipt_sha256") != GENERATION_ZERO_RECEIPT_SHA256
        or receipt.get("checkpoint_sha256") != GENERATION_ZERO_CHECKPOINT_SHA256
        or receipt.get("receipt_sha256") != GENERATION_ZERO_RECEIPT_SHA256
        or receipt.get("external_authority") != dict(binding)
    ):
        raise TDG8SuccessorRuntimeError("generation-zero receipt binding differs")
    try:
        expected_receipt = hlt15.Proto17HLT15Runtime(
            root,
            derived["genesis_spec"],
            authority_binding=binding,
        )._receipt(checkpoint)
    except hlt15.Proto17HLT15Error as exc:
        raise TDG8SuccessorRuntimeError("generation-zero runtime authority differs") from exc
    if expected_receipt != receipt:
        raise TDG8SuccessorRuntimeError("generation-zero receipt differs from HLT15")

    try:
        identities = source_identities_from_pref26_result(pref26_record, auth_evidence)
    except Proto18Auth1InputError as exc:
        raise TDG8SuccessorRuntimeError("PREF26 source identities differ") from exc
    evidence = Pref27Evidence(
        receipt=receipt,
        checkpoint=stored_checkpoint,
        genesis_spec=derived["genesis_spec"],
        authority_binding=dict(binding),
        leaves=ordered,
        tree_sha256=PREF27_TREE_SHA256,
    )
    return evidence, identities


def authenticate_and_reconstruct(root: Path) -> ReconstructedGR0MemberSet:
    """Authenticate PREF28 and restore the immutable generation-one members."""
    repository = _repository_root(root)
    try:
        config_raw = read_nofollow(
            repository, pref28.CONFIG_PATH, "PREF28 config"
        )
        pref28.bind_terminal_store(config_raw, repository)
        evidence, identities = _pref27_evidence(repository)
        store = HLT16CampaignStore(repository / SOURCE_RELATIVE)
        checkpoint = store.authenticated_checkpoint_at_generation(1)
    except (
        HLT16CampaignStoreError,
        HLT16StateStoreError,
        pref28.PREF28BinderError,
    ) as exc:
        raise TDG8SuccessorRuntimeError("sealed PREF28 predecessor differs") from exc
    if (
        checkpoint.generation != 1
        or checkpoint.sha256 != GENERATION_ONE_CHECKPOINT_SHA256
        or checkpoint.protocol != "FGC-2-SF1-PROTO18"
        or checkpoint.campaign_id != pref28.CAMPAIGN_ID
        or checkpoint.event != 23
        or checkpoint.target
        != {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise TDG8SuccessorRuntimeError("historical generation-one anchor differs")

    try:
        shells = build_static_gr0_shells(repository)
        members: dict[str, object] = {}
        for key in MEMBER_KEYS:
            member = restore_member_with_overlay(store, checkpoint, shells[key], key=key)
            _check_member(key, member, evidence)
            members[key] = member
        return ReconstructedGR0MemberSet(
            evidence=evidence,
            members=members,
            source_identity=identities,
        )
    except (
        HLT16CampaignRuntimeError,
        Proto19GR0StaticFactoryError,
        ValueError,
        KeyError,
    ) as exc:
        if isinstance(exc, TDG8SuccessorRuntimeError):
            raise
        raise TDG8SuccessorRuntimeError("generation-one member restoration differs") from exc


def construct_successor_plan(root: Path, authorization_commit: str) -> ProgressionPlan:
    """Construct the fixed fresh-campaign PROTO18 event-23 plan."""
    repository = _repository_root(root)
    plan = replace(
        construct_first_event(
            repository, authorization_commit=authorization_commit
        ),
        campaign_id=CAMPAIGN_ID,
    )
    if (
        plan.protocol_artifact_id != "FGC-2-SF1-PROTO18"
        or plan.campaign_id != CAMPAIGN_ID
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
        or plan.common_event_index != 23
        or plan.start_time.rational != "23/16"
        or plan.target_time.rational != "3/2"
        or tuple(member.member_key for member in plan.members) != MEMBER_KEYS
    ):
        raise TDG8SuccessorRuntimeError("successor progression plan differs")
    return plan


def _sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise TDG8SuccessorRuntimeError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as exc:
        raise TDG8SuccessorRuntimeError(f"{label} is not SHA-256") from exc
    return value


def _destination(root: Path, destination: Path) -> Path:
    supplied = Path(destination)
    if not supplied.is_absolute():
        supplied = root / supplied
    supplied = Path(os.path.abspath(os.fspath(supplied)))
    expected = root / DESTINATION_RELATIVE
    if supplied != expected:
        raise TDG8SuccessorRuntimeError("successor destination differs")
    return supplied


def _validate_source(source: ReconstructedGR0MemberSet) -> None:
    if not isinstance(source, ReconstructedGR0MemberSet):
        raise TDG8SuccessorRuntimeError("successor source type differs")
    evidence = source.evidence
    if (
        evidence.tree_sha256 != PREF27_TREE_SHA256
        or evidence.checkpoint.get("checkpoint_sha256")
        != GENERATION_ZERO_CHECKPOINT_SHA256
        or evidence.receipt.get("receipt_sha256") != GENERATION_ZERO_RECEIPT_SHA256
        or evidence.genesis_spec.get("namespace") != SOURCE_RELATIVE.as_posix()
        or tuple(source.members) != MEMBER_KEYS
        or source.branch != "GR-0"
        or source.amplitude != "3"
        or source.restart_time.hex() != "0x1.7000000000000p+0"
        or source.target_time.hex() != "0x1.8000000000000p+0"
    ):
        raise TDG8SuccessorRuntimeError("successor source identity differs")
    for key in MEMBER_KEYS:
        _check_member(key, source.members[key], evidence)


def _verify_staged_prefix(
    store: Path,
    source: ReconstructedGR0MemberSet,
    materialized: Mapping[str, Any],
) -> None:
    prefix = SOURCE_RELATIVE.as_posix() + "/"
    expected: dict[str, str] = {}
    for row in source.evidence.leaves:
        path = str(row["path"])
        if not path.startswith(prefix):
            raise TDG8SuccessorRuntimeError("PREF27 staged path differs")
        expected[path.removeprefix(prefix)] = str(row["raw_sha256"])
    observed: dict[str, str] = {}
    directories: set[str] = set()
    pending = [store]
    while pending:
        directory = pending.pop()
        try:
            entries = sorted(os.scandir(directory), key=lambda entry: entry.name)
        except OSError as exc:
            raise TDG8SuccessorRuntimeError("private HLT15 stage is unreadable") from exc
        for entry in entries:
            path = Path(entry.path)
            relative = path.relative_to(store).as_posix()
            info = entry.stat(follow_symlinks=False)
            if stat.S_ISLNK(info.st_mode):
                raise TDG8SuccessorRuntimeError("private HLT15 stage contains a symlink")
            if stat.S_ISDIR(info.st_mode):
                directories.add(relative)
                pending.append(path)
            elif stat.S_ISREG(info.st_mode) and info.st_nlink == 1:
                observed[relative] = sha256(path.read_bytes()).hexdigest()
            else:
                raise TDG8SuccessorRuntimeError("private HLT15 stage is unsafe")
    if directories != {"checkpoints", "receipts", "states"} or observed != expected:
        raise TDG8SuccessorRuntimeError("private HLT15 eight-leaf prefix differs")
    if (
        materialized.get("checkpoint") != source.evidence.checkpoint
        or materialized.get("receipt") != source.evidence.receipt
    ):
        raise TDG8SuccessorRuntimeError("private HLT15 materialization differs")


def install_generation_zero_prefix(
    root: Path,
    destination: Path,
    source: ReconstructedGR0MemberSet,
    plan_sha256: str,
) -> Mapping[str, Any]:
    """Install only the exact PREF27 base into the absent successor namespace.

    HLT15 first derives and verifies the eight leaves inside a private
    repository on the destination filesystem.  Darwin ``RENAME_EXCL`` then
    adopts that private store as one atomic, no-clobber operation.  A target
    that already exists is always foreign and is never inspected or resumed.
    """
    repository = _repository_root(root)
    target = _destination(repository, destination)
    plan_digest = _sha256(plan_sha256, "successor plan")
    try:
        _validate_source(source)
    except (ValueError, KeyError, TypeError) as exc:
        if isinstance(exc, TDG8SuccessorRuntimeError):
            raise
        raise TDG8SuccessorRuntimeError("successor source validation differs") from exc

    try:
        checked_target = hlt15._safe_parent_chain(
            repository, DESTINATION_RELATIVE.as_posix()
        )
        if checked_target != target:
            raise TDG8SuccessorRuntimeError("successor destination resolution differs")
        hlt15._lstat_absent(target)
    except hlt15.Proto17HLT15Error as exc:
        raise TDG8SuccessorRuntimeError(
            "successor destination is unsafe or already exists"
        ) from exc

    stage_repository = Path(
        tempfile.mkdtemp(
            prefix=f".tdg8-successor-{plan_digest[:12]}-",
            dir=target.parent,
        )
    )
    adopted = False
    try:
        runtime = hlt15.Proto17HLT15Runtime(
            stage_repository,
            source.evidence.genesis_spec,
            authority_binding=source.evidence.authority_binding,
        )
        materialized = runtime.materialize_generation_zero()
        staged_store = stage_repository / SOURCE_RELATIVE
        _verify_staged_prefix(staged_store, source, materialized)
        try:
            hlt15._rename_exclusive(staged_store, target)
        except hlt15.Proto17HLT15Error as exc:
            raise TDG8SuccessorRuntimeError(
                "successor destination appeared or exclusive install failed"
            ) from exc
        adopted = True
        parent = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
        return materialized
    except hlt15.Proto17HLT15Error as exc:
        raise TDG8SuccessorRuntimeError("private HLT15 materialization failed") from exc
    finally:
        # This path is created by this invocation with mode 0700.  Once the
        # store is adopted it contains only empty HLT15 parent directories;
        # before adoption it contains only this invocation's private copy.
        if stage_repository.exists():
            shutil.rmtree(stage_repository)
        if adopted and not target.exists():
            raise TDG8SuccessorRuntimeError("successor destination vanished after adoption")


__all__ = [
    "CAMPAIGN_ID",
    "DESTINATION_RELATIVE",
    "SOURCE_RELATIVE",
    "TDG8SuccessorRuntimeError",
    "authenticate_and_reconstruct",
    "construct_successor_plan",
    "install_generation_zero_prefix",
]
