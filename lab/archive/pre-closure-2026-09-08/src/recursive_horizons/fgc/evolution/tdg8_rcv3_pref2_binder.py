"""Independent terminal binder for the recovered TDG8 RCV3 event attempt.

PREF2 performs one read-only binding operation over two independently stabilized
snapshots (four non-following low-level tree scans), reconstructs the complete
content-addressed lineage from the authenticated generation-nine boundary
through the exact generation-ten recovery and terminal generation 29, and
safely decodes every persisted array.  It imports no recovery, runner, retry,
or evolution decision code.  Ordinary verification consumes only the compact
tracked result.
"""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib
from typing import Any, Iterable, Mapping

from . import tdg8_rcv3_pref1_binder as pref1


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV3-PREF2"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "outcome_neutral_rcv3_temporal_retry_exhaustion_terminal_binder"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-pref2.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-pref2.json"
OWNER_DOCUMENT = "docs/fgc-tdg8-rcv3-pref2.md"
STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"

AUTHORITY_COMMIT = "98f2e4f20c3c09cf052cf5a8a228d31478842792"
AUTHORITY_BLOBS = (
    (
        "configs/fgc/fgc-1-tdg8-rcv3-rec1-auth1.toml",
        "dab950eb11f8e883222524f6ccbf13702776caf7e8833bbc87fcfae7fa3a4996",
    ),
    (
        "results/fgc-1-tdg8-rcv3-rec1-auth1.json",
        "456ddfb53304025bf53eb8bb58efe0d60ed151eb682a40fd0d8778d9bfda27ae",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg8_rcv3_rec1_authority.py",
        "1988623897e1b7a05d07e120fdf7cd39be6b4a0b7d58be88919efab15d938d54",
    ),
    (
        "scripts/run_fgc_tdg8_rcv3_rec1_event.py",
        "00cc24163d0eb37ca8ac44614b360565f127bb8629ecbd2fba2d008ec8e25a5d",
    ),
)
PREF1_COMMIT = "46abf80753a53e13cf7206fb76d90dc28e4be242"
PREF1_PATH = "src/recursive_horizons/fgc/evolution/tdg8_rcv3_pref1_binder.py"
PREF1_SHA256 = "98caf18141ca83abaf2559ed7eeaa58d42fce4cc1e796fea4fd6171fab01c733"

CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
INTERNAL_AUTHORIZATION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
EVENT = 23
TARGET = {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"}
MEMBER_KEY = "RK4-2049"
DESCRIPTOR_SHA256 = "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
PHYSICAL_STATE_SHA256 = "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"

ANCHOR = {
    "generation": 9,
    "checkpoint_sha256": "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56",
    "checkpoint_raw_sha256": "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846",
    "journal_sequence": 10,
    "journal_sha256": "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff",
}
RECOVERY = {
    "generation": 10,
    "checkpoint_sha256": "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187",
    "checkpoint_raw_sha256": "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b",
    "checkpoint_byte_count": 186_497,
    "parent_sha256": ANCHOR["checkpoint_sha256"],
    "journal_sequence": 12,
    "journal_sha256": "d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4",
    "retry_count": 3,
    "pending_cap_binary64_hex": "0x1.aaa9612df8000p-12",
}
TERMINAL = {
    "generation": 29,
    "checkpoint_sha256": "13eb3ffcbf0a30e588a0edf31867d98c14724f4542cb8272d400a0f52f8bc856",
    "checkpoint_raw_sha256": "e1a012dfac1423d7a70a456aa582f64593564e297e096055b497976ebe5ea43d",
    "checkpoint_byte_count": 1_063_710,
    "parent_sha256": "8e7784d74d2f94538ae2c2cbb8c2f79538d5275943d6dd6f36294612991f4779",
    "journal_sequence": 50,
    "journal_sha256": "29aee78763893fa99e82f8b4b6d2918a4f4a091eb37bad7ad03707a5c33d26da",
    "journal_raw_sha256": "b2b9e71c90be5ebbd30e361bbaccd51898cff5f365a656ee2c4bbe5dde32a530",
    "journal_byte_count": 23_598,
    "rejection_sequence": 49,
    "rejection_sha256": "928cfd1baddf6f4a148e6e6e7277328e50eaeccbf6a45a04c52e9b20cf9117d1",
    "rejection_raw_sha256": "a0bf45404f66c487a76d35b79fa7e879ca4831b0ce9856bb344756bd4560e27a",
    "rejection_byte_count": 23_448,
    "terminal_lock_sha256": "298ec9bd0b57ddc807b8305da2eb1e79a9133c39ce57500ff08eb7f5730922d8",
    "terminal_lock_raw_sha256": "a81c6babfadb699611486e835338abd7ed19016ad075264692fea8c63ad584fd",
    "terminal_lock_byte_count": 297,
    "classification": "invalid",
    "disposition": "invalid_terminal",
    "owner": "temporal",
    "reason": "temporal_retry_exhausted",
    "retry_count": 22,
    "attempted_macro_step_binary64_hex": "0x1.aaa9600000000p-30",
    "retry_step_binary64_hex": "0x1.aaa9600000000p-31",
    "cursor_sha256": "586f7050ee53d4ffcd5cae37f8ad78f173140ab97016178706ff5b48b284fbfa",
    "member_map_sha256": "c81af6f77419ec052416827c2a8aa6f8cb3b631889e13ce6a0cdd78f92df48cd",
}
STORE = {
    "path": STORE_PATH,
    "leaf_count": 115,
    "byte_count": 25_734_805,
    "manifest_sha256": "46b57bc38f3bbfbec70c4b6b11a089cb799dc73135af0d9dac0e3dc289b952d4",
    "checkpoint_count": 30,
    "journal_count": 51,
    "state_count": 17,
    "payload_count": 11,
    "decoded_array_count": 99,
    "rejection_count_after_anchor": 20,
    "cursor_transition_count_after_anchor": 19,
    "terminal_record_count": 1,
}
LOCK_PATHS = (
    "locks/.active-write.lock.hlt16-quarantine-278f0c5aa4cfd86661e8bbaef44e340e39592c9b0206cc8b4dc216f387f29605",
    "locks/.active-write.lock.hlt16-quarantine-8507692d2d970adf895818e81ccc095fad2c657c570aafda22fe5986cab6367d",
    "locks/.active-write.lock.hlt16-quarantine-e2770fba1a988123dc347f8cc8e69bc03eb67b947d205525298ecc1415934f92",
    "locks/terminal.lock",
    "locks/writer.guard",
)

NONCLAIMS = (
    "PREF2 binds a finite numerical/instrumental terminal; it does not diagnose or prescribe a successor remedy.",
    "Temporal retry exhaustion before one accepted post-restart macro-step is not a GR-0, FGC-QR, gradient, mechanism, or physical obstruction.",
    "The terminal preserves the accepted physical state and accumulated temporal debit; no common event, trapping calibration, candidate trajectory, activation, DEF1 result, retained-EFT result, transition, or physical result was earned.",
)
SCOPE = {
    "live_store_read_only": True,
    "compact_verifier_store_blind": True,
    "complete_lineage_reconstructed": True,
    "recovery_runtime_imported": False,
    "evolution_runtime_imported": False,
    "PDE_proposal_executed_by_binder": False,
    "accepted_state_mutated_by_binder": False,
    "candidate_branch_opened": False,
    "successor_remedy_selected": False,
}
CLAIMS = {
    "REC1_authority_commit_bound": True,
    "generation9_anchor_bound": True,
    "exact_generation10_recovery_bound": True,
    "complete_generation9_to_terminal_lineage_bound": True,
    "all_persisted_states_and_payloads_decoded": True,
    "accepted_physical_state_preserved_through_terminal": True,
    "temporal_retry_exhaustion_terminal_bound": True,
    "terminal_classification_numerical_not_physical": True,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "candidate_trajectory_read": False,
    "mechanism_result_earned": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}

_CHECKPOINT_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.json\Z")
_JOURNAL_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.journal\Z")
_STATE_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.json\Z")
_PAYLOAD_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.npz\Z")
_CHANNELS = tuple(f"{kind}:{name}" for kind in ("u", "p", "q") for name in (
    "alpha", "v", "lambda", "R", "phi", "chi",
))


class TDG8RCV3PREF2Error(ValueError):
    """The compact contract or authenticated terminal store differs."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


def _stop(stop_id: str, detail: str) -> None:
    raise TDG8RCV3PREF2Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        _stop("PREF2_CANONICAL_DRIFT", "value is not canonical JSON")
        raise AssertionError from error


def canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        _stop("PREF2_COMPACT_DRIFT", "result is not canonical JSON")
        raise AssertionError from error


def _digest(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(key)
        answer[key] = value
    return answer


def _json(raw: bytes, label: str, *, pretty: bool = False) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"), object_pairs_hook=_duplicates,
            parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("PREF2_JSON_DRIFT", f"{label} is malformed")
        raise AssertionError from error
    expected = canonical_result(value) if pretty else _canonical(value)
    if not isinstance(value, dict) or raw != expected:
        _stop("PREF2_JSON_DRIFT", f"{label} is noncanonical")
    return value


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        _stop("PREF2_CONFIG_DRIFT", "PREF2 config is malformed")
        raise AssertionError from error
    expected_keys = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "store", "anchor",
        "recovery", "terminal", "scope", "claims",
    }
    expected_predecessor = {
        "authority_commit": AUTHORITY_COMMIT,
        "authority_blobs": [
            {"path": path, "sha256": digest} for path, digest in AUTHORITY_BLOBS
        ],
        "pref1_commit": PREF1_COMMIT,
        "pref1_path": PREF1_PATH,
        "pref1_sha256": PREF1_SHA256,
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        _stop("PREF2_CONFIG_DRIFT", "PREF2 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != list(NONCLAIMS)
        or value["predecessor"] != expected_predecessor
        or value["store"] != STORE
        or value["anchor"] != ANCHOR
        or value["recovery"] != RECOVERY
        or value["terminal"] != TERMINAL
        or value["scope"] != SCOPE
        or value["claims"] != CLAIMS
    ):
        _stop("PREF2_CONFIG_DRIFT", "PREF2 exact contract differs")
    return value


def _git(repository: Path, *args: str) -> bytes:
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    try:
        return subprocess.run(
            ["git", "--no-replace-objects", "--no-optional-locks", "-c",
             "core.fsmonitor=false", "-c", "core.untrackedCache=false", *args],
            cwd=repository, check=True, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=environment,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        _stop("PREF2_AUTHORITY_DRIFT", "Git authority lookup failed")
        raise AssertionError from error


def _bind_predecessor(repository: Path) -> dict[str, Any]:
    for commit in (AUTHORITY_COMMIT, PREF1_COMMIT):
        if _git(repository, "rev-parse", "--verify", f"{commit}^{{commit}}") != f"{commit}\n".encode("ascii"):
            _stop("PREF2_AUTHORITY_DRIFT", "predecessor commit differs")
        _git(repository, "merge-base", "--is-ancestor", commit, "HEAD")
    for path, expected in AUTHORITY_BLOBS + ((PREF1_PATH, PREF1_SHA256),):
        commit = PREF1_COMMIT if path == PREF1_PATH else AUTHORITY_COMMIT
        tree = _git(repository, "ls-tree", "-z", commit, "--", path)
        rows = [row for row in tree.split(b"\0") if row]
        if len(rows) != 1 or not rows[0].startswith(b"100644 blob "):
            _stop("PREF2_AUTHORITY_DRIFT", f"authority blob differs: {path}")
        if sha256(_git(repository, "show", f"{commit}:{path}")).hexdigest() != expected:
            _stop("PREF2_AUTHORITY_DRIFT", f"authority hash differs: {path}")
    result_raw = _git(
        repository, "show",
        f"{AUTHORITY_COMMIT}:results/fgc-1-tdg8-rcv3-rec1-auth1.json",
    )
    result = _json(result_raw, "REC1 authority result")
    payload = result.get("artifact_payload")
    if (
        result.get("artifact_id") != "FGC-1-TDG8-RCV3-REC1-AUTH1"
        or result.get("classification")
        != "premise_only_exact_recovery_then_same_event_continuation_authority"
        or not isinstance(payload, Mapping)
        or payload.get("predicted_recovery", {}).get("checkpoint_sha256")
        != RECOVERY["checkpoint_sha256"]
        or payload.get("claims", {}).get("bounded_event23_continuation_authorized") is not True
        or payload.get("claims", {}).get("event24_completed") is not False
        or payload.get("claims", {}).get("candidate_execution_authorized") is not False
    ):
        _stop("PREF2_AUTHORITY_DRIFT", "REC1 authority result meaning differs")
    return {
        "authority_commit": AUTHORITY_COMMIT,
        "authority_blobs": [
            {"path": path, "sha256": digest} for path, digest in AUTHORITY_BLOBS
        ],
        "pref1_commit": PREF1_COMMIT,
        "pref1_path": PREF1_PATH,
        "pref1_sha256": PREF1_SHA256,
    }


def _wrap_pref1(label: str, function: Any, *args: Any) -> Any:
    try:
        return function(*args)
    except pref1.TDG8RCV3PREF1Error as error:
        _stop("PREF2_INDEPENDENT_DECODE_DRIFT", f"{label}: {error.stop_id}")
        raise AssertionError from error


def _validate_tree_grammar(snapshot: Any) -> None:
    if snapshot.directories != pref1.STORE_DIRECTORIES:
        _stop("PREF2_TREE_DRIFT", "store directory grammar differs")
    if tuple(sorted(path for path in snapshot.leaves if path.startswith("locks/"))) != LOCK_PATHS:
        _stop("PREF2_TREE_DRIFT", "terminal lock inventory differs")
    for relative, leaf in snapshot.leaves.items():
        parts = Path(relative).parts
        if len(parts) != 2 or parts[0] not in pref1.STORE_DIRECTORIES:
            _stop("PREF2_TREE_DRIFT", f"nested or foreign leaf appears: {relative}")
        directory, name = parts
        valid = (
            (directory == "checkpoints" and _CHECKPOINT_NAME.fullmatch(name) is not None)
            or (directory == "journal" and _JOURNAL_NAME.fullmatch(name) is not None)
            or (directory == "states" and _STATE_NAME.fullmatch(name) is not None)
            or (directory == "payloads" and _PAYLOAD_NAME.fullmatch(name) is not None)
            or (directory == "receipts" and name == "generation-zero.json")
            or (directory == "locks" and relative in LOCK_PATHS)
        )
        if not valid:
            _stop("PREF2_TREE_DRIFT", f"store leaf grammar differs: {relative}")
        if relative == "locks/writer.guard" and leaf.raw != b"":
            _stop("PREF2_TREE_DRIFT", "writer guard has bytes")
        if relative.startswith("locks/.active-write.lock.hlt16-quarantine-"):
            expected = relative.rsplit("-", 1)[-1]
            if leaf.raw_sha256 != expected:
                _stop("PREF2_TREE_DRIFT", "quarantined lease identity differs")


def _address(value: Mapping[str, Any], field: str, expected: str, label: str) -> None:
    if value.get(field) != expected:
        _stop("PREF2_HASH_DRIFT", f"{label} embedded address differs")
    body = dict(value)
    body.pop(field, None)
    if _digest(body) != expected:
        _stop("PREF2_HASH_DRIFT", f"{label} content address differs")


def _common_checkpoint(checkpoint: Mapping[str, Any], generation: int) -> None:
    if (
        checkpoint.get("generation") != generation
        or checkpoint.get("authorization_commit") != INTERNAL_AUTHORIZATION_COMMIT
        or checkpoint.get("plan_sha256") != PLAN_SHA256
        or checkpoint.get("protocol") != TARGET_PROTOCOL
        or checkpoint.get("campaign_id") != CAMPAIGN_ID
        or checkpoint.get("event") != EVENT
        or checkpoint.get("target") != TARGET
        or not isinstance(checkpoint.get("members"), Mapping)
        or set(checkpoint["members"]) != set(pref1.MEMBER_KEYS)
    ):
        _stop("PREF2_LINEAGE_DRIFT", f"generation {generation} common fields differ")


def _rejection_evidence(record: Mapping[str, Any], retry_count: int) -> Mapping[str, Any]:
    payload = record.get("payload")
    if not isinstance(payload, Mapping):
        _stop("PREF2_LINEAGE_DRIFT", "rejection payload differs")
    evidence = payload.get("evidence")
    if not isinstance(evidence, Mapping):
        _stop("PREF2_LINEAGE_DRIFT", "rejection evidence differs")
    admissions = evidence.get("channel_admissions")
    if not isinstance(admissions, list) or len(admissions) != 18:
        _stop("PREF2_LINEAGE_DRIFT", "TDG6 channel inventory differs")
    channels = [row.get("channel") for row in admissions if isinstance(row, Mapping)]
    failed = [row.get("channel") for row in admissions if isinstance(row, Mapping) and row.get("admission_passed") is False]
    if (
        tuple(channels) != _CHANNELS
        or evidence.get("failed_channels") != failed
        or not failed
        or payload.get("member_key") != MEMBER_KEY
        or payload.get("predecessor_descriptor_sha256") != DESCRIPTOR_SHA256
        or evidence.get("event_type") != "rejected_TDG6_temporal_admission"
        or evidence.get("classification_is_numerical_not_physical") is not True
        or evidence.get("accepted_state_bitwise_preserved") is not True
        or evidence.get("accumulated_debit_unchanged") is not True
        or evidence.get("real_monitor_and_causal_ledgers_preserved") is not True
        or evidence.get("real_tracer_ledger_preserved") is not True
        or evidence.get("initial_state_sha256") != PHYSICAL_STATE_SHA256
        or evidence.get("retry_count_for_current_macro_step") != retry_count
        or evidence.get("cumulative_temporal_retry_count") != retry_count
        or evidence.get("previous_step_index") != 362
        or evidence.get("previous_transaction_serial") != 1810
    ):
        _stop("PREF2_LINEAGE_DRIFT", f"retry {retry_count} evidence differs")
    return evidence


def _validate_accepted_state_preservation(
    member: Mapping[str, Any],
    anchor_member: Mapping[str, Any],
    generation: int,
) -> None:
    """Compare persisted accepted-state accounting without trusting retry evidence."""
    ledger = member.get("ledger")
    anchor_ledger = anchor_member.get("ledger")
    cursor = member.get("cursor")
    anchor_cursor = anchor_member.get("cursor")
    if not all(
        isinstance(value, Mapping)
        for value in (ledger, anchor_ledger, cursor, anchor_cursor)
    ):
        _stop("PREF2_STATE_DRIFT", f"generation {generation} state accounting differs")
    if (
        member.get("descriptor_sha256") != anchor_member.get("descriptor_sha256")
        or cursor.get("accepted_state_sha256")
        != anchor_cursor.get("accepted_state_sha256")
        or cursor.get("accepted_boundary_time")
        != anchor_cursor.get("accepted_boundary_time")
        or cursor.get("previous_step_index")
        != anchor_cursor.get("previous_step_index")
        or cursor.get("previous_transaction_serial")
        != anchor_cursor.get("previous_transaction_serial")
        or cursor.get("committed_common_event_index")
        != anchor_cursor.get("committed_common_event_index")
        or ledger.get("accumulated_debit_vector_hex")
        != anchor_ledger.get("accumulated_debit_vector_hex")
        or ledger.get("accepted_macro_step_count")
        != anchor_ledger.get("accepted_macro_step_count")
        or ledger.get("last_accepted_time_hex")
        != anchor_ledger.get("last_accepted_time_hex")
        or ledger.get("last_accepted_macro_step_temporal_retry_count")
        != anchor_ledger.get("last_accepted_macro_step_temporal_retry_count")
    ):
        _stop("PREF2_STATE_DRIFT", f"generation {generation} accepted state changed")


def _validate_lineage(snapshot: Any) -> dict[str, Any]:
    checkpoints = _wrap_pref1(
        "checkpoint decode", pref1._indexed, snapshot.leaves,
        "checkpoints", _CHECKPOINT_NAME, "checkpoint_sha256",
    )
    journals = _wrap_pref1(
        "journal decode", pref1._indexed, snapshot.leaves,
        "journal", _JOURNAL_NAME, "record_sha256",
    )
    if set(checkpoints) != set(range(30)) or set(journals) != set(range(51)):
        _stop("PREF2_LINEAGE_DRIFT", "checkpoint or journal chain has a gap/fork")

    early_leaves = {
        path: leaf for path, leaf in snapshot.leaves.items()
        if (
            (path.startswith("checkpoints/") and int(Path(path).name[:20]) <= 9)
            or (path.startswith("journal/") and int(Path(path).name[:20]) <= 10)
        )
    }
    early_checkpoints, _ = _wrap_pref1(
        "generation-nine anchor", pref1._validate_lineage, early_leaves,
    )
    descriptors = _wrap_pref1(
        "state descriptors", pref1._descriptor_objects,
        snapshot.leaves, early_checkpoints,
    )
    payloads = _wrap_pref1(
        "payload arrays", pref1._load_payloads, snapshot.leaves, descriptors,
    )
    members = _wrap_pref1(
        "anchor members", pref1._member_evidence,
        early_checkpoints[9], descriptors, payloads,
    )
    if members != list(pref1.EXPECTED_MEMBERS):
        _stop("PREF2_STATE_DRIFT", "generation-nine member evidence differs")

    parent_record = "0" * 64
    for sequence in range(51):
        record = journals[sequence]
        if (
            record.get("schema") != "FGC-1-HLT16-campaign-journal-v1"
            or record.get("sequence") != sequence
            or record.get("campaign_id") != CAMPAIGN_ID
            or record.get("event") != EVENT
            or record.get("previous_record_sha256") != parent_record
        ):
            _stop("PREF2_LINEAGE_DRIFT", f"journal {sequence} relation differs")
        parent_record = str(record["record_sha256"])

    anchor = checkpoints[9]
    if (
        anchor.get("checkpoint_sha256") != ANCHOR["checkpoint_sha256"]
        or anchor.get("journal_sequence") != ANCHOR["journal_sequence"]
        or anchor.get("journal_tip_sha256") != ANCHOR["journal_sha256"]
    ):
        _stop("PREF2_LINEAGE_DRIFT", "generation-nine anchor differs")

    anchor_members = anchor["members"]
    parent_checkpoint = ANCHOR["checkpoint_sha256"]
    previous_cursor = anchor_members[MEMBER_KEY]["cursor"]
    for generation in range(10, 29):
        rejection_sequence = 2 * generation - 9
        transition_sequence = rejection_sequence + 1
        rejection = journals[rejection_sequence]
        transition = journals[transition_sequence]
        retry_count = generation - 7
        checkpoint = checkpoints[generation]
        _common_checkpoint(checkpoint, generation)
        if (
            checkpoint.get("parent_sha256") != parent_checkpoint
            or checkpoint.get("journal_sequence") != transition_sequence
            or checkpoint.get("journal_tip_sha256") != transition.get("record_sha256")
            or checkpoint.get("disposition") != "nonterminal"
            or checkpoint.get("terminal") is not None
            or rejection.get("generation") != generation
            or rejection.get("kind") != "tdg6_rejection"
            or transition.get("generation") != generation
            or transition.get("kind") != "cursor_transition"
        ):
            _stop("PREF2_LINEAGE_DRIFT", f"generation {generation} transition differs")
        evidence = _rejection_evidence(rejection, retry_count)
        rejection_payload = rejection["payload"]
        transition_payload = transition.get("payload")
        if not isinstance(transition_payload, Mapping):
            _stop("PREF2_LINEAGE_DRIFT", "cursor transition payload differs")
        successor = transition_payload.get("successor_cursor")
        member = checkpoint["members"][MEMBER_KEY]
        if not isinstance(successor, Mapping) or not isinstance(member, Mapping):
            _stop("PREF2_LINEAGE_DRIFT", "successor member differs")
        successor_sha = str(transition_payload.get("successor_cursor_sha256"))
        _address(successor, "cursor_chain_sha256", successor_sha, "successor cursor")
        retry = successor.get("retry_successor_payload_or_none")
        if not isinstance(retry, Mapping):
            _stop("PREF2_LINEAGE_DRIFT", "retry successor payload differs")
        ledger = member.get("ledger")
        if not isinstance(ledger, Mapping):
            _stop("PREF2_LINEAGE_DRIFT", "member ledger differs")
        _validate_accepted_state_preservation(
            member, anchor_members[MEMBER_KEY], generation,
        )
        if (
            rejection_payload.get("predecessor_cursor_sha256")
            != previous_cursor.get("cursor_chain_sha256")
            or transition_payload.get("predecessor_cursor_sha256")
            != previous_cursor.get("cursor_chain_sha256")
            or transition_payload.get("rejection_record_sha256")
            != rejection.get("record_sha256")
            or successor.get("cursor_chain_parent_sha256")
            != previous_cursor.get("cursor_chain_sha256")
            or successor.get("journal_tip_sha256") != rejection.get("record_sha256")
            or successor.get("accepted_state_sha256") != DESCRIPTOR_SHA256
            or successor != member.get("cursor")
            or member.get("descriptor_sha256") != DESCRIPTOR_SHA256
            or member.get("pending_owner") != "temporal"
            or retry.get("retry_count") != retry_count
            or retry.get("TDG6_rejection_evidence") != evidence
            or retry.get("durable_rejection_record_sha256")
            != rejection.get("record_sha256")
            or ledger.get("cumulative_temporal_retry_count") != retry_count
            or ledger.get("current_macro_step_temporal_retry_count") != retry_count
            or _digest(ledger) != successor.get("TDG6_ledger_sha256")
        ):
            _stop("PREF2_LINEAGE_DRIFT", f"generation {generation} retry relation differs")
        for key in pref1.MEMBER_KEYS:
            if key != MEMBER_KEY and checkpoint["members"][key] != anchor_members[key]:
                _stop("PREF2_STATE_DRIFT", f"unaffected member changed: {key}")
            if checkpoint["members"][key].get("descriptor_sha256") != anchor_members[key].get("descriptor_sha256"):
                _stop("PREF2_STATE_DRIFT", f"accepted descriptor changed: {key}")
        if generation == RECOVERY["generation"]:
            recovery_path = (
                f"checkpoints/{generation:020d}-{RECOVERY['checkpoint_sha256']}.json"
            )
            leaf = snapshot.leaves.get(recovery_path)
            if (
                checkpoint.get("checkpoint_sha256") != RECOVERY["checkpoint_sha256"]
                or checkpoint.get("parent_sha256") != RECOVERY["parent_sha256"]
                or checkpoint.get("journal_tip_sha256") != RECOVERY["journal_sha256"]
                or leaf is None
                or leaf.raw_sha256 != RECOVERY["checkpoint_raw_sha256"]
                or leaf.byte_count != RECOVERY["checkpoint_byte_count"]
                or member.get("pending_cap_hex") != RECOVERY["pending_cap_binary64_hex"]
            ):
                _stop("PREF2_RECOVERY_DRIFT", "exact predicted generation ten differs")
        previous_cursor = successor
        parent_checkpoint = str(checkpoint["checkpoint_sha256"])

    terminal_checkpoint = checkpoints[29]
    _common_checkpoint(terminal_checkpoint, 29)
    terminal_member = terminal_checkpoint["members"][MEMBER_KEY]
    if not isinstance(terminal_member, Mapping):
        _stop("PREF2_STATE_DRIFT", "generation 29 member differs")
    _validate_accepted_state_preservation(
        terminal_member, anchor_members[MEMBER_KEY], 29,
    )
    rejection = journals[49]
    terminal_record = journals[50]
    evidence = _rejection_evidence(rejection, TERMINAL["retry_count"])
    terminal_payload = terminal_record.get("payload")
    if not isinstance(terminal_payload, Mapping):
        _stop("PREF2_TERMINAL_DRIFT", "terminal payload differs")
    expected_terminal_summary = {
        "classification": TERMINAL["classification"],
        "event": EVENT,
        "member_key": MEMBER_KEY,
        "owner": TERMINAL["owner"],
        "reason": TERMINAL["reason"],
        "target": TARGET,
    }
    if (
        rejection.get("generation") != 29
        or rejection.get("kind") != "tdg6_rejection"
        or rejection["payload"].get("predecessor_cursor_sha256")
        != previous_cursor.get("cursor_chain_sha256")
        or terminal_record.get("generation") != 29
        or terminal_record.get("kind") != "terminal_lock"
        or terminal_record.get("previous_record_sha256") != rejection.get("record_sha256")
        or terminal_payload.get("classification") != TERMINAL["classification"]
        or terminal_payload.get("owner") != TERMINAL["owner"]
        or terminal_payload.get("reason") != TERMINAL["reason"]
        or terminal_payload.get("member_key") != MEMBER_KEY
        or terminal_payload.get("cursor_sha256") != previous_cursor.get("cursor_chain_sha256")
        or terminal_payload.get("descriptor_sha256") != DESCRIPTOR_SHA256
        or terminal_payload.get("rejection_record_sha256") != rejection.get("record_sha256")
        or terminal_payload.get("evidence") != evidence
        or float(evidence.get("attempted_macro_step_size")).hex()
        != TERMINAL["attempted_macro_step_binary64_hex"]
        or float(evidence.get("retry_step_size")).hex()
        != TERMINAL["retry_step_binary64_hex"]
        or terminal_checkpoint.get("parent_sha256") != parent_checkpoint
        or terminal_checkpoint.get("journal_sequence") != 50
        or terminal_checkpoint.get("journal_tip_sha256") != terminal_record.get("record_sha256")
        or terminal_checkpoint.get("disposition") != TERMINAL["disposition"]
        or terminal_checkpoint.get("terminal") != expected_terminal_summary
        or terminal_checkpoint.get("members") != checkpoints[28].get("members")
        or _digest(terminal_checkpoint.get("members")) != TERMINAL["member_map_sha256"]
    ):
        _stop("PREF2_TERMINAL_DRIFT", "terminal lineage or meaning differs")

    for directory, sequence, digest, raw_digest, byte_count in (
        ("checkpoints", 29, TERMINAL["checkpoint_sha256"], TERMINAL["checkpoint_raw_sha256"], TERMINAL["checkpoint_byte_count"]),
        ("journal", 49, TERMINAL["rejection_sha256"], TERMINAL["rejection_raw_sha256"], TERMINAL["rejection_byte_count"]),
        ("journal", 50, TERMINAL["journal_sha256"], TERMINAL["journal_raw_sha256"], TERMINAL["journal_byte_count"]),
    ):
        suffix = "json" if directory == "checkpoints" else "journal"
        path = f"{directory}/{sequence:020d}-{digest}.{suffix}"
        leaf = snapshot.leaves.get(path)
        if leaf is None or leaf.raw_sha256 != raw_digest or leaf.byte_count != byte_count:
            _stop("PREF2_HASH_DRIFT", f"terminal raw object differs: {path}")

    lock_leaf = snapshot.leaves["locks/terminal.lock"]
    lock = _json(lock_leaf.raw, "terminal lock")
    _address(lock, "lock_sha256", TERMINAL["terminal_lock_sha256"], "terminal lock")
    if (
        lock_leaf.raw_sha256 != TERMINAL["terminal_lock_raw_sha256"]
        or lock_leaf.byte_count != TERMINAL["terminal_lock_byte_count"]
        or lock.get("schema") != "FGC-1-HLT16-terminal-lock-v1"
        or lock.get("checkpoint_sha256") != TERMINAL["checkpoint_sha256"]
        or lock.get("journal_tip_sha256") != TERMINAL["journal_sha256"]
    ):
        _stop("PREF2_TERMINAL_DRIFT", "terminal lock differs")

    return {
        "checkpoint_count": len(checkpoints),
        "journal_count": len(journals),
        "state_count": len(snapshot.leaves) - len([
            path for path in snapshot.leaves if not path.startswith("states/")
        ]),
        "payload_count": len(payloads),
        "decoded_array_count": sum(len(arrays) for arrays in payloads.values()),
        "rejection_count_after_anchor": sum(
            journals[n].get("kind") == "tdg6_rejection" for n in range(11, 51)
        ),
        "cursor_transition_count_after_anchor": sum(
            journals[n].get("kind") == "cursor_transition" for n in range(11, 51)
        ),
        "terminal_record_count": sum(
            journals[n].get("kind") == "terminal_lock" for n in range(11, 51)
        ),
    }


def expected_evidence(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "predecessor": dict(config["predecessor"]),
        "store": dict(config["store"]),
        "anchor": dict(config["anchor"]),
        "recovery": dict(config["recovery"]),
        "terminal": dict(config["terminal"]),
    }


def bind_terminal_store(config_raw: bytes, repository: Path) -> dict[str, Any]:
    _config(config_raw)
    predecessor = _bind_predecessor(repository)
    try:
        snapshot = pref1._snapshot(repository, STORE_PATH)
    except pref1.TDG8RCV3PREF1Error as error:
        _stop("PREF2_PATH_UNSAFE", error.stop_id)
        raise AssertionError from error
    if (
        len(snapshot.leaves) != STORE["leaf_count"]
        or snapshot.byte_count() != STORE["byte_count"]
        or snapshot.manifest_sha256() != STORE["manifest_sha256"]
    ):
        _stop("PREF2_TREE_DRIFT", "terminal store manifest differs")
    _validate_tree_grammar(snapshot)
    observed = _validate_lineage(snapshot)
    for key in (
        "checkpoint_count", "journal_count", "state_count", "payload_count",
        "decoded_array_count", "rejection_count_after_anchor",
        "cursor_transition_count_after_anchor", "terminal_record_count",
    ):
        if observed[key] != STORE[key]:
            _stop("PREF2_TREE_DRIFT", f"terminal store count differs: {key}")
    # Close the read window with a second independently stabilized snapshot.
    try:
        if snapshot != pref1._snapshot(repository, STORE_PATH):
            _stop("PREF2_PATH_UNSAFE", "terminal store changed during binding")
    except pref1.TDG8RCV3PREF1Error as error:
        _stop("PREF2_PATH_UNSAFE", error.stop_id)
        raise AssertionError from error
    evidence = expected_evidence(config_raw)
    if evidence["predecessor"] != predecessor:
        _stop("PREF2_AUTHORITY_DRIFT", "predecessor/config evidence differs")
    return evidence


def _result_from_config(config_raw: bytes, config: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            "terminal_evidence": expected_evidence(config_raw),
            "scope": dict(config["scope"]),
            "claims": dict(config["claims"]),
            "nonclaims": list(config["nonclaims"]),
        },
    }


def build_pref2_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    config = _config(config_raw)
    if bind_terminal_store(config_raw, repository) != expected_evidence(config_raw):
        _stop("PREF2_TERMINAL_DRIFT", "live/config terminal binding differs")
    return _result_from_config(config_raw, config)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    value = _json(result_raw, "PREF2 compact result", pretty=True)
    if value != _result_from_config(config_raw, config):
        _stop("PREF2_COMPACT_DRIFT", "PREF2 compact result differs")
    return value


__all__ = [
    "ARTIFACT_ID", "CONFIG_PATH", "RESULT_PATH", "STORE_PATH",
    "TDG8RCV3PREF2Error", "bind_terminal_store", "build_pref2_result",
    "canonical_result", "expected_evidence", "validate_compact_result",
]
