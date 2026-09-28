"""Independent read-only binder for the installed TDG8 RCV3 projection.

The bootstrap writer is deliberately outside this module's import graph.  This
binder reconstructs the terminal source tree, the projected generation-nine
tree, the external installation receipt, the complete HLT16 lineage, and every
persisted NumPy payload directly from bytes.  It exposes no evolution, retry,
candidate, or publication entry point.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import re
import stat
import struct
import subprocess
import tomllib
from typing import Any, Iterable, Mapping, Sequence
import zipfile

import numpy as np


ARTIFACT_ID = "FGC-1-TDG8-RCV3-PREF1"
SCHEMA_VERSION = 1
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "read_only_generation9_recovery_projection_binder"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-pref1.json"

SOURCE_STORE = "runs/fgc-2-sf1/proto19/calibration"
DESTINATION_WRAPPER = "runs/fgc-2-sf1/tdg8-rcv3"
DESTINATION_STORE = f"{DESTINATION_WRAPPER}/calibration"
RECEIPT_PATH = f"{DESTINATION_WRAPPER}/bootstrap-receipt.json"

AUTHORITY_COMMIT = "ec24cd25bc2d226e4bb26188c10d742bd9ff3b1f"
FRZ1_RESULT_PATH = "results/fgc-1-tdg8-rcv3-frz1.json"
FRZ1_RESULT_SHA256 = "4693cba6c3f990ce774ed04998392b1aff49c0a956d85e81646ac6d5ec91c9c9"
RUNTIME_PATH = "src/recursive_horizons/fgc/evolution/tdg8_rcv3_fork_runtime.py"
RUNTIME_SHA256 = "9c7e40dd76d126fa4ce708025483d0c15a65e9f36a80c5919335bd53c73dca36"
BOOTSTRAP_PATH = "scripts/bootstrap_fgc_tdg8_rcv3.py"
BOOTSTRAP_SHA256 = "0e3bbfc4f23508c8a5a03764c54c56d60a3e6e9c7f2eeec53bc9e1a9f78b614a"

PROJECTION_ID = "FGC-1-TDG8-RCV3-PROJECTION-1"
RECEIPT_SHA256 = "59e928830c213bd05d798bbc89b5994667d53a8b3c9ecbad3bf9023cf756dcba"
RECEIPT_RAW_SHA256 = "ed6b9c7ff72f18e8843ccfb0992d39dda7ee199e57cda8c17912d6d81004bd8a"

CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
INTERNAL_AUTHORIZATION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
ANCHOR_GENERATION = 9
ANCHOR_CHECKPOINT_SHA256 = "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
ANCHOR_CHECKPOINT_RAW_SHA256 = "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846"
ANCHOR_JOURNAL_SEQUENCE = 10
ANCHOR_JOURNAL_SHA256 = "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
ANCHOR_JOURNAL_RAW_SHA256 = "b50dc39ebd8bbdd3729d40d9e9ae4b223ff2ef19389d4ca74f9a316bac72202f"
TERMINAL_GENERATION = 10
TERMINAL_CHECKPOINT_SHA256 = "291eb466d5cb566b1f03822bf109c0c9e0f6ff930e9be496d26309b53b5005ea"
TERMINAL_CHECKPOINT_RAW_SHA256 = "5062b48d9e85fdb702fa3ed87e63fec78cb500baaeccf0b252d49acb83604054"
TERMINAL_JOURNAL_SEQUENCE = 11
TERMINAL_JOURNAL_SHA256 = "394c803df15f0b73827425dec4c914a1dd998aad8de23d8a8ecf3d637e1cf340"
TERMINAL_JOURNAL_RAW_SHA256 = "649dd7a21798b41196b99212c6ef86476d42065426058b21c6247aba80e0f5be"
TERMINAL_LOCK_SHA256 = "3bfb54c0e23d0940d52418402cdd24e8444a843398f3bbfb194d32eac05e8606"
TERMINAL_LOCK_RAW_SHA256 = "86eac601504de31f846bc54d958430c946d645e2a3ca93d7914aa2495b30ee0c"

SOURCE_LEAF_COUNT = 59
SOURCE_BYTE_COUNT = 6_342_040
SOURCE_TREE_SHA256 = "b076b5233e5e2ee9898444f6cade495e2e0ed8d24697348f640485104d2cf145"
PROJECTED_LEAF_COUNT = 50
PROJECTED_BYTE_COUNT = 6_201_938
PROJECTED_TREE_SHA256 = "d84e85c3dd29ac1554111ed00b37289d2b89273e7a7c2c47e3d21f5b9fda9e0b"
MEMBER_MAP_SHA256 = "4a7cf5daadbda2e3d085d355583645a5de21eff91cff89ad10aaddb1c2d3ff32"

STORE_DIRECTORIES = (
    "checkpoints", "journal", "locks", "payloads", "receipts", "states",
)
WRAPPER_DIRECTORIES = ("calibration",) + tuple(
    f"calibration/{name}" for name in STORE_DIRECTORIES
)
MEMBER_KEYS = (
    "RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-16385",
    "SSPRK3-4097", "SSPRK3-8193",
)
ARRAY_NAMES = (
    "u", "p", "q", "grid_coordinates", "tracer_labels",
    "tracer_positions", "tracer_proper_times", "event_proper_times",
    "event_fields",
)
EXCLUDED_PATHS = (
    f"checkpoints/{TERMINAL_GENERATION:020d}-{TERMINAL_CHECKPOINT_SHA256}.json",
    f"journal/{TERMINAL_JOURNAL_SEQUENCE:020d}-{TERMINAL_JOURNAL_SHA256}.journal",
    "locks/.active-write.lock.hlt16-quarantine-1d1ebcfbcb57ecf91adb191815acd57c1137609cd50f35fa98e7d18d2c91086c",
    "locks/.active-write.lock.hlt16-quarantine-2f76eeac85746761d38617804da02f0d520676a86a243a5e29b2a9f0fdfa71e9",
    "locks/.active-write.lock.hlt16-quarantine-6c409d46d6c049cced6a26d9ad733911070dad9d53e966674a7fdce4ba372946",
    "locks/.active-write.lock.hlt16-quarantine-f474863d36b1c0312ff677560f1385c2913add7353bde694f08e69161bcef5ec",
    "locks/bootstrap.guard", "locks/terminal.lock", "locks/writer.guard",
)

EXPECTED_AUTHORITY = {
    "authority_commit": AUTHORITY_COMMIT,
    "frz1_result_path": FRZ1_RESULT_PATH,
    "frz1_result_sha256": FRZ1_RESULT_SHA256,
    "runtime_path": RUNTIME_PATH,
    "runtime_sha256": RUNTIME_SHA256,
    "bootstrap_path": BOOTSTRAP_PATH,
    "bootstrap_sha256": BOOTSTRAP_SHA256,
}
EXPECTED_SOURCE = {
    "store_path": SOURCE_STORE,
    "leaf_count": SOURCE_LEAF_COUNT,
    "byte_count": SOURCE_BYTE_COUNT,
    "manifest_sha256": SOURCE_TREE_SHA256,
}
EXPECTED_SOURCE_TERMINAL = {
    "generation": TERMINAL_GENERATION,
    "checkpoint_sha256": TERMINAL_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": TERMINAL_CHECKPOINT_RAW_SHA256,
    "journal_sequence": TERMINAL_JOURNAL_SEQUENCE,
    "journal_sha256": TERMINAL_JOURNAL_SHA256,
    "journal_raw_sha256": TERMINAL_JOURNAL_RAW_SHA256,
    "terminal_lock_sha256": TERMINAL_LOCK_SHA256,
    "terminal_lock_raw_sha256": TERMINAL_LOCK_RAW_SHA256,
    "disposition": "invalid_terminal",
}
EXPECTED_PROJECTION = {
    "projection_id": PROJECTION_ID,
    "wrapper_path": DESTINATION_WRAPPER,
    "store_path": DESTINATION_STORE,
    "leaf_count": PROJECTED_LEAF_COUNT,
    "byte_count": PROJECTED_BYTE_COUNT,
    "manifest_sha256": PROJECTED_TREE_SHA256,
    "excluded_paths": list(EXCLUDED_PATHS),
}
EXPECTED_RECEIPT = {
    "path": RECEIPT_PATH,
    "receipt_sha256": RECEIPT_SHA256,
    "raw_sha256": RECEIPT_RAW_SHA256,
    "byte_count": 3070,
}
EXPECTED_ANCHOR = {
    "generation": ANCHOR_GENERATION,
    "checkpoint_sha256": ANCHOR_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": ANCHOR_CHECKPOINT_RAW_SHA256,
    "journal_sequence": ANCHOR_JOURNAL_SEQUENCE,
    "journal_sha256": ANCHOR_JOURNAL_SHA256,
    "journal_raw_sha256": ANCHOR_JOURNAL_RAW_SHA256,
    "campaign_id": CAMPAIGN_ID,
    "authorization_commit": INTERNAL_AUTHORIZATION_COMMIT,
    "plan_sha256": PLAN_SHA256,
    "event": 23,
    "target_rational": "3/2",
    "target_binary64_hex": "0x1.8000000000000p+0",
    "disposition": "nonterminal",
    "member_map_sha256": MEMBER_MAP_SHA256,
}
EXPECTED_MEMBERS = (
    {
        "member_key": "RK4-2049", "descriptor_sha256": "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4",
        "semantic_sha256": "0e48e0e2396ea7d9ef79d1641e6944a2dfdaa3716bb65e0b882b97156ff9372f",
        "raw_archive_sha256": "f47a4f533eb8e0ed9f8ea337d520fa6f4e0bb6c1b27a4ef8023ee78bc830f5a4",
        "cursor_sha256": "8f17e3efd76df310c97229672363609881c32cd063cd2386d6086387181a95c4",
        "ledger_sha256": "9f4003090de1fab36bf0c8aa833d6fdbc6c3cee84d431f2cf33a0bd248c6c3b2",
        "accepted_time_hex": "0x1.78554de5a30e0p+0", "mode": "RETRY_PENDING",
        "pending_owner": "temporal", "pending_cap_hex": "0x1.aaa9612df8000p-11",
        "retry_count": 2, "cumulative_retry_count": 2, "current_retry_count": 2,
        "physical_state_sha256": "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a",
        "u_bytes_sha256": "856fc327d3163d038f9fdbfc660968509dd7360faa0bb9d158721eb3049c2fe4",
        "p_bytes_sha256": "32939115b86938cf8a813bbf9a50a8c90bbaceaee161fa65695deb435d1243b7",
        "q_bytes_sha256": "48ec07bcace3935173c06e7b41915cacb17f14e1c9d3500522d4eb59a205c4a3",
    },
    {
        "member_key": "RK4-4097", "descriptor_sha256": "4f269d78a3f1fb25a2e28ddd14c0d6d95838b918d1bdae6c258339eede16c20c",
        "semantic_sha256": "febfa2b7a09885b7720308eeff2c1a2e69238cd391f7ddabcf396cf332dd8c30",
        "raw_archive_sha256": "588ccfe5b27910b369ff652f700c1aaf2000f38d3a72ddb1d47910a3eaca45bd",
        "cursor_sha256": "bc5cde83729395c5ac679ccb5e6fb087f2499b69b175a4c9840e4928e03cd34d",
        "ledger_sha256": "d8d46f7bdf3d8fc158aa6cd6bf5b387fcade13848ab2356c4bef491f0a6768ca",
        "accepted_time_hex": "0x1.7000000000000p+0", "mode": "FRESH_READY",
        "pending_owner": "", "pending_cap_hex": "", "retry_count": 0,
        "cumulative_retry_count": 0, "current_retry_count": 0,
        "physical_state_sha256": "9d27ca221fbac2a5b3b86a8d9aef49a8b8ac451b255b80529563b8cb7f8b73a9",
        "u_bytes_sha256": "0baf5adf0f80ba5f6f6af7c0dcf2590e7ae8f8580f443455db6e3032479fdc9c",
        "p_bytes_sha256": "aba9b66c1e30b0e20427aab4aeb5149a8a558e6607059e9bcb95f8972efd7e3f",
        "q_bytes_sha256": "df0f5590dec9ace9694ec9fb586fa0bbad0c1fef3b5a0545fa6a8ad5aecd3382",
    },
    {
        "member_key": "RK4-8193", "descriptor_sha256": "a5fbfae456a93318317a77cc6eb9ddc690e0552aa6374afeb5f3f9180a229d4b",
        "semantic_sha256": "37f65d160740815e102b9af121d0373fecf110a53043aeb95f59ff5f82c9c82f",
        "raw_archive_sha256": "21fcb06603ff42e725b91b26b40e87909a946c2d4d40f9c6cb914f2b265ecae0",
        "cursor_sha256": "ae4dafa295b44f0ac15c5b139858be9adba5d17485ac2a6d4ee1062fbea0725b",
        "ledger_sha256": "d8d46f7bdf3d8fc158aa6cd6bf5b387fcade13848ab2356c4bef491f0a6768ca",
        "accepted_time_hex": "0x1.7000000000000p+0", "mode": "FRESH_READY",
        "pending_owner": "", "pending_cap_hex": "", "retry_count": 0,
        "cumulative_retry_count": 0, "current_retry_count": 0,
        "physical_state_sha256": "d5b4cfe25199862fff5be7d524cc1c3a6030cc634d268cedb626c3ced7220f4a",
        "u_bytes_sha256": "44e26799d3e86ecf74ba94554223d1a9feecd3c36000f850be3148f7a7dd8386",
        "p_bytes_sha256": "64e454928b74d79070163cd2eede548d0991279ae5e58b7752b8472e786e16b0",
        "q_bytes_sha256": "19bab9efaa70ab8dbed51fe5039e9118c44605f206d07dc850a164b14e14a6b7",
    },
    {
        "member_key": "SSPRK3-16385", "descriptor_sha256": "69954e7e27672f35a9987657a602e69b2635b320d48061d180ea21e80b469234",
        "semantic_sha256": "b21e6d3da3785fa2d4be4c5869c5b05aba17bbe168c657d2cd2d70192c4d2205",
        "raw_archive_sha256": "72b49ffc72fc70185e1f611e273302dbf97b3ddc693fb54c3710d58dbecda5e9",
        "cursor_sha256": "9234ff83cca1f8d3e3a53634fb95428add21c1145888f7f1af9166c15a0cffa5",
        "ledger_sha256": "d8d46f7bdf3d8fc158aa6cd6bf5b387fcade13848ab2356c4bef491f0a6768ca",
        "accepted_time_hex": "0x1.7000000000000p+0", "mode": "FRESH_READY",
        "pending_owner": "", "pending_cap_hex": "", "retry_count": 0,
        "cumulative_retry_count": 0, "current_retry_count": 0,
        "physical_state_sha256": "91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d",
        "u_bytes_sha256": "d0ac42c0ad74a0a9d700731992a1504d7d9d067418d6bf5d7ebf351a91ed1ce6",
        "p_bytes_sha256": "df42a973169995b4edc181f3a8e2a28cdd0570b27522705198c4009ce36c1184",
        "q_bytes_sha256": "428f5f8b97a314c3bbdfcd4f94cf56fad781a78528ce8408b5f3e09c11bd4435",
    },
    {
        "member_key": "SSPRK3-4097", "descriptor_sha256": "886f39f2081145c59ebac9567719b5aed3ed6894acf2a7c797612259336fbc49",
        "semantic_sha256": "063ff123069409894b58529bcf7bd9282fcafff65a05791790e5b17b4dc26c5d",
        "raw_archive_sha256": "cea48b83168250b969d45f953001c559c2ef44b634e984b2952df36221bfb1d3",
        "cursor_sha256": "feb710c582a159fb614d9d588ad64ad1f7f3adf483ae817d643857bfc991cc71",
        "ledger_sha256": "d8d46f7bdf3d8fc158aa6cd6bf5b387fcade13848ab2356c4bef491f0a6768ca",
        "accepted_time_hex": "0x1.7000000000000p+0", "mode": "FRESH_READY",
        "pending_owner": "", "pending_cap_hex": "", "retry_count": 0,
        "cumulative_retry_count": 0, "current_retry_count": 0,
        "physical_state_sha256": "ab96602c0f05633375710efee2be01be84e24c32a31d07d6ef23390788a72162",
        "u_bytes_sha256": "3e3368a23f281ad62af81401f53fc6799d187d28d4f126b475793d25def925d1",
        "p_bytes_sha256": "7bb726d4a33a451ada46012c51928f0d97c47923bf19ad082f128d6ea1db9c97",
        "q_bytes_sha256": "62e994a54f617d37e58d4a93e76724d2dd576763a4a0fe0ad89e56e747756429",
    },
    {
        "member_key": "SSPRK3-8193", "descriptor_sha256": "9c5e1faf39bef93c3cefd8f83d47439851c4dada69d9b14ae1ed416dced7873a",
        "semantic_sha256": "eec392e02d2cf6b777e308885f4047faa7968377969499e8e4bd72d94238ea5f",
        "raw_archive_sha256": "8120b8787feaf356841f4d2318cb46a34e8fbbf15db8ac6122bc3fa78d8a9c0e",
        "cursor_sha256": "178b748157fad8802dcb498ae7d46cc07dcc1a8dd06fbcdf847fdb363b15af05",
        "ledger_sha256": "d8d46f7bdf3d8fc158aa6cd6bf5b387fcade13848ab2356c4bef491f0a6768ca",
        "accepted_time_hex": "0x1.7000000000000p+0", "mode": "FRESH_READY",
        "pending_owner": "", "pending_cap_hex": "", "retry_count": 0,
        "cumulative_retry_count": 0, "current_retry_count": 0,
        "physical_state_sha256": "c46abc3f0f72fa2a13c1fee7bb329b94698cc2c26003406ae3eb1592f11b648c",
        "u_bytes_sha256": "ed56d358550c22b3817c17f4d358220f8dce7006037705cb533c154de92c2855",
        "p_bytes_sha256": "d9ceaefbfa6df6e840085461d65eef1194007d0e138c59d172d7d8a8ee75e923",
        "q_bytes_sha256": "1f82429ba9321a4777f219a169d61ce8ddf0d1fd7b3f52b9c8fac0469b00b8c5",
    },
)
EXPECTED_SCOPE = {
    "source_store_read_only": True,
    "destination_store_read_only": True,
    "bootstrap_reexecuted": False,
    "PDE_proposal_execution": False,
    "accepted_state_mutation": False,
    "bounded_event_execution_authorized": False,
    "candidate_branch_opened": False,
}
EXPECTED_CLAIMS = {
    "source_terminal_store_unchanged": True,
    "installed_projection_authenticated": True,
    "generation9_restart_boundary_authenticated": True,
    "six_member_payloads_decoded_and_hash_bound": True,
    "retry_pending_depth2_and_half_cap_preserved": True,
    "source_or_destination_mutated": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}
EXPECTED_NONCLAIMS = [
    "PREF1 authenticates a copied finite-state prefix; it does not authorize or execute the next PDE proposal.",
    "The internal HLT16 campaign identity is intentionally unchanged; the external projection ID and committed receipt distinguish the recovery fork.",
    "A restartable generation-nine boundary is not a completed common event, GR-0 calibration, candidate trajectory, activation, trapping, DEF1 result, retained-EFT result, transition, or physical mechanism.",
]

_CHECKPOINT_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.json\Z")
_JOURNAL_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.journal\Z")
_STATE_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.json\Z")
_PAYLOAD_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.npz\Z")
_MAX_LEAF_BYTES = 64 * 1024 * 1024
_MAX_TREE_BYTES = 128 * 1024 * 1024
_MAX_ARRAY_BYTES = 48 * 1024 * 1024


class TDG8RCV3PREF1Error(ValueError):
    """The compact contract or either authenticated tree differs."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


def _stop(stop_id: str, detail: str) -> None:
    raise TDG8RCV3PREF1Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        _stop("PREF1_CANONICAL_DRIFT", "value is not canonical JSON")
        raise AssertionError from error


def canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        _stop("PREF1_COMPACT_DRIFT", "result is not canonical JSON")
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
        _stop("PREF1_JSON_DRIFT", f"{label} is malformed")
        raise AssertionError from error
    expected = canonical_result(value) if pretty else _canonical(value)
    if not isinstance(value, dict) or raw != expected:
        _stop("PREF1_JSON_DRIFT", f"{label} is noncanonical")
    return value


def _object_address(value: Mapping[str, Any], field: str, address: str, label: str) -> None:
    if value.get(field) != address:
        _stop("PREF1_HASH_DRIFT", f"{label} embedded address differs")
    body = dict(value)
    body.pop(field, None)
    if _digest(body) != address:
        _stop("PREF1_HASH_DRIFT", f"{label} content address differs")


def _sha(value: object, label: str, length: int = 64) -> str:
    if not isinstance(value, str) or len(value) != length or value.lower() != value:
        _stop("PREF1_HASH_DRIFT", f"{label} is not a lowercase digest")
    try:
        int(value, 16)
    except ValueError as error:
        _stop("PREF1_HASH_DRIFT", f"{label} is not hexadecimal")
        raise AssertionError from error
    return value


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        _stop("PREF1_CONFIG_DRIFT", "PREF1 config is malformed")
        raise AssertionError from error
    expected_keys = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "authority", "source", "source_terminal",
        "projection", "receipt", "anchor", "members", "validation", "scope", "claims",
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        _stop("PREF1_CONFIG_DRIFT", "PREF1 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != EXPECTED_NONCLAIMS
        or value["authority"] != EXPECTED_AUTHORITY
        or value["source"] != EXPECTED_SOURCE
        or value["source_terminal"] != EXPECTED_SOURCE_TERMINAL
        or value["projection"] != EXPECTED_PROJECTION
        or value["receipt"] != EXPECTED_RECEIPT
        or value["anchor"] != EXPECTED_ANCHOR
        or tuple(value["members"]) != EXPECTED_MEMBERS
        or value["validation"] != {
            "checkpoint_count": 10, "journal_count": 11,
            "legacy_state_count": 6, "runtime_descriptor_count": 11,
            "payload_count": 11, "array_count_per_payload": 9,
            "destination_lock_leaf_count": 0,
        }
        or value["scope"] != EXPECTED_SCOPE
        or value["claims"] != EXPECTED_CLAIMS
    ):
        _stop("PREF1_CONFIG_DRIFT", "PREF1 exact contract differs")
    return value


@dataclass(frozen=True, slots=True)
class _Leaf:
    path: str
    byte_count: int
    raw_sha256: str
    raw: bytes

    def manifest_row(self) -> dict[str, object]:
        return {"path": self.path, "byte_count": self.byte_count, "sha256": self.raw_sha256}


@dataclass(frozen=True, slots=True)
class _Snapshot:
    leaves: Mapping[str, _Leaf]
    directories: tuple[str, ...]

    def rows(self) -> list[dict[str, object]]:
        return [self.leaves[path].manifest_row() for path in sorted(self.leaves)]

    def manifest_sha256(self) -> str:
        return sha256(_canonical(self.rows())).hexdigest()

    def byte_count(self) -> int:
        return sum(leaf.byte_count for leaf in self.leaves.values())


def _directory_flags() -> int:
    return os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)


def _open_chain(root: Path, parts: Sequence[str], label: str) -> int:
    try:
        before = root.lstat()
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
            _stop("PREF1_PATH_UNSAFE", f"{label} root is unsafe")
        current = os.open(root, _directory_flags())
        active = os.fstat(current)
        if (before.st_dev, before.st_ino) != (active.st_dev, active.st_ino):
            os.close(current)
            _stop("PREF1_PATH_UNSAFE", f"{label} root raced")
    except OSError as error:
        _stop("PREF1_PATH_UNSAFE", f"{label} root cannot be opened")
        raise AssertionError from error
    try:
        for component in parts:
            if component in {"", ".", ".."} or "/" in component or "\\" in component:
                _stop("PREF1_PATH_UNSAFE", f"{label} component is unsafe")
            item = os.stat(component, dir_fd=current, follow_symlinks=False)
            if stat.S_ISLNK(item.st_mode) or not stat.S_ISDIR(item.st_mode):
                _stop("PREF1_PATH_UNSAFE", f"{label} ancestor is unsafe")
            child = os.open(component, _directory_flags(), dir_fd=current)
            active = os.fstat(child)
            if (item.st_dev, item.st_ino) != (active.st_dev, active.st_ino):
                os.close(child)
                _stop("PREF1_PATH_UNSAFE", f"{label} ancestor raced")
            os.close(current)
            current = child
        return current
    except BaseException:
        os.close(current)
        raise


def _read_leaf(parent: int, name: str, before: os.stat_result, label: str) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1 or before.st_size < 0 or before.st_size > _MAX_LEAF_BYTES
    ):
        _stop("PREF1_PATH_UNSAFE", f"{label} is not a bounded single-linked file")
    descriptor = -1
    try:
        descriptor = os.open(name, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent)
        active = os.fstat(descriptor)
        identity = (
            before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
            before.st_ctime_ns, before.st_nlink,
        )
        if identity != (
            active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns,
            active.st_ctime_ns, active.st_nlink,
        ):
            _stop("PREF1_PATH_UNSAFE", f"{label} raced before read")
        remaining = before.st_size
        chunks: list[bytes] = []
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("PREF1_PATH_UNSAFE", f"{label} ended early")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF1_PATH_UNSAFE", f"{label} grew during read")
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        path_after = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if identity != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
            after.st_ctime_ns, after.st_nlink,
        ) or (path_after.st_dev, path_after.st_ino) != (before.st_dev, before.st_ino):
            _stop("PREF1_PATH_UNSAFE", f"{label} changed during read")
        return raw
    except OSError as error:
        _stop("PREF1_PATH_UNSAFE", f"{label} cannot be read safely")
        raise AssertionError from error
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _snapshot_once(repository: Path, relative: str) -> _Snapshot:
    path = Path(relative)
    if path.is_absolute() or path.as_posix() != relative or any(part in {"", ".", ".."} for part in path.parts):
        _stop("PREF1_PATH_UNSAFE", "tree path is unsafe")
    root = _open_chain(repository, path.parts, f"tree {relative}")
    leaves: dict[str, _Leaf] = {}
    directories: list[str] = []
    total = 0

    def descend(directory: int, prefix: str) -> None:
        nonlocal total
        before_directory = os.fstat(directory)
        try:
            with os.scandir(directory) as scan:
                entries = sorted(scan, key=lambda item: item.name)
        except OSError as error:
            _stop("PREF1_PATH_UNSAFE", f"tree cannot be scanned: {prefix or '.'}")
            raise AssertionError from error
        for entry in entries:
            name = entry.name
            if name in {"", ".", ".."} or "/" in name or "\\" in name:
                _stop("PREF1_PATH_UNSAFE", "tree name is unsafe")
            child_path = f"{prefix}/{name}" if prefix else name
            item = os.stat(name, dir_fd=directory, follow_symlinks=False)
            if stat.S_ISLNK(item.st_mode):
                _stop("PREF1_PATH_UNSAFE", f"symlink appears: {child_path}")
            if stat.S_ISREG(item.st_mode):
                raw = _read_leaf(directory, name, item, child_path)
                total += len(raw)
                if total > _MAX_TREE_BYTES:
                    _stop("PREF1_PATH_UNSAFE", "tree exceeds its byte budget")
                leaves[child_path] = _Leaf(child_path, len(raw), sha256(raw).hexdigest(), raw)
                continue
            if not stat.S_ISDIR(item.st_mode):
                _stop("PREF1_PATH_UNSAFE", f"special node appears: {child_path}")
            child = os.open(name, _directory_flags(), dir_fd=directory)
            try:
                active = os.fstat(child)
                if (item.st_dev, item.st_ino) != (active.st_dev, active.st_ino):
                    _stop("PREF1_PATH_UNSAFE", f"directory raced: {child_path}")
                directories.append(child_path)
                descend(child, child_path)
            finally:
                os.close(child)
        after_directory = os.fstat(directory)
        if (
            before_directory.st_dev, before_directory.st_ino,
            before_directory.st_mtime_ns, before_directory.st_ctime_ns,
        ) != (
            after_directory.st_dev, after_directory.st_ino,
            after_directory.st_mtime_ns, after_directory.st_ctime_ns,
        ):
            _stop("PREF1_PATH_UNSAFE", f"directory changed: {prefix or '.'}")

    try:
        descend(root, "")
    finally:
        os.close(root)
    return _Snapshot(dict(sorted(leaves.items())), tuple(sorted(directories)))


def _snapshot(repository: Path, relative: str) -> _Snapshot:
    first = _snapshot_once(repository, relative)
    second = _snapshot_once(repository, relative)
    if first != second:
        _stop("PREF1_PATH_UNSAFE", f"tree changed during authentication: {relative}")
    return first


def _manifest(snapshot: _Snapshot, expected: Mapping[str, Any], label: str) -> None:
    if (
        len(snapshot.leaves) != expected["leaf_count"]
        or snapshot.byte_count() != expected["byte_count"]
        or snapshot.manifest_sha256() != expected["manifest_sha256"]
    ):
        _stop("PREF1_TREE_DRIFT", f"{label} manifest differs")


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
        _stop("PREF1_AUTHORITY_DRIFT", "Git authority lookup failed")
        raise AssertionError from error


def _bind_authority(repository: Path) -> None:
    resolved = _git(repository, "rev-parse", "--verify", f"{AUTHORITY_COMMIT}^{{commit}}")
    if resolved != f"{AUTHORITY_COMMIT}\n".encode("ascii"):
        _stop("PREF1_AUTHORITY_DRIFT", "installation authority commit differs")
    _git(repository, "merge-base", "--is-ancestor", AUTHORITY_COMMIT, "HEAD")
    for relative, expected_hash in (
        (FRZ1_RESULT_PATH, FRZ1_RESULT_SHA256),
        (RUNTIME_PATH, RUNTIME_SHA256),
        (BOOTSTRAP_PATH, BOOTSTRAP_SHA256),
    ):
        tree = _git(repository, "ls-tree", "-z", AUTHORITY_COMMIT, "--", relative)
        entries = [entry for entry in tree.split(b"\0") if entry]
        if len(entries) != 1 or not entries[0].startswith(b"100644 blob "):
            _stop("PREF1_AUTHORITY_DRIFT", f"authority blob differs: {relative}")
        raw = _git(repository, "show", f"{AUTHORITY_COMMIT}:{relative}")
        if sha256(raw).hexdigest() != expected_hash:
            _stop("PREF1_AUTHORITY_DRIFT", f"authority hash differs: {relative}")
        if relative == FRZ1_RESULT_PATH:
            result = _json(raw, "FRZ1 result", pretty=True)
            payload = result.get("artifact_payload")
            if result.get("artifact_id") != "FGC-1-TDG8-RCV3-FRZ1" or not isinstance(payload, Mapping):
                _stop("PREF1_AUTHORITY_DRIFT", "FRZ1 result identity differs")
            if payload.get("claims", {}).get("projection_authenticated") is not False:
                _stop("PREF1_AUTHORITY_DRIFT", "FRZ1 already promoted the projection")


def _validate_tree_grammar(snapshot: _Snapshot, *, source: bool) -> None:
    if snapshot.directories != STORE_DIRECTORIES:
        _stop("PREF1_TREE_DRIFT", "store directory grammar differs")
    for relative, leaf in snapshot.leaves.items():
        parts = Path(relative).parts
        if len(parts) != 2 or parts[0] not in STORE_DIRECTORIES:
            _stop("PREF1_TREE_DRIFT", f"nested or foreign leaf appears: {relative}")
        directory, name = parts
        valid = False
        if directory == "checkpoints":
            valid = _CHECKPOINT_NAME.fullmatch(name) is not None
        elif directory == "journal":
            valid = _JOURNAL_NAME.fullmatch(name) is not None
        elif directory == "states":
            valid = _STATE_NAME.fullmatch(name) is not None
        elif directory == "payloads":
            valid = _PAYLOAD_NAME.fullmatch(name) is not None
        elif directory == "receipts":
            valid = name == "generation-zero.json"
        elif directory == "locks" and source:
            valid = relative in EXCLUDED_PATHS
            if name in {"bootstrap.guard", "writer.guard"} and leaf.raw != b"":
                _stop("PREF1_TREE_DRIFT", f"source guard has bytes: {name}")
        if not valid:
            _stop("PREF1_TREE_DRIFT", f"store leaf grammar differs: {relative}")
    if not source and any(path.startswith("locks/") for path in snapshot.leaves):
        _stop("PREF1_TREE_DRIFT", "projected lock directory is not empty")


def _indexed(
    leaves: Mapping[str, _Leaf], directory: str, pattern: re.Pattern[str], field: str,
) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for relative in sorted(path for path in leaves if path.startswith(f"{directory}/")):
        match = pattern.fullmatch(Path(relative).name)
        if match is None:
            _stop("PREF1_TREE_DRIFT", f"{directory} filename differs")
        number = int(match.group("number"))
        address = match.group("digest")
        if number in result:
            _stop("PREF1_TREE_DRIFT", f"{directory} fork appears at {number}")
        value = _json(leaves[relative].raw, f"{directory} {number}")
        _object_address(value, field, address, f"{directory} {number}")
        result[number] = value
    return result


def _validate_lineage(
    leaves: Mapping[str, _Leaf],
) -> tuple[dict[int, dict[str, Any]], dict[int, dict[str, Any]]]:
    checkpoints = _indexed(leaves, "checkpoints", _CHECKPOINT_NAME, "checkpoint_sha256")
    journals = _indexed(leaves, "journal", _JOURNAL_NAME, "record_sha256")
    if set(checkpoints) != set(range(10)) or set(journals) != set(range(11)):
        _stop("PREF1_LINEAGE_DRIFT", "checkpoint or journal chain has a gap/fork")
    genesis = checkpoints[0]
    if (
        genesis.get("campaign_generation") != 0
        or genesis.get("parent_checkpoint_sha256") != "0" * 64
        or genesis.get("terminal_lock") is not False
        or genesis.get("disposition") != "nonterminal"
    ):
        _stop("PREF1_LINEAGE_DRIFT", "generation-zero checkpoint differs")
    parent_checkpoint = str(genesis["checkpoint_sha256"])
    for generation in range(1, 10):
        checkpoint = checkpoints[generation]
        if (
            checkpoint.get("generation") != generation
            or checkpoint.get("parent_sha256") != parent_checkpoint
            or checkpoint.get("authorization_commit") != INTERNAL_AUTHORIZATION_COMMIT
            or checkpoint.get("plan_sha256") != PLAN_SHA256
            or checkpoint.get("protocol") != TARGET_PROTOCOL
            or checkpoint.get("campaign_id") != CAMPAIGN_ID
            or checkpoint.get("event") != 23
            or checkpoint.get("target") != {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"}
            or checkpoint.get("disposition") != "nonterminal"
            or checkpoint.get("terminal") is not None
            or not isinstance(checkpoint.get("members"), Mapping)
            or set(checkpoint["members"]) != set(MEMBER_KEYS)
        ):
            _stop("PREF1_LINEAGE_DRIFT", f"generation {generation} relation differs")
        sequence = checkpoint.get("journal_sequence")
        if not isinstance(sequence, int) or sequence not in journals or checkpoint.get("journal_tip_sha256") != journals[sequence].get("record_sha256"):
            _stop("PREF1_LINEAGE_DRIFT", f"generation {generation} journal tip differs")
        for member_key, member in checkpoint["members"].items():
            if not isinstance(member, Mapping) or not isinstance(member.get("cursor"), Mapping) or not isinstance(member.get("ledger"), Mapping):
                _stop("PREF1_LINEAGE_DRIFT", f"generation {generation} {member_key} differs")
            cursor = member["cursor"]
            _object_address(cursor, "cursor_chain_sha256", str(cursor.get("cursor_chain_sha256")), f"{member_key} cursor")
            if (
                member.get("descriptor_sha256") != cursor.get("accepted_state_sha256")
                or _digest(member["ledger"]) != cursor.get("TDG6_ledger_sha256")
                or member["ledger"].get("last_accepted_time_hex") != cursor.get("accepted_boundary_time", {}).get("binary64_hex")
            ):
                _stop("PREF1_LINEAGE_DRIFT", f"generation {generation} {member_key} state/ledger differs")
        parent_checkpoint = str(checkpoint["checkpoint_sha256"])
    parent_record = "0" * 64
    for sequence in range(11):
        record = journals[sequence]
        if (
            record.get("schema") != "FGC-1-HLT16-campaign-journal-v1"
            or record.get("sequence") != sequence
            or record.get("campaign_id") != CAMPAIGN_ID
            or record.get("event") != 23
            or record.get("previous_record_sha256") != parent_record
        ):
            _stop("PREF1_LINEAGE_DRIFT", f"journal {sequence} relation differs")
        parent_record = str(record["record_sha256"])
    anchor_path = f"checkpoints/{ANCHOR_GENERATION:020d}-{ANCHOR_CHECKPOINT_SHA256}.json"
    journal_path = f"journal/{ANCHOR_JOURNAL_SEQUENCE:020d}-{ANCHOR_JOURNAL_SHA256}.journal"
    if (
        leaves[anchor_path].raw_sha256 != ANCHOR_CHECKPOINT_RAW_SHA256
        or leaves[journal_path].raw_sha256 != ANCHOR_JOURNAL_RAW_SHA256
    ):
        _stop("PREF1_HASH_DRIFT", "generation-nine raw anchor differs")
    anchor = checkpoints[ANCHOR_GENERATION]
    if (
        anchor.get("checkpoint_sha256") != ANCHOR_CHECKPOINT_SHA256
        or anchor.get("journal_sequence") != ANCHOR_JOURNAL_SEQUENCE
        or anchor.get("journal_tip_sha256") != ANCHOR_JOURNAL_SHA256
        or _digest(anchor.get("members")) != MEMBER_MAP_SHA256
    ):
        _stop("PREF1_LINEAGE_DRIFT", "generation-nine anchor differs")
    return checkpoints, journals


def _descriptor_objects(
    leaves: Mapping[str, _Leaf], checkpoints: Mapping[int, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    observed = {
        match.group("digest") for relative in leaves if relative.startswith("states/")
        for match in [_STATE_NAME.fullmatch(Path(relative).name)] if match is not None
    }
    genesis_states = checkpoints[0].get("states")
    if not isinstance(genesis_states, Mapping) or set(genesis_states) != set(MEMBER_KEYS):
        _stop("PREF1_DESCRIPTOR_DRIFT", "generation-zero state set differs")
    legacy = {_sha(value, "legacy state") for value in genesis_states.values()}
    runtime = {
        _sha(checkpoints[generation]["members"][key]["descriptor_sha256"], "runtime descriptor")
        for generation in range(1, 10) for key in MEMBER_KEYS
    }
    if observed != legacy | runtime or legacy & runtime or len(legacy) != 6 or len(runtime) != 11:
        _stop("PREF1_DESCRIPTOR_DRIFT", "state descriptor orphan/omission appears")
    for address in legacy:
        leaf = leaves[f"states/{address}.json"]
        _json(leaf.raw, f"legacy state {address}")
        if leaf.raw_sha256 != address:
            _stop("PREF1_DESCRIPTOR_DRIFT", "legacy state raw address differs")
    descriptors: dict[str, dict[str, Any]] = {}
    for address in runtime:
        value = _json(leaves[f"states/{address}.json"].raw, f"runtime descriptor {address}")
        if set(value) != {"schema", "semantic_sha256", "raw_archive_sha256", "arrays", "metadata", "descriptor_sha256"} or value.get("schema") != "FGC-1-HLT16-member-state-v1":
            _stop("PREF1_DESCRIPTOR_DRIFT", "runtime descriptor schema differs")
        _object_address(value, "descriptor_sha256", address, "runtime descriptor")
        _sha(value.get("semantic_sha256"), "semantic payload")
        _sha(value.get("raw_archive_sha256"), "raw archive")
        descriptors[address] = value
    return descriptors


def _decode_npz(raw: bytes, descriptor: Mapping[str, Any], label: str) -> dict[str, np.ndarray]:
    try:
        archive = zipfile.ZipFile(BytesIO(raw), "r")
        infos = archive.infolist()
    except (OSError, zipfile.BadZipFile) as error:
        _stop("PREF1_PAYLOAD_DRIFT", f"{label} ZIP is malformed")
        raise AssertionError from error
    arrays: dict[str, np.ndarray] = {}
    try:
        expected_names = tuple(f"{name}.npy" for name in ARRAY_NAMES)
        if tuple(info.filename for info in infos) != expected_names or len({info.filename for info in infos}) != len(infos):
            _stop("PREF1_PAYLOAD_DRIFT", f"{label} ZIP inventory/order differs")
        total = 0
        for info in infos:
            mode = (info.external_attr >> 16) & 0o170000
            if (
                info.is_dir() or info.flag_bits & 0x1 or info.file_size <= 0
                or info.file_size > _MAX_ARRAY_BYTES + (1 << 20)
                or info.compress_size < 0 or mode not in {0, stat.S_IFREG}
                or Path(info.filename).name != info.filename
            ):
                _stop("PREF1_PAYLOAD_DRIFT", f"{label} ZIP member is unsafe")
            total += info.file_size
        if total > _MAX_TREE_BYTES:
            _stop("PREF1_PAYLOAD_DRIFT", f"{label} ZIP expansion exceeds budget")
        for name, info in zip(ARRAY_NAMES, infos, strict=True):
            encoded = archive.read(info)
            stream = BytesIO(encoded)
            try:
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    shape, fortran_order, dtype = np.lib.format.read_array_header_1_0(
                        stream,
                    )
                elif version == (2, 0):
                    shape, fortran_order, dtype = np.lib.format.read_array_header_2_0(
                        stream,
                    )
                else:
                    _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY version differs")
                element_count = 1
                for extent in shape:
                    if not isinstance(extent, int) or extent < 0:
                        _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY shape is unsafe")
                    element_count *= extent
                    if element_count * 8 > _MAX_ARRAY_BYTES:
                        _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY shape exceeds budget")
                expected_bytes = element_count * dtype.itemsize
                if (
                    not shape or fortran_order or dtype.str != "<f8"
                    or expected_bytes > _MAX_ARRAY_BYTES
                    or len(encoded) - stream.tell() != expected_bytes
                ):
                    _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY header differs")
                stream.seek(0)
                value = np.load(stream, allow_pickle=False)
            except (OSError, ValueError, EOFError) as error:
                _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY member is malformed")
                raise AssertionError from error
            if stream.tell() != len(encoded):
                _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY member has trailing bytes")
            if (
                not isinstance(value, np.ndarray) or value.dtype.str != "<f8"
                or not value.flags.c_contiguous or not value.shape
                or value.nbytes > _MAX_ARRAY_BYTES or not np.isfinite(value).all()
            ):
                _stop("PREF1_PAYLOAD_DRIFT", f"{label} NPY dtype/layout/value differs")
            arrays[name] = value.copy(order="C")
    except (RuntimeError, zipfile.BadZipFile, KeyError) as error:
        _stop("PREF1_PAYLOAD_DRIFT", f"{label} ZIP member cannot be decoded")
        raise AssertionError from error
    finally:
        archive.close()
    manifest = [
        {
            "name": name, "dtype": "<f8", "shape": list(arrays[name].shape),
            "order": "C", "byte_count": int(arrays[name].nbytes),
            "bytes_sha256": sha256(arrays[name].tobytes(order="C")).hexdigest(),
        }
        for name in ARRAY_NAMES
    ]
    if manifest != descriptor.get("arrays") or _digest(manifest) != descriptor.get("semantic_sha256"):
        _stop("PREF1_PAYLOAD_DRIFT", f"{label} manifest differs")
    return arrays


def _load_payloads(
    leaves: Mapping[str, _Leaf], descriptors: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, np.ndarray]]:
    expected = {str(value["semantic_sha256"]) for value in descriptors.values()}
    observed = {
        match.group("digest") for relative in leaves if relative.startswith("payloads/")
        for match in [_PAYLOAD_NAME.fullmatch(Path(relative).name)] if match is not None
    }
    if observed != expected or len(observed) != 11:
        _stop("PREF1_PAYLOAD_DRIFT", "payload orphan/omission appears")
    payloads: dict[str, dict[str, np.ndarray]] = {}
    expected_raw: dict[str, str] = {}
    for descriptor in descriptors.values():
        semantic = str(descriptor["semantic_sha256"])
        raw_hash = str(descriptor["raw_archive_sha256"])
        if semantic in expected_raw and expected_raw[semantic] != raw_hash:
            _stop("PREF1_PAYLOAD_DRIFT", "shared semantic payload conflicts")
        expected_raw[semantic] = raw_hash
        if semantic in payloads:
            continue
        leaf = leaves[f"payloads/{semantic}.npz"]
        if leaf.raw_sha256 != raw_hash:
            _stop("PREF1_PAYLOAD_DRIFT", f"payload raw hash differs: {semantic}")
        payloads[semantic] = _decode_npz(leaf.raw, descriptor, f"payload {semantic}")
    return payloads


def _physical_state_sha256(*arrays: np.ndarray) -> str:
    digestor = sha256()
    for value in arrays:
        normalized = np.asarray(value, dtype="<f8", order="C")
        digestor.update(struct.pack("<I", normalized.ndim))
        digestor.update(struct.pack(f"<{normalized.ndim}Q", *normalized.shape))
        digestor.update(normalized.tobytes(order="C"))
    return digestor.hexdigest()


def _member_evidence(
    checkpoint: Mapping[str, Any], descriptors: Mapping[str, Mapping[str, Any]],
    payloads: Mapping[str, Mapping[str, np.ndarray]],
) -> list[dict[str, Any]]:
    observed: list[dict[str, Any]] = []
    for expected in EXPECTED_MEMBERS:
        key = str(expected["member_key"])
        member = checkpoint["members"].get(key)
        if not isinstance(member, Mapping):
            _stop("PREF1_MEMBER_DRIFT", f"member is absent: {key}")
        cursor = member.get("cursor")
        ledger = member.get("ledger")
        if not isinstance(cursor, Mapping) or not isinstance(ledger, Mapping):
            _stop("PREF1_MEMBER_DRIFT", f"member cursor/ledger differs: {key}")
        descriptor_id = str(member.get("descriptor_sha256"))
        descriptor = descriptors.get(descriptor_id)
        if descriptor is None:
            _stop("PREF1_MEMBER_DRIFT", f"member descriptor differs: {key}")
        semantic = str(descriptor["semantic_sha256"])
        arrays = payloads[semantic]
        retry = cursor.get("retry_successor_payload_or_none")
        retry_count = 0
        if retry is not None:
            if not isinstance(retry, Mapping):
                _stop("PREF1_RETRY_DRIFT", f"retry payload differs: {key}")
            retry_count = retry.get("retry_count")
        evidence = {
            "member_key": key,
            "descriptor_sha256": descriptor_id,
            "semantic_sha256": semantic,
            "raw_archive_sha256": descriptor["raw_archive_sha256"],
            "cursor_sha256": cursor.get("cursor_chain_sha256"),
            "ledger_sha256": _digest(ledger),
            "accepted_time_hex": cursor.get("accepted_boundary_time", {}).get("binary64_hex"),
            "mode": cursor.get("mode"),
            "pending_owner": member.get("pending_owner") or "",
            "pending_cap_hex": member.get("pending_cap_hex") or "",
            "retry_count": retry_count,
            "cumulative_retry_count": ledger.get("cumulative_temporal_retry_count"),
            "current_retry_count": ledger.get("current_macro_step_temporal_retry_count"),
            "physical_state_sha256": _physical_state_sha256(arrays["u"], arrays["p"], arrays["q"]),
            "u_bytes_sha256": sha256(arrays["u"].tobytes(order="C")).hexdigest(),
            "p_bytes_sha256": sha256(arrays["p"].tobytes(order="C")).hexdigest(),
            "q_bytes_sha256": sha256(arrays["q"].tobytes(order="C")).hexdigest(),
        }
        if evidence != expected:
            _stop("PREF1_MEMBER_DRIFT", f"generation-nine member differs: {key}")
        if key == "RK4-2049":
            prefix = retry.get("complete_rejection_prefix") if isinstance(retry, Mapping) else None
            if (
                cursor.get("mode") != "RETRY_PENDING"
                or member.get("pending_owner") != "temporal"
                or member.get("pending_cap_hex") != "0x1.aaa9612df8000p-11"
                or retry.get("half_cap") != "0x1.aaa9612df8000p-11"
                or retry.get("retry_count") != 2
                or not isinstance(prefix, list) or len(prefix) != 2
                or retry.get("accepted_state_sha256") != descriptor_id
                or retry.get("accepted_boundary_time") != cursor.get("accepted_boundary_time")
                or retry.get("TDG6_rejection_evidence", {}).get("initial_state_sha256") != evidence["physical_state_sha256"]
            ):
                _stop("PREF1_RETRY_DRIFT", "depth-two persisted retry differs")
        elif retry is not None or cursor.get("mode") != "FRESH_READY":
            _stop("PREF1_RETRY_DRIFT", f"unexpected retry state: {key}")
        observed.append(evidence)
    return observed


def _validate_source_terminal(snapshot: _Snapshot) -> None:
    checkpoint_path = f"checkpoints/{TERMINAL_GENERATION:020d}-{TERMINAL_CHECKPOINT_SHA256}.json"
    journal_path = f"journal/{TERMINAL_JOURNAL_SEQUENCE:020d}-{TERMINAL_JOURNAL_SHA256}.journal"
    lock_path = "locks/terminal.lock"
    if (
        snapshot.leaves[checkpoint_path].raw_sha256 != TERMINAL_CHECKPOINT_RAW_SHA256
        or snapshot.leaves[journal_path].raw_sha256 != TERMINAL_JOURNAL_RAW_SHA256
        or snapshot.leaves[lock_path].raw_sha256 != TERMINAL_LOCK_RAW_SHA256
    ):
        _stop("PREF1_SOURCE_DRIFT", "source terminal raw bytes differ")
    checkpoint = _json(snapshot.leaves[checkpoint_path].raw, "source terminal checkpoint")
    journal = _json(snapshot.leaves[journal_path].raw, "source terminal journal")
    lock = _json(snapshot.leaves[lock_path].raw, "source terminal lock")
    _object_address(checkpoint, "checkpoint_sha256", TERMINAL_CHECKPOINT_SHA256, "source terminal checkpoint")
    _object_address(journal, "record_sha256", TERMINAL_JOURNAL_SHA256, "source terminal journal")
    _object_address(lock, "lock_sha256", TERMINAL_LOCK_SHA256, "source terminal lock")
    if (
        checkpoint.get("generation") != TERMINAL_GENERATION
        or checkpoint.get("journal_sequence") != TERMINAL_JOURNAL_SEQUENCE
        or checkpoint.get("journal_tip_sha256") != TERMINAL_JOURNAL_SHA256
        or checkpoint.get("disposition") != "invalid_terminal"
        or journal.get("sequence") != TERMINAL_JOURNAL_SEQUENCE
        or journal.get("previous_record_sha256") != ANCHOR_JOURNAL_SHA256
        or journal.get("kind") != "terminal_lock"
        or lock.get("checkpoint_sha256") != TERMINAL_CHECKPOINT_SHA256
        or lock.get("journal_tip_sha256") != TERMINAL_JOURNAL_SHA256
    ):
        _stop("PREF1_SOURCE_DRIFT", "source terminal relation differs")


def _expected_receipt() -> dict[str, Any]:
    body: dict[str, Any] = {
        "schema": "FGC-1-TDG8-RCV3-external-recovery-fork-v1",
        "projection_id": PROJECTION_ID,
        "kind": "externally_identified_recovery_fork_projection",
        "source_terminal_store": SOURCE_STORE,
        "wrapper_root": DESTINATION_WRAPPER,
        "projected_store": DESTINATION_STORE,
        "external_install_authority": {
            "authority_commit_sha": AUTHORITY_COMMIT,
            "frz1_result_sha256": FRZ1_RESULT_SHA256,
            "runtime_sha256": RUNTIME_SHA256,
            "bootstrap_script_sha256": BOOTSTRAP_SHA256,
        },
        "source_terminal": {
            "generation": TERMINAL_GENERATION,
            "checkpoint_sha256": TERMINAL_CHECKPOINT_SHA256,
            "journal_sequence": TERMINAL_JOURNAL_SEQUENCE,
            "journal_sha256": TERMINAL_JOURNAL_SHA256,
            "disposition": "invalid_terminal",
            "terminal_lock_authenticated": True,
            "writer_absent": True,
            "suffix_absent": True,
        },
        "selected_lineage_anchor": {
            "generation": ANCHOR_GENERATION,
            "checkpoint_sha256": ANCHOR_CHECKPOINT_SHA256,
            "journal_sequence": ANCHOR_JOURNAL_SEQUENCE,
            "journal_sha256": ANCHOR_JOURNAL_SHA256,
            "disposition": "nonterminal", "event": 23,
            "target": {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"},
        },
        "internal_hlt16_identity": {
            "campaign_id": CAMPAIGN_ID,
            "authorization_commit": INTERNAL_AUTHORIZATION_COMMIT,
            "plan_sha256": PLAN_SHA256,
            "byte_identical_to_anchor": True,
            "new_internal_campaign_created": False,
        },
        "tree_projection": {
            "source_leaf_count": SOURCE_LEAF_COUNT,
            "source_tree_sha256": SOURCE_TREE_SHA256,
            "projected_leaf_count": PROJECTED_LEAF_COUNT,
            "projected_tree_sha256": PROJECTED_TREE_SHA256,
            "excluded_paths": list(EXCLUDED_PATHS),
        },
        "claims": {
            "source_store_mutated": False,
            "accepted_physical_state_advanced": False,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "SGBL_or_FGCQR_opened": False,
            "mechanism_or_physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
    }
    return {**body, "receipt_sha256": sha256(_canonical(body)).hexdigest()}


def _bind_receipt(wrapper: _Snapshot) -> None:
    leaf = wrapper.leaves.get("bootstrap-receipt.json")
    if leaf is None or leaf.byte_count != 3070 or leaf.raw_sha256 != RECEIPT_RAW_SHA256:
        _stop("PREF1_RECEIPT_DRIFT", "bootstrap receipt raw identity differs")
    value = _json(leaf.raw, "bootstrap receipt")
    _object_address(value, "receipt_sha256", RECEIPT_SHA256, "bootstrap receipt")
    if value != _expected_receipt():
        _stop("PREF1_RECEIPT_DRIFT", "bootstrap receipt contract differs")


def expected_anchor(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "authority": dict(config["authority"]),
        "source": dict(config["source"]),
        "source_terminal": dict(config["source_terminal"]),
        "projection": dict(config["projection"]),
        "receipt": dict(config["receipt"]),
        "anchor": dict(config["anchor"]),
        "members": [dict(item) for item in config["members"]],
        "validation": dict(config["validation"]),
    }


def bind_installed_projection(config_raw: bytes, repository: Path) -> dict[str, Any]:
    """Authenticate source, projection, receipt, lineage, and all payloads."""
    _config(config_raw)
    repository = Path(repository)
    _bind_authority(repository)
    source = _snapshot(repository, SOURCE_STORE)
    wrapper = _snapshot(repository, DESTINATION_WRAPPER)
    if source.directories != STORE_DIRECTORIES or wrapper.directories != WRAPPER_DIRECTORIES:
        _stop("PREF1_TREE_DRIFT", "source or wrapper directory grammar differs")
    source_leaves = source.leaves
    destination_leaves = {
        path.removeprefix("calibration/"): _Leaf(
            path.removeprefix("calibration/"), leaf.byte_count,
            leaf.raw_sha256, leaf.raw,
        )
        for path, leaf in wrapper.leaves.items() if path.startswith("calibration/")
    }
    destination = _Snapshot(destination_leaves, STORE_DIRECTORIES)
    _validate_tree_grammar(source, source=True)
    _validate_tree_grammar(destination, source=False)
    _manifest(source, EXPECTED_SOURCE, "source terminal")
    _manifest(destination, EXPECTED_PROJECTION, "installed projection")
    expected_paths = set(source_leaves) - set(EXCLUDED_PATHS)
    if set(destination.leaves) != expected_paths:
        _stop("PREF1_PROJECTION_DRIFT", "projected path set differs")
    for path in sorted(expected_paths):
        source_leaf = source_leaves[path]
        destination_leaf = destination.leaves[path]
        if (
            source_leaf.byte_count != destination_leaf.byte_count
            or source_leaf.raw_sha256 != destination_leaf.raw_sha256
            or source_leaf.raw != destination_leaf.raw
        ):
            _stop("PREF1_PROJECTION_DRIFT", f"projected bytes differ: {path}")
    if set(EXCLUDED_PATHS) & set(destination.leaves):
        _stop("PREF1_PROJECTION_DRIFT", "terminal/excluded suffix entered projection")
    if set(wrapper.leaves) != {"bootstrap-receipt.json"} | {
        f"calibration/{path}" for path in destination.leaves
    }:
        _stop("PREF1_TREE_DRIFT", "wrapper has a foreign leaf")
    _validate_source_terminal(source)
    _bind_receipt(wrapper)
    checkpoints, _journals = _validate_lineage(destination.leaves)
    descriptors = _descriptor_objects(destination.leaves, checkpoints)
    payloads = _load_payloads(destination.leaves, descriptors)
    members = _member_evidence(checkpoints[ANCHOR_GENERATION], descriptors, payloads)
    if members != list(EXPECTED_MEMBERS):
        _stop("PREF1_MEMBER_DRIFT", "six-member evidence ordering differs")
    # Close the read window by re-reading both trees.  This does not invoke the
    # installer or HLT16 decision code.
    if source != _snapshot(repository, SOURCE_STORE) or wrapper != _snapshot(repository, DESTINATION_WRAPPER):
        _stop("PREF1_PATH_UNSAFE", "source or projection changed during binding")
    return expected_anchor(config_raw)


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
            "live_anchor": expected_anchor(config_raw),
            "scope": dict(config["scope"]),
            "claims": dict(config["claims"]),
            "nonclaims": list(config["nonclaims"]),
        },
    }


def build_pref1_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    config = _config(config_raw)
    if bind_installed_projection(config_raw, repository) != expected_anchor(config_raw):
        _stop("PREF1_PROJECTION_DRIFT", "live/config projection binding differs")
    return _result_from_config(config_raw, config)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    value = _json(result_raw, "PREF1 compact result", pretty=True)
    if value != _result_from_config(config_raw, config):
        _stop("PREF1_COMPACT_DRIFT", "PREF1 compact result differs")
    return value


__all__ = [
    "ARTIFACT_ID", "CONFIG_PATH", "RESULT_PATH", "SOURCE_STORE",
    "DESTINATION_STORE", "TDG8RCV3PREF1Error", "bind_installed_projection",
    "build_pref1_result", "canonical_result", "expected_anchor",
    "validate_compact_result",
]
