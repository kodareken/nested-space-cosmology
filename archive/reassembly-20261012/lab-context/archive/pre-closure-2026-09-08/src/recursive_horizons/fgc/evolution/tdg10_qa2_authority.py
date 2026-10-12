"""Prospective authority for two-width exact complete-C robustness.

TDG10-QA2 may authorize one diagnostic, non-accepted measurement of the
immutable retry-4 and retry-5 GR-0 proposals under SSPRK3-on-inherited-SBP4
with the already-qualified exact radius-free complete-C interval owner.  It
does not authorize a campaign-state advance, a new protocol ID, diagnostic
endpoint serialization, or restoration of retry 3.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
import platform
from pathlib import Path
import re
import stat
import sys
from typing import Mapping, NoReturn

from . import tdg10_qa1_authority as qa1
from . import tdg9_ti2_authority as ti2
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG10-QA2-FRZ1"
CLASSIFICATION = (
    "premise_only_retries4_5_ssprk3_sbp4_exact_complete_C_two_width_robustness"
)
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa2-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa2-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg10-qa2-frz1.md"
AUTHORITY_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg10_qa2_authority.py"
)
RUNNER_PATH = "scripts/run_fgc_tdg10_qa2.py"
EXACT_RUNTIME_PATH = qa1.EXACT_RUNTIME_PATH
EXACT_ADMISSION_PATH = qa1.EXACT_ADMISSION_PATH
TI2_RUNNER_PATH = qa1.TI2_RUNNER_PATH
QA1_PREF1_ARTIFACT_ID = "FGC-1-TDG10-QA1-PREF1"
QA1_PREF1_RESULT_PATH = "results/fgc-1-tdg10-qa1-pref1.json"
QA1_PREF1_RESULT_SHA256 = (
    "3a107bcede479df4326aac5e194beefa1042860a8be498d7dc5d0519fdbc8603"
)
OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/tdg10-qa2/retries4-5-ssprk3-sbp4-exact-complete-c"
)
STAGING_PREFIX = ".retries4-5-ssprk3-sbp4-exact-complete-c-qa2.stage-"
BASE_COMMIT = "a7915925a23ceaff6a491b8d08c572a70d242985"

STORE_PATH = qa1.STORE_PATH
SEALED_STORE_LEAF_COUNT = qa1.SEALED_STORE_LEAF_COUNT
SEALED_STORE_SNAPSHOT_SHA256 = qa1.SEALED_STORE_SNAPSHOT_SHA256
MEMBER_KEY = qa1.MEMBER_KEY
POINT_COUNT = qa1.POINT_COUNT
OWNED_ROW_COUNT = qa1.OWNED_ROW_COUNT
ACCEPTED_TIME_HEX = qa1.ACCEPTED_TIME_HEX
PHYSICAL_STATE_SHA256 = qa1.PHYSICAL_STATE_SHA256
MEMBER_DESCRIPTOR_SHA256 = qa1.MEMBER_DESCRIPTOR_SHA256
COORDINATES_SHA256 = qa1.COORDINATES_SHA256
CURSOR_MODE = qa1.CURSOR_MODE
PENDING_OWNER = qa1.PENDING_OWNER
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
CHANNEL_COUNT = 18
TABLEAU_SELECTOR = qa1.TABLEAU_SELECTOR
TABLEAU_RUNTIME_SELECTOR = qa1.TABLEAU_RUNTIME_SELECTOR
ACTUAL_SPATIAL_OPERATOR = qa1.ACTUAL_SPATIAL_OPERATOR
INTERVAL_OWNER = qa1.INTERVAL_OWNER
PRIMARY_REFINEMENT_DEPTH = qa1.PRIMARY_REFINEMENT_DEPTH
PER_CHANNEL_MAXIMUM_CANDIDATES_D01 = qa1.PER_CHANNEL_MAXIMUM_CANDIDATES_D01
PER_CHANNEL_MAXIMUM_CANDIDATES_D12 = qa1.PER_CHANNEL_MAXIMUM_CANDIDATES_D12
SELECTED_RETRIES = (4, 5)
FORBIDDEN_RETRY3_GENERATION = 9
LEDGER_OWNER = ti2.SEMANTICS["ledger_owner"]
TRACER_OWNER = ti2.SEMANTICS["tracer_owner"]
TRANSACTION_OWNER = ti2.SEMANTICS["transaction_owner"]
RHS_OWNER = ti2.SEMANTICS["RHS_owner"]
PROJECTOR_OWNER = ti2.SEMANTICS["projector_owner"]

WIDTHS = (
    {
        "retry": 4,
        "prior_retry_count": 3,
        "generation": 10,
        "forbidden_predecessor_generation": 11,
        "checkpoint_sha256": (
            "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187"
        ),
        "checkpoint_raw_sha256": (
            "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b"
        ),
        "checkpoint_journal_sequence": 12,
        "checkpoint_journal_sha256": (
            "d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4"
        ),
        "checkpoint_journal_raw_sha256": (
            "dd02233dd5e4cf28becbc4fe29259d0879005d806b70493aa2281255081eb219"
        ),
        "historical_rejection_sequence": 13,
        "historical_rejection_sha256": (
            "1bd21553f3b6e0af613997da094a46d4dcd966a0e4346692bae62b05046862ca"
        ),
        "historical_rejection_raw_sha256": (
            "26542403a17b05664d5a447b88860186d5f554b1f4d555624493ed40b941694c"
        ),
        "attempted_width_hex": "0x1.aaa9612df8000p-12",
        "transaction_sha256": (
            "abab435c8222de0a86368e9562a2949578de38fe822c6032319410b49e4c6f6c"
        ),
        "historical_journal_sha256": (
            "fc68ce01980cb7a66f82de59bbd59d5e53ace08537b6b0dbfe0be1e070e6912f"
        ),
    },
    {
        "retry": 5,
        "prior_retry_count": 4,
        "generation": 11,
        "forbidden_predecessor_generation": 12,
        "checkpoint_sha256": (
            "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8"
        ),
        "checkpoint_raw_sha256": (
            "786bdcb835e043c65ddd33fad3a2af0e3d33e768dc26a9f9f00d2757b0da3ea5"
        ),
        "checkpoint_journal_sequence": 14,
        "checkpoint_journal_sha256": (
            "0810018134df958b8dab32bd9c0ecd7ebfc4376fbc87869a8bee7dd02eb1e199"
        ),
        "checkpoint_journal_raw_sha256": (
            "a1882d904091591d7fa3ccf3eab04aa748161115f97c15abd57e535e6fa01f3d"
        ),
        "historical_rejection_sequence": 15,
        "historical_rejection_sha256": (
            "475014d918952adef032d71b81073773364a37fa795491e739aae609f77e8537"
        ),
        "historical_rejection_raw_sha256": (
            "1997342415f74f330ce6b45ca6fd550ad0e5e674b85148fea5047ff82ac3f152"
        ),
        "attempted_width_hex": "0x1.aaa9612df0000p-13",
        "transaction_sha256": (
            "69fd069034fdd3e2a80a5f3fa989b639b6659105a8ff3e8742e32e78739fd461"
        ),
        "historical_journal_sha256": (
            "b602d12ab00757090734ac63042f172a372154d51b41c7d286c7c0a758b62221"
        ),
    },
)
REPLAYS = tuple(
    dict(item) for item in ti2.REPLAYS if int(item["retry"]) in SELECTED_RETRIES
)
ENVIRONMENT = dict(qa1.ENVIRONMENT)
IMPLEMENTATION_INVENTORY = (
    ("runner_path", RUNNER_PATH),
    ("authority_path", AUTHORITY_PATH),
    ("exact_runtime_path", EXACT_RUNTIME_PATH),
    ("exact_admission_path", EXACT_ADMISSION_PATH),
    ("TI2_runner_path", TI2_RUNNER_PATH),
)
IMPLEMENTATION_DIGEST_KEYS = tuple(
    key.replace("_path", "_sha256") for key, _path in IMPLEMENTATION_INVENTORY
)
RUNTIME_BYTE_PATHS = (
    CONFIG_PATH,
    RESULT_PATH,
    RUNNER_PATH,
    AUTHORITY_PATH,
)
DELTA_PATHS = tuple(
    sorted(
        (
            "Makefile",
            "README.md",
            CONFIG_PATH,
            "docs/claim-ledger.md",
            "docs/fgc-runtime-matrix.md",
            OWNER_DOCUMENT,
            "docs/research-roadmap.md",
            RESULT_PATH,
            "results/README.md",
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg10_qa2_frz1.py",
            RUNNER_PATH,
            AUTHORITY_PATH,
            "tests/test_check_repo_tdg10_qa2_frz1.py",
            "tests/test_fgc_tdg10_qa2_authority.py",
            "tests/test_fgc_tdg10_qa2_runner.py",
        )
    )
)
ALLOWED_UNTRACKED = qa1.ALLOWED_UNTRACKED
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


class QA2AuthorityError(RuntimeError):
    """The prospective QA2 image differs from its frozen robustness authority."""


def _fail(message: str) -> NoReturn:
    raise QA2AuthorityError(message)


def canonical_pretty(value: object) -> bytes:
    return qa1.canonical_pretty(value)


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def _git(root: Path, *arguments: str) -> bytes:
    try:
        return qa1._git(root, *arguments)
    except qa1.QA1AuthorityError as exc:
        raise QA2AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def read_leaf(root: Path, relative: str) -> bytes:
    """Read one regular leaf without following any parent or leaf symlink."""

    try:
        return qa1.read_leaf(root, relative)
    except qa1.QA1AuthorityError as exc:
        raise QA2AuthorityError(str(exc)) from exc


def _working(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    return (
        _paths(root, "diff", "--name-only", "-z", "--"),
        _paths(root, "diff", "--cached", "--name-only", "-z", "--"),
        _paths(root, "ls-files", "--others", "--exclude-standard", "-z", "--"),
    )


def _require_no_replace_refs(root: Path) -> None:
    if _git(root, "for-each-ref", "--format=%(refname)", "refs/replace").strip():
        _fail("Git replace refs are present")


def _environment() -> dict[str, str]:
    import numpy as np

    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
    }
    if observed != ENVIRONMENT:
        _fail("execution environment differs")
    return observed


def width_spec(retry: int) -> dict[str, object]:
    for item in WIDTHS:
        if int(item["retry"]) == int(retry):
            return dict(item)
    _fail(f"QA2 width is not selected: {retry}")


def _require_frozen_replays() -> None:
    selected = tuple(
        dict(item) for item in ti2.REPLAYS if int(item["retry"]) in SELECTED_RETRIES
    )
    if tuple(int(item["retry"]) for item in selected) != SELECTED_RETRIES:
        _fail("frozen TI2 retry-4/5 replays differ")
    if selected != REPLAYS:
        _fail("QA2 replay binding differs from TI2")
    for replay, width in zip(selected, WIDTHS, strict=True):
        if (
            int(replay["retry"]) != width["retry"]
            or int(replay["prior_retry_count"]) != width["prior_retry_count"]
            or int(replay["predecessor_generation"]) != width["generation"]
            or int(replay["predecessor_generation"])
            == width["forbidden_predecessor_generation"]
            or int(replay["predecessor_generation"]) == FORBIDDEN_RETRY3_GENERATION
            or str(replay["checkpoint_sha256"]) != width["checkpoint_sha256"]
            or str(replay["checkpoint_raw_sha256"]) != width["checkpoint_raw_sha256"]
            or int(replay["journal_sequence"]) != width["historical_rejection_sequence"]
            or str(replay["journal_sha256"]) != width["historical_rejection_sha256"]
            or str(replay["journal_raw_sha256"])
            != width["historical_rejection_raw_sha256"]
            or str(replay["attempted_width_hex"]) != width["attempted_width_hex"]
            or str(replay["transaction_sha256"]) != width["transaction_sha256"]
        ):
            _fail(f"frozen retry-{width['retry']} replay differs")


def _width_block(width: Mapping[str, object]) -> dict[str, object]:
    return {
        "retry": width["retry"],
        "prior_retry_count": width["prior_retry_count"],
        "generation": width["generation"],
        "forbidden_predecessor_generation": width["forbidden_predecessor_generation"],
        "checkpoint_sha256": width["checkpoint_sha256"],
        "checkpoint_raw_sha256": width["checkpoint_raw_sha256"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "checkpoint_journal_sha256": width["checkpoint_journal_sha256"],
        "checkpoint_journal_raw_sha256": width["checkpoint_journal_raw_sha256"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "historical_rejection_sha256": width["historical_rejection_sha256"],
        "historical_rejection_raw_sha256": width["historical_rejection_raw_sha256"],
        "attempted_width_binary64_hex": width["attempted_width_hex"],
        "transaction_sha256": width["transaction_sha256"],
        "historical_journal_sha256": width["historical_journal_sha256"],
    }


def _predecessor() -> dict[str, object]:
    retry_4, retry_5 = WIDTHS
    return {
        "store": STORE_PATH,
        "store_leaf_count": SEALED_STORE_LEAF_COUNT,
        "store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        "member_key": MEMBER_KEY,
        "point_count": POINT_COUNT,
        "owned_row_count": OWNED_ROW_COUNT,
        "selected_retries": [4, 5],
        "accepted_time_binary64_hex": ACCEPTED_TIME_HEX,
        "physical_state_sha256": PHYSICAL_STATE_SHA256,
        "member_descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "coordinates_sha256": COORDINATES_SHA256,
        "cursor_mode": CURSOR_MODE,
        "pending_owner": PENDING_OWNER,
        "ledger_owner": LEDGER_OWNER,
        "tracer_owner": TRACER_OWNER,
        "transaction_owner": TRANSACTION_OWNER,
        "RHS_owner": RHS_OWNER,
        "projector_owner": PROJECTOR_OWNER,
        "immediate_predecessor_artifact_id": QA1_PREF1_ARTIFACT_ID,
        "immediate_predecessor_commit": BASE_COMMIT,
        "immediate_predecessor_result_path": QA1_PREF1_RESULT_PATH,
        "immediate_predecessor_result_sha256": QA1_PREF1_RESULT_SHA256,
        "immediate_predecessor_licenses_qa2_method_design": True,
        "immediate_predecessor_licenses_old_member_adoption": False,
        "retry_4": _width_block(retry_4),
        "retry_5": _width_block(retry_5),
    }


def _formulation() -> dict[str, object]:
    return {
        "source_member_method": "RK4",
        "source_member_integrator": (
            "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4"
        ),
        "tableau_selector": TABLEAU_SELECTOR,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "actual_spatial_operator": ACTUAL_SPATIAL_OPERATOR,
        "interval_owner": INTERVAL_OWNER,
        "channel_order_owner": "TDG6_COMPLETE_STATE_CHANNELS",
        "channel_count": CHANNEL_COUNT,
        "complete_admission_is_all_of": True,
        "minimum_observed_order": "3/2",
        "order_squared_multiplier": 8,
        "equality_passes": True,
        "outer_subintervals": 2,
        "finest_subintervals": 4,
        "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
        "per_channel_maximum_candidates_D01": PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "per_channel_maximum_candidates_D12": PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "independent_route_required": True,
        "absolute_tolerance_used": False,
        "physical_signal_used_for_normalization": False,
        "declared_cubic_is_exact_PDE_history": False,
        "pde_trajectory_order_certified": False,
        "both_widths_pass_licenses_production_method_design": True,
        "any_width_nonpass_rejects_this_exact_remedy_on_tested_neighborhood": True,
    }


def _transaction() -> dict[str, object]:
    return {
        "width_count": 2,
        "shadow_path_count_per_width": 7,
        "shadow_proposal_count_per_width": 7,
        "SSPRK3_records_per_proposal": 4,
        "maximum_stage_and_endpoint_records_per_width": 28,
        "maximum_stage_and_endpoint_records": 56,
        "historical_store_read_only": True,
        "output_namespace_must_be_absent": True,
        "overwrite_forbidden": True,
        "diagnostic_fine_endpoint_serialized": False,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "campaign_state_write_authorized": False,
        "PDE_state_commit_authorized": False,
        "fine_path_commit_authorized": False,
        "independent_postrun_binder_required": True,
        "later_production_method_requires_separate_freeze": True,
    }


def _exclusions() -> dict[str, object]:
    return {
        "retry_3_authorized": False,
        "fourth_width_authorized": False,
        "new_grid_authorized": False,
        "threshold_change_authorized": False,
        "historical_module_mutation_authorized": False,
        "diagnostic_endpoint_serialization_authorized": False,
        "endpoint_transplant_authorized": False,
        "old_member_adoption_authorized": False,
        "common_event_authorized": False,
        "GR0_calibration_authorized": False,
        "production_SSPRK3_comparator_earned": False,
        "independent_method_agreement_earned": False,
        "production_method_earned": False,
        "width_robustness_passed": False,
        "RA1_authorized": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "candidate_execution_authorized": False,
        "mechanism_result_authorized": False,
        "physical_result_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "external_publication_authorized": False,
        "push_authorized": False,
    }


def _require_implementation_block(value: object) -> dict[str, object]:
    """Validate the exact path schema and config-supplied hashes. No self-hash constants."""

    if not isinstance(value, Mapping):
        _fail("QA2 implementation block is required")
    expected_keys = []
    for key, path in IMPLEMENTATION_INVENTORY:
        expected_keys.append(key)
        expected_keys.append(key.replace("_path", "_sha256"))
    if set(value) != set(expected_keys):
        _fail("QA2 implementation inventory keys differ")
    answer: dict[str, object] = {}
    for key, path in IMPLEMENTATION_INVENTORY:
        digest_key = key.replace("_path", "_sha256")
        declared = value.get(key)
        digest = value.get(digest_key)
        if declared != path:
            _fail(f"QA2 implementation path differs: {key}")
        if not isinstance(digest, str) or not _SHA256.fullmatch(digest):
            _fail(f"QA2 implementation hash is malformed: {digest_key}")
        answer[key] = path
        answer[digest_key] = digest
    return answer


def _expected_config(implementation: Mapping[str, object]) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "base_commit": BASE_COMMIT,
        "output_namespace": OUTPUT_NAMESPACE,
        "predecessor": _predecessor(),
        "formulation": _formulation(),
        "transaction": _transaction(),
        "environment": dict(ENVIRONMENT),
        "exclusions": _exclusions(),
        "implementation": dict(implementation),
    }


def parse_config(raw: bytes) -> dict[str, object]:
    import tomllib

    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise QA2AuthorityError("invalid QA2 TOML") from exc
    _require_frozen_replays()
    implementation = _require_implementation_block(value.get("implementation"))
    expected = _expected_config(implementation)
    if value != expected:
        _fail("QA2 config differs from frozen contract")
    predecessor = value["predecessor"]
    retry_4 = predecessor["retry_4"]
    retry_5 = predecessor["retry_5"]
    if (
        predecessor["selected_retries"] != [4, 5]
        or predecessor["immediate_predecessor_artifact_id"] != QA1_PREF1_ARTIFACT_ID
        or predecessor["immediate_predecessor_commit"] != BASE_COMMIT
        or predecessor["immediate_predecessor_result_path"] != QA1_PREF1_RESULT_PATH
        or predecessor["immediate_predecessor_result_sha256"] != QA1_PREF1_RESULT_SHA256
        or predecessor["immediate_predecessor_licenses_qa2_method_design"] is not True
        or predecessor["immediate_predecessor_licenses_old_member_adoption"] is not False
        or retry_4["generation"] != 10
        or retry_4["generation"] == retry_4["forbidden_predecessor_generation"]
        or retry_4["checkpoint_journal_sequence"] != 12
        or retry_4["historical_rejection_sequence"] != 13
        or retry_4["retry"] != 4
        or retry_5["generation"] != 11
        or retry_5["generation"] == retry_5["forbidden_predecessor_generation"]
        or retry_5["checkpoint_journal_sequence"] != 14
        or retry_5["historical_rejection_sequence"] != 15
        or retry_5["retry"] != 5
        or predecessor["cursor_mode"] != CURSOR_MODE
        or predecessor["pending_owner"] != PENDING_OWNER
        or predecessor["ledger_owner"] != LEDGER_OWNER
        or predecessor["tracer_owner"] != TRACER_OWNER
        or predecessor["transaction_owner"] != TRANSACTION_OWNER
    ):
        _fail("QA2 retry-4/5 predecessor boundary differs")
    if value["target_protocol"] != TARGET_PROTOCOL:
        _fail("QA2 target protocol differs")
    if value["project_version"] != PROJECT_VERSION:
        _fail("QA2 project version differs")
    if any(value["exclusions"].values()) or any(
        value["transaction"][name] is not False
        for name in (
            "diagnostic_fine_endpoint_serialized",
            "diagnostic_fine_endpoint_is_accepted_state",
            "campaign_state_write_authorized",
            "PDE_state_commit_authorized",
            "fine_path_commit_authorized",
        )
    ):
        _fail("QA2 authorized a state, campaign, candidate, or physics flag")
    return value


def require_output_absent(root: Path) -> None:
    target = root / OUTPUT_NAMESPACE
    current = root
    for part in target.relative_to(root).parent.parts:
        current /= part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("QA2 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("QA2 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in target.parent.iterdir()):
        _fail("QA2 staging namespace already exists")


def filled_implementation_pairs(
    config: Mapping[str, object],
) -> tuple[tuple[str, str], ...]:
    """Return the required (path, sha256) pairs from the parsed implementation block."""

    block = _require_implementation_block(config.get("implementation"))
    return tuple(
        (path, str(block[key.replace("_path", "_sha256")]))
        for key, path in IMPLEMENTATION_INVENTORY
    )


def verify_implementation_inventory(
    root: Path,
    config: Mapping[str, object],
    *,
    commit: str | None,
) -> tuple[tuple[str, str], ...]:
    pairs = filled_implementation_pairs(config)
    if len(pairs) != len(IMPLEMENTATION_INVENTORY):
        _fail("QA2 implementation inventory is incomplete")
    for path, digest in pairs:
        raw = (
            read_leaf(root, path)
            if commit is None
            else _git(root, "show", f"{commit}:{path}")
        )
        if sha256(raw).hexdigest() != digest:
            _fail(f"QA2 implementation binding differs: {path}")
    return pairs


def _verify_qa1_pref1_result(root: Path) -> None:
    live = read_leaf(root, QA1_PREF1_RESULT_PATH)
    if sha256(live).hexdigest() != QA1_PREF1_RESULT_SHA256:
        _fail("live QA1-PREF1 compact result differs")
    committed = _git(root, "show", f"{BASE_COMMIT}:{QA1_PREF1_RESULT_PATH}")
    if sha256(committed).hexdigest() != QA1_PREF1_RESULT_SHA256:
        _fail("committed QA1-PREF1 compact result differs")
    if committed != live:
        _fail("QA1-PREF1 compact result drifted from the base commit")


def _verify_predecessor_leaves(root: Path) -> None:
    for width in WIDTHS:
        generation = int(width["generation"])
        checkpoint = (
            f"{STORE_PATH}/checkpoints/"
            f"{generation:020d}-{width['checkpoint_sha256']}.json"
        )
        journal_tip = (
            f"{STORE_PATH}/journal/"
            f"{int(width['checkpoint_journal_sequence']):020d}-"
            f"{width['checkpoint_journal_sha256']}.journal"
        )
        rejection = (
            f"{STORE_PATH}/journal/"
            f"{int(width['historical_rejection_sequence']):020d}-"
            f"{width['historical_rejection_sha256']}.journal"
        )
        for relative, digest in (
            (checkpoint, width["checkpoint_raw_sha256"]),
            (journal_tip, width["checkpoint_journal_raw_sha256"]),
            (rejection, width["historical_rejection_raw_sha256"]),
        ):
            if sha256(read_leaf(root, relative)).hexdigest() != digest:
                _fail(f"predecessor identity differs: {relative}")
        checkpoint_raw = json.loads(read_leaf(root, checkpoint))
        if (
            checkpoint_raw.get("generation") != generation
            or checkpoint_raw.get("checkpoint_sha256") != width["checkpoint_sha256"]
            or checkpoint_raw.get("journal_sequence")
            != width["checkpoint_journal_sequence"]
            or checkpoint_raw.get("journal_tip_sha256")
            != width["checkpoint_journal_sha256"]
        ):
            _fail(
                f"generation-{generation} checkpoint tip is not journal sequence "
                f"{width['checkpoint_journal_sequence']}"
            )
        journal_raw = json.loads(read_leaf(root, journal_tip))
        if (
            journal_raw.get("sequence") != width["checkpoint_journal_sequence"]
            or journal_raw.get("record_sha256") != width["checkpoint_journal_sha256"]
            or journal_raw.get("kind") == "tdg6_rejection"
        ):
            _fail(
                f"checkpoint journal sequence {width['checkpoint_journal_sequence']} "
                "is not the accepted generation tip"
            )
        rejection_raw = json.loads(read_leaf(root, rejection))
        if (
            rejection_raw.get("sequence") != width["historical_rejection_sequence"]
            or rejection_raw.get("record_sha256")
            != width["historical_rejection_sha256"]
            or rejection_raw.get("kind") != "tdg6_rejection"
        ):
            _fail(
                f"sequence {width['historical_rejection_sequence']} is not the "
                f"historical retry-{width['retry']} rejection"
            )


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return qa1._snapshot_store(root)
    except qa1.QA1AuthorityError as exc:
        raise QA2AuthorityError("sealed store cannot be authenticated") from exc


def _expected_fingerprint_receipt(width: Mapping[str, object]) -> dict[str, object]:
    return {
        "retry": width["retry"],
        "predecessor_generation": width["generation"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "transaction_sha256": width["transaction_sha256"],
        "state_sha256": PHYSICAL_STATE_SHA256,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "shadow_proposal_constructed": False,
    }


def expected_raw_replay_receipt(retry: int) -> dict[str, object]:
    width = width_spec(retry)
    return {
        "member_key": MEMBER_KEY,
        "retry": width["retry"],
        "predecessor_generation": width["generation"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "attempted_width_hex": width["attempted_width_hex"],
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": width["transaction_sha256"],
        "historical_journal_sha256": width["historical_journal_sha256"],
        "ledger_owner": LEDGER_OWNER,
        "tracer_owner": TRACER_OWNER,
        "transaction_owner": TRANSACTION_OWNER,
    }


def real_predecessor_fingerprints(root: Path) -> tuple[dict[str, object], ...]:
    """Restore only retry-4/5 RK4-2049 predecessors and fingerprint without proposals."""

    from scripts import run_fgc_tdg9_ti2 as runner

    from . import tdg9_ar1_authority as ar1
    from .hlt16_campaign_store import HLT16CampaignStore
    from .proto19_gr0_static_factory import build_static_gr0_shells

    _require_frozen_replays()
    repository = root.resolve()
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    receipts: list[dict[str, object]] = []
    for replay, width in zip(REPLAYS, WIDTHS, strict=True):
        if int(replay["predecessor_generation"]) == FORBIDDEN_RETRY3_GENERATION:
            _fail("QA2 must not restore generation 9")
        if int(replay["predecessor_generation"]) == width[
            "forbidden_predecessor_generation"
        ]:
            _fail(
                f"retry-{width['retry']} must not restore generation "
                f"{width['forbidden_predecessor_generation']}"
            )
        try:
            restored = runner._restore_replay(repository, store, shells, replay)
            fingerprint = restored.fingerprint
        except Exception as exc:
            raise QA2AuthorityError(
                f"retry-{width['retry']} predecessor fingerprint failed"
            ) from exc
        if (
            fingerprint["transaction_sha256"] != width["transaction_sha256"]
            or fingerprint["state_sha256"] != PHYSICAL_STATE_SHA256
            or fingerprint["accepted_time_hex"] != ACCEPTED_TIME_HEX
            or fingerprint["descriptor_sha256"] != MEMBER_DESCRIPTOR_SHA256
        ):
            _fail(f"retry-{width['retry']} predecessor fingerprint differs")
        receipts.append(
            {
                "retry": width["retry"],
                "predecessor_generation": width["generation"],
                "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
                "historical_rejection_sequence": width["historical_rejection_sequence"],
                "transaction_sha256": fingerprint["transaction_sha256"],
                "state_sha256": fingerprint["state_sha256"],
                "accepted_time_hex": fingerprint["accepted_time_hex"],
                "descriptor_sha256": fingerprint["descriptor_sha256"],
                "shadow_proposal_constructed": False,
            }
        )
    return tuple(receipts)


def expected_compact(config_raw: bytes) -> dict[str, object]:
    parsed = parse_config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            **_expected_config(parsed["implementation"]),
            "store_identity_observed_at_prelaunch": {
                "leaf_count": SEALED_STORE_LEAF_COUNT,
                "sha256": SEALED_STORE_SNAPSHOT_SHA256,
            },
            "QA2_output_namespace_absent_at_prelaunch": True,
            "replace_refs_absent": True,
            "restored_predecessor_fingerprints": [
                _expected_fingerprint_receipt(width) for width in WIDTHS
            ],
            "shadow_executed": False,
            "shadow_proposals_constructed": 0,
            "diagnostic_qualification_only": True,
            "diagnostic_fine_endpoint_serialized": False,
            "new_protocol_id_authorized": False,
            "target_protocol": TARGET_PROTOCOL,
            "retry4_restores_generation": 10,
            "retry4_never_restores_generation_11": True,
            "checkpoint_journal_is_sequence_12": True,
            "historical_retry4_rejection_is_sequence_13": True,
            "retry5_restores_generation": 11,
            "retry5_never_restores_generation_12": True,
            "checkpoint_journal_is_sequence_14": True,
            "historical_retry5_rejection_is_sequence_15": True,
            "width_robustness_passed": False,
            "production_method_earned": False,
            "QA1_licenses_QA2_method_design_never_old_member_adoption": True,
            "state_advance_authorized": False,
            "campaign_state_write_authorized": False,
            "PDE_state_commit_authorized": False,
            "fine_path_commit_authorized": False,
            "candidate_execution_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
        },
    }


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    *,
    store_snapshot: tuple[int, str],
) -> dict[str, object]:
    repository = root.resolve()
    config = parse_config(config_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not the sealed QA2 base commit")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted(
            (
                *unstaged,
                *(item for item in untracked if item not in ALLOWED_UNTRACKED),
            )
        )
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective QA2 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    verify_implementation_inventory(repository, config, commit=None)
    _verify_qa1_pref1_result(repository)
    _verify_predecessor_leaves(repository)
    sealed = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256)
    if store_snapshot != sealed:
        _fail("sealed store differs at QA2 prelaunch")
    require_output_absent(repository)
    fingerprints = list(real_predecessor_fingerprints(repository))
    if fingerprints != [_expected_fingerprint_receipt(width) for width in WIDTHS]:
        _fail("restored predecessor fingerprints differ")
    value = expected_compact(config_raw)
    return {
        **value,
        "artifact_payload": {
            **value["artifact_payload"],
            "store_identity_observed_at_prelaunch": {
                "leaf_count": store_snapshot[0],
                "sha256": store_snapshot[1],
            },
            "restored_predecessor_fingerprints": fingerprints,
        },
    }


def _reject_compact_json_constant(token: str) -> NoReturn:
    _fail(f"nonfinite QA2 compact JSON constant: {token}")


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    expected = expected_compact(config_raw)
    try:
        result = json.loads(
            result_raw,
            object_pairs_hook=_unique,
            parse_constant=_reject_compact_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QA2AuthorityError("invalid QA2 compact JSON") from exc
    if result != expected or result_raw != canonical_pretty(result):
        _fail("QA2 compact result differs")
    payload = result["artifact_payload"]
    if payload["shadow_executed"] is not False or payload[
        "shadow_proposals_constructed"
    ] != 0:
        _fail("QA2 compact result claimed shadow work")
    if payload["state_advance_authorized"] is not False:
        _fail("QA2 compact result authorized a state advance")
    if payload["width_robustness_passed"] is not False:
        _fail("QA2 compact result claimed width robustness")
    if payload["production_method_earned"] is not False:
        _fail("QA2 compact result earned a production method")
    if payload["diagnostic_fine_endpoint_serialized"] is not False:
        _fail("QA2 compact result serialized a diagnostic endpoint")
    for forbidden in (
        "complete_admission_passed",
        "failed_channels",
        "diagnostic_payload_sha256",
        "fine_endpoint_descriptor_sha256",
    ):
        if forbidden in payload:
            _fail("QA2 compact result claimed a diagnostic outcome")
    return result


def _require_bytes_match_commit(root: Path, commit: str, relative: str) -> None:
    if _git(root, "show", f"{commit}:{relative}") != read_leaf(root, relative):
        _fail(f"working {relative} differs from the authority commit")


@dataclass(frozen=True, slots=True)
class QA2Authority:
    authority_commit: str
    execution_authorized: bool = True
    diagnostic_qualification_only: bool = True
    tableau_selector: str = TABLEAU_SELECTOR
    actual_spatial_operator: str = ACTUAL_SPATIAL_OPERATOR
    selected_retries: tuple[int, int] = SELECTED_RETRIES
    target_protocol: str = TARGET_PROTOCOL
    new_protocol_id_authorized: bool = False
    state_advance_authorized: bool = False
    campaign_state_write_authorized: bool = False
    PDE_state_commit_authorized: bool = False
    fine_path_commit_authorized: bool = False
    diagnostic_fine_endpoint_serialized: bool = False
    candidate_execution_authorized: bool = False
    mechanism_result_authorized: bool = False
    physical_result_authorized: bool = False
    production_SSPRK3_comparator_earned: bool = False
    independent_method_agreement_earned: bool = False
    production_method_earned: bool = False
    width_robustness_passed: bool = False


def authorize(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> QA2Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("QA2 authority commit is malformed")
    config = validate_compact(config_raw, result_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = (
        _git(repository, "rev-list", "--parents", "-n", "1", authority_commit)
        .decode()
        .split()[1:]
    )
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("QA2 authority is not the exact direct successor")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("QA2 authority working image is not clean")
    delta = tuple(
        sorted(
            _paths(
                repository,
                "diff",
                "--name-only",
                "-z",
                BASE_COMMIT,
                authority_commit,
                "--",
            )
        )
    )
    if delta != DELTA_PATHS:
        _fail("committed QA2 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed QA2 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed QA2 compact result differs")
    for relative in RUNTIME_BYTE_PATHS:
        _require_bytes_match_commit(repository, authority_commit, relative)
    verify_implementation_inventory(
        repository, config["artifact_payload"], commit=authority_commit
    )
    _verify_qa1_pref1_result(repository)
    _verify_predecessor_leaves(repository)
    if _snapshot_store(repository) != (
        SEALED_STORE_LEAF_COUNT,
        SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("sealed store differs before QA2 authorization")
    if list(real_predecessor_fingerprints(repository)) != [
        _expected_fingerprint_receipt(width) for width in WIDTHS
    ]:
        _fail("restored predecessor fingerprints differ before authorization")
    return QA2Authority(authority_commit=authority_commit)


__all__ = [
    "ARTIFACT_ID",
    "AUTHORITY_PATH",
    "BASE_COMMIT",
    "CHANNEL_ORDER",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "DELTA_PATHS",
    "FORBIDDEN_RETRY3_GENERATION",
    "IMPLEMENTATION_INVENTORY",
    "MEMBER_KEY",
    "OUTPUT_NAMESPACE",
    "PROJECT_VERSION",
    "QA1_PREF1_ARTIFACT_ID",
    "QA1_PREF1_RESULT_PATH",
    "QA1_PREF1_RESULT_SHA256",
    "QA2Authority",
    "QA2AuthorityError",
    "REPLAYS",
    "RESULT_PATH",
    "RUNNER_PATH",
    "SELECTED_RETRIES",
    "STAGING_PREFIX",
    "STORE_PATH",
    "TARGET_PROTOCOL",
    "WIDTHS",
    "authorize",
    "build_prelaunch",
    "canonical_pretty",
    "expected_compact",
    "expected_raw_replay_receipt",
    "filled_implementation_pairs",
    "parse_config",
    "read_leaf",
    "real_predecessor_fingerprints",
    "require_output_absent",
    "validate_compact",
    "verify_implementation_inventory",
    "width_spec",
]
