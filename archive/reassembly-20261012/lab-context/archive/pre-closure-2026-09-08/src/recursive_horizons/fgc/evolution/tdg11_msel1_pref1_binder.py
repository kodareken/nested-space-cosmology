"""Independent TDG11 result binding, with a raw/store/Git-blind compact route.

Only ``bind_raw_result`` restores source data or constructs shadow families.
It replays the frozen proposer and guards as the implementation under test,
then uses separate reconstruction, localization, and protocol reductions.
No diagnostic selection, classification, or publication function is imported.
"""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha1, sha256
import json
import os
from pathlib import Path
import signal
import stat
import time
import tomllib

from recursive_horizons import evidence_io as io
from . import tdg11_msel1_contract as freeze
from . import tdg11_msel1_pref1_protocol as protocol


ARTIFACT_ID = "FGC-1-TDG11-MSEL1-PREF1"
CONFIG_PATH = "configs/fgc/fgc-1-tdg11-msel1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg11-msel1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg11-msel1-pref1.md"
AUTHORITY_COMMIT = "a04265794fe0a3266e34b50fcf034b9dcb462263"
AUTHORITY_PARENT = "dbfe92df40bf18fd287663e0b9997a812fe4ccb2"
FREEZE_CONFIG_SHA256 = (
    "4d468b3569b34eea77fe5613bfb733c1db05321118ae2dad266474317997b717"
)
FREEZE_RESULT_SHA256 = (
    "69acb67a760e0a8e43223e13e2fa6ed1c0234b7e26e9ea93d4f8a5366250a561"
)
SOURCE_TREE_SHA256 = "d33a52b6ca096884b826748484ac5671a8c3fa775a5ec1d48c978b8f319cd0c8"
SOURCE_TREE_BLOB_COUNT = 211
INDEPENDENT_COMPLETE_C_ID = "tdg11_msel1_pref1_independent_complete_c_v1"
RAW_COMPLETE_C_ID = "tdg11_rational_complete_c_dual_route_v1"

# Hashes observed after the single authorized diagnostic published both leaves.
RAW_MANIFEST_SHA256 = "26737b5dd3d437ec4e7da4844667565e2b9e255afb769afe4c73f052b2977911"
RAW_TERMINAL_SHA256 = "29c1f3a33757bd83c5a2b19d6b18e873c8fd0b5b4efae72865894adaefa5a221"

IMPLEMENTATION_PATHS = tuple(
    sorted(
        [
            "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_binder.py",
            "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_protocol.py",
            "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_reconstruction.py",
            "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_localization.py",
            "scripts/reproduce_fgc_tdg11_msel1_pref1.py",
            "tests/test_fgc_tdg11_msel1_pref1_binder.py",
            "tests/test_fgc_tdg11_msel1_pref1_protocol.py",
            "tests/test_fgc_tdg11_msel1_pref1_reconstruction.py",
            "tests/test_fgc_tdg11_msel1_pref1_localization.py",
        ]
    )
)
TASK_PATHS = frozenset(IMPLEMENTATION_PATHS) | {
    CONFIG_PATH,
    RESULT_PATH,
    OWNER_DOCUMENT,
    "README.md",
    "docs/active-code-map.md",
    "docs/claim-ledger.md",
    "docs/research-roadmap.md",
    "docs/fgc-runtime-matrix.md",
    "results/README.md",
    "configs/fgc/artifact-catalog.json",
    "scripts/build_artifact_catalog.py",
    "mk/current-foundation.mk",
    "mk/closed-live-targets.mk",
    "scripts/check_repo.py",
    "scripts/repo_checks/core.py",
    "scripts/repo_checks/catalog.py",
    "scripts/repo_checks/temporal_selection.py",
    "tests/test_check_repo_tdg11_msel1_pref1.py",
    "tests/test_phase_minus1_artifact_catalog.py",
    "tests/test_phase_minus1_make_routing.py",
}
LIMITS = {
    "maximum_wall_seconds": 14400,
    "maximum_raw_bytes": 16 * 1024 * 1024,
    "maximum_compact_bytes": 16 * 1024 * 1024,
    "maximum_guarded_families": 6,
    "maximum_shadow_proposals": 42,
    "maximum_stage_and_endpoint_records": 210,
    "maximum_shadow_rhs_calls": 252,
    "maximum_rational_bits": 32768,
}


class TDG11PREF1Error(RuntimeError):
    """The independent certificate could not be established."""


class TDG11PREF1ResourceStop(TDG11PREF1Error):
    """Binding exhausted its own prospective resource allowance."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise TDG11PREF1Error(message)


def _same(actual: object, expected: object, label: str) -> None:
    try:
        protocol.strict_equal(actual, expected, label)
    except protocol.PREF1ProtocolError as error:
        raise TDG11PREF1Error(str(error)) from error


def _root(repository: Path) -> Path:
    path = Path(repository)
    _require(
        path.is_absolute() and path.resolve(strict=True) == path,
        "repository path is not absolute and canonical",
    )
    _require(stat.S_ISDIR(path.lstat().st_mode), "repository root is not a directory")
    return path


def _freeze_config(raw: bytes) -> dict:
    _require(
        type(raw) is bytes and sha256(raw).hexdigest() == FREEZE_CONFIG_SHA256,
        "frozen config bytes differ",
    )
    return freeze.parse_config(raw)


def _measurement_pins() -> dict[str, str]:
    for value in (RAW_MANIFEST_SHA256, RAW_TERMINAL_SHA256):
        _require(
            type(value) is str,
            "the authorized measurement has not yet supplied raw pins",
        )
        protocol.digest(value, "published raw leaf")
    return {
        "manifest_sha256": RAW_MANIFEST_SHA256,
        "terminal_sha256": RAW_TERMINAL_SHA256,
    }


def expected_config(freeze_raw: bytes, implementation: list[dict[str, str]]) -> dict:
    source = _freeze_config(freeze_raw)
    _require(
        type(implementation) is list
        and len(implementation) == len(IMPLEMENTATION_PATHS),
        "PREF1 implementation inventory differs",
    )
    for item, path in zip(implementation, IMPLEMENTATION_PATHS, strict=True):
        protocol.exact_keys(item, {"path", "sha256"}, "PREF1 implementation")
        _same(item["path"], path, "implementation path")
        protocol.digest(item["sha256"], "implementation hash")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "authority_parent": AUTHORITY_PARENT,
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_result_sha256": FREEZE_RESULT_SHA256,
        "source_tree_sha256": SOURCE_TREE_SHA256,
        "source_tree_blob_count": SOURCE_TREE_BLOB_COUNT,
        "candidate_precedence": list(freeze.CANDIDATES),
        "selected_retries": [3, 4, 5],
        "restored_generations": [9, 10, 11],
        "complete_state_channels": list(freeze.TDG6_COMPLETE_STATE_CHANNELS),
        "raw": _measurement_pins(),
        "limits": dict(LIMITS),
        "environment": deepcopy(source["environment"]),
        "implementation": deepcopy(implementation),
    }


def parse_config(raw: bytes, freeze_raw: bytes) -> dict:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeError, tomllib.TOMLDecodeError) as error:
        raise TDG11PREF1Error("PREF1 config is not valid TOML") from error
    expected = expected_config(freeze_raw, value.get("implementation"))
    _same(value, expected, "PREF1 config")
    _require(render_config(expected) == raw, "PREF1 config is not canonical TOML")
    return expected


def render_config(value: dict) -> bytes:
    """Deterministic TOML for the closed scalar/table/table-array schema."""

    def literal(item):
        _require(type(item) in (str, int, bool, list), "unsupported config scalar")
        return json.dumps(item, ensure_ascii=True, separators=(", ", ": "))

    tables = {"raw", "limits", "environment", "implementation"}
    lines = [
        f"{key} = {literal(item)}" for key, item in value.items() if key not in tables
    ]
    for table in ("raw", "limits", "environment"):
        lines.extend(("", f"[{table}]"))
        lines.extend(
            f"{key} = {literal(item)}" for key, item in sorted(value[table].items())
        )
    for item in value["implementation"]:
        lines.extend(
            (
                "",
                "[[implementation]]",
                f"path = {literal(item['path'])}",
                f"sha256 = {literal(item['sha256'])}",
            )
        )
    return ("\n".join(lines) + "\n").encode("utf-8")


def _implementation_records(root: Path) -> list[dict[str, str]]:
    return [
        {"path": path, "sha256": sha256(io.read_regular_file(root, path)).hexdigest()}
        for path in IMPLEMENTATION_PATHS
    ]


def emit_config_bytes(repository: Path) -> bytes:
    root = _root(repository)
    return render_config(
        expected_config(
            io.read_regular_file(root, freeze.CONFIG_PATH),
            _implementation_records(root),
        )
    )


def _tracked_inputs(config_raw: bytes, root: Path) -> dict:
    freeze_raw = io.read_regular_file(root, freeze.CONFIG_PATH)
    source = _freeze_config(freeze_raw)
    config = parse_config(config_raw, freeze_raw)
    _same(_implementation_records(root), config["implementation"], "PREF1 code closure")
    frozen_result = io.read_regular_file(root, freeze.RESULT_PATH)
    _require(
        sha256(frozen_result).hexdigest() == FREEZE_RESULT_SHA256,
        "frozen compact result bytes differ",
    )
    # This route is compact and does not authorize or restore anything.
    freeze.validate_compact(
        freeze_raw, io.load_canonical_json(frozen_result.removesuffix(b"\n"))
    )
    for item in source["implementation"]:
        observed = sha256(io.read_regular_file(root, item["path"])).hexdigest()
        _same(observed, item["sha256"], "frozen implementation bytes")
    for path, expected in freeze.PINNED_REFERENCES.items():
        _same(
            sha256(io.read_regular_file(root, path)).hexdigest(),
            expected,
            "frozen compact predecessor/origin reference",
        )
    return config


def require_live_config(config_raw: bytes, repository: Path) -> None:
    """Bind a live construction to the on-disk config, not just cached bytes."""
    try:
        current = io.read_regular_file(_root(repository), CONFIG_PATH)
    except io.EvidenceIOError as error:
        raise TDG11PREF1Error("live PREF1 config became unreadable") from error
    _require(current == config_raw, "live PREF1 config changed during binding")


def _git(root: Path, operation: io.GitOperation) -> object:
    try:
        return io.git_read(root, operation)
    except io.EvidenceIOError as error:
        raise TDG11PREF1Error("immutable authority query failed") from error


def _source_listing_bytes(entries: tuple) -> bytes:
    """Reconstruct the exact nul-delimited ls-tree record, in Git tree order."""
    records = []
    for item in entries:
        _require(
            isinstance(item, io.TreeEntry)
            and item.kind == "blob"
            and item.mode in {"100644", "100755"},
            "nonregular frozen source object",
        )
        records.append(
            f"{item.mode} {item.kind} {item.object_id}\t{item.path}\0".encode("utf-8")
        )
    return b"".join(records)


def _authenticate_authority(root: Path) -> dict:
    """Authenticate the immutable freeze; deliberately do not require it at HEAD."""
    _same(
        _git(root, io.ResolveCommit(AUTHORITY_COMMIT)),
        AUTHORITY_COMMIT,
        "authority identity",
    )
    _require(
        _git(root, io.CommitParents(AUTHORITY_COMMIT)) == (AUTHORITY_PARENT,),
        "freeze is not its declared direct successor",
    )
    delta = _git(root, io.InspectDelta(AUTHORITY_COMMIT, AUTHORITY_PARENT))
    _require(
        tuple(sorted(item.path for item in delta)) == freeze.DELTA_PATHS
        and all(item.status in {"A", "M"} for item in delta),
        "freeze authority delta differs",
    )
    for path, expected in (
        (freeze.CONFIG_PATH, FREEZE_CONFIG_SHA256),
        (freeze.RESULT_PATH, FREEZE_RESULT_SHA256),
    ):
        blob = _git(root, io.ReadBlob(AUTHORITY_COMMIT, path))
        _same(sha256(blob).hexdigest(), expected, "authority compact object")
        _require(
            blob == io.read_regular_file(root, path),
            "live freeze differs from Git object",
        )
    image = _git(root, io.InspectTree(AUTHORITY_COMMIT, "src"))
    _require(len(image) == SOURCE_TREE_BLOB_COUNT, "frozen source object count differs")
    _same(
        sha256(_source_listing_bytes(image)).hexdigest(),
        SOURCE_TREE_SHA256,
        "frozen source tree commitment",
    )
    for entry in image:
        raw = io.read_regular_file(root, entry.path)
        git_blob = sha1(
            b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
        ).hexdigest()
        _same(git_blob, entry.object_id, f"instrument-under-test source {entry.path}")
        mode = "100755" if (root / entry.path).lstat().st_mode & 0o111 else "100644"
        _same(mode, entry.mode, "instrument-under-test source mode")
    pending = _git(root, io.InspectWorktree())
    _require(
        not pending.hidden_index_paths,
        "hidden Git index flags are not binding authority",
    )
    _require(
        all(item.path in TASK_PATHS for item in pending.changes),
        "unrelated pending files overlap the live binding image",
    )
    return {
        "commit": AUTHORITY_COMMIT,
        "sole_parent": AUTHORITY_PARENT,
        "exact_delta_paths": list(freeze.DELTA_PATHS),
        "source_tree_sha256": SOURCE_TREE_SHA256,
        "source_tree_blob_count": SOURCE_TREE_BLOB_COUNT,
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_result_sha256": FREEZE_RESULT_SHA256,
        "frozen_source_objects_match_live_instrument": True,
    }


def _authority_receipt() -> dict:
    """Fixed receipt shape for compact checking; it performs no Git query."""
    return {
        "commit": AUTHORITY_COMMIT,
        "sole_parent": AUTHORITY_PARENT,
        "exact_delta_paths": list(freeze.DELTA_PATHS),
        "source_tree_sha256": SOURCE_TREE_SHA256,
        "source_tree_blob_count": SOURCE_TREE_BLOB_COUNT,
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_result_sha256": FREEZE_RESULT_SHA256,
        "frozen_source_objects_match_live_instrument": True,
    }


def _directory_identity(path: Path) -> tuple:
    info = path.lstat()
    _require(
        stat.S_ISDIR(info.st_mode) and path.resolve(strict=True) == path,
        "raw namespace is not a canonical regular directory",
    )
    return info.st_dev, info.st_ino, info.st_mode, info.st_mtime_ns, info.st_ctime_ns


def _raw_snapshot(root: Path, config: dict) -> tuple[dict[str, bytes], dict, dict]:
    namespace = root / freeze.OUTPUT_NAMESPACE
    before = _directory_identity(namespace)
    with os.scandir(namespace) as entries:
        names = {entry.name for entry in entries}
    _require(
        names == {"manifest.json", "terminal.json"},
        "raw namespace is partial or has extra leaves",
    )
    blobs = {
        name: io.read_regular_file(
            root,
            f"{freeze.OUTPUT_NAMESPACE}/{name}",
            max_bytes=LIMITS["maximum_raw_bytes"],
        )
        for name in sorted(names)
    }
    _require(
        sum(len(raw) for raw in blobs.values()) <= LIMITS["maximum_raw_bytes"],
        "raw bundle exceeds the frozen bound",
    )
    with os.scandir(namespace) as entries:
        _require(
            {entry.name for entry in entries} == names, "raw namespace names changed"
        )
    _require(
        _directory_identity(namespace) == before, "raw namespace changed during capture"
    )
    for leaf, key in (
        ("manifest.json", "manifest_sha256"),
        ("terminal.json", "terminal_sha256"),
    ):
        _same(sha256(blobs[leaf]).hexdigest(), config["raw"][key], f"published {leaf}")
    manifest = io.load_canonical_json(blobs["manifest.json"])
    terminal = io.load_canonical_json(blobs["terminal.json"])
    protocol.validate_raw_manifest(
        manifest,
        authority_commit=AUTHORITY_COMMIT,
        config_sha256=FREEZE_CONFIG_SHA256,
        freeze_sha256=FREEZE_RESULT_SHA256,
        terminal_sha256=config["raw"]["terminal_sha256"],
    )
    protocol.validate_raw_terminal(
        terminal,
        authority_commit=AUTHORITY_COMMIT,
        config_sha256=FREEZE_CONFIG_SHA256,
        freeze_sha256=FREEZE_RESULT_SHA256,
        environment=config["environment"],
    )
    return blobs, manifest, terminal


def _expected_manifest(config: dict) -> dict:
    return {
        "schema": protocol.RAW_SCHEMA,
        "runner_id": protocol.RUNNER_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "config_sha256": FREEZE_CONFIG_SHA256,
        "freeze_sha256": FREEZE_RESULT_SHA256,
        "terminal_sha256": config["raw"]["terminal_sha256"],
        "leaf_names": ["manifest.json", "terminal.json"],
        "endpoint_serialized": False,
        "historical_store_written": False,
    }


def _conclusion(selection: dict) -> dict:
    classification = selection["classification"]
    selected = classification == protocol.SELECTED
    if selected:
        branch = "candidate_qualified_only_for_separate_TDG11_IMP1"
    elif classification == protocol.NONPASS:
        branch = "tested_numerical_branch_rejected_requires_new_solver_formulation_plan"
    else:
        branch = "bounded_inconclusive_requires_new_prospective_freeze"
    return {
        "branch": branch,
        "selected_candidate": selection["selected_candidate"],
        "licenses_only_separate_TDG11_IMP1": selected,
        "production_method_selected_or_implemented": False,
        "diagnostic_endpoint_may_be_adopted": False,
        "campaign_may_be_resumed": False,
        "GR0_calibration_completed": False,
        "SGBL_candidate_or_holdout_execution_opened": False,
        "FGC_QR_action_rejected": False,
        "global_PDE_error_enclosure": False,
        "physical_claim": False,
        "publication_authorized": False,
    }


def _channel_record(candidate: str, channel: str, data: object, localize) -> dict:
    """Build a record solely from fresh independent finite-data reductions."""
    rounds = protocol.three_magnitudes(
        protocol.wire(data.accumulation_bounds), "independent accumulation bounds"
    )
    if candidate == freeze.CANDIDATES[2]:
        e0, e1, e2 = protocol.three_magnitudes(
            protocol.wire(data.embedded_defect_bounds), "independent embedded defects"
        )
        record = protocol.wire(
            {
                "channel": channel,
                "estimators": (e0, e1, e2),
                "roundoff_bounds": rounds,
                "contraction_01": protocol.decision_from_bounds(e0, e0, e1, e1),
                "contraction_12": protocol.decision_from_bounds(e1, e1, e2, e2),
                "public_fine_debit": e2 + rounds[2],
                "global_PDE_enclosure": False,
            }
        )
    else:
        rows = (
            data.corrected_rows if candidate == freeze.CANDIDATES[0] else data.raw_rows
        )
        evidence = localize(
            rows,
            expected_row_count=freeze.OWNED_ROW_COUNT,
            maximum_candidates_D01=freeze.RESOURCES[
                "per_channel_maximum_candidates_D01"
            ],
            maximum_candidates_D12=freeze.RESOURCES[
                "per_channel_maximum_candidates_D12"
            ],
            refinement_depth=freeze.RESOURCES["primary_refinement_depth"],
        )
        exact = protocol.wire(evidence)
        declared_wire = _raw_channel({"exact_complete_C": exact}, candidate)[
            "exact_complete_C"
        ]
        validated = protocol.validate_complete_c(
            declared_wire, expected_rows=freeze.OWNED_ROW_COUNT
        )
        debit = validated["d12"][1]
        record = {"channel": channel, "exact_complete_C": exact}
        if candidate == freeze.CANDIDATES[0]:
            record["roundoff_bounds"] = protocol.wire(rounds)
            debit += rounds[2]
        record["public_fine_debit"] = protocol.wire(debit)
    protocol.classify_channel(
        _raw_channel(record, candidate), candidate, expected_rows=freeze.OWNED_ROW_COUNT
    )
    return record


def _raw_channel(value: dict, candidate: str) -> dict:
    """Map only an adapter provenance ID when rebuilding the raw commitment.

    Saved PREF1 records retain the independent adapter ID. Numeric values,
    classifications, counts, and stream hashes are never translated.
    """
    _require(type(value) is dict, "independent channel must be a record")
    result = deepcopy(value)
    if candidate != freeze.CANDIDATES[2]:
        evidence = result.get("exact_complete_C")
        _require(type(evidence) is dict, "independent complete-C evidence is absent")
        _same(
            evidence.get("evaluator_id"),
            INDEPENDENT_COMPLETE_C_ID,
            "independent adapter identity",
        )
        evidence["evaluator_id"] = RAW_COMPLETE_C_ID
    return result


def _raw_widths(value: object) -> list[dict]:
    _require(type(value) is list and len(value) == 3, "independent width slots differ")
    result = deepcopy(value)
    for width in result:
        protocol.exact_keys(
            width,
            {
                "retry",
                "generation",
                "width_hex",
                "restored",
                "families",
                "baseline",
                "candidates",
            },
            "independent width",
        )
        protocol.exact_keys(
            width["candidates"], set(freeze.CANDIDATES), "independent candidate slots"
        )
        pairs = [("baseline", width["baseline"]), *width["candidates"].items()]
        for candidate, group in pairs:
            protocol.exact_keys(
                group, {"status", "channels", "stop"}, "independent group"
            )
            _require(
                type(group["channels"]) is list, "independent channels are not a list"
            )
            group["channels"] = [
                _raw_channel(record, candidate) for record in group["channels"]
            ]
    return result


def _group_pairs(width: dict, arithmetic: str) -> list[tuple[str, dict]]:
    if arithmetic == freeze.ORIGINAL_ARITHMETIC_ID:
        return [
            ("baseline", width["baseline"]),
            (freeze.CANDIDATES[0], width["candidates"][freeze.CANDIDATES[0]]),
            (freeze.CANDIDATES[2], width["candidates"][freeze.CANDIDATES[2]]),
        ]
    return [(freeze.CANDIDATES[1], width["candidates"][freeze.CANDIDATES[1]])]


def _counts_from_family(family: object) -> dict[str, int]:
    return {
        "accepted_source_prechecks": family.accepted_source_precheck_count,
        "stage_and_endpoint_records": family.stage_record_count,
        "rhs_calls": family.rhs_call_count,
    }


def _stop_scope(widths: list[dict]) -> list[dict]:
    """Deterministic disclosure of which stopped work is not a numerical replay."""
    rows = []
    for arithmetic in (freeze.ORIGINAL_ARITHMETIC_ID, freeze.COMPENSATED_ARITHMETIC_ID):
        for width in widths:
            family = width["families"][arithmetic]
            if family is None:
                continue
            if family["status"] != "complete":
                rows.append(
                    {
                        "retry": width["retry"],
                        "arithmetic_id": arithmetic,
                        "scope": "guarded_family",
                        "owner": family["stop"]["owner"],
                        "code": family["stop"]["code"],
                        "verification": "premise_reproduced"
                        if family["status"] == "premise_stop"
                        else "resource_receipt_authenticated_without_reexecution",
                        "scientific_nonpass_inferred_from_stop": False,
                    }
                )
                continue
            for candidate, group in _group_pairs(width, arithmetic):
                if group["status"] not in {"premise_stop", "resource_stop"}:
                    continue
                stop = group["stop"]
                rows.append(
                    {
                        "retry": width["retry"],
                        "arithmetic_id": arithmetic,
                        "scope": candidate,
                        "owner": stop["owner"],
                        "code": stop["code"],
                        "verification": "localizer_stop_reproduced"
                        if stop["owner"] == "exact_localizer"
                        else "resource_receipt_authenticated_without_reexecution",
                        "scientific_nonpass_inferred_from_stop": False,
                    }
                )
    return rows


def _replay_accounting(widths: list[dict], static_shells: int) -> dict:
    """Counts follow replayed complete/premise families, not raw resource work."""
    complete = 0
    premise_stops = 0
    resource_stops = 0
    counts = {
        "accepted_source_prechecks": 0,
        "stage_and_endpoint_records": 0,
        "rhs_calls": 0,
    }
    channel_counts = {"baseline": 0, **dict.fromkeys(freeze.CANDIDATES, 0)}
    for width in widths:
        for arithmetic, family in width["families"].items():
            if family is None:
                continue
            if family["status"] == "resource_stop":
                resource_stops += 1
                continue
            if family["status"] == "complete":
                complete += 1
            else:
                premise_stops += 1
            for key in counts:
                counts[key] += family[key]
            for candidate, group in _group_pairs(width, arithmetic):
                channel_counts[candidate] += len(group["channels"])
    return {
        "static_shells_rebuilt": static_shells,
        "restored_retries": [
            width["retry"] for width in widths if width["restored"] is not None
        ],
        "complete_guarded_families_recomputed": complete,
        "complete_shadow_proposals_recomputed": 7 * complete,
        "guarded_premise_stops_reproduced": premise_stops,
        "guarded_resource_stops_authenticated_only": resource_stops,
        "replayed_work": counts,
        "channel_records_recomputed": channel_counts,
        "typed_stop_scope": _stop_scope(widths),
        "independent_reconstruction_and_localization": True,
        "independent_complete_c_evaluator_id": INDEPENDENT_COMPLETE_C_ID,
        "raw_commitment_complete_c_evaluator_id": RAW_COMPLETE_C_ID,
        "runner_classification_or_publication_imported": False,
        "complete_channel_record_bytes_match_raw_commitment": True,
        "source_store_unchanged": True,
        "raw_leaves_unchanged": True,
        "accepted_endpoint_serialized_or_adopted": False,
    }


def _compact_payload(
    config_raw: bytes, config: dict, raw_terminal: dict, widths: list[dict]
) -> dict:
    header = {
        key: deepcopy(value) for key, value in raw_terminal.items() if key != "widths"
    }
    selection = protocol.reduce_widths(
        _raw_widths(widths), expected_rows=freeze.OWNED_ROW_COUNT
    )
    return {
        "artifact_id": ARTIFACT_ID,
        "schema_version": 1,
        "classification": selection["classification"],
        "selected_candidate": selection["selected_candidate"],
        "config_sha256": sha256(config_raw).hexdigest(),
        "authority": _authority_receipt(),
        "raw_commitment": {**config["raw"], "terminal_header": header},
        "independent_replay": {
            "widths": widths,
            **_replay_accounting(widths, header["accounting"]["static_shells"]),
        },
        "conclusion": _conclusion(selection),
        "nonclaims": freeze.nonclaims(),
    }


def validate_compact_result(
    config_raw: bytes, result_raw: bytes, repository: Path
) -> dict:
    """Validate tracked bytes and the complete raw commitment without live data.

    There is no Git query, environment observation, raw namespace traversal,
    store access, shadow construction, or localizer execution on this route.
    """
    root = _root(repository)
    config = _tracked_inputs(config_raw, root)
    _require(
        type(result_raw) is bytes
        and len(result_raw) <= LIMITS["maximum_compact_bytes"],
        "compact result byte bound differs",
    )
    result = io.load_canonical_json(result_raw.removesuffix(b"\n"))
    protocol.exact_keys(
        result,
        {
            "artifact_id",
            "schema_version",
            "classification",
            "selected_candidate",
            "config_sha256",
            "authority",
            "raw_commitment",
            "independent_replay",
            "conclusion",
            "nonclaims",
        },
        "compact result",
    )
    commitment = protocol.exact_keys(
        result["raw_commitment"],
        {"manifest_sha256", "terminal_sha256", "terminal_header"},
        "raw commitment",
    )
    for key in ("manifest_sha256", "terminal_sha256"):
        _same(commitment[key], config["raw"][key], f"raw commitment {key}")
    _require(
        type(result["independent_replay"]) is dict, "independent replay is not a record"
    )
    widths = result["independent_replay"].get("widths")
    raw_widths = _raw_widths(widths)
    header = commitment["terminal_header"]
    _require(
        type(header) is dict and "widths" not in header,
        "raw header embeds duplicate widths",
    )
    terminal = {**header, "widths": raw_widths}
    protocol.validate_raw_terminal(
        terminal,
        authority_commit=AUTHORITY_COMMIT,
        config_sha256=FREEZE_CONFIG_SHA256,
        freeze_sha256=FREEZE_RESULT_SHA256,
        environment=config["environment"],
    )
    _same(
        sha256(io.canonical_json_bytes(terminal)).hexdigest(),
        config["raw"]["terminal_sha256"],
        "reconstructed full raw terminal commitment",
    )
    _same(
        sha256(io.canonical_json_bytes(_expected_manifest(config))).hexdigest(),
        config["raw"]["manifest_sha256"],
        "reconstructed raw manifest commitment",
    )
    expected = _compact_payload(config_raw, config, terminal, widths)
    _same(result, expected, "independent compact certificate")
    return result


def _replay_independently(root: Path, terminal: dict, progress) -> list[dict]:
    # These are deliberately live-only imports. The proposer/RHS/guards are
    # the frozen instrument under test; the reductions below are independent.
    from . import tdg11_msel1_authority as authority
    from .numerical_engine import PRIMARY_METHOD
    from .proto19_gr0_static_factory import build_static_gr0_shells
    from .tdg11_msel1_runtime import build_guarded_shadow_family, TDG11ShadowPremiseStop
    from .tdg11_msel1_pref1_reconstruction import (
        validate_recorded_family_independently,
        reconstruct_channel_independently,
    )
    from .tdg11_msel1_pref1_localization import (
        assess_complete_c_independently,
        IndependentCompleteCResourceExhausted,
        IndependentCompleteCRouteDisagreement,
    )

    widths = deepcopy(terminal["widths"])
    templates = None
    if terminal["accounting"]["static_shells"] == 6:
        templates = build_static_gr0_shells(root)
        _require(len(templates) == 6, "independent static factory differs")
    restored = {}
    for width in widths:
        if width["restored"] is None:
            continue
        _require(templates is not None, "restoration lacks the static factory")
        item = authority.restore_predecessor(root, width["retry"], templates=templates)
        receipt = {
            "checkpoint_sha256": item.checkpoint_sha256,
            "descriptor_sha256": item.descriptor_sha256,
            "state_sha256": item.fingerprint["state_sha256"],
            "fingerprint_sha256": sha256(
                io.canonical_json_bytes(protocol.wire(item.fingerprint))
            ).hexdigest(),
        }
        _same(receipt, width["restored"], "independent restoration receipt")
        width["restored"] = receipt
        restored[width["retry"]] = item

    for arithmetic in (freeze.ORIGINAL_ARITHMETIC_ID, freeze.COMPENSATED_ARITHMETIC_ID):
        for width, spec in zip(widths, freeze.replay_specs(), strict=True):
            raw_family = width["families"][arithmetic]
            if raw_family is None:
                continue
            if raw_family["status"] == "resource_stop":
                # Replaying a wall-clock interruption would be a new timing
                # experiment. Its receipt is authenticated, never a nonpass.
                progress(
                    {
                        "retry": spec["retry"],
                        "arithmetic": arithmetic,
                        "phase": "raw_resource_receipt_authenticated_only",
                    }
                )
                continue
            item = restored[spec["retry"]]
            member = item.member
            progress(
                {
                    "retry": spec["retry"],
                    "arithmetic": arithmetic,
                    "phase": "independent_guarded_family",
                }
            )
            try:
                family = build_guarded_shadow_family(
                    method=PRIMARY_METHOD,
                    arithmetic_id=arithmetic,
                    time=member.time,
                    step_size=float.fromhex(spec["width_hex"]),
                    state=member.state,
                    rhs=member.operator,
                    projector=member.projector,
                    transaction=member.transaction,
                    tracers=member.tracers,
                    coordinates=member.initial.grid.coordinates,
                    previous_step_index=member.step_index,
                    previous_transaction_serial=member.transaction_serial,
                )
            except TDG11ShadowPremiseStop as error:
                if isinstance(error.cause, TDG11PREF1ResourceStop):
                    raise error.cause
                _same(
                    raw_family["status"], "premise_stop", "independent guarded premise"
                )
                _same(
                    raw_family["stop"]["owner"],
                    "guarded_shadow",
                    "independent premise owner",
                )
                code = (
                    "source_retry"
                    if error.source_retry is not None
                    else type(error.cause).__name__
                )
                _same(raw_family["stop"]["code"], code, "independent premise cause")
                counts = _counts_from_family(error)
                _same(
                    counts,
                    {key: raw_family[key] for key in counts},
                    "independent premise prefix",
                )
                _same(
                    protocol.wire(
                        authority._member_fingerprint(member, item.descriptor_sha256)
                    ),
                    protocol.wire(item.fingerprint),
                    "restored member after premise replay",
                )
                progress(
                    {
                        "retry": spec["retry"],
                        "arithmetic": arithmetic,
                        "phase": "guarded_premise_reproduced",
                        "code": code,
                        "cause": str(error.cause)[:640],
                    }
                )
                continue
            _same(
                raw_family["status"], "complete", "guarded premise was not reproduced"
            )
            paths = tuple(
                tuple(attempt.proposal for attempt in path.attempts)
                for path in family.paths
            )
            recorded = validate_recorded_family_independently(
                paths,
                method=PRIMARY_METHOD,
                arithmetic_id=arithmetic,
                owned_row_count=freeze.OWNED_ROW_COUNT,
            )
            counts = _counts_from_family(family)
            new_family = {
                "status": "complete",
                "family_sha256": recorded.family_sha256,
                **counts,
            }
            _same(
                new_family,
                raw_family,
                "independent recorded family and work accounting",
            )
            width["families"][arithmetic] = new_family
            for candidate, group in _group_pairs(width, arithmetic):
                raw_channels = group["channels"]
                rebuilt = []
                for index, raw_channel in enumerate(raw_channels):
                    name = freeze.TDG6_COMPLETE_STATE_CHANNELS[index]
                    data = reconstruct_channel_independently(recorded, name)
                    record = _channel_record(
                        candidate, name, data, assess_complete_c_independently
                    )
                    _same(
                        _raw_channel(record, candidate),
                        raw_channel,
                        "independent channel record",
                    )
                    rebuilt.append(record)
                    progress(
                        {
                            "retry": spec["retry"],
                            "candidate": candidate,
                            "channel": name,
                            "phase": "independent_channel_matched",
                        }
                    )
                group["channels"] = rebuilt
                if (
                    group["status"] == "resource_stop"
                    and group["stop"]["owner"] == "exact_localizer"
                ):
                    _require(
                        candidate != freeze.CANDIDATES[2]
                        and len(rebuilt) < len(freeze.TDG6_COMPLETE_STATE_CHANNELS),
                        "impossible localizer stop scope",
                    )
                    name = freeze.TDG6_COMPLETE_STATE_CHANNELS[len(rebuilt)]
                    data = reconstruct_channel_independently(recorded, name)
                    try:
                        _channel_record(
                            candidate, name, data, assess_complete_c_independently
                        )
                    except (
                        IndependentCompleteCResourceExhausted,
                        IndependentCompleteCRouteDisagreement,
                    ) as error:
                        _same(
                            error.evidence.reason,
                            group["stop"]["code"],
                            "independent localizer stop",
                        )
                    else:
                        raise TDG11PREF1Error(
                            "the raw localizer stop was not reproduced"
                        )
            _same(
                protocol.wire(
                    authority._member_fingerprint(member, item.descriptor_sha256)
                ),
                protocol.wire(item.fingerprint),
                "restored member after complete replay",
            )
    _same(
        _raw_widths(widths),
        terminal["widths"],
        "all independently reconstructed width records",
    )
    return widths


def bind_raw_result(config_raw: bytes, repository: Path, *, progress=None) -> dict:
    """Explicit read-only reconstruction; returns a result but publishes nothing."""
    root = _root(repository)
    require_live_config(config_raw, root)
    config = _tracked_inputs(config_raw, root)
    authority_receipt = _authenticate_authority(root)
    _same(authority_receipt, _authority_receipt(), "live authority receipt")
    blobs, _manifest, terminal = _raw_snapshot(root, config)
    from . import tdg11_msel1_authority as authority

    _same(
        authority.observe_environment(),
        config["environment"],
        "independent replay environment",
    )
    source_before = authority.snapshot_store(root)
    _require(
        source_before == (freeze.STORE_LEAF_COUNT, freeze.STORE_SHA256),
        "source store is not the frozen 115-leaf image",
    )
    authority.inspect_predecessors(root)
    started = time.monotonic()
    callback = (lambda _event: None) if progress is None else progress

    def report(event):
        if time.monotonic() - started > LIMITS["maximum_wall_seconds"]:
            raise TDG11PREF1ResourceStop("independent replay wall budget exhausted")
        callback(event)

    def alarm(_signum, _frame):
        raise TDG11PREF1ResourceStop("independent replay wall budget exhausted")

    _require(
        signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0),
        "another process alarm owns the timer",
    )
    previous_handler = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, LIMITS["maximum_wall_seconds"])
    try:
        widths = _replay_independently(root, terminal, report)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous_handler)
    _require(
        authority.snapshot_store(root) == source_before,
        "source store changed during binding",
    )
    _require(
        _raw_snapshot(root, config)[0] == blobs,
        "published raw leaves changed during binding",
    )
    _same(
        _authenticate_authority(root),
        authority_receipt,
        "authority after independent replay",
    )
    _same(
        authority.observe_environment(),
        config["environment"],
        "environment after independent replay",
    )
    require_live_config(config_raw, root)
    _tracked_inputs(config_raw, root)
    result = _compact_payload(config_raw, config, terminal, widths)
    raw = io.canonical_json_bytes(result) + b"\n"
    validate_compact_result(config_raw, raw, root)
    return result


def verify_compact(repository: Path) -> dict:
    root = _root(repository)
    return validate_compact_result(
        io.read_regular_file(root, CONFIG_PATH),
        io.read_regular_file(root, RESULT_PATH),
        root,
    )
