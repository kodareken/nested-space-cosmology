"""Artifact-catalog and generated result-index validation."""

from __future__ import annotations

from ._shared import *  # noqa: F401,F403
from ._shared import (
    Path,
    REPOSITORY,
    _probe_unique_regular_leaf,
    _read_unique_regular_bytes,
    _reject_duplicate_json_pairs,
    json,
    re,
    sha256,
    subprocess,
    sys,
    tomllib,
)


_CATALOG_PATH = "configs/fgc/artifact-catalog.json"
_CATALOG_SCHEMA = "FGC-artifact-catalog-v2"
_CATALOG_CURRENT = "FGC-1-PRO20-EV1-PREF1"
_CATALOG_NEXT = "FGC-1-PRO20-EV1-RSRC1"
_RSRC1_CLASSIFICATION = "prospective_resource_isolation_core_no_authority"
_RSRC1_SOURCES = {
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py",
}
_RSRC1_REPRODUCERS = {"scripts/run_fgc_pro20_rsrc1_child.py"}
_RSRC1_TESTS = {
    "tests/test_fgc_pro20_rsrc1_attempt.py",
    "tests/test_fgc_pro20_rsrc1_isolation.py",
    "tests/test_fgc_pro20_rsrc1_member.py",
    "tests/test_fgc_pro20_rsrc1_seed.py",
}
_PRO20_FRZ1 = "FGC-1-PRO20-EV1-FRZ1"
_REC1_PREF1 = "FGC-1-HLT17-SRCQ1-REC1-PREF1"
_PRO20_SOURCES = {
    "src/recursive_horizons/fgc/evolution/pro20_ev1_authority.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_protocol.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_runtime.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_store.py",
}
_PRO20_TESTS = {
    "tests/test_fgc_pro20_ev1_authority.py",
    "tests/test_fgc_pro20_ev1_protocol.py",
    "tests/test_fgc_pro20_ev1_runtime.py",
    "tests/test_fgc_pro20_ev1_store.py",
}
_PRO20_REPRODUCERS = {
    "scripts/reproduce_fgc_pro20_ev1_frz1.py",
    "scripts/run_fgc_pro20_ev1.py",
}
_PRO20_OWNER = "docs/fgc-pro20-ev1-frz1.md"
_PRO20_CONFIG = "configs/fgc/fgc-1-pro20-ev1-frz1.toml"
_IMP1_ARTIFACT = "FGC-1-TDG11-IMP1"
_REC1_FRZ1 = "FGC-1-HLT17-SRCQ1-REC1-FRZ1"
_PARALLEL_INSTRUMENTS = {
    "FGC-1-SGB1-CTL1": (
        "partial_outcome_blind_sgb_source_and_synthetic_runtime_instruments",
        {"src/recursive_horizons/fgc/sgb1_ctl1_constraints.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_source.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_family.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_center.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_principal.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_interval_health.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_runtime.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_admission.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_initial_health.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_controls.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_cone.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_continuity.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_family_principal.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_cell_admission.py",
         "src/recursive_horizons/fgc/sgb1_ctl1_local_symmetrizer.py"},
        {"docs/fgc-sgb1-ctl1.md",
         "docs/fgc-sgb1-ctl1-admission.md",
         "docs/fgc-sgb1-ctl1-initial-health.md",
         "docs/fgc-sgb1-ctl1-cone.md",
         "docs/fgc-sgb1-ctl1-continuity.md",
         "docs/fgc-sgb1-ctl1-trap-refinement.md",
         "docs/fgc-sgb1-ctl1-family-principal.md",
         "docs/fgc-sgb1-ctl1-cell-admission.md",
         "docs/fgc-sgb1-ctl1-local-symmetrizer.md"},
    ),
    "FGC-1-SGB1-CTL1-TRAP-REFINEMENT2": (
        "finite_product_box_no_trap_refinement_nonpass_no_health",
        {"src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement2.py"},
        {"docs/fgc-sgb1-ctl1-trap-refinement2.md"},
    ),
    "FGC-1-SGB1-CTL1-TRAP-TAYLOR": (
        "order4_taylor_no_trap_enclosure_nonpass_no_health",
        {"src/recursive_horizons/fgc/sgb1_ctl1_trap_taylor.py"},
        {"docs/fgc-sgb1-ctl1-trap-taylor.md"},
    ),
    "FGC-1-SGB1-CTL1-TRAP-BARRIER": (
        "barrier_nonpass_not_trajectory_trapping",
        {"src/recursive_horizons/fgc/sgb1_ctl1_trap_barrier.py"},
        {"docs/fgc-sgb1-ctl1-trap-barrier.md"},
    ),
    "FGC-1-SGB1-CTL1-SOL1": (
        "orbit_local_picard_lindelof_interval_inconclusive_no_health",
        {"src/recursive_horizons/fgc/sgb1_ctl1_sol1.py"},
        {"docs/fgc-sgb1-ctl1-sol1.md"},
    ),
    "FGC-1-DEF1-STAB1": (
        "partial_outcome_blind_error_map_and_margin_algebra",
        {"src/recursive_horizons/fgc/def1_stab1.py",
         "src/recursive_horizons/fgc/def1_geometry_error.py",
         "src/recursive_horizons/fgc/def1_stab1_providers.py",
         "src/recursive_horizons/fgc/def1_stab1_qualification.py",
         "src/recursive_horizons/fgc/def1_stab1_freeze_contract.py"},
        {"docs/fgc-def1-stab1.md", "docs/fgc-def1-stab1-geometry.md"},
    ),
    "FGC-1-TDG11-C1R1": (
        "partial_exact_c1_representation_no_production_authority",
        {"src/recursive_horizons/fgc/evolution/tdg11_c1r1_ring.py",
         "src/recursive_horizons/fgc/evolution/tdg11_c1r1_enclosure.py",
         "src/recursive_horizons/fgc/evolution/tdg11_c1r1_runtime.py"},
        {"docs/fgc-tdg11-c1r1.md"},
    ),
    "FGC-1-HLT17-MON17": (
        "partial_c1r1_finite_state_runtime_and_atomic_store",
        {"src/recursive_horizons/fgc/evolution/hlt17_admission_runtime.py",
         "src/recursive_horizons/fgc/evolution/hlt17_runtime_member.py",
         "src/recursive_horizons/fgc/evolution/hlt17_imp1_cursor.py",
         "src/recursive_horizons/fgc/evolution/hlt17_imp1_bridge.py",
         "src/recursive_horizons/fgc/evolution/hlt17_member_codec.py",
         "src/recursive_horizons/fgc/evolution/hlt17_member_checkpoint.py",
         "src/recursive_horizons/fgc/evolution/protocol_v19.py",
         "src/recursive_horizons/fgc/evolution/hlt17_campaign_store.py"},
        {"docs/fgc-hlt17-mon17.md"},
    ),
}
_QA2_PREF1 = "FGC-1-TDG10-QA2-PREF1"
_MSEL1_FRZ1 = "FGC-1-TDG11-MSEL1-FRZ1"
_PREF1_ARTIFACT = "FGC-1-TDG11-MSEL1-PREF1"
_PREF1_SOURCES = (
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_binder.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_protocol.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_reconstruction.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_localization.py",
)
_PREF1_CONFIG = "configs/fgc/fgc-1-tdg11-msel1-pref1.toml"
_PREF1_RESULT = "results/fgc-1-tdg11-msel1-pref1.json"
_PREF1_OWNER = "docs/fgc-tdg11-msel1-pref1.md"
_PREF1_REPRODUCER = "scripts/reproduce_fgc_tdg11_msel1_pref1.py"
_PREF1_TESTS = (
    "tests/test_fgc_tdg11_msel1_pref1_binder.py",
    "tests/test_fgc_tdg11_msel1_pref1_protocol.py",
    "tests/test_fgc_tdg11_msel1_pref1_reconstruction.py",
    "tests/test_fgc_tdg11_msel1_pref1_localization.py",
    "tests/test_check_repo_tdg11_msel1_pref1.py",
)
_IMP1_SOURCES = (
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_enclosure.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_ledger.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_runtime.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_qualification.py",
)
_IMP1_CONFIG = "configs/fgc/fgc-1-tdg11-imp1.toml"
_IMP1_RESULT = "results/fgc-1-tdg11-imp1.json"
_IMP1_OWNER = "docs/fgc-tdg11-imp1.md"
_IMP1_REPRODUCER = "scripts/reproduce_fgc_tdg11_imp1.py"
_IMP1_TESTS = (
    "tests/test_fgc_tdg11_imp1_enclosure.py",
    "tests/test_fgc_tdg11_imp1_ledger.py",
    "tests/test_fgc_tdg11_imp1_runtime.py",
    "tests/test_fgc_tdg11_imp1_retry_limits.py",
    "tests/test_fgc_tdg11_imp1_reproduction.py",
    "tests/test_check_repo_tdg11_imp1.py",
)
_DEF1_FRZ1 = "FGC-1-DEF1-STAB1-FRZ1"
_DEF1_FRZ1_CONFIG = "configs/fgc/fgc-1-def1-stab1-frz1.toml"
_DEF1_FRZ1_RESULT = "results/fgc-1-def1-stab1-frz1.json"
_DEF1_FRZ1_OWNER = "docs/fgc-def1-stab1-frz1.md"
_DEF1_FRZ1_REPRODUCER = "scripts/reproduce_fgc_def1_stab1_frz1.py"
_DEF1_FRZ1_SOURCE = "src/recursive_horizons/fgc/def1_stab1_frz1_certificate.py"
_DEF1_FRZ1_TESTS = (
    "tests/test_fgc_def1_stab1_frz1_certificate.py",
    "tests/test_check_repo_def1_stab1_frz1.py",
)
_DEF1_FRZ1_CLASSIFICATION = (
    "candidate_blind_conversion_error_map_instrument_freeze_no_trajectory"
)
_DEF1_PREF1 = "FGC-1-DEF1-STAB1-PREF1"
_DEF1_PREF1_CONFIG = "configs/fgc/fgc-1-def1-stab1-pref1.toml"
_DEF1_PREF1_RESULT = "results/fgc-1-def1-stab1-pref1.json"
_DEF1_PREF1_OWNER = "docs/fgc-def1-stab1-pref1.md"
_DEF1_PREF1_REPRODUCER = "scripts/reproduce_fgc_def1_stab1_pref1.py"
_DEF1_PREF1_SOURCE = "src/recursive_horizons/fgc/def1_stab1_pref1_binder.py"
_DEF1_PREF1_TESTS = (
    "tests/test_fgc_def1_stab1_pref1_binder.py",
    "tests/test_check_repo_def1_stab1_pref1.py",
)
_DEF1_PREF1_CLASSIFICATION = (
    "independently_bound_candidate_blind_error_map_readiness_only_no_trajectory"
)
_SOL1_FRZ1 = "FGC-1-SGB1-CTL1-SOL1-FRZ1"
_SOL1_FRZ1_CONFIG = "configs/fgc/fgc-1-sgb1-ctl1-sol1-frz1.toml"
_SOL1_FRZ1_RESULT = "results/fgc-1-sgb1-ctl1-sol1-frz1.json"
_SOL1_FRZ1_OWNER = "docs/fgc-sgb1-ctl1-sol1-frz1.md"
_SOL1_FRZ1_REPRODUCER = "scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py"
_SOL1_FRZ1_SOURCE = "src/recursive_horizons/fgc/sgb1_ctl1_sol1_frz1_certificate.py"
_SOL1_FRZ1_TESTS = (
    "tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py",
    "tests/test_check_repo_sgb1_ctl1_sol1_frz1.py",
)
_SOL1_FRZ1_CLASSIFICATION = (
    "prospective_sol1_frozen_nominal_interval_inconclusive_no_health"
)
_SOL1_PREF1 = "FGC-1-SGB1-CTL1-SOL1-PREF1"
_SOL1_PREF1_CONFIG = "configs/fgc/fgc-1-sgb1-ctl1-sol1-pref1.toml"
_SOL1_PREF1_RESULT = "results/fgc-1-sgb1-ctl1-sol1-pref1.json"
_SOL1_PREF1_OWNER = "docs/fgc-sgb1-ctl1-sol1-pref1.md"
_SOL1_PREF1_REPRODUCER = "scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py"
_SOL1_PREF1_SOURCE = "src/recursive_horizons/fgc/sgb1_ctl1_sol1_pref1_binder.py"
_SOL1_PREF1_TESTS = (
    "tests/test_fgc_sgb1_ctl1_sol1_pref1_binder.py",
    "tests/test_check_repo_sgb1_ctl1_sol1_pref1.py",
)
_SOL1_PREF1_CLASSIFICATION = (
    "independently_bound_sol1_nominal_interval_inconclusive_no_health"
)
_TRANSIENT_TEST_ID = "FGC-1-PRO19-SID3-REAL-STORE-PREFLIGHT"
_TRANSIENT_TEST_PATH = (
    "archive/historical-tests/test_fgc_pro19_sid3_real_store_preflight.py"
)
_TRANSIENT_TEST_SHA256 = (
    "09d3515947530384063ad83c11faaea25ef7d43d527f23b15b50597067e63b3e"
)
_TDG8_OWNER = "docs/fgc-tdg8-frz1.md"
_CATALOG_STATUSES = {
    "compact_active",
    "current_frontier",
    "historical_authority",
    "historical_runner_disabled",
    "prospective_authority",
    "superseded_but_hash_bound",
    "transient_test_archived",
}
_CATALOG_KEYS = {
    "artifact_id",
    "classification",
    "compact_verify_target_or_none",
    "conclusion_or_nonclaim",
    "config_paths",
    "explicit_live_target_or_none",
    "family",
    "may_execute_state_change",
    "owner_documents",
    "predecessor_ids",
    "raw_dependencies",
    "reproducer_paths",
    "result_paths",
    "source_paths",
    "status",
    "successor_ids",
    "test_paths",
}
_CLASSIFICATION_KEYS = {"classification", "heading_artifact_id", "path", "sha256"}
_PATH_LIST_KEYS = {
    "config_paths",
    "owner_documents",
    "reproducer_paths",
    "source_paths",
    "test_paths",
}
_EXISTENCE_PATH_KEYS = _PATH_LIST_KEYS
_HEADING_ID = re.compile(
    r"(?m)^#\s+(?P<id>(?:FGC-[12]-[A-Z0-9]+(?:-[A-Z0-9]+)*"
    r"|GMF-[A-Z0-9]+(?:-[A-Z0-9]+)*"
    r"|EC-[A-Z0-9]+(?:-[A-Z0-9]+)*))"
)
_FGC_HEADING_ID = re.compile(r"\AFGC-[12]-[A-Z0-9]+(?:-[A-Z0-9]+)*\Z")
_MAKE_LIST_BODY = re.compile(
    r"(?P<name>[A-Z0-9_]+) := \\\n(?P<body>(?:\t[^\n]+(?: \\)?\n)+)"
)
_MAKEFILE_TARGET = re.compile(r"(?m)^([A-Za-z0-9_.-]+):")
_STATE_PREFIXES = ("run-", "resume-", "recover-")
_CLOSED_COMPACT_TARGETS = {
    "run-fgc-pro13-calibration": "verify-fgc-pro13-result",
    "run-fgc-pro14-calibration": "verify-fgc-pro14-result",
    "recover-fgc-pro19-sid1": "fgc-pro19-sid2-pref1",
    "resume-fgc-pro19-sid3-event1": "fgc-pro19-pref28",
    "run-fgc-pro19-event1": "fgc-pro19-pref28",
    "resume-fgc-pro19-event1": "fgc-pro19-pref28",
    "run-fgc-tdg8-successor-event1": "verify-fgc-tdg8-rcv3-pref2",
    "recover-fgc-tdg8-rcv1": "verify-fgc-tdg8-rcv3-pref2",
    "resume-fgc-tdg8-rcv1": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv2": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv3": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv3-rec1": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg9-ar1": "verify-fgc-tdg9-ar1-pref1",
    "run-fgc-tdg9-loc1": "verify-fgc-tdg9-loc2-pref2",
    "run-fgc-tdg9-loc2": "verify-fgc-tdg9-loc2-pref2",
    "run-fgc-tdg9-ti1": "verify-fgc-tdg9-ti2-pref1",
    "run-fgc-tdg9-ti2": "verify-fgc-tdg9-ti2-pref1",
    "run-fgc-tdg9-ac1": "verify-fgc-tdg9-ac1-pref1",
    "run-fgc-tdg9-ur1": "verify-fgc-tdg9-ur1-pref1",
    "run-fgc-tdg10-qa1": "verify-fgc-tdg10-qa1-pref1",
    "run-fgc-tdg10-qa2": "verify-fgc-tdg10-qa2-pref1",
    "run-fgc-tdg11-msel1": "fgc-tdg11-msel1-pref1",
    "run-fgc-hlt17-srcq1-rec1": "fgc-hlt17-srcq1-rec1-pref1",
    "run-fgc-pro20-ev1": "fgc-pro20-ev1-pref1",
}


def _catalog_path(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} is not a path")
    path = Path(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"{label} is unsafe")
    return value


def _catalog_top_id(path: Path) -> str | None:
    try:
        raw = path.read_bytes()
        if path.suffix == ".json":
            value = json.loads(raw, object_pairs_hook=_reject_duplicate_json_pairs)
        elif path.suffix == ".toml":
            value = tomllib.loads(raw.decode("utf-8"))
        else:
            return None
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        tomllib.TOMLDecodeError,
        ValueError,
    ):
        return None
    artifact_id = value.get("artifact_id") if isinstance(value, dict) else None
    return artifact_id if isinstance(artifact_id, str) and artifact_id else None


def _parse_make_list(source: str, name: str) -> tuple[str, ...]:
    for matched in _MAKE_LIST_BODY.finditer(source):
        if matched.group("name") != name:
            continue
        return tuple(
            line.strip().removesuffix(" \\")
            for line in matched.group("body").splitlines()
            if line.strip()
        )
    raise ValueError(f"Make list {name} is absent")


def _tracked_paths() -> tuple[str, ...]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=REPOSITORY,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    chunks = completed.stdout.split(b"\0")
    if chunks[-1] != b"":
        raise ValueError("Git path stream is not NUL terminated")
    return tuple(sorted(item.decode("utf-8") for item in chunks[:-1]))


def _catalog_coverage(tracked: tuple[str, ...]) -> dict[str, dict[str, set[str]]]:
    coverage: dict[str, dict[str, set[str]]] = {}
    for relative in tracked:
        path = REPOSITORY / relative
        if relative == _CATALOG_PATH:
            continue
        if relative.startswith("configs/fgc/") and relative.endswith((".toml", ".json")):
            key = "config_paths"
        elif relative.startswith("results/") and relative.endswith(".json"):
            key = "result_paths"
        else:
            continue
        artifact_id = _catalog_top_id(path)
        if artifact_id is None:
            continue
        record = coverage.setdefault(
            artifact_id, {"config_paths": set(), "result_paths": set()}
        )
        record[key].add(relative)
    return coverage


def _heading_owner_documents(tracked: tuple[str, ...]) -> dict[str, set[str]]:
    owners: dict[str, set[str]] = {}
    for relative in tracked:
        if not relative.startswith("docs/") or not relative.endswith(".md"):
            continue
        text = (REPOSITORY / relative).read_text(encoding="utf-8")
        matched = _HEADING_ID.search(text)
        if matched is None:
            continue
        heading = matched.group("id")
        owners.setdefault(heading, set()).add(relative)
    return owners


def check_artifact_catalog() -> list[str]:
    failures: list[str] = []
    try:
        raw = _read_unique_regular_bytes(_CATALOG_PATH, "artifact catalog")
        catalog = json.loads(
            raw.decode("ascii"), object_pairs_hook=_reject_duplicate_json_pairs
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        return [f"cannot inspect artifact catalog: {exc}"]
    if (
        not isinstance(catalog, dict)
        or set(catalog)
        != {"schema", "allowed_statuses", "frontier", "explicit_classifications", "artifacts"}
        or catalog.get("schema") != _CATALOG_SCHEMA
        or set(catalog.get("allowed_statuses", ())) != _CATALOG_STATUSES
        or list(catalog.get("allowed_statuses", ())) != sorted(_CATALOG_STATUSES)
    ):
        return ["artifact catalog top-level contract differs"]
    frontier = catalog.get("frontier")
    if (
        not isinstance(frontier, dict)
        or frontier.get("current_artifact_id") != _CATALOG_CURRENT
        or frontier.get("next_design_artifact_id") != _CATALOG_NEXT
        or frontier.get("state_advance_authorized") is not False
    ):
        failures.append("artifact catalog frontier differs")
    artifacts = catalog.get("artifacts")
    classifications = catalog.get("explicit_classifications")
    if not isinstance(artifacts, list) or not artifacts:
        return ["artifact catalog entries are absent"]
    if not isinstance(classifications, list):
        return ["artifact catalog explicit classifications are absent"]
    ids: list[str] = []
    by_id: dict[str, dict[str, object]] = {}
    status_counts: dict[str, int] = {name: 0 for name in _CATALOG_STATUSES}
    classified_paths: set[str] = set()
    for index, entry in enumerate(classifications):
        if not isinstance(entry, dict) or set(entry) != _CLASSIFICATION_KEYS:
            failures.append(f"artifact catalog classification {index} schema differs")
            continue
        try:
            relative = _catalog_path(entry.get("path"), label=f"classification.{index}")
        except ValueError as exc:
            failures.append(str(exc))
            continue
        if relative in classified_paths:
            failures.append(f"artifact catalog classification path repeats: {relative}")
        classified_paths.add(relative)
        probe_error: OSError | ValueError | None = None
        try:
            present = _probe_unique_regular_leaf(
                relative, f"catalog classification {relative}"
            )
        except (OSError, ValueError) as exc:
            present = False
            probe_error = exc
        if probe_error is not None:
            failures.append(
                f"artifact catalog classified path is unsafe: {relative}: {probe_error}"
            )
        elif not present:
            failures.append(f"artifact catalog classified path is absent: {relative}")
        else:
            try:
                digest = sha256(
                    _read_unique_regular_bytes(
                        relative, f"catalog classification {relative}"
                    )
                ).hexdigest()
            except (OSError, ValueError) as exc:
                failures.append(
                    f"artifact catalog classified path is unsafe: {relative}: {exc}"
                )
            else:
                if digest != entry.get("sha256"):
                    failures.append(
                        f"artifact catalog classified path hash differs: {relative}"
                    )
        heading = entry.get("heading_artifact_id")
        if heading is not None and not isinstance(heading, str):
            failures.append(f"artifact catalog classification heading differs: {relative}")
        if not isinstance(entry.get("classification"), str) or not entry["classification"]:
            failures.append(
                f"artifact catalog classification label differs: {relative}"
            )
    if [item.get("path") for item in classifications if isinstance(item, dict)] != sorted(
        classified_paths
    ):
        failures.append("artifact catalog classifications are not unique and sorted")
    for index, entry in enumerate(artifacts):
        if not isinstance(entry, dict) or set(entry) != _CATALOG_KEYS:
            failures.append(f"artifact catalog entry {index} schema differs")
            continue
        artifact_id = entry.get("artifact_id")
        if not isinstance(artifact_id, str) or not artifact_id:
            failures.append(f"artifact catalog entry {index} ID differs")
            continue
        if artifact_id.endswith(":") or artifact_id in {"Claim", "Closed", "FGC", "FGC-1:"}:
            failures.append(f"artifact catalog ID is not a strict artifact: {artifact_id}")
        ids.append(artifact_id)
        by_id[artifact_id] = entry
        status = entry.get("status")
        if status not in _CATALOG_STATUSES:
            failures.append(f"artifact catalog status differs: {artifact_id}")
        else:
            status_counts[status] += 1
        if not isinstance(entry.get("family"), str) or not entry["family"]:
            failures.append(f"artifact catalog family differs: {artifact_id}")
        execution_flag = entry.get("may_execute_state_change")
        if not isinstance(execution_flag, bool):
            failures.append(f"artifact catalog execution flag differs: {artifact_id}")
        elif execution_flag and status != "prospective_authority":
            failures.append(
                f"artifact catalog nonprospective execution flag differs: {artifact_id}"
            )
        if not isinstance(entry.get("conclusion_or_nonclaim"), str) or not entry[
            "conclusion_or_nonclaim"
        ]:
            failures.append(f"artifact catalog conclusion differs: {artifact_id}")
        classification = entry.get("classification")
        if classification is not None and not isinstance(classification, str):
            failures.append(f"artifact catalog classification differs: {artifact_id}")
        for optional_key in (
            "compact_verify_target_or_none",
            "explicit_live_target_or_none",
        ):
            value = entry.get(optional_key)
            if value is not None and not isinstance(value, str):
                failures.append(f"artifact catalog {optional_key} differs: {artifact_id}")
        for key in _PATH_LIST_KEYS:
            values = entry.get(key)
            if not isinstance(values, list) or values != sorted(set(values)):
                failures.append(f"artifact catalog {key} order differs: {artifact_id}")
                continue
            for value in values:
                try:
                    relative = _catalog_path(value, label=f"{artifact_id}.{key}")
                except ValueError as exc:
                    failures.append(str(exc))
                    continue
                if key in _EXISTENCE_PATH_KEYS:
                    try:
                        present = _probe_unique_regular_leaf(
                            relative, f"{artifact_id}.{key}"
                        )
                    except (OSError, ValueError) as exc:
                        failures.append(
                            f"artifact catalog path is unsafe: {relative}: {exc}"
                        )
                    else:
                        if not present:
                            failures.append(
                                f"artifact catalog path is absent: {relative}"
                            )
        raw_dependencies = entry.get("raw_dependencies")
        if not isinstance(raw_dependencies, list) or raw_dependencies != sorted(
            set(raw_dependencies)
        ):
            failures.append(f"artifact catalog raw_dependencies order differs: {artifact_id}")
        else:
            for value in raw_dependencies:
                try:
                    _catalog_path(value, label=f"{artifact_id}.raw")
                except ValueError as exc:
                    failures.append(str(exc))
        result_paths = entry.get("result_paths")
        if not isinstance(result_paths, list):
            failures.append(f"artifact catalog results differ: {artifact_id}")
        else:
            observed = []
            for result in result_paths:
                if not isinstance(result, dict) or set(result) != {"path", "sha256"}:
                    failures.append(f"artifact catalog result schema differs: {artifact_id}")
                    continue
                try:
                    relative = _catalog_path(
                        result.get("path"), label=f"{artifact_id}.result"
                    )
                except ValueError as exc:
                    failures.append(str(exc))
                    continue
                observed.append(relative)
                try:
                    result_bytes = _read_unique_regular_bytes(
                        relative, f"{artifact_id}.result"
                    )
                except FileNotFoundError:
                    failures.append(f"artifact catalog result is absent: {relative}")
                except (OSError, ValueError) as exc:
                    failures.append(
                        f"artifact catalog result is unsafe: {relative}: {exc}"
                    )
                else:
                    if sha256(result_bytes).hexdigest() != result.get("sha256"):
                        failures.append(
                            f"artifact catalog result hash differs: {relative}"
                        )
            if observed != sorted(set(observed)):
                failures.append(f"artifact catalog result order differs: {artifact_id}")
        for edge_key in ("predecessor_ids", "successor_ids"):
            values = entry.get(edge_key)
            if not isinstance(values, list) or values != sorted(set(values)):
                failures.append(f"artifact catalog {edge_key} order differs: {artifact_id}")
    if ids != sorted(set(ids)):
        failures.append("artifact catalog IDs are not unique and sorted")
    current = by_id.get(_CATALOG_CURRENT)
    imp1 = by_id.get(_IMP1_ARTIFACT)
    successor = by_id.get(_CATALOG_NEXT)
    qa2 = by_id.get(_QA2_PREF1)
    frozen = by_id.get(_MSEL1_FRZ1)
    pref1 = by_id.get(_PREF1_ARTIFACT)
    frz1 = by_id.get(_PRO20_FRZ1)
    rec1 = by_id.get(_REC1_PREF1)
    def1_frz1 = by_id.get(_DEF1_FRZ1)
    def1_pref1 = by_id.get(_DEF1_PREF1)
    sol1_frz1 = by_id.get(_SOL1_FRZ1)
    sol1_pref1 = by_id.get(_SOL1_PREF1)
    if (
        current is None
        or current.get("status") != "current_frontier"
        or current.get("may_execute_state_change") is not False
        or current.get("explicit_live_target_or_none") is not None
        or current.get("compact_verify_target_or_none") != "fgc-pro20-ev1-pref1"
        or current.get("classification")
        != "independently_bound_resource_exhausted_no_common_event"
        or _PRO20_FRZ1 not in current.get("predecessor_ids", ())
        or _CATALOG_NEXT not in current.get("successor_ids", ())
        or successor is None
        or successor.get("status") != "prospective_authority"
        or successor.get("predecessor_ids") != [_CATALOG_CURRENT]
        or successor.get("may_execute_state_change") is not False
        or successor.get("explicit_live_target_or_none") is not None
        or successor.get("compact_verify_target_or_none") is not None
        or successor.get("config_paths") != []
        or successor.get("result_paths") != []
        or successor.get("classification") != _RSRC1_CLASSIFICATION
        or set(successor.get("source_paths", ())) != _RSRC1_SOURCES
        or set(successor.get("reproducer_paths", ())) != _RSRC1_REPRODUCERS
        or set(successor.get("test_paths", ())) != _RSRC1_TESTS
        or successor.get("raw_dependencies") != []
        or "docs/fgc-pro20-ev1-rsrc1.md" not in successor.get("owner_documents", ())
        or "no-store per-member seed and attempt" not in str(
            successor.get("conclusion_or_nonclaim", "")
        ).lower()
        or "RSRC1-FRZ1" not in str(successor.get("conclusion_or_nonclaim", ""))
        or frz1 is None
        or frz1.get("status") != "historical_authority"
        or frz1.get("may_execute_state_change") is not False
        or frz1.get("explicit_live_target_or_none") is not None
        or frz1.get("compact_verify_target_or_none") != "fgc-pro20-ev1-frz1"
        or frz1.get("config_paths") != [_PRO20_CONFIG]
        or not _PRO20_SOURCES.issubset(set(frz1.get("source_paths", ())))
        or not _PRO20_TESTS.issubset(set(frz1.get("test_paths", ())))
        or set(frz1.get("reproducer_paths", ())) != _PRO20_REPRODUCERS
        or _PRO20_OWNER not in frz1.get("owner_documents", ())
        or rec1 is None
        or rec1.get("status") != "compact_active"
        or rec1.get("classification")
        != "independently_bound_all_six_source_c1r1_qualified_no_state_advance"
        or rec1.get("compact_verify_target_or_none") != "fgc-hlt17-srcq1-rec1-pref1"
        or _REC1_FRZ1 not in rec1.get("predecessor_ids", ())
        or def1_frz1 is None
        or def1_frz1.get("status") != "compact_active"
        or def1_frz1.get("classification") != _DEF1_FRZ1_CLASSIFICATION
        or def1_frz1.get("may_execute_state_change") is not False
        or def1_frz1.get("explicit_live_target_or_none") is not None
        or def1_frz1.get("compact_verify_target_or_none") != "fgc-def1-stab1-frz1"
        or def1_frz1.get("raw_dependencies") != []
        or "FGC-1-DEF1-STAB1" not in def1_frz1.get("predecessor_ids", ())
        or _DEF1_PREF1 not in def1_frz1.get("successor_ids", ())
        or def1_pref1 is None
        or def1_pref1.get("status") != "compact_active"
        or def1_pref1.get("classification") != _DEF1_PREF1_CLASSIFICATION
        or def1_pref1.get("may_execute_state_change") is not False
        or def1_pref1.get("explicit_live_target_or_none") is not None
        or def1_pref1.get("compact_verify_target_or_none") != "fgc-def1-stab1-pref1"
        or def1_pref1.get("raw_dependencies") != []
        or _DEF1_FRZ1 not in def1_pref1.get("predecessor_ids", ())
        or sol1_frz1 is None
        or sol1_frz1.get("status") != "compact_active"
        or sol1_frz1.get("classification") != _SOL1_FRZ1_CLASSIFICATION
        or sol1_frz1.get("may_execute_state_change") is not False
        or sol1_frz1.get("explicit_live_target_or_none") is not None
        or sol1_frz1.get("compact_verify_target_or_none")
        != "fgc-sgb1-ctl1-sol1-frz1"
        or sol1_frz1.get("raw_dependencies") != []
        or "FGC-1-SGB1-CTL1-SOL1" not in sol1_frz1.get("predecessor_ids", ())
        or _SOL1_PREF1 not in sol1_frz1.get("successor_ids", ())
        or sol1_pref1 is None
        or sol1_pref1.get("status") != "compact_active"
        or sol1_pref1.get("classification") != _SOL1_PREF1_CLASSIFICATION
        or sol1_pref1.get("may_execute_state_change") is not False
        or sol1_pref1.get("explicit_live_target_or_none") is not None
        or sol1_pref1.get("compact_verify_target_or_none")
        != "fgc-sgb1-ctl1-sol1-pref1"
        or sol1_pref1.get("raw_dependencies") != []
        or _SOL1_FRZ1 not in sol1_pref1.get("predecessor_ids", ())
        or status_counts["current_frontier"] != 1
        or status_counts["prospective_authority"] != 1 + len(_PARALLEL_INSTRUMENTS)
        or [
            entry.get("artifact_id")
            for entry in by_id.values()
            if entry.get("may_execute_state_change") is True
        ]
        != []
        or [
            entry.get("artifact_id")
            for entry in by_id.values()
            if entry.get("explicit_live_target_or_none") is not None
        ]
        != []
    ):
        failures.append("artifact catalog current-to-next edge differs")
    for artifact_id, (classification, sources, owners) in _PARALLEL_INSTRUMENTS.items():
        entry = by_id.get(artifact_id, {})
        if (
            entry.get("status") != "prospective_authority"
            or entry.get("classification") != classification
            or set(entry.get("source_paths", ())) != sources
            or set(entry.get("owner_documents", ())) != owners
            or entry.get("may_execute_state_change") is not False
            or entry.get("compact_verify_target_or_none") is not None
            or entry.get("explicit_live_target_or_none") is not None
            or any(entry.get(key) for key in ("config_paths", "result_paths", "raw_dependencies", "reproducer_paths"))
        ):
            failures.append(f"artifact catalog partial instrument scope differs: {artifact_id}")
    if (
        pref1 is None
        or pref1.get("status") != "compact_active"
        or pref1.get("may_execute_state_change") is not False
        or pref1.get("explicit_live_target_or_none") is not None
        or pref1.get("compact_verify_target_or_none") != "fgc-tdg11-msel1-pref1"
        or _MSEL1_FRZ1 not in pref1.get("predecessor_ids", ())
        or _IMP1_ARTIFACT not in pref1.get("successor_ids", ())
    ):
        failures.append("artifact catalog compact-active PREF1 routing differs")
    if (
        qa2 is None
        or qa2.get("status") == "current_frontier"
        or "exact-C remedy" not in str(qa2.get("conclusion_or_nonclaim", ""))
    ):
        failures.append("artifact catalog QA2 exact-C rejection is not distinct")
    if (
        frozen is None
        or frozen.get("status") != "historical_authority"
        or frozen.get("may_execute_state_change") is not False
        or frozen.get("explicit_live_target_or_none") is not None
        or frozen.get("compact_verify_target_or_none") != "fgc-tdg11-msel1-frz1"
        or _PREF1_ARTIFACT not in frozen.get("successor_ids", ())
    ):
        failures.append("artifact catalog historical MSEL1 freeze routing differs")
    tdg8 = by_id.get("FGC-1-TDG8-FRZ1")
    if tdg8 is None or _TDG8_OWNER not in tdg8.get("owner_documents", ()):
        failures.append("artifact catalog does not link the TDG8 owner document")
    transient = by_id.get(_TRANSIENT_TEST_ID)
    if (
        transient is None
        or transient.get("status") != "transient_test_archived"
        or _TRANSIENT_TEST_PATH not in transient.get("test_paths", ())
        or status_counts["transient_test_archived"] != 1
    ):
        failures.append("artifact catalog transient archive record differs")
    else:
        try:
            transient_bytes = _read_unique_regular_bytes(
                _TRANSIENT_TEST_PATH, "artifact catalog transient archive"
            )
        except (OSError, ValueError) as exc:
            failures.append(f"artifact catalog transient archive is unsafe: {exc}")
        else:
            if sha256(transient_bytes).hexdigest() != _TRANSIENT_TEST_SHA256:
                failures.append("artifact catalog transient archive hash differs")
    make_source = "\n".join(
        (REPOSITORY / relative).read_text(encoding="utf-8")
        for relative in (
            "Makefile",
            "mk/historical-certificates.mk",
            "mk/current-foundation.mk",
            "mk/closed-live-targets.mk",
        )
    )
    try:
        closed = _parse_make_list(make_source, "CLOSED_HISTORICAL_STATE_TARGETS")
    except ValueError as exc:
        failures.append(str(exc))
        closed = ()
    if set(closed) != set(_CLOSED_COMPACT_TARGETS):
        failures.append("closed-operation compact-verifier map differs")
    disabled = tuple(
        artifact_id
        for artifact_id, entry in by_id.items()
        if entry.get("status") == "historical_runner_disabled"
    )
    if disabled != tuple(sorted(closed)) or set(disabled) != set(closed):
        failures.append("artifact catalog disabled operations differ from closed targets")
    for target in closed:
        entry = by_id.get(target)
        if (
            entry is None
            or entry.get("status") != "historical_runner_disabled"
            or entry.get("explicit_live_target_or_none") is not None
            or entry.get("compact_verify_target_or_none")
            != _CLOSED_COMPACT_TARGETS.get(target)
            or entry.get("may_execute_state_change") is not False
        ):
            failures.append(f"artifact catalog disabled operation differs: {target}")
    make_targets = set(_MAKEFILE_TARGET.findall(make_source)) | set(closed)
    for artifact_id, entry in by_id.items():
        compact_target = entry.get("compact_verify_target_or_none")
        if compact_target is not None and compact_target not in make_targets:
            failures.append(
                f"artifact catalog compact target is absent: {artifact_id}: "
                f"{compact_target}"
            )
        live_target = entry.get("explicit_live_target_or_none")
        if live_target is None:
            if entry.get("may_execute_state_change") is True:
                failures.append(
                    f"artifact catalog executable authority lacks a live target: "
                    f"{artifact_id}"
                )
            continue
        if (
            live_target not in make_targets
            or not live_target.startswith(_STATE_PREFIXES)
            or live_target in closed
            or entry.get("status") != "prospective_authority"
            or entry.get("may_execute_state_change") is not True
        ):
            failures.append(
                f"artifact catalog live target authority differs: {artifact_id}: "
                f"{live_target}"
            )
    state_changing = sorted(
        name for name in make_targets if name.startswith(_STATE_PREFIXES)
    )
    for target in state_changing:
        if target in closed:
            continue
        owners = [
            artifact_id
            for artifact_id, entry in by_id.items()
            if entry.get("status") == "prospective_authority"
            and entry.get("may_execute_state_change") is True
            and entry.get("explicit_live_target_or_none") == target
        ]
        if len(owners) != 1:
            failures.append(
                "state-changing Make target lacks one active prospective authority: "
                f"{target}: {owners}"
            )
    for artifact_id, entry in by_id.items():
        if (
            artifact_id != _CATALOG_CURRENT
            and entry.get("status") == "current_frontier"
        ):
            failures.append(f"artifact catalog marks non-frontier current: {artifact_id}")
        if entry.get("status") == "current_frontier" and any(
            Path(path).name.startswith("run_")
            for path in entry.get("reproducer_paths", ())
            if isinstance(path, str)
        ):
            failures.append("existing runner marked the current frontier executable")
        for predecessor in entry.get("predecessor_ids", ()):
            if predecessor in by_id and artifact_id not in by_id[predecessor].get(
                "successor_ids", ()
            ):
                failures.append(
                    f"artifact catalog predecessor/successor edge differs: "
                    f"{predecessor} -> {artifact_id}"
                )
    tracked = _tracked_paths()
    for entry, label, config, result, owner, reproducer, sources, tests in (
        (
            pref1, "PREF1", _PREF1_CONFIG, _PREF1_RESULT, _PREF1_OWNER,
            _PREF1_REPRODUCER, _PREF1_SOURCES, _PREF1_TESTS,
        ),
        (
            imp1, "IMP1", _IMP1_CONFIG, _IMP1_RESULT, _IMP1_OWNER,
            _IMP1_REPRODUCER, _IMP1_SOURCES, _IMP1_TESTS,
        ),
        (
            def1_frz1, "DEF1-STAB1-FRZ1", _DEF1_FRZ1_CONFIG, _DEF1_FRZ1_RESULT,
            _DEF1_FRZ1_OWNER, _DEF1_FRZ1_REPRODUCER, {_DEF1_FRZ1_SOURCE},
            _DEF1_FRZ1_TESTS,
        ),
        (
            def1_pref1, "DEF1-STAB1-PREF1", _DEF1_PREF1_CONFIG, _DEF1_PREF1_RESULT,
            _DEF1_PREF1_OWNER, _DEF1_PREF1_REPRODUCER, {_DEF1_PREF1_SOURCE},
            _DEF1_PREF1_TESTS,
        ),
        (
            sol1_frz1, "SGB1-CTL1-SOL1-FRZ1", _SOL1_FRZ1_CONFIG,
            _SOL1_FRZ1_RESULT, _SOL1_FRZ1_OWNER, _SOL1_FRZ1_REPRODUCER,
            {_SOL1_FRZ1_SOURCE}, _SOL1_FRZ1_TESTS,
        ),
        (
            sol1_pref1, "SGB1-CTL1-SOL1-PREF1", _SOL1_PREF1_CONFIG,
            _SOL1_PREF1_RESULT, _SOL1_PREF1_OWNER, _SOL1_PREF1_REPRODUCER,
            {_SOL1_PREF1_SOURCE}, _SOL1_PREF1_TESTS,
        ),
    ):
        if (
            entry is None
            or config not in entry.get("config_paths", ())
            or result not in {
                item.get("path")
                for item in entry.get("result_paths", ())
                if isinstance(item, dict)
            }
            or owner not in entry.get("owner_documents", ())
            or reproducer not in entry.get("reproducer_paths", ())
            or set(entry.get("source_paths", ())) != set(sources)
            or any(path not in entry.get("test_paths", ()) for path in tests)
        ):
            failures.append(
                f"artifact catalog {label} source/config/result/owner/test association differs"
            )
    coverage = _catalog_coverage(tracked)
    for artifact_id, expected in coverage.items():
        entry = by_id.get(artifact_id)
        if entry is None:
            failures.append(f"artifact catalog coverage is absent: {artifact_id}")
            continue
        if set(entry.get("config_paths", ())) != expected["config_paths"]:
            failures.append(f"artifact catalog config coverage differs: {artifact_id}")
        observed_results = {
            item.get("path")
            for item in entry.get("result_paths", ())
            if isinstance(item, dict)
        }
        if observed_results != expected["result_paths"]:
            failures.append(f"artifact catalog result coverage differs: {artifact_id}")
    covered_configs: set[str] = set()
    covered_results: set[str] = set()
    for entry in by_id.values():
        covered_configs.update(
            path for path in entry.get("config_paths", ()) if isinstance(path, str)
        )
        covered_results.update(
            item.get("path")
            for item in entry.get("result_paths", ())
            if isinstance(item, dict)
        )
    for relative in tracked:
        if relative == _CATALOG_PATH:
            continue
        if relative.startswith("configs/fgc/") and relative.endswith((".toml", ".json")):
            if relative in covered_configs:
                continue
            if relative not in classified_paths:
                failures.append(f"artifact catalog config is unclassified: {relative}")
        elif relative.startswith("results/") and relative.endswith(".json"):
            if relative in covered_results:
                continue
            if relative not in classified_paths:
                failures.append(f"artifact catalog result is unclassified: {relative}")
    heading_owners = _heading_owner_documents(tracked)
    linked_owners: set[str] = set()
    for entry in by_id.values():
        linked_owners.update(
            path for path in entry.get("owner_documents", ()) if isinstance(path, str)
        )
    for heading, paths in heading_owners.items():
        if _FGC_HEADING_ID.fullmatch(heading):
            entry = by_id.get(heading)
            if entry is None:
                failures.append(f"artifact catalog heading artifact is absent: {heading}")
                continue
            missing = paths - set(entry.get("owner_documents", ()))
            if missing:
                failures.append(
                    f"artifact catalog heading owner is unlinked: {sorted(missing)}"
                )
            continue
        for path in paths:
            if path not in classified_paths and path not in linked_owners:
                failures.append(
                    f"artifact catalog heading document is unclassified: {path}"
                )
    completed = subprocess.run(
        [sys.executable, "scripts/build_artifact_catalog.py", "--check"],
        cwd=REPOSITORY,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        failures.append(
            "artifact catalog/result index regeneration differs"
            + (f": {detail}" if detail else "")
        )
    return failures


def _check_def1_stab1_frz1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.def1_stab1_frz1_certificate import verify_compact

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"DEF1-STAB1-FRZ1 compact freeze differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_def1_stab1_pref1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.def1_stab1_pref1_binder import verify_compact

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"DEF1-STAB1-PREF1 compact binder differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_sgb1_ctl1_sol1_frz1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.sgb1_ctl1_sol1_frz1_certificate import verify_compact

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"SGB1-CTL1-SOL1-FRZ1 compact freeze differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


def _check_sgb1_ctl1_sol1_pref1_result() -> list[str]:
    previous = list(sys.path)
    try:
        sys.path.insert(0, str(REPOSITORY / "src"))
        from recursive_horizons.fgc.sgb1_ctl1_sol1_pref1_binder import verify_compact

        verify_compact(REPOSITORY)
    except (OSError, ValueError, RuntimeError) as error:
        return [f"SGB1-CTL1-SOL1-PREF1 compact binder differs: {error}"]
    finally:
        sys.path[:] = previous
    return []


__all__ = tuple(name for name in globals() if not name.startswith("__"))
