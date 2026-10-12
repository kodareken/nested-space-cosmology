#!/usr/bin/env python3
"""Read-only authentication of the sealed spherical feedback episode.

``verify_saved`` on the current producer compares frozen hashes with the
working tree. The working-tree Galerkin file is not the producer pin, so
that check stops. This consumer recovers the producer Galerkin bytes from
an explicit commit, authenticates the episode's frozen dependency context
and saved arrays, and parses the historical import graph before deciding
whether that source may be loaded. It does not import the historical
module and does not generate episode files.
"""
from __future__ import annotations

import ast
import hashlib
import io
import json
import os
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

from recursive_horizons.provenance import resolve_pinned_source_bytes

LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
EPISODE_JSON = LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"
EPISODE_NPZ = LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"

SCHEMA = "NSC-SPHERICAL-FEEDBACK-EPISODE-v1"
CONSUMER_SCHEMA = "NSC-SPHERICAL-FEEDBACK-EPISODE-HISTORICAL-READ-v1"
HISTORICAL_COMMIT = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"
HISTORICAL_GALERKIN_SHA256 = (
    "85d1f3dbbd82a38fe14d2405ad9782af6de68095052ff522ac6602e5d50a14ab"
)
CURRENT_GALERKIN_SHA256 = (
    "aa8a64f593e8f0d931eeee407d944b72568fdf4f2c53436048142c3473ec4911"
)

# Frozen names stored by the episode producer, in record order.
DEPENDENCIES = (
    ("v5_json", "lab/results/development/nsc-spherical-coupling-refinement-v5.json", "scientific_record"),
    ("v5_npz", "lab/results/development/nsc-spherical-coupling-refinement-v5.npz", "scientific_record"),
    ("regional_json", "lab/results/development/nsc-regional-energy-exchange-v1.json", "scientific_record"),
    ("v1_json", "lab/results/development/nsc-spherical-coupling-control-v1.json", "scientific_record"),
    ("v2_json", "lab/results/development/nsc-spherical-coupling-control-v2.json", "scientific_record"),
    ("v3_json", "lab/results/development/nsc-spherical-coupling-refinement-v3.json", "scientific_record"),
    ("galerkin", "lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py", "scientific_source"),
    ("regional_module", "lab/src/recursive_horizons/nsc_regional_energy_exchange.py", "scientific_source"),
    ("feedback_action", "lab/src/recursive_horizons/nsc_spherical_feedback_action.py", "scientific_source"),
    ("coupling", "lab/src/recursive_horizons/nsc_spherical_coupling.py", "scientific_source"),
    ("conformal_source", "lab/src/recursive_horizons/nsc_conformal_adm_source.py", "scientific_source"),
    ("manuscript_md", "paper/nested-space-cosmology.md", "manuscript"),
    ("manuscript_pdf", "paper/nested-space-cosmology.pdf", "manuscript"),
    ("companion_pdf", "paper/local-incoming-gate-draft.pdf", "manuscript"),
    ("driver", "lab/scripts/derive_nsc_spherical_feedback_episode.py", "producer"),
)

PIN_TO_MODULE = {
    "galerkin": "nsc_spherical_galerkin_coupling",
    "regional_module": "nsc_regional_energy_exchange",
    "feedback_action": "nsc_spherical_feedback_action",
    "coupling": "nsc_spherical_coupling",
    "conformal_source": "nsc_conformal_adm_source",
}
MODULE_TO_PIN = {module: pin for pin, module in PIN_TO_MODULE.items()}

SETUP_COPIES = (
    ("nf512_initial_r", "nf512_geometry_r"),
    ("nf512_initial_Q", "nf512_geometry_Q"),
    ("nf512_initial_chi", "nf512_geometry_chi"),
    ("nf512_initial_phi0", "nf512_columns_phi0"),
    ("nf512_initial_phi1", "nf512_columns_phi1"),
    ("nf256_initial_r", "geometry_r"),
    ("nf256_initial_phi0", "columns_phi0"),
)
REPLAY_STATE = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
REPLAY_RATES = (
    "Q_dot", "r_dot", "chi_dot", "p_Q_dot", "p_r_dot", "p_chi_dot", "phi0_dot", "phi1_dot",
)

PRESERVED = tuple(path for _name, path, _role in DEPENDENCIES) + (
    "lab/results/development/nsc-spherical-cauchy-data-v1.json",
    "lab/results/development/nsc-spherical-feedback-episode-v1.json",
    "lab/results/development/nsc-spherical-feedback-episode-v1.npz",
)


def _sha256_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def _sha256_path(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _snapshot():
    return {relative: _sha256_path(ROOT / relative) for relative in dict.fromkeys(PRESERVED)}


def _assignment_names(target):
    if isinstance(target, ast.Name):
        return {target.id}
    if isinstance(target, (ast.Tuple, ast.List)):
        found = set()
        for item in target.elts:
            found.update(_assignment_names(item))
        return found
    return set()


def _bound_names(tree):
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                names.update(_assignment_names(target))
        elif isinstance(node, ast.AnnAssign):
            names.update(_assignment_names(node.target))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.asname or alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if alias.name != "*":
                    names.add(alias.asname or alias.name)
    return names


def _module_body_imports(tree):
    rows = []
    calls = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            rows.append({
                "kind": "import",
                "level": 0,
                "module": None,
                "names": [alias.name for alias in node.names],
            })
        elif isinstance(node, ast.ImportFrom):
            rows.append({
                "kind": "from",
                "level": node.level,
                "module": node.module,
                "names": [alias.name for alias in node.names],
            })
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            calls.append(ast.unparse(node.value.func))
    return rows, calls


def inspect_historical_dependencies(pinned_sources):
    """Parse pinned source bytes. This does not import or execute them.

    The returned gate is the decision that has to exist before any loader
    of the historical Galerkin module. This consumer stops at the gate.
    """
    if "galerkin" not in pinned_sources or "driver" not in pinned_sources:
        raise ValueError("historical inspection requires pinned Galerkin and producer bytes")
    graphs = {}
    for pin, module in PIN_TO_MODULE.items():
        tree = ast.parse(pinned_sources[pin])
        rows, calls = _module_body_imports(tree)
        graphs[module] = {"rows": rows, "calls": calls, "names": _bound_names(tree)}
    galerkin = graphs["nsc_spherical_galerkin_coupling"]
    direct = [row for row in galerkin["rows"] if row["kind"] == "from" and row["level"] == 1]
    missing = []
    for row in direct:
        target = graphs.get(row["module"])
        if target is None:
            missing.append(row["module"])
            continue
        absent = [name for name in row["names"] if name not in target["names"]]
        if absent:
            missing.append(row["module"] + ":" + ",".join(absent))
    seen = set()
    queue = ["nsc_spherical_galerkin_coupling"]
    unpinned = []
    authenticated = []
    while queue:
        module = queue.pop(0)
        if module in seen:
            continue
        seen.add(module)
        graph = graphs.get(module)
        if graph is None:
            unpinned.append(module)
            continue
        authenticated.append(module)
        for row in graph["rows"]:
            if row["kind"] == "from" and row["level"] >= 1 and row["module"]:
                queue.append(row["module"])
    driver_tree = ast.parse(pinned_sources["driver"])
    driver_names = _bound_names(driver_tree)
    import_allowed = not missing and not unpinned and not galerkin["calls"]
    return {
        "inspected_before_import": True,
        "historical_direct_modules": sorted(row["module"] for row in direct),
        "missing_names": missing,
        "authenticated_import_modules": sorted(set(authenticated)),
        "unpinned_transitive_modules": sorted(set(unpinned)),
        "module_level_calls": list(galerkin["calls"]),
        "import_allowed": import_allowed,
        "producer_defines_run_episode": "run_episode" in driver_names,
        "producer_defines_reassess_draft": "reassess_draft" in driver_names,
    }


def _finite(array):
    values = np.asarray(array)
    if values.size == 0 or not np.isfinite(values).all():
        return False
    return True


def _authenticate_saved_arrays(record, episode_npz, v5_npz):
    if _sha256_bytes(episode_npz) != record.get("npz_sha256"):
        raise AssertionError("episode payload hash does not match npz_sha256")
    with np.load(io.BytesIO(episode_npz), allow_pickle=False) as episode, np.load(
        io.BytesIO(v5_npz), allow_pickle=False
    ) as saved:
        for episode_key, saved_key in SETUP_COPIES:
            if episode_key not in episode.files or saved_key not in saved.files:
                raise AssertionError(f"saved setup array missing: {episode_key}")
            if not np.array_equal(episode[episode_key], saved[saved_key]):
                raise AssertionError(f"setup does not copy {saved_key}")
        for name in REPLAY_STATE:
            key = "replay_" + name
            if key not in episode.files or not _finite(episode[key]):
                raise AssertionError(f"saved replay array missing: {key}")
        for name in REPLAY_RATES:
            key = "replay_rate_" + name
            if key not in episode.files or not _finite(episode[key]):
                raise AssertionError(f"saved replay rate missing: {key}")
        if "replay_time" not in episode.files:
            raise AssertionError("saved replay time missing")
        replay_time = float(episode["replay_time"])
        declared = record.get("replay") or {}
        if declared.get("stored") is not True or abs(replay_time - float(declared.get("time"))) > 1e-12:
            raise AssertionError("saved replay time does not match the record")
        primary = "nf512_dt_0_0005_"
        run = (record.get("runs") or {}).get("nf512_dt_0_0005")
        if run is None or primary + "time" not in episode.files or primary + "proper_max" not in episode.files:
            raise AssertionError("primary saved series missing")
        attained = float(episode[primary + "time"][-1])
        if abs(attained - float(run["attained_T"])) > 1e-12:
            raise AssertionError("attained T does not match the series")
        proper = episode[primary + "proper_max"]
        if abs(float(proper[0]) - float(run["initial"]["proper_max"])) > 1e-12:
            raise AssertionError("initial proper max does not match the series")
        if abs(float(proper[-1]) - float(run["final"]["proper_max"])) > 1e-12:
            raise AssertionError("final proper max does not match the series")
        if "verdict" not in episode.files or str(episode["verdict"]) != record.get("verdict"):
            raise AssertionError("saved verdict array does not match the record")
    return True


def authenticate():
    """Authenticate the sealed episode in memory and return the report.

    Historical Galerkin bytes come only from ``resolve_pinned_source_bytes``.
    A working-tree mismatch stays visible as current-code incompatibility.
    """
    repository = ROOT
    before = _snapshot()
    if not EPISODE_JSON.is_file() or not EPISODE_NPZ.is_file():
        raise FileNotFoundError("episode evidence is missing")
    episode_json = EPISODE_JSON.read_bytes()
    episode_npz = EPISODE_NPZ.read_bytes()
    record = json.loads(episode_json)
    if record.get("schema") != SCHEMA:
        raise AssertionError("schema mismatch")
    if record.get("frozen_bytes_unchanged") is not True:
        raise AssertionError("episode does not declare frozen bytes unchanged")
    before_hashes = record.get("hashes_before")
    after_hashes = record.get("hashes_after")
    expected = {name for name, _path, _role in DEPENDENCIES}
    if set(before_hashes or {}) != expected or set(after_hashes or {}) != expected:
        raise AssertionError("episode dependency context is not the frozen set")
    if before_hashes["galerkin"] != HISTORICAL_GALERKIN_SHA256 or after_hashes["galerkin"] != HISTORICAL_GALERKIN_SHA256:
        raise AssertionError("episode Galerkin pin is not the named historical bytes")
    if before_hashes["driver"] == after_hashes["driver"]:
        raise AssertionError("driver before/after hashes were collapsed")

    entries = []
    pinned_sources = {}
    for name, relative, role in DEPENDENCIES:
        declared_before = before_hashes[name]
        declared_after = after_hashes[name]
        if name != "driver" and declared_before != declared_after:
            raise AssertionError(f"{name} changed between hashes_before and hashes_after")
        working = repository / relative
        working_sha = _sha256_path(working) if working.is_file() else None
        if role in {"scientific_source", "producer"}:
            expected_sha = declared_after
            raw = resolve_pinned_source_bytes(
                repository, relative, expected_sha, commit=HISTORICAL_COMMIT,
            )
            if _sha256_bytes(raw) != expected_sha:
                raise AssertionError(f"pinned bytes do not match {name}")
            pinned_sources[name] = raw
            origin = "resolve_pinned_source_bytes"
        else:
            if working_sha != declared_after:
                raise AssertionError(f"working-tree record does not match {name}")
            origin = "working-tree"
        entries.append({
            "name": name,
            "path": relative,
            "role": role,
            "declared_before": declared_before,
            "declared_after": declared_after,
            "authenticated_sha256": declared_after if role == "producer" else declared_before,
            "origin": origin,
            "commit": HISTORICAL_COMMIT if origin == "resolve_pinned_source_bytes" else None,
            "working_tree_sha256": working_sha,
            "working_tree_matches_historical": working_sha == declared_after,
        })

    galerkin_raw = pinned_sources["galerkin"]
    if _sha256_bytes(galerkin_raw) != HISTORICAL_GALERKIN_SHA256:
        raise AssertionError("recovered Galerkin bytes are not the named scientific pin")
    inspection = inspect_historical_dependencies(pinned_sources)
    if inspection["inspected_before_import"] is not True:
        raise AssertionError("dependencies were not inspected")
    if inspection["missing_names"]:
        raise AssertionError(
            "historical import names are absent: " + ", ".join(inspection["missing_names"])
        )
    if inspection["module_level_calls"]:
        raise AssertionError("historical Galerkin module has module-level calls")
    if not inspection["unpinned_transitive_modules"]:
        raise AssertionError("import gate opened without an authenticated closure")
    # Unpinned transitive imports would bind to whatever is importable now.
    # The historical module and the episode generator stay unloaded.
    historical_source_imported = False
    producer_invoked = False
    current_galerkin = next(item for item in entries if item["name"] == "galerkin")
    current_bytes = (repository / current_galerkin["path"]).read_bytes()
    current_rows, _current_calls = _module_body_imports(ast.parse(current_bytes))
    current_direct = sorted(
        row["module"] for row in current_rows if row["kind"] == "from" and row["level"] == 1
    )
    current_sha = _sha256_bytes(current_bytes)
    if current_sha != current_galerkin["working_tree_sha256"]:
        raise AssertionError("Galerkin working-tree hash changed while authenticating")
    compatible = current_sha == HISTORICAL_GALERKIN_SHA256
    v5_relative = next(path for name, path, _role in DEPENDENCIES if name == "v5_npz")
    saved_arrays = _authenticate_saved_arrays(
        record, episode_npz, (repository / v5_relative).read_bytes(),
    )
    file_payload = len(episode_json) + len(episode_npz)
    if int(record.get("payload_bytes")) != file_payload:
        raise AssertionError("payload_bytes does not match the sealed files")
    if record.get("payload_within_64MiB") is not True or file_payload > 64 * 1024 * 1024:
        raise AssertionError("payload limit was not kept")
    after = _snapshot()
    if after != before:
        changed = sorted(name for name in before if before[name] != after.get(name))
        raise AssertionError("historical consumer changed sealed bytes: " + ", ".join(changed))
    return {
        "schema": CONSUMER_SCHEMA,
        "historical_commit": HISTORICAL_COMMIT,
        "historical_galerkin_sha256": HISTORICAL_GALERKIN_SHA256,
        "historical_galerkin_bytes": len(galerkin_raw),
        "source_origin": "resolve_pinned_source_bytes",
        "current_galerkin_sha256": current_sha,
        "current_code_compatible": compatible,
        "current_direct_modules": current_direct,
        "current_code_replay_claimed": False,
        "current_code_limit": (
            "working-tree Galerkin bytes are a different pin from the episode producer; "
            "current verify_saved stops on that pin, and this consumer does not "
            "substitute a current-code step for the saved replay"
        ),
        "historical_source_imported": historical_source_imported,
        "producer_invoked": producer_invoked,
        "generation_invoked": False,
        "saved_arrays_authenticated": saved_arrays,
        "dependencies": entries,
        "dependency_inspection": inspection,
        "driver_before_sha256": before_hashes["driver"],
        "driver_after_sha256": after_hashes["driver"],
        "driver_before_matches_commit": before_hashes["driver"] == _sha256_bytes(pinned_sources["driver"]),
        "sealed_bytes_unchanged": True,
        "episode_json_sha256": _sha256_bytes(episode_json),
        "episode_npz_sha256": _sha256_bytes(episode_npz),
    }


def main(argv=None):
    del argv
    report = authenticate()
    public = {
        key: report[key]
        for key in (
            "schema",
            "historical_commit",
            "historical_galerkin_sha256",
            "historical_galerkin_bytes",
            "current_galerkin_sha256",
            "current_code_compatible",
            "current_direct_modules",
            "current_code_replay_claimed",
            "historical_source_imported",
            "producer_invoked",
            "generation_invoked",
            "saved_arrays_authenticated",
            "dependency_inspection",
            "driver_before_sha256",
            "driver_after_sha256",
            "sealed_bytes_unchanged",
            "episode_json_sha256",
            "episode_npz_sha256",
        )
    }
    print(json.dumps(public, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
