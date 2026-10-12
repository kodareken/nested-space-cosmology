"""Read-only authentication and no-adoption controls for the HLT15 GEN1 store."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
from typing import Any, Mapping

from . import proto17_hlt15_runtime as hlt15_runtime
from .proto17_hlt15_runtime import Proto17HLT15Error, Proto17HLT15Runtime
from .proto17_pure_construction import Proto17ConstructionError, build_genesis, canonical, digest, validate_checkpoint


class Proto18Pref27Error(ValueError):
    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


ROOT_RELATIVE = "runs/fgc-2-sf1/proto17/calibration"
AUTH_RESULT = "results/fgc-1-pro18-auth1.json"
AUTH_CONFIG = "configs/fgc/fgc-1-pro18-auth1.toml"
EXPECTED_CLAIMS = {
    "generation_zero_materialized": True, "state_advanced": False,
    "calibration_completed": False, "candidate_execution_authorized": False,
    "physical_transition_claim_authorized": False,
}


def _stop(stop_id: str, detail: str) -> None:
    raise Proto18Pref27Error(stop_id, detail)


def validate_pref27_config(config: Mapping[str, Any]) -> None:
    """Reject scope/claim promotion before touching the real GEN1 store."""
    if not isinstance(config, Mapping) or (config.get("schema_version"), config.get("artifact_id"), config.get("project_version"), config.get("target_protocol")) != (1, "FGC-1-PRO18-PREF27", "0.11.0", "FGC-2-SF1-PROTO17"):
        _stop("PREF27_CONTRACT_DRIFT", "PREF27 identity differs")
    scope = config.get("scope"); predecessor = config.get("predecessor"); claims = config.get("claims")
    if (not isinstance(scope, Mapping) or any(scope.get(key) is not False for key in ("real_store_write_access", "holdout_root_observed", "state_advance", "candidate_execution"))
            or not isinstance(predecessor, Mapping) or predecessor != {"authority_artifact": "FGC-1-PRO18-AUTH1", "authority_result": AUTH_RESULT, "calibration_root": ROOT_RELATIVE, "expected_leaf_count": 8, "expected_state_count": 6, "absent_or_partial_hlt15_store_is_rejected": True, "raw_predecessors_not_reopened": True}
            or not isinstance(claims, Mapping) or claims != {"HLT15_generation_zero_store_authenticated": True, "production_namespace_reuse_refusal_verified": True, "foreign_store_adoption_refusal_verified": True, "all_four_HLT14_production_duties_completed": True, "state_advanced": False, "fresh_GR0_dynamic_calibration_authorized": False, "candidate_execution_authorized": False, "physical_transition_claim_authorized": False}):
        _stop("PREF27_CONTRACT_DRIFT", "PREF27 scope, predecessor, or claims differ")


def _json(raw: bytes, *, label: str, indented: bool = False) -> dict[str, Any]:
    def duplicate(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in pairs:
            if key in answer: raise ValueError(key)
            answer[key] = value
        return answer
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=duplicate)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("PREF27_STORE_JSON_INVALID", f"{label} is malformed")
    expected = ((json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()
                if indented else canonical(value))
    if not isinstance(value, dict) or expected != raw:
        _stop("PREF27_STORE_JSON_INVALID", f"{label} is noncanonical")
    return value


def _safe_file(root: Path, relative: str) -> bytes:
    target = root / relative
    if Path(relative).is_absolute() or ".." in Path(relative).parts:
        _stop("PREF27_STORE_PATH_UNSAFE", "unsafe relative path")
    current = root
    for part in Path(relative).parts[:-1]:
        current /= part
        try: item = current.lstat()
        except OSError as error: _stop("PREF27_STORE_PATH_UNSAFE", f"missing ancestor: {relative}")
        if stat.S_ISLNK(item.st_mode) or not stat.S_ISDIR(item.st_mode):
            _stop("PREF27_STORE_PATH_UNSAFE", f"unsafe ancestor: {relative}")
    try: item = target.lstat()
    except OSError: _stop("PREF27_STORE_PATH_UNSAFE", f"missing leaf: {relative}")
    if stat.S_ISLNK(item.st_mode) or not stat.S_ISREG(item.st_mode) or item.st_nlink != 1:
        _stop("PREF27_STORE_PATH_UNSAFE", f"unsafe leaf: {relative}")
    return target.read_bytes()


def _nonfollowing_entries(root: Path, *, label: str) -> list[tuple[Path, os.stat_result]]:
    """Return a lexical tree inventory without following any directory entry."""
    try:
        root_info = root.stat(follow_symlinks=False)
    except OSError as error:
        _stop("PREF27_STORE_PATH_UNSAFE", f"{label} cannot be statted")
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        _stop("PREF27_STORE_PATH_UNSAFE", f"{label} root is unsafe")
    pending = [root]
    entries: list[tuple[Path, os.stat_result]] = []
    while pending:
        directory = pending.pop()
        try:
            with os.scandir(directory) as iterator:
                children = sorted(iterator, key=lambda entry: entry.name, reverse=True)
        except OSError as error:
            _stop("PREF27_STORE_PATH_UNSAFE", f"{label} cannot be scanned")
        for child in children:
            path = Path(child.path)
            try:
                info = child.stat(follow_symlinks=False)
            except OSError as error:
                _stop("PREF27_STORE_PATH_UNSAFE", f"{label} entry cannot be statted")
            if stat.S_ISLNK(info.st_mode):
                _stop("PREF27_STORE_PATH_UNSAFE", f"symlink appears in {label}")
            if stat.S_ISDIR(info.st_mode):
                pending.append(path)
            elif stat.S_ISREG(info.st_mode):
                if info.st_nlink != 1:
                    _stop("PREF27_STORE_PATH_UNSAFE", f"hard-linked file appears in {label}")
            else:
                _stop("PREF27_STORE_PATH_UNSAFE", f"nonregular entry appears in {label}")
            entries.append((path, info))
    return entries


def _git(root: Path, *args: str) -> bytes:
    try:
        return subprocess.run(["git", *args], cwd=root, check=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        _stop("PREF27_AUTHORITY_GIT_INVALID", "Git authority lookup failed")


def _auth_record(root: Path, receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    binding = receipt.get("external_authority")
    if not isinstance(binding, Mapping) or set(binding) != {
        "artifact_id", "authorization_commit", "authorization_result_path",
        "authorization_result_sha256", "genesis_spec_sha256", "source_closure_sha256",
        "numerical_environment_sha256", "pref26_auth1_input_evidence_sha256",
    } or binding.get("artifact_id") != "FGC-1-PRO18-AUTH1" or binding.get("authorization_result_path") != AUTH_RESULT:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "receipt external authority differs")
    commit = binding["authorization_commit"]
    if not isinstance(commit, str) or len(commit) != 40 or commit.lower() != commit:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "receipt commit differs")
    try: int(commit, 16)
    except ValueError: _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "receipt commit differs")
    _git(root, "rev-parse", "--verify", f"{commit}^{{commit}}")
    try:
        subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=root, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (OSError, subprocess.CalledProcessError): _stop("PREF27_AUTHORITY_GIT_INVALID", "AUTH1 commit is not an ancestor")
    tree = _git(root, "ls-tree", "-r", "-z", commit, "--", AUTH_RESULT).split(b"\0")
    entries = [item for item in tree if item]
    if len(entries) != 1 or not entries[0].startswith(b"100644 blob "):
        _stop("PREF27_AUTHORITY_GIT_INVALID", "AUTH1 result is not a tracked regular blob")
    raw = _git(root, "show", f"{commit}:{AUTH_RESULT}")
    if sha256(raw).hexdigest() != binding["authorization_result_sha256"]:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "AUTH1 blob hash differs")
    live = _safe_file(root, AUTH_RESULT)
    if live != raw:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "live AUTH1 differs from recorded authority")
    record = _json(raw, label="AUTH1 result", indented=True)
    payload = record.get("artifact_payload")
    if record.get("artifact_id") != "FGC-1-PRO18-AUTH1" or not isinstance(payload, Mapping):
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "AUTH1 identity differs")
    for name in ("genesis_spec_sha256", "source_closure_sha256"):
        if payload.get(name) != binding[name]:
            _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", f"AUTH1 {name} differs")
    environment = payload.get("numerical_environment")
    if not isinstance(environment, Mapping) or environment.get("canonical_sha256") != binding["numerical_environment_sha256"]:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "AUTH1 environment differs")
    if payload.get("auth1_input_evidence_sha256") != binding["pref26_auth1_input_evidence_sha256"]:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "AUTH1 PREF26 evidence differs")
    return {"binding": dict(binding), "record": record, "payload": dict(payload)}


@dataclass(frozen=True, slots=True)
class Pref27Evidence:
    receipt: Mapping[str, Any]
    checkpoint: Mapping[str, Any]
    genesis_spec: Mapping[str, Any]
    authority_binding: Mapping[str, Any]
    leaves: tuple[Mapping[str, str], ...]
    tree_sha256: str


def verify_generation_zero_store(root: Path) -> Pref27Evidence:
    """Authenticate the eight-leaf, nonterminal GEN1 store without mutation."""
    root = Path(root)
    try: root_info = root.lstat()
    except OSError: _stop("PREF27_STORE_PATH_UNSAFE", "repository root is absent")
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        _stop("PREF27_STORE_PATH_UNSAFE", "repository root is unsafe")
    store = root / ROOT_RELATIVE
    try: info = store.lstat()
    except OSError: _stop("PREF27_STORE_PATH_UNSAFE", "calibration root is absent")
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode): _stop("PREF27_STORE_PATH_UNSAFE", "calibration root is unsafe")
    if any(item.name.startswith(".calibration.hlt15-stage-") for item in store.parent.iterdir()):
        _stop("PREF27_STORE_TREE_DRIFT", "stale HLT15 staging residue is present")
    receipt_raw = _safe_file(root, f"{ROOT_RELATIVE}/receipts/generation-zero.json")
    receipt = _json(receipt_raw, label="GEN1 receipt")
    bare = dict(receipt); receipt_digest = bare.pop("receipt_sha256", None)
    if receipt_digest != digest(bare) or receipt.get("claims") != EXPECTED_CLAIMS:
        _stop("PREF27_RECEIPT_DRIFT", "GEN1 receipt digest or claims differ")
    authority = _auth_record(root, receipt)
    spec = authority["payload"].get("genesis_spec")
    if not isinstance(spec, Mapping) or digest(spec) != authority["binding"]["genesis_spec_sha256"]:
        _stop("PREF27_AUTHORITY_RECEIPT_DRIFT", "AUTH1 GenesisSpec differs")
    try:
        reconstructed = build_genesis(spec)
    except (Proto17ConstructionError, TypeError, ValueError):
        _stop("PREF27_CHECKPOINT_DRIFT", "AUTH1 GenesisSpec does not reconstruct")
    checkpoint_hash = reconstructed["checkpoint"]["checkpoint_sha256"]
    checkpoint_rel = f"{ROOT_RELATIVE}/checkpoints/00000000000000000000-{checkpoint_hash}.json"
    checkpoint_raw = _safe_file(root, checkpoint_rel)
    checkpoint = _json(checkpoint_raw, label="GEN1 checkpoint")
    try:
        if validate_checkpoint(checkpoint) != reconstructed["checkpoint"]:
            _stop("PREF27_CHECKPOINT_DRIFT", "stored checkpoint differs from AUTH1 reconstruction")
    except Proto17ConstructionError:
        _stop("PREF27_CHECKPOINT_DRIFT", "stored checkpoint is invalid")
    if (receipt.get("checkpoint_sha256") != checkpoint_hash
            or receipt.get("state_object_set_sha256") != checkpoint["accepted_state_set_sha256"]
            or receipt.get("namespace_relative_path") != ROOT_RELATIVE):
        _stop("PREF27_RECEIPT_DRIFT", "receipt/checkpoint cross-binding differs")
    if _runtime(root, reconstructed["genesis_spec"], authority["binding"])._receipt(reconstructed["checkpoint"]) != receipt:
        _stop("PREF27_RECEIPT_DRIFT", "receipt differs from exact reconstructed runtime receipt")
    expected_paths = {checkpoint_rel, f"{ROOT_RELATIVE}/receipts/generation-zero.json"}
    leaves: list[dict[str, str]] = []
    for state_hash, state in reconstructed["store_plan"]["state_objects"].items():
        rel = f"{ROOT_RELATIVE}/states/{state_hash}.json"; expected_paths.add(rel)
        raw = _safe_file(root, rel)
        if sha256(raw).hexdigest() != state_hash or _json(raw, label=f"state {state_hash}") != state:
            _stop("PREF27_STATE_DRIFT", f"state differs: {state_hash}")
        leaves.append({"path": rel, "raw_sha256": sha256(raw).hexdigest(), "content_sha256": state_hash})
    observed: set[str] = set(); directories: set[str] = set()
    for path, item in _nonfollowing_entries(store, label="generation-zero store"):
        rel = path.relative_to(root).as_posix()
        if stat.S_ISREG(item.st_mode): observed.add(rel)
        elif stat.S_ISDIR(item.st_mode): directories.add(rel)
    if observed != expected_paths or directories != {
        f"{ROOT_RELATIVE}/checkpoints", f"{ROOT_RELATIVE}/receipts", f"{ROOT_RELATIVE}/states",
    }:
        _stop("PREF27_STORE_TREE_DRIFT", "store leaf grammar differs")
    leaves.extend([
        {"path": checkpoint_rel, "raw_sha256": sha256(checkpoint_raw).hexdigest(), "content_sha256": checkpoint_hash},
        {"path": f"{ROOT_RELATIVE}/receipts/generation-zero.json", "raw_sha256": sha256(receipt_raw).hexdigest(), "content_sha256": receipt_digest},
    ])
    ordered = tuple(sorted(leaves, key=lambda x: x["path"]))
    return Pref27Evidence(receipt=receipt, checkpoint=checkpoint, genesis_spec=reconstructed["genesis_spec"],
                           authority_binding=authority["binding"], leaves=ordered, tree_sha256=digest(list(ordered)))


def _runtime(root: Path, spec: Mapping[str, Any], binding: Mapping[str, Any]) -> Proto17HLT15Runtime:
    return Proto17HLT15Runtime(root, spec, authority_binding=binding)


def _temporary_controls(repository_root: Path, evidence: Pref27Evidence) -> Mapping[str, Any]:
    """Exercise no-overwrite, foreign, malformed, and race cases off-store only."""
    controls: dict[str, bool] = {}
    with tempfile.TemporaryDirectory(prefix="fgc-pref27-") as directory:
        base = Path(directory)
        def target(name: str) -> Path:
            return base / name / ROOT_RELATIVE
        def assert_refusal(name: str, populate) -> None:
            repo = base / name; path = target(name); path.parent.mkdir(parents=True, exist_ok=True); populate(path)
            before = _temp_tree(path)
            try: _runtime(repo, evidence.genesis_spec, evidence.authority_binding).materialize_generation_zero()
            except Proto17HLT15Error as error:
                if str(error) != "HLT15 refuses overwrite or resume of a namespace":
                    _stop("PREF27_CONTROL_FAILED", f"{name} wrong refusal: {error}")
            else: _stop("PREF27_CONTROL_FAILED", f"{name} was adopted")
            if before != _temp_tree(path): _stop("PREF27_CONTROL_FAILED", f"{name} mutated")
            if any(item.name.startswith(".calibration.hlt15-stage-") for item in path.parent.iterdir()): _stop("PREF27_CONTROL_FAILED", f"{name} left staging")
            controls[name] = True
        assert_refusal("partial", lambda p: (p.mkdir(), (p / "partial").write_text("x")))
        assert_refusal("file", lambda p: p.write_text("foreign"))
        assert_refusal("extra", lambda p: (p.mkdir(), (p / "extra").write_text("x")))
        assert_refusal("symlink", lambda p: os.symlink(base, p))
        # A byte-identical copy of the actual eight-leaf store is still a
        # pre-existing namespace and therefore must be refused, never resumed.
        copied_repo = base / "canonical-copy"; copied_target = target("canonical-copy")
        copied_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(Path(repository_root) / ROOT_RELATIVE, copied_target)
        before_copy = _temp_tree(copied_target)
        try: _runtime(copied_repo, evidence.genesis_spec, evidence.authority_binding).materialize_generation_zero()
        except Proto17HLT15Error as error:
            if str(error) != "HLT15 refuses overwrite or resume of a namespace": _stop("PREF27_CONTROL_FAILED", "copy wrong refusal")
        else: _stop("PREF27_CONTROL_FAILED", "canonical copied store was reused")
        if before_copy != _temp_tree(copied_target): _stop("PREF27_CONTROL_FAILED", "canonical copied store mutated")
        if any(item.name.startswith(".calibration.hlt15-stage-") for item in copied_target.parent.iterdir()): _stop("PREF27_CONTROL_FAILED", "canonical copied store left staging")
        controls["byte_identical_actual_store_copy"] = True
        controls["byte_identical_actual_store_copy_no_staging_residue"] = True
        # An internally valid but different GenesisSpec is still foreign and
        # must not be adopted by the original authority.
        repo = base / "foreign"; repo.mkdir(); foreign = deepcopy(dict(evidence.genesis_spec)); foreign["campaign_id"] = "FOREIGN-GENESIS"
        derived = build_genesis(foreign, _validate_derived=False)
        foreign_spec = derived["genesis_spec"]; foreign_binding = dict(evidence.authority_binding); foreign_binding["genesis_spec_sha256"] = digest(foreign_spec)
        _runtime(repo, foreign_spec, foreign_binding).materialize_generation_zero()
        before = _temp_tree(target("foreign"))
        try: _runtime(repo, evidence.genesis_spec, evidence.authority_binding).materialize_generation_zero()
        except Proto17HLT15Error as error:
            if str(error) != "HLT15 refuses overwrite or resume of a namespace": _stop("PREF27_CONTROL_FAILED", "foreign wrong refusal")
        else: _stop("PREF27_CONTROL_FAILED", "foreign valid store was adopted")
        if before != _temp_tree(target("foreign")): _stop("PREF27_CONTROL_FAILED", "foreign store mutated")
        if any(item.name.startswith(".calibration.hlt15-stage-") for item in target("foreign").parent.iterdir()): _stop("PREF27_CONTROL_FAILED", "foreign store left staging")
        controls["internally_self_consistent_foreign_store"] = True
        controls["internally_self_consistent_foreign_store_no_staging_residue"] = True
        # Deterministically model an arriving foreign target exactly between
        # private staging and Darwin's real RENAME_EXCL call.  The wrapper only
        # injects the arrival; the unmodified no-clobber primitive must reject.
        repo = base / "rename-race"; repo.mkdir(); arriving = target("rename-race")
        original_rename = hlt15_runtime._rename_exclusive
        def arriving_foreign(source: Path, destination: Path) -> None:
            destination.mkdir()
            (destination / "foreign-sentinel").write_bytes(b"PREF27-arriving-foreign")
            original_rename(source, destination)
        hlt15_runtime._rename_exclusive = arriving_foreign
        try:
            try: _runtime(repo, evidence.genesis_spec, evidence.authority_binding).materialize_generation_zero()
            except Proto17HLT15Error as error:
                if str(error) != "HLT15 namespace appeared during staging": _stop("PREF27_CONTROL_FAILED", f"race wrong refusal: {error}")
            else: _stop("PREF27_CONTROL_FAILED", "arriving foreign root was overwritten")
        finally:
            hlt15_runtime._rename_exclusive = original_rename
        if _temp_tree(arriving) != digest([("foreign-sentinel", sha256(b"PREF27-arriving-foreign").hexdigest())]):
            _stop("PREF27_CONTROL_FAILED", "arriving foreign root changed")
        if any(item.name.startswith(".calibration.hlt15-stage-") for item in arriving.parent.iterdir()):
            _stop("PREF27_CONTROL_FAILED", "arriving-root race left staging")
        controls["arriving_foreign_root_at_exclusive_rename"] = True
        controls["arriving_foreign_root_no_staging_residue"] = True
    return controls


def _temp_tree(path: Path) -> str:
    if path.is_symlink(): return "symlink:" + os.readlink(path)
    if path.is_file(): return "file:" + sha256(path.read_bytes()).hexdigest()
    rows: list[tuple[str, str]] = []
    for item, info in _nonfollowing_entries(path, label="temporary control store"):
        if stat.S_ISREG(info.st_mode):
            rows.append((item.relative_to(path).as_posix(), sha256(item.read_bytes()).hexdigest()))
    return digest(rows)


def build_pref27(root: Path, config: Mapping[str, Any] | None = None, *, config_sha256: str | None = None,
                 reproducer_sha256: str | None = None) -> dict[str, Any]:
    if config is not None: validate_pref27_config(config)
    before = verify_generation_zero_store(root)
    controls = _temporary_controls(root, before)
    after = verify_generation_zero_store(root)
    if before.tree_sha256 != after.tree_sha256:
        _stop("PREF27_REAL_STORE_MUTATED", "real store changed during temporary controls")
    implementation = {
        "src/recursive_horizons/fgc/evolution/proto18_pref27_binder.py": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scripts/reproduce_fgc_pro18_pref27.py": reproducer_sha256,
    }
    return {"schema_version": 1, "project_version": "0.11.0", "artifact_id": "FGC-1-PRO18-PREF27",
            "classification": "completed_read_only_generation_zero_store_binder", "source_config_sha256": config_sha256,
            "implementation_sha256": implementation,
            "gate_status": {"state_advance_authorized": False, "candidate_execution_authorized": False},
            "nonclaims": {"holdout_opened": False, "state_advanced": False, "physical_transition_claim": False},
            "artifact_payload": {
        "authority_binding": before.authority_binding,
        "receipt_sha256": before.receipt["receipt_sha256"],
        "checkpoint_sha256": before.checkpoint["checkpoint_sha256"],
        "store_leaves": list(before.leaves), "store_tree_sha256": before.tree_sha256,
        "temporary_refusal_controls": controls,
        "claims": {"HLT15_generation_zero_store_authenticated": True,
                   "production_namespace_reuse_refusal_verified": True,
                   "foreign_store_adoption_refusal_verified": True,
                   "all_four_HLT14_production_duties_completed": True,
                   "state_advanced": False, "fresh_GR0_dynamic_calibration_authorized": False,
                   "candidate_execution_authorized": False, "physical_transition_claim_authorized": False},
    }}
