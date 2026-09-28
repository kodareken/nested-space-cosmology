"""Construct the non-self-authenticating PRO18 AUTH1 genesis authority.

AUTH1 is deliberately a content artifact.  It binds already compact PREF26
evidence, pins the prospective production implementation and derives one
calibration-only generation-zero ``GenesisSpec``.  It does *not* authenticate
the resulting file through Git, import the evolution implementation, observe
either future namespace, or write a runtime store.  Those are HLT15/PREF27
responsibilities.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .proto17_pure_construction import (
    GENESIS_INPUT_FIELDS,
    MEMBER_KEYS,
    Proto17ConstructionError,
    build_genesis,
    canonical,
    digest,
)
from .proto18_auth1_inputs import Proto18Auth1InputError, build_calibration_genesis


class Proto18AuthorityError(ValueError):
    """A required AUTH1 construction premise is false."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


DIRECT_ROLE_PATHS = {
    "runtime_module": "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py",
    "adapter_module": "src/recursive_horizons/fgc/evolution/proto18_auth1_inputs.py",
    "runner": "scripts/run_fgc_proto17_hlt15_gen1.py",
}
REQUIRED_CLOSURE_PATHS = frozenset({
    "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py",
    "src/recursive_horizons/fgc/evolution/proto18_auth1_inputs.py",
    "scripts/run_fgc_proto17_hlt15_gen1.py",
    "src/recursive_horizons/fgc/evolution/proto18_production_inputs.py",
    "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py",
    "src/recursive_horizons/fgc/evolution/protocol_v17.py",
    "src/recursive_horizons/fgc/evolution/numerical_engine.py",
    "src/recursive_horizons/fgc/evolution/proto18_authority.py",
})
PREF26_COMMIT = "ab3dd2613f70028388e0658c873c33b4225233c7"
PREF26_RESULT = "results/fgc-1-pro18-pref26.json"
AUTH1_CONFIG = "configs/fgc/fgc-1-pro18-auth1.toml"
RAW_SOURCE_CONFIG = "configs/fgc/fgc-1-pro18-frz1.toml"

_EXPECTED_TRUE_CLAIMS = frozenset({
    "canonical_payload_constructed",
    "PREF26_evidence_bound",
    "complete_calibration_GenesisSpec_constructed",
    "calibration_launch_closure_pinned",
})
_EXPECTED_FALSE_CLAIMS = frozenset({
    "committed_external_Git_live_import_validation_completed",
    "HLT15_authorized",
    "future_output_roots_observed",
    "namespace_created",
    "pretrajectory_operation_authorized",
    "fresh_GR0_dynamic_calibration_completed",
    "candidate_execution_authorized",
    "physical_transition_claim_authorized",
})


def derive_campaign_id(*, pref26_evidence_sha256: str, run_plan_sha256: str) -> str:
    """Derive the one PROTO17 calibration identity from frozen inputs.

    The identifier is not reused from an old campaign and avoids ``:`` because
    cursor attempt IDs reserve that delimiter.  The short suffix is a readable
    projection of a full canonical derivation recorded in the AUTH1 payload.
    """
    material = {
        "target_protocol": "FGC-2-SF1-PROTO17",
        "branch": "GR-0",
        "amplitude": "3",
        "namespace": "runs/fgc-2-sf1/proto17/calibration",
        "pref26_auth1_input_evidence_sha256": _sha(pref26_evidence_sha256, label="PREF26 evidence"),
        "run_plan_sha256": _sha(run_plan_sha256, label="run plan"),
    }
    return f"FGC-2-SF1-PROTO17-GR0-A3-CAL-{digest(material)[:20]}"


def _canonical_json(raw: bytes, *, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result
    try:
        parsed = json.loads(raw.decode("utf-8"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", f"{label} is malformed") from error
    # PREF26 defines its compact evidence surface as sorted, indented JSON
    # (rather than the pure constructor's compact JSON), so authenticate that
    # owner format exactly before consuming it.
    expected = (json.dumps(parsed, sort_keys=True, indent=2, ensure_ascii=True,
                           allow_nan=False) + "\n").encode("utf-8")
    if not isinstance(parsed, dict) or expected != raw:
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", f"{label} is not canonical JSON")
    return parsed


def _sha(value: object, *, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", f"{label} is not lowercase SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", f"{label} is not hexadecimal") from error
    return value


def _relative_path(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"{label} is absent")
    path = Path(value)
    if path.is_absolute() or "\\" in value or any(part in {"", ".", ".."} for part in path.parts):
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"{label} is unsafe")
    return value


def _module_for_path(relative: str) -> str | None:
    if not relative.startswith("src/") or not relative.endswith(".py"):
        return None
    return relative.removeprefix("src/").removesuffix(".py").replace("/", ".")


def _path_for_module(module: str) -> str | None:
    if not module.startswith("recursive_horizons."):
        return None
    return f"src/{module.replace('.', '/')}.py"


def _local_imports(relative: str, source: bytes) -> tuple[str, ...]:
    """Return project-local import targets without importing the module."""
    try:
        tree = ast.parse(source, filename=relative)
    except SyntaxError as error:
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"{relative} does not parse") from error
    current = _module_for_path(relative)
    package = current.rsplit(".", 1)[0] if current and "." in current else ""
    found: set[str] = set()
    for node in ast.walk(tree):
        module: str | None = None
        if isinstance(node, ast.Import):
            for alias in node.names:
                target = _path_for_module(alias.name)
                if target:
                    found.add(target)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if not package:
                    continue
                pieces = package.split(".")
                if node.level > len(pieces) + 1:
                    continue
                base = ".".join(pieces[: len(pieces) - node.level + 1])
                module = f"{base}.{node.module}" if node.module else base
            else:
                module = node.module
            if module:
                target = _path_for_module(module)
                if target:
                    found.add(target)
                # ``from .foo import X`` points directly at foo.  ``from .
                # import foo`` needs the imported names checked as modules.
                if node.module is None:
                    for alias in node.names:
                        sibling = _path_for_module(f"{module}.{alias.name}")
                        if sibling:
                            found.add(sibling)
    return tuple(sorted(found))


def derive_source_closure(root: Path, direct_paths: Mapping[str, str]) -> tuple[str, ...]:
    """Dependency-first ordered closure for declared project-local sources."""
    root = Path(root)
    ordered: list[str] = []
    visiting: set[str] = set()
    visited: set[str] = set()

    def walk(relative: str) -> None:
        relative = _relative_path(relative, label="source closure path")
        if relative in visited:
            return
        if relative in visiting:
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "project-local import cycle")
        candidate = root / relative
        if not candidate.is_file() or candidate.is_symlink():
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"declared source is absent: {relative}")
        visiting.add(relative)
        for child in _local_imports(relative, candidate.read_bytes()):
            child_path = root / child
            if child_path.exists():
                walk(child)
            elif child.startswith("src/recursive_horizons/"):
                raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"local import is absent: {child}")
        visiting.remove(relative)
        visited.add(relative)
        ordered.append(relative)

    if set(direct_paths) != set(DIRECT_ROLE_PATHS):
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "direct role vocabulary differs")
    for role in ("runtime_module", "adapter_module", "runner"):
        if direct_paths[role] != DIRECT_ROLE_PATHS[role]:
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"{role} path differs")
        walk(direct_paths[role])
    return tuple(ordered)


def _load_pref26_evidence(root: Path, config: Mapping[str, Any]) -> Mapping[str, Any]:
    predecessor = config.get("pref26_predecessor")
    if not isinstance(predecessor, Mapping):
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 predecessor declaration is absent")
    required = {"commit", "result", "result_sha256", "auth1_input_evidence_sha256"}
    if set(predecessor) != required or predecessor.get("commit") != PREF26_COMMIT or predecessor.get("result") != PREF26_RESULT:
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 predecessor identity differs")
    path = root / PREF26_RESULT
    if not path.is_file() or path.is_symlink():
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 compact result is absent")
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != _sha(predecessor["result_sha256"], label="PREF26 result"):
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 compact result hash differs")
    result = _canonical_json(raw, label="PREF26 compact result")
    payload = result.get("artifact_payload")
    if result.get("artifact_id") != "FGC-1-PRO18-PREF26" or not isinstance(payload, Mapping):
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 result identity differs")
    evidence = payload.get("AUTH1_input_evidence")
    if not isinstance(evidence, Mapping):
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 AUTH1 input evidence is absent")
    evidence = dict(evidence)
    expected = _sha(predecessor["auth1_input_evidence_sha256"], label="PREF26 AUTH1 evidence")
    if evidence.get("AUTH1_input_evidence_sha256") != expected:
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 AUTH1 evidence self digest differs")
    bare = dict(evidence); bare.pop("AUTH1_input_evidence_sha256", None)
    if digest(bare) != expected:
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 AUTH1 evidence content differs")
    if (evidence.get("authority_identity"), evidence.get("target_protocol"),
            evidence.get("final_PROTO17_GenesisSpec_construction_deferred_to_AUTH1"),
            evidence.get("future_runtime_adapter_and_runner_source_hashes_present")) != (
        "FGC-1-PRO18-AUTH1", "FGC-2-SF1-PROTO17", True, False,
    ):
        raise Proto18AuthorityError("AUTH1_PREF26_EVIDENCE_INVALID", "PREF26 evidence scope differs")
    return evidence


def _claims(config: Mapping[str, Any]) -> dict[str, bool]:
    claims = config.get("claims")
    if not isinstance(claims, Mapping) or set(claims) != _EXPECTED_TRUE_CLAIMS | _EXPECTED_FALSE_CLAIMS:
        raise Proto18AuthorityError("AUTH1_CLAIM_DRIFT", "AUTH1 claim vocabulary differs")
    if any(claims[name] is not True for name in _EXPECTED_TRUE_CLAIMS) or any(
        claims[name] is not False for name in _EXPECTED_FALSE_CLAIMS
    ):
        raise Proto18AuthorityError("AUTH1_CLAIM_DRIFT", "AUTH1 claim values differ")
    return dict(claims)


def _hash_bound_file(root: Path, declaration: Mapping[str, Any], *, label: str) -> str:
    if not isinstance(declaration, Mapping) or set(declaration) != {"path", "sha256"}:
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", f"{label} declaration differs")
    relative = _relative_path(declaration["path"], label=label)
    expected = _sha(declaration["sha256"], label=label)
    path = root / relative
    if not path.is_file() or path.is_symlink() or sha256(path.read_bytes()).hexdigest() != expected:
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", f"{label} bytes differ")
    return expected


def _environment(config: Mapping[str, Any]) -> Mapping[str, Any]:
    environment = config.get("numerical_environment")
    if not isinstance(environment, Mapping) or set(environment) != {"contract", "canonical_sha256"}:
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "environment declaration differs")
    contract = environment["contract"]
    expected_keys = {
        "contract_id", "implementation", "python_version", "python_build", "python_compiler",
        "numpy_version", "numpy_build_configuration", "blas_lapack", "operating_system",
        "system", "machine", "processor", "platform", "byteorder",
    }
    if not isinstance(contract, Mapping) or set(contract) != expected_keys:
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "environment contract vocabulary differs")
    if (contract.get("implementation"), contract.get("python_version"), contract.get("numpy_version"),
            contract.get("blas_lapack"), contract.get("system"), contract.get("machine")) != (
        "CPython", "3.14.3", "2.5.1", "Accelerate", "Darwin", "arm64",
    ):
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "evolution environment identity differs")
    expected = _sha(environment["canonical_sha256"], label="environment")
    if digest(dict(contract)) != expected:
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "environment canonical digest differs")
    return {"contract": dict(contract), "canonical_sha256": expected}


def _source_pins(root: Path, config: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    pins = config.get("source_closure")
    if not isinstance(pins, list) or not pins:
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "source closure is absent")
    declared: list[dict[str, Any]] = []
    direct: dict[str, str] = {}
    for entry in pins:
        if not isinstance(entry, Mapping) or set(entry) != {"path", "sha256", "role"}:
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "source pin schema differs")
        relative = _relative_path(entry["path"], label="source pin")
        role = entry["role"]
        if role not in {"runtime_module", "adapter_module", "runner", "dependency"}:
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "source pin role differs")
        if role != "dependency":
            if role in direct or relative != DIRECT_ROLE_PATHS[role]:
                raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "direct source role differs")
            direct[role] = relative
        expected = _sha(entry["sha256"], label=f"source pin {relative}")
        candidate = root / relative
        if not candidate.is_file() or candidate.is_symlink() or sha256(candidate.read_bytes()).hexdigest() != expected:
            raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", f"source bytes differ: {relative}")
        declared.append({"path": relative, "sha256": expected, "role": role})
    if set(direct) != set(DIRECT_ROLE_PATHS) or len({pin["path"] for pin in declared}) != len(declared):
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "direct source closure is incomplete or duplicated")
    observed = derive_source_closure(root, direct)
    if tuple(pin["path"] for pin in declared) != observed:
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "declared source closure differs from AST closure")
    if not REQUIRED_CLOSURE_PATHS.issubset(set(observed)):
        raise Proto18AuthorityError("AUTH1_SOURCE_CLOSURE_DRIFT", "required production closure is incomplete")
    return tuple(declared)


def build_authority(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    """Build the content-only AUTH1 result and exact legal GenesisSpec."""
    if not isinstance(config, Mapping) or (config.get("schema_version"), config.get("artifact_id"), config.get("target_protocol")) != (
        1, "FGC-1-PRO18-AUTH1", "FGC-2-SF1-PROTO17",
    ):
        raise Proto18AuthorityError("AUTH1_CONTRACT_DRIFT", "AUTH1 identity differs")
    evidence = _load_pref26_evidence(Path(root), config)
    claims = _claims(config)
    environment = _environment(config)
    source_pins = _source_pins(Path(root), config)
    runtime_bindings = config.get("runtime_bindings")
    if (not isinstance(runtime_bindings, Mapping)
            or set(runtime_bindings) != {"raw_source_config"}):
        raise Proto18AuthorityError(
            "AUTH1_GENESIS_INPUT_DRIFT", "runtime input binding vocabulary differs",
        )
    raw_source_config_sha256 = _hash_bound_file(
        Path(root), runtime_bindings["raw_source_config"], label="raw source config",
    )
    if runtime_bindings["raw_source_config"]["path"] != RAW_SOURCE_CONFIG:
        raise Proto18AuthorityError(
            "AUTH1_GENESIS_INPUT_DRIFT", "raw source config path differs",
        )
    authority_config_path = Path(root) / AUTH1_CONFIG
    if (not authority_config_path.is_file() or authority_config_path.is_symlink()):
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", "AUTH1 config is absent")
    authority_config_sha256 = sha256(authority_config_path.read_bytes()).hexdigest()
    bindings = config.get("genesis_bindings")
    required_bindings = {
        "protocol_config", "protocol_freeze_result", "hlt13_result", "run_plan",
        "campaign_id", "namespace",
    }
    if not isinstance(bindings, Mapping) or set(bindings) != required_bindings:
        raise Proto18AuthorityError("AUTH1_GENESIS_INPUT_DRIFT", "genesis binding vocabulary differs")
    # This strings-only check intentionally does not observe either path.
    if bindings["namespace"] != "runs/fgc-2-sf1/proto17/calibration":
        raise Proto18AuthorityError("AUTH1_HOLDOUT_NAMESPACE_FORBIDDEN", "only calibration namespace may appear")
    hashes = {
        name: _hash_bound_file(Path(root), bindings[name], label=name)
        for name in ("protocol_config", "protocol_freeze_result", "hlt13_result", "run_plan")
    }
    expected_campaign_id = derive_campaign_id(
        pref26_evidence_sha256=str(evidence["AUTH1_input_evidence_sha256"]),
        run_plan_sha256=hashes["run_plan"],
    )
    if bindings["campaign_id"] != expected_campaign_id:
        raise Proto18AuthorityError("AUTH1_GENESIS_INPUT_DRIFT", "campaign identity derivation differs")
    direct_hashes = {pin["role"]: pin["sha256"] for pin in source_pins if pin["role"] != "dependency"}
    try:
        constructed = build_calibration_genesis(
            evidence,
            campaign_id=bindings["campaign_id"],
            namespace=bindings["namespace"],
            source_hashes={
                "protocol_config_sha256": hashes["protocol_config"],
                "protocol_freeze_result_sha256": hashes["protocol_freeze_result"],
                "hlt13_result_sha256": hashes["hlt13_result"],
                "run_plan_sha256": hashes["run_plan"],
                "runtime_module_sha256": direct_hashes["runtime_module"],
                "adapter_module_sha256": direct_hashes["adapter_module"],
                "runner_sha256": direct_hashes["runner"],
            },
            numerical_environment=environment["contract"],
        )
    except (Proto18Auth1InputError, Proto17ConstructionError, TypeError, ValueError) as error:
        raise Proto18AuthorityError("AUTH1_GENESIS_INPUT_DRIFT", "GenesisSpec construction failed") from error
    genesis = constructed.genesis_spec
    result = {
        "artifact_id": "FGC-1-PRO18-AUTH1",
        "artifact_payload": {
            "pref26_binding": {
                "historical_commit": PREF26_COMMIT,
                "compact_result": PREF26_RESULT,
                "compact_result_sha256": config["pref26_predecessor"]["result_sha256"],
                "AUTH1_input_evidence_sha256": evidence["AUTH1_input_evidence_sha256"],
            },
            "authority_inputs": {
                "authority_config": {
                    "path": AUTH1_CONFIG,
                    "sha256": authority_config_sha256,
                },
                "raw_source_config": {
                    "path": RAW_SOURCE_CONFIG,
                    "sha256": raw_source_config_sha256,
                },
            },
            "genesis_spec": genesis,
            "genesis_spec_sha256": digest(genesis),
            "auth1_input_evidence": evidence,
            "auth1_input_evidence_sha256": evidence["AUTH1_input_evidence_sha256"],
            "source_hashes": {
                "protocol_config_sha256": hashes["protocol_config"],
                "protocol_freeze_result_sha256": hashes["protocol_freeze_result"],
                "hlt13_result_sha256": hashes["hlt13_result"],
                "run_plan_sha256": hashes["run_plan"],
                "runtime_module_sha256": direct_hashes["runtime_module"],
                "adapter_module_sha256": direct_hashes["adapter_module"],
                "runner_sha256": direct_hashes["runner"],
            },
            "campaign_id": expected_campaign_id,
            "source_closure": list(source_pins),
            "source_closure_sha256": digest(list(source_pins)),
            "numerical_environment": environment,
            "claims": claims,
        },
    }
    # The content payload must remain free of the external launch tuple.
    forbidden = {"authorization_commit", "authorization_result_path", "authorization_result_sha256"}
    def contains(value: object) -> bool:
        if isinstance(value, Mapping):
            return any(key in forbidden or contains(item) for key, item in value.items())
        return isinstance(value, list) and any(contains(item) for item in value)
    if contains(result):
        raise AssertionError("AUTH1 payload accidentally self-authenticates")
    return result


def verify_authority_payload(root: Path, config: Mapping[str, Any], supplied: Mapping[str, Any]) -> dict[str, Any]:
    """Reconstruct AUTH1 and reject every payload or GenesisSpec mutation.

    This is intentionally content-only verification: it does not read Git,
    import production modules, touch raw predecessors, or inspect a future
    namespace.  HLT15 later authenticates a committed external tuple.
    """
    expected = build_authority(root, config)
    if not isinstance(supplied, Mapping) or set(supplied) != set(expected):
        raise Proto18AuthorityError("AUTH1_CLAIM_DRIFT", "authority result schema differs")
    if supplied.get("artifact_id") != expected["artifact_id"]:
        raise Proto18AuthorityError("AUTH1_CLAIM_DRIFT", "authority result identity differs")
    payload = supplied.get("artifact_payload")
    expected_payload = expected["artifact_payload"]
    if not isinstance(payload, Mapping) or set(payload) != set(expected_payload):
        raise Proto18AuthorityError("AUTH1_CLAIM_DRIFT", "authority payload vocabulary differs")
    for field in GENESIS_INPUT_FIELDS:
        if payload.get("genesis_spec", {}).get(field) != expected_payload["genesis_spec"][field]:
            raise Proto18AuthorityError("AUTH1_GENESIS_INPUT_DRIFT", f"GenesisSpec field differs: {field}")
    if dict(payload) != expected_payload:
        raise Proto18AuthorityError("AUTH1_HASH_DRIFT", "AUTH1 canonical payload differs")
    return expected


__all__ = [
    "AUTH1_CONFIG", "DIRECT_ROLE_PATHS", "PREF26_COMMIT", "PREF26_RESULT",
    "RAW_SOURCE_CONFIG", "Proto18AuthorityError",
    "build_authority", "derive_campaign_id", "derive_source_closure", "verify_authority_payload",
]
