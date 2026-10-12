#!/usr/bin/env python3
"""Reproduce the compact PROTO16 common-event/trusted-genesis freeze."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import tomllib
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from recursive_horizons.fgc.evolution.protocol_v16 import (
    FALSE_CLAIMS,
    TRUE_CLAIMS,
    validate_sf1_protocol_v16,
)

ARTIFACT_ID = "FGC-1-PRO16-FRZ1"
CONFIG = ROOT / "configs/fgc/fgc-1-pro16-frz1.toml"
OUTPUT = ROOT / "results/fgc-1-pro16-frz1.json"
DOC = ROOT / "docs/fgc-pro16-frz1.md"
LINEAGE = {
    "configs/fgc/fgc-2-sf1-protocol-v15.toml": "f1f1849d85b6cbe1eb61267392a63ca11389f551df127033ddf905abb38fdc3b",
    "configs/fgc/fgc-1-pro15-frz1.toml": "1a45af49cca92701a79cf103fd6795ffebf6353af372bff779cd4aae90f71676",
    "results/fgc-1-pro15-frz1.json": "99d24425cb5d40da2e248605265b83bcdc82a2126cfa4b5bce4239c6ef30635c",
    "results/fgc-1-hlt13-mon13.json": "de55f8d55a1db0043945d3553c7d2070703809b82c1f7141bf923c036b9df5b7",
    "src/recursive_horizons/fgc/evolution/proto15_runtime.py": "73f07ae26e90642bc6e5158d9a85bfa499d8082c13e082fc8f19e968d8e624bc",
}


def sha(p: Path) -> str:
    return sha256(p.read_bytes()).hexdigest()


def canonical(v: Any) -> str:
    return (
        json.dumps(v, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        + "\n"
    )


def reject(pairs):
    d = {}
    for k, v in pairs:
        if k in d:
            raise ValueError("duplicate JSON key")
        d[k] = v
    return d


def load_config(path=CONFIG):
    with path.open("rb") as f:
        v = tomllib.load(f)
    required = {
        "schema_version",
        "artifact_id",
        "project_version",
        "protocol_config",
        "protocol_module",
        "scope",
        "immutable_lineage",
        "future_authority",
        "namespace_precondition",
        "proof_contract",
        "claims",
    }
    if (
        set(v) != required
        or (v["schema_version"], v["artifact_id"], v["project_version"])
        != (1, ARTIFACT_ID, "0.11.0")
        or v["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v16.toml"
        or v["protocol_module"]
        != "src/recursive_horizons/fgc/evolution/protocol_v16.py"
    ):
        raise ValueError("PROTO16 freeze config identity differs")
    if v["scope"] != {
        "target_protocol": "FGC-2-SF1-PROTO16",
        "predecessor_protocol": "FGC-2-SF1-PROTO15",
        "successor_runtime_owner": "FGC-1-HLT14-MON14",
        "role": "premise_only_common_event_commit_and_external_trusted_genesis_contract_certificate",
        "PROTO16_trajectory_read": False,
        "historical_raw_campaign_or_checkpoint_opened": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "mechanism_question_answered": False,
    }:
        raise ValueError("scope differs")
    if (
        v["claims"].keys() != (TRUE_CLAIMS | FALSE_CLAIMS)
        or any(v["claims"][k] is not True for k in TRUE_CLAIMS)
        or any(v["claims"][k] is not False for k in FALSE_CLAIMS)
    ):
        raise ValueError("claims differ")
    if (
        any(x is not True for x in v["proof_contract"].values())
        or v["future_authority"].get("authority_artifact_id") != "FGC-1-HLT14-MON14"
        or v["future_authority"].get("authority_result_is_not_yet_created") is not True
    ):
        raise ValueError("future authority differs")
    return v


def record(config_path=CONFIG):
    c = load_config(config_path)
    commit = c["immutable_lineage"]["evidence_checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=ROOT
    ).returncode:
        raise ValueError("sealed HLT13 commit is not ancestor")
    records = []
    for rel, digest in LINEAGE.items():
        b = subprocess.run(
            ["git", "show", f"{commit}:{rel}"],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            check=True,
        ).stdout
        if sha256(b).hexdigest() != digest or sha(ROOT / rel) != digest:
            raise ValueError(f"lineage differs: {rel}")
        records.append(
            {"path": rel, "sha256": digest, "matches_sealed_commit_and_worktree": True}
        )
    with (ROOT / c["protocol_config"]).open("rb") as f:
        protocol = tomllib.load(f)
    validate_sf1_protocol_v16(protocol)
    for rel in c["namespace_precondition"].values():
        if isinstance(rel, str) and (ROOT / rel).exists():
            raise ValueError("PROTO16 namespace exists")
    hlt = json.loads(
        (ROOT / "results/fgc-1-hlt13-mon13.json").read_text(), object_pairs_hook=reject
    )
    if (
        hlt.get("gate_status", {}).get("PROTO15_runtime_synthetic_qualification_passed")
        is not True
    ):
        raise ValueError("HLT13 not qualified")
    mutations = {}
    for name, edit in {
        "event_kind": lambda x: x["event_commit"].__setitem__(
            "journal_record_kind", "OTHER"
        ),
        "receipt_cycle": lambda x: x["event_commit"].__setitem__(
            "receipt_payload_required_fields", ["common_event_receipt_sha256"]
        ),
        "cadence": lambda x: x["event_commit"].__setitem__(
            "frozen_common_event_cadence", "1/8"
        ),
        "claim": lambda x: x["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
        "genesis": lambda x: x["trusted_genesis"].__setitem__(
            "authority_artifact_id", "OTHER"
        ),
        "launch_tuple": lambda x: x["trusted_genesis"].__setitem__(
            "launch_authority_tuple_required_fields", []
        ),
        "genesis_self_hash": lambda x: x["trusted_genesis"].__setitem__(
            "genesis_spec_required_fields", ["genesis_spec_sha256"]
        ),
        "hidden_member_descriptors": lambda x: x["trusted_genesis"].__setitem__(
            "member_descriptor_required_fields", []
        ),
    }.items():
        a = deepcopy(protocol)
        edit(a)
        try:
            validate_sf1_protocol_v16(a)
        except (TypeError, ValueError):
            mutations[name] = True
        else:
            mutations[name] = False
    if not all(mutations.values()):
        raise ValueError("mutation control failed")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": "0.11.0",
        "classification": "prospective_common_event_commit_and_external_trusted_genesis_protocol_freeze_without_runtime_or_trajectory",
        "generated_by": "scripts/reproduce_fgc_pro16_frz1.py",
        "derivation_document": "docs/fgc-pro16-frz1.md",
        "derivation_document_sha256": sha(DOC),
        "source_config_sha256": {
            "configs/fgc/fgc-1-pro16-frz1.toml": sha(config_path),
            "configs/fgc/fgc-2-sf1-protocol-v16.toml": sha(ROOT / c["protocol_config"]),
        },
        "implementation_sha256": {
            "scripts/reproduce_fgc_pro16_frz1.py": sha(Path(__file__)),
            "src/recursive_horizons/fgc/evolution/protocol_v16.py": sha(
                ROOT / c["protocol_module"]
            ),
        },
        "scope_bindings": c["scope"],
        "artifact_payload": {
            "immutable_lineage": {
                "evidence_checkpoint_commit": commit,
                "checkpoint_is_ancestor_of_HEAD": True,
                "tracked_compact_records": records,
            },
            "event_commit_contract": protocol["event_commit"],
            "trusted_genesis_contract": protocol["trusted_genesis"],
            "restart_descriptor_contract": protocol["restart_inputs"],
            "namespace_precondition": {
                "both_new_namespaces_absent": True,
                "freeze_created_no_namespace": True,
                "records": [
                    {
                        "path": c["namespace_precondition"]["calibration_output_root"],
                        "absent_before_freeze": True,
                    },
                    {
                        "path": c["namespace_precondition"]["holdout_output_root"],
                        "absent_before_freeze": True,
                    },
                ],
            },
            "mutation_controls": mutations,
            "decision": {
                "PROTO16_frozen": True,
                "successor_runtime_owner": "FGC-1-HLT14-MON14",
                "successor_runtime_implemented": False,
                "pretrajectory_authorized": False,
                "calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
        },
        "gate_status": c["claims"],
        "nonclaims": {k: False for k in FALSE_CLAIMS},
    }


def verify(path=OUTPUT):
    src = path.read_text()
    v = json.loads(src, object_pairs_hook=reject)
    if src != canonical(v) or v != record():
        raise ValueError("PROTO16 result differs")


def atomic(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(payload)
            f.flush()
            os.fsync(f.fileno())
        Path(tmp).replace(path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    p.add_argument("--output", type=Path, default=OUTPUT)
    a = p.parse_args()
    if a.verify:
        verify(a.output)
        print(f"verified {a.output}")
    else:
        atomic(a.output, canonical(record()).encode())
        print(f"wrote {a.output}")
