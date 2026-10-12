#!/usr/bin/env python3
"""Reproduce the independent compact PROTO16 transition theorem."""

from __future__ import annotations
import argparse, json, os, subprocess, tempfile, tomllib
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto16_transition_theorem import audit

CONFIG = ROOT / "configs/fgc/fgc-1-pro16-pref24.toml"
OUTPUT = ROOT / "results/fgc-1-pro16-pref24.json"
DOC = ROOT / "docs/fgc-pro16-pref24.md"
LINEAGE = {
    "configs/fgc/fgc-2-sf1-protocol-v16.toml": "3b0ab8663ed0568f03841167a2db096419b0978a71b534eefa9f9bb9a6787c75",
    "configs/fgc/fgc-1-pro16-frz1.toml": "5c296d158e02b40faf576be3ff3261f9ccf76f9c6f01025c489e6b0aba5aaea4",
    "results/fgc-1-pro16-frz1.json": "f8e14deda1638bda6e98831ae1ea80c8d1fdffb23957e183acd8b57a154ce43f",
    "src/recursive_horizons/fgc/evolution/protocol_v16.py": "c1755e917649a9e1d81ff5fe250847d5b9bc9e273d2c3cc6dba689119d6b0318",
    "results/fgc-1-hlt13-mon13.json": "de55f8d55a1db0043945d3553c7d2070703809b82c1f7141bf923c036b9df5b7",
    "src/recursive_horizons/fgc/evolution/proto15_runtime.py": "73f07ae26e90642bc6e5158d9a85bfa499d8082c13e082fc8f19e968d8e624bc",
}
TRUE = {
    "PROTO16_receipt_hash_cycle_removed",
    "PROTO16_exact_successor_and_genesis_derivation_not_frozen",
    "PROTO17_successor_freeze_design_authorized",
}


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def canon(x):
    return (
        json.dumps(x, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        + "\n"
    )


def reject(pairs):
    out = {}
    for k, v in pairs:
        if k in out:
            raise ValueError("duplicate key")
        out[k] = v
    return out


def record():
    with CONFIG.open("rb") as f:
        c = tomllib.load(f)
    if (
        c["immutable_lineage"]["checkpoint_commit"]
        != "8a1bb20ba43484685dad7eff39e300d61504b01c"
        or c["immutable_lineage"].get("proto15_runtime_sha256")
        != LINEAGE["src/recursive_horizons/fgc/evolution/proto15_runtime.py"]
        or subprocess.run(
            [
                "git",
                "merge-base",
                "--is-ancestor",
                c["immutable_lineage"]["checkpoint_commit"],
                "HEAD",
            ],
            cwd=ROOT,
        ).returncode
    ):
        raise ValueError("lineage differs")
    rows = []
    for path, d in LINEAGE.items():
        old = subprocess.run(
            ["git", "show", f"{c['immutable_lineage']['checkpoint_commit']}:{path}"],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            check=True,
        ).stdout
        if sha256(old).hexdigest() != d or sha(ROOT / path) != d:
            raise ValueError("sealed byte drift")
        rows.append({"path": path, "sha256": d, "matches_commit_and_worktree": True})
    with (ROOT / c["protocol_config"]).open("rb") as handle:
        protocol = tomllib.load(handle)
    t = audit(protocol)
    if (
        not t["receipt"]["receipt_acyclic"]
        or t["successor_schema"]["exact_derivation_frozen"]
        or t["genesis_schema"]["exact_derivation_frozen"]
    ):
        raise ValueError("theorem failed")
    claims = c["claims"]
    if any(claims[x] is not True for x in TRUE) or any(
        v is not False for k, v in claims.items() if k not in TRUE
    ):
        raise ValueError("claim boundary")
    return {
        "schema_version": 1,
        "artifact_id": "FGC-1-PRO16-PREF24",
        "project_version": "0.11.0",
        "classification": "independent_PROTO16_hash_cycle_and_under_specification_diagnosis_without_runtime_or_trajectory",
        "generated_by": "scripts/reproduce_fgc_pro16_pref24.py",
        "derivation_document": "docs/fgc-pro16-pref24.md",
        "derivation_document_sha256": sha(DOC),
        "source_config_sha256": {"configs/fgc/fgc-1-pro16-pref24.toml": sha(CONFIG)},
        "implementation_sha256": {
            "scripts/reproduce_fgc_pro16_pref24.py": sha(Path(__file__)),
            "src/recursive_horizons/fgc/evolution/proto16_transition_theorem.py": sha(
                ROOT
                / "src/recursive_horizons/fgc/evolution/proto16_transition_theorem.py"
            ),
        },
        "scope_bindings": c["scope"],
        "artifact_payload": {
            "immutable_lineage": {
                "checkpoint_commit": c["immutable_lineage"]["checkpoint_commit"],
                "records": rows,
            },
            "static_completeness_audit": t,
            "decision": {
                "PROTO17_successor_freeze_design_authorized": True,
                "HLT14_implementation_and_synthetic_qualification_authorized": False,
                "runtime_implemented": False,
                "pretrajectory_authorized": False,
                "calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
        },
        "gate_status": claims,
        "nonclaims": {k: False for k in claims if k not in TRUE},
    }


def verify(path=OUTPUT):
    s = path.read_text()
    v = json.loads(s, object_pairs_hook=reject)
    if s != canon(v) or v != record():
        raise ValueError("canonical result differs")


def write(p, data):
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=p.parent)
    with os.fdopen(fd, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    Path(tmp).replace(p)


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--verify", action="store_true")
    a.add_argument("--output", type=Path, default=OUTPUT)
    x = a.parse_args()
    if x.verify:
        verify(x.output)
        print(f"verified {x.output}")
    else:
        write(x.output, canon(record()).encode())
        print(f"wrote {x.output}")
