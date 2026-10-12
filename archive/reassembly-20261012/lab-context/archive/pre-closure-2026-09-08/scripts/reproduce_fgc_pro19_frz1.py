#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, sys, tomllib
from hashlib import sha256
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto19_progression_freeze import build_freeze
CONFIG = ROOT / "configs/fgc/fgc-1-pro19-frz1.toml"; RESULT = ROOT / "results/fgc-1-pro19-frz1.json"
def canonical(v): return (json.dumps(v, sort_keys=True, indent=2) + "\n").encode()
def main():
 p=argparse.ArgumentParser(); g=p.add_mutually_exclusive_group(required=True); g.add_argument("--write",action="store_true"); g.add_argument("--verify",action="store_true"); a=p.parse_args()
 c=tomllib.loads(CONFIG.read_text()); payload=build_freeze(ROOT,c)
 payload.update({"schema_version":1,"project_version":"0.11.0","source_config_sha256":sha256(CONFIG.read_bytes()).hexdigest()})
 if a.write: RESULT.write_bytes(canonical(payload))
 elif RESULT.read_bytes()!=canonical(payload): raise SystemExit("PROTO19 freeze result differs")
 print(json.dumps({"artifact_id":"FGC-1-PRO19-FRZ1","verified":a.verify},sort_keys=True))
if __name__ == "__main__": main()
