#!/usr/bin/env python3
"""Reproduce the read-only PREF27 GEN1 store binder at its fixed result path."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto18_pref27_binder import build_pref27  # noqa: E402

RESULT = ROOT / "results/fgc-1-pro18-pref27.json"
CONFIG = ROOT / "configs/fgc/fgc-1-pro18-pref27.toml"

def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()

def read_canonical(path: Path) -> dict[str, object]:
    def duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in pairs:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("PREF27 result is malformed") from error
    if not isinstance(value, dict) or canonical(value) != raw:
        raise ValueError("PREF27 result is noncanonical")
    return value

def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    raw_config = CONFIG.read_bytes()
    expected = build_pref27(
        ROOT, tomllib.loads(raw_config.decode("utf-8")),
        config_sha256=sha256(raw_config).hexdigest(),
        reproducer_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    if args.write:
        temporary = RESULT.with_name(f".{RESULT.name}.tmp")
        if temporary.exists():
            raise SystemExit("stale PREF27 temporary result exists")
        with temporary.open("xb") as handle:
            handle.write(canonical(expected))
        temporary.replace(RESULT)
    else:
        observed = read_canonical(RESULT)
        if observed != expected: raise SystemExit("PREF27 result differs from reconstructed evidence")
    print(json.dumps({"artifact_id": "FGC-1-PRO18-PREF27", "verified": bool(args.verify)}, sort_keys=True)); return 0

if __name__ == "__main__": raise SystemExit(main())
