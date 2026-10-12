#!/usr/bin/env python3
"""Build or verify the compact, non-executing HLT16 MON16 record."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.hlt16_mon16_artifact import build_mon16  # noqa: E402


CONFIG = ROOT / "configs/fgc/fgc-1-hlt16-mon16.toml"
RESULT = ROOT / "results/fgc-1-hlt16-mon16.json"


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()


def build() -> dict[str, object]:
    config_raw = CONFIG.read_bytes()
    result = build_mon16(ROOT, tomllib.loads(config_raw.decode("utf-8")))
    result.update({
        "schema_version": 1,
        "project_version": "0.11.0",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "implementation_sha256": {
            "src/recursive_horizons/fgc/evolution/hlt16_mon16_artifact.py": sha256(
                (ROOT / "src/recursive_horizons/fgc/evolution/hlt16_mon16_artifact.py").read_bytes()
            ).hexdigest(),
            "scripts/reproduce_fgc_hlt16_mon16.py": sha256(Path(__file__).read_bytes()).hexdigest(),
        },
    })
    return result


def refresh_config_inventory() -> None:
    """Replace only visible pending inventory hashes with observed digests."""
    raw = CONFIG.read_text(encoding="utf-8")
    parsed = tomllib.loads(raw)
    observed = build_mon16(ROOT, parsed)["artifact_payload"]["inventory"]
    replacements = {
        item["path"]: item["observed_sha256"]
        for item in observed
        if item["binding_status"] == "pending_coordinator_refresh"
    }
    lines = raw.splitlines()
    current_path: str | None = None
    changed: set[str] = set()
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('path = "') and stripped.endswith('"'):
            current_path = stripped[8:-1]
            continue
        if (
            current_path in replacements
            and stripped == 'sha256 = "PENDING_COORDINATOR_REFRESH"'
        ):
            lines[index] = f'sha256 = "{replacements[current_path]}"'
            changed.add(current_path)
            current_path = None
    if changed != set(replacements):
        missing = sorted(set(replacements) - changed)
        raise SystemExit(f"MON16 pending inventory refresh was incomplete: {missing}")
    temporary = CONFIG.with_name(f".{CONFIG.name}.tmp")
    if temporary.exists():
        raise SystemExit("stale MON16 temporary configuration exists")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(CONFIG)


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--refresh-config", action="store_true")
    args = parser.parse_args()
    if args.refresh_config:
        refresh_config_inventory()
        print(json.dumps({"artifact_id": "FGC-1-HLT16-MON16", "config_refreshed": True}, sort_keys=True))
        return 0
    result = build()
    if args.write:
        temporary = RESULT.with_name(f".{RESULT.name}.tmp")
        if temporary.exists():
            raise SystemExit("stale MON16 temporary result exists")
        temporary.write_bytes(canonical(result))
        temporary.replace(RESULT)
    elif not RESULT.is_file() or RESULT.read_bytes() != canonical(result):
        raise SystemExit("MON16 result differs from reconstructed compact inventory")
    print(json.dumps({"artifact_id": "FGC-1-HLT16-MON16", "verified": args.verify}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
