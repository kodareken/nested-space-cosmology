#!/usr/bin/env python3
"""Bind a completed three-receipt comparison into the NSC-P001 register.

The comparison and its input receipts are checked before the register is
replaced. Failure leaves the existing register bytes in place. No evolution
is run and no endpoint value is invented here.

    python scripts/run_python.py scripts/register_nsc_atlas_fixed_endpoint.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

from analyze_nsc_atlas_fixed_endpoint import compare
from check_prediction_register import (
    REGISTER, assert_comparison_matches, check, check_register, endpoint_receipts,
    finite_amplitude_block, load_comparison, local_file, read_json, require,
)

ROOT = Path(__file__).resolve().parents[1]


def repo_relative(root, path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    full = candidate.resolve()
    require(full.is_file() and full.is_relative_to(root),
            "comparison must be an existing repository file")
    relative = full.relative_to(root).as_posix()
    local_file(root, relative)
    require(relative != REGISTER, "comparison must not be the prediction register")
    return relative


def endpoint_rows(found):
    rows = []
    for binding, record in found.values():
        require(type(record.get("level")) is int, "receipt level must be an integer")
        rows.append(dict(level=record["level"], path=binding["path"], bytes=binding["bytes"],
                         sha256=binding["sha256"]))
    rows.sort(key=lambda row: row["level"])
    return rows


def install_extension(entry, block):
    require(isinstance(entry.get("observational_test"), dict), "observational test missing")
    entry["observational_test"]["SI_scale_selected"] = True
    if "finite_amplitude" in entry:
        entry["finite_amplitude"] = block
        return
    rebuilt = {}
    for key, value in entry.items():
        if key == "evidence":
            rebuilt["finite_amplitude"] = block
        rebuilt[key] = value
    if "finite_amplitude" not in rebuilt:
        rebuilt["finite_amplitude"] = block
    entry.clear()
    entry.update(rebuilt)


def commit_register(path, original, text, root):
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
        return check(root)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        if path.exists() and path.read_bytes() != original:
            path.write_bytes(original)
        raise
    finally:
        if temporary.exists():
            temporary.unlink()


def register(root=ROOT, comparison=None):
    require(comparison is not None, "comparison path required")
    root = Path(root).resolve()
    register_path = local_file(root, REGISTER)
    original = register_path.read_bytes()
    current = read_json(original)
    check_register(root, current)
    relative = repo_relative(root, comparison)
    raw = local_file(root, relative).read_bytes()
    binding = dict(path=relative, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    comparison_record, found = load_comparison(root, binding)
    rows = endpoint_rows(found)
    assert_comparison_matches(comparison_record, compare(endpoint_receipts(found, rows)))
    updated = read_json(original)
    install_extension(updated["entries"][0], finite_amplitude_block(comparison_record, binding, rows))
    text = json.dumps(updated, indent=2, allow_nan=False) + "\n"
    check_register(root, read_json(text))
    confirmed = commit_register(register_path, original, text, root)
    return dict(prediction="NSC-P001", extension="nsc-p001-finite-amplitude-v1",
                prediction_status=updated["entries"][0]["finite_amplitude"]["prediction_status"],
                endpoints=3, comparison=relative, physical_evolutions=0, scientific_replays=0,
                register_check=confirmed)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    try:
        result = register(args.root, args.comparison)
    except (ValueError, OSError, KeyError, TypeError, IndexError) as error:
        print("prediction register update failed: " + str(error))
        return 1
    print("prediction register updated: " + json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
