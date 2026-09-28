#!/usr/bin/env python3
"""Explicitly register intended lab edits while preserving import provenance."""
import argparse
import json
import subprocess

from check_lab_snapshot import ROOT, digest, safe_path


def update(root, names, parent_commit):
    target = root / "docs/lab-snapshot.json"
    manifest = json.loads(target.read_text())
    entries = {row["path"]: row for row in manifest["files"]}
    omitted = {"lab/" + path for path in manifest["intentionally_omitted"]}
    for name in names:
        if (not name.startswith("lab/") or name.startswith("lab/.source-history/")
                or name in omitted):
            raise ValueError("only explicit current lab files may be registered: " + name)
        path = safe_path(root, name)
        if not path.is_file():
            raise ValueError("a current file is required: " + name)
        previous = entries.get(name)
        checksum = digest(path)
        if previous is not None and previous["sha256"] == checksum:
            continue
        row = dict(previous) if previous else {"path": name}
        row["bytes"] = path.stat().st_size
        row["sha256"] = checksum
        row["change"] = "explicit successor after consolidation; original import fields preserved"
        row["changed_from_commit"] = parent_commit
        entries[name] = row
    manifest["files"] = [entries[name] for name in sorted(entries)]
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="explicit repository-relative lab file paths")
    args = parser.parse_args()
    parent = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    update(ROOT, args.paths, parent)
    print("Updated explicit lab file bindings. Rebuild publication provenance and run relevant checks.")


if __name__ == "__main__":
    main()
