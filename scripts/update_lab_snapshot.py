#!/usr/bin/env python3
"""Retired import updater; current successor registration has one root owner.

The original update function is retained for isolated historical fixture tests.
The active command cannot rewrite docs/lab-snapshot.json after reassembly.
Use scripts/reassemble_repository.py --record-successors with explicit current
repository-relative paths instead.
"""
import argparse
import json

from check_lab_snapshot import digest, safe_path


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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="legacy arguments; no files will be changed")
    parser.parse_args(argv)
    parser.error("retired: use python scripts/reassemble_repository.py --record-successors PATH [PATH ...]; "
                 "the original laboratory import receipt is preserved")


if __name__ == "__main__":
    main()
