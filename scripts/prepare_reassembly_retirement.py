#!/usr/bin/env python3
"""Preserve and authenticate old-tree retirement targets; never delete them.

Default creates a receipt and any missing recovery copies. --check rechecks the
complete resident inventory and recovery bytes; --check-recovery works after
the explicitly named trees have been retired. All Git invocations are read-only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PINNED = "d273ce6c86cacdf35d57c97292b270fac5fa6713"
TARGETS = (
    "BlackHoles-Infinity-w2", "BlackHoles-Infinity",
    ".atlas-v2-preview-20261010", ".atlas-v3-preview-20261010",
    ".atlas-v4-preview-20261010", ".atlas-v5-preview-20261010",
    ".atlas-v5-preview-b-20261010", ".atlas-v6-preview-20261010",
)
META = ".git/worktrees/BlackHoles-Infinity-w2"
BASE = Path("archive/reassembly-20261012/retirement-additional")
RECEIPT = Path("results/reassembly/retirement.json")
MIGRATION = Path("results/reassembly/manifest.json")


def confined(relative: str | Path) -> Path:
    candidate = ROOT / relative
    resolved = candidate.resolve()
    if not resolved.is_relative_to(ROOT.resolve()) or resolved == ROOT.resolve():
        raise ValueError(f"Path leaves named repository: {relative}")
    # Reject intermediate indirection as well as final symlinks/junctions.
    part = ROOT
    for component in candidate.relative_to(ROOT).parts:
        part /= component
        if part.is_symlink() or part.is_junction():
            raise ValueError(f"Filesystem indirection: {part}")
    return resolved


def digest(path: Path) -> tuple[int, str]:
    before = path.stat()
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"File changed during authentication: {path}")
    return after.st_size, sha


def inventory(relative: str) -> tuple[list[Path], list[str]]:
    start = confined(relative)
    if not start.is_dir():
        raise ValueError(f"Missing named target directory: {relative}")
    files, directories = [], [relative]
    for directory, names, filenames in os.walk(start, followlinks=False):
        names.sort()
        filenames.sort()
        for name in names:
            child = Path(directory) / name
            confined(child.relative_to(ROOT))
            directories.append(child.relative_to(ROOT).as_posix())
        for name in filenames:
            child = confined((Path(directory) / name).relative_to(ROOT))
            if not child.is_file():
                raise ValueError(f"Not a resident regular file: {child}")
            files.append(child)
    return files, directories


def git(*args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True)
    if result.returncode:
        raise ValueError(result.stderr.decode("utf-8", errors="replace"))
    return result.stdout


def commit_tree() -> dict[str, str]:
    actual = git("rev-parse", "--verify", PINNED + "^{commit}").decode().strip()
    if actual != PINNED:
        raise ValueError("Pinned commit mismatch")
    result = {}
    for entry in git("ls-tree", "-r", "-z", PINNED).split(b"\0"):
        if not entry:
            continue
        meta, path = entry.split(b"\t", 1)
        mode, kind, oid = meta.split()
        if kind == b"blob" and mode in (b"100644", b"100755"):
            result[path.decode("utf-8")] = oid.decode("ascii")
    return result


class BlobReader:
    def __init__(self) -> None:
        self.process = subprocess.Popen(
            ["git", "cat-file", "--batch"], cwd=ROOT,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.cache: dict[str, tuple[int, str]] = {}

    def digest(self, oid: str) -> tuple[int, str]:
        if oid in self.cache:
            return self.cache[oid]
        proc = self.process
        assert proc.stdin is not None and proc.stdout is not None
        proc.stdin.write(oid.encode("ascii") + b"\n")
        proc.stdin.flush()
        header = proc.stdout.readline().split()
        if len(header) != 3 or header[0].decode() != oid or header[1] != b"blob":
            raise ValueError(f"Missing/non-blob Git object: {oid}")
        size = int(header[2])
        remaining, sha = size, hashlib.sha256()
        while remaining:
            chunk = proc.stdout.read(min(1024 * 1024, remaining))
            if not chunk:
                raise ValueError("Truncated Git blob stream")
            sha.update(chunk)
            remaining -= len(chunk)
        if proc.stdout.read(1) != b"\n":
            raise ValueError("Malformed Git batch boundary")
        self.cache[oid] = (size, sha.hexdigest())
        return self.cache[oid]

    def close(self) -> None:
        assert self.process.stdin is not None
        self.process.stdin.close()
        code = self.process.wait()
        if code:
            raise ValueError("Git batch reader failed")


def write_zip(relative: Path, members: dict[str, Path]) -> dict:
    """Create once; on a repeat authenticate all members without overwriting."""
    target = confined(relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = {}
    for member, source in sorted(members.items()):
        size, sha = digest(source)
        expected[member] = {"bytes": size, "sha256": sha,
                            "source_mtime_ns": source.stat().st_mtime_ns,
                            "source_mode": source.stat().st_mode}
    if not target.exists():
        with zipfile.ZipFile(target, "x", zipfile.ZIP_DEFLATED, compresslevel=6) as output:
            for member, source in sorted(members.items()):
                info = zipfile.ZipInfo(member, date_time=(1980, 1, 1, 0, 0, 0))
                info.create_system = 3
                info.external_attr = (expected[member]["source_mode"] & 0xFFFF) << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                with source.open("rb") as incoming, output.open(info, "w", force_zip64=True) as outgoing:
                    shutil.copyfileobj(incoming, outgoing, length=1024*1024)
    with zipfile.ZipFile(target) as packed:
        if len(packed.namelist()) != len(expected) or set(packed.namelist()) != set(expected):
            raise ValueError(f"Archive inventory drift: {relative}")
        for member, entry in expected.items():
            with packed.open(member) as stream:
                sha = hashlib.file_digest(stream, "sha256").hexdigest()
            if packed.getinfo(member).file_size != entry["bytes"] or sha != entry["sha256"]:
                raise ValueError(f"Archive content drift: {relative}:{member}")
    size, sha = digest(target)
    return {"path": relative.as_posix(), "bytes": size, "sha256": sha,
            "members": expected}


class RecoveryVerifier:
    def __init__(self, blobs: BlobReader) -> None:
        self.blobs = blobs
        self.files: dict[str, tuple[int, str]] = {}
        self.archives: dict[str, zipfile.ZipFile] = {}
        self.members: dict[tuple[str, str], tuple[int, str]] = {}

    def file(self, relative: str) -> tuple[int, str]:
        if relative not in self.files:
            self.files[relative] = digest(confined(relative))
        return self.files[relative]

    def member(self, archive: str, member: str) -> tuple[int, str]:
        if Path(member).is_absolute() or ".." in Path(member).parts:
            raise ValueError("Unsafe archive member")
        key = (archive, member)
        if key not in self.members:
            if archive not in self.archives:
                self.archives[archive] = zipfile.ZipFile(confined(archive))
                names = self.archives[archive].namelist()
                if len(names) != len(set(names)):
                    raise ValueError("Duplicate archive member names")
            packed = self.archives[archive]
            with packed.open(member) as stream:
                sha = hashlib.file_digest(stream, "sha256").hexdigest()
            self.members[key] = (packed.getinfo(member).file_size, sha)
        return self.members[key]

    def target(self, recovery: dict) -> tuple[int, str]:
        if recovery["kind"] == "git_blob":
            return self.blobs.digest(recovery["blob"])
        if recovery["kind"] == "file":
            return self.file(recovery["path"])
        if recovery["kind"] == "zip_member":
            return self.member(recovery["archive"], recovery["member"])
        raise ValueError("Unknown recovery kind")

    def close(self) -> None:
        for archive in self.archives.values():
            archive.close()


def prepare() -> dict:
    destination = confined(RECEIPT)
    if destination.exists():
        raise ValueError("Receipt already exists; use --check rather than replace its inventory.")
    migration = json.loads(confined(MIGRATION).read_text(encoding="utf-8"))
    candidates = defaultdict(list)
    for row in migration["rows"]:
        candidates[row["source"]].append(row)
    tree, blobs = commit_tree(), BlobReader()
    verifier = RecoveryVerifier(blobs)
    rows, directory_rows, additional = [], [], []
    try:
        metadata_files, metadata_directories = inventory(META)
        pointer = confined("BlackHoles-Infinity-w2/.git")
        if not pointer.is_file():
            raise ValueError("Worktree .git pointer is absent")
        metadata_members = {p.relative_to(ROOT).as_posix(): p for p in metadata_files + [pointer]}
        metadata = write_zip(BASE / "worktree-metadata.zip", metadata_members)
        archive_receipts = {metadata["path"]: metadata}
        for target in TARGETS:
            files, directories = inventory(target)
            directory_rows.append({
                "path": target, "absolute_path": str(confined(target)),
                "files": len(files), "directories": directories,
            })
            for source in files:
                relative = source.relative_to(ROOT).as_posix()
                size, sha = digest(source)
                expected = (size, sha)
                recovery = None
                if relative == "BlackHoles-Infinity-w2/.git":
                    recovery = {"kind": "zip_member", "archive": metadata["path"], "member": relative}
                logical = source.relative_to(confined(target)).as_posix()
                if recovery is None and target == "BlackHoles-Infinity-w2" and logical in tree:
                    oid = tree[logical]
                    if blobs.digest(oid) == expected:
                        recovery = {"kind": "git_blob", "commit": PINNED,
                                    "path": logical, "blob": oid}
                if recovery is None:
                    for previous in candidates[relative]:
                        candidate = confined(previous["destination"])
                        if candidate.is_file() and verifier.file(previous["destination"]) == expected:
                            recovery = {"kind": "file", "path": previous["destination"]}
                            break
                        archive, member = previous.get("storage_archive"), previous.get("archive_member")
                        if archive and member and verifier.member(archive, member) == expected:
                            recovery = {"kind": "zip_member", "archive": archive, "member": member}
                            break
                if recovery is None:
                    options = []
                    if target == "BlackHoles-Infinity-w2":
                        options.append(Path("archive/reassembly-20261012/w2/files") / logical)
                        canonical = Path(logical)
                        if canonical.parts[0] == "lab":
                            canonical = Path(*canonical.parts[1:])
                        options.append(canonical)
                    elif target == "BlackHoles-Infinity":
                        canonical = Path(logical)
                        if canonical.parts[0] == "lab":
                            canonical = Path(*canonical.parts[1:])
                        options.append(canonical)
                    for option in options:
                        path = confined(option)
                        if path.is_file() and verifier.file(option.as_posix()) == expected:
                            recovery = {"kind": "file", "path": option.as_posix()}
                            break
                if recovery is None:
                    retained = BASE / "files" / relative
                    retained_path = confined(retained)
                    retained_path.parent.mkdir(parents=True, exist_ok=True)
                    if not retained_path.exists():
                        # Exclusive creation prevents accidental overwriting during retries.
                        with source.open("rb") as incoming, retained_path.open("xb") as outgoing:
                            shutil.copyfileobj(incoming, outgoing, length=1024*1024)
                        shutil.copystat(source, retained_path)
                    if digest(retained_path) != expected:
                        raise ValueError(f"Additional preservation mismatch: {relative}")
                    recovery = {"kind": "file", "path": retained.as_posix()}
                    additional.append(retained.as_posix())
                if verifier.target(recovery) != expected:
                    raise ValueError(f"Recovery target mismatch: {relative}")
                rows.append({"source": relative, "bytes": size, "sha256": sha,
                             "recovery": recovery})
                if len(rows) % 1000 == 0:
                    print(json.dumps({"authenticated_files": len(rows), "additional": len(additional)}), flush=True)
        # Pack every ignored file selected as a recovery target. Its copied bytes
        # must remain recoverable even in a clean clone that omits ignored paths.
        file_targets = sorted({r["recovery"]["path"] for r in rows if r["recovery"]["kind"] == "file"})
        ignored_run = subprocess.run(
            ["git", "check-ignore", "--no-index", "--stdin", "-z"],
            input=("\0".join(file_targets) + "\0").encode(), cwd=ROOT, capture_output=True)
        if ignored_run.returncode not in (0, 1):
            raise ValueError(ignored_run.stderr.decode(errors="replace"))
        ignored = {p.decode("utf-8") for p in ignored_run.stdout.split(b"\0") if p}
        if ignored:
            packed = write_zip(BASE / "ignored-recovery.zip", {p: confined(p) for p in sorted(ignored)})
            archive_receipts[packed["path"]] = packed
            for row in rows:
                recovery = row["recovery"]
                if recovery["kind"] == "file" and recovery["path"] in ignored:
                    row["recovery"] = {"kind": "zip_member", "archive": packed["path"],
                                       "member": recovery["path"]}
        for entry in migration.get("storage_archives", []):
            if any(r["recovery"].get("archive") == entry["path"] for r in rows):
                size, sha = digest(confined(entry["path"]))
                if (size, sha) != (entry["bytes"], entry["sha256"]):
                    raise ValueError("Migration storage archive drift")
                archive_receipts[entry["path"]] = entry
        receipt = {
            "schema": "NSC-RETIREMENT-PRESERVATION-v1",
            "action": "authentication and preservation only; no deletion or Git metadata writes",
            "canonical_root": str(ROOT.resolve()),
            "pinned_worktree_commit": PINNED,
            "commit_references": git("for-each-ref", "--contains", PINNED, "--format=%(refname)").decode().splitlines(),
            "producer": "scripts/prepare_reassembly_retirement.py",
            "producer_sha256": digest(Path(__file__))[1],
            "migration_manifest_at_preparation": {
                "path": MIGRATION.as_posix(), "bytes": digest(confined(MIGRATION))[0],
                "sha256": digest(confined(MIGRATION))[1],
            },
            "targets": directory_rows,
            "metadata": {"original_directory": META, "directories": metadata_directories,
                         "archive": metadata["path"], "pointer_member": "BlackHoles-Infinity-w2/.git"},
            "storage_archives": sorted(archive_receipts.values(), key=lambda x: x["path"]),
            "additional_copied_files": additional,
            "rows": rows,
            "counts": {
                "resident_files": len(rows), "resident_bytes": sum(r["bytes"] for r in rows),
                "recovery_kinds": dict(Counter(r["recovery"]["kind"] for r in rows)),
                "additional_copied_files": len(additional),
                "metadata_files_including_pointer": len(metadata_members),
                "empty_retirement_targets": sum(t["files"] == 0 for t in directory_rows),
            },
            "restoration": [
                "For git_blob rows, read git cat-file blob BLOB; verify the recorded byte count and SHA256.",
                "For file rows, copy the recorded canonical-relative recovery path.",
                "For zip_member rows, extract only the named member to its recorded source path.",
                "The worktree metadata archive preserves original bytes, including stale gitdir pointers; repair path registration separately if restoring.",
                "Run --check immediately before recycling; run --check-recovery afterwards. Neither mode deletes anything.",
            ],
        }
        check(receipt, source_inventory=True)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write("\n")
        return receipt["counts"]
    finally:
        verifier.close()
        blobs.close()


def check(receipt: dict, *, source_inventory: bool) -> dict:
    tree, blobs = commit_tree(), BlobReader()
    verifier = RecoveryVerifier(blobs)
    try:
        if receipt["pinned_worktree_commit"] != PINNED:
            raise ValueError("Wrong pinned recovery commit")
        if [t["path"] for t in receipt["targets"]] != list(TARGETS):
            raise ValueError("Changed retirement target set")
        if len({r["source"] for r in receipt["rows"]}) != len(receipt["rows"]):
            raise ValueError("Duplicate resident inventory")
        for archive in receipt["storage_archives"]:
            if digest(confined(archive["path"])) != (archive["bytes"], archive["sha256"]):
                raise ValueError("Changed recovery archive")
            if isinstance(archive.get("members"), dict):
                with zipfile.ZipFile(confined(archive["path"])) as packed:
                    names = packed.namelist()
                    if len(names) != len(set(names)) or set(names) != set(archive["members"]):
                        raise ValueError("Recovery archive inventory mismatch")
                for member, entry in archive["members"].items():
                    if verifier.member(archive["path"], member) != (entry["bytes"], entry["sha256"]):
                        raise ValueError("Recovery archive member mismatch")
        if source_inventory:
            seen = []
            for target in receipt["targets"]:
                files, directories = inventory(target["path"])
                if directories != target["directories"] or len(files) != target["files"]:
                    raise ValueError(f"Source directory inventory drift: {target['path']}")
                seen.extend(p.relative_to(ROOT).as_posix() for p in files)
            if set(seen) != {r["source"] for r in receipt["rows"]}:
                raise ValueError("Source resident file inventory drift")
            metadata_files, metadata_directories = inventory(META)
            if metadata_directories != receipt["metadata"]["directories"]:
                raise ValueError("Worktree metadata directory drift")
            meta_receipt = next(a for a in receipt["storage_archives"]
                                if a["path"] == receipt["metadata"]["archive"])
            expected_meta = {m for m in meta_receipt["members"] if m.startswith(META + "/")}
            if {p.relative_to(ROOT).as_posix() for p in metadata_files} != expected_meta:
                raise ValueError("Worktree metadata resident inventory drift")
            for path in metadata_files:
                entry = meta_receipt["members"][path.relative_to(ROOT).as_posix()]
                if digest(path) != (entry["bytes"], entry["sha256"]):
                    raise ValueError("Worktree metadata byte drift")
        for row in receipt["rows"]:
            expected = (row["bytes"], row["sha256"])
            recovery = row["recovery"]
            if recovery["kind"] == "git_blob":
                if recovery["commit"] != PINNED or tree.get(recovery["path"]) != recovery["blob"]:
                    raise ValueError("Git recovery path does not name the pinned tree blob")
            if verifier.target(recovery) != expected:
                raise ValueError(f"Recovery bytes mismatch: {row['source']}")
            if source_inventory and digest(confined(row["source"])) != expected:
                raise ValueError(f"Resident source drift: {row['source']}")
        return {"authenticated_resident_rows": len(receipt["rows"]),
                "source_inventory_checked": source_inventory,
                "metadata_archives_authenticated": len(receipt["storage_archives"])}
    finally:
        verifier.close()
        blobs.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--check-recovery", action="store_true")
    args = parser.parse_args()
    if args.check or args.check_recovery:
        record = json.loads(confined(RECEIPT).read_text(encoding="utf-8"))
        result = check(record, source_inventory=args.check)
    else:
        result = prepare()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(json.dumps({"error": type(error).__name__, "detail": str(error)}), flush=True)
        raise
