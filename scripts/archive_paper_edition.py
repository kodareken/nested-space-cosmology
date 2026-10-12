"""Preserve the original foundation before the living manuscript replaces it."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "paper/archive/foundation-v0.26.0"
SOURCES = (
    "paper/nested-space-cosmology.md", "paper/nested-space-cosmology.pdf",
    "paper/metadata.json", "paper/build-manifest.json", "paper/claim-bindings.json",
    "paper/catalog.json", "scripts/build_paper.py",
)

def main():
    if DEST.exists():
        raise FileExistsError("edition already archived; authenticate it instead of overwriting")
    baseline = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    expected = {
        "paper/nested-space-cosmology.md": "509f61e7cec3e0a294bcb5f492b805c6ae0216eb1ae1b6fab0a7eb29ec329ff3",
        "paper/nested-space-cosmology.pdf": "159925ae0f785f0b734ce0a18499dde77d10a0692b805cb4615a531df2246faf",
    }
    snapshots = {name: (ROOT/name).read_bytes() for name in SOURCES}
    for name, digest in expected.items():
        if hashlib.sha256(snapshots[name]).hexdigest() != digest:
            raise ValueError(f"unexpected foundation bytes: {name}")
    files = []
    for name, data in snapshots.items():
        target = DEST/name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as handle:
            handle.write(data)
        if target.read_bytes() != data or (ROOT/name).read_bytes() != data:
            raise ValueError(f"archive verification failed: {name}")
        files.append({"original_path": name, "storage_path": target.relative_to(ROOT).as_posix(),
                      "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
    manifest = {"schema": "NSC-PAPER-ARCHIVE-v1", "id": "foundation", "version": "0.26.0",
                "source_commit": baseline, "files": files}
    (DEST/"archive-manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"archived": str(DEST), "verified_files": len(files), "source_commit": baseline}))

if __name__ == "__main__":
    main()
