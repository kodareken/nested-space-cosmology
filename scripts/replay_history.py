#!/usr/bin/env python3
"""Explicit, isolated replay of preserved editions; never a current-source fallback."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "04a0cb9378977772243ebf7c9be3c74182b030f2"
CONTEXT = ROOT / ".build/replay" / BASELINE[:12]
PREFIXES = ("src", "scripts", "tests", "docs", "paper", "results",
            "lab/src", "lab/scripts", "lab/tests", "lab/docs", "lab/.source-history")
BUNDLES = ["lab/results/development/nsc-cosmological-atlas" + ("-v" + str(v) if v > 1 else "-v1")
           for v in range(1, 7)]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True)
    if result.returncode:
        raise ValueError("historical Git object unavailable: " + result.stderr.decode(errors="replace"))
    return result.stdout


def safe_member(name: str, context: Path) -> Path:
    path = (context / name).resolve()
    if not path.is_relative_to(context.resolve()):
        raise ValueError("historical archive path escapes context")
    return path


def current_artifact(name: str) -> Path:
    if name.startswith("lab/results/"):
        return ROOT / name.removeprefix("lab/")
    if name == "paper/nested-space-cosmology.md" or name == "paper/nested-space-cosmology.pdf":
        return ROOT / "paper/archive/foundation-v0.26.0" / name
    return ROOT / name


def preserved_artifact(name: str) -> bool:
    return name.startswith("lab/results/") or (
        name.startswith("paper/") and not name.startswith(("paper/figures/", "paper/pending/"))
        and name not in {"paper/metadata.json", "paper/catalog.json", "paper/claim-bindings.json", "paper/build-manifest.json"})


def install_archive(data: bytes, context: Path) -> dict[str, str]:
    """Creation/resume only: changed context bytes fail instead of being repaired."""
    hashes = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:") as archive:
        for member in archive:
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError("non-file historical archive member: " + member.name)
            target = safe_member(member.name, context)
            content = archive.extractfile(member).read()
            if preserved_artifact(member.name):
                actual = current_artifact(member.name)
                if not actual.is_file() or sha(actual.read_bytes()) != sha(content):
                    raise ValueError("preserved evidence changed or missing: " + str(actual))
            if target.exists():
                if not target.is_file() or target.read_bytes() != content:
                    raise ValueError("historical context changed: " + str(target))
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
            hashes[member.name] = sha(content)
    return hashes


def prepare(context: Path = CONTEXT) -> dict:
    data = git("archive", "--format=tar", BASELINE, *PREFIXES, *BUNDLES)
    hashes = install_archive(data, context)
    receipt = {"schema": "NSC-EXPLICIT-REPLAY-v1", "source_commit": BASELINE,
               "purpose": "Historical context only; never substitutes for a current source check.",
               "files": hashes}
    path = context / "replay-context.json"
    encoded = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode()
    if path.exists() and path.read_bytes() != encoded:
        raise ValueError("historical context receipt changed")
    path.write_bytes(encoded)
    return {"context": str(context), "source_commit": BASELINE, "files": len(hashes)}


def environment(context: Path, laboratory: bool = False) -> dict[str, str]:
    env = os.environ.copy()
    owner = context / "lab" if laboratory else context
    paths = [str(owner / "src"), str(owner / "scripts")]
    # Preserve dependency lookup (some Windows venv redirectors need this),
    # but never import live project source during historical replay.
    paths.extend(p for p in env.get("PYTHONPATH", "").split(os.pathsep) if "site-packages" in p)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    env["GIT_DIR"] = str(ROOT / ".git")
    env["GIT_WORK_TREE"] = str(context)
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(context / "lab/.source-history/objects")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def run(kind: str, context: Path = CONTEXT) -> dict:
    info = prepare(context)
    commands = []
    if kind in {"atlas", "all"}:
        for version in range(1, 7):
            suffix = "" if version == 1 else "_v" + str(version)
            commands.append((True, f"scripts/derive_nsc_cosmological_atlas{suffix}.py", ["--check"]))
    if kind in {"finite", "all"}:
        commands.append((False, "scripts/verify_finite_regeneration.py", []))
    if kind in {"companion", "all"}:
        commands.append((False, "scripts/verify_local_gate_draft.py", []))
    if kind in {"tests", "all"}:
        commands.append((False, "-m", ["pytest", "-q", "tests/test_publication.py"]))
    results = []
    for lab, script, args in commands:
        owner = context / "lab" if lab else context
        result = subprocess.run([sys.executable, script, *args], cwd=owner,
                                env=environment(context, lab), capture_output=True, text=True)
        if result.returncode:
            raise ValueError(f"historical replay failed: {script}\n{result.stdout}\n{result.stderr}")
        results.append({"script": script, "result": result.stdout.strip() if script == "-m" else json.loads(result.stdout)})
    return {**info, "checks": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["prepare", "atlas", "finite", "companion", "tests", "all"])
    args = parser.parse_args()
    print(json.dumps(prepare() if args.kind == "prepare" else run(args.kind), indent=2))


if __name__ == "__main__":
    main()
