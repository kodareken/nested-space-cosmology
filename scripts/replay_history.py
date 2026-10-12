#!/usr/bin/env python3
"""Explicit, isolated replay of preserved editions; never a current-source fallback."""
from __future__ import annotations

import argparse
import hashlib
import io
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
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
V7_COMMIT = "7fee894ad5a4da2ba378f724e45e926401f82670"
V7_CONTEXT = ROOT / ".build/replay" / V7_COMMIT[:12]
V7_BUNDLE = "results/development/nsc-cosmological-atlas-v7"


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


def install_archive(data: bytes, context: Path, *, retained_paths=None) -> dict[str, str]:
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
            retained = (preserved_artifact(member.name) if retained_paths is None
                        else member.name in retained_paths)
            if retained:
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


def prepare_v7(context: Path = V7_CONTEXT) -> dict:
    """Authenticate v7 evidence and materialize its original root source layout.

    Current scientific owners may evolve. Their frozen bytes come only from
    V7_COMMIT; retained result payloads must still match that commit in ROOT.
    No renderer, package installation or scientific calculation is selected.
    """
    manifest_name = V7_BUNDLE + "/manifest.json"
    manifest = json.loads(git("show", V7_COMMIT + ":" + manifest_name))
    source_paths = set(manifest["source_bindings"])
    bundle_paths = {V7_BUNDLE + "/" + name for name in manifest["outputs"]}
    bundle_paths.add(manifest_name)
    # Full src imports remain isolated even where the producer binds only the
    # modules it uses directly. The source-bound scripts and owners are exact.
    paths = sorted({"src", "pyproject.toml", *source_paths, *bundle_paths})
    retained = bundle_paths | {name for name in source_paths
                               if name.startswith(("results/", "paper/atlas-working-draft"))}
    data = git("archive", "--format=tar", V7_COMMIT, *paths)
    hashes = install_archive(data, context, retained_paths=retained)
    if not set(paths).difference({"src"}) <= set(hashes):
        raise ValueError("v7 historical source inventory is incomplete")
    receipt = {"schema": "NSC-EXPLICIT-REPLAY-v1", "source_commit": V7_COMMIT,
               "purpose": "Explicit v7 context; never substitutes for a current source check.",
               "files": hashes}
    path = context / "replay-context.json"
    encoded = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode()
    if path.exists() and path.read_bytes() != encoded:
        raise ValueError("historical context receipt changed")
    path.write_bytes(encoded)
    return {"context": str(context), "source_commit": V7_COMMIT, "files": len(hashes),
            "evidence_authenticated": len(retained), "scientific_replay": False}


def require_v7_profile(context: Path) -> dict:
    """Reject an incompatible renderer explicitly; never choose a fallback."""
    expected = json.loads((context / V7_BUNDLE / "manifest.json").read_bytes())["runtime"]
    actual = {"python": platform.python_version()}
    for name in expected.keys() - {"python"}:
        try:
            actual[name] = version(name)
        except PackageNotFoundError:
            actual[name] = "unavailable"
    if platform.system() != "Windows" or actual != expected:
        raise ValueError("v7 exact image replay requires Windows and its recorded profile "
                         + json.dumps(expected, sort_keys=True)
                         + "; actual platform=" + platform.system() + " runtime="
                         + json.dumps(actual, sort_keys=True)
                         + ". Use prepare-v7 for evidence authentication; no fallback was selected.")
    return {"platform": platform.system(), **actual}


def run_v7(context: Path = V7_CONTEXT) -> dict:
    """Run the untouched original v7 checker, including its exact image checks."""
    info = prepare_v7(context)
    profile = require_v7_profile(context)
    script = "scripts/derive_nsc_cosmological_atlas_v7.py"
    env = environment(context)
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(ROOT / "archive/source-history/objects")
    env["PYTHONNOUSERSITE"] = "1"
    result = subprocess.run([sys.executable, script, "--check"], cwd=context,
                            env=env, capture_output=True, text=True)
    if result.returncode:
        raise ValueError(f"v7 historical replay failed: {script}\n{result.stdout}\n{result.stderr}")
    return {**info, "scientific_replay": True, "runtime": profile,
            "checks": [{"script": script, "result": json.loads(result.stdout)}]}


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
    parser.add_argument("kind", choices=["prepare", "atlas", "finite", "companion", "tests", "all",
                                        "prepare-v7", "atlas-v7"],
                        help="v7 commands are explicit; legacy all retains its separate v1-v6 profile")
    args = parser.parse_args()
    if args.kind == "prepare-v7":
        result = prepare_v7()
    elif args.kind == "atlas-v7":
        result = run_v7()
    else:
        result = prepare() if args.kind == "prepare" else run(args.kind)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
