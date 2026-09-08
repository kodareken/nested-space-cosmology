#!/usr/bin/env python3
"""Validate the curated repository's public boundary and claim routing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
RESULT_MANIFEST = ROOT / "results" / "manifest.json"
PROVENANCE = ROOT / "PUBLICATION-PROVENANCE.json"
PAPER_MANIFEST = ROOT / "paper" / "build-manifest.json"
IGNORED_PARTS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    ".build",
}
EXCLUDED_ROOTS = {"archive", "collected-data", "configs", "mk", "runs"}
REQUIRED = {
    "README.md",
    "THEORY.md",
    "AGENTS.md",
    "CONTRIBUTING.md",
    "CITATION.cff",
    "LICENSE",
    "SECURITY.md",
    "docs/current-result.md",
    "docs/prior-art-and-open-claim.md",
    "docs/reproducing.md",
    "docs/nsc-closure-verification-2026-09-07.md",
    "docs/nsc-regulated-recursion.md",
    "docs/nsc-radial-spectrum.md",
    "docs/nsc-geometric-chain.md",
    "paper/nested-space-cosmology.md",
    "paper/nested-space-cosmology.pdf",
    "paper/build-manifest.json",
    "results/manifest.json",
    "PUBLICATION-PROVENANCE.json",
}
ALLOWED_CATEGORIES = {
    "imported_benchmark_reproduction",
    "repository_derived_exact_identity",
    "repository_derived_numerical_result",
    "diagnostic_nonpass",
    "superseded_candidate",
    "current_frontier",
}
FORBIDDEN_LITERAL = {
    "/Users" + "/admin": "private absolute user path",
    ".co" + "dex": "private Codex state",
    ".ag" + "ents": "private agent state",
    "DouglasMac" + ".local": "private host name",
    "bh-infinity" + "@localhost": "private commit identity",
}
SECRET_PATTERNS = {
    "GitHub token": re.compile(r"gh[pousr]_[A-Za-z0-9_]{20,}"),
    "OpenAI key": re.compile(r"sk-(?:proj-)?[A-Za-z0-9_-]{20,}"),
    "AWS access key": re.compile(r"AKIA[0-9A-Z]{16}"),
    "private key": re.compile(r"BEGIN (?:RSA|OPENSSH|EC|PGP) PRIVATE KEY"),
}
LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def public_files() -> list[Path]:
    paths = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if any(
            part in IGNORED_PARTS or part.endswith(".egg-info")
            for part in relative.parts
        ):
            continue
        if path.is_symlink():
            paths.append(path)
        elif path.is_file() and path.name != ".DS_Store" and path.suffix != ".pyc":
            paths.append(path)
    return sorted(paths, key=lambda item: str(item.relative_to(ROOT)))


def check_manifest(errors: list[str]) -> None:
    try:
        manifest = json.loads(RESULT_MANIFEST.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"cannot read result manifest: {exc}")
        return
    steps = manifest.get("steps", [])
    if manifest.get("schema") != "NSC-PUBLIC-RESULT-MANIFEST-v1":
        errors.append("unexpected result manifest schema")
    if len(steps) != 64 or manifest.get("result_count") != 64:
        errors.append("result manifest does not contain exactly 64 steps")
    if manifest.get("historical_result_count") != 58:
        errors.append("historical result count must remain 58")
    if manifest.get("source_commit") != "ff2cf2722b966589b98a61accdbb6cee819a58c7":
        errors.append("result manifest has the wrong historical source commit")
    if (
        manifest.get("follow_up_source_commit")
        != "5f38712ca01ddd71e715fd265088925a73369aba"
    ):
        errors.append("result manifest has the wrong follow-up source commit")
    if manifest.get("follow_up_artifact_id") != "NSC-2-ZETA1-UNIT-CLOSURE-CHECK":
        errors.append("follow-up identity is incorrect")
    if (
        manifest.get("nsc3_source_commit")
        != "445d5b069adea8b5641384b0ee02f07fe9cfcc0a"
    ):
        errors.append("result manifest has the wrong nsc-3 source commit")
    expected_scoped = [
        "results/nsc-2-zeta1-unit-closure-check.json",
        "results/nsc-3-regulated-recursion.json",
        "results/nsc-3-radial-spectrum.json",
        "results/nsc-3-geometric-chain.json",
        "results/nsc-3-threshold-response.json",
        "results/nsc-3-boundary-response.json",
    ]
    if manifest.get("scoped_follow_up_outputs") != expected_scoped:
        errors.append("scoped follow-up outputs are incorrect")
    outputs: set[str] = set()
    current = []
    for step in steps:
        output = str(step.get("output", ""))
        generator = str(step.get("generator", ""))
        if output in outputs:
            errors.append(f"duplicate result output: {output}")
        missing = sorted(set(step.get("dependencies", [])) - outputs)
        if missing:
            errors.append(f"non-topological dependencies at {output}: {missing}")
        outputs.add(output)
        category = step.get("category")
        if category not in ALLOWED_CATEGORIES:
            errors.append(f"invalid category for {output}: {category}")
        if category == "current_frontier":
            current.append(step.get("artifact_id"))
        for relative, field in (
            (output, "output_sha256"),
            (generator, "generator_sha256"),
        ):
            path = ROOT / relative
            if not path.is_file():
                errors.append(f"manifest path is absent: {relative}")
            elif sha256(path) != step.get(field):
                errors.append(f"manifest hash mismatch: {relative}")
        for auxiliary in step.get("auxiliary_inputs", []):
            relative = str(auxiliary.get("path", ""))
            path = ROOT / relative
            if not path.is_file() or sha256(path) != auxiliary.get("sha256"):
                errors.append(f"auxiliary provenance input mismatch: {relative}")
    if current != ["NSC-2-ZETA1-RECURSION-MAP"]:
        errors.append(f"unexpected current frontier: {current}")
    if manifest.get("frontier_artifact_id") != "NSC-2-ZETA1-RECURSION-MAP":
        errors.append("frontier identity is incorrect")
    result_files = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "results").glob("nsc-*.json")
    }
    if result_files != outputs:
        errors.append("tracked NSC result files differ from the manifest closure")


def check_provenance(errors: list[str], files: list[Path]) -> None:
    try:
        value = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"cannot read publication provenance: {exc}")
        return
    if value.get("source_commit") != "ff2cf2722b966589b98a61accdbb6cee819a58c7":
        errors.append("publication provenance has the wrong source commit")
    declared = {entry["path"]: entry for entry in value.get("included_paths", [])}
    actual = {
        str(path.relative_to(ROOT))
        for path in files
        if path != PROVENANCE
    }
    if set(declared) != actual:
        errors.append("publication provenance path set differs from the checkout")
    for relative, entry in declared.items():
        path = ROOT / relative
        if path.is_file() and (
            sha256(path) != entry.get("sha256") or path.stat().st_size != entry.get("bytes")
        ):
            errors.append(f"publication provenance mismatch: {relative}")


def check_paper(errors: list[str]) -> None:
    try:
        value = json.loads(PAPER_MANIFEST.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"cannot read paper build manifest: {exc}")
        return
    for path_field, hash_field in (
        ("source", "source_sha256"),
        ("metadata", "metadata_sha256"),
        ("output", "pdf_sha256"),
    ):
        path = ROOT / str(value.get(path_field, ""))
        if not path.is_file() or sha256(path) != value.get(hash_field):
            errors.append(f"paper build manifest mismatch: {path_field}")


def check_markdown_links(errors: list[str], path: Path, text: str) -> None:
    for raw_target in LINK_RE.findall(text):
        target = raw_target.strip().strip("<>").split()[0]
        if target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target = unquote(target.split("#", 1)[0])
        if not target:
            continue
        if target.startswith("/"):
            errors.append(f"absolute Markdown link in {path.relative_to(ROOT)}: {target}")
            continue
        resolved = (path.parent / target).resolve()
        try:
            resolved.relative_to(ROOT.resolve())
        except ValueError:
            errors.append(f"Markdown link escapes repository in {path.relative_to(ROOT)}")
            continue
        if not resolved.exists():
            errors.append(f"broken Markdown link in {path.relative_to(ROOT)}: {target}")


def check_public_boundary(errors: list[str], files: list[Path]) -> None:
    relative_files = {str(path.relative_to(ROOT)) for path in files}
    for required in sorted(REQUIRED - relative_files):
        errors.append(f"required public file is absent: {required}")
    for excluded in sorted(EXCLUDED_ROOTS & {path.parts[0] for path in map(Path, relative_files)}):
        errors.append(f"excluded historical root is present: {excluded}")
    if "PLAN.md" in relative_files:
        errors.append("private PLAN.md must not be published")
    if "paper/recursive-horizons.pdf" in relative_files:
        errors.append("superseded Recursive Horizons PDF must not be published")
    for path in files:
        if path.is_symlink():
            errors.append(f"symlink is not allowed in the public core: {path.relative_to(ROOT)}")
            continue
        if path.suffix.lower() not in {".md", ".py", ".json", ".toml", ".txt", ".yml", ".yaml", ".cff", ""}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for literal, label in FORBIDDEN_LITERAL.items():
            if literal in text:
                errors.append(f"{label} appears in {path.relative_to(ROOT)}")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                errors.append(f"possible {label} appears in {path.relative_to(ROOT)}")
        if path.suffix.lower() == ".md":
            check_markdown_links(errors, path, text)


def main() -> int:
    errors: list[str] = []
    files = public_files()
    check_public_boundary(errors, files)
    check_manifest(errors)
    check_paper(errors)
    check_provenance(errors, files)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"publication check failed with {len(errors)} finding(s)", file=sys.stderr)
        return 1
    print(f"publication check passed: {len(files)} curated files, 64 result steps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
