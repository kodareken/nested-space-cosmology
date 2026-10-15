#!/usr/bin/env python3
"""Validate the curated repository's public boundary and claim routing."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from types import ModuleType, SimpleNamespace
from pathlib import Path
import re
import subprocess
import sys
import tomllib
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
    # Match the publication producer: agent routing is not scientific evidence.
    ".pi",
    # Machine-local agent state. Keep in step with
    # build_publication_provenance.IGNORED_PARTS, or the two disagree on
    # the path set and the provenance record can never validate.
    "node_modules",
    ".co" + "dex",
    "subagents",
    ".claude",
    ".cursor",
    ".research-reset-preservation",
    "BlackHoles-Infinity",
    "BlackHoles-Infinity-w2",
}
#: Host-local agent configuration files, skipped by name.
IGNORED_NAMES = {
    "settings.json",
}
EXCLUDED_ROOTS = {"collected-data", "configs", "mk", "runs"}
REQUIRED = {
    'AGENTS.md',
    'CITATION.cff',
    'CONTRIBUTING.md',
    'LICENSE',
    'PUBLICATION-PROVENANCE.json',
    'README.md',
    'SECURITY.md',
    'THEORY.md',
    'docs/current-result.md',
    'docs/instructions.md',
    'docs/nsc-radial-spectrum.md',
    'docs/nsc-regulated-recursion.md',
    'docs/papers.md',
    'docs/prior-art-and-open-claim.md',
    'docs/reproducing.md',
    'paper/build-manifest.json',
    'paper/catalog.json',
    'paper/local-gate-draft-manifest.json',
    'paper/local-incoming-gate-draft.pdf',
    'paper/nested-space-cosmology.md',
    'paper/nested-space-cosmology.pdf',
    'results/manifest.json',
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
GITHUB_MARKDOWN_ENTRYPOINTS = {
    "README.md",
    "THEORY.md",
    "AGENTS.md",
    "CONTRIBUTING.md",
    "paper/nested-space-cosmology.md",
    "docs/current-result.md",
    "docs/development-update-2026-09-10.md",
    "docs/nsc-background-projection.md",
    "docs/nsc-child-metric-backreaction.md",
    "docs/nsc-constraint-complete-neck.md",
    "docs/nsc-coupled-ctp-metric-evolution.md",
    "docs/nsc-mode-resolved-cauchy-state.md",
    "docs/nsc-landau-cauchy-isometry.md",
    "docs/nsc-joint-history-state-bvp.md",
    "docs/nsc-general-ks-same-action-history.md",
    "docs/nsc-general-ks-reference.md",
    "docs/nsc-general-ks-local-history.md",
    "docs/nsc-tilted-landau-interface.md",
    "docs/nsc-extended-tilted-history-gate.md",
    "docs/nsc-transmitting-boundary-selection.md",
    "docs/nsc-weyl-endpoint-match.md",
    "docs/nsc-full-extended-history-gate.md",
    "docs/nsc-transmitting-ctp-variation.md",
    "docs/nsc-transmitting-metric-pullback.md",
    "docs/nsc-extended-action-completion-gate.md",
    "docs/nsc-dirac-endpoint-jets.md",
    "docs/nsc-transmitting-boundary-remainder.md",
    "docs/nsc-physical-jet-extended-gate.md",
    "docs/nsc-transmitting-dirac-domain.md",
    "docs/nsc-hamiltonian-trace-representation.md",
    "docs/nsc-common-time-bulk-split.md",
    "docs/nsc-transmitting-cross-resolvent.md",
    "docs/nsc-transmitting-ctp-resolvent.md",
    "docs/nsc-pg-lll-preparation.md",
    "docs/nsc-massive-signed-preparation.md",
    "docs/nsc-horizon-paired-pg-map.md",
    "docs/nsc-compact-matched-restart.md",
    "results/README.md",
}
UNSUPPORTED_GITHUB_MATH = {
    "operatorname macro": re.compile(r"\\operatorname\*?\{"),
    "parenthesized inline delimiter": re.compile(r"\\\(|\\\)"),
    "bracket display delimiter": re.compile(r"(?m)^\\\[$|^\\\]$"),
}


ARCHIVE_DOI = "10.5281/zenodo.23090015"


def retained_receipts() -> tuple[dict | None, dict[str, dict]]:
    path = ROOT / "docs/lab-snapshot.json"
    if not path.is_file():
        return None, {}
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    if snapshot.get("schema") != "NSC-CONSOLIDATED-LAB-v2":
        return None, {}
    if (snapshot.get("canonical_root") != "."
            or snapshot.get("provenance", {}).get("archived_import_and_reassembly_doi") != ARCHIVE_DOI):
        raise ValueError("invalid retained-file receipt identity")
    receipts = {}
    for row in snapshot.get("files", []):
        name = row["path"]
        _safe_relative(name)
        if (name in receipts or not re.fullmatch("[0-9a-f]{64}", row.get("sha256", ""))
                or not isinstance(row.get("bytes"), int) or row["bytes"] < 0):
            raise ValueError("invalid or repeated retained-file receipt: " + name)
        receipts[name] = row
    return snapshot, receipts


def _safe_relative(relative: str) -> Path:
    path = Path(relative)
    if (path.is_absolute() or path.drive or ".." in path.parts or not path.parts
            or not (ROOT / path).resolve().is_relative_to(ROOT.resolve())):
        raise ValueError("unsafe publication evidence path: " + str(relative))
    return path


def _git_bytes(commit: str, relative: str) -> bytes:
    """Read one exact pinned object without materializing a historical checkout."""
    _safe_relative(relative)
    if not re.fullmatch("[0-9a-f]{40}", commit):
        raise ValueError("historical source is not pinned")
    key = (str(ROOT), commit, relative)
    if key not in _PINNED_BYTES:
        environment = os.environ.copy()
        objects = ROOT / "archive/source-history/objects"
        if objects.is_dir():
            inherited = environment.get("GIT_ALTERNATE_OBJECT_DIRECTORIES", "")
            environment["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(objects) + (os.pathsep + inherited if inherited else "")
        result = subprocess.run(["git", "show", f"{commit}:{relative}"], cwd=ROOT,
                                env=environment, capture_output=True)
        if result.returncode:
            raise FileNotFoundError("pinned historical object unavailable: " + relative)
        _PINNED_BYTES[key] = result.stdout
    return _PINNED_BYTES[key]


_PINNED_BYTES = {}


class _PinnedPath:
    """The frozen checker's read-only path interface over one explicit Git context."""
    def __init__(self, commit, relative=""):
        self.commit, self.relative = commit, relative

    def __truediv__(self, relative):
        return type(self)(self.commit, (Path(self.relative) / relative).as_posix())

    def read_bytes(self):
        return _git_bytes(self.commit, self.relative)

    def read_text(self, encoding="utf-8"):
        return self.read_bytes().decode(encoding)

    def is_file(self):
        try:
            self.read_bytes()
            return True
        except FileNotFoundError:
            return False

    def stat(self):
        return SimpleNamespace(st_size=len(self.read_bytes()))

    def __str__(self):
        return f"pinned:{self.commit}:{self.relative}"


def historical_checker(commit):
    source = _git_bytes(commit, "scripts/reproduce_public_results.py")
    module = ModuleType("nsc_pinned_release_checker")
    module.__file__ = str(ROOT / "scripts/reproduce_public_results.py")
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    module.ROOT = _PinnedPath(commit)
    return module


def _baseline() -> str:
    snapshot, _ = retained_receipts()
    if snapshot is None:
        from replay_history import BASELINE
        return BASELINE
    return snapshot["provenance"]["historical_replay_baseline"]


def authenticate_frozen_path(relative, expected, *, result=False, context=None):
    """Local frozen outputs remain strict; disposed outputs need DOI/hash receipts."""
    path = ROOT / _safe_relative(relative)
    if path.is_file() and sha256(path) == expected:
        return "local"
    if result:
        if path.is_file():
            raise ValueError("retained frozen evidence changed: " + relative)
        value = json.loads(PROVENANCE.read_text(encoding="utf-8"))
        rows = {row["path"]: row for row in value.get("evidence_locations", [])}
        row = rows.get(relative, {})
        if (row.get("status") != "archived" or row.get("sha256") != expected
                or row.get("doi") != ARCHIVE_DOI):
            raise ValueError("missing frozen evidence has no archival receipt: " + relative)
    raw = _git_bytes(context or _baseline(), relative)
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("pinned historical evidence differs: " + relative)
    return "archived" if result else "pinned_source"


def check_development_snapshot() -> set[str]:
    """Authenticate the later snapshot without extending the frozen release."""
    snapshot_receipt, _ = retained_receipts()
    if snapshot_receipt is not None:
        return check_archived_development_snapshot()
    from reproduce_public_results import validate_authenticated_inputs
    snapshot = json.loads((ROOT / "results/development-snapshot.json").read_text())
    if snapshot.get("schema") != "NSC-DEVELOPMENT-SNAPSHOT-v1":
        raise ValueError("unexpected development snapshot schema")
    if not re.fullmatch(r"[0-9a-f]{40}", snapshot.get("source_commit", "")):
        raise ValueError("development source is not pinned")
    if snapshot["release_manifest_sha256"] != sha256(RESULT_MANIFEST):
        raise ValueError("development snapshot's preprint baseline changed")
    sources = {}
    for entry in snapshot["source_files"]:
        relative = entry["path"]
        path = Path(relative)
        if (path.is_absolute() or ".." in path.parts or relative in sources
                or not path.parts or path.parts[0] not in {"docs", "src", "scripts", "tests", "results"}):
            raise ValueError("invalid development source path")
        target = ROOT / path
        if target.stat().st_size != entry["bytes"] or sha256(target) != entry["sha256"]:
            raise ValueError(f"development source hash mismatch: {relative}")
        sources[relative] = entry
    outputs = set()
    payloads = set()
    for record in snapshot["records"]:
        output = record["output"]
        if (output in outputs or output not in sources
                or not output.startswith("results/development/")
                or record["generator"] not in sources or record["arguments"] != ["--check"]):
            raise ValueError("invalid development record declaration")
        value = json.loads((ROOT / output).read_text())
        if record["comparison"] != value["comparison"]:
            raise ValueError(f"development comparison policy changed: {output}")
        validate_authenticated_inputs(ROOT, value)
        declared_payloads = [value.get("payload")]
        auxiliary = value.get("auxiliary_payloads", [])
        if not isinstance(auxiliary, list) or any(not isinstance(p, dict) or not p.get("path") for p in auxiliary):
            raise ValueError("invalid auxiliary payload declaration")
        declared_payloads.extend(auxiliary)
        if value.get("schema") == "NSC-COMPACT-MATCHED-RESTART-v1":
            for key in ("scattering_payload", "seed_generations_payload"):
                if not isinstance(value.get(key), dict) or not value[key].get("path"):
                    raise ValueError(f"missing versioned restart payload: {key}")
                declared_payloads.append(value[key])
        for payload in declared_payloads:
            if isinstance(payload, dict) and payload.get("path"):
                relative = payload["path"]
                if (relative in payloads or relative not in sources
                        or not relative.startswith("results/development/artifacts/")):
                    raise ValueError("invalid development payload declaration")
                target = ROOT / relative
                if target.stat().st_size != payload.get("bytes") or sha256(target) != payload.get("sha256"):
                    raise ValueError(f"development payload mismatch: {relative}")
                payloads.add(relative)
        outputs.add(output)
    if outputs | payloads != {path for path in sources if path.startswith("results/")}:
        raise ValueError("development record index differs from imported results")
    return outputs


def check_archived_development_snapshot() -> set[str]:
    baseline = _baseline()
    path = ROOT / "results/development-snapshot.json"
    if path.read_bytes() != _git_bytes(baseline, "results/development-snapshot.json"):
        raise ValueError("immutable development snapshot differs from its pinned context")
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    if snapshot.get("schema") != "NSC-DEVELOPMENT-SNAPSHOT-v1":
        raise ValueError("unexpected development snapshot schema")
    if snapshot["release_manifest_sha256"] != sha256(RESULT_MANIFEST):
        raise ValueError("development release baseline changed")
    legacy = historical_checker(baseline)
    sources = {}
    for row in snapshot["source_files"]:
        relative = row["path"]
        _safe_relative(relative)
        if relative in sources or Path(relative).parts[0] not in {"docs", "src", "scripts", "tests", "results"}:
            raise ValueError("invalid development source path")
        raw = _git_bytes(baseline, relative)
        if len(raw) != row["bytes"] or hashlib.sha256(raw).hexdigest() != row["sha256"]:
            raise ValueError("pinned development source differs: " + relative)
        if relative.startswith("results/"):
            authenticate_frozen_path(relative, row["sha256"], result=True, context=baseline)
        sources[relative] = row
    outputs, payloads, local = set(), set(), set()
    for record in snapshot["records"]:
        output = record["output"]
        if (output in outputs or output not in sources or not output.startswith("results/development/")
                or record["generator"] not in sources or record["arguments"] != ["--check"]):
            raise ValueError("invalid development record declaration")
        value = json.loads(_git_bytes(baseline, output))
        if record["comparison"] != value["comparison"]:
            raise ValueError("development comparison policy changed: " + output)
        legacy.validate_authenticated_inputs(legacy.ROOT, value)
        declared = [value.get("payload"), *value.get("auxiliary_payloads", [])]
        if value.get("schema") == "NSC-COMPACT-MATCHED-RESTART-v1":
            declared += [value.get("scattering_payload"), value.get("seed_generations_payload")]
        for payload in declared:
            if isinstance(payload, dict) and payload.get("path"):
                relative = payload["path"]
                if relative in payloads or relative not in sources or not relative.startswith("results/development/artifacts/"):
                    raise ValueError("invalid development payload declaration")
                raw = _git_bytes(baseline, relative)
                if len(raw) != payload.get("bytes") or hashlib.sha256(raw).hexdigest() != payload.get("sha256"):
                    raise ValueError("pinned development payload differs: " + relative)
                payloads.add(relative)
        outputs.add(output)
        if (ROOT / output).is_file():
            local.add(output)
    if outputs | payloads != {p for p in sources if p.startswith("results/")}:
        raise ValueError("development record closure differs")
    return local


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def public_files() -> list[Path]:
    paths = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        # Historical source contexts have their own authenticated receipts.
        if relative.parts[0] in {"archive", "lab"}:
            continue
        if any(
            part in IGNORED_PARTS or part.endswith(".egg-info")
            for part in relative.parts
        ):
            continue
        if path.is_symlink():
            paths.append(path)
        elif path.is_file() and path.name != ".DS_Store" and path.suffix != ".pyc":
            if path.name in IGNORED_NAMES and relative.parts[0] in {".pi", ".co" + "dex"}:
                continue
            paths.append(path)
    return sorted(paths, key=lambda item: item.relative_to(ROOT).as_posix())


def check_manifest(errors: list[str]) -> None:
    try:
        manifest = json.loads(RESULT_MANIFEST.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"cannot read result manifest: {exc}")
        return
    steps = manifest.get("steps", [])
    if manifest.get("schema") != "NSC-PUBLIC-RESULT-MANIFEST-v1":
        errors.append("unexpected result manifest schema")
    release = json.loads((ROOT / "results/release-spec.json").read_text(encoding="utf-8"))
    expected_count = 58 + len(release["scoped_follow_ups"])
    if len(steps) != expected_count or manifest.get("result_count") != expected_count:
        errors.append("result manifest count differs from release specification")
    try:
        validate_historical_release(manifest)
    except Exception as exc:
        errors.append(f"result closure validation failed: {exc}")
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
    expected_scoped = {row["output"] for row in release["scoped_follow_ups"]}
    if set(manifest.get("scoped_follow_up_outputs", [])) != expected_scoped:
        errors.append("scoped follow-up outputs differ from release specification")
    paper = (ROOT / "paper/nested-space-cosmology.md").read_text(encoding="utf-8")
    targets = release.get("current_targets", [])
    required_targets = {"target-constants", "target-continuation", "target-dark-sector",
                        "target-measurement", "target-antimatter", "target-recursion"}
    if {target.get("claim_id") for target in targets} != required_targets or len(targets) != 6:
        errors.append("release must declare the six physical targets")
    artifact_ids = {step.get("artifact_id") for step in steps}
    for target in targets:
        if target.get("status") != "open_target" or not target.get("required_physical_result"):
            errors.append(f"invalid physical target status: {target.get('claim_id')}")
        if not target.get("supporting_artifact_ids") or not set(target["supporting_artifact_ids"]) <= artifact_ids:
            errors.append(f"physical target has absent evidence: {target.get('claim_id')}")
        if f"<!-- nsc-claim:{target['claim_id']} -->" not in paper:
            errors.append(f"paper lacks target claim anchor: {target['claim_id']}")
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
            try:
                authenticate_frozen_path(relative, step.get(field), result=(relative == output))
            except (OSError, ValueError, KeyError, TypeError) as error:
                errors.append("manifest evidence: " + str(error))
        for auxiliary in step.get("auxiliary_inputs", []):
            relative = str(auxiliary.get("path", ""))
            try:
                authenticate_frozen_path(relative, auxiliary.get("sha256"), result=relative.startswith("results/"))
            except (OSError, ValueError, KeyError, TypeError) as error:
                errors.append("auxiliary provenance: " + str(error))
    if current != [manifest.get("frontier_artifact_id")]:
        errors.append(f"historical frontier differs from its declared manifest identity: {current}")
    result_files = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "results").glob("nsc-*.json")
    }
    development = ROOT / "results" / "development"
    if development.is_dir():
        result_files.update(
            path.relative_to(ROOT).as_posix()
            for path in development.glob("*.json")
        )
    try:
        development_outputs = check_development_snapshot()
        if development_outputs & outputs:
            errors.append("development snapshot overlaps the frozen release")
    except Exception as exc:
        errors.append(f"development snapshot validation failed: {exc}")
        development_outputs = set()
    receipt_outputs = set()
    try:
        snapshot, receipts = retained_receipts()
        if snapshot is None:
            raise ValueError("live publication checks require NSC-CONSOLIDATED-LAB-v2 receipts")
        for relative in result_files & set(receipts):
            row = receipts[relative]
            if sha256(ROOT / relative) != row["sha256"] or (ROOT / relative).stat().st_size != row["bytes"]:
                errors.append("retained result differs from its receipt: " + relative)
            receipt_outputs.add(relative)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("retained result receipts: " + str(error))
    local_outputs = {relative for relative in outputs if (ROOT / relative).is_file()}
    if result_files != local_outputs | development_outputs | receipt_outputs:
        errors.append("tracked NSC result files differ from authenticated local/archived closure")


def validate_historical_release(manifest: dict) -> None:
    """Run the original release contract on pinned Git bytes, without old checkouts."""
    baseline = _baseline()
    if RESULT_MANIFEST.read_bytes() != _git_bytes(baseline, "results/manifest.json"):
        raise ValueError("immutable release manifest differs from its preserved edition")
    legacy = historical_checker(baseline)
    legacy.validate_checkout(manifest)


def check_provenance(errors: list[str], files: list[Path]) -> None:
    try:
        value = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"cannot read publication provenance: {exc}")
        return
    if value.get("source_commit") != "ff2cf2722b966589b98a61accdbb6cee819a58c7":
        errors.append("publication provenance has the wrong source commit")
    try:
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        metadata = json.loads((ROOT / "paper/metadata.json").read_text())
        catalog = json.loads((ROOT / "paper/catalog.json").read_text())
        foundation = next(row for row in catalog["documents"] if row["id"] == "foundation")
        archived = json.loads((ROOT / foundation["metadata"]).read_text())
        if value.get("public_version") != project["project"]["version"]:
            errors.append("publication provenance has the wrong current package version")
        if value.get("foundation_version") != archived["version"]:
            errors.append("publication provenance has the wrong frozen foundation version")
        if value.get("manuscript_revision") != metadata["revision"]:
            errors.append("publication provenance has the wrong manuscript revision")
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("publication version routing: " + str(error))
    try:
        from build_publication_provenance import evidence_index
        expected = evidence_index(ROOT)
        rows = value.get("evidence_locations", [])
        locations = {row["path"]: row for row in rows}
        if len(rows) != len(locations) or set(locations) != set(expected):
            raise ValueError("publication evidence location set differs from frozen declarations")
        for relative, identity in expected.items():
            row = locations[relative]
            if row.get("sha256") != identity["sha256"] or row.get("bytes") != identity["bytes"]:
                raise ValueError("publication evidence identity differs: " + relative)
            disposition = row.get("archive_disposition")
            if disposition is not None and (
                    disposition.get("sha256") != identity["sha256"]
                    or disposition.get("doi") != ARCHIVE_DOI):
                raise ValueError("publication archival disposition differs: " + relative)
            status = authenticate_frozen_path(relative, identity["sha256"], result=True)
            if row.get("status") != status:
                raise ValueError("publication local/archive status differs: " + relative)
        if value.get("lab_snapshot_sha256") != sha256(ROOT / "docs/lab-snapshot.json"):
            raise ValueError("publication retained-file receipt hash differs")
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("publication evidence locations: " + str(error))
    declared = {entry["path"]: entry for entry in value.get("included_paths", [])}
    actual = {
        path.relative_to(ROOT).as_posix()
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
    for item in value.get("scientific_build_inputs", []):
        try:
            path = checked_path(item["path"])
            if sha256(path) != item["sha256"]:
                errors.append("paper scientific input mismatch: " + item["path"])
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append("paper scientific input: " + str(error))


def checked_path(relative: str) -> Path:
    path = Path(relative)
    candidate = ROOT / path
    if (path.is_absolute() or ".." in path.parts
            or not candidate.resolve().is_relative_to(ROOT.resolve())
            or not candidate.is_file()):
        raise ValueError("missing or unsafe publication path: " + str(relative))
    return candidate


def archive_context(row: dict) -> dict[str, Path]:
    """Authenticate exact historical bytes while retaining their logical paths."""
    receipt = json.loads(checked_path(row["archive_manifest"]).read_text(encoding="utf-8"))
    if (receipt.get("schema") != "NSC-PAPER-ARCHIVE-v1"
            or receipt.get("id") != row["id"]
            or receipt.get("version") != row["version"]
            or not re.fullmatch("[0-9a-f]{40}", receipt.get("source_commit", ""))):
        raise ValueError("invalid foundation archive identity")
    context = {}
    for item in receipt.get("files", []):
        original = item["original_path"]
        logical = Path(original)
        if logical.is_absolute() or ".." in logical.parts or original in context:
            raise ValueError("invalid or repeated archive logical path")
        stored = checked_path(item["storage_path"])
        if sha256(stored) != item["sha256"]:
            raise ValueError("archived bytes differ: " + original)
        context[original] = stored
    for original in ("paper/nested-space-cosmology.md", "paper/nested-space-cosmology.pdf",
                     "paper/metadata.json", "paper/build-manifest.json",
                     "paper/claim-bindings.json", "scripts/build_paper.py"):
        if original not in context:
            raise ValueError("incomplete foundation archive: " + original)
    build = json.loads(context["paper/build-manifest.json"].read_text(encoding="utf-8"))
    for path_key, hash_key in (("source", "source_sha256"), ("metadata", "metadata_sha256"),
                               ("output", "pdf_sha256")):
        if build[path_key] not in context or sha256(context[build[path_key]]) != build[hash_key]:
            raise ValueError("historical build binding differs: " + path_key)
    for item in build.get("scientific_build_inputs", []):
        path = context.get(item["path"])
        if path is None:
            path = checked_path(item["path"])
        if sha256(path) != item["sha256"]:
            raise ValueError("historical scientific build input differs: " + item["path"])
    return context


def check_paper_catalog(errors: list[str]) -> None:
    """One living main manuscript; prior science remains an authenticated archive."""
    try:
        catalog = json.loads((ROOT / "paper/catalog.json").read_text(encoding="utf-8"))
        if catalog.get("schema") != "NSC-PAPER-CATALOG-v3":
            raise ValueError("unexpected catalog schema")
        rows = catalog["documents"]
        documents = {row["id"]: row for row in rows}
        if len(rows) != 3 or set(documents) != {"main", "foundation", "local-gate"}:
            raise ValueError("main manuscript and both historical science roles are required")
        if (catalog.get("current_document") != "main"
                or catalog.get("current_status") != "WORKING_DRAFT"
                or [row["id"] for row in rows if row.get("active")] != ["main"]):
            raise ValueError("exactly one active main manuscript is required")
        main = documents["main"]
        if (main.get("role") != "main_manuscript" or main.get("status") != "WORKING_DRAFT"
                or main.get("submission_ready") is not False
                or main.get("author_review") != "pending"):
            raise ValueError("main manuscript must remain an unsubmitted working draft")
        if (main["source"] != "paper/nested-space-cosmology.md"
                or main["pdf"] != "paper/nested-space-cosmology.pdf"
                or main["manifest"] != "paper/build-manifest.json"):
            raise ValueError("canonical main manuscript paths differ")
        build = json.loads(checked_path(main["manifest"]).read_text(encoding="utf-8"))
        metadata = json.loads(checked_path(main["metadata"]).read_text(encoding="utf-8"))
        for field, build_field, hash_field in (("source", "source", "source_sha256"),
                                              ("pdf", "output", "pdf_sha256"),
                                              ("metadata", "metadata", "metadata_sha256")):
            if main[field] != build[build_field] or sha256(checked_path(main[field])) != build[hash_field]:
                raise ValueError("current manuscript requires a matching build: " + field)
        if (main["version"] != build["version"] or main["version"] != metadata["version"]
                or main["date"] != metadata["date_iso"]):
            raise ValueError("current manuscript revision metadata differs")
        foundation = documents["foundation"]
        if foundation.get("role") != "archived_foundation" or foundation.get("active") is not False:
            raise ValueError("foundation is a historical edition")
        context = archive_context(foundation)
        for field, logical in (("source", "paper/nested-space-cosmology.md"),
                               ("pdf", "paper/nested-space-cosmology.pdf"),
                               ("manifest", "paper/build-manifest.json"),
                               ("metadata", "paper/metadata.json")):
            if checked_path(foundation[field]) != context[logical]:
                raise ValueError("foundation catalog differs from archive: " + field)
        companion = documents["local-gate"]
        if (companion.get("role") != "archived_supporting_result"
                or companion.get("active") is not False
                or companion.get("status") != "MEASURED_FINITE_REALIZATION"
                or companion.get("submission_ready") is not False):
            raise ValueError("finite companion is a scoped historical result")
        for edition in [companion, *companion.get("historical_editions", [])]:
            prior = json.loads(checked_path(edition["manifest"]).read_text(encoding="utf-8"))
            for relative, digest in ((edition["source"], prior["inputs"][edition["source"]]),
                                     (edition["pdf"], prior["pdf_sha256"]),
                                     (edition["source_archive"], prior["source_archive_sha256"])):
                if (ROOT / relative).is_file():
                    if sha256(checked_path(relative)) != digest:
                        raise ValueError("historical companion source, archive or PDF differs")
                else:
                    authenticate_frozen_path(relative, digest, result=True)
        if not companion.get("historical_editions"):
            raise ValueError("historical companion edition must remain available")
        if (catalog.get("finite_result_status") != "MEASURED_FINITE_REALIZATION"
                or catalog.get("physical_local_gate_status") != "OPEN"):
            raise ValueError("finite evidence and open physical application are distinct")
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append(f"paper catalog: {error}")


def check_markdown_links(errors: list[str], path: Path, text: str) -> None:
    for match in LINK_RE.finditer(text):
        raw_target = match.group(1)
        # Plain-text arXiv extractions can join a set notation and equation
        # number as ``\subseteq[n](14)``. This is not a Markdown hyperlink.
        if (raw_target.isdecimal()
                and re.search(r"\\subseteq\s*$", text[:match.start()])
                and re.match(r"[ \t]*(?:\r?\n|$)", text[match.end():])):
            continue
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


def historical_content() -> set[str]:
    """Recognize exact migrated provenance, never exempt current entrypoints."""
    candidates = {}
    snapshot, receipts = retained_receipts()
    if snapshot is not None:
        candidates.update({name: row["import_sha256"] for name, row in receipts.items()
                           if row.get("sha256") == row.get("import_sha256")
                           and row.get("bytes") == row.get("import_bytes")
                           and not name.startswith(("archive/", "."))})
    elif not (ROOT / "results/release-spec.json").is_file():
        # Compatibility for small, isolated historical-link test fixtures only.
        # A real checkout never routes through the retired migration manifest.
        fixture = ROOT / "results/reassembly/manifest.json"
        if fixture.is_file():
            record = json.loads(fixture.read_text(encoding="utf-8"))
            candidates.update({row["destination"]: row["sha256"] for row in record["rows"]
                               if row["origin"] in {"lab", "w2"}
                               and not row["destination"].startswith(("archive/", "."))})
    for manifest in (ROOT / "results/development").glob("nsc-cosmological-atlas*/manifest.json"):
        record = json.loads(manifest.read_text(encoding="utf-8"))
        candidates.update({path: digest for path, digest in record.get("source_bindings", {}).items()
                           if path.startswith("paper/atlas-working-draft") and path.endswith(".md")})
    for manifest in (ROOT / "paper/archive").glob("*/archive-manifest.json"):
        record = json.loads(manifest.read_text(encoding="utf-8"))
        candidates.update({row["storage_path"]: row["sha256"] for row in record["files"]
                           if row["storage_path"].endswith(".md")})
    # These pre-existing historical notes were already in the publication root,
    # rather than imported by the layout migration. Bind their original bytes.
    from replay_history import BASELINE
    environment = os.environ.copy()
    objects = ROOT / "archive/source-history/objects"
    if objects.is_dir():
        environment["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(objects)
    for relative in ("docs/PLAN.md", "lit/cosmo_v3.md"):
        if not (ROOT / relative).is_file():
            continue
        # This global note inventory is independent of current-record reads.
        preserved = subprocess.run(["git", "show", f"{BASELINE}:{relative}"],
                                   cwd=ROOT, env=environment, capture_output=True,
                                   timeout=30)
        if preserved.returncode == 0:
            candidates[relative] = hashlib.sha256(preserved.stdout).hexdigest()
    return {path for path, digest in candidates.items()
            if path not in GITHUB_MARKDOWN_ENTRYPOINTS
            and (ROOT / path).is_file() and sha256(ROOT / path) == digest}


def check_public_boundary(errors: list[str], files: list[Path]) -> None:
    historical = historical_content()
    relative_files = {path.relative_to(ROOT).as_posix() for path in files}
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
        if path.suffix.lower() not in {".md", ".py", ".cpp", ".json", ".toml", ".txt", ".yml", ".yaml", ".cff", ""}:
            continue
        # A .gitignore exists to name the paths that must stay out of the
        # tree. Naming a private directory there is the fix, not the leak.
        if path.name == ".gitignore":
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        relative = path.relative_to(ROOT).as_posix()
        for literal, label in FORBIDDEN_LITERAL.items():
            if literal in text and relative not in historical:
                errors.append(f"{label} appears in {path.relative_to(ROOT)}")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                errors.append(f"possible {label} appears in {path.relative_to(ROOT)}")
        if relative in GITHUB_MARKDOWN_ENTRYPOINTS:
            for label, pattern in UNSUPPORTED_GITHUB_MATH.items():
                if pattern.search(text):
                    errors.append(f"unsupported GitHub math {label} in {relative}")
        # Immutable laboratory source copies retain their original relative
        # citation paths. Their bytes/inventory are checked by the nested
        # snapshot verifier; the reader-facing proof in docs/ is link-checked.
        if (path.suffix.lower() == ".md"
                and not relative.startswith("paper/nested-quality-evidence/")
                and relative not in historical):
            check_markdown_links(errors, path, text)


def check_bulk_response_history(errors: list[str]) -> None:
    """Select bulk-v1's named original source context, never a failure fallback."""
    try:
        from replay_history import check_bulk_v1
        check_bulk_v1()
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("bulk response original-context evidence: " + str(error))


def main() -> int:
    errors: list[str] = []
    files = public_files()
    check_public_boundary(errors, files)
    check_manifest(errors)
    check_paper(errors)
    check_paper_catalog(errors)
    # A compact analytic successor authenticates its stored prerequisites;
    # publication checks never rerun the predecessor field integrations.
    try:
        from record_nsc_atlas_observer_inference import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("observer inference evidence: " + str(error))
    try:
        from record_nsc_atlas_scalar_wall_consistency import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("scalar/wall consistency evidence: " + str(error))
    try:
        from record_nsc_atlas_characteristic_initialization import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("characteristic initialization evidence: " + str(error))
    try:
        from record_nsc_atlas_source_mass import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("source-leg mass evidence: " + str(error))
    try:
        from record_nsc_atlas_scalar_normal_force import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("scalar normal-force evidence: " + str(error))
    try:
        from record_nsc_atlas_cap_wall_response import check, RELATIVE_OUTPUT
        check(ROOT / RELATIVE_OUTPUT)
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append("cap/wall response evidence: " + str(error))
    check_bulk_response_history(errors)
    check_provenance(errors, files)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"publication check failed with {len(errors)} finding(s)", file=sys.stderr)
        return 1
    count = json.loads(RESULT_MANIFEST.read_text())["result_count"]
    print(f"publication check passed: {len(files)} curated files, {count} result steps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
