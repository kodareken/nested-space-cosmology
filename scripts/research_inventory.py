"""Inventory the preserved research baseline without running scientific code.

Git plumbing supplies immutable tracked objects; a hash-checked preservation
receipt overlays the original dirty files. Current working files are never used
as baseline evidence. Git object IDs identify unchanged payloads; document and
dirty-file SHA-256s identify text variants. Reading/hashing is not verification
of the physics. Outputs are deterministic and add no machine-specific root
paths. Historical absolute references are retained and labelled. This
inventory's own outputs and runtime directories are excluded.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess


SCHEMA = "NSC-RESEARCH-INVENTORY-v1"
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("lab/results/research-reset")
EXCLUDED_PARTS = {"__pycache__", ".pytest_cache", ".ruff_cache", "node_modules",
                  "pytest-of-swedi", ".venv", ".co" + "dex", ".build"}
EXCLUDED_PREFIXES = ("lab/results/research-reset/", "lab/results/nested-atlas/",
                     ".pi/subagents/", "build/", "dist/")
TEXT_SUFFIXES = {".md", ".tex", ".bib", ".py", ".cpp", ".toml", ".json"}
MARKDOWN_LINK = re.compile(r"\[[^\]\n]*\]\(([^)\n]+)\)")
PATH_TOKEN = re.compile(r"(?:lab/)?(?:src/|scripts/|tests/|results/|docs/|paper/|archive/)[A-Za-z0-9_.+/@-]+")
FILE_TOKEN = re.compile(r"(?<![\w./-])([A-Za-z0-9_-]+\.(?:py|json|npz|md|tex))(?![\w.])")
ARXIV = re.compile(r"(?:arxiv\.org/(?:abs|html|pdf)/|arXiv[:\s]+)([A-Za-z-]+/\d{7}|\d{4}\.\d{4,5})(?:v\d+)?", re.I)


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def excluded(path):
    return (any(part in EXCLUDED_PARTS for part in PurePosixPath(path).parts)
            or path.startswith(EXCLUDED_PREFIXES)
            or path in {"scripts/research_inventory.py", "tests/test_research_inventory.py",
                        "lab/docs/research-inventory.md"})


def canonical_stem(path):
    name = PurePosixPath(path).stem.replace("_", "-")
    name = re.sub(r"^(?:test|derive|assess|check|run|map|build)-", "", name)
    return re.sub(r"-v\d+$", "", name)


def result_family(path):
    """Index every result root; a directory's numeric rows belong to its family."""
    if path.startswith(("lab/archive/", "archive/")) and "/results/" in path:
        prefix, remainder = path.split("/results/", 1)
        # Preserve archive identity; never blend a retired family with live files.
        return prefix + "/" + result_family("results/" + remainder)
    for root in ("lab/results/development/", "lab/results/", "results/development/", "results/"):
        if path.startswith(root):
            parts = path[len(root):].split("/")
            # This exploratory import contains several independent calculations.
            if parts[0] == "child-neck-v1" and len(parts) > 1:
                return root + "child-neck/" + canonical_stem(parts[1])
            if parts[0] == "artifacts" and len(parts) > 1:
                return root + "artifacts/" + canonical_stem(parts[1])
            return root + canonical_stem(parts[0])
    return None


def is_document(path):
    return PurePosixPath(path).suffix.lower() in {".md", ".tex"} and (
        path.startswith(("docs/", "lab/docs/", "paper/", "lab/paper/", "lab/archive/", "archive/"))
        or path in {"AGENTS.md", "README.md", "lab/AGENTS.md", "lab/README.md", "lab/PLAN.md", "lab/handover.md"})


def kind(path):
    if "/.source-history/" in path:
        return "pinned_source_history"
    if path.startswith(("lab/archive/", "archive/")):
        return "historical_archive"
    if path.startswith(("paper/", "lab/paper/")):
        return "publication"
    if result_family(path):
        return "result"
    if is_document(path):
        return "document"
    if "/tests/" in "/" + path:
        return "test"
    if "/scripts/" in "/" + path:
        return "driver"
    if "/src/" in "/" + path:
        return "source"
    return "support"


def strand(path):
    name = path.lower().replace("_", "-")
    if name.startswith(("lab/archive/", "archive/")):
        return "historical_archive"
    if "child-neck" in name or "child/neck" in name:
        return "child_neck"
    if any(word in name for word in ("influence-", "matryoshka", "cosmic-band", "cosmic-shell", "constant-table", "centenary", "programme-index")):
        return "w2_influence_interpretations"
    if any(word in name for word in ("nsc-wall-", "nsc-march-", "nsc-t30-", "nsc-charged-", "nsc-tilt-", "nsc-two-ends-", "nsc-oneway-", "nsc-isotropy-", "nsc-fringe-", "nsc-fee-", "nsc-extinction-", "nsc-interior-", "nsc-continuation-", "nsc-trap-", "nsc-time-awareness", "nsc-lowband-", "nsc-q-")):
        return "w2_radial_extensions"
    if "nsc-rn" in name:
        return "rn_radial"
    if "nsc-discovery" in name:
        return "finite_discovery"
    if any(word in name for word in ("nested-parent-child", "nested-qualities", "energy-transfer", "local-observer", "regulated-recursion", "geometric-chain", "finite-turnover", "imbalance-turnover")):
        return "finite_nesting"
    if any(word in name for word in ("spherical-conformal", "conformal-local", "conformal-memory")):
        return "conformal_realization"
    if any(word in name for word in ("spherical", "regeneration", "galerkin", "coupled-local", "conformal-adm", "regional-", "geometry-graded", "local-boundary")):
        return "spherical_action"
    if any(word in name for word in ("nsc-ks", "-ks-", "nsc-027", "incoming", "middle-source", "prepared-history", "endpoint", "local-gate", "collar", "uv-", "subgap", "extended-gate", "completion-gate", "history-gate")):
        return "paused_incoming_gate"
    if any(word in name for word in ("torsion", "neutrino", "charged-sector", "nuclear", "baryon", "proton", "neutron", "electron", "particle", "dark-matter", "higgs", "standard-model")):
        return "matter_torsion"
    if any(word in name for word in ("casimir", "determinant", "spectral", "quantum", "fermion", "dirac", "anomaly", "flow-", "renormal", "holonomy", "heat-kernel", "measure", "vacuum", "landau", "lll-", "jost", "retarded", "unitary", "resolvent", "ctp", "normalization", "gauge", "zeta", "radial-spectrum", "finite-terms", "chiral", "threshold", "tail-", "weyl", "massless", "source", "state-regulator", "hamiltonian", "boundary-response", "reference-band", "spatial-reference", "unruh")):
        return "spectral_quantum"
    if any(word in name for word in ("compact", "neck", "horizon", "child", "parent", "warped", "gravity", "curvature", "black-universe", "black-hole", "bounce", "clock", "lorentz", "pg-", "geometry", "geometric", "scale-binding", "mmp", "angular-stress", "shape-response", "observable-bridge")):
        return "geometry_transition"
    if path.startswith(("paper/", "lab/paper/")):
        return "publication"
    if any(word in name for word in ("claim", "audit", "prior-art", "arxiv", "readme", "instructions", "plan", "handover", "code-map", "provenance", "snapshot", "storage", "consolidation", "reproduc")):
        return "governance_provenance"
    return "other_research"


class GitObjects:
    """Unfiltered Git plumbing is required for binary object identities."""
    def __init__(self, root):
        self.process = subprocess.Popen(["git", "cat-file", "--batch"], cwd=root,
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    def read(self, oid):
        self.process.stdin.write((oid + "\n").encode("ascii"))
        self.process.stdin.flush()
        header = self.process.stdout.readline().decode().split()
        if len(header) != 3 or header[1] != "blob":
            raise ValueError(f"not an available blob: {oid}")
        size = int(header[2])
        data = self.process.stdout.read(size)
        if len(data) != size or self.process.stdout.read(1) != b"\n":
            raise ValueError("truncated Git object stream")
        return data

    def close(self):
        self.process.stdin.close()
        self.process.stdout.close()
        self.process.wait()


def baseline(root, preservation):
    receipt_bytes = (preservation / "receipt.json").read_bytes()
    receipt = json.loads(receipt_bytes)
    head = receipt["head"]
    tree = subprocess.check_output(["git", "ls-tree", "-rlz", head], cwd=root)
    rows = {}
    omitted = []
    for entry in tree.split(b"\0"):
        if not entry:
            continue
        metadata, rawpath = entry.split(b"\t", 1)
        mode, objtype, oid, size = metadata.decode().split()
        path = rawpath.decode("utf-8")
        if objtype != "blob":
            continue
        if excluded(path):
            omitted.append(path)
            continue
        rows[path] = {"path": path, "bytes": int(size), "git_blob": oid, "mode": mode,
                      "origin": "committed", "kind": kind(path), "strand": strand(path)}
    for item in receipt["files"]:
        path = item["path"].replace("\\", "/")
        copy = preservation / "files" / path
        data = copy.read_bytes()
        if digest(data) != item["sha256"] or len(data) != item["bytes"]:
            raise ValueError(f"preserved input has drifted: {receipt['repository']}/{path}")
        if excluded(path):
            omitted.append(path)
            continue
        previous = rows.get(path)
        rows[path] = {"path": path, "bytes": len(data), "sha256": item["sha256"],
                      "baseline_git_blob": previous["git_blob"] if previous else None,
                      "origin": "modified" if previous else "untracked", "kind": kind(path),
                      "strand": strand(path)}
    # Deletions have no backup payload and must not be resurrected from HEAD.
    status = (preservation / "status.bin").read_bytes().split(b"\0")
    deleted = []
    for entry in status:
        if len(entry) > 3 and b"D" in entry[:2]:
            path = entry[3:].decode("utf-8").replace("\\", "/")
            rows.pop(path, None)
            deleted.append(path)
    summary = {"repository": receipt["repository"], "head": head, "branch": receipt["branch"],
               "preservation_receipt_sha256": digest(receipt_bytes),
               "preserved_dirty_files_verified": len(receipt["files"]),
               "excluded_runtime_paths": sorted(set(omitted)), "deleted_at_baseline": sorted(deleted),
               "inventory_files": len(rows), "inventory_bytes": sum(row["bytes"] for row in rows.values())}
    return rows, summary


def resolve_link(source, raw, paths):
    target = raw.strip().strip("<>").split(' "', 1)[0]
    if target.startswith(("https://", "http://", "mailto:", "doi:", "arxiv:")):
        return {"target": target, "status": "external_reference"}
    if target.startswith("#"):
        return {"target": target, "status": "local_anchor_unchecked"}
    target = target.split("#", 1)[0].split("?", 1)[0]
    if not target:
        return {"target": raw, "status": "empty"}
    if target.startswith(("/", "~", "file:")) or re.match(r"^[A-Za-z]:", target):
        return {"target": raw, "status": "historical_or_host_absolute"}
    target = target.replace("\\", "/")
    resolved = posixpath.normpath(posixpath.join(posixpath.dirname(source), target))
    if resolved in paths:
        return {"target": target, "resolved": resolved, "status": "resolved"}
    if any(path.startswith(resolved.rstrip("/") + "/") for path in paths):
        return {"target": target, "resolved": resolved, "status": "resolved_directory"}
    return {"target": target, "resolved": resolved, "status": "missing_in_preserved_baseline"}


def extract_document(path, data):
    text = data.decode("utf-8", errors="replace")
    lines = text.splitlines()
    headings = [{"line": i + 1, "text": line.lstrip("# ")} for i, line in enumerate(lines) if line.startswith("#")]
    evidence = []
    for i, line in enumerate(lines):
        if re.search(r"\b(scope|supersed\w*|withdraw\w*|not (?:a |an )?(?:proof|claim|certif\w*)|remains? open|status|does not|cannot|unresolved|paused)\b", line, re.I):
            evidence.append({"line": i + 1, "text": line[:500]})
    return {"sha256": digest(data), "bytes": len(data), "lines": len(lines),
            "title": headings[0]["text"] if headings else PurePosixPath(path).name,
            "headings": headings, "scope_and_status_excerpts": evidence,
            "arxiv_ids": sorted(set(ARXIV.findall(text))),
            "raw_links": MARKDOWN_LINK.findall(text), "path_tokens": sorted(set(PATH_TOKEN.findall(text))),
            "file_tokens": sorted(set(FILE_TOKEN.findall(text))),
            "full_text_scanned": True, "equations_independently_verified": False}


def record_metadata(data):
    try:
        obj = json.loads(data)
    except (ValueError, UnicodeDecodeError):
        return {"parse": "invalid_or_non_utf8_json"}
    if not isinstance(obj, dict):
        return {"parse": "json_non_object"}
    fields = {}
    for key in ("schema", "api_version", "status", "verdict", "scope", "claim_status", "supersedes", "producer_commit", "source_commit"):
        value = obj.get(key)
        if isinstance(value, (str, int, bool)):
            fields[key] = str(value)[:600]
    refs = set()
    def walk(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if isinstance(key, str) and "/" in key and PurePosixPath(key).suffix in TEXT_SUFFIXES | {".npz"}:
                    refs.add(key)
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)
        elif isinstance(value, str) and len(value) < 300 and re.search(r"\.(?:json|npz|py|md|tex|toml|bib|cpp)(?:$|#)", value):
            if not value.startswith(("http", "sha")) and "\n" not in value:
                refs.add(value)
    walk(obj)
    return {"parse": "json", "declared": fields, "references": sorted(refs)}


def build(root, preservation, review_path):
    review = json.loads(review_path.read_bytes())
    corpora, snapshots, files, families, link_issues = {}, [], [], {}, []
    family_references = defaultdict(set)
    repository_paths = {}
    objects = GitObjects(root)
    cache = {}
    try:
        for repository in ("BlackHoles-Infinity", "BlackHoles-Infinity-w2"):
            saved = preservation / repository
            rows, snapshot = baseline(root, saved)
            snapshots.append(snapshot)
            paths = set(rows)
            repository_paths[repository] = paths
            bystem = defaultdict(list)
            byname = defaultdict(list)
            for path in paths:
                bystem[canonical_stem(path).removeprefix("nsc-")].append(path)
                byname[PurePosixPath(path).name].append(path)
            def read(row):
                if row["origin"] != "committed":
                    return (saved / "files" / row["path"]).read_bytes()
                return objects.read(row["git_blob"])
            for path, row in sorted(rows.items()):
                files.append({"repository": repository, **row})
                family = result_family(path)
                if family:
                    fid = repository + ":" + family
                    fam = families.setdefault(fid, {"id": fid, "repository": repository, "family": family,
                        "strand": strand(path), "files": [], "json_records": [], "note_owners": [],
                        "disposition": "unresolved", "review": "inventoried_only"})
                    fam["files"].append(path)
                    if path.endswith(".json"):
                        # Read every result JSON, including its producer/source references.
                        key = row.get("git_blob", row.get("sha256"))
                        if key not in cache:
                            cache[key] = record_metadata(read(row))
                        metadata = cache[key]
                        references = metadata.get("references", [])
                        family_references[fid].update(references)
                        fam["json_records"].append({"path": path,
                            **{name: value for name, value in metadata.items() if name != "references"},
                            "reference_count": len(references),
                            "references_sha256": digest(encoded(references))})
                if not is_document(path):
                    continue
                doc = extract_document(path, read(row))
                sha = doc["sha256"]
                location = {"repository": repository, "path": path, "origin": row["origin"], "strand": row["strand"]}
                references = [resolve_link(path, link, paths) for link in doc.pop("raw_links")]
                for link in references:
                    if link["status"] in {"missing_in_preserved_baseline", "historical_or_host_absolute"}:
                        link_issues.append({"repository": repository, "source": path, **link})
                linked = {link["resolved"] for link in references if link["status"] == "resolved"}
                # Code-style paths are evidence hints, not Markdown link successes.
                hints = set()
                for token in doc.pop("path_tokens"):
                    options = [token, ("lab/" + token) if path.startswith("lab/") else token,
                               posixpath.normpath(posixpath.join(posixpath.dirname(path), token))]
                    hints.update(option for option in options if option in paths)
                for token in doc.pop("file_tokens"):
                    candidates = byname.get(token, [])
                    local = [candidate for candidate in candidates if candidate.startswith("lab/") == path.startswith("lab/")]
                    hints.update(local)
                hints.discard(path)
                names = set(bystem.get(canonical_stem(path).removeprefix("nsc-"), [])) - {path}
                location["links"] = references
                location["related_files"] = [{"path": item, "kind": rows[item]["kind"],
                    "relationship": "explicit_link" if item in linked else "path_mention" if item in hints else "basename_candidate"}
                    for item in sorted(linked | hints | names)]
                reviewed = path in review.get("reviewed_owner_paths", {}).get(repository, [])
                location["review"] = "owner_scope_reviewed" if reviewed else "inventoried_only"
                location["disposition"] = disposition(path, row, review, reviewed)
                if path in review.get("document_overrides", {}):
                    location["review_note"] = review["document_overrides"][path]
                if sha not in corpora:
                    corpora[sha] = {**doc, "locations": []}
                corpora[sha]["locations"].append(location)
                for linked_path in linked | hints | names:
                    rid = result_family(linked_path)
                    if rid:
                        fid = repository + ":" + rid
                        fam = families.setdefault(fid, {"id": fid, "repository": repository, "family": rid,
                            "strand": strand(linked_path), "files": [], "json_records": [], "note_owners": [],
                            "disposition": "unresolved", "review": "inventoried_only"})
                        fam["note_owners"].append({"path": path, "sha256": sha,
                            "owner_candidate": PurePosixPath(path).name not in {"AGENTS.md", "README.md", "active-code-map.md", "claim-ledger.md", "current-result.md", "instructions.md", "PLAN.md", "handover.md", "arxiv-studies.md", "prior-art-and-open-claim.md"},
                            "review": location["review"], "relationship": "explicit_link" if linked_path in linked else "path_mention" if linked_path in hints else "basename_candidate"})
            snapshot["counts_by_kind"] = dict(sorted(Counter(row["kind"] for row in rows.values()).items()))
            snapshot["counts_by_origin"] = dict(sorted(Counter(row["origin"] for row in rows.values()).items()))
    finally:
        objects.close()
    document_checksums = {(loc["repository"], loc["path"]): doc["sha256"]
                          for doc in corpora.values() for loc in doc["locations"]}
    for fam in families.values():
        fam["files"] = sorted(set(fam["files"]))
        fam["note_owners"] = list({(row["path"], row["relationship"]): row for row in fam["note_owners"]}.values())
        fam["note_owners"].sort(key=lambda row: (row["path"], row["relationship"]))
        owners = [row for row in fam["note_owners"] if row["owner_candidate"]]
        reviewed = [row for row in owners if row["review"] == "owner_scope_reviewed" and row["relationship"] != "basename_candidate"]
        fam["review"] = "owner_scope_reviewed" if reviewed else "inventoried_only"
        policy = review["strands"][fam["strand"]]
        fam["disposition"] = policy["default_disposition"]
        fam["admitted_to_atlas_as_physical_result"] = False
        fam["equations_independently_verified"] = False
        fam["owner_resolution"] = "linked_or_mentioned" if any(row["relationship"] != "basename_candidate" for row in owners) else "basename_candidates_only" if owners else "no_owner_identified"
        fam["family_role"] = ("historical_archive_family" if fam["strand"] == "historical_archive" else
            "payload_support_family" if "/artifacts/" in fam["family"] else
            "support_without_json_record" if not fam["json_records"] else
            "historical_campaign_record" if fam["strand"] == "paused_incoming_gate" and re.search(r"(?:trial|iterate|proposal|stall|blocked|search)", fam["family"]) else
            "scientific_record_family")
        source_paths, result_dependencies, unresolved = set(), set(), set()
        paths = repository_paths[fam["repository"]]
        for raw in family_references[fam["id"]]:
            target = raw.split("#", 1)[0].replace("\\", "/")
            options = (["lab/" + target, target] if fam["family"].startswith("lab/") else [target])
            if fam["strand"] == "historical_archive":
                options.insert(0, fam["family"].split("/results/", 1)[0] + "/" + target)
            resolved = next((candidate for candidate in options if candidate in paths), None)
            if resolved:
                dependency = result_family(resolved)
                if dependency and dependency != fam["family"]:
                    result_dependencies.add(dependency)
                elif not dependency:
                    source_paths.add(resolved)
            else:
                unresolved.add(raw)
        fam["dependency_summary"] = {
            "source_test_note_paths": sorted(source_paths),
            "referenced_result_families": sorted(result_dependencies),
            "unresolved_or_historical_references": sorted(unresolved),
            "boundary": "Exact per-record reference sets remain in the identified source JSON; this index condenses resolved numeric payload references to their result family. No deletion is authorised by this graph."}
        # A record's explicit note dependency is stronger than a matching name,
        # even if that note never linked back to this particular successor.
        for note in sorted(source_paths):
            if is_document(note) and PurePosixPath(note).name not in {
                    "AGENTS.md", "README.md", "active-code-map.md", "claim-ledger.md",
                    "current-result.md", "instructions.md", "PLAN.md", "handover.md"}:
                if not any(row["path"] == note and row["relationship"] != "basename_candidate" for row in fam["note_owners"]):
                    checksum = document_checksums.get((fam["repository"], note))
                    reviewed_note = note in review.get("reviewed_owner_paths", {}).get(fam["repository"], [])
                    fam["note_owners"].append({"path": note, "sha256": checksum,
                        "owner_candidate": True, "relationship": "record_bound_note_dependency",
                        "review": "owner_scope_reviewed" if reviewed_note else "inventoried_only"})
                    if fam["owner_resolution"] in {"no_owner_identified", "basename_candidates_only"}:
                        fam["owner_resolution"] = "record_bound_note_dependency"
                    if reviewed_note:
                        fam["review"] = "owner_scope_reviewed"
        fam["note_owners"].sort(key=lambda row: (row["path"], row["relationship"]))
        fam["audit_owner"] = {
            "note": "lab/docs/research-inventory.md",
            "strand_disposition": "lab/results/research-reset/strand-dispositions.json#/strands/" + fam["strand"],
            "responsibility": "Resolve scientific-note ownership and applicability before any new compute or physical atlas admission; preserve the baseline evidence."}
    corpus = sorted(corpora.values(), key=lambda doc: doc["sha256"])
    for doc in corpus:
        doc["locations"].sort(key=lambda loc: (loc["repository"], loc["path"]))
    duplicates = [{"sha256": doc["sha256"], "locations": [{"repository": loc["repository"], "path": loc["path"]} for loc in doc["locations"]]}
                  for doc in corpus if len(doc["locations"]) > 1]
    documents_by_location = {(loc["repository"], loc["path"]): doc["sha256"] for doc in corpus for loc in doc["locations"]}
    alias_audit = []
    for (repository, path), checksum in sorted(documents_by_location.items()):
        if path.startswith("docs/") and (repository, "lab/" + path) in documents_by_location:
            alias_audit.append({"repository": repository, "publication_path": path, "lab_path": "lab/" + path,
                "byte_identical": checksum == documents_by_location[(repository, "lab/" + path)]})
    worktree_variants = []
    for (repository, path), checksum in sorted(documents_by_location.items()):
        other = documents_by_location.get(("BlackHoles-Infinity-w2", path))
        if repository == "BlackHoles-Infinity" and other is not None and other != checksum:
            worktree_variants.append({"path": path, "canonical_sha256": checksum, "w2_sha256": other})
    family_list = sorted(families.values(), key=lambda row: row["id"])
    unknown_owners = [{"repository": fam["repository"], "family": fam["family"], "strand": fam["strand"],
                       "role": fam["family_role"], "source_test_note_paths": fam["dependency_summary"]["source_test_note_paths"]}
                      for fam in family_list if fam["owner_resolution"] == "no_owner_identified"]
    summary = {"schema": SCHEMA, "baseline": preservation.name, "physics_computations_performed": False,
        "source_snapshots": snapshots, "files": len(files), "document_locations": sum(len(doc["locations"]) for doc in corpus),
        "distinct_document_contents": len(corpus), "duplicate_content_groups": len(duplicates),
        "reviewed_owner_locations": sum(loc["review"] == "owner_scope_reviewed" for doc in corpus for loc in doc["locations"]),
        "result_families": len(family_list), "json_results_read": sum(len(fam["json_records"]) for fam in family_list),
        "families_by_strand": dict(sorted(Counter(fam["strand"] for fam in family_list).items())),
        "families_by_owner_resolution": dict(sorted(Counter(fam["owner_resolution"] for fam in family_list).items())),
        "families_by_role": dict(sorted(Counter(fam["family_role"] for fam in family_list).items())),
        "unowned_families_by_role": dict(sorted(Counter(fam["family_role"] for fam in family_list if fam["owner_resolution"] == "no_owner_identified").items())),
        "families_with_reviewed_owner_scope": sum(fam["review"] == "owner_scope_reviewed" for fam in family_list),
        "link_issues": dict(sorted(Counter(row["status"] for row in link_issues).items())),
        "review_boundary": "All inventoried document bytes were read and extracted; named owning notes received scope review. Other derivations, numerical payloads and every equation remain independently unverified. File identity is not scientific validation.",
        "baseline_boundary": "Tracked Git blobs at preserved HEADs plus hash-checked dirty copies; excludes concurrent implementation and runtime paths. Same-content documents remain separate path aliases; historical archive is not current proof."}
    header = {"schema": SCHEMA, "baseline": preservation.name}
    outputs = {
        "baseline-snapshots.json": {**header, "snapshots": snapshots},
        "file-inventory.json": {**header, "files": files},
        "document-corpus.json": {**header, "documents": corpus, "byte_identical_groups": duplicates,
                                 "publication_lab_alias_audit": alias_audit, "same_path_worktree_variants": worktree_variants},
        "result-families.json": {**header, "families": family_list},
        "unowned-families.json": {**header, "families": unknown_owners,
            "boundary": "Unowned means no explicit scientific-note link or mention was resolved in this preserved baseline. Payload support, historical campaigns and scientific records are separated; every item remains preserved and inactive until ownership is reviewed."},
        "link-audit.json": {**header, "issues": sorted(link_issues, key=lambda row: (row["repository"], row["source"], row["target"]))},
        "inventory-summary.json": summary,
    }
    return outputs


def disposition(path, row, review, reviewed):
    if path.startswith(("lab/archive/", "archive/")):
        return "retired_from_active_path_preserve"
    if path in review.get("document_overrides", {}):
        return review["document_overrides"][path]["disposition"]
    if row["origin"] == "untracked":
        return "exploratory_preserve_unreviewed"
    return review["strands"][row["strand"]]["default_disposition"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preservation", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--check", action="store_true", help="Rebuild in memory and compare; never write.")
    args = parser.parse_args(argv)
    destination = args.root / OUTPUT
    review_path = destination / "strand-dispositions.json"
    outputs = build(args.root, args.preservation, review_path)
    mismatches = []
    for name, value in outputs.items():
        path = destination / name
        data = encoded(value)
        if args.check:
            if not path.is_file() or path.read_bytes() != data:
                mismatches.append(name)
        else:
            destination.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    print(json.dumps({"check": args.check, "mismatches": mismatches,
                      **{key: outputs["inventory-summary.json"][key] for key in ("files", "distinct_document_contents", "result_families", "json_results_read", "reviewed_owner_locations")}}, indent=2))
    return bool(mismatches)


if __name__ == "__main__":
    raise SystemExit(main())
