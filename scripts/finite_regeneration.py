"""Authenticate, build and package the focused finite-regeneration article.

Scientific numerical inputs come only from an explicit immutable commit.
Recorded historical source bytes are authenticated separately; a current
guarded producer never impersonates its original hash. No simulation runs.
"""
from __future__ import annotations

import ast
from datetime import datetime, timezone
import gzip
from hashlib import sha1, sha256
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("paper/finite-regeneration")
EVIDENCE = SOURCE / "evidence/snapshot.json"
PDF = Path("paper/finite-regeneration.pdf")
ARCHIVE = Path("paper/finite-regeneration-arxiv.tar.gz")
MANIFEST = Path("paper/finite-regeneration-manifest.json")
METADATA = SOURCE / "arxiv-metadata.json"
DRYRUN = SOURCE / "arxiv-dryrun.json"
COMPILER = ROOT / ".build/tectonic-0.17.0/tectonic"
HISTORICAL = "c2d5fb6f02722d9f69ec70a7145ea3f2287266df"
OLDER = ("fa6855d60655bdfa7eee9ed7ad8d9193639a0aa8", "5f10ecd365843d1616e50eb16a20d7acd8377e2c")
RECORDS = (
    "nsc-spherical-conformal-episode-v1", "nsc-spherical-conformal-frames-v1",
    "nsc-spherical-conformal-episode-v2", "nsc-spherical-conformal-transport-v2",
    "nsc-spherical-conformal-clock-v2", "nsc-spherical-conformal-curvature-v1",
    "nsc-spherical-conformal-source-controls-v1", "nsc-conformal-local-response-v2",
    "nsc-conformal-memory-control-v1", "nsc-conformal-local-response-feasibility-v1",
    "nsc-spherical-coupling-refinement-v5",
)
HASH_GROUPS = {"source_bindings", "source_bindings_before", "source_bindings_after",
    "source_hashes", "source_hashes_before", "source_hashes_after", "hashes_before", "hashes_after",
    "predecessor_hashes", "source_trajectory_hashes", "bound_sources", "code_identities", "production_bindings"}
ALIASES = {
    "v5_json": "nsc-spherical-coupling-refinement-v5.json", "v5_npz": "nsc-spherical-coupling-refinement-v5.npz",
    "legacy_episode_json": "nsc-spherical-feedback-episode-v1.json", "legacy_episode_npz": "nsc-spherical-feedback-episode-v1.npz",
    "v1_json": "nsc-spherical-conformal-episode-v1.json", "v1_npz": "nsc-spherical-conformal-episode-v1.npz",
    "v1_source_json": "nsc-spherical-conformal-episode-v1.json", "v1_source_npz": "nsc-spherical-conformal-episode-v1.npz",
    "v2_source_json": "nsc-spherical-conformal-episode-v2.json", "v2_source_npz": "nsc-spherical-conformal-episode-v2.npz",
    "episode_v1": "nsc-spherical-conformal-episode-v1.json", "episode_v2": "nsc-spherical-conformal-episode-v2.json",
    "payload_v2": "nsc-spherical-conformal-episode-v2.npz", "v1_frames_json": "nsc-spherical-conformal-frames-v1.json",
    "v1_frames_npz": "nsc-spherical-conformal-frames-v1.npz",
    "episode_npz": "nsc-spherical-conformal-episode-v1.npz", "episode_json": "nsc-spherical-conformal-episode-v1.json",
    "manuscript_pdf": "paper/nested-space-cosmology.pdf", "manuscript_md": "paper/nested-space-cosmology.md",
    "companion_pdf": "paper/local-incoming-gate-draft.pdf",
}


def digest(data):
    return sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def read_git(commit, path):
    return git("show", f"{commit}:{path}")


def blob_id(data):
    return sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def safe_name(name):
    p = PurePosixPath(name)
    if not name or p.is_absolute() or ".." in p.parts or "\\" in name or ":" in name or p.as_posix() != name:
        raise ValueError("unsafe artifact path: " + name)
    return p


def protected_files():
    """The catalog may change; historical source, PDFs and manifests may not."""
    paths = ("paper/nested-space-cosmology.md", "paper/nested-space-cosmology.pdf",
             "paper/build-manifest.json", "paper/metadata.json", "paper/local-incoming-gate-draft.pdf",
             "paper/local-incoming-gate-draft-source.tar.gz", "paper/local-gate-draft-manifest.json",
             "scripts/local_gate_draft.py")
    result = {name: digest((ROOT / name).read_bytes()) for name in paths}
    for folder in ("paper/local-gate-draft", "paper/local-gate-evidence", "paper/nested-quality-evidence"):
        for path in (ROOT / folder).rglob("*"):
            if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc":
                result[path.relative_to(ROOT).as_posix()] = digest(path.read_bytes())
    return result


class Graph:
    """A content-authenticated record/source graph, with explicit byte origins."""
    def __init__(self, science):
        if not re.fullmatch("[0-9a-f]{40}", science):
            raise ValueError("full immutable science commit required")
        git("cat-file", "-e", science + "^{commit}")
        self.science = science
        self.refs = tuple(dict.fromkeys((science, HISTORICAL, *OLDER)))
        self.index = {}
        self.small_hashes = {}
        self.bytes = {}
        self.nodes = {}
        self.edges = []
        self.pending = []
        for ref in self.refs:
            tree = {}
            for line in git("ls-tree", "-rl", ref).decode().splitlines():
                header, path = line.split("\t", 1)
                mode, kind, oid, size = header.split()
                if kind == "blob":
                    tree[path] = (oid, int(size))
            self.index[ref] = tree

    def data(self, ref, path):
        oid, _ = self.index[ref][path]
        if oid not in self.bytes:
            self.bytes[oid] = git("cat-file", "blob", oid)
        return self.bytes[oid]

    def small_index(self):
        if self.small_hashes:
            return
        # Hash each shared blob once; large numerical arrays are resolved by path.
        for ref in self.refs:
            for path, (oid, size) in self.index[ref].items():
                if size <= 1024 * 1024 and Path(path).suffix in (".py", ".json", ".md", ".toml", ".cpp", ".tex", ".bib"):
                    key = digest(self.data(ref, path))
                    self.small_hashes.setdefault(key, (ref, path))

    def resolve(self, hint, expected, context):
        candidates = []
        if "/" in hint and not hint.startswith("/"):
            candidates += [hint, "lab/" + hint]
        basename = ALIASES.get(hint, hint)
        if "/" in basename:
            candidates.append(basename)
        if basename.endswith((".json", ".npz")):
            candidates += ["lab/results/development/" + basename]
        if hint in ("npz", "payload_npz", "payload_sha256"):
            candidates.append(str(Path(context).with_suffix(".npz")))
        for ref in self.refs:
            for path in candidates:
                if path in self.index[ref] and digest(self.data(ref, path)) == expected:
                    if Path(path).suffix in (".npz", ".json") and ref != self.science:
                        raise ValueError("numerical evidence must be resident at the science pin: " + path)
                    return ref, path
        self.small_index()
        if expected in self.small_hashes:
            ref, path = self.small_hashes[expected]
            if Path(path).suffix in (".npz", ".json") and ref != self.science:
                raise ValueError("historical numerical substitution refused: " + path)
            return ref, path
        if "npz" in hint or "payload" in hint:
            for target, (oid, _) in self.index[self.science].items():
                if target.endswith(".npz") and digest(self.data(self.science, target)) == expected:
                    return self.science, target
        raise ValueError(f"unresolved declared dependency {hint}={expected} in {context}")

    def add(self, ref, path, expected=None, parent=None, label="declared"):
        data = self.data(ref, path)
        actual = digest(data)
        if expected is not None and expected != actual:
            raise ValueError("evidence hash mismatch: " + path)
        key = f"{path}@{actual}"
        if parent:
            self.edges.append({"from": parent, "to": key, "role": label})
        if key not in self.nodes:
            oid, size = self.index[ref][path]
            self.nodes[key] = {"path": path, "commit": ref, "git_blob": oid, "sha256": actual, "bytes": size}
            if not label.startswith("preservation:"):
                self.pending.append((key, ref, path, data))
        return key

    def record_edges(self, key, ref, path, record):
        for group in HASH_GROUPS:
            values = record.get(group)
            if not isinstance(values, dict):
                continue
            for hint, expected in values.items():
                if not isinstance(expected, str) or not re.fullmatch("[0-9a-f]{64}", expected):
                    continue
                origin, target = self.resolve(hint, expected, path)
                preservation = hint in ("legacy_episode_json", "legacy_episode_npz", "manuscript_pdf", "manuscript_md", "companion_pdf")
                role = ("preservation:" if preservation else "") + group + ":" + hint
                self.add(origin, target, expected, key, role)
        for field in ("source_code_sha256", "reducer_sha256", "consumer_helpers_sha256", "frozen_feasibility_module_sha256"):
            expected = record.get(field)
            if isinstance(expected, str) and re.fullmatch("[0-9a-f]{64}", expected):
                origin, target = self.resolve(field, expected, path)
                self.add(origin, target, expected, key, field)
        expected = record.get("payload_sha256")
        if isinstance(expected, dict):
            expected = expected.get("npz")
        if expected is None:
            expected = record.get("payload_npz_sha256")
        if isinstance(expected, str) and re.fullmatch("[0-9a-f]{64}", expected):
            origin, target = self.resolve("npz", expected, path)
            self.add(origin, target, expected, key, "numerical-payload")

    def imports(self, key, ref, path, data):
        if not path.endswith(".py"):
            return
        try:
            parsed = ast.parse(data.decode())
        except (SyntaxError, UnicodeDecodeError):
            return
        package = "lab/src/recursive_horizons/" if path.startswith("lab/") else "src/recursive_horizons/"
        scripts = "lab/scripts/" if path.startswith("lab/") else "scripts/"
        for node in ast.walk(parsed):
            names = []
            if isinstance(node, ast.ImportFrom):
                if node.level == 1 and node.module:
                    names = [package + node.module.replace(".", "/") + ".py"]
                elif node.module and node.module.startswith("recursive_horizons."):
                    names = [package + node.module.split(".", 1)[1].replace(".", "/") + ".py"]
                elif node.module == "recursive_horizons":
                    names = [package + item.name + ".py" for item in node.names]
            elif isinstance(node, ast.Import):
                for item in node.names:
                    if item.name.startswith("recursive_horizons."):
                        names.append(package + item.name.split(".", 1)[1].replace(".", "/") + ".py")
                    elif item.name.startswith("derive_nsc_"):
                        names.append(scripts + item.name + ".py")
            for name in names:
                if name in self.index[ref]:
                    self.add(ref, name, parent=key, label="local-import")
                elif name in self.index[self.science]:
                    self.add(self.science, name, parent=key, label="local-import-at-science-pin")

    def build(self):
        seeds = []
        for name in RECORDS:
            path = "lab/results/development/" + name + ".json"
            seeds.append(self.add(self.science, path))
        # All source-cache objects and original pin declarations remain intact.
        for path in self.index[self.science]:
            if path.startswith("lab/.source-history/"):
                self.add(self.science, path, label="preserved-source-cache")
        # The reused finite operator construction is a separate preserved strand.
        for path in self.index[self.science]:
            if path.startswith("paper/nested-quality-evidence/") or path == "docs/nsc-nested-qualities.md":
                self.add(self.science, path, label="finite-operator-background")
        while self.pending:
            key, ref, path, data = self.pending.pop()
            if path.endswith(".json"):
                value = json.loads(data)
                if isinstance(value, dict):
                    self.record_edges(key, ref, path, value)
            self.imports(key, ref, path, data)
        return {"schema": "NSC-FINITE-REGENERATION-EVIDENCE-v1", "science_commit": self.science,
            "historical_source_commit": HISTORICAL,
            "complete_declared_scientific_dependency_graph": True,
            "numerical_history_substitution": False, "simulation_reexecuted": False,
            "scope": "Pinned finite-realization records, numerical payloads and declared producers/local imports; preserved finite-operator background and replay cache; byte authentication is not continuum certification",
            "roots": seeds, "files": [self.nodes[k] for k in sorted(self.nodes)],
            "edges": sorted(self.edges, key=lambda item: (item["from"], item["to"], item["role"]))}


def authenticate_snapshot(snapshot):
    if snapshot.get("complete_declared_scientific_dependency_graph") is not True:
        raise ValueError("incomplete declared scientific graph")
    for item in snapshot["files"]:
        data = read_git(item["commit"], item["path"])
        if digest(data) != item["sha256"] or len(data) != item["bytes"] or blob_id(data) != item["git_blob"]:
            raise ValueError("scientific graph authentication failed: " + item["path"])
    return snapshot


def source_files():
    folder = ROOT / SOURCE
    files = {name: (folder / name).read_bytes() for name in ("main.tex", "references.bib")}
    if (folder / "result_macros.tex").exists():
        files["result_macros.tex"] = (folder / "result_macros.tex").read_bytes()
    main = files["main.tex"].decode("ascii")
    if any(value < 32 and value not in (9, 10) for data in files.values() for value in data):
        raise ValueError("control bytes in manuscript source")
    if re.search(r"\\(?:write18|openout|read|immediate|catcode)\b", main):
        raise ValueError("unsupported external TeX operation")
    figures = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", main)
    if len(figures) != 4 or len(set(figures)) != 4:
        raise ValueError("exactly four declared vector figures are required")
    for name in figures:
        safe_name(name)
        if not re.fullmatch(r"figures/[a-z-]+\.pdf", name):
            raise ValueError("unexpected figure dependency")
        files[name] = (folder / name).read_bytes()
    inputs = re.findall(r"\\input\{([^}]+)\}", main)
    if any(name not in ("result_macros", "result_macros.tex") for name in inputs):
        raise ValueError("unexpected TeX input graph")
    return files


def compile_twice(files):
    if not COMPILER.is_file():
        raise ValueError("pinned Tectonic binary is unavailable")
    env = dict(os.environ, SOURCE_DATE_EPOCH="0", FORCE_SOURCE_DATE="1", TZ="UTC")
    products = []
    version = subprocess.check_output([str(COMPILER), "--version"], text=True).strip()
    for _ in range(2):
        with tempfile.TemporaryDirectory(prefix="nsc-finite-article-") as directory:
            folder = Path(directory)
            for name, data in files.items():
                destination = folder / safe_name(name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(data)
            command = [str(COMPILER), "--untrusted", "--only-cached", "--keep-intermediates", "--keep-logs", "main.tex"]
            result = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True)
            log = (folder / "main.log").read_text(errors="replace") if (folder / "main.log").exists() else ""
            if result.returncode:
                raise ValueError("local TeX build failed:\n" + result.stdout[-3000:] + result.stderr[-2500:] + log[-1500:])
            if re.search(r"undefined references|undefined citations|Citation .* undefined|Reference .* undefined|Rerun to get", log):
                raise ValueError("unresolved TeX references")
            if re.search(r"Overfull \\[hv]box", log):
                raise ValueError("overfull TeX box:\n" + "\n".join(re.findall(r"Overfull[^\n]*\n(?:[^\n]*\n){0,4}", log)))
            pdf, bbl = (folder / "main.pdf").read_bytes(), (folder / "main.bbl").read_bytes()
            products.append((pdf, bbl, log))
    if products[0][:2] != products[1][:2]:
        raise ValueError("two clean PDF/bibliography builds differ")
    cache = Path.home() / "Library/Caches/TectonicProject.Tectonic"
    bundle_records = []
    for path in sorted((cache / "bundles/hashes").glob("*.tar")):
        bundle_records.append({"cache_entry": path.name, "resource_digest": path.read_text().strip(),
                               "entry_sha256": digest(path.read_bytes())})
    return products[0][0], products[0][1], products[0][2], {
        "engine": version, "binary_sha256": digest(COMPILER.read_bytes()), "metadata_epoch": 0,
        "untrusted": True, "only_cached_packages": True, "clean_builds": 2,
        "standard_resource_bundle": bundle_records,
        "pdf_and_bibliography_byte_identical": True, "arxiv_TeX_Live_2025_processor_verified": False}


def packed_source(files):
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode="wb", filename="", mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.USTAR_FORMAT) as tar:
            for name, data in sorted(files.items()):
                safe_name(name)
                if not re.fullmatch(r"(?:[A-Za-z0-9_-]+\.(?:tex|bib|bbl)|figures/[a-z-]+\.pdf)", name):
                    raise ValueError("unexpected archive member: " + name)
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                tar.addfile(info, io.BytesIO(data))
    return output.getvalue()


def metadata(main, pages):
    text = main.decode("ascii")
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", text, re.S).group(1)
    abstract = " ".join(abstract.split())
    title = re.search(r"\\title\{(.*?)\}\s*\\author", text, re.S).group(1)
    title = " ".join(title.replace("\\\\", " ").split())
    authors = re.search(r"\\author\{([^}]+)\}", text).group(1)
    if authors != "Douglas Ek" or not title.isascii() or not abstract.isascii() or len(abstract) > 1920:
        raise ValueError("invalid arXiv title/author/ASCII abstract metadata")
    if not re.search(r"(?i)(significant|substantial).*AI assistance", text, re.S):
        raise ValueError("significant AI assistance disclosure is missing")
    fields = {"title": title, "authors": ["Douglas Ek"], "abstract": abstract,
              "comments": f"{pages} pages, 4 figures. Finite numerical realization; author review pending.",
              "primary_category_proposed": "gr-qc", "cross_lists_proposed": [],
              "license_proposed": "arXiv non-exclusive distribution license",
              "license_url": "https://arxiv.org/licenses/nonexclusive-distrib/1.0/",
              "human_author_review_pending": True, "submission_performed": False}
    check = {"schema": "NSC-FINITE-ARXIV-DRYRUN-v1", "local_checks": {
        "ASCII_metadata": True, "abstract_characters": len(abstract), "abstract_limit": 1920,
        "sole_human_author": True, "AI_assistance_disclosed": True, "metadata_matches_source": True},
        "official_requirements_checked": "2026-10-01", "requirement_sources": {
            "metadata": "https://info.arxiv.org/help/prep.html", "TeX": "https://info.arxiv.org/help/submit_tex.html",
            "processor": "https://info.arxiv.org/help/faq/texlive.html", "license": "https://info.arxiv.org/help/license/index.html"},
        "arxiv_default_documented": "TeX Live 2025", "arxiv_processor_verified": False,
        "local_compiler_only": True, "author_review": "pending", "license_choice": "proposed; human approval pending",
        "submission_performed": False, "account_or_UI_accessed": False}
    return fields, check


def pdf_checks(path):
    info = subprocess.check_output(["pdfinfo", str(path)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)$", info, re.M).group(1))
    font_rows = subprocess.check_output(["pdffonts", str(path)], text=True).splitlines()[2:]
    if not font_rows or any(row.split()[-5] != "yes" for row in font_rows if row.strip()):
        raise ValueError("PDF contains unembedded fonts")
    text = subprocess.check_output(["pdftotext", str(path), "-"], text=True)
    if "Douglas Ek" not in text or len(text.strip()) < 1000:
        raise ValueError("PDF is not searchable or does not contain the author")
    return pages, {"embedded_fonts": True, "searchable_text": True, "pages": pages}


def build(science):
    before = protected_files()
    snapshot = Graph(science).build()
    authenticate_snapshot(snapshot)
    folder = ROOT / SOURCE
    (folder / "evidence").mkdir(parents=True, exist_ok=True)
    (ROOT / EVIDENCE).write_text(json.dumps(snapshot, indent=2) + "\n")
    files = source_files()
    pdf, bbl, log, compiler = compile_twice(files)
    (ROOT / PDF).write_bytes(pdf)
    archive_files = dict(files, **{"main.bbl": bbl})
    archive = packed_source(archive_files)
    if archive != packed_source(archive_files):
        raise ValueError("archive packing is nondeterministic")
    (ROOT / ARCHIVE).write_bytes(archive)
    pages, quality = pdf_checks(ROOT / PDF)
    fields, dryrun = metadata(files["main.tex"], pages)
    (ROOT / METADATA).write_text(json.dumps(fields, indent=2) + "\n")
    (ROOT / DRYRUN).write_text(json.dumps(dryrun, indent=2) + "\n")
    inputs = {(SOURCE / name).as_posix(): digest(data) for name, data in files.items()}
    for path in (EVIDENCE, METADATA, DRYRUN, SOURCE / "figure-manifest.json"):
        inputs[path.as_posix()] = digest((ROOT / path).read_bytes())
    inputs["scripts/finite_regeneration.py"] = digest(Path(__file__).read_bytes())
    inputs["scripts/finite_regeneration_figures.py"] = digest((ROOT / "scripts/finite_regeneration_figures.py").read_bytes())
    manifest = {"schema": "NSC-FINITE-REGENERATION-MANIFEST-v1", "status": "MEASURED_FINITE_REALIZATION",
        "submission_ready": False, "author_review": "pending", "science_commit": science,
        "title": fields["title"], "authors": fields["authors"], "pages": pages, "inputs": inputs,
        "pdf_sha256": digest(pdf), "source_archive_sha256": digest(archive), "bbl_sha256": digest(bbl),
        "archive_members": sorted(archive_files), "compiler": compiler, "pdf_quality": quality,
        "complete_declared_scientific_dependency_graph": True, "arxiv_processor_verified": False,
        "historical_artifacts_unchanged": before == protected_files(), "simulation_reexecuted": False}
    if not manifest["historical_artifacts_unchanged"] or source_files() != files:
        raise ValueError("historical artifact or frozen article changed during build")
    (ROOT / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    log_folder = ROOT / ".build/finite-regeneration"
    log_folder.mkdir(parents=True, exist_ok=True)
    (log_folder / "main.log").write_text(log)
    return manifest


def verify():
    manifest = json.loads((ROOT / MANIFEST).read_text())
    if manifest.get("status") != "MEASURED_FINITE_REALIZATION" or manifest.get("submission_ready") is not False:
        raise ValueError("invalid finite publication scope")
    snapshot = authenticate_snapshot(json.loads((ROOT / EVIDENCE).read_text()))
    if snapshot["science_commit"] != manifest["science_commit"]:
        raise ValueError("science pin differs")
    for name, expected in manifest["inputs"].items():
        safe_name(name)
        if digest((ROOT / name).read_bytes()) != expected:
            raise ValueError("publication input differs: " + name)
    if digest((ROOT / PDF).read_bytes()) != manifest["pdf_sha256"]:
        raise ValueError("new PDF hash differs")
    archive = (ROOT / ARCHIVE).read_bytes()
    if digest(archive) != manifest["source_archive_sha256"]:
        raise ValueError("source archive hash differs")
    files = source_files()
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
        members = tar.getmembers()
        if sorted(item.name for item in members) != manifest["archive_members"] or any(not item.isfile() for item in members):
            raise ValueError("archive inventory differs")
        for item in members:
            content = tar.extractfile(item).read()
            expected = manifest["bbl_sha256"] if item.name == "main.bbl" else digest(files[item.name])
            if digest(content) != expected or item.mtime != 0 or item.uid != 0 or item.gid != 0:
                raise ValueError("archive content or deterministic metadata differs")
    fields, dryrun = metadata(files["main.tex"], manifest["pages"])
    if fields != json.loads((ROOT / METADATA).read_text()) or dryrun != json.loads((ROOT / DRYRUN).read_text()):
        raise ValueError("metadata dry run differs from source")
    return {"status": "MEASURED_FINITE_REALIZATION", "artifact_integrity": "PASS", "pages": manifest["pages"],
            "science_commit": manifest["science_commit"], "pdf_sha256": manifest["pdf_sha256"],
            "submission_ready": False, "arxiv_processor_verified": False}
