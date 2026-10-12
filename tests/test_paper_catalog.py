"""One current manuscript and authenticated historical editions."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("publication_catalog_check", ROOT / "scripts/check_publication.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@pytest.fixture
def catalog_tree(tmp_path, monkeypatch):
    def put(path, data):
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data if isinstance(data, bytes) else json.dumps(data).encode())
        return hashlib.sha256(target.read_bytes()).hexdigest()

    main = {"id": "main", "role": "main_manuscript", "active": True,
            "version": "working-2026-10-12", "date": "2026-10-12",
            "status": "WORKING_DRAFT", "submission_ready": False, "author_review": "pending",
            "source": "paper/nested-space-cosmology.md", "pdf": "paper/nested-space-cosmology.pdf",
            "metadata": "paper/metadata.json", "manifest": "paper/build-manifest.json"}
    source_hash = put(main["source"], b"current scientific argument")
    pdf_hash = put(main["pdf"], b"current rendered paper")
    metadata_hash = put(main["metadata"], {"version": main["version"], "date_iso": main["date"]})
    put(main["manifest"], {"source": main["source"], "output": main["pdf"], "metadata": main["metadata"],
                          "version": main["version"], "source_sha256": source_hash,
                          "pdf_sha256": pdf_hash, "metadata_sha256": metadata_hash})
    prefix = "paper/archive/foundation-v0.26.0/"
    files = []
    old = {}
    for original, content in [
        ("paper/nested-space-cosmology.md", b"old argument"),
        ("paper/nested-space-cosmology.pdf", b"old rendered paper"),
        ("paper/metadata.json", {"version": "0.26.0"}),
        ("paper/claim-bindings.json", {"claims": {}}),
        ("scripts/build_paper.py", b"old renderer"),
    ]:
        digest = put(prefix + original, content)
        old[original] = digest
        files.append({"original_path": original, "storage_path": prefix + original, "sha256": digest})
    old_build = {"source": "paper/nested-space-cosmology.md", "output": "paper/nested-space-cosmology.pdf",
                 "metadata": "paper/metadata.json", "source_sha256": old["paper/nested-space-cosmology.md"],
                 "pdf_sha256": old["paper/nested-space-cosmology.pdf"], "metadata_sha256": old["paper/metadata.json"],
                 "scientific_build_inputs": [{"path": "scripts/build_paper.py", "sha256": old["scripts/build_paper.py"]}]}
    digest = put(prefix + "paper/build-manifest.json", old_build)
    files.append({"original_path": "paper/build-manifest.json",
                  "storage_path": prefix + "paper/build-manifest.json", "sha256": digest})
    receipt = prefix + "archive-manifest.json"
    put(receipt, {"schema": "NSC-PAPER-ARCHIVE-v1", "id": "foundation", "version": "0.26.0",
                  "source_commit": "1" * 40, "files": files})
    foundation = {"id": "foundation", "role": "archived_foundation", "active": False,
                  "version": "0.26.0", "archive_manifest": receipt,
                  "source": prefix + "paper/nested-space-cosmology.md",
                  "pdf": prefix + "paper/nested-space-cosmology.pdf",
                  "metadata": prefix + "paper/metadata.json", "manifest": prefix + "paper/build-manifest.json"}

    def companion(name):
        source, pdf, manifest, archive = [f"paper/{name}{suffix}" for suffix in (".tex", ".pdf", ".json", ".tar.gz")]
        source_hash, pdf_hash, archive_hash = put(source, b"source"), put(pdf, b"pdf"), put(archive, b"archive")
        put(manifest, {"inputs": {source: source_hash}, "pdf_sha256": pdf_hash, "source_archive_sha256": archive_hash})
        return {"source": source, "pdf": pdf, "manifest": manifest, "source_archive": archive}

    finite = dict(companion("finite"), id="local-gate", role="archived_supporting_result", active=False,
                  status="MEASURED_FINITE_REALIZATION", submission_ready=False,
                  historical_editions=[companion("gate")])
    catalog = {"schema": "NSC-PAPER-CATALOG-v3", "current_document": "main", "current_status": "WORKING_DRAFT",
               "documents": [main, foundation, finite], "finite_result_status": "MEASURED_FINITE_REALIZATION",
               "physical_local_gate_status": "OPEN"}
    put("paper/catalog.json", catalog)
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    return tmp_path, catalog, put


def errors():
    found = []
    checker.check_paper_catalog(found)
    return found


def test_one_main_and_historical_evidence_authenticate(catalog_tree):
    assert errors() == []


@pytest.mark.parametrize("field", ["source", "pdf", "metadata"])
def test_changed_current_content_requires_rebuild(catalog_tree, field):
    root, catalog, _ = catalog_tree
    (root / catalog["documents"][0][field]).write_bytes(b"changed")
    assert errors()


def test_archived_original_cannot_resolve_to_current_pdf(catalog_tree):
    _, catalog, put = catalog_tree
    catalog["documents"][1]["pdf"] = catalog["documents"][0]["pdf"]
    put("paper/catalog.json", catalog)
    assert errors()


def test_archived_source_tampering_is_rejected(catalog_tree):
    root, catalog, _ = catalog_tree
    (root / catalog["documents"][1]["source"]).write_bytes(b"changed history")
    assert errors()


def test_second_active_paper_is_rejected(catalog_tree):
    _, catalog, put = catalog_tree
    catalog["documents"][2]["active"] = True
    put("paper/catalog.json", catalog)
    assert errors()


def test_old_finite_result_does_not_close_physical_application(catalog_tree):
    _, catalog, put = catalog_tree
    catalog["physical_local_gate_status"] = "CLOSED"
    put("paper/catalog.json", catalog)
    assert errors()


def test_historical_companion_cannot_disappear(catalog_tree):
    _, catalog, put = catalog_tree
    catalog["documents"][2]["historical_editions"] = []
    put("paper/catalog.json", catalog)
    assert errors()


def test_unsubmitted_working_draft_cannot_claim_submission(catalog_tree):
    _, catalog, put = catalog_tree
    catalog["documents"][0]["submission_ready"] = True
    put("paper/catalog.json", catalog)
    assert errors()
