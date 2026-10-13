"""Historical link contexts do not weaken live-document or secret checks."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("history_boundary_check", ROOT / "scripts/check_publication.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def prepare(tmp_path, monkeypatch, name, data, *, legacy=False):
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    if legacy:
        manifest = tmp_path / "results/reassembly/manifest.json"
        manifest.parent.mkdir(parents=True)
        manifest.write_text(json.dumps({"rows": [{"origin": "lab", "destination": name,
                                                  "sha256": digest}]}))
    else:
        manifest = tmp_path / "docs/lab-snapshot.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps({
            "schema": "NSC-CONSOLIDATED-LAB-v2", "canonical_root": ".",
            "provenance": {"archived_import_and_reassembly_doi": checker.ARCHIVE_DOI},
            "files": [{"path": name, "sha256": digest, "bytes": len(data),
                       "import_sha256": digest, "import_bytes": len(data)}],
        }))
        release = tmp_path / "results/release-spec.json"
        release.parent.mkdir(parents=True, exist_ok=True)
        release.write_text("{}")
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    monkeypatch.setattr(checker, "REQUIRED", set())
    return path


def findings(path):
    found = []
    checker.check_public_boundary(found, [path])
    return found


def test_exact_old_note_keeps_original_logical_link_context(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"[Old record](../archived-input/record.json)")
    assert not (tmp_path / "results/reassembly/manifest.json").exists()
    assert findings(path) == []
    path.write_bytes(path.read_bytes() + b"\nchanged")
    assert any("broken Markdown link" in item for item in findings(path))


def test_live_entrypoint_never_inherits_historical_exemption(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "README.md", b"[Broken route](lab/no-longer-here.md)")
    assert any("broken Markdown link" in item for item in findings(path))


def test_exact_historical_bytes_still_receive_secret_scan(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"BEGIN RSA " + b"PRIVATE KEY")
    assert any("possible private key" in item for item in findings(path))


def test_rebound_current_receipt_does_not_exempt_changed_note(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"[Old record](missing.json)")
    path.write_bytes(path.read_bytes() + b"\nchanged")
    manifest = tmp_path / "docs/lab-snapshot.json"
    value = json.loads(manifest.read_text())
    value["files"][0].update(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             bytes=path.stat().st_size)
    manifest.write_text(json.dumps(value))
    assert any("broken Markdown link" in item for item in findings(path))


def test_v2_receipt_requires_exact_archive_identity(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"historical note")
    manifest = tmp_path / "docs/lab-snapshot.json"
    value = json.loads(manifest.read_text())
    value["provenance"]["archived_import_and_reassembly_doi"] = "wrong archive"
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="retained-file receipt identity"):
        findings(path)


def test_isolated_legacy_fixture_retains_original_context(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md",
                   b"[Old record](../archived-input/record.json)", legacy=True)
    assert findings(path) == []


def test_extracted_latex_equation_number_is_not_a_link(tmp_path, monkeypatch):
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    found = []
    checker.check_markdown_links(found, tmp_path / "source.md", r"\mathcal{S}\subseteq[n](14)" + "\n")
    assert found == []
    checker.check_markdown_links(found, tmp_path / "source.md", "[An actual numbered link](14)\n")
    assert found == ["broken Markdown link in source.md: 14"]


def restored_capture(tmp_path, monkeypatch, *, name=None, data=None):
    """A reviewed immutable capture has its own original identity, not a blanket scope."""
    name = name or "results/development/reviewed-control/record.json"
    data = data or json.dumps({"source_hashes": {
        "/Users/" + "admin/original/producer.py": "a" * 64,
    }}).encode()
    path = prepare(tmp_path, monkeypatch, name, data)
    receipt = tmp_path / "docs/lab-snapshot.json"
    value = json.loads(receipt.read_text())
    row = value["files"][0]
    row.pop("import_sha256")
    row.pop("import_bytes")
    row["restored_from_commit"] = "1" * 40
    receipt.write_text(json.dumps(value))
    identity = {"source_commit": row["restored_from_commit"],
                "sha256": row["sha256"], "bytes": row["bytes"]}
    monkeypatch.setattr(checker, "RESTORED_CREATION_METADATA", {name: identity})
    monkeypatch.setattr(checker, "_git_bytes", lambda commit, relative: data)
    return path, receipt


def test_reviewed_restored_creator_metadata_keeps_original_path_labels(tmp_path, monkeypatch):
    path, _ = restored_capture(tmp_path, monkeypatch)
    assert findings(path) == []
    assert path.relative_to(tmp_path).as_posix() in checker.historical_content()


def test_modified_restored_bytes_cannot_refresh_away_original_identity(tmp_path, monkeypatch):
    path, receipt = restored_capture(tmp_path, monkeypatch)
    path.write_bytes(path.read_bytes() + b"\nchanged")
    value = json.loads(receipt.read_text())
    value["files"][0].update(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             bytes=path.stat().st_size)
    receipt.write_text(json.dumps(value))
    assert any("private absolute user path" in item for item in findings(path))


@pytest.mark.parametrize("context", ["HEAD", "../outside", "2" * 40])
def test_restored_creator_metadata_rejects_unreviewed_source_ref(tmp_path, monkeypatch, context):
    path, receipt = restored_capture(tmp_path, monkeypatch)
    value = json.loads(receipt.read_text())
    value["files"][0]["restored_from_commit"] = context
    receipt.write_text(json.dumps(value))
    monkeypatch.setattr(checker, "_git_bytes", lambda *_: pytest.fail("unsafe ref reached Git"))
    assert any("private absolute user path" in item for item in findings(path))


def test_restored_creator_metadata_requires_original_blob_identity(tmp_path, monkeypatch):
    _, _ = restored_capture(tmp_path, monkeypatch)
    monkeypatch.setattr(checker, "_git_bytes", lambda *_: b"changed original bytes")
    with pytest.raises(ValueError, match="creator metadata differs from original bytes"):
        checker.historical_content()


def test_reviewed_restored_metadata_still_receives_secret_scan(tmp_path, monkeypatch):
    data = (b"/Users/" + b"admin/original/path\nBEGIN RSA " + b"PRIVATE KEY")
    path, _ = restored_capture(tmp_path, monkeypatch, data=data)
    assert any("possible private key" in item for item in findings(path))


def test_current_entrypoint_cannot_acquire_restored_capture_scope(tmp_path, monkeypatch):
    path, _ = restored_capture(tmp_path, monkeypatch, name="README.md")
    assert any("private absolute user path" in item for item in findings(path))


def test_actual_reviewed_capture_matches_its_local_and_original_identities():
    snapshot, receipts = checker.retained_receipts()
    assert snapshot is not None
    expected = {relative: identity["sha256"]
                for relative, identity in checker.RESTORED_CREATION_METADATA.items()}
    assert checker.restored_creation_metadata(receipts) == expected
