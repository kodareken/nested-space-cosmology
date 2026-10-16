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


@pytest.mark.parametrize("expression", [
    r"$F(T+\operatorname{atanh}r)$",
    "$$\n" + r"\operatorname*{arg\,min}_x f(x)" + "\n$$",
])
def test_mathjax_operator_names_are_valid_in_live_markdown(tmp_path, monkeypatch, expression):
    path = prepare(tmp_path, monkeypatch, "README.md", expression.encode())
    assert findings(path) == []


def test_markdown_delimiter_guard_remains_active(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "README.md", br"\(x=1\)")
    assert findings(path) == ["unsupported GitHub math parenthesized inline delimiter in README.md"]


def test_result_inventory_reports_unregistered_and_missing_evidence():
    found = {"results/kept.json", "results/unregistered.json"}
    authenticated = {"results/kept.json", "results/missing.json"}
    errors = []
    checker.check_result_inventory(errors, found, authenticated)
    assert errors == [
        "tracked NSC result files differ from authenticated local/archived closure"
        "; unregistered: ['results/unregistered.json']; absent: ['results/missing.json']"
    ]
    errors = []
    checker.check_result_inventory(errors, found, found)
    assert errors == []


def test_portable_result_summaries_require_the_same_inventory_receipts(tmp_path, monkeypatch):
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    records = {"results/nsc-release.json", "results/development/audit.json",
               "results/summaries/audit.json"}
    for name in (*records, ".build/results/raw.json", "results/development/nested/record.json"):
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")
    assert checker.result_inventory_files() == records
    errors = []
    checker.check_result_inventory(errors, records, records-{"results/summaries/audit.json"})
    assert len(errors) == 1 and "unregistered: ['results/summaries/audit.json']" in errors[0]


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


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_display_equation_arguments_are_not_links_but_adjacent_links_are(tmp_path, monkeypatch, newline):
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    (tmp_path / "present.md").write_text("present")
    text = newline.join(("[Existing note](present.md)", "$$",
                         r"\int [K_0(z+g_1+B)+G](q,V)\,dq", "$$",
                         "[Missing note](missing.md)", ""))
    found = []
    checker.check_markdown_links(found, tmp_path / "source.md", text)
    assert found == ["broken Markdown link in source.md: missing.md"]


def test_unclosed_math_opener_does_not_hide_a_document_link(tmp_path, monkeypatch):
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    found = []
    checker.check_markdown_links(found, tmp_path / "source.md",
                                 "$$\nx=1\n[Missing note](missing.md)\n")
    assert found == ["broken Markdown link in source.md: missing.md"]


def current_record(tmp_path, monkeypatch, *, name=None, data=None):
    """Current normalized records receive current checks, including after restoration."""
    name = name or "results/development/nsc-discovery-vacuum-control-v1/record.json"
    data = data or json.dumps({"source_hashes": {"scripts/producer.py": "a" * 64}}).encode()
    path = prepare(tmp_path, monkeypatch, name, data)
    receipt = tmp_path / "docs/lab-snapshot.json"
    value = json.loads(receipt.read_text())
    row = value["files"][0]
    row.pop("import_sha256")
    row.pop("import_bytes")
    row["restored_from_commit"] = "4390d23901b666659c7c419ba3fdf8c9a08a9794"
    receipt.write_text(json.dumps(value))
    return path, receipt


def test_normalized_relative_record_is_current_and_accepted(tmp_path, monkeypatch):
    path, _ = current_record(tmp_path, monkeypatch)
    assert findings(path) == []
    assert path.relative_to(tmp_path).as_posix() not in checker.historical_content()


@pytest.mark.parametrize("context", ["HEAD", "../outside", "4390d23901b666659c7c419ba3fdf8c9a08a9794"])
def test_restoration_and_refreshed_receipt_never_exempt_private_paths(tmp_path, monkeypatch, context):
    path, receipt = current_record(tmp_path, monkeypatch)
    path.write_bytes(json.dumps({"source_hashes": {
        "/Users/" + "admin/original/producer.py": "a" * 64,
    }}).encode())
    value = json.loads(receipt.read_text())
    value["files"][0].update(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             bytes=path.stat().st_size, restored_from_commit=context)
    receipt.write_text(json.dumps(value))
    monkeypatch.setattr(checker, "_git_bytes", lambda *_: pytest.fail("current record read from Git"))
    assert any("private absolute user path" in item for item in findings(path))


def test_normalized_record_still_receives_secret_scan(tmp_path, monkeypatch):
    data = json.dumps({"source_hashes": {"scripts/producer.py": "a" * 64},
                       "diagnostic": "BEGIN RSA " + "PRIVATE KEY"}).encode()
    path, _ = current_record(tmp_path, monkeypatch, data=data)
    assert any("possible private key" in item for item in findings(path))


def test_current_entrypoint_has_no_restoration_exception(tmp_path, monkeypatch):
    path, _ = current_record(tmp_path, monkeypatch, name="README.md",
                             data=b"/Users/" + b"admin/original/path")
    assert any("private absolute user path" in item for item in findings(path))


def test_actual_normalized_vacuum_record_receives_current_checks(monkeypatch):
    relative = "results/development/nsc-discovery-vacuum-control-v1/record.json"
    path = ROOT / relative
    text = path.read_text(encoding="utf-8")
    assert "/Users/" + "admin" not in text
    record = json.loads(text)
    for name in (*record["source_hashes"], *record["input_hashes"]):
        assert not Path(name).is_absolute() and not Path(name).drive
        assert ".." not in Path(name).parts
    monkeypatch.setattr(checker, "REQUIRED", set())
    assert relative not in checker.historical_content()
    assert findings(path) == []
