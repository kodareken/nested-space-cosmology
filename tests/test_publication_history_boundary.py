"""Historical link contexts do not weaken live-document or secret checks."""
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("history_boundary_check", ROOT / "scripts/check_publication.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def prepare(tmp_path, monkeypatch, name, data):
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    manifest = tmp_path / "results/reassembly/manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({"rows": [{"origin": "lab", "destination": name,
                                              "sha256": hashlib.sha256(data).hexdigest()}]}))
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    monkeypatch.setattr(checker, "REQUIRED", set())
    return path


def findings(path):
    found = []
    checker.check_public_boundary(found, [path])
    return found


def test_exact_old_note_keeps_original_logical_link_context(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"[Old record](../archived-input/record.json)")
    assert findings(path) == []
    path.write_bytes(path.read_bytes() + b"\nchanged")
    assert any("broken Markdown link" in item for item in findings(path))


def test_live_entrypoint_never_inherits_historical_exemption(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "README.md", b"[Broken route](lab/no-longer-here.md)")
    assert any("broken Markdown link" in item for item in findings(path))


def test_exact_historical_bytes_still_receive_secret_scan(tmp_path, monkeypatch):
    path = prepare(tmp_path, monkeypatch, "docs/old-note.md", b"BEGIN RSA " + b"PRIVATE KEY")
    assert any("possible private key" in item for item in findings(path))


def test_extracted_latex_equation_number_is_not_a_link(tmp_path, monkeypatch):
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    found = []
    checker.check_markdown_links(found, tmp_path / "source.md", r"\mathcal{S}\subseteq[n](14)" + "\n")
    assert found == []
    checker.check_markdown_links(found, tmp_path / "source.md", "[An actual numbered link](14)\n")
    assert found == ["broken Markdown link in source.md: 14"]
