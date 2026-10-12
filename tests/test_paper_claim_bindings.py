"""Claims in the main manuscript must name their actual finite evidence."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("main_claim_check", ROOT / "scripts/check_claim_bindings.py")
claims = importlib.util.module_from_spec(spec)
spec.loader.exec_module(claims)


def fixture(tmp_path, *, link=True):
    (tmp_path / "paper").mkdir()
    (tmp_path / "results").mkdir()
    (tmp_path / "results/wave.json").write_text('{"domain": "finite"}')
    source = "## Observer response\n\n<!-- nsc-claim:atlas-field-response -->\nFinite conditional response.\n"
    if link:
        source += "[Record](../results/wave.json)\n"
    (tmp_path / "paper/nested-space-cosmology.md").write_text(source)
    (tmp_path / "paper/claim-bindings.json").write_text(json.dumps({
        "paper": "paper/nested-space-cosmology.md", "status_values": {"measured": "finite result"},
        "claims": {"atlas-field-response": {"status": "measured"}}, "findings": []}))


def test_current_measured_claim_resolves_direct_record(tmp_path):
    fixture(tmp_path)
    assert claims.main(["--root", str(tmp_path)]) == 0


def test_reported_result_without_evidence_link_fails(tmp_path):
    fixture(tmp_path, link=False)
    assert claims.main(["--root", str(tmp_path)]) != 0


def test_missing_evidence_is_not_a_valid_claim(tmp_path):
    fixture(tmp_path)
    (tmp_path / "results/wave.json").unlink()
    assert claims.main(["--root", str(tmp_path)]) != 0


def test_duplicate_marker_is_not_accepted(tmp_path):
    fixture(tmp_path)
    path = tmp_path / "paper/nested-space-cosmology.md"
    path.write_text(path.read_text() + "\n<!-- nsc-claim:atlas-field-response -->\n")
    with pytest.raises(ValueError, match="duplicate"):
        claims.main(["--root", str(tmp_path)])
