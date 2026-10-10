"""Guard the inventory's provenance and scope boundaries, not the cosmology."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


SPEC = importlib.util.spec_from_file_location(
    "research_inventory", Path(__file__).resolve().parents[1] / "scripts/research_inventory.py")
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


def git(root, *args):
    return subprocess.check_output(["git", "-c", "user.name=Inventory test",
                                   "-c", "user.email=inventory@example.invalid", *args], cwd=root)


def preserved_fixture(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q")
    (root / "docs").mkdir()
    (root / "docs/owner.md").write_text("# Original owner\n")
    (root / "docs/deleted.md").write_text("# Retired\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "fixture")
    head = git(root, "rev-parse", "HEAD").decode().strip()
    saved = tmp_path / "saved"
    (saved / "files/docs").mkdir(parents=True)
    data = b"# Preserved dirty owner\n"
    (saved / "files/docs/owner.md").write_bytes(data)
    receipt = {"repository": "fixture", "head": head, "branch": "fixture",
               "files": [{"path": "docs/owner.md", "bytes": len(data), "sha256": inventory.digest(data)}]}
    (saved / "receipt.json").write_bytes(inventory.encoded(receipt))
    (saved / "status.bin").write_bytes(b" M docs/owner.md\0 D docs/deleted.md\0")
    return root, saved, data


def test_baseline_uses_preservation_not_live_files_and_honours_deletion(tmp_path):
    root, saved, data = preserved_fixture(tmp_path)
    (root / "docs/owner.md").write_text("# Concurrent implementation\n")
    (root / "docs/new-work.md").write_text("# Not in baseline\n")
    rows, snapshot = inventory.baseline(root, saved)
    assert rows["docs/owner.md"]["sha256"] == inventory.digest(data)
    assert rows["docs/owner.md"]["origin"] == "modified"
    assert "docs/deleted.md" not in rows
    assert "docs/new-work.md" not in rows
    assert snapshot["deleted_at_baseline"] == ["docs/deleted.md"]


def test_corrupt_preserved_input_refused(tmp_path):
    root, saved, _ = preserved_fixture(tmp_path)
    (saved / "files/docs/owner.md").write_text("corrupt")
    with pytest.raises(ValueError, match="preserved input has drifted"):
        inventory.baseline(root, saved)


def test_git_object_reader_ignores_worktree_drift(tmp_path):
    root, saved, _ = preserved_fixture(tmp_path)
    head = json.loads((saved / "receipt.json").read_bytes())["head"]
    oid = git(root, "rev-parse", head + ":docs/owner.md").decode().strip()
    expected = git(root, "show", head + ":docs/owner.md")
    (root / "docs/owner.md").write_text("different")
    objects = inventory.GitObjects(root)
    try:
        assert objects.read(oid) == expected
        assert expected != (root / "docs/owner.md").read_bytes()
    finally:
        objects.close()


def test_document_identity_is_not_a_scientific_verdict():
    doc = inventory.extract_document("lab/docs/test.md", b"# Example\nStatus: PROVED\nNot a global horizon.\n")
    assert doc["full_text_scanned"] is True
    assert doc["equations_independently_verified"] is False
    assert doc["scope_and_status_excerpts"][0]["line"] == 2


def test_links_keep_external_historical_missing_and_resolved_distinct():
    paths = {"lab/results/example.json", "lab/docs/owner.md"}
    source = "lab/docs/owner.md"
    assert inventory.resolve_link(source, "../results/example.json", paths)["status"] == "resolved"
    assert inventory.resolve_link(source, "../results/missing.json", paths)["status"] == "missing_in_preserved_baseline"
    assert inventory.resolve_link(source, "/old/mac/evidence.json", paths)["status"] == "historical_or_host_absolute"
    assert inventory.resolve_link(source, "https://arxiv.org/abs/1234.56789", paths)["status"] == "external_reference"


def test_family_grouping_preserves_distinct_child_neck_calculations():
    assert inventory.result_family("lab/results/development/nsc-example-v2/points/12.json") == "lab/results/development/nsc-example"
    assert inventory.result_family("lab/results/development/nsc-example-v1.json") == "lab/results/development/nsc-example"
    assert inventory.result_family("lab/results/development/child-neck-v1/israel_matching.json") != inventory.result_family("lab/results/development/child-neck-v1/neck_flux.json")
    assert inventory.result_family("lab/archive/old/results/development/nsc-example-v1.json") == "lab/archive/old/results/development/nsc-example"
    assert inventory.result_family("lab/results/development/artifacts/nsc-example-v2/point-1.npz") == "lab/results/development/artifacts/nsc-example"


def test_inventory_cannot_hash_itself_or_runtime_cache():
    assert inventory.excluded("lab/results/research-reset/inventory-summary.json")
    assert inventory.excluded("lab/pytest-of-swedi/pytest-1/result.json")
    assert not inventory.excluded("lab/.source-history/manifest.json")
