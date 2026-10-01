"""Focused-publication byte, archive, source-graph and metadata contracts."""
from copy import deepcopy
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import finite_regeneration as article


def test_finite_publication_authenticates_exact_declared_inputs():
    result = article.verify()
    assert result["artifact_integrity"] == "PASS"
    assert result["status"] == "MEASURED_FINITE_REALIZATION"
    assert result["submission_ready"] is False
    assert result["arxiv_processor_verified"] is False
    snapshot = json.loads((article.ROOT / article.EVIDENCE).read_text())
    assert snapshot["numerical_history_substitution"] is False
    assert snapshot["simulation_reexecuted"] is False
    assert all(row["commit"] == snapshot["science_commit"] for row in snapshot["files"] if row["path"].endswith((".json", ".npz")))
    assert sum(row["path"].startswith("lab/.source-history/objects/") for row in snapshot["files"]) == 25


def test_snapshot_rejects_mutated_scientific_content_or_blob_identity():
    snapshot = json.loads((article.ROOT / article.EVIDENCE).read_text())
    corrupted = deepcopy(snapshot)
    corrupted["files"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="authentication failed"):
        article.authenticate_snapshot(corrupted)
    corrupted = deepcopy(snapshot)
    corrupted["files"][0]["git_blob"] = "0" * 40
    with pytest.raises(ValueError, match="authentication failed"):
        article.authenticate_snapshot(corrupted)


def test_source_archive_is_deterministic_and_only_contains_used_tex_assets():
    files = article.source_files()
    files["main.bbl"] = b"deterministic test bibliography\n"
    first = article.packed_source(files)
    assert article.packed_source(dict(reversed(list(files.items())))) == first
    with tarfile.open(fileobj=io.BytesIO(first), mode="r:gz") as tar:
        members = tar.getmembers()
        assert {item.name for item in members} == set(files)
        assert all(item.isfile() and item.mtime == 0 and item.uid == 0 and item.gid == 0 for item in members)
        assert all(not item.name.endswith((".npz", ".json", ".aux", ".log")) for item in members)
    with pytest.raises(ValueError, match="unexpected archive member"):
        article.packed_source({"trajectory.npz": b"never duplicate scientific arrays"})
    with pytest.raises(ValueError, match="unsafe artifact path"):
        article.packed_source({"../main.tex": b"escape"})


def test_metadata_is_ascii_short_single_human_author_and_pending_review():
    manifest = json.loads((article.ROOT / article.MANIFEST).read_text())
    fields, checks = article.metadata(article.source_files()["main.tex"], manifest["pages"])
    assert fields["authors"] == ["Douglas Ek"]
    assert fields["abstract"].isascii()
    assert len(fields["abstract"]) <= 1920
    assert fields["human_author_review_pending"] is True
    assert fields["submission_performed"] is False
    assert checks["local_checks"]["AI_assistance_disclosed"] is True
    assert checks["arxiv_processor_verified"] is False
    main = article.source_files()["main.tex"]
    altered = main.replace(b"\\author{Douglas Ek}", b"\\author{Douglas Ek and AI}")
    with pytest.raises(ValueError, match="metadata"):
        article.metadata(altered, manifest["pages"])


def test_compiler_contract_is_local_cached_untrusted_and_reproducible():
    manifest = json.loads((article.ROOT / article.MANIFEST).read_text())
    compiler = manifest["compiler"]
    assert compiler["engine"] == "Tectonic 0.17.0"
    assert compiler["only_cached_packages"] is True
    assert compiler["untrusted"] is True
    assert compiler["clean_builds"] == 2
    assert compiler["pdf_and_bibliography_byte_identical"] is True
    assert compiler["metadata_epoch"] == 0
    assert compiler["arxiv_TeX_Live_2025_processor_verified"] is False
    assert compiler["standard_resource_bundle"]
    assert manifest["historical_artifacts_unchanged"] is True
