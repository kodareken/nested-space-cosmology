"""Original-source authentication never substitutes missing numerical data."""

import hashlib
from pathlib import Path

import pytest

from recursive_horizons.nsc_frozen_sources import (
    authenticate_source, frozen_source_context, selected_source_sha256,
    is_implementation_source, selected_source_context, source_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
COMMIT = "329ef809e466a86ed4a290ee32c19485f85add9d"
SOURCE = "scripts/derive_nsc_evolved_incoming_state.py"
EXPECTED = "f33ee6bb7854feb424aefcd5063105f6c1222ff082a570c70ff36b841d411939"


def test_original_selected_source_matches_saved_producer():
    raw = authenticate_source(ROOT, SOURCE, EXPECTED, context=COMMIT)
    assert hashlib.sha256(raw).hexdigest() == EXPECTED
    assert source_sha256(ROOT, SOURCE, context=COMMIT) == EXPECTED
    with pytest.raises(RuntimeError, match="frozen source drift"):
        authenticate_source(ROOT, SOURCE, "0" * 64, context=COMMIT)


def test_current_authentication_rejects_missing_and_changed_source(tmp_path):
    relative = "scripts/owner.py"
    with pytest.raises(FileNotFoundError, match="current source unavailable"):
        authenticate_source(tmp_path, relative, EXPECTED)
    path = tmp_path / relative
    path.parent.mkdir()
    path.write_bytes(b"current source\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert authenticate_source(tmp_path, relative, digest) == path.read_bytes()
    path.write_bytes(b"changed current source\n")
    with pytest.raises(RuntimeError, match="current source drift"):
        authenticate_source(tmp_path, relative, digest)


@pytest.mark.parametrize("relative", [
    "../scripts/source.py", "/scripts/source.py", "C:/scripts/source.py",
    r"scripts\source.py", "scripts/../source.py", "scripts//source.py", ".",
])
def test_unsafe_paths_are_rejected_before_lookup(tmp_path, relative):
    with pytest.raises(ValueError, match="repository-relative"):
        authenticate_source(tmp_path, relative, EXPECTED, context=COMMIT)


def test_full_immutable_context_is_required(tmp_path):
    with pytest.raises(ValueError, match="full immutable source commit"):
        source_sha256(tmp_path, SOURCE, context="329ef809")
    with pytest.raises(ValueError, match="full immutable source commit"):
        with frozen_source_context("main"):
            pytest.fail("an unpinned context was selected")


def test_selected_context_restores_on_nested_and_exceptional_exit():
    assert selected_source_context() is None
    with frozen_source_context(COMMIT):
        assert selected_source_context() == COMMIT
        assert selected_source_sha256(ROOT, SOURCE) == EXPECTED
        with pytest.raises(RuntimeError, match="frozen source unavailable"):
            with frozen_source_context("0" * 40):
                assert selected_source_context() == "0" * 40
                selected_source_sha256(ROOT, SOURCE)
        assert selected_source_sha256(ROOT, SOURCE) == EXPECTED
        assert selected_source_context() == COMMIT
    assert selected_source_context() is None
    assert selected_source_sha256(ROOT, "pyproject.toml") == source_sha256(ROOT, "pyproject.toml")


def test_default_current_lookup_ignores_an_explicit_scope(tmp_path):
    with frozen_source_context(COMMIT):
        with pytest.raises(FileNotFoundError, match="current source unavailable"):
            source_sha256(tmp_path, SOURCE, context=None)
        with pytest.raises(FileNotFoundError, match="current source unavailable"):
            authenticate_source(tmp_path, SOURCE, EXPECTED, context=None)
    with pytest.raises(FileNotFoundError, match="current source unavailable"):
        selected_source_sha256(tmp_path, SOURCE)


def test_numerical_outputs_stay_local_inside_a_source_scope(tmp_path):
    relative = "results/development/retained.json"
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_bytes(b'{"value": 7}\n')
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with frozen_source_context(COMMIT):
        assert selected_source_sha256(tmp_path, relative) == digest
        with pytest.raises(ValueError, match="implementation source"):
            authenticate_source(tmp_path, relative, digest, context=COMMIT)
        path.unlink()
        with pytest.raises(FileNotFoundError, match="current source unavailable"):
            selected_source_sha256(tmp_path, relative)


@pytest.mark.parametrize("relative", [
    "docs/nsc-pg-tail-moment-audit.json",
    "docs/nsc-pg-retained-join-audit.json",
    "paper/local-gate-evidence/results/development/nsc-control.json",
    "paper/foundation-evidence/results/reference.json",
    "paper/local-gate-evidence/receipt.md",
    "paper/local-gate-evidence/src/owner.py",
    "lab/docs/nsc-pg-tail-moment-audit.json",
    "tests/fixtures/results/numerical.json",
    "paper/catalog.json",
])
def test_numerical_receipts_and_evidence_subtrees_never_select_git(tmp_path, relative):
    assert not is_implementation_source(relative)
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    raw = b'{"value": 7}\n'
    path.write_bytes(raw)
    expected = hashlib.sha256(raw).hexdigest()
    with frozen_source_context(COMMIT):
        assert selected_source_sha256(tmp_path, relative) == expected
        with pytest.raises(ValueError, match="implementation source"):
            authenticate_source(tmp_path, relative, expected, context=COMMIT)
        with pytest.raises(ValueError, match="implementation source"):
            source_sha256(tmp_path, relative, context=COMMIT)
        path.write_bytes(raw + b"changed\n")
        with pytest.raises(RuntimeError, match="current source drift"):
            authenticate_source(tmp_path, relative, expected)
        path.unlink()
        with pytest.raises(FileNotFoundError, match="current source unavailable"):
            selected_source_sha256(tmp_path, relative)


def test_public_source_predicate_is_strict_and_does_not_select_a_scope():
    assert is_implementation_source("scripts/original.py")
    assert is_implementation_source("docs/owner.md")
    assert is_implementation_source("pyproject.toml")
    assert not is_implementation_source("results/owner.py")
    with pytest.raises(ValueError, match="repository-relative"):
        is_implementation_source("../docs/owner.md")
    assert selected_source_context() is None
