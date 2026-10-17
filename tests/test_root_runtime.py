"""The relocated runtime selects one root package and project-local dependencies."""
import os
from pathlib import Path
import sys

import pytest

import run_python


def test_runtime_discards_stale_source_paths(tmp_path, monkeypatch):
    validation = tmp_path / ".venv/validation"
    (validation / "Lib/site-packages").mkdir(parents=True)
    (validation / "pyvenv.cfg").write_text("executable = " + sys.executable + "\n")
    monkeypatch.setenv("PYTHONPATH", "old-worktree/lab/src")
    executable, packages, env = run_python.runtime(tmp_path)
    assert executable == Path(sys.executable)
    assert packages == validation / "Lib/site-packages"
    assert env["PYTHONPATH"].split(os.pathsep) == [str((tmp_path / "src").resolve()),
        str((tmp_path / "scripts").resolve()), str(packages.resolve())]
    assert env["PYTHONNOUSERSITE"] == "1"


def test_missing_validation_environment_is_explicit(tmp_path):
    with pytest.raises(RuntimeError, match="missing existing validation configuration"):
        run_python.runtime(tmp_path)


def test_complete_project_runtime_takes_precedence_and_isolates_sources(tmp_path, monkeypatch):
    prefix = tmp_path / ".venv/science-runtime"
    (prefix / "Lib/site-packages").mkdir(parents=True)
    (prefix / "conda-meta").mkdir()
    (prefix / "conda-meta/history").write_text("project-owned\n")
    (prefix / "python.exe").touch()
    monkeypatch.setenv("PYTHONHOME", "other-python")
    monkeypatch.setenv("PYTHONPATH", "other-project")
    monkeypatch.setenv("PATH", "existing-path")
    executable, packages, env = run_python.runtime(tmp_path)
    assert executable == prefix / "python.exe"
    assert packages == prefix / "Lib/site-packages"
    assert "PYTHONHOME" not in env
    assert env["PYTHONPATH"].split(os.pathsep) == [str((tmp_path / "src").resolve()),
        str((tmp_path / "scripts").resolve()), str(packages.resolve())]
    assert env["CONDA_PREFIX"] == str(prefix.resolve())
    assert env["PATH"].split(os.pathsep)[:3] == [str((prefix / p).resolve())
        for p in ("", "Library/bin", "Scripts")]
    assert env["PATH"].endswith("existing-path")
    assert os.environ["PYTHONHOME"] == "other-python"


def test_partial_project_runtime_is_not_silently_replaced(tmp_path):
    (tmp_path / ".venv/science-runtime").mkdir(parents=True)
    with pytest.raises(RuntimeError, match="incomplete project scientific environment"):
        run_python.runtime(tmp_path)
