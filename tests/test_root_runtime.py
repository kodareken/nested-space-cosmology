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
