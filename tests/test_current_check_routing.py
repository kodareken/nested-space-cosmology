"""The legacy command must use consolidated owners without old Git objects."""
import importlib.util
from pathlib import Path


def test_current_checks_use_existing_owners_and_propagate_failure(monkeypatch):
    path = Path(__file__).resolve().parents[1] / "lab/scripts/check_current.py"
    spec = importlib.util.spec_from_file_location("current_check_routing", path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    calls = []

    def invoke(command, **kwargs):
        calls.append(Path(command[-1]).name)
        return 0

    monkeypatch.setattr(checker.subprocess, "call", invoke)
    assert checker.main([]) == 0
    assert calls == ["check_lab_snapshot.py", "check_publication.py"]
    calls.clear()

    def fail(command, **kwargs):
        calls.append(Path(command[-1]).name)
        return 3

    monkeypatch.setattr(checker.subprocess, "call", fail)
    assert checker.main([]) == 3
    assert calls == ["check_lab_snapshot.py"]


def test_public_envelope_ignores_local_agent_routing_but_keeps_research(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    (tmp_path / ".pi/npm").mkdir(parents=True)
    (tmp_path / ".pi/MATH_SYSTEM.md").write_text("local routing")
    (tmp_path / ".pi/npm/.gitignore").write_text("*")
    (tmp_path / "THEORY.md").write_text("scientific narrative")
    for name, method in (("build_publication_provenance.py", "included_files"),
                         ("check_publication.py", "public_files")):
        spec = importlib.util.spec_from_file_location(name[:-3], root / "scripts" / name)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        monkeypatch.setattr(module, "ROOT", tmp_path)
        names = {p.relative_to(tmp_path).as_posix() for p in getattr(module, method)()}
        assert names == {"THEORY.md"}
