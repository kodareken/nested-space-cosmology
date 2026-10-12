"""Isolated subprocess controls for the generic PRO19 execution closure."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/proto19_execution_closure.py"


def _git(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ("git", "-C", str(root), *arguments),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


ISOLATED_DRIVER = r"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

module_path = Path(sys.argv[1])
root = Path(sys.argv[2])
mode = sys.argv[3]
spec = importlib.util.spec_from_file_location("proto19_execution_closure_isolated", module_path)
if spec is None or spec.loader is None:
    raise RuntimeError("closure module cannot be loaded")
ec = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ec
spec.loader.exec_module(ec)

def git(*arguments):
    return subprocess.run(
        ("git", "-C", str(root), *arguments),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()

failure_modes = {
    "duplicate",
    "forbidden_module",
    "forbidden_script",
    "hardlink",
    "ignored_shadow",
    "later_source",
    "mutation",
    "outside_module",
    "path_escape",
    "postverify_mutation",
    "preloaded_outside",
    "symlink",
    "symlink_shadow",
    "untracked_shadow",
    "wrong_head",
}

try:
    implementation = git("rev-parse", "HEAD")
    record = ec.build_execution_closure_record(
        root,
        implementation_commit=implementation,
        files={
            "data/input.txt": "input",
            "src/fixture/helper.py": "python",
            "src/fixture/runtime.py": "python",
        },
        imports={
            "fixture.helper": "src/fixture/helper.py",
            "fixture.runtime": "src/fixture/runtime.py",
        },
        namespaces={"fixture": "src/fixture"},
        forbidden_module_prefixes=("fixture.candidate",),
        forbidden_script_prefixes=("scripts/run_candidate",),
        authority_delta_paths=("configs/authority.toml",),
    )
    if mode == "path_escape":
        mapping = record.to_mapping()
        mapping["files"][0]["path"] = "../escape.py"
        ec.parse_execution_closure_record(mapping)
        raise AssertionError("path escape was accepted")

    (root / "configs").mkdir(exist_ok=True)
    (root / "configs/authority.toml").write_text(
        "schema = 'fixture-authority'\n", encoding="utf-8"
    )
    if mode == "later_source":
        (root / "src/fixture/later.py").write_text("LATER = True\n", encoding="utf-8")
    git("add", ".")
    git("commit", "-qm", "authority")
    authority = git("rev-parse", "HEAD")

    if mode == "mutation":
        (root / "src/fixture/helper.py").write_text("VALUE = 99\n", encoding="utf-8")
    elif mode == "symlink":
        target = root / "src/fixture/helper.py"
        target.unlink()
        os.symlink(root / "data/input.txt", target)
    elif mode == "hardlink":
        os.link(root / "src/fixture/helper.py", root / "hardlink-owner.txt")
    elif mode == "symlink_shadow":
        os.symlink(root / "src/fixture", root / "fixture")
    elif mode == "untracked_shadow":
        (root / "src/fixture/rogue.py").write_text("ROGUE = True\n", encoding="utf-8")
    elif mode == "ignored_shadow":
        (root / "fixture").mkdir()
        (root / "fixture/runtime.py").write_text("RESULT = -1\n", encoding="utf-8")
    elif mode == "wrong_head":
        (root / "docs.txt").write_text("later authority drift\n", encoding="utf-8")
        git("add", "docs.txt")
        git("commit", "-qm", "wrong head")
    elif mode == "preloaded_outside":
        outside_path = root / "src/fixture/extra.py"
        outside_spec = importlib.util.spec_from_file_location("preloaded_rogue", outside_path)
        outside = importlib.util.module_from_spec(outside_spec)
        sys.modules[outside_spec.name] = outside
        outside_spec.loader.exec_module(outside)

    verified = ec.verify_execution_closure(
        root, record.to_mapping(), authority_commit=authority
    )
    if verified.record.implementation_commit != implementation:
        raise AssertionError("implementation commit was not preserved")
    if verified.authority_commit != authority or implementation == authority:
        raise AssertionError("two-commit authority topology was not preserved")
    if verified.record.interpreter.flags == ():
        raise AssertionError("isolated interpreter flags were not bound")
    if mode == "postverify_mutation":
        with ec.ImportOriginGuard(verified) as guard:
            (root / "src/fixture/helper.py").write_text("VALUE = 88\n", encoding="utf-8")
            guard.import_module("fixture.runtime")
    elif mode == "forbidden_script":
        with ec.ImportOriginGuard(verified) as guard:
            guard.assert_script_allowed(root / "scripts/run_candidate_resume.py")
    else:
        marker = root / "executed-marker"
        with ec.ImportOriginGuard(verified) as guard:
            if tuple(sys.path[:2]) != verified.import_path:
                raise AssertionError("src-before-root order differs")
            runtime = guard.import_module("fixture.runtime")
            if runtime.RESULT != 42:
                raise AssertionError("verified module result differs")
            guard.audit()
        if any(name == "fixture" or name.startswith("fixture.") for name in sys.modules):
            raise AssertionError("guarded repository modules leaked after cleanup")
        if tuple(sys.path) != record.interpreter.base_sys_path:
            raise AssertionError("isolated base sys.path was not restored")
        if marker.exists():
            raise AssertionError("a refused module executed before the guard stopped it")
    if mode in failure_modes:
        raise AssertionError(f"attack was accepted: {mode}")
    print(json.dumps({
        "state": "accepted",
        "implementation": implementation,
        "authority": authority,
        "record_sha256": record.canonical_sha256,
    }, sort_keys=True))
except ec.ExecutionClosureError as exc:
    if mode not in failure_modes:
        raise
    print("REFUSED:" + str(exc))
"""


class Proto19ExecutionClosureTests(unittest.TestCase):
    def make_repository(self, mode: str) -> tuple[TemporaryDirectory[str], Path, str]:
        directory = TemporaryDirectory()
        root = Path(directory.name)
        _git(root, "init", "-q")
        _git(root, "config", "user.email", "test@example.invalid")
        _git(root, "config", "user.name", "Test")
        _write(root / ".gitignore", "/fixture/\n")
        _write(root / "data/input.txt", "bound input\n")
        _write(root / "src/fixture/helper.py", "VALUE = 41\n")
        if mode == "forbidden_module":
            runtime = "from .candidate import payload\nRESULT = payload.VALUE\n"
        elif mode == "outside_module":
            runtime = "from . import extra\nRESULT = extra.VALUE\n"
        else:
            runtime = "from .helper import VALUE\nRESULT = VALUE + 1\n"
        _write(root / "src/fixture/runtime.py", runtime)
        _write(
            root / "src/fixture/candidate/payload.py",
            "from pathlib import Path\n(Path(__file__).parents[3] / 'executed-marker').write_text('bad')\nVALUE = -1\n",
        )
        _write(
            root / "src/fixture/extra.py",
            "from pathlib import Path\n(Path(__file__).parents[2] / 'executed-marker').write_text('bad')\nVALUE = -1\n",
        )
        _write(
            root / "scripts/run_candidate_resume.py",
            "raise RuntimeError('forbidden')\n",
        )
        if mode == "duplicate":
            _write(root / "fixture/runtime.py", "RESULT = -2\n")
            # The root fixture is normally ignored; force-add it to model a
            # tracked duplicate already present at implementation commit A.
            _git(root, "add", "-f", "fixture/runtime.py")
        _git(root, "add", ".")
        _git(root, "commit", "-qm", "implementation")
        return directory, root, _git(root, "rev-parse", "HEAD")

    def run_mode(self, mode: str) -> subprocess.CompletedProcess[str]:
        directory, root, _implementation = self.make_repository(mode)
        self.addCleanup(directory.cleanup)
        return subprocess.run(
            (
                sys.executable,
                "-I",
                "-B",
                "-c",
                ISOLATED_DRIVER,
                str(MODULE_PATH),
                str(root),
                mode,
            ),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_attached_descendant_authority_and_guarded_import_succeed(self) -> None:
        result = self.run_mode("success")
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["state"], "accepted")
        self.assertNotEqual(payload["implementation"], payload["authority"])
        self.assertEqual(len(payload["record_sha256"]), 64)

    def test_mutation_link_and_shadow_attacks_fail_closed(self) -> None:
        for mode in (
            "mutation",
            "symlink",
            "symlink_shadow",
            "hardlink",
            "untracked_shadow",
            "ignored_shadow",
            "duplicate",
        ):
            with self.subTest(mode=mode):
                result = self.run_mode(mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(result.stdout.startswith("REFUSED:"), result.stdout)

    def test_commit_path_and_preimport_attacks_fail_closed(self) -> None:
        for mode in ("later_source", "wrong_head", "path_escape", "preloaded_outside"):
            with self.subTest(mode=mode):
                result = self.run_mode(mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(result.stdout.startswith("REFUSED:"), result.stdout)

    def test_forbidden_and_outside_closure_imports_never_execute(self) -> None:
        for mode in ("forbidden_module", "outside_module", "forbidden_script"):
            with self.subTest(mode=mode):
                result = self.run_mode(mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(result.stdout.startswith("REFUSED:"), result.stdout)

    def test_postverification_mutation_is_rechecked_before_import(self) -> None:
        result = self.run_mode("postverify_mutation")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.stdout.startswith("REFUSED:"), result.stdout)

    def test_module_loads_without_numpy_or_project_package_imports(self) -> None:
        code = r"""
import importlib.util
import json
import sys
before = set(sys.modules)
spec = importlib.util.spec_from_file_location("closure_import_purity", sys.argv[1])
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
added = set(sys.modules) - before
print(json.dumps({
    "numpy": any(name == "numpy" or name.startswith("numpy.") for name in added),
    "project": any(name == "recursive_horizons" or name.startswith("recursive_horizons.") for name in added),
    "isolated": sys.flags.isolated,
    "dont_write_bytecode": sys.flags.dont_write_bytecode,
}, sort_keys=True))
"""
        result = subprocess.run(
            (sys.executable, "-I", "-B", "-c", code, str(MODULE_PATH)),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        payload = json.loads(result.stdout)
        self.assertEqual(
            payload,
            {
                "dont_write_bytecode": 1,
                "isolated": 1,
                "numpy": False,
                "project": False,
            },
        )

    def test_nonisolated_interpreter_is_refused(self) -> None:
        code = r"""
import importlib.util
import sys
spec = importlib.util.spec_from_file_location("closure_nonisolated", sys.argv[1])
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
try:
    module.current_interpreter_identity()
except module.ExecutionClosureError as exc:
    print("REFUSED:" + str(exc))
else:
    raise AssertionError("nonisolated interpreter was accepted")
"""
        result = subprocess.run(
            (sys.executable, "-B", "-c", code, str(MODULE_PATH)),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.assertTrue(
            result.stdout.startswith("REFUSED:bootstrap requires python -I -B")
        )


if __name__ == "__main__":
    unittest.main()
