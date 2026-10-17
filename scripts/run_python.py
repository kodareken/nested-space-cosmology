#!/usr/bin/env python3
"""Run root commands in the project-owned Windows scientific environment.

Prefer the complete local Conda prefix. The older validation-package fallback
remains available for checkouts without that prefix. Neither route installs
packages or changes the user's other Python installations.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def runtime(root: Path = ROOT) -> tuple[Path, Path, dict[str, str]]:
    isolated = root / ".venv/science-runtime"
    if isolated.is_dir():
        executable = isolated / "python.exe"
        site_packages = isolated / "Lib/site-packages"
        if not (executable.is_file() and site_packages.is_dir()
                and (isolated / "conda-meta/history").is_file()):
            raise RuntimeError("incomplete project scientific environment: " + str(isolated))
        env = source_environment(root, site_packages)
        env.pop("PYTHONHOME", None)
        env["CONDA_PREFIX"] = str(isolated.resolve())
        env["PATH"] = os.pathsep.join(str((isolated / part).resolve())
                                    for part in ("", "Library/bin", "Scripts")) + os.pathsep + env.get("PATH", "")
        return executable, site_packages, env
    validation = root / ".venv/validation"
    cfg = validation / "pyvenv.cfg"
    if not cfg.is_file():
        raise RuntimeError("missing existing validation configuration: " + str(cfg))
    values = dict(line.split("=", 1) for line in cfg.read_text(encoding="utf-8").splitlines() if "=" in line)
    values = {key.strip(): value.strip() for key, value in values.items()}
    executable = Path(values["executable"])
    site_packages = validation / "Lib/site-packages"
    if not executable.is_file() or not site_packages.is_dir():
        raise RuntimeError("the configured base interpreter or validation packages are unavailable")
    return executable, site_packages, source_environment(root, site_packages)


def source_environment(root: Path, site_packages: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(str(path.resolve()) for path in
                                       (root / "src", root / "scripts", site_packages))
    env["PYTHONNOUSERSITE"] = "1"
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(root / "archive/source-history/objects")
    return env


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("Usage: python scripts/run_python.py <script.py | -m module | -c code>", file=sys.stderr)
        return 2
    executable, site_packages, env = runtime()
    if args == ["--runtime-info"]:
        print(json.dumps({"executable": str(executable), "site_packages": str(site_packages),
                          "source": str(ROOT / "src"), "user_site_enabled": False}, indent=2))
        return 0
    return subprocess.call([str(executable), *args], cwd=ROOT, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
