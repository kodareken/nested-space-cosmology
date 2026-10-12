#!/usr/bin/env python3
"""Run root commands with the existing validation packages and base interpreter.

The moved Windows virtual-environment redirector cannot locate pyvenv.cfg on
this host. This launcher reads that configuration directly, uses its declared
base executable, and supplies only canonical root sources and local packages.
It does not install packages, alter permissions, or use a competing checkout.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def runtime(root: Path = ROOT) -> tuple[Path, Path, dict[str, str]]:
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
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(str(path.resolve()) for path in
                                       (root / "src", root / "scripts", site_packages))
    env["PYTHONNOUSERSITE"] = "1"
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(root / "archive/source-history/objects")
    return executable, site_packages, env


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
