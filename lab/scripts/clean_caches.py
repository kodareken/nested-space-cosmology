#!/usr/bin/env python3
"""Move confirmed Python/Ruff/macOS cache files to recoverable Trash."""
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    command = shutil.which("trashit")
    if command is None:
        raise SystemExit("No recoverable trashit command; no files removed.")
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0"))
    targets = []
    for base in (ROOT, ROOT / "src", ROOT / "scripts", ROOT / "tests"):
        paths = list(base.iterdir()) if base == ROOT else list(base.rglob("*"))
        for path in paths:
            if path.is_symlink() or str(path.relative_to(ROOT)) in tracked:
                continue
            if path.name in {"__pycache__", ".ruff_cache", ".pytest_cache"} or path.name == ".DS_Store":
                if any(p.startswith(str(path.relative_to(ROOT)) + "/") for p in tracked):
                    continue
                targets.append(path)
    unique = [p for p in sorted(set(targets)) if not any(q in p.parents for q in targets)]
    if unique:
        subprocess.run([command, *map(str, unique)], check=True)
    print(f"Moved {len(unique)} disposable cache paths to Trash; scientific files preserved.")


if __name__ == "__main__":
    main()
