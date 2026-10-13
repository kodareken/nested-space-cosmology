"""Current Atlas owner bindings, separate from frozen scientific receipts.

Only declared source paths are normalized. Numerical inputs and generated
outputs use their original local paths and are never routed through this map.
"""
from __future__ import annotations

import hashlib
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping


# Ultimate owners from FINAL-MERGE.tsv, with the current Atlas index owning
# archived edition explanations and the retained v7 test owning earlier
# producer coverage. The living paper consumes this evidence.
_CURRENT_OWNERS = {
    "docs/child-neck/balanced-layers-note.md": "docs/child-neck/child-and-neck-status-note.md",
    "docs/child-neck/israel-matching-note.md": "docs/child-neck/child-and-neck-status-note.md",
    "docs/nsc-atlas-connection-observer.md": "docs/nsc-atlas-causal-observers.md",
    "docs/nsc-atlas-junction-sources.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-atlas-parent-gluing.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-atlas-parent-initial-data.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-clock-horizon.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v2.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v3.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v4.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v5.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v6.md": "docs/nsc-cosmological-atlas.md",
    "docs/nsc-cosmological-atlas-v7.md": "docs/nsc-cosmological-atlas-v8.md",
    "docs/nsc-vacuum-charge-matching.md": "docs/prior-art-and-open-claim.md",
    "paper/atlas-working-draft.md": "docs/nsc-cosmological-atlas.md",
    "paper/atlas-working-draft-v3.md": "docs/nsc-cosmological-atlas.md",
    "paper/atlas-working-draft-v4.md": "docs/nsc-cosmological-atlas.md",
    "paper/atlas-working-draft-v5.md": "docs/nsc-cosmological-atlas.md",
    "paper/atlas-working-draft-v6.md": "docs/nsc-cosmological-atlas.md",
    "tests/test_nsc_cosmological_atlas_v2.py": "tests/test_nsc_cosmological_atlas_v7.py",
    "tests/test_nsc_cosmological_atlas_v3.py": "tests/test_nsc_cosmological_atlas_v7.py",
    "tests/test_nsc_cosmological_atlas_v4.py": "tests/test_nsc_cosmological_atlas_v7.py",
    "tests/test_nsc_cosmological_atlas_v5.py": "tests/test_nsc_cosmological_atlas_v7.py",
    "tests/test_nsc_cosmological_atlas_v6.py": "tests/test_nsc_cosmological_atlas_v7.py",
}
_LIVING_PAPER = {"paper/nested-space-cosmology.md", "paper/nested-space-cosmology.pdf"}


def _relative_source(root: Path, source: str | Path) -> str:
    path = Path(source)
    if path.is_absolute():
        path = path.relative_to(root)
    name = path.as_posix()
    parts = PurePosixPath(name).parts
    if not parts or name in {"", "."} or any(part in {".", ".."} for part in parts):
        raise ValueError("Atlas source binding must use a repository-relative file")
    if name in _LIVING_PAPER:
        raise ValueError("the living paper is an Atlas consumer, not a source binding")
    return _CURRENT_OWNERS.get(name, name)


def current_source_bindings(root: Path, sources: Iterable[str | Path], *,
                            inherited: Mapping[str, str] | None = None) -> dict[str, str]:
    """Hash current retained owners and name the exact files whose bytes were read.

    Inherited keys are gathered as paths and rehashed, rather than carrying an
    old digest forward. Missing files and changed numerical inputs remain hard
    failures. Frozen producer hashes belong to explicit historical replay.
    """
    root = Path(root).resolve()
    names = {_relative_source(root, source) for source in sources}
    names.update(_relative_source(root, source) for source in (inherited or {}))
    names.add(_relative_source(root, Path(__file__).resolve()))
    bindings = {}
    for name in sorted(names):
        path = root / name
        path.resolve(strict=True).relative_to(root)
        bindings[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return bindings
