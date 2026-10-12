"""Test-only extraction of the immutable eight-leaf PREF27 GEN0 fixture.

The production calibration namespace grows after HLT16 starts.  Historical
GEN0 tests must therefore select the eight leaves already hash-bound by
PREF27, not copy the mutable successor tree or relax the production verifier.
"""
from __future__ import annotations

import json
from pathlib import Path
import shutil
from unittest.mock import patch

from recursive_horizons.fgc.evolution import proto18_pref27_binder as pref27
from recursive_horizons.fgc.evolution import proto19_progression_inputs as inputs


STORE_RELATIVE = Path("runs/fgc-2-sf1/proto17/calibration")


def sealed_gen0_leaves(root: Path) -> tuple[str, ...]:
    result = json.loads((root / "results/fgc-1-pro18-pref27.json").read_text("utf-8"))
    payload = result.get("artifact_payload")
    rows = payload.get("store_leaves") if isinstance(payload, dict) else None
    if (
        result.get("artifact_id") != "FGC-1-PRO18-PREF27"
        or not isinstance(rows, list)
        or len(rows) != 8
    ):
        raise AssertionError("sealed PREF27 fixture inventory differs")
    prefix = STORE_RELATIVE.as_posix() + "/"
    paths = tuple(sorted(str(row.get("path")) for row in rows if isinstance(row, dict)))
    if len(paths) != 8 or any(not path.startswith(prefix) for path in paths):
        raise AssertionError("sealed PREF27 fixture path differs")
    return paths


def copy_sealed_gen0_store(root: Path, destination: Path) -> None:
    prefix = STORE_RELATIVE.as_posix() + "/"
    destination.mkdir(parents=True, exist_ok=False)
    for path in sealed_gen0_leaves(root):
        relative = Path(path.removeprefix(prefix))
        source = root / path
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def reconstruct_sealed_gen0_members(root: Path):
    """Run the real GEN0 verifier over only PREF27's immutable leaf view."""
    expected_files = set(sealed_gen0_leaves(root))
    expected_directories = {
        f"{STORE_RELATIVE.as_posix()}/checkpoints",
        f"{STORE_RELATIVE.as_posix()}/receipts",
        f"{STORE_RELATIVE.as_posix()}/states",
    }
    inventory = pref27._nonfollowing_entries

    def sealed_inventory(store: Path, *, label: str):
        selected = []
        for path, info in inventory(store, label=label):
            relative = path.relative_to(root).as_posix()
            if relative in expected_files or relative in expected_directories:
                selected.append((path, info))
        return selected

    with patch.object(pref27, "_nonfollowing_entries", side_effect=sealed_inventory):
        evidence = pref27.verify_generation_zero_store(root)
    with patch.object(inputs, "verify_generation_zero_store", return_value=evidence):
        return inputs.reconstruct_gr0_members(root)


__all__ = ["copy_sealed_gen0_store", "reconstruct_sealed_gen0_members"]
