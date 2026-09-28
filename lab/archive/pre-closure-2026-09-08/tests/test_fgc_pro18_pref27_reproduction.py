from __future__ import annotations

from hashlib import sha256
import importlib.util
import subprocess
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_STORE = ROOT / "runs/fgc-2-sf1/proto17/calibration"


def reproducer_module():
    path = ROOT / "scripts/reproduce_fgc_pro18_pref27.py"
    spec = importlib.util.spec_from_file_location("pref27_reproducer_test", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tree_digest(path: Path) -> tuple[tuple[str, str], ...]:
    return tuple(
        (entry.relative_to(path).as_posix(), sha256(entry.read_bytes()).hexdigest())
        for entry in sorted(path.rglob("*")) if entry.is_file()
    )


class Pref27ReproductionTests(unittest.TestCase):
    def test_duplicate_and_noncanonical_stored_results_are_rejected(self) -> None:
        module = reproducer_module()
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            path.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_canonical(path)
            path.write_text('{"a":1}\n', encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_canonical(path)

    def test_fixed_pref27_result_reproduces_without_advancing_canonical_store(self) -> None:
        if not CANONICAL_STORE.is_dir():  # pragma: no cover - concrete production prerequisite
            self.skipTest("HLT15 canonical generation-zero store is absent")
        before = tree_digest(CANONICAL_STORE)
        completed = subprocess.run(
            [sys.executable, "scripts/reproduce_fgc_pro18_pref27.py", "--verify"],
            cwd=ROOT, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn('"artifact_id": "FGC-1-PRO18-PREF27"', completed.stdout)
        self.assertTrue((ROOT / "results/fgc-1-pro18-pref27.json").is_file())
        self.assertEqual(tree_digest(CANONICAL_STORE), before)


if __name__ == "__main__":
    unittest.main()
