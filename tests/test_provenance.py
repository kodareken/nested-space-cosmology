import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from recursive_horizons.provenance import resolve_pinned_source_bytes, resolve_source


class ProvenanceTests(unittest.TestCase):
    def test_explicit_review_replay_authenticates_exact_commit(self):
        root = Path(__file__).resolve().parents[1]
        path = "lab/tests/test_nsc_spherical_null_expansion.py"
        expected = "a7588bfac6bbe7049c2b64d0df7539b497e2cb8f9273dc2606797eb30b8d18aa"
        commit = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"
        raw = resolve_pinned_source_bytes(root, path, expected, commit=commit)
        self.assertEqual(hashlib.sha256(raw).hexdigest(), expected)
        with self.assertRaisesRegex(RuntimeError, "current source unavailable"):
            resolve_source(root, path, expected)
        with self.assertRaisesRegex(RuntimeError, "pinned source drift"):
            resolve_pinned_source_bytes(root, path, "0" * 64, commit=commit)
        with self.assertRaisesRegex(RuntimeError, "unavailable"):
            resolve_pinned_source_bytes(root, path, expected, commit="0" * 40)
        with self.assertRaisesRegex(ValueError, "implementation source"):
            resolve_pinned_source_bytes(root, "lab/results/review.json", expected, commit=commit)

    def test_explicit_cache_pin_replays_without_rebinding_current_source(self):
        root = Path(__file__).resolve().parents[1]
        manifest = json.loads((root / "archive/source-history/manifest.json").read_text())
        pin = manifest["pins"][-1]
        raw = resolve_pinned_source_bytes(
            root, pin["path"], pin["sha256"], commit=pin["commit"])
        self.assertEqual(hashlib.sha256(raw).hexdigest(), pin["sha256"])

    def test_retired_receipts_never_substitute_current_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            data = b"historical scientific source\n"
            digest = hashlib.sha256(data).hexdigest()
            retired = "archive/pre-closure-2026-09-08"
            stored = root / retired / "scripts/original.py"
            stored.parent.mkdir(parents=True)
            stored.write_bytes(data)
            (root / retired / "manifest.json").write_text(json.dumps({"entries": [{
                "original_path": "scripts/original.py", "sha256": digest,
                "storage_path": f"{retired}/scripts/original.py"}]}))
            with self.assertRaisesRegex(RuntimeError, "current source unavailable"):
                resolve_source(root, "scripts/original.py", digest)
            migration = root / "results/reassembly/manifest.json"
            migration.parent.mkdir(parents=True)
            migration.write_text(json.dumps({"rows": [{
                "origin": "lab", "original_path": "scripts/original.py",
                "destination": f"{retired}/scripts/original.py"}]}))
            with self.assertRaisesRegex(RuntimeError, "current source unavailable"):
                resolve_source(root, "scripts/original.py", digest)
            active = root / "scripts/original.py"
            active.parent.mkdir()
            active.write_bytes(data)
            self.assertEqual(resolve_source(root, "scripts/original.py", digest), active)
            active.write_bytes(b"changed")
            with self.assertRaisesRegex(RuntimeError, "source drift"):
                resolve_source(root, "scripts/original.py", digest)
            with self.assertRaises(ValueError):
                resolve_source(root, "../original.py", digest)


if __name__ == "__main__":
    unittest.main()
