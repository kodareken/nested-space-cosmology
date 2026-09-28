import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from recursive_horizons.provenance import ARCHIVE, resolve_source


class ProvenanceTests(unittest.TestCase):
    def test_archive_preserves_identity_and_active_drift_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            data = b"historical scientific source\n"
            digest = hashlib.sha256(data).hexdigest()
            stored = root / ARCHIVE / "scripts/original.py"
            stored.parent.mkdir(parents=True)
            stored.write_bytes(data)
            (root / ARCHIVE / "manifest.json").write_text(json.dumps({"entries": [{
                "original_path": "scripts/original.py", "sha256": digest,
                "storage_path": f"{ARCHIVE}/scripts/original.py"}]}))
            self.assertEqual(resolve_source(root, "scripts/original.py", digest), stored)
            active = root / "scripts/original.py"
            active.parent.mkdir()
            active.write_bytes(b"changed")
            with self.assertRaisesRegex(RuntimeError, "source drift"):
                resolve_source(root, "scripts/original.py", digest)
            with self.assertRaises(ValueError):
                resolve_source(root, "../original.py", digest)


if __name__ == "__main__":
    unittest.main()
