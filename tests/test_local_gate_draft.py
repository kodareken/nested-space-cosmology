"""Adversarial integrity checks for the public OPEN draft (standard library)."""
import io
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import local_gate_draft as draft


class LocalGateDraftTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        shutil.copytree(ROOT / draft.EVIDENCE, self.root / draft.EVIDENCE)
        shutil.copytree(ROOT / draft.SOURCE, self.root / draft.SOURCE)

    def test_original_snapshot_and_generated_values(self):
        self.assertEqual(draft.authenticate_snapshot(self.root)['status'], 'OPEN')
        self.assertIn('main.tex', draft.source_files(self.root))

    def test_corrupted_original_record_rejected(self):
        path = self.root / draft.EVIDENCE / 'results/development/nsc-ks-gate-budget-v4.json'
        path.write_bytes(path.read_bytes() + b' ')
        with self.assertRaisesRegex(ValueError, 'evidence mismatch'):
            draft.authenticate_snapshot(self.root)

    def test_closed_claim_rejected(self):
        path = self.root / draft.SOURCE / 'result_macros.tex'
        path.write_text(path.read_text().replace('"verdict":"OPEN"', '"verdict":"EXISTENCE"'))
        with self.assertRaisesRegex(ValueError, 'closed certificate'):
            draft.source_files(self.root)

    def test_stale_values_rejected(self):
        path = self.root / draft.SOURCE / 'evidence_values.tex'
        path.write_text(path.read_text().replace('9.64093', '0.00000'))
        with self.assertRaisesRegex(ValueError, 'evidence values'):
            draft.source_files(self.root)

    def test_unsafe_paths_rejected(self):
        for name in ('../secret', '/etc/passwd', 'C:/secret', 'x\\y', '.hidden', 'a/../b'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                draft.safe_path(self.root, name)

    def test_forged_git_blob_rejected(self):
        path = self.root / draft.EVIDENCE / 'snapshot.json'
        record = json.loads(path.read_text())
        record['files'][0]['git_blob'] = 'a' * 40
        path.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, 'Git blob'):
            draft.authenticate_snapshot(self.root)

    def test_source_archive_is_deterministic_and_plain_files_only(self):
        files = draft.source_files(self.root)
        first = draft.packed_source(files)
        self.assertEqual(first, draft.packed_source(dict(reversed(list(files.items())))))
        with tarfile.open(fileobj=io.BytesIO(first), mode='r:gz') as tar:
            self.assertEqual(set(tar.getnames()), set(files))
            self.assertTrue(all(m.isfile() and m.mtime == 0 for m in tar.getmembers()))
        with self.assertRaises(ValueError):
            draft.packed_source({'../main.tex': b'bad'})


if __name__ == '__main__':
    unittest.main()
