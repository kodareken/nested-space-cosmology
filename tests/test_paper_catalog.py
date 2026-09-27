"""Prevent either manuscript disappearing or being silently relabelled."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('publication_catalog_check',ROOT/'scripts/check_publication.py')
checker=importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class PaperCatalogTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.catalog=json.loads((ROOT/'paper/catalog.json').read_text())
        files={'paper/catalog.json'}
        for row in self.catalog['documents']:
            files.update(row[k] for k in ['source','pdf','manifest'])
        for name in files:
            p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,p)
        self.patch=patch.object(checker,'ROOT',self.root)
        self.patch.start();self.addCleanup(self.patch.stop)

    def errors(self):
        errors=[];checker.check_paper_catalog(errors);return errors

    def save_catalog(self):
        (self.root/'paper/catalog.json').write_text(json.dumps(self.catalog))

    def test_both_papers_and_builds_are_bound(self):
        self.assertEqual(self.errors(),[])

    def test_missing_companion_is_detected(self):
        (self.root/self.catalog['documents'][1]['pdf']).unlink()
        self.assertTrue(self.errors())

    def test_replaced_foundation_is_detected(self):
        (self.root/self.catalog['documents'][0]['pdf']).write_bytes(b'new unrelated PDF')
        self.assertTrue(self.errors())

    def test_swapped_or_duplicate_roles_are_detected(self):
        self.catalog['documents'][1]['role']='foundational_manuscript';self.save_catalog()
        self.assertTrue(self.errors())
        self.catalog['documents'][1]=self.catalog['documents'][0];self.save_catalog()
        self.assertTrue(self.errors())

    def test_false_closed_companion_is_detected(self):
        self.catalog['documents'][1]['status']='EXISTENCE';self.save_catalog()
        self.assertTrue(self.errors())

    def test_changed_source_requires_a_rebuilt_manifest(self):
        p=self.root/self.catalog['documents'][1]['source'];p.write_bytes(p.read_bytes()+b'\nchanged')
        self.assertTrue(self.errors())


if __name__=='__main__':unittest.main()
