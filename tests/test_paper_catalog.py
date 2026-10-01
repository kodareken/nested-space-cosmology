"""Prevent either manuscript disappearing or being silently relabelled."""
import importlib.util
import copy
import hashlib
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
            if 'source_archive' in row:
                files.add(row['source_archive'])
            for prior in row.get('historical_editions', []):
                files.update(prior[k] for k in ['source', 'pdf', 'manifest', 'source_archive'])
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

    def measured_catalog(self):
        """A fixture for an edition transition, retaining the historical bound artifact."""
        active = self.catalog['documents'][1]
        prior = copy.deepcopy(active.get('historical_editions', [active])[0])
        archive = prior.get('source_archive', 'paper/local-incoming-gate-draft-source.tar.gz')
        shutil.copyfile(ROOT / archive, self.root / archive)
        source, pdf, manifest = ('paper/new-edition/main.tex', 'paper/new-edition.pdf',
                                 'paper/new-edition-manifest.json')
        (self.root / source).parent.mkdir(parents=True)
        (self.root / source).write_bytes(b'new measured edition source')
        (self.root / pdf).write_bytes(b'new measured edition PDF')
        current_archive = 'paper/new-edition-source.tar.gz'
        (self.root / current_archive).write_bytes(b'new measured edition source archive')
        digest = lambda path: hashlib.sha256((self.root / path).read_bytes()).hexdigest()
        (self.root / manifest).write_text(json.dumps({
            'schema': 'NSC-FINITE-REGENERATION-MANIFEST-v1',
            'status': 'MEASURED_FINITE_REALIZATION', 'submission_ready': False,
            'author_review': 'pending', 'science_commit': '1' * 40,
            'pdf_sha256': digest(pdf), 'source_archive_sha256': digest(current_archive),
            'inputs': {source: digest(source)},
        }))
        prior['source_archive'] = archive
        active.update(role='focused_companion', status='MEASURED_FINITE_REALIZATION',
                      author_review='pending', source=source, pdf=pdf, manifest=manifest,
                      source_archive=current_archive, historical_editions=[prior])
        self.catalog.update(schema='NSC-PAPER-CATALOG-v2',
                            current_status='MEASURED_FINITE_REALIZATION',
                            finite_result_status='MEASURED_FINITE_REALIZATION',
                            physical_local_gate_status='OPEN')
        self.save_catalog()

    def test_measured_edition_keeps_two_roles_and_separate_open_application(self):
        self.measured_catalog()
        self.assertEqual(self.errors(), [])
        self.catalog['physical_local_gate_status'] = 'CLOSED'
        self.save_catalog()
        self.assertTrue(self.errors())

    def test_measured_edition_cannot_drop_or_overwrite_its_history(self):
        self.measured_catalog()
        active = self.catalog['documents'][1]
        history = active.pop('historical_editions')
        self.save_catalog()
        self.assertTrue(self.errors())
        active['historical_editions'] = history
        self.save_catalog()
        (self.root / history[0]['pdf']).write_bytes(b'overwritten historical edition')
        self.assertTrue(self.errors())

    def test_measured_label_requires_its_own_built_source_and_review_domain(self):
        self.measured_catalog()
        self.catalog['documents'][1]['author_review'] = 'approved'
        self.save_catalog()
        self.assertTrue(self.errors())


if __name__=='__main__':unittest.main()
