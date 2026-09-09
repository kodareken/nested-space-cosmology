"""New source-fixture boundary and preservation of the published 91-record chain."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import tomllib
import unittest

from test_publication import ROOT,load_script

PIN='161028d52ef206a4bc99a25bf97ea467c30e7f16'
NAMES=('gauge-source','spherical-action','curvature-eft','spectral-endpoint','child-state',
       'massless-reference','angular-stress','unruh-state','state-regulator')


class V070ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer=load_script('reproduce_public_results.py')
        cls.importer=load_script('stage_release_import.py')
        cls.manifest=json.loads((ROOT/'results/manifest.json').read_text())
        cls.release=json.loads((ROOT/'results/release-spec.json').read_text())

    def test_current_version_and_historical_step_preservation(self):
        old=json.loads(subprocess.check_output(['git','show','v0.6.0:results/manifest.json'],cwd=ROOT))
        self.assertEqual(old['steps'],self.manifest['steps'][:91])
        self.assertEqual(100,len(self.manifest['steps']))
        self.assertEqual(58,self.manifest['historical_result_count'])
        self.assertEqual(42,len(self.release['scoped_follow_ups']))
        self.assertEqual(PIN,self.release['source_commit'])
        for version in (self.release['release_version'],self.manifest['release_version'],
                        json.loads((ROOT/'paper/metadata.json').read_text())['version'],
                        tomllib.loads((ROOT/'pyproject.toml').read_text())['project']['version']):
            self.assertEqual('0.7.0',version)
        self.assertIn('version: 0.7.0',(ROOT/'CITATION.cff').read_text())

    def test_new_policies_preserve_the_recorded_tolerances(self):
        self.assertEqual(['results/development/'+n+'.json' for n in NAMES],
                         [r['output'] for r in self.manifest['steps'][91:]])
        for step in self.manifest['steps'][91:]:
            record=json.loads((ROOT/step['output']).read_text())
            policy=step['comparison_policy'];original=record['comparison']
            self.assertEqual('all_fields',policy['kind'])
            self.assertEqual(original['float_atol'],policy['absolute_tolerance'])
            self.assertEqual(original['float_rtol'],policy['relative_tolerance'])
            self.assertEqual(record['schema'],step['identity_policy']['schema'])
            self.assertEqual(record['status'],step['identity_policy']['status'])
            self.assertEqual(PIN,step['follow_up_source_commit'])
            self.assertEqual(1,len(step['paper_claim_ids']))
            self.reproducer.validate_authenticated_inputs(ROOT,record,step,use_auxiliary=False)

    def test_fixture_preserves_the_original_configuration_and_native_source(self):
        step=self.manifest['steps'][97]  # angular-stress, order98
        self.assertEqual('NSC-28-ANGULAR-STRESS',step['artifact_id'])
        item=step['source_fixtures'][0]
        self.assertEqual('pyproject.toml',item['target'])
        self.assertEqual('tests/fixtures/source-checkpoints/'+PIN+'/pyproject.toml',item['path'])
        source=ROOT/item['path']
        record=json.loads((ROOT/step['output']).read_text())
        self.assertEqual(record['source_hashes']['pyproject.toml'],hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertNotEqual(source.read_bytes(),(ROOT/'pyproject.toml').read_bytes())
        native='src/recursive_horizons/_angular_transport.cpp'
        self.assertIn(native,step['source_dependencies'])
        self.assertEqual(record['source_hashes'][native],hashlib.sha256((ROOT/native).read_bytes()).hexdigest())

    def test_materialized_fixture_tampering_is_rejected(self):
        step=self.manifest['steps'][97];item=step['source_fixtures'][0]
        with tempfile.TemporaryDirectory() as temporary:
            work=Path(temporary);fixture=work/item['path'];fixture.parent.mkdir(parents=True)
            shutil.copyfile(ROOT/item['path'],fixture)
            shutil.copyfile(fixture,work/item['target'])
            self.reproducer.source_fixture_lookup(work,step)
            (work/item['target']).write_text('changed = true\n')
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.source_fixture_lookup(work,step)

    def test_fixture_paths_and_manifest_overrides_are_bounded(self):
        step=self.manifest['steps'][97]
        for target in ('../pyproject.toml','scripts/check_nsc_unruh_state.py','README.md'):
            changed=copy.deepcopy(step);changed['source_fixtures'][0]['target']=target
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.source_fixture_lookup(ROOT,changed)
        manifest=copy.deepcopy(self.manifest)
        manifest['steps'][97]['source_fixtures'][0]['sha256']='0'*64
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_checkout(manifest)
        for path in ('pyproject.toml','tests/private.toml','src/other.cpp'):
            with self.assertRaises(ValueError):self.importer.safe_relative(path)


if __name__=='__main__':unittest.main()
