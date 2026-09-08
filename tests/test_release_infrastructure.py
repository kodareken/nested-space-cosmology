"""Negative controls for the public importer and authenticated reproduction DAG."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import unittest

from test_publication import ROOT, load_script


class ReleaseInfrastructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer = load_script('reproduce_public_results.py')
        cls.importer = load_script('stage_release_import.py')
        cls.release = json.loads((ROOT / 'results/release-spec.json').read_text())
        cls.manifest = cls.reproducer.load_manifest()
        cls.steps = {row['artifact_id']: row for row in cls.manifest['steps']}

    def test_all_64_existing_scientific_files_and_pinned_imports_preserved(self):
        preserved = [
            entry
            for key, value in self.release.items()
            if key.startswith('preserved_') and key.endswith('_scientific_files')
            for entry in value
        ]
        for entry in (preserved + self.release['import_files']
                      + self.release['retained_byte_identical_files']):
            self.assertEqual(entry['sha256'], hashlib.sha256((ROOT / entry['path']).read_bytes()).hexdigest())
        result_paths = {
            entry['path']
            for entry in preserved + self.release['import_files']
            if entry['path'].startswith('results/') and entry['path'].endswith('.json')
        }
        self.assertEqual({row['output'] for row in self.manifest['steps']}, result_paths)

    def test_explicit_smooth_module_generator_and_full_source_closure(self):
        smooth = self.steps['NSC-4-SMOOTH-GEOMETRY']
        self.assertEqual('scripts/check_nsc_smooth_geometry.py', smooth['generator'])
        self.assertEqual({'results/nsc-3-geometric-chain.json',
                          'results/nsc-3-regulated-recursion.json',
                          'results/nsc-3-radial-spectrum.json'}, set(smooth['dependencies']))
        for row in self.manifest['steps'][58:]:
            self.assertEqual(set(row['source_dependencies']), set(row['source_dependency_hashes']))
        self.assertIn('tests/test_nsc_smooth_geometry.py', smooth['source_dependencies'])
        self.assertIn('docs/nsc-smooth-geometry.md', smooth['source_dependencies'])
        self.assertIn('src/recursive_horizons/nsc_geometric_chain.py', smooth['source_dependencies'])

    def test_current_physical_targets_are_explicit_and_open(self):
        targets = self.release['current_targets']
        self.assertEqual(6, len(targets))
        for target in targets:
            self.assertEqual('open_target', target['status'])
            self.assertTrue(target['required_physical_result'])
            self.assertLessEqual(set(target['supporting_artifact_ids']), set(self.steps))

    def test_plateau_schema_scope_and_false_prediction_are_required(self):
        step = self.steps['NSC-5-PLATEAU-CONDITIONS']
        value = json.loads((ROOT / step['output']).read_text())
        self.reproducer.validate_identity(value, step)
        self.assertNotIn('terminal', value)
        self.assertNotIn('classification', value)
        self.assertNotIn('nonclaims', value)
        for field in ('schema', 'scope', 'artifact_id'):
            wrong = copy.deepcopy(value)
            wrong[field] += ' altered'
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_identity(wrong, step)
        wrong = copy.deepcopy(value)
        wrong['observational_audit']['prediction_claim'] = True
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_identity(wrong, step)

    def test_all_hash_schemas_reject_tampering_before_normalization(self):
        for artifact in ('NSC-4-SHAPE-RESPONSE', 'NSC-4-LORENTZIAN-TRANSPORT',
                         'NSC-6-ENERGY-TRANSFER', 'NSC-4-COVARIANT-MEASURE',
                         'NSC-8-CHIRAL-BOUNDARY', 'NSC-8-FINITE-TERMS',
                         'NSC-9-COVARIANT-SOURCE', 'NSC-10-MEASURE-NORMALIZATION',
                         'NSC-10-INFLUENCE'):
            step = self.steps[artifact]
            original = json.loads((ROOT / step['output']).read_text())
            locations = self.reproducer.validate_authenticated_inputs(ROOT, original, step)
            normalized = self.reproducer.normalize_dynamic_hashes(original, locations)
            self.assertEqual(original['source_hashes'], normalized['source_hashes'])
            if 'comparison_function_sha256' in original:
                self.assertEqual(original['comparison_function_sha256'], normalized['comparison_function_sha256'])
            for location, _, _, _ in self.reproducer.hash_bindings(original):
                wrong = copy.deepcopy(original)
                parent = wrong
                for part in location[:-1]:
                    parent = parent[part]
                parent[location[-1]] = '0' * 64
                with self.assertRaises(self.reproducer.ReproductionError):
                    self.reproducer.validate_authenticated_inputs(ROOT, wrong, step)

    def test_comparator_hash_is_never_ignored(self):
        step = self.steps['NSC-4-COVARIANT-MEASURE']
        value = json.loads((ROOT / step['output']).read_text())
        value['comparison_function_sha256'] = '0' * 64
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_authenticated_inputs(ROOT, value, step)

    def test_normalization_requires_declared_regenerated_results(self):
        step = self.steps['NSC-4-SHAPE-RESPONSE']
        value = json.loads((ROOT / step['output']).read_text())
        locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
        self.assertEqual(2, len(locations))
        unchanged = self.reproducer.normalize_dynamic_hashes(value, set())
        self.assertEqual(value, unchanged)
        source = next(iter(value['source_hashes']))
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.normalize_dynamic_hashes(value, {('source_hashes', source)})
        auxiliary_step = self.steps['NSC-2-ZETA1-UNIT-CLOSURE-CHECK']
        auxiliary = json.loads((ROOT / auxiliary_step['output']).read_text())
        locations = self.reproducer.validate_authenticated_inputs(ROOT, auxiliary, auxiliary_step, use_auxiliary=False)
        self.assertEqual(set(), locations)

    def test_missing_result_dependency_and_changed_source_closure_fail(self):
        step = copy.deepcopy(self.steps['NSC-4-SHAPE-RESPONSE'])
        value = json.loads((ROOT / step['output']).read_text())
        step['dependencies'] = []
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_authenticated_inputs(ROOT, value, step)
        manifest = copy.deepcopy(self.manifest)
        row = next(item for item in manifest['steps'] if item['artifact_id'] == 'NSC-4-SHAPE-RESPONSE')
        relative = next(iter(row['source_dependency_hashes']))
        row['source_dependency_hashes'][relative] = '0' * 64
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_checkout(manifest)

    def test_malformed_hash_mappings_fail(self):
        for value in ({'authenticated_input_hashes': []},
                      {'input_hashes': {'results/example.json': {'sha256': '0' * 64}}},
                      {'source_hashes': {'../outside.py': '0' * 64}},
                      {'authenticated_inputs': [{'path': 'results/example.json'}]}):
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_authenticated_inputs(ROOT, value)

    def test_declared_numeric_tolerances_have_negative_controls(self):
        specifications = {
            'NSC-4-SHAPE-RESPONSE': (1e-6, 2e-7),
            'NSC-5-PLATEAU-CONDITIONS': (2e-10, 3e-12),
            'NSC-7-OBSERVABLE-BRIDGE': (2e-10, 3e-12),
            'NSC-6-VACUUM-WORK': (1e-7, 5e-10),
            'NSC-8-FINITE-TERMS': (2e-8, 2e-6),
            'NSC-8-CHIRAL-BOUNDARY': (1e-8, 1e-8),
            'NSC-9-COVARIANT-SOURCE': (2e-8, 2e-7),
            'NSC-10-MEASURE-NORMALIZATION': (2e-8, 2e-7),
            'NSC-10-INFLUENCE': (2e-8, 2e-7),
            'NSC-11-RESPONSE-MATCHING': (2e-8, 2e-8),
        }
        for artifact, (relative, absolute) in specifications.items():
            policy = self.steps[artifact]['comparison_policy']
            self.assertEqual((relative, absolute), (policy['relative_tolerance'], policy['absolute_tolerance']))
            for expected, increment in ((1000.0, relative * 1000.0), (0.0, absolute)):
                self.reproducer.compare_portable(expected, expected + increment / 2,
                    path='/control', relative_tolerance=relative,
                    absolute_tolerance=absolute, compare_numbers=True)
                with self.assertRaises(self.reproducer.ReproductionError):
                    self.reproducer.compare_portable(expected, expected + increment * 4,
                        path='/control', relative_tolerance=relative,
                        absolute_tolerance=absolute, compare_numbers=True)

    def test_chiral_identity_requires_its_schema_classification_and_exact_gate(self):
        step = self.steps['NSC-8-CHIRAL-BOUNDARY']
        original = json.loads((ROOT / step['output']).read_text())
        self.assertNotIn('artifact_id', original)
        self.assertNotIn('terminal', original)
        self.reproducer.validate_identity(original, step)
        for field in ('schema', 'classification'):
            wrong = copy.deepcopy(original)
            wrong[field] += ' changed'
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_identity(wrong, step)
        for replacement in (False, 1, None):
            wrong = copy.deepcopy(original)
            wrong['gate'][next(iter(wrong['gate']))] = replacement
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_identity(wrong, step)
        wrong = copy.deepcopy(original)
        del wrong['gate'][next(iter(wrong['gate']))]
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.validate_identity(wrong, step)

    def test_finite_terms_document_hashes_stay_exact_and_only_three_results_normalize(self):
        step = self.steps['NSC-8-FINITE-TERMS']
        value = json.loads((ROOT / step['output']).read_text())
        locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
        self.assertEqual(3, len(locations))
        normalized = self.reproducer.normalize_dynamic_hashes(value, locations)
        for path, digest in value['authenticated_input_hashes'].items():
            actual = normalized['authenticated_input_hashes'][path]
            self.assertEqual('<validated-generated-result>' if path.startswith('results/') else digest, actual)
        for actual in (4.0, True, 5):
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.compare_portable(4, actual, path='/rank', relative_tolerance=2e-8,
                                                 absolute_tolerance=2e-6, compare_numbers=True)

    def test_tolerance_widening_and_graph_changes_cannot_override_release_spec(self):
        for mutation in ('comparison_policy', 'dependencies', 'identity_policy'):
            manifest = copy.deepcopy(self.manifest)
            row = next(item for item in manifest['steps'] if item['artifact_id'] == 'NSC-4-SHAPE-RESPONSE')
            row[mutation] = {} if mutation != 'dependencies' else []
            with self.assertRaises(self.reproducer.ReproductionError):
                self.reproducer.validate_checkout(manifest)

    def test_nsc9_through_nsc11_are_appended_terminal_all_field_records(self):
        outputs = [row['output'] for row in self.manifest['steps']]
        self.assertEqual(outputs[77:81], [
            'results/nsc-9-covariant-source.json',
            'results/nsc-10-measure-normalization.json',
            'results/nsc-10-influence.json',
            'results/nsc-11-response-matching.json',
        ])
        self.assertEqual(outputs[-4:], [
            'results/development/compact-interaction.json',
            'results/development/torsion-uv-map.json',
            'results/development/flow-compatibility.json',
            'results/development/charged-sector.json',
        ])
        self.assertEqual(85, len(self.manifest['steps']))
        self.assertEqual(27, len(self.release['scoped_follow_ups']))
        self.assertEqual('95b96be312feb667377cdbc3bbfe453697a458dd', self.release['source_commit'])
        self.assertEqual('3a747cc17e33a6a3d6cc58634eaa40dd69e30a26',
                         self.steps['NSC-9-COVARIANT-SOURCE']['follow_up_source_commit'])
        for artifact in ('NSC-9-COVARIANT-SOURCE', 'NSC-10-MEASURE-NORMALIZATION', 'NSC-10-INFLUENCE', 'NSC-11-RESPONSE-MATCHING'):
            step = self.steps[artifact]
            value = json.loads((ROOT / step['output']).read_text())
            self.reproducer.validate_identity(value, step)
            self.assertTrue(value['terminal'])
            self.assertEqual('all_fields', step['comparison_policy']['kind'])
            absolute = 2e-8 if artifact == 'NSC-11-RESPONSE-MATCHING' else 2e-7
            self.assertEqual((2e-8, absolute), (step['comparison_policy']['relative_tolerance'],
                                            step['comparison_policy']['absolute_tolerance']))
            locations = self.reproducer.validate_authenticated_inputs(ROOT, value, step)
            normalized = self.reproducer.normalize_dynamic_hashes(value, locations)
            self.assertEqual(value['source_hashes'], normalized['source_hashes'])
            for path, digest in value['input_hashes'].items():
                actual = normalized['input_hashes'][path]
                self.assertEqual('<validated-generated-result>' if path.startswith('results/') else digest, actual)

    def test_importer_rejects_root_copy_private_and_escape_paths(self):
        for path in ('.', 'README.md', '../results/x.json', '/tmp/x.json',
                     'archive/results/x.json', 'runs/raw.json', '.private/key.json',
                     'docs/.private.md', 'scripts/tool.sh'):
            with self.assertRaises(ValueError):
                self.importer.safe_relative(path)
        self.assertEqual(
            'results/development/charged-sector.json',
            self.importer.safe_relative('results/development/charged-sector.json'),
        )


if __name__ == '__main__':
    unittest.main()
