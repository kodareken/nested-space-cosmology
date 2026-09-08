"""Conditioning controls for the smooth-profile convergence-order diagnostic."""
from __future__ import annotations

import copy
import json
import math
import unittest

from test_publication import ROOT, load_script


class DerivedPortabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer = load_script('reproduce_public_results.py')
        cls.expected = json.loads((ROOT / 'results/nsc-4-smooth-geometry.json').read_text())
        manifest = cls.reproducer.load_manifest()
        cls.policy = next(row['comparison_policy'] for row in manifest['steps']
                          if row['artifact_id'] == 'NSC-4-SMOOTH-GEOMETRY')

    def compare(self, actual, *, expected=None):
        expected = self.expected if expected is None else expected
        overrides = self.reproducer.derived_numeric_overrides(expected, actual, self.policy)
        self.reproducer.compare_portable(expected, actual, path='',
            relative_tolerance=self.policy['relative_tolerance'],
            absolute_tolerance=self.policy['absolute_tolerance'],
            compare_numbers=True, numeric_overrides=overrides)
        return overrides

    def consistent_error_change(self, changed_error):
        actual = copy.deepcopy(self.expected)
        family = actual['gap_families'][0]
        family['lattice'][2]['minimum_sampled_bloch_gap'] = family['continuum']['band_edge'] + changed_error
        family['absolute_errors'] = [abs(row['minimum_sampled_bloch_gap'] - family['continuum']['band_edge'])
                                     for row in family['lattice']]
        errors = family['absolute_errors']
        family['observed_orders'] = [math.log2(left / right) for left, right in zip(errors, errors[1:])]
        return actual

    def test_consistent_linux_scale_roundoff_is_accepted_with_exact_scope(self):
        # Reconstruct a gap perturbation at the scale reported by Linux CI.
        observed_linux_order = 1.9996935208105302
        family = self.expected['gap_families'][0]
        changed_error = family['absolute_errors'][1] / 2**observed_linux_order
        actual = self.consistent_error_change(changed_error)
        self.assertGreater(abs(actual['gap_families'][0]['observed_orders'][1]
                               - family['observed_orders'][1]), 3e-8)
        self.assertLess(abs(actual['gap_families'][0]['absolute_errors'][2]
                            - family['absolute_errors'][2]), 2e-12)
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.compare_portable(self.expected, actual, path='',
                relative_tolerance=1e-8, absolute_tolerance=1e-8, compare_numbers=True)
        overrides = self.compare(actual)
        allowed = {f'/gap_families/{family}/{field}/{index}'
                   for family in range(3)
                   for field, count in (('absolute_errors', 3), ('observed_orders', 2))
                   for index in range(count)}
        self.assertEqual(allowed, set(overrides))
        for pointer, tolerances in overrides.items():
            if '/absolute_errors/' in pointer:
                self.assertEqual((0.0, 2e-12), tolerances)

    def test_order_corruption_is_rejected_even_inside_the_propagated_interval(self):
        actual = copy.deepcopy(self.expected)
        actual['gap_families'][0]['observed_orders'][1] += 1e-7
        with self.assertRaisesRegex(self.reproducer.ReproductionError, 'order does not match'):
            self.compare(actual)

    def test_raw_error_beyond_budget_is_rejected_despite_consistent_order(self):
        error = self.expected['gap_families'][0]['absolute_errors'][2]
        actual = self.consistent_error_change(error + 3e-12)
        with self.assertRaisesRegex(self.reproducer.ReproductionError, 'raw error exceeds'):
            self.compare(actual)

    def test_gap_and_error_identity_must_hold_in_each_record(self):
        actual = copy.deepcopy(self.expected)
        actual['gap_families'][0]['lattice'][2]['minimum_sampled_bloch_gap'] += 1e-13
        with self.assertRaisesRegex(self.reproducer.ReproductionError, 'actual error does not match'):
            self.compare(actual)
        with self.assertRaisesRegex(self.reproducer.ReproductionError, 'expected error does not match'):
            self.compare(self.expected, expected=actual)

    def test_scientific_flags_and_unrelated_numbers_stay_strict(self):
        actual = copy.deepcopy(self.expected)
        actual['nonclaims']['physical_metric_stability_proved'] = True
        with self.assertRaises(self.reproducer.ReproductionError):
            self.compare(actual)
        actual = copy.deepcopy(self.expected)
        actual['gap_families'][0]['radius'] += 1e-4
        with self.assertRaises(self.reproducer.ReproductionError):
            self.compare(actual)

    def test_unresolved_nonfinite_and_boolean_errors_are_rejected(self):
        for error in (2e-12, 1e-12, float('nan'), float('inf'), True):
            actual = copy.deepcopy(self.expected)
            actual['gap_families'][0]['absolute_errors'][2] = error
            with self.assertRaises(self.reproducer.ReproductionError):
                self.compare(actual)

    def test_identity_and_interval_rounding_allowances_are_compatible(self):
        error = self.expected['gap_families'][0]['absolute_errors'][2]
        actual = self.consistent_error_change(error + 1.99e-12)
        order = actual['gap_families'][0]['observed_orders'][1]
        actual['gap_families'][0]['observed_orders'][1] += 8 * math.ulp(order)
        self.compare(actual)


if __name__ == '__main__':
    unittest.main()
