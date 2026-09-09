"""Condition a cancellation diagnostic while preserving its raw numerical tests."""
import copy
import json
import unittest

from test_publication import ROOT,load_script


class DerivativeResidualPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reproducer=load_script('reproduce_public_results.py')
        cls.record=json.loads((ROOT/'results/development/spectral-endpoint.json').read_text())
        release=json.loads((ROOT/'results/release-spec.json').read_text())
        cls.policy=next(r['comparison_policy'] for r in release['scoped_follow_ups']
                        if r['output']=='results/development/spectral-endpoint.json')

    def compare(self,actual):
        overrides=self.reproducer.derived_numeric_overrides(self.record,actual,self.policy)
        self.reproducer.compare_portable(self.record,actual,path='',compare_numbers=True,
            absolute_tolerance=self.policy['absolute_tolerance'],
            relative_tolerance=self.policy['relative_tolerance'],numeric_overrides=overrides)
        return overrides

    @staticmethod
    def perturb_raw_derivatives(record,change):
        row=record['numerical_controls'][3]
        row['central_difference_steps']=[x+change for x in row['central_difference_steps']]
        a,b=row['central_difference_steps']
        row['extrapolated_log_derivative']=(4*b-a)/3
        row['independent_derivative_error']=row['extrapolated_log_derivative']-row['compact_18_26_34'][2]['d_log_radius']

    def test_small_raw_change_propagates_to_the_residual(self):
        actual=copy.deepcopy(self.record)
        # Synthetic conditioning control of the magnitude observed in CI;
        # this is not substituted for an actual scientific output.
        self.perturb_raw_derivatives(actual,3e-7)
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.compare_portable(self.record,actual,path='',compare_numbers=True,
                absolute_tolerance=3e-9,relative_tolerance=3e-8)
        overrides=self.compare(actual)
        self.assertEqual({f'/numerical_controls/{i}/independent_derivative_error' for i in range(4)},set(overrides))

    def test_false_residual_is_rejected_even_inside_the_old_absolute_tolerance(self):
        actual=copy.deepcopy(self.record)
        actual['numerical_controls'][3]['independent_derivative_error']+=1e-10
        with self.assertRaisesRegex(self.reproducer.ReproductionError,'residual identity'):
            self.compare(actual)

    def test_false_richardson_extrapolation_is_rejected(self):
        actual=copy.deepcopy(self.record);row=actual['numerical_controls'][3]
        row['extrapolated_log_derivative']+=3e-7
        row['independent_derivative_error']=row['extrapolated_log_derivative']-row['compact_18_26_34'][2]['d_log_radius']
        with self.assertRaisesRegex(self.reproducer.ReproductionError,'Richardson identity'):
            self.compare(actual)

    def test_original_derivative_accuracy_is_still_required(self):
        actual=copy.deepcopy(self.record);self.perturb_raw_derivatives(actual,2e-4)
        with self.assertRaisesRegex(self.reproducer.ReproductionError,'derivative accuracy'):
            self.compare(actual)

    def test_no_raw_physical_quantity_gets_an_override(self):
        actual=copy.deepcopy(self.record)
        actual['endpoint_comparisons'][0]['full']['d_log_radius']+=1e-5
        with self.assertRaises(self.reproducer.ReproductionError):self.compare(actual)

    def test_unchanged_source_contract_and_nonfinite_rejection(self):
        self.compare(self.record)
        bad=copy.deepcopy(self.policy)
        bad['derived_quantities'][0]['derivative_relative_accuracy']=3e-6
        with self.assertRaises(self.reproducer.ReproductionError):
            self.reproducer.derived_numeric_overrides(self.record,self.record,bad)
        actual=copy.deepcopy(self.record)
        actual['numerical_controls'][3]['central_difference_steps'][0]=float('nan')
        with self.assertRaises(self.reproducer.ReproductionError):self.compare(actual)


if __name__=='__main__':unittest.main()
