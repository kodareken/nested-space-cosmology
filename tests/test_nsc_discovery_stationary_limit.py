"""Independent finite-limit evidence checks; old producers remain unchanged."""
import json
from pathlib import Path
import tempfile

import numpy as np
import pytest

import assess_nsc_discovery_stationary_limit as limit


@pytest.fixture
def tmp_path():
    with tempfile.TemporaryDirectory(prefix='nsc-stationary-limit-test-',dir='/tmp') as path:
        yield Path(path)


@pytest.fixture(scope='module')
def report():
    return limit.assess()


def test_limit_lemma_ap_weights_copy_factor_and_length():
    coefficients=limit.coupling.locked_coefficients()
    base=limit.analytic_limits(coefficients)
    assert base['massless_base_field_energy']==pytest.approx(7*np.pi/2)
    assert base['auxiliary_energy_limit']==pytest.approx(13.400474772648874)
    assert base['source_eta_limit']==pytest.approx(1.2187153141973204)
    longer=limit.analytic_limits(coefficients,length=16.)
    assert longer['massless_base_field_energy']==pytest.approx(base['massless_base_field_energy']/2)
    assert longer['auxiliary_energy_limit']==pytest.approx(base['auxiliary_energy_limit']/2)
    assert longer['source_eta_limit']==pytest.approx(base['source_eta_limit'])
    assert longer['rQ_limit']==pytest.approx(base['rQ_limit']/2)
    assert longer['proper_period_length_limit']==pytest.approx(base['proper_period_length_limit'])
    doubled=limit.analytic_limits(coefficients,multiplicity=8)
    assert doubled['massless_base_field_energy']==pytest.approx(2*base['massless_base_field_energy'])
    assert doubled['source_eta_limit']==pytest.approx(base['source_eta_limit']/2)
    assert not base['new_spectral_trial_performed']


def test_leading_radial_and_auxiliary_signs_algebraically():
    # Resolved-mode arithmetic, no operator eigenvalue/source trial.
    x=np.arange(512)*8/512;k=2*np.pi;c=np.cos(k*x)
    a=.001;s2=a*a*k*k/6
    R1=-a*c/3
    first_radial=3*k*k*R1+a*k*k*c
    np.testing.assert_allclose(first_radial,0,atol=1e-16)
    assert s2+np.mean(a*k*k*c*R1)==pytest.approx(0,abs=1e-18)
    leading_chi=-12*c/a
    np.testing.assert_allclose(s2*leading_chi,-2*a*k*k*c,atol=1e-15)
    assert np.mean(s2*leading_chi**2)==pytest.approx(12*k*k)


def test_saved_limit_matches_but_force_and_original_replay_remain_unsatisfied(report):
    assert report['binding']['producing_commit']==limit.PRODUCING_COMMIT
    assert report['binding']['envelope_and_payload_authenticated']
    assert report['strict_original_replay']['old_checker_exit_code']==1
    assert report['strict_original_replay']['expected_failure_observed']
    assert report['strict_original_replay']['numeric_identity_unresolved']
    assert not report['strict_original_replay']['strict_original_numeric_replay_passed']
    match=report['selected_saved_match'];theory=report['analytic_limits']
    assert match['s_over_abs_a']==pytest.approx(theory['s_over_abs_a_limit'],rel=1e-6)
    assert match['r_times_abs_a']==pytest.approx(theory['r_times_abs_a_limit'],rel=2e-6)
    assert match['rQ_mean']==pytest.approx(theory['rQ_limit'],rel=2e-6)
    assert match['eta']==pytest.approx(theory['source_eta_limit'],rel=2e-6)
    assert report['actual_force']['actual_retained_pQ_max']==pytest.approx(32.95694394768389,abs=1e-8)
    assert report['actual_force']['actual_retained_pQ_vs_saved_max_gap']==0
    assert not report['actual_force']['full_force_satisfied']
    assert not report['scope']['finite_branch_nonexistence_proven']
    assert not report['scope']['holding_or_stability_claim']
    assert report['cpu_seconds']<=10


def test_stable_einstein_extraction_preserves_original_action_and_force_failure(report):
    affine=report['einstein_affine_conditioning']
    assert affine['well_conditioned_initial_stable_vs_ordinary_max']<1e-9
    assert affine['selected_stable_vs_stored_ordinary_max']>1e-8
    assert affine['original_normal_denominator']>0
    assert affine['stable_normal_denominator']>0
    assert affine['stable_denominator_identity_defect']==pytest.approx(0,abs=1e-10)
    assert affine['ordinary_affine_at_actual_radius_gap_max']>.2
    assert max(affine['ordinary_affine_terms_max'])>1e7
    assert affine['coefficient_isolation']['isolated_source_force_max']==0
    assert not affine['coefficient_isolation']['original_action_coefficients_changed']
    assert not affine['actual_candidate_replaced_with_stable_projection']
    assert not affine['indicators_are_certified_bounds']


def test_reordering_failure_is_recorded_with_indicators_not_relabelled(report):
    strict=report['strict_original_replay']
    assert strict['producer_order_gap_max']==0
    assert strict['checker_order_gap_max']>strict['historical_absolute_gradient_threshold']
    assert strict['checker_order_gap_max']==pytest.approx(2.0510802431483732e-10,rel=.1)
    assert len(strict['eps_weighted_dot_indicators'])==7
    assert all(value>=0 for value in strict['eps_weighted_dot_indicators'])
    assert not strict['indicator_is_certified_bound']
    assert strict['failure_message_not_replaced_with_success']


def test_creation_only_check_preserves_failure_outcome(report,tmp_path):
    path=tmp_path/'assessment.json'
    limit.write_record(report,path)
    old=path.read_bytes()
    checked=limit.check_record(path)
    assert checked['assessment_record_authenticated']
    assert not checked['strict_original_numeric_replay_passed']
    assert checked['numeric_identity_unresolved']
    assert checked['bytes_written']==0 and not checked['numerical_model_evaluated']
    assert path.read_bytes()==old
    with pytest.raises(FileExistsError):limit.write_record(report,path)
    changed=json.loads(path.read_text())
    changed['strict_original_replay']['old_checker_exit_code']=0
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError,match='erased'):limit.check_record(path)
    with pytest.raises(ValueError):limit.output_path(limit.INPUT)


def test_evidence_and_old_checker_bytes_unchanged_after_consumption(report):
    for name,digest in report['input_hashes'].items():
        assert limit.sha256(limit.repository_path(name))==digest
    for name,digest in report['source_hashes'].items():
        assert limit.sha256(limit.repository_path(name))==digest
    assert not report['scope']['old_evidence_or_checker_modified']


def test_consumer_does_not_call_eigensolvers_or_initial_radius_solver(monkeypatch):
    import scipy.linalg
    def forbidden(*args, **kwargs):
        raise AssertionError('no new source eigensolve or geometry solve is authorized')
    monkeypatch.setattr(np.linalg,'eigh',forbidden)
    monkeypatch.setattr(np.linalg,'eigvalsh',forbidden)
    monkeypatch.setattr(scipy.linalg,'eigh',forbidden)
    monkeypatch.setattr(limit.galerkin,'solve_initial_radius',forbidden)
    measured=limit.assess()
    assert measured['actual_force']['new_eigenstate_or_source_trial'] is False
    assert measured['cpu_seconds']<=limit.CPU_BUDGET
