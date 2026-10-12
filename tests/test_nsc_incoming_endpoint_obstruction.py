"""Finite operator identities, actual guards and analytic-proof receipt scope."""
import importlib.util
import json
from pathlib import Path
import pytest
from fractions import Fraction
import recursive_horizons.nsc_incoming_endpoint_obstruction as M
ROOT=Path(__file__).resolve().parents[1]


def script():
    spec=importlib.util.spec_from_file_location('endpoint_proof',ROOT/'scripts/derive_nsc_incoming_endpoint_obstruction.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_scalar_projector_recurrence_and_endpoint_sign():
    proof=M.symbolic_proof()
    assert proof['checked_scalar_residual_count']==40
    assert proof['scalar']['energy_degree_of_G']==1
    assert proof['recurrence']['orders_0_through_3_residuals']==['0']*4
    assert set(proof['recurrence']['remaining_coefficients'])==set(map(str,range(4,9)))
    assert proof['endpoint']['original_limit']=='-W0*a1**2*ell/(4*r1**2)'
    assert proof['endpoint']['two_mean_necessary_conditions']==[{'Amean':'0','Rmean':'0'}]
    assert not proof['flat_shift_zero_transfer']['nonzero_transfer_claimed']
    assert proof['linearized_bridge']['physical_linearized_error_order'].startswith('X-deltaP4=O(E^-4)')


def test_actual_positive_bump_and_nonzero_mass_angular_guards():
    module=script();records,_,_=module.load_inputs()
    _,compensator,pilot,inventory,restart=records
    channel=inventory['channels'][14];config=restart['scattering_provenance']['config'];pulse=pilot['control']['support']
    guard=M.rational_guard(channel,pulse,config,compensator)
    upper=guard['limit_upper_strict']
    assert Fraction(upper['numerator'],upper['denominator'])==-Fraction(2,375)
    with pytest.raises(ValueError,match='angular'):
        M.rational_guard(dict(channel,angular_eigenvalue=0.),pulse,config,compensator)
    with pytest.raises(ValueError,match='positive compact'):
        M.rational_guard(channel,dict(pulse,axial_halfwidth=.01),config,compensator)


def test_receipt_replay_never_calls_physical_or_old_tail_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('physical/scientific producer entered algebraic proof')
    import scipy.integrate
    import recursive_horizons.nsc_retarded_radial_response as radial
    import recursive_horizons.nsc_incoming_vacuum_tail_bound as tail
    import recursive_horizons.nsc_upstream_compensator as comp
    monkeypatch.setattr(scipy.integrate,'solve_ivp',forbidden)
    monkeypatch.setattr(radial,'short_response',forbidden)
    monkeypatch.setattr(tail.IncomingVacuumTailBoundGeometry,'__init__',forbidden)
    monkeypatch.setattr(comp,'operator_response',forbidden)
    module=script();record=module.make_record()
    assert record==json.loads((ROOT/module.OUTPUT).read_text())
    assert record['conclusion']['finite_energy_threshold'] is None
    assert record['conclusion']['numerical_Big_O_constant'] is None
    assert not record['conclusion']['fixed_C0_linear_matching_in_this_class']
    assert record['conclusion']['finite_pair_calibration_can_still_hold']
    assert record['scope']['field_radial_source_or_old_producer_runs']==0
    assert not record['scope']['undeclared_Planck_cutoff_assumed']
    assert not record['scope']['broader_NSC_architecture_excluded']
    assert record['analytic_uniform_proof']['uniform_orders_are_analytic_not_numerically_estimated']
