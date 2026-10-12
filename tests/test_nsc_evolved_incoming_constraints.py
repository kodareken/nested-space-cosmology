import json
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import block_diag

from recursive_horizons.nsc_evolved_incoming_constraints import (
    source_column_matter,surface_geometry_response,evaluate_source_fixed_history)
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_prepared_history_jets import restriction_covariance_tangent
from recursive_horizons.nsc_common_subtracted_ks_source import (
    physical_source_kernel,raw_ks_vertex_coefficients,SpectrumDeclaration)
from recursive_horizons.nsc_transmitting_dirac_domain import S1, S2, S3


ROOT = Path(__file__).resolve().parents[1]
PARAM = dict(mass=np.pi/2, angular=np.sqrt(5), axial_scale=np.sqrt(3*np.pi/2-4), radius=np.sqrt(2), multiplicity=7.)


def columns():
    rng = np.random.default_rng(240918)
    f = (rng.normal(size=(3,2,6))+1j*rng.normal(size=(3,2,6)))/7
    fz = (rng.normal(size=f.shape)+1j*rng.normal(size=f.shape))/5
    df = (rng.normal(size=(2,*f.shape))+1j*rng.normal(size=(2,*f.shape)))/9
    dfz = (rng.normal(size=df.shape)+1j*rng.normal(size=df.shape))/8
    C = block_diag(*(source_covariance(E,.238,3.97,np.pi/2) for E in (.4,1.3)))
    return f,fz,C,df,dfz


def test_stationary_limit_has_owned_raw_vertex_signs_and_no_extra_measure():
    f,_,C,df,_ = columns()
    E = np.repeat([.4,1.3],3)
    fz = -1j*E*f
    dfz = -1j*E*df
    result = source_column_matter(f,fz,C,df,dfz,**PARAM)
    expected = np.zeros((len(f),2),complex)
    for j,energy in enumerate((.4,1.3)):
        block = slice(3*j,3*j+3)
        rho = f[:,:,block]@C[block,block]@f[:,:,block].swapaxes(-1,-2).conj()
        VN = -PARAM['mass']*S1+PARAM['angular']/PARAM['radius']*S2-energy/PARAM['axial_scale']*S3
        expected[:,0] -= PARAM['multiplicity']*np.einsum('zij,ji->z',rho,VN)
        expected[:,1] -= PARAM['multiplicity']*energy*np.trace(rho,axis1=-2,axis2=-1)
    np.testing.assert_allclose(result['action_gradient'],expected.real,atol=3e-13,rtol=0)


def test_retarded_source_tangent_matches_full_coherent_restriction_and_difference():
    f,fz,C,df,dfz = columns()
    result = source_column_matter(f,fz,C,df,dfz,**PARAM)
    for direction in range(2):
        full = restriction_covariance_tangent(f.reshape(-1,6),df[direction].reshape(-1,6),C,np.zeros_like(C))
        diagonal = np.array([full[2*z:2*z+2,2*z:2*z+2] for z in range(len(f))])
        np.testing.assert_allclose(result['density_tangent'][direction],diagonal,atol=3e-15,rtol=0)
        h = 2e-5
        plus = source_column_matter(f+h*df[direction],fz+h*dfz[direction],C,df,dfz,**PARAM)
        minus = source_column_matter(f-h*df[direction],fz-h*dfz[direction],C,df,dfz,**PARAM)
        np.testing.assert_allclose((plus['action_gradient']-minus['action_gradient'])/(2*h),
                                   result['action_gradient_tangent'][direction],atol=3e-10,rtol=0)
    diagonal_source = source_column_matter(f,fz,np.diag(np.diag(C)),df,dfz,**PARAM)
    assert np.max(abs(result['action_gradient_tangent']-diagonal_source['action_gradient_tangent'])) > 1e-6


def test_nonstationary_momentum_cannot_be_replaced_by_source_energy():
    f,fz,C,df,dfz = columns()
    actual = source_column_matter(f,fz,C,df,dfz,**PARAM)
    energy = np.repeat([.4,1.3],3)
    frozen_label = source_column_matter(f,-1j*energy*f,C,df,-1j*energy*df,**PARAM)
    assert np.max(abs(actual['action_gradient']-frozen_label['action_gradient'])) > .01
    assert np.max(abs(actual['action_gradient_tangent']-frozen_label['action_gradient_tangent'])) > .01


def test_nonstationary_value_matches_owned_near_diagonal_source_kernel():
    f,fz,C,df,dfz = columns()
    weights = np.array([.3,.7])
    rootw = np.repeat(np.sqrt(weights/(2*np.pi)),3)
    actual = source_column_matter(f*rootw,fz*rootw,C,df*rootw,dfz*rootw,**PARAM)
    unflatten = lambda v: v.reshape(len(v),2,2,3).transpose(0,2,1,3)
    fields, derivatives = unflatten(f),unflatten(fz)
    kernel = physical_source_kernel(fields,fields,derivatives,derivatives,
        np.stack((C[:3,:3],C[3:,3:])),np.array([.4,1.3]),weights,separation=0.,
        declaration=SpectrumDeclaration('control_only','independent finite-column contraction test'))
    v = raw_ks_vertex_coefficients(np.array([1.,0.,PARAM['axial_scale'],PARAM['radius']]),
                                  PARAM['mass'],PARAM['angular'],envelopes=np.ones(4))
    expected = -PARAM['multiplicity']*(np.einsum('bij,zji->zb',v.multiplication[:2],kernel.kernel)
        -1j*np.einsum('bij,zji->zb',v.momentum[:2],kernel.separation_derivative))
    np.testing.assert_allclose(actual['action_gradient'],expected.real,atol=3e-14,rtol=0)


def owned_coefficients():
    coefficients = json.loads((ROOT/'results/development/nsc-incoming-surface-coefficients.json').read_text())['coefficients']
    branch = json.loads((ROOT/'results/development/nsc-incoming-surface-regular-branch.json').read_text())['branch']
    mid = lambda x: float(np.mean(x))
    return {'a':mid(branch['intrinsic_intervals']['a']), 'r':mid(branch['intrinsic_intervals']['r']),
            'Hr':mid(branch['intrinsic_intervals']['Hr']), 'A0':mid(branch['imported_unchanged']['A0']),
            'A1':mid(branch['intervals']['A1_total']), 'd':mid(branch['imported_unchanged']['d']),
            'e':mid(branch['imported_unchanged']['e']), 'C':mid(coefficients['C_interval']),
            'D':np.r_[0.,[mid(coefficients['polynomial_intervals']['D'][str(j)]) for j in range(1,5)]],
            'F':np.array([mid(coefficients['polynomial_intervals']['F'][str(j)]) for j in range(3)])}


def test_geometry_response_keeps_old_constants_out_and_has_complete_derivative():
    c = owned_coefficients()
    z = np.array([-.1,.2])
    basis = np.stack((z**3,3*z*z,6*z,np.full_like(z,6),1+z,np.ones_like(z)),axis=-1)
    slots = .02*basis
    tangent = basis[None]
    result = surface_geometry_response(slots,tangent,c)
    h = 1e-6
    plus = surface_geometry_response(slots+h*basis,tangent,c)['action_gradient_change']
    minus = surface_geometry_response(slots-h*basis,tangent,c)['action_gradient_change']
    np.testing.assert_allclose((plus-minus)/(2*h),result['action_gradient_tangent'][0],atol=3e-10,rtol=0)
    zero = surface_geometry_response(np.zeros_like(slots),tangent,c)
    np.testing.assert_array_equal(zero['action_gradient_change'],0)
    bad = {**c,'D':c['D'].copy()};bad['D'][0]=1.
    with pytest.raises(ValueError,match=r'D\[0\]=0'):
        surface_geometry_response(slots,tangent,bad)


def test_fresh_history_entry_delegates_owned_evolution_with_fixed_source_inputs(monkeypatch):
    spec=importlib.util.spec_from_file_location('evolved_constraint_control',ROOT/'scripts/derive_nsc_evolved_incoming_constraints.py')
    control=importlib.util.module_from_spec(spec);spec.loader.exec_module(control)
    records,arrays,meta,channel,source,coeff=control.load_inputs()
    state=control.restore_saved_control(arrays,meta,source,'center',records[0]['control']['amplitude'])
    x,phi,data,_=control.STATE.RETARDED.inputs()
    owner=control.STATE.FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item())
    provider=control.STATE.CachedCompatibleIncomingMetric(x,control.STATE.family(.001,meta['support']))
    _,expected_tangent,expected_incoming=control.STATE.incident(owner,phi,data,1)
    calls=[]
    def checked_factory(given_owner,times,given_provider,fields,tangents,incoming,given_source):
        assert given_owner is owner and given_provider is provider and given_source is source
        np.testing.assert_array_equal(fields,phi)
        np.testing.assert_array_equal(tangents,expected_tangent)
        for time in (times[:-1]+times[1:])/2:
            actual,delta=incoming(time);expected,dexpected=expected_incoming(time)
            np.testing.assert_array_equal(actual,expected)
            np.testing.assert_array_equal(delta,dexpected)
        calls.append(True)
        return state  # authenticated output of this same actual evolution, no new run
    import recursive_horizons.nsc_evolved_incoming_state as owner_module
    monkeypatch.setattr(owner_module,'evolve_incoming_state',checked_factory)
    Dz=control.make_interp_spline(state.z,np.eye(len(state.z)),k=3).derivative()(state.z)
    returned,result=evaluate_source_fixed_history(owner,arrays['times'],provider,phi,source,Dz,
        records[1]['baseline']['action_gradient_approximant'],coeff,channel,1)
    assert returned is state and calls==[True]
    assert np.max(abs(result['evolved_matter_tangent']))>1e-5
    assert result['reference_state_tangent_subtracted'] is False
    assert result['physical_constraint_status']=='OPEN'
    with pytest.raises(TypeError,match='not incoming C0'):
        evaluate_source_fixed_history(owner,arrays['times'],provider,phi,np.eye(2)/2,Dz,
            records[1]['baseline']['action_gradient_approximant'],coeff,channel,1)


def test_integrated_receipt_resolves_state_terms_without_claiming_a_physical_root():
    record=json.loads((ROOT/'results/development/nsc-evolved-incoming-constraints.json').read_text())
    error=record['residuals']['complete_discrete_constraint_derivative']
    assert error<record['control']['derivative_control_tolerance']
    assert min(record['omitted_term_effects'].values())>10*error
    assert record['frozen_inputs_rejected']==['bare_C0','legacy_matter']
    assert record['scope']['physical_constraint_status']=='OPEN'
    assert record['scope']['full_constraint_residual_certified'] is False
    assert record['error_budget']['changed_history_complete_source_error'] is None
