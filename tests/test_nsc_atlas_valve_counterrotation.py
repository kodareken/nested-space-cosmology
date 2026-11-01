"""Public same-action counterrotation data, constraints and authenticated load."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import recursive_horizons.nsc_atlas_valve_counterrotation as counterrotation
from recursive_horizons.nsc_atlas_kick_york import tensor_divergence
from recursive_horizons.nsc_atlas_s3_frame import arbitrary_points_values
from recursive_horizons.nsc_atlas_s3_spectral import harmonic_values
from recursive_horizons.nsc_atlas_valve_counterrotation import (
    CounterrotatingData,original_points,rotated_points,solve_counterrotating_initial_data,
)
from recursive_horizons.nsc_atlas_valve_neck import solve_necks


@pytest.fixture(scope='module')
def finite_data():
    # One actual lambda1 continuation and independent full validation,
    # shared by all tests; no dependence on any archived raw cache.
    return solve_counterrotating_initial_data()


def _points():
    points=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
                     [.5,.5,.5,.5],[.3,-.4,.5,.7],[-.2,.6,.3,-.7]],dtype=float)
    return points/np.linalg.norm(points,axis=1)[:,None]


def test_zero_kick_recovers_one_mode_neck_without_legacy_source():
    zero=solve_counterrotating_initial_data(max_j=4,lambda_path=(0.,))
    neck=solve_necks(40,(1.,))[0];points=_points()
    np.testing.assert_allclose(zero.values(points),neck.values(points[:,0]),atol=4e-14,rtol=0)
    np.testing.assert_allclose(zero.chi(points),zero.A*points[:,0],atol=0,rtol=0)
    assert zero.kick_amplitude==0 and zero.metadata['coefficient_change_from_neck']==0
    accounting=zero.metadata['stages'][0]['accounting']
    assert accounting['normal_Pi_energy']==accounting['integrated_KijKij']==0
    assert zero.chi(np.array([[0.,2**-.5,0.,2**-.5]]))[0]==0


def test_new_source_and_york_divergence_match_analytic_single_mode(finite_data):
    data=finite_data;points=_points();x,y,z,w=points.T;p=data.parent
    B=np.sqrt(12)*p.H**1.5/np.pi
    np.testing.assert_allclose(data.chi(points),data.A*x,atol=0,rtol=0)
    np.testing.assert_allclose(data.f(points),B*y*w,atol=1e-18,rtol=3e-16)
    np.testing.assert_allclose(original_points(rotated_points(points)),points,atol=0,rtol=0)
    basis=data.york_basis
    div_coeff=tensor_divergence(basis,data.tensor_coefficients,r=p.H)
    divergence=arbitrary_points_values(points,data.tensor_labels)@div_coeff.T
    analytic=-p.kappa*data.f(points)[:,None]*p.H*data.A*np.column_stack([w,-z,y])
    np.testing.assert_allclose(divergence,analytic,atol=7e-18,rtol=1e-13)
    tensor=data.tensor(points)
    np.testing.assert_allclose(tensor,tensor.swapaxes(1,2),atol=2e-18,rtol=0)
    np.testing.assert_allclose(np.trace(tensor,axis1=1,axis2=2),0,atol=2e-18)
    # Consequential explicit component of the independently derived cubic seed.
    normalized_A00=(13*w*w*z/60-2*w*x*y/15+11*x*x*z/60)
    np.testing.assert_allclose(tensor[:,0,0],p.kappa*data.A*B*normalized_A00,atol=7e-18,rtol=1e-12)


def test_lambda_one_full_constraints_positive_bound_and_energy(finite_data):
    data=finite_data;metadata=data.metadata;validation=metadata['independent_validation']
    assert data.kick_amplitude==1 and metadata['branch_completed']
    assert metadata['stages'][-1]['analytic_finite_polynomial_psi_lower_bound']>.45
    assert validation['sampled_psi_min']>.69
    assert validation['physical_Hamiltonian_RMS_over_H_squared']<1e-4
    assert validation['physical_momentum_RMS_over_H_squared']<2e-13
    assert validation['full_conformal_residual_RMS']>1e3*metadata['stages'][-1]['projected_residual_norm']
    # Independent local reconstruction includes the correct scalar source,
    # full vacuum and inverse-power shear/normal-momentum contribution.
    points=_points();Y=harmonic_values(rotated_points(points),data.labels)
    eigen=np.array([l*(l+2) for l,_,_ in data.labels]);psi=data.values(points)
    lap=-Y@(eigen*data.coefficients);p=data.parent;z=p.mass_squared/p.H**2
    a=p.kappa*data.A**2;u=points[:,0];tensor2=data.squared_norm(points);f=data.f(points)
    residual=-8*lap+(6-a*(1-u*u))*psi-(6+a*z*u*u)*psi**5-(tensor2+p.kappa*f*f)/p.H**2*psi**-7
    assert np.max(abs(residual))<1.5e-4
    accounting=validation['accounting'];stage=metadata['stages'][-1]['accounting']
    assert accounting['normal_Pi_energy']>0 and accounting['integrated_KijKij']>0
    assert accounting['total_proper_matter_energy']==pytest.approx(sum(accounting[key] for key in (
        'normal_Pi_energy','spatial_gradient_energy','scalar_mass_energy','vacuum_energy')),rel=2e-14)
    assert accounting['shear_constraint_energy']==pytest.approx(accounting['integrated_KijKij']/(2*p.kappa),rel=2e-14)
    for key,value in accounting.items():
        assert value==pytest.approx(stage[key],rel=3e-12,abs=1e-13)
    assert not data.coefficients.flags.writeable and not data.tensor_coefficients.flags.writeable


def test_signed_initial_tensor_momentum_and_countercirculation(finite_data):
    data=finite_data;points=_points();psi=data.values(points);Abar=data.tensor(points)
    Pi=data.kick_amplitude*data.f(points)*psi**-6
    Kcov=data.kick_amplitude*psi[:,None,None]**-2*Abar
    Korth=Kcov*psi[:,None,None]**-4
    np.testing.assert_allclose(Korth,data.kick_amplitude*psi[:,None,None]**-6*Abar,atol=2e-18)
    assert np.max(abs(Pi))>0 and np.max(abs(Korth))>0
    # The public family deliberately admits lambda>=0. Sign reversal below
    # is initial-data algebra at the SAME psi, not a new negative solver run.
    minusPi=-data.kick_amplitude*data.f(points)*psi**-6
    minusK=-data.kick_amplitude*psi[:,None,None]**-2*Abar
    np.testing.assert_array_equal(minusPi,-Pi);np.testing.assert_array_equal(minusK,-Kcov)
    np.testing.assert_array_equal(minusPi**2,Pi**2)
    np.testing.assert_array_equal(np.sum(minusK**2,axis=(1,2)),np.sum(Kcov**2,axis=(1,2)))
    angle=2*np.pi*np.arange(64)/64;radius=2**-.5
    integrals=[]
    for sign in (1,-1):
        circle=np.column_stack([radius*np.cos(angle),radius*np.sin(angle),np.zeros_like(angle),
                                np.full_like(angle,sign*radius)])
        dx=-radius*np.sin(angle)
        integrand=-data.kick_amplitude*data.values(circle)**-6*data.f(circle)*data.A*dx
        integrals.append(float(2*np.pi*np.mean(integrand)))
    assert integrals[0]>0 and integrals[1]<0
    np.testing.assert_allclose(integrals[0],-integrals[1],rtol=2e-13,atol=1e-16)


def _write_authenticated_record(directory,data):
    # A fresh temporary declared root contains the current owner-source bytes
    # and newly generated state cache. This tests authentication without
    # relying on .build files, network, Git or archived producer output in CI.
    source_name='src/recursive_horizons/nsc_atlas_valve_counterrotation.py'
    owner=directory/source_name;owner.parent.mkdir(parents=True)
    owner.write_bytes(Path(counterrotation.__file__).read_bytes())
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    pins={source_name:sha(owner)};metadata=data.metadata;raw=directory/'initial-data.npz'
    np.savez_compressed(raw,coefficients=data.coefficients,labels=np.asarray(data.labels),A=data.A,
        kick_amplitude=data.kick_amplitude,tensor_coefficients=data.tensor_coefficients,
        tensor_lrs=np.array([label[:3] for label in data.tensor_labels]),
        tensor_kind=np.array([label[3] for label in data.tensor_labels]),
        vacuum=data.parent.vacuum,mass_squared=data.parent.mass_squared,
        offset=data.parent.offset,kappa=data.parent.kappa,
        metadata=json.dumps(dict(preparation_id='one-mode-counterrotation-v1',source_sha256=pins,data_metadata=metadata)))
    record=dict(schema='nsc-atlas-valve-counterrotation-v1',repository_root_relative_to_record='.',
        provenance=dict(source_sha256_before=pins,source_sha256_after=pins,source_unchanged=True),
        raw_cache=dict(path=raw.name,sha256=sha(raw)),metadata=metadata,
        preparation=dict(A=data.A,kick_amplitude=data.kick_amplitude))
    output=directory/'record.json';output.write_text(json.dumps(record),encoding='utf-8')
    return output,raw,owner,record


def test_authenticated_temporary_roundtrip_keeps_actual_fields(tmp_path,finite_data):
    output,_,_,_=_write_authenticated_record(tmp_path,finite_data)
    loaded=CounterrotatingData.from_record(output);points=_points()
    for name in ('values','tensor','f','chi'):
        np.testing.assert_array_equal(getattr(loaded,name)(points),getattr(finite_data,name)(points))
    np.testing.assert_array_equal(loaded.coefficients,finite_data.coefficients)
    np.testing.assert_array_equal(loaded.tensor_coefficients,finite_data.tensor_coefficients)
    assert loaded.metadata==finite_data.metadata


@pytest.mark.parametrize('target',['source','cache','record_preparation'])
def test_authenticated_load_rejects_changed_source_or_state(tmp_path,finite_data,target):
    output,raw,owner,record=_write_authenticated_record(tmp_path,finite_data)
    if target=='source':owner.write_bytes(owner.read_bytes()+b'\n# changed\n')
    elif target=='cache':raw.write_bytes(raw.read_bytes()+b'changed')
    else:
        record['preparation']['kick_amplitude']=0.
        output.write_text(json.dumps(record),encoding='utf-8')
    with pytest.raises(ValueError,match='SHA mismatch|parameter mismatch'):
        CounterrotatingData.from_record(output)
