"""Same-state continuation, physical receiver, flux and composition controls."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_handoff import (
    Diamond,C,diamond_advance,discrete_symplectic_balance,
    inherited_reference,inherited_second_leg,geometry,receiver_record,calculation,flux_balance)
from recursive_horizons.nsc_atlas_horizon_response import solve_source


@pytest.fixture(scope='module')
def result():return calculation()


def test_zero_potential_recovers_two_arbitrary_null_traces():
    u=np.linspace(0,1,41);v=np.linspace(0,1,41)
    f=np.c_[np.sin(u),u*u];g=np.c_[np.sin(2*v),v*v*v]
    solution=diamond_advance(f,g,np.zeros((40,40)),1/40**2)
    exact=f[:,None,:]+g[None,:,:]
    np.testing.assert_allclose(solution,exact,atol=2e-14,rtol=0)


def test_null_corner_compatibility_and_zero_source_balance():
    f=np.zeros((9,2));g=f.copy();g[0,0]=.1
    with pytest.raises(ValueError,match='corner'):
        diamond_advance(f,g,np.ones((8,8)),.01)
    zeros=diamond_advance(f,f,np.ones((8,8)),.01)
    assert np.max(abs(zeros))==0
    assert discrete_symplectic_balance(zeros)['relative_residual']==0
    axes=np.linspace(0,1,9)
    assert flux_balance(axes,axes,zeros[...,0],zeros[...,1])['relative_residual']==0


def test_same_potential_discrete_green_identity_is_exact():
    u=np.linspace(0,.2,33);v=np.linspace(0,.2,33)
    f=np.c_[np.cos(7*u),np.cos(5*u)]
    g=np.c_[np.cos(5*v),np.cos(7*v)]
    U,V=np.meshgrid(u[:-1]+.2/64,v[:-1]+.2/64,indexing='ij')
    potential=140+10*U+20*V
    psi=diamond_advance(f,g,potential,(.2/32)**2)
    assert discrete_symplectic_balance(psi)['relative_residual']<1e-11


def test_parent_horizon_trace_fixes_second_leg_in_the_declared_sector():
    source=solve_source(endpoint_u=2.);domain=Diamond.inherited(source.parent)
    _,v=domain.axes(64);g,ratio=inherited_second_leg(source,domain,v)
    direct=inherited_reference(source,np.full_like(v,domain.u0),v)[0]
    np.testing.assert_allclose(g,direct,atol=2e-16,rtol=2e-14)
    np.testing.assert_allclose(ratio[0],[1,1],atol=2e-15)
    assert domain.eta_corner+domain.width<math.pi/2


def test_receiver_density_agrees_with_independent_embedding_gradient(result):
    source=solve_source(endpoint_u=2.);domain=Diamond.inherited(source.parent)
    u,v=domain.u0+domain.width,domain.v0+domain.width
    psi,du,dv=inherited_reference(source,np.array(u),np.array(v))
    n=np.array(result['domain']['receiver_angular_axis']);chi=domain.chi_receiver
    X=np.r_[math.sin(chi)*n,math.cos(chi)]
    eta,_,a,_=geometry(source,u,v);w=math.asinh(math.tan(float(eta)));z=source.state(w)
    q1,p1,q2,p2=z[2:6]
    field=C[0]*q1*X[0]+C[1]*q2*X[1]*X[3]
    pi=C[0]*p1*X[0]+C[1]*p2*X[1]*X[3]
    ambient=np.array([C[0]*q1,C[1]*q2*X[3],0.,C[1]*q2*X[1]])
    grad=ambient-X*np.dot(X,ambient)
    rho=(pi*pi+np.dot(grad,grad)/float(a*a)+source.parent.mass_squared*field*field)/2
    record=receiver_record(source,domain,psi,du,dv,n)
    assert record['field_over_epsilon']==pytest.approx(field,abs=2e-16)
    assert record['density_over_epsilon_squared_vacuum']==pytest.approx(rho/source.parent.vacuum,rel=2e-13)


def test_refinement_contracts_field_phase_normal_jet_and_noether_errors(result):
    rows=result['resolutions']
    for key in ('relative_RMS_field_error','horizon_normal_derivative_relative_RMS'):
        for a,b in zip(rows,rows[1:]):assert 3.7<a[key]/b[key]<4.3
    for a,b in zip(rows,rows[1:]):
        assert 3.5<a['noether_flux']['relative_residual']/b['noether_flux']['relative_residual']<4.5
        assert 3.7<abs(a['receiver_wrapped_phase_error']/b['receiver_wrapped_phase_error'])<4.3
    assert rows[-1]['relative_RMS_field_error']<.0014
    assert all(row['discrete_symplectic_flux']['relative_residual']<1e-10 for row in rows)


def test_horizon_energy_is_the_previous_parent_loss_and_composition_never_resets(result):
    assert result['horizon_energy_connection']['relative_volume_flux_residual']<1e-10
    assert result['horizon_energy_connection']['finite_grid_horizon_flux_relative_difference']<.0006
    assert result['composition']['max_field_difference']==0
    assert result['composition']['no_new_source_reset'] is True


def test_horizon_trace_alone_is_not_general_exterior_state_selection(result):
    control=result['one_trace_control']
    assert control['horizon_trace_max_difference']==0
    assert control['interior_relative_RMS_change']>.01
    assert abs(control['receiver_field_change_over_epsilon'])>1e-6


def test_saved_source_bindings():
    root=Path(__file__).resolve().parents[1]
    path=root/'results/summaries/nsc-atlas-characteristic-handoff-v1.json'
    if not path.exists():pytest.skip('record produced after source freeze')
    record=json.loads(path.read_text(encoding='utf-8'))
    for name,digest in record['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
