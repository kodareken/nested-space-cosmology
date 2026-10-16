"""Bounded common-value, frame, vacuum and independent energy controls."""
from fractions import Fraction as Q

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_fixed_hybrid import fixed_cap_hybrid_reference
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice
from test_nsc_atlas_reference_wall_gauge import GEOMETRY
from test_nsc_atlas_cap_fixed_reference import exact_operators

CONTROLS=ControlBudget(maximum_evaluations=1500000,cpu_seconds=45.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed hybrid controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def fields(budget):
    basis=FixedCapBasis([0,1],budget=budget)
    return basis,np.array([.001,.002,.003]),np.array([.003,-.002,.001]),Q(1,2)


def test_callback_matches_direct_port_and_independent_zero_alpha_energy(budget):
    basis,c,d,Z=fields(budget);seen=[]
    inputs=dict(T_x=1.2,radius_rate=.001,incoming_w_plus=.03)
    direct=fixed_cap_reference(basis,c,d,Z,alpha=0.,budget=budget,**inputs)
    def port(packet):
        with pytest.raises(TypeError): packet['radius']=0.
        seen.append(dict(packet))
        return inputs
    coupled=fixed_cap_reference(basis,c,d,Z,alpha=0.,budget=budget,wall_coupler=port)
    assert len(seen)==1
    for key in ('mass','lapse','mass_rate','physical_weak_flux','weak_balance_defect'):
        assert coupled[key]==direct[key]
    for key in direct['rates']:
        np.testing.assert_array_equal(coupled['rates'][key],direct['rates'][key])
    M,K,e=exact_operators(basis,Z,budget)
    R=np.tanh(float(Z));A=1/np.cosh(float(Z))**2
    stiffness=K-(A/R)*np.outer(e,e)
    assert direct['mass']==pytest.approx((c@stiffness@c+d@M@d)/4,abs=2e-12)
    np.testing.assert_allclose(direct['rates']['chi'],1.2*d,rtol=2e-8,atol=2e-10)
    expected=1.2*np.linalg.solve(M,-stiffness@c+e*(A*R*.03-e@d))
    np.testing.assert_allclose(direct['rates']['canonical_density'],expected,rtol=3e-8,atol=3e-9)


def test_quiet_material_stage_recovers_existing_vacuum_jets_and_clocks(budget):
    gauge=ReferenceWallGauge(GEOMETRY,budget=budget).sample(.4)
    Z=Q(float(np.arctanh(.45)));basis=FixedCapBasis([0,2*Z],budget=budget)
    zero=np.zeros(basis.dimension);seen=[]
    vac=vacuum_slice(gauge['W'],.4,budget=budget)
    def exterior(common):
        seen.append(dict(common))
        assert common['c']==0.
        return dict(r_V=float(vac['r_x']),c_V=0.,log_C_V=float(vac['log_C_x']))
    result=fixed_cap_hybrid_reference(basis,zero,zero,Z,alpha=.002,exterior_wall=exterior,
        C_wall=np.exp(gauge['L']),eta_plus=gauge['eta_plus'],eta_plus_x=gauge['eta_plus_x'],budget=budget)
    assert len(seen)==1
    for key,name in (('r','r_x_total'),('L','L_x_total'),('T','T_x'),('s','s_x')):
        assert result['wall_rates'][key]==pytest.approx(gauge[name],abs=4e-13)
    assert all(value==0. for value in result['residual'].values())
    assert result['physical_mass_rate']==pytest.approx(1.5*result['common']['r']**2*result['wall_rates']['r'])
    assert all(np.all(value==0.) for value in result['cap']['rates'].values())


def test_nonzero_stage_preserves_common_value_and_exposes_actual_trace_residuals(budget):
    basis,c,d,Z=fields(budget);untouched=(c.copy(),d.copy());seen=[]
    def exterior(common):
        with pytest.raises(TypeError): common['c']=1.
        seen.append(dict(common))
        return dict(r_V=.11,c_V=.006,log_C_V=.003)
    start=budget.evaluations
    value=fixed_cap_hybrid_reference(basis,c,d,Z,alpha=.02,exterior_wall=exterior,
        C_wall=1.,eta_plus=0.,eta_plus_x=.01,budget=budget)
    assert value['charged_evaluations']==budget.evaluations-start
    assert len(seen)==1
    cap=value['cap'];wall=cap['wall'];hybrid=value['hybrid'];res=value['residual']
    R=cap['radius'];tx,rx=value['wall_rates']['T'],value['wall_rates']['r']
    e=basis.values(0,.5,budget=budget)
    f=(e['value']@c)/R;scale=np.sqrt(.02/6)
    assert seen[0]['c']==pytest.approx(R*np.sinh(scale*f)/scale,abs=3e-16)
    assert cap['T_x']==tx and cap['radius_rate']==rx
    assert (wall['p_weak']+wall['h_weak'])/R==pytest.approx(hybrid['w_plus'])
    assert (wall['p_trace']-wall['h_trace'])/R==pytest.approx(hybrid['Pi']-hybrid['Phi'])
    assert res['scalar_tangent']==pytest.approx(res['scalar_tangent_factored'],abs=3e-15)
    assert res['conformal_wall_flux']==pytest.approx(res['conformal_wall_flux_factored'],abs=3e-15)
    assert abs(res['scalar_tangent'])>1e-6
    for pair in ('actual','matched'):
        assert res['temporal_mass_'+pair]-res['weak_flux_'+pair]==pytest.approx(cap['weak_balance_defect']/tx,abs=1e-17)
        assert res['physical_temporal_mass_'+pair]==.02*res['temporal_mass_'+pair]
    vacuum_removed=(value['physical_mass_rate']-rx*value['physical_mass_density'])/tx
    matched=wall['A']**2*np.exp(cap['lapse'])*R**2*float(hybrid['Pi'])*float(hybrid['Phi'])/2
    assert vacuum_removed-.02*matched==pytest.approx(res['physical_temporal_mass_matched'],abs=2e-16)
    # Differentiate the reconstructed field along the actual coefficient
    # and moving-radius direction, without another metric/gradient solve.
    step=1e-7;trace=[];invariant=[]
    for sign in (-1,1):
        radius=R+sign*step*rx;z=np.arctanh(radius)
        row=basis.values(0,z,budget=budget)
        fnew=(row['value']@(c+sign*step*cap['rates']['chi']))/radius
        trace.append(fnew);invariant.append(radius*np.sinh(scale*fnew)/scale)
    assert (trace[1]-trace[0])/(2*step)==pytest.approx(value['cap_f_x'],abs=2e-9)
    assert (invariant[1]-invariant[0])/(2*step)==pytest.approx(value['cap_c_x'],abs=2e-9)
    np.testing.assert_array_equal(c,untouched[0]);np.testing.assert_array_equal(d,untouched[1])
    budget.check(4000)


def test_invalid_port_stops_without_a_second_or_prescribed_wall_input(budget):
    basis,c,d,Z=fields(budget)
    for port in (lambda _:None,False):
        with pytest.raises(ValueError,match='mutually exclusive'):
            fixed_cap_reference(basis,c,d,Z,alpha=0.,budget=budget,wall_coupler=port,
                                T_x=1.,radius_rate=0.,incoming_w_plus=0.)
    zero=np.zeros(basis.dimension)
    with pytest.raises(ValueError,match='exactly the three'):
        fixed_cap_reference(basis,zero,zero,Z,alpha=0.,budget=budget,
                            wall_coupler=lambda _:dict(T_x=1.))
    with pytest.raises(ValueError,match='three outgoing'):
        fixed_cap_hybrid_reference(basis,zero,zero,Z,alpha=0.,budget=budget,
            exterior_wall=lambda _:dict(r_V=.1,c_V=0.),C_wall=1.,eta_plus=0.,eta_plus_x=0.)
