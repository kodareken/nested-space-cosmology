"""Known ODEs, failure accounting, and one manufactured common update."""
import copy
import math
from fractions import Fraction as Q

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_fixed_cap_stage import quiet_fixed_state,coupled_fixed_stage
from recursive_horizons.nsc_atlas_fixed_cap_time import (
    advance_fixed_ssprk3,fixed_hermite,fixed_time_residual,
)
from recursive_horizons.nsc_atlas_coupled_evolution import CoupledEvolutionError
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from test_nsc_atlas_reference_wall_gauge import GEOMETRY

CONTROLS=ControlBudget(maximum_evaluations=3000000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed time controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def setup(budget):
    basis=FixedCapBasis([0,1],budget=budget)
    state=quiet_fixed_state(basis,exterior_cells=4)
    for i,group in enumerate(state.values(),1):
        for j,key in enumerate(group,1): group[key]=np.full_like(group[key],i+j/8)
    return state,dict(basis=basis,exterior_cells=4,budget=budget)


def scaled(state,value):
    return {g:{k:np.asarray(v)*value for k,v in row.items()} for g,row in state.items()}


def compare(a,b,atol=2e-13):
    for g in a:
        for k in a[g]: np.testing.assert_allclose(a[g][k],b[g][k],rtol=2e-14,atol=atol)


def test_known_cubic_solution_dense_derivative_and_exact_cache_keys(budget):
    weights,kw=setup(budget);left=.25;right=.5;mid=np.nextafter(.375,1.)
    def p(x): return x**3+2*x*x-x+3
    def dp(x): return 3*x*x+4*x-1
    seen=[]
    def rhs(state,*,x,budget):
        seen.append(x);budget.check(200)
        return dict(rates=scaled(weights,dp(x)))
    step=advance_fixed_ssprk3(scaled(weights,p(left)),left,mid,right,stage=rhs,**kw)
    assert seen==[left,right,mid,right] and step['stage_calls']==4 and step['retries']==0
    compare(step['accepted_state'],scaled(weights,p(right)))
    for x in (left,.3125,mid,.4375,right):
        row=fixed_hermite(step,x,**kw)
        compare(row['state'],scaled(weights,p(x)))
        compare(row['derivative'],scaled(weights,dp(x)))
    residual=fixed_time_residual(step,.3125,stage=rhs,**kw)
    compare(residual['residual'],scaled(weights,0.))
    assert not residual['uniform_error_bound']
    with pytest.raises(ValueError,match='accepted interval'):
        fixed_hermite(step,np.nextafter(right,1.),**kw)


def test_linear_solution_and_accepted_endpoint_slope_are_independent_of_predictor(budget):
    initial,kw=setup(budget);saved=copy.deepcopy(initial);lam=-3.;h=1/64
    def rhs(state,*,x,budget):
        budget.check(200);return dict(rates=scaled(state,lam))
    step=advance_fixed_ssprk3(initial,0.,h/2,h,stage=rhs,**kw)
    z=lam*h;p3=1+z+z*z/2+z**3/6
    compare(step['accepted_state'],scaled(initial,p3))
    compare(step['rates'][3],scaled(initial,lam*p3))
    assert abs(step['rates'][3]['wall']['r']-step['rates'][1]['wall']['r'])>1e-3
    # Independent exponential solution and Taylor remainder, not a second mesh.
    bound=math.exp(abs(z))*abs(z)**4/24
    for g in initial:
        for k,v in initial[g].items():
            assert np.max(abs(step['accepted_state'][g][k]-v*math.exp(z)))<=np.max(abs(v))*bound
    t=Q(1,3);x=float(t)*h
    row=fixed_time_residual(step,x,stage=rhs,**kw)
    # Scalar Hermite polynomial computed independently from endpoint basis.
    s=float(t)
    H=(2*s**3-3*s*s+1)+(s**3-2*s*s+s)*z+(-2*s**3+3*s*s)*p3+(s**3-s*s)*z*p3
    Hx=((6*s*s-6*s)+(3*s*s-4*s+1)*z+(-6*s*s+6*s)*p3+(3*s*s-2*s)*z*p3)/h
    compare(row['residual'],scaled(initial,Hx-lam*H),atol=3e-13)
    compare(initial,saved,atol=0.)


def test_callback_failure_and_mutation_do_not_change_accepted_input(budget):
    initial,kw=setup(budget);saved=copy.deepcopy(initial);calls=[]
    def failed(state,*,x,budget):
        calls.append(x)
        state['wall']['r'][...]=999.
        if len(calls)==3: raise RuntimeError('declared stop')
        return dict(rates=scaled(saved,0.))
    with pytest.raises(CoupledEvolutionError) as caught:
        advance_fixed_ssprk3(initial,0.,.005,.01,stage=failed,**kw)
    d=caught.value.diagnostic
    assert d['stage_calls']==3 and d['completed_stages']==2 and d['phase']=='midpoint_stage'
    compare(d['accepted_state'],saved,atol=0.);compare(initial,saved,atol=0.)
    with pytest.raises(CoupledEvolutionError) as caught:
        advance_fixed_ssprk3(initial,0.,.004,.01,stage=failed,**kw)
    assert caught.value.diagnostic['stage_calls']==0
    def malformed(state,*,x,budget):
        state['cap']['chi']=np.zeros(2);return dict(rates=state)
    with pytest.raises(CoupledEvolutionError) as caught:
        advance_fixed_ssprk3(initial,0.,.005,.01,stage=malformed,**kw)
    assert caught.value.diagnostic['stage_calls']==1


def test_one_manufactured_nonzero_common_step_and_midpoint_residual(budget):
    # This smooth cubic source tests integration, not the retained preparation.
    left=.4;h=1/1024;right=left+h;mid=left+h/2
    Z=Q(float(np.arctanh(.45)))
    basis=FixedCapBasis([0,2*Z],budget=budget)
    initial=quiet_fixed_state(basis,exterior_cells=4)
    gauge=ReferenceWallGauge(GEOMETRY,budget=budget)
    def stage(state,*,x,budget):
        t=x-left
        source=dict(r=x,log_C=0.,c=t**3,r_x=1.,log_C_x=0.,c_x=3*t*t,r_xx=0.)
        reference=gauge.sample(x)
        row=coupled_fixed_stage(state,basis=basis,x=x,eta=np.linspace(0,1,5),
            source=source,gauge=reference,alpha=.002,
            tube=dict(r=[.2,.8],L=[-2.,2.],c=[-4.,4.]),budget=budget)
        row['test_W_x']=reference['W_x']
        return row
    kw=dict(basis=basis,exterior_cells=4,budget=budget)
    step=advance_fixed_ssprk3(initial,left,mid,right,stage=stage,**kw)
    row=fixed_time_residual(step,mid,stage=stage,**kw)
    assert np.max(abs(step['accepted_state']['cap']['canonical_density']))>0.
    assert np.max(abs(step['accepted_state']['departure']['c']))>0.
    assert np.max(abs(row['residual']['cap']['canonical_density']))>0.
    assert len(step['stages'])==4 and row['stage_calls']==1
    for actual in [*step['stages'],row['stage']]:
        assert actual['physical_wall']['r']==actual['connection']['cap']['radius']
        for key in ('r','L','c'):
            assert actual['wall_kinematic_defect'][key]==pytest.approx(
                -actual['test_W_x']*actual['wall_incoming_jet_defect'][key],abs=3e-13)
        assert actual['source_constraint'] is not None
        connection=actual['connection'];defect=connection['residual']
        R=actual['physical_wall']['r'];rootD=math.sqrt(float(connection['hybrid']['D']))
        assert defect['conformal_wall_flux']==pytest.approx(
            R*rootD*defect['scalar_tangent'],abs=2e-20)
        for suffix in ('actual','matched'):
            assert defect['physical_temporal_mass_'+suffix]==pytest.approx(
                .002*defect['temporal_mass_'+suffix],rel=2e-15,abs=1e-30)
    # The derived wall scalar is absent from the stored state. Check its
    # moving-trace derivative separately against an independent difference.
    def scalar(path,reference):
        R=float(reference['r']+path['state']['wall']['r']);Z=float(np.arctanh(R))
        shape=basis.values(0,Z/float(2*Q(float(np.arctanh(.45)))),budget=budget)
        coeff=path['state']['cap']['chi'];f=float(shape['value']@coeff)/R
        k=math.sqrt(.002/6);u=math.sinh(k*f)/k
        Rx=float(reference['r_x_total']+path['derivative']['wall']['r'])
        chix=float(shape['value']@path['derivative']['cap']['chi']
                   +(shape['first']@coeff)*Rx/(1-R*R))
        return R*u,math.cosh(k*f)*chix+(u-f*math.cosh(k*f))*Rx
    path=fixed_hermite(step,mid,**kw)
    _,analytic=scalar(path,gauge.sample(mid))
    delta=h/4096
    minus=scalar(fixed_hermite(step,mid-delta,**kw),gauge.sample(mid-delta))[0]
    plus=scalar(fixed_hermite(step,mid+delta,**kw),gauge.sample(mid+delta))[0]
    assert analytic==pytest.approx((plus-minus)/(2*delta),abs=2e-13)
    # At the same state, inserting the stage rates recovers its actual c_x.
    stage_path=dict(state=path['state'],derivative=row['stage']['rates'])
    assert scalar(stage_path,gauge.sample(mid))[1]==pytest.approx(
        row['stage']['connection']['cap_c_x'],abs=2e-20)
    for group in row['residual'].values():
        for v in group.values(): assert np.all(np.isfinite(v))
