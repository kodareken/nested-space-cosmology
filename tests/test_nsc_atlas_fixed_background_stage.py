"""Physical-clock history hookup and one common step; no prepared source."""
import copy
from fractions import Fraction as Q

import numpy as np
from numpy.polynomial import Polynomial
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_transfer import _transfer_coefficients
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_fixed_cap_stage import coupled_fixed_stage
from recursive_horizons.nsc_atlas_fixed_cap_time import advance_fixed_ssprk3,fixed_time_residual
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from test_nsc_atlas_fixed_background_reference import history
from test_nsc_atlas_fixed_cap_stage import case
from test_nsc_atlas_reference_wall_gauge import GEOMETRY

CONTROLS=ControlBudget(maximum_evaluations=8000000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('common history controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def manufactured_case(budget):
    state,kw=case(budget,active=True)
    state['wall']['T']=np.asarray(1/32)
    bg=history(Polynomial([0,0,0,1/256]),knots=(-4.,-2.,0.,2.))
    return state,kw,bg


def test_total_polynomial_and_independent_common_mass_direction(budget):
    state,kw,bg=manufactured_case(budget)
    affine=coupled_fixed_stage(state,background=bg,**kw)
    T=Q(float(kw['gauge']['T']+state['wall']['T']));a=Q(1,256)
    assert affine['connection']['cap']['background']['time']==float(T)
    def coeff(power):
        # Transfer consumes powers of the local cell coordinate, z/h.
        h=kw['basis'].exact_knots[1]
        local=tuple(value*h**p for p,value in enumerate(power))
        return np.array(_transfer_coefficients(kw['basis'],[local],budget),dtype=float)
    cb=coeff([0,6*a*T*T,0,2*a,0,0])
    kb=coeff([0,12*a*T,0,0,0,0])
    bb=coeff([0,12*a,0,0,0,0])
    total_state=copy.deepcopy(state)
    total_state['cap']['chi']+=cb;total_state['cap']['canonical_density']+=kb
    total=coupled_fixed_stage(total_state,**kw)
    for key in ('r','L','T','s','c'):
        assert affine['physical_wall'][key]==pytest.approx(total['physical_wall'][key],abs=2e-12)
    for group in ('departure','wall'):
        for key in affine['rates'][group]:
            np.testing.assert_allclose(affine['rates'][group][key],total['rates'][group][key],rtol=0,atol=2e-10)
    tx=affine['wall_rates']['T'];rx=affine['wall_rates']['r']
    Z=Q(float(np.arctanh(affine['physical_wall']['r'])))
    gram=kw['basis'].gram(Z,budget=budget)
    for key,rate in (('chi',kb),('canonical_density',bb)):
        gap=affine['rates']['cap'][key]+tx*rate-total['rates']['cap'][key]
        # A clipped basis can magnify individual extension coefficients.
        # Compare the actual cap field and wall trace in their own norms.
        assert np.sqrt(gram.squared_norm(gap,budget=budget))<2e-9
        assert abs(gram.wall['value']@gap[gram.wall['indices']])<2e-9
    cap=affine['connection']['cap'];R=cap['radius']
    assert cap['mass']==pytest.approx(total['connection']['cap']['mass'],abs=2e-11)
    assert abs(tx*cap['gradient']['T'])>2e-6
    h=1e-5;mass=[]
    for sign in (-1,1):
        delta=sign*h
        c=state['cap']['chi']+delta*affine['rates']['cap']['chi']
        d=state['cap']['canonical_density']+delta*affine['rates']['cap']['canonical_density']
        shifted=fixed_cap_reference(kw['basis'],c,d,Q(float(np.arctanh(R+delta*rx))),
            alpha=kw['alpha'],T_x=tx,radius_rate=rx,incoming_w_plus=0.,
            background=bg,background_time=float(T)+delta*tx,budget=budget)
        mass.append(shifted['mass'])
    difference=(mass[1]-mass[0])/(2*h)
    assert difference==pytest.approx(cap['mass_rate'],abs=2e-8)
    physical=1.5*R*R*rx+kw['alpha']*difference
    assert physical==pytest.approx(affine['connection']['physical_mass_rate'],abs=4e-11)
    print('common mass direction:',difference,cap['mass_rate'],'explicit clock term:',tx*cap['gradient']['T'])


def test_one_common_history_step_keeps_stage_clock_and_mass_remainder(budget):
    initial,kw,bg=manufactured_case(budget);saved=copy.deepcopy(initial)
    left=kw.pop('x');kw.pop('gauge');kw.pop('budget')
    gauge=ReferenceWallGauge(GEOMETRY,budget=budget)
    h=1/16384;right=left+h;mid=left+h/2;seen=[]
    def stage(state,*,x,budget):
        reference=gauge.sample(x)
        row=coupled_fixed_stage(state,x=x,gauge=reference,budget=budget,background=bg,**kw)
        T=float(reference['T']+state['wall']['T']);cap=row['connection']['cap']
        assert cap['background']['time']==T
        assert abs(T-reference['T'])>.02
        assert cap['T_x']==row['wall_rates']['T']
        expected=(cap['T_x']*cap['gradient']['T']+
                  cap['gradient']['chi']@row['rates']['cap']['chi']+
                  cap['gradient']['canonical_density']@row['rates']['cap']['canonical_density']+
                  cap['moving_density'])
        assert cap['mass_rate']==pytest.approx(expected,abs=2e-15)
        assert abs(cap['remainder_identity_defect'])<2e-12
        seen.append((x,T))
        return row
    opts=dict(basis=kw['basis'],exterior_cells=4,stage=stage,budget=budget)
    step=advance_fixed_ssprk3(initial,left,mid,right,**opts)
    residual=fixed_time_residual(step,mid,**opts)
    assert len(seen)==5 and step['retries']==0
    assert seen[1][0]==seen[3][0] and seen[1][1]!=seen[3][1]
    assert not residual['uniform_error_bound']
    for group in residual['residual'].values():
        for value in group.values(): assert np.all(np.isfinite(value))
    for group in initial:
        for key in initial[group]: np.testing.assert_array_equal(initial[group][key],saved[group][key])
    print('common history midpoint residuals:',{g:max(float(np.max(abs(v))) for v in row.values())
          for g,row in residual['residual'].items()})


def test_constant_history_clock_origin_preserves_the_common_physical_state(budget):
    state,kw,bg=manufactured_case(budget);data=bg.inputs;shift=.25
    shifted=RetainedCapBackground(data['labels']+shift,data['values'],data['first'],data['second'])
    original=coupled_fixed_stage(state,background=bg,**kw)
    moved=coupled_fixed_stage(state,background=shifted,background_clock_offset=shift,**kw)
    assert moved['connection']['cap']['background']['time']==original['connection']['cap']['background']['time']+shift
    for key in original['physical_wall']:
        assert moved['physical_wall'][key]==pytest.approx(original['physical_wall'][key],rel=0.,abs=2e-13)
    for group in original['rates']:
        for key in original['rates'][group]:
            np.testing.assert_allclose(moved['rates'][group][key],original['rates'][group][key],rtol=2e-11,atol=2e-10)
    assert moved['connection']['cap']['mass_rate']==pytest.approx(original['connection']['cap']['mass_rate'],rel=0.,abs=2e-12)
    with pytest.raises(ValueError,match='requires a retained background'):
        coupled_fixed_stage(state,background_clock_offset=shift,**kw)
    with pytest.raises(ValueError,match='background clock offset'):
        coupled_fixed_stage(state,background=bg,background_clock_offset=float('nan'),**kw)
