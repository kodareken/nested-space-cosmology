"""One common exterior/cap state, with independent boundary discrepancies."""
import copy
from fractions import Fraction as Q

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_fixed_cap_stage import (
    quiet_fixed_state,fixed_state,exterior_departure_stage,coupled_fixed_stage,
)
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_gradient_evolution import coupled_stage as old_stage
from test_nsc_atlas_reference_wall_gauge import GEOMETRY
from test_nsc_atlas_gradient_evolution import stage_options,initial_state

CONTROLS=ControlBudget(maximum_evaluations=2000000,cpu_seconds=60.)
TUBE=dict(r=[.2,.8],L=[-2.,2.],c=[-4.,4.])


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed common stage controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def case(budget,active=False):
    Z=Q(float(np.arctanh(.45)));basis=FixedCapBasis([0,2*Z],budget=budget)
    state=quiet_fixed_state(basis,exterior_cells=4)
    gauge=ReferenceWallGauge(GEOMETRY,budget=budget).sample(.4)
    source=dict(r=.4,log_C=0.,c=0.,r_x=1.,log_C_x=0.,c_x=.03 if active else 0.,r_xx=0.)
    if active:
        state['cap']['chi']=np.array([.001,.002,.003])
        state['cap']['canonical_density']=np.array([.003,-.002,.001])
        state['departure']['r']=np.array([.0001,.0002,.0001])
        state['departure']['L']=np.array([.0002,.0003,.0001])
        state['departure']['c']=np.array([.0004,.0008,.0007])
    return state,dict(basis=basis,x=.4,eta=np.linspace(0,1,5),source=source,gauge=gauge,
                      alpha=.002,tube=TUBE,budget=budget)


def test_exterior_orchestration_reproduces_frozen_stage_fields(budget):
    kw=stage_options(budget,active=True,alpha=.01,n=4);state=initial_state(4,active=True)
    old=old_stage(state,**kw)
    row=exterior_departure_stage(state['departure'],old['physical_wall'],x=kw['x'],eta=kw['eta'],
        source=kw['source'],W=kw['gauge']['W'],W_x=kw['gauge']['W_x'],alpha=kw['alpha'],
        tube=kw['tube'],budget=budget)
    for new,prior in ((row['exterior'],old['exterior']),
                      (row['longitudinal'],old['exterior_longitudinal']),
                      (row['gradients'],old['departure_gradients']),
                      (row['rates'],old['rates']['departure'])):
        for key in new: np.testing.assert_array_equal(new[key],prior[key])
    for key in ('r_x','c_x','log_C_x','r_Ubar','c_Ubar','log_C_Ubar'):
        np.testing.assert_array_equal(row['prefix'][key],old['exterior_prefix'][key])


def test_quiet_common_stage_recovers_all_vacuum_rates(budget):
    state,kw=case(budget);before=budget.evaluations
    value=coupled_fixed_stage(state,**kw)
    assert value['charged_evaluations']==budget.evaluations-before
    for group in value['rates'].values():
        for row in group.values(): np.testing.assert_allclose(row,0.,rtol=0.,atol=4e-13)
    assert value['source_constraint']==0.
    assert max(abs(v) for v in value['wall_kinematic_defect'].values())<4e-13
    assert max(abs(v) for v in value['wall_incoming_jet_defect'].values())<4e-13
    assert abs(value['coordinate_radius_difference'])<2e-16
    assert 'courant_rate' not in value
    assert value['exterior_characteristic_courant_rate']>0.


def test_active_common_stage_uses_one_gradient_and_separates_the_two_tangents(budget):
    state,kw=case(budget,True);saved=copy.deepcopy(state);seen=[]
    def derivative(eta,fields,*,budget):
        assert not eta.flags.writeable and not fields.flags.writeable
        seen.append(fields.copy());budget.check(fields.size)
        return np.stack((.001*(1+eta),-.002*eta,.003*(1-eta)))
    value=coupled_fixed_stage(state,departure_gradient=derivative,**kw)
    assert len(seen)==1
    W,Wx=kw['gauge']['W'],kw['gauge']['W_x'];eta=kw['eta'];ext=value['exterior_prefix']
    for key,jet,delta,source_rate in (
            ('r','r_Ubar','delta_r_x',0.),('L','log_C_Ubar','delta_log_C_x',0.),
            ('c','c_Ubar','delta_c_x',.03)):
        G=value['departure_gradients'][key]
        base=0. if key=='c' else ext['vacuum'][jet]
        np.testing.assert_allclose(value['exterior_first_jets'][key],base+G/W,rtol=2e-15,atol=2e-15)
        np.testing.assert_allclose(value['rates']['departure'][key],
            (ext[delta]-source_rate+eta*(Wx/W)*G)[1:-1],rtol=2e-15,atol=2e-15)
        assert value['exterior'][key][0]==kw['source']['log_C' if key=='L' else key]
        assert value['exterior'][key][-1]==value['physical_wall'][key]
        assert value['wall_kinematic_defect'][key]==pytest.approx(
            -Wx*value['wall_incoming_jet_defect'][key],abs=2e-14)
    difference=value['shared_wall_update_defect']['c']-value['wall_kinematic_defect']['c']
    assert difference==pytest.approx(value['connection']['residual']['conformal_wall_flux'],abs=2e-15)
    assert abs(difference)>1e-6
    assert value['shared_wall_rates']['c']==value['connection']['cap_c_x']
    assert 'c' not in value['rates']['wall']
    assert value['source_constraint']==pytest.approx(.002*.03**2/(2*.4),abs=2e-18)
    for group in state:
        for key in state[group]: np.testing.assert_array_equal(state[group][key],saved[group][key])


def test_missing_source_curvature_is_unavailable_and_bad_domains_stop(budget):
    state,kw=case(budget)
    source=dict(kw['source']);source.pop('r_xx')
    row=coupled_fixed_stage(state,**dict(kw,source=source))
    assert row['source_constraint'] is None
    bad=copy.deepcopy(state);bad['cap']={'Ip':np.zeros(3),'Iq':np.zeros(3)}
    with pytest.raises(ValueError,match='fixed cap fields'): fixed_state(bad,4,kw['basis'])
    bad=copy.deepcopy(state);bad['wall']['r']=1.
    with pytest.raises(ValueError,match='0<r<1'): coupled_fixed_stage(bad,**kw)
    bad_gauge=dict(kw['gauge'],eta_plus=1.)
    with pytest.raises(ValueError,match='compatible boost'):
        coupled_fixed_stage(state,**dict(kw,gauge=bad_gauge))
    with pytest.raises(ValueError,match='departure gradient'):
        coupled_fixed_stage(state,departure_gradient=None,**kw)
    narrow=FixedCapBasis([0,Q(1,4)],budget=budget)
    with pytest.raises(ValueError,match='wall cut'):
        coupled_fixed_stage(state,**dict(kw,basis=narrow))
