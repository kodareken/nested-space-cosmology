"""Exact family closure, operator controls, and producer-free directed replay."""
import importlib.util
import json
from pathlib import Path
import sympy as sp
import recursive_horizons.nsc_incoming_surface_coefficients as owner
ROOT=Path(__file__).resolve().parents[1]

def read(name):return json.loads((ROOT/f'results/development/nsc-incoming-{name}.json').read_text())

def test_full_two_dimensional_local_C_and_independent_reference_identity():
    proof=owner.polynomial_proof(read('shift-coefficient'))
    assert set(proof['identities'].values())=={'0'}
    local=owner.local_surface_polynomials();a,r,Ha,Hr,A2,R2,A3,R3,w,A,CW,CE,Cg,Cq,hq,q=local['symbols']
    assert sp.factor(sum(local['C_local_channels'].values())-64*sp.pi*CW/(3*a)+(hq+sp.log(r))/(15*sp.pi*a))==0
    for name in ('light/LLL_geometry_bar','compact/einstein_bulk','light/WZ_Euler','light/WZ_boxR'):
        assert local['C_local_channels'][name]==0

def test_homogeneous_polynomial_preserves_fixed_raw_second_third_derivatives():
    ref=owner.homogeneous_reference_energy();m,L,p,W,a,r,D,Ha,Hr,A2,R2,A3,R3,w=ref['symbols']
    energy=sum(ref['integrated_energy_orders'].values())
    polynomial=sum(value*w**j for j,value in ref['D_reference_coefficients'].items())
    assert sp.factor(energy.subs(Hr,Hr+w/r)-energy-polynomial)==0
    assert sp.factor(polynomial.subs(L,0))==0
    assert ref['odd_order_full_line_residual']==ref['reused_Cu_equation_residual']=='0'
    local=owner.local_surface_polynomials();a,r,Ha,Hr,A2,R2,A3,R3,w,*_=local['symbols']
    for name,energy in local['homogeneous_lapse_channels'].items():
        polynomial=sum(value*w**j for j,value in local['D_local_coefficients'][name].items())
        assert sp.factor(energy.subs(Hr,Hr+w/r)-energy-polynomial)==0

def test_new_w_point_energy_and_original_density_bridge():
    result=owner.new_algebra_controls(read('local-constraints')['locked_inputs'])
    assert result['new_w_point_energy_maximum']<3e-13
    assert result['local_density_maximum']<3e-13
    assert result['new_reference_point_count']==4
    assert result['new_reference_momentum_integrations']==0

def test_replay_never_enters_old_or_new_symbolic_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('symbolic/probe producer entered replay')
    for name in ('homogeneous_reference_energy','local_density_expressions','local_surface_polynomials','polynomial_proof','new_algebra_controls'):
        monkeypatch.setattr(owner,name,forbidden)
    import recursive_horizons.nsc_incoming_lapse_coefficient as cu
    import recursive_horizons.nsc_incoming_shift_coefficient as cv
    import recursive_horizons.nsc_spatial_reference_symbol as ref
    monkeypatch.setattr(cu,'closed_lapse_coefficient',forbidden);monkeypatch.setattr(cv,'closed_shift_coefficient',forbidden)
    monkeypatch.setattr(ref,'reference_projector',forbidden)
    spec=importlib.util.spec_from_file_location('surface_coefficients_replay',ROOT/'scripts/derive_nsc_incoming_surface_coefficients.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    prior=read('surface-coefficients')
    assert script.make_record(ROOT/prior['payload']['path'])==prior
    assert prior['coefficients']['polynomial_intervals']['F']['0']==read('shift-coefficient')['coefficient']['total_coefficient_interval']
    assert prior['coefficients']['C_interval'][0]>0
    assert prior['stored_controls']['maximum']<3e-11
    assert not prior['scope']['w_profile_integrated']
