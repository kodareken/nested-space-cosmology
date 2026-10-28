"""Independent action, matching and evolving-source controls."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_inherited_matter import exact_initial_geometry
from recursive_horizons.nsc_atlas_parent_scale_selector import (
    VolumeParent, calculation, evolve, matched_coefficients)


def test_reduction_measure_inverse_metrics_and_wrapped_brane():
    d, s, V0 = sp.symbols('d s V0', positive=True)
    # Determinant weight, external and internal inverse-metric contractions.
    measure = sp.exp(-d*s)
    top_form = measure*sp.exp(4*d*s)*sp.exp(-2*d*s)/V0
    scalar_kinetic = measure*sp.exp(d*s)
    wall_measure = sp.exp(-3*d*s/2)*sp.exp(d*s)
    assert sp.simplify(top_form-sp.exp(d*s)/V0) == 0
    assert sp.simplify(scalar_kinetic-1) == 0
    assert sp.simplify(wall_measure-sp.exp(-d*s/2)) == 0
    # Jordan reduction R -> R -2d box(s)-d(d+1)(ds)^2;
    # Einstein Weyl term adds -3d^2/2; integrated box adds +2d^2.
    coefficient = sp.simplify(-d*(d+1)+2*d*d-3*d*d/2)
    assert sp.simplify(coefficient+d*(d+2)/2) == 0


def test_minimal_positive_vacuum_and_general_stiffness_identity():
    d, s, rho = sp.symbols('d s rho', positive=True)
    A, B = sp.symbols('A B')
    W = sp.Function('W')(s)
    V = A*sp.exp(-d*s)-B*sp.exp(-(d+2)*s)+W
    match = sp.solve([V.subs(s, 0)-rho, sp.diff(V, s).subs(s, 0)], (A, B))
    ss = sp.diff(V, s, 2).subs(s, 0).subs(match)
    target = -d*(d+2)*rho+(sp.diff(W,s,2)+2*(d+1)*sp.diff(W,s)+d*(d+2)*W).subs(s,0)
    assert sp.simplify(ss-target) == 0
    K = d*(d+2)/2
    assert sp.simplify((-d*(d+2)*rho/K)/(rho/3)) == -6


def test_ancestor_state_changes_the_stability_requirement():
    d, s, E = sp.symbols('d s E', positive=True)
    T, rho = sp.symbols('T rho', positive=True)
    B = d*(rho+T+E/2)/2
    A = B+rho-T
    V = A*sp.exp(-d*s)-B*sp.exp(-(d+2)*s)+T*sp.exp(-2*d*s)+E*sp.exp(-d*s/2)
    assert sp.simplify(V.subs(s,0)-rho-E) == 0
    assert sp.simplify(sp.diff(V,s).subs(s,0)) == 0
    ss = sp.diff(V,s,2).subs(s,0)
    frozen = {rho:sp.Rational(9,32), E:sp.Rational(3,8)}
    assert sp.simplify(32*ss.subs(frozen).subs(T,sp.Rational(1,2))/d) == -2*d-62
    assert sp.simplify(32*ss.subs(frozen).subs(T,sp.Rational(25,32))/d) == 7*d-80
    assert matched_coefficients(11,2)['initial_stiffness'] < 0
    c = matched_coefficients(12,2)
    assert (c['A'],c['B'],c['T'],c['initial_stiffness'],c['kinetic']) == (7,F(15,2),F(25,32),F(3,2),84)
    assert c['required_UV_bare'] == F(215,32)


def test_exact_energy_exchange_and_generalized_corrected_action():
    N,s,v,chi,pi,H,w,K,d,V,Vs = sp.symbols('N s v chi pi H w K d V Vs', real=True)
    E = (pi*pi+w*w*chi*chi)/2
    dots = {N:H,s:v,v:-3*H*v+(-Vs+d*w*w*chi*chi/2)/K,
            chi:pi,pi:-3*H*pi-w*w*chi,w:-d*v*w/2,V:Vs*v}
    dt = lambda f: sp.expand(sum(sp.diff(f,x)*dx for x,dx in dots.items()))
    total = K*v*v/2+V+E
    assert sp.simplify(dt(total)+3*H*(K*v*v+pi*pi)) == 0
    D = sp.exp(3*N)*E/w
    S = sp.exp(3*N)*chi*pi
    c = (3*H-d*v/2)/(2*w)
    assert sp.simplify(dt(D)+c*dt(S)) == 0
    # I=D+cS then I_dot=c_dot*S follows without a fixed frequency assumption.
    assert sp.simplify(D**2-S**2-sp.exp(6*N)*(pi*pi-w*w*chi*chi)**2/(4*w*w)) == 0


def test_exact_full_source_initial_density_and_average_force():
    m = VolumeParent()
    initial = np.array([0.,0.,0.,math.sqrt(.75)/400,0.])
    dy = m.scalar_rhs(0.,initial,400)
    assert dy[0]**2 == pytest.approx(float(exact_initial_geometry()['H02']),rel=2e-14)
    # Initial phase is coherent: exact force differs from phase average.
    assert dy[2] > 0
    assert m.averaged_rhs(0.,[0.,0.,0.])[2] == 0
    c = matched_coefficients(12,2)
    assert m.potential(0) == float(c['vacuum'])


def test_finite_evolution_extrema_clock_and_averaging_controls():
    r = calculation()
    v = next(x for x in r['vacuum_extrema'] if x['kind']=='minimum')
    assert v['vacuum'] > 0 and v['stiffness'] > 0
    assert abs(VolumeParent().potential(v['s'],1)) < 1e-10
    assert r['late_vacuum_Lref0_squared'] == pytest.approx(9.7022370093,rel=1e-9)
    assert r['volume_field_clock_difference'] < 1e-9
    assert abs(r['averaged_evolution']['energy_ledger_residual']) < 1e-9
    controls = r['resolved_scalar_controls']
    assert controls[1]['absolute_volume_field_difference_from_average'] < .4*controls[0]['absolute_volume_field_difference_from_average']
    for row in controls:
        assert abs(row['corrected_action_ledger_relative']) < 1e-8
        lo, hi = row['action_ratio_enclosure']
        assert lo < row['action_ratio'] < hi
    assert not r['SI_scale_selected'] and not r['Lambda_obs_derived']
    assert r['required_wall_volume_normal_derivative_jump'] == pytest.approx(-1/70)
    assert r['averaged_basin_entry']['barrier_energy_margin'] > .13
    assert r['averaged_basin_entry']['expanding_branch_H_squared_floor'] > .07
    assert json.loads(json.dumps(r))['averaged_basin_entry']['left_of_barrier']


def test_finite_candidate_rejects_inappropriate_cycle_or_clock_budget():
    with pytest.raises(ValueError): evolve(endpoint=20)
    with pytest.raises(ValueError): evolve(endpoint=4,omega=1e9)
    with pytest.raises(ValueError): matched_coefficients(0,1)


def test_finite_solution_tolerance_controls_its_claimed_digits():
    loose,_ = evolve()
    tight,_ = evolve(rtol=2e-12,atol=2e-14)
    assert abs(loose['s']-tight['s']) < 2e-10
    assert abs(loose['vacuum_density']-tight['vacuum_density']) < 2e-10
