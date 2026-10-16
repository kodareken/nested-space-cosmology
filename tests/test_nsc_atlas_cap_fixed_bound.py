"""Independent polynomial hull identities and constrained-metric limits."""
from fractions import Fraction as Q
from math import comb

import numpy as np
import pytest
from flint import arb,ctx

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_bound import _ball,_bernstein_controls,fixed_metric_bound

CONTROLS=ControlBudget(maximum_evaluations=500_000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed metric controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def linear_coefficients(basis,slope):
    result=np.zeros(basis.dimension)
    result[0]=float(slope*(basis.exact_knots[1]-basis.exact_knots[0]))
    for node in range(1,basis.cells+1):
        result[4*node-3]=float(slope*basis.exact_knots[node])
        result[4*node-2]=float(slope*basis.exact_scales[node])
    return result


def test_exact_bernstein_expansion_reconstructs_cut_polynomial(budget):
    power=(Q(1,7),Q(-3,5),Q(4,3),Q(-2,9),Q(1,11),Q(-5,4))
    for cut in (Q(1),Q(1,4),Q(3,4)):
        controls=_bernstein_controls(power,cut)
        reconstructed=[Q(0)]*6
        for j,b in enumerate(controls):
            for k in range(6-j):
                reconstructed[j+k]+=b*comb(5,j)*comb(5-j,k)*(-1)**k
        assert tuple(reconstructed)==tuple(v*cut**i for i,v in enumerate(power))
        budget.check(1024)


def test_quiet_cap_has_geometric_metric_floor_and_zero_mass(budget):
    basis=FixedCapBasis([0,Q(1,4),Q(1,2)],budget=budget)
    zero=np.zeros(basis.dimension);Z=Q(3,8)
    result=fixed_metric_bound(basis,zero,zero,Z,alpha=.5,budget=budget)
    with ctx.workprec(160):
        expected=1/_ball(Z).cosh()**2
        assert _ball(result['A_lower_float_down'])<=expected
        assert expected-_ball(result['A_lower_float_down'])<arb('1e-14')
    assert result['mass_upper']==0.
    assert result['interval_count']==2 and result['wall_cut_fraction']==.5


def test_geometric_linear_field_mass_is_below_source_envelope(budget):
    basis=FixedCapBasis([0,Q(1,8),Q(1,4),Q(1,2)],budget=budget)
    c=linear_coefficients(basis,Q(1,2));d=linear_coefficients(basis,Q(3,4))
    c0,d0=c.copy(),d.copy();Z=Q(3,8)
    result=fixed_metric_bound(basis,c,d,Z,alpha=0.,budget=budget)
    with ctx.workprec(160):
        z=_ball(Z);r=z.tanh();A0=1/z.cosh()**2
        # Independent integration: 4N=integral(chi_z^2+k^2)dz-A0*chi_W^2/R.
        mass=(_ball(Q(1,4))*z+_ball(Q(9,16))*z**3/3-A0*_ball(Q(1,4))*z*z/r)/4
        assert mass>0 and mass<_ball(result['mass_upper'])
        assert _ball(result['A_lower_float_down'])<=A0
    np.testing.assert_array_equal(c,c0);np.testing.assert_array_equal(d,d0)


def test_finite_coupling_and_cut_domain_have_positive_conservative_margin(budget):
    basis=FixedCapBasis([0,Q(1,4),Q(1,2)],budget=budget)
    c=linear_coefficients(basis,Q(1,8));d=linear_coefficients(basis,Q(1,4))
    small=fixed_metric_bound(basis,c,d,Q(5,16),alpha=.001,budget=budget)
    large=fixed_metric_bound(basis,c,d,Q(7,16),alpha=.001,budget=budget)
    assert 0<large['A_lower_float_down']<small['A_lower_float_down']<1
    assert 0<small['mass_upper']<large['mass_upper']
    assert small['wall_cut_fraction']==.25 and large['wall_cut_fraction']==.75


def test_certificate_rejects_unresolved_margin_and_bad_inputs(budget):
    basis=FixedCapBasis([0,Q(1,2)],budget=budget)
    c=linear_coefficients(basis,Q(1));zero=np.zeros(basis.dimension)
    with pytest.raises(ValueError,match='physical sign is undetermined'):
        fixed_metric_bound(basis,c,zero,Q(1,4),alpha=10.,budget=budget)
    for Z in (Q(1,16),Q(7,16)):
        with pytest.raises(ValueError,match='wall cut'):
            fixed_metric_bound(basis,zero,zero,Z,alpha=0.,budget=budget)
    with pytest.raises(ValueError,match='dimensions'):
        fixed_metric_bound(basis,[0],zero,Q(1,4),alpha=0.,budget=budget)
    with pytest.raises(ValueError):
        fixed_metric_bound(basis,c,zero,Q(1,4),alpha=-1.,budget=budget)
    with pytest.raises(ValueError):
        fixed_metric_bound(basis,c,zero,Q(1,4),alpha=np.nan,budget=budget)
