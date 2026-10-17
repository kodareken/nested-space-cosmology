"""Positive comparison controls without an original-source application."""
from fractions import Fraction as Q

from flint import arb,ctx
import pytest

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_geometric_departure_bound as m

CONTROLS=ControlBudget(maximum_evaluations=12000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('geometric departure controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_triangular_inverse_against_back_substitution(budget):
    matrix=[[0,Q(1,4),0],[0,0,Q(1,3)],[0,0,Q(1,2)]]
    result=m.positive_resolvent(matrix,[1,2,3],[1,1,1],budget=budget)
    assert result['solution']==[2,4,6]
    assert result['weighted_norm']==Q(1,2)


def test_cycle_against_exact_geometric_series(budget):
    # Each three turns of the cycle supplies another factor 1/24.
    matrix=[[0,Q(1,2),0],[0,0,Q(1,3)],[Q(1,4),0,0]]
    result=m.positive_resolvent(matrix,[1,0,0],[1,1,1],budget=budget)
    assert result['solution']==[Q(24,23),Q(2,23),Q(6,23)]
    assert result['weighted_norm']==Q(1,2)


def test_quiet_forcing(budget):
    result=m.positive_resolvent([[Q(1,4),0,0],[Q(1,3),0,0],[0,Q(1,5),0]],
                               [0,0,0],[1,2,3],budget=budget)
    assert result['solution']==[0,0,0]


@pytest.mark.parametrize('entry,match',[(1,'norm'),(-1,'nonnegative')])
def test_noncontractive_or_negative_matrix_stops(entry,match,budget):
    with pytest.raises(ValueError,match=match):
        m.positive_resolvent([[entry,0,0],[0,0,0],[0,0,0]],[1,1,1],[1,1,1],budget=budget)


def test_coefficient_extraction_preserves_affine_forcing(budget):
    budget.check(1024)
    with ctx.workprec(160):
        c={k:arb(1) for k in ('dt','W','Wxhi','Wxlo','Upmin','Upmax','Vpmin','Vpmax',
                              'Amin','Hmin','rx0hi','Sw','drt','dSw','u_wall','nc','bbar',
                              'mu_coefficient','N','vmax','pabs','k','Cmax','kinetic')}
        c.update(rlo=m._ball(Q(1,4)),rhi=m._ball(Q(1,2)),Rlo=m._ball(Q(1,2)),
                 R0lo=m._ball(Q(1,2)),R0hi=m._ball(Q(1,2)))
        c['source']={k:m._ball(Q(1,10)) for k in ('r','log_C','log_C_x','r_times_r_x_minus_x')}
        zero=m.affine_image(c,(0,0,0),0)
        assert all(v==0 for v in zero.values())
        f={k:m._round_upper(v) for k,v in m.affine_image(c,(0,0,0),1).items()}
        cols=[{k:m._round_upper(v) for k,v in m.affine_image(c,tuple(int(i==j) for i in range(3)),0).items()}
              for j in range(3)]
        vector,scale=(Q(2,3),Q(3,5),Q(5,7)),Q(1,11)
        direct=m.affine_image(c,vector,scale)
        for k in m.OUTPUTS:
            upper=f[k]*scale+sum(cols[j][k]*vector[j] for j in range(3))
            assert m.dyadic(direct[k].upper())<=upper,k
        # Scaling forcing affects source and explicit-alpha slots together;
        # it does not replace the frozen scalar or metric coefficients.
        assert all(v>=0 for v in f.values())


def test_budget_stops_before_frozen_context(budget):
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.geometric_departure_bound({}, {}, {}, budget=Spent())
