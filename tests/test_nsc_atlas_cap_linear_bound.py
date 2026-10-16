"""Independent monomial integrals and a two-cell polynomial control."""
from fractions import Fraction as Q
import sympy as s
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_linear_bound import frozen_cap_bound

CONTROLS=ControlBudget(maximum_evaluations=60000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('linear bound controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_one_cell_bound_is_independent_of_the_odd_quintic_basis(budget):
    basis=FixedCapBasis([0,1],budget=budget)
    row=frozen_cap_bound(basis,Q(1,2),budget=budget)
    assert row['bound_squared']==3444 and row['integer_bound']==59
    assert row['M']==row['M'].T and row['K']==row['K'].T
    for n in range(1,4):
        assert row['M'][:n,:n].det()>0 and row['K'][:n,:n].det()>0


def test_two_cell_matrices_integrate_one_global_linear_field(budget):
    basis=FixedCapBasis([0,Q(1,2),1],budget=budget);Z=Q(3,4)
    row=frozen_cap_bound(basis,Z,budget=budget)
    # On the odd cell, z is represented by (h,h,h). On the next full
    # quintic cell, endpoint values and scaled derivatives represent z.
    c=s.Matrix([s.Rational(1,2)]*3+[0,0,1,s.Rational(1,2)])
    assert (c.T*row['M']*c)[0]==s.Rational(9,64) # integral_0^.75 z^2
    assert (c.T*row['K']*c)[0]==s.Rational(3,4)
    assert (row['e'].T*c)[0]==s.Rational(3,4)
    assert row['integer_bound']**2>=row['bound_squared']
    assert (row['integer_bound']-1)**2<row['bound_squared']
    with pytest.raises(ValueError,match='wall cut'):
        frozen_cap_bound(basis,Q(1,2),budget=budget)
