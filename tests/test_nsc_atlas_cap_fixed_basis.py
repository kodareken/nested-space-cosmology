"""Exact local polynomials, independent Gram integrals and initial data."""
from fractions import Fraction as Q

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


CONTROLS=ControlBudget(maximum_evaluations=2_000_000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed local basis controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def dense(lower):
    n=lower.shape[1];result=np.zeros((n,n))
    for d in range(len(lower)):
        for j in range(n-d): result[j+d,j]=result[j,j+d]=lower[d,j]
    return result


def exact_piece(basis,cell,coefficients):
    return tuple(sum(row[j]*coefficients[int(i)] for j,i in enumerate(basis.indices(cell)))
                 for row in basis.exact_shapes(cell))


def evaluate_exact(poly,x,derivative=0,h=Q(1)):
    for _ in range(derivative): poly=tuple(p*poly[p]/h for p in range(1,len(poly)))
    return sum(c*x**p for p,c in enumerate(poly))


def profile(kind='piecewise'):
    labels=[Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)]
    def jet(u):
        if kind=='global': return u**5-u**3+u*u+2*u,5*u**4-3*u*u+2*u+2,20*u**3-6*u+2
        x=max(u-Q(1,4),Q(0))
        return x**3+x**5,3*x*x+5*x**4,6*x+20*x**3
    rows=[jet(u) for u in labels]
    return RetainedCapBackground([float(u) for u in labels],
        *[[float(row[i]) for row in rows] for i in range(3)])


def test_minimal_cell_indices_share_only_the_internal_node(budget):
    one=FixedCapBasis([0,Q(1,2)],budget=budget)
    two=FixedCapBasis([0,Q(1,4),Q(1,2)],budget=budget)
    assert one.dimension==3 and list(one.indices(0))==[0,1,2]
    assert two.dimension==7 and list(two.indices(1))==[1,2,3,4,5,6]


def test_independent_bernstein_source_matches_known_polynomial(budget):
    from check_nsc_atlas_cap_fixed_transfer import BernsteinSource
    source=BernsteinSource(profile('global'),budget)
    for u in (Q(0),Q(1,7),Q(1,4),Q(4,7),Q(1)):
        assert source.evaluate(u,budget)==(u**5-u**3+u*u+2*u,
                                           5*u**4-3*u*u+2*u+2,20*u**3-6*u+2)
    with pytest.raises(ValueError,match='exact source domain'):
        source.evaluate(Q(-1,7),budget)
    with pytest.raises(ValueError,match='exact source domain'):
        source.evaluate(Q(8,7),budget)


def test_independent_transfer_comparison_covers_piecewise_and_quiet_source(budget):
    from check_nsc_atlas_cap_fixed_transfer import BernsteinSource, compare_transfer
    background=profile()
    source=BernsteinSource(background,budget)
    for u in (Q(-1,7),Q(0),Q(1,4),Q(1,3),Q(7,8),Q(1)):
        x=max(u-Q(1,4),Q(0))
        assert source.evaluate(u,budget)==(x**3+x**5,3*x*x+5*x**4,6*x+20*x**3)
    transfer=transfer_initial_background(background,Q(1,8),Q(7,8),budget=budget)
    result=compare_transfer(background,transfer,budget)
    assert result['evaluation_points']==6*transfer['basis'].cells+2
    assert result['exact_polynomial_identities']==4*transfer['basis'].cells
    assert result['exact_point_equalities']==4*result['evaluation_points']
    assert set(result['boundary_targets'])=={'axis','wall'}


def test_basis_is_c1_and_axis_polynomial_is_odd_for_every_state(budget):
    basis=FixedCapBasis([0,Q(1,8),Q(1,4),Q(1,2)],budget=budget)
    c=tuple(Q(i-3,7) for i in range(basis.dimension))
    pieces=[exact_piece(basis,i,c) for i in range(basis.cells)]
    assert all(pieces[0][p]==0 for p in (0,2,4))
    for i in range(basis.cells-1):
        hl=basis.exact_knots[i+1]-basis.exact_knots[i]
        hr=basis.exact_knots[i+2]-basis.exact_knots[i+1]
        assert evaluate_exact(pieces[i],Q(1))==evaluate_exact(pieces[i+1],Q(0))
        assert evaluate_exact(pieces[i],Q(1),1,hl)==evaluate_exact(pieces[i+1],Q(0),1,hr)
    assert basis.dimension==11


def test_local_values_and_first_two_derivatives_match_exact_polynomials(budget):
    basis=FixedCapBasis([0,Q(1,7),Q(2,7),Q(5,7)],budget=budget)
    c=tuple(Q(i+1,11) for i in range(basis.dimension))
    for cell,x in enumerate((Q(1,9),Q(2,5),Q(3,7))):
        row=basis.values(cell,float(x),budget=budget)
        poly=exact_piece(basis,cell,c)
        h=basis.exact_knots[cell+1]-basis.exact_knots[cell]
        local=np.array([float(c[int(i)]) for i in row['indices']])
        actual=[row[key]@local for key in ('value','first','second')]
        expected=[float(evaluate_exact(poly,x,order,h)) for order in range(3)]
        np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=2e-13)


def test_clipped_gram_and_solve_match_independent_exact_power_integrals(budget):
    basis=FixedCapBasis([0,Q(1,8),Q(1,4),Q(1,2)],budget=budget)
    Z=Q(3,8);n=basis.dimension
    exact=[[Q(0) for _ in range(n)] for _ in range(n)]
    for cell in range(basis.cells):
        low,high=basis.exact_knots[cell:cell+2];h=high-low
        cut=min(Q(1),(Z-low)/h)
        shapes=basis.exact_shapes(cell);idx=basis.indices(cell)
        for a,ga in enumerate(idx):
            for b,gb in enumerate(idx):
                exact[int(ga)][int(gb)]+=h*sum(shapes[p][a]*shapes[q][b]*cut**(p+q+1)/Q(p+q+1)
                    for p in range(6) for q in range(6))
    budget.check(12_000)
    operator=basis.gram(Z,budget=budget)
    expected=np.array(exact,dtype=float)
    np.testing.assert_allclose(dense(operator.lower),expected,rtol=5e-14,atol=2e-17)
    true=np.arange(1,n+1,dtype=float)/13
    actual=operator.solve(expected@true,budget=budget)
    np.testing.assert_allclose(actual,true,rtol=2e-8,atol=2e-9)
    assert np.isfinite(operator.condition) and operator.condition>1


def test_same_fixed_fields_are_evaluated_without_a_wall_dependent_basis(budget):
    basis=FixedCapBasis([0,Q(1,4),Q(3,4)],budget=budget)
    before=basis.values(1,.3,budget=budget)
    first=basis.gram(Q(1,2),budget=budget)
    second=basis.gram(Q(9,16),budget=budget)
    after=basis.values(1,.3,budget=budget)
    for key in before: np.testing.assert_array_equal(before[key],after[key])
    assert first.Z!=second.Z
    assert not np.array_equal(first.lower,second.lower)


def test_transfer_reconstructs_a_global_polynomial_and_axis_limits_exactly(budget):
    bg=profile('global');T,Z=Q(1,2),Q(1,4)
    result=transfer_initial_background(bg,T,Z,budget=budget)
    basis=result['basis']
    def F(u): return u**5-u**3+u*u+2*u
    def Fp(u): return 5*u**4-3*u*u+2*u+2
    for cell in range(basis.cells):
        low,high=basis.exact_knots[cell:cell+2]
        x=Q(1,3);z=low+(high-low)*x
        if z>Z: continue
        assert evaluate_exact(result['exact_polynomials']['chi'][cell],x)==F(T+z)-F(T-z)
        assert evaluate_exact(result['exact_polynomials']['canonical_density'][cell],x)==Fp(T+z)-Fp(T-z)
    h=basis.exact_knots[1]
    assert result['exact_chi'][0]/h==2*Fp(T)
    assert all(result['exact_polynomials']['canonical_density'][0][p]==0 for p in (0,2,4,5))


def test_piecewise_initial_momentum_needs_c1_not_c2_and_support_is_not_source_extrapolation(budget):
    bg=profile();T,Z=Q(1,8),Q(7,8)
    result=transfer_initial_background(bg,T,Z,budget=budget)
    basis=result['basis']
    assert T+basis.exact_knots[-1]>Q(1)
    assert (Z-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2])==Q(1,2)
    # The activation knot z=1/8 has zero value/slope but a momentum
    # second-derivative jump of 6. Its C1 representation is exact.
    knot=basis.exact_knots.index(Q(1,8))
    right=result['exact_polynomials']['canonical_density'][knot]
    h=basis.exact_knots[knot+1]-basis.exact_knots[knot]
    assert right[0]==right[1]==0
    assert 2*right[2]/(h*h)==6
    for field,key in (('chi','exact_chi'),('canonical_density','exact_canonical_density')):
        for i in range(basis.cells):
            assert exact_piece(basis,i,result[key])==result['exact_polynomials'][field][i]
    with pytest.raises(ValueError,match='finite labels'):
        bg.reduced_jets(float(T),float(basis.exact_knots[-1]),budget=budget)


def test_conversion_bound_covers_rounded_coefficients_in_the_exact_space(budget):
    result=transfer_initial_background(profile(),.123456789,.3333333333,budget=budget)
    basis=result['basis']
    for name,key in (('chi','exact_chi'),('canonical_density','exact_canonical_density')):
        errors=tuple(Q(float(v))-e for v,e in zip(result[name],result[key]))
        bounds=result['coefficient_conversion_bounds'][name]
        for cell in range(basis.cells):
            poly=exact_piece(basis,cell,errors)
            h=basis.exact_knots[cell+1]-basis.exact_knots[cell]
            for x in (Q(0),Q(1,7),Q(1,2),Q(1)):
                assert abs(evaluate_exact(poly,x))<=Q(bounds['value'])
                assert abs(evaluate_exact(poly,x,1,h))<=Q(bounds['first_derivative'])
    assert result['numerical_evaluation_enclosure'] is False
    assert result['physical_evolution'] is False


def test_quiet_transfer_and_copied_data_remain_unchanged(budget):
    bg=RetainedCapBackground([0.,1.],[0.,0.],[0.,0.],[0.,0.])
    before=bg.inputs
    result=transfer_initial_background(bg,.2,.3,budget=budget)
    assert np.all(result['chi']==0) and np.all(result['canonical_density']==0)
    result['chi'][0]=1
    for key in before: np.testing.assert_array_equal(bg.inputs[key],before[key])
    knots=result['basis'].knots;knots[-1]=99
    assert result['basis'].knots[-1]!=99


def test_incompatible_axis_parity_and_exact_domain_escape_are_rejected(budget):
    # At T=1/4, the merely C2 cubic activation has quadratic momentum
    # at the axis. It cannot be smoothed into an odd polynomial silently.
    with pytest.raises(ValueError,match='axis jets'):
        transfer_initial_background(profile(),Q(1,4),Q(1,4),budget=budget)
    bg=RetainedCapBackground([0.,1.],[0.,0.],[0.,0.],[0.,0.])
    assert .1+.9==1.
    with pytest.raises(ValueError,match='exact supplied'):
        transfer_initial_background(bg,.1,.9,budget=budget)


def test_basis_domain_budget_and_solve_dimensions_are_bounded(budget):
    for knots in ([.1,.2],[0.,0.],[0.,1.,.5],[0.,1e-324]):
        with pytest.raises(ValueError): FixedCapBasis(knots,budget=budget)
    basis=FixedCapBasis([0.,.5],budget=budget)
    for Z in (.01,.49,.6):
        with pytest.raises(ValueError,match='wall cut'): basis.gram(Z,budget=budget)
    with pytest.raises(ValueError,match='dimension'):
        basis.power_coefficients([0.],budget=budget)
    with pytest.raises(ValueError,match='coordinate'):
        basis.values(0,1.1,budget=budget)
    operator=basis.gram(.25,budget=budget)
    with pytest.raises(ValueError,match='dimension'): operator.solve([1.],budget=budget)
    with pytest.raises(RuntimeError,match='budget'):
        FixedCapBasis([0.,.5],budget=Budget(maximum_evaluations=1,cpu_seconds=1.))
