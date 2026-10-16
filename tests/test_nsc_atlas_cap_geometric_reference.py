"""Regular fixed-cap controls; no retained fields or physical evolution."""
import copy

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from numpy.polynomial.legendre import leggauss

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_geometric_reference as model

CONTROLS = ControlBudget(maximum_evaluations=200000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("geometric cap controls:", CONTROLS.budget.evaluations,
          "coefficient/manufactured values;", CONTROLS.active_cpu_seconds,
          "active CPU seconds")


def evaluate(budget, **kwargs):
    args = dict(z_nodes=np.array([0., .07, .21, .5]),
                chi=np.array([0., .03, .12, .13]),
                velocity=np.array([0., -.02, .07, .08]),
                T_x=.7, incoming_w_plus=.4)
    args.update(kwargs)
    return model.geometric_cap_reference(**args, budget=budget)


@pytest.mark.parametrize("nodes", [[0., .6], [0., .03, .11, .39, .6]])
def test_exact_time_linear_wave_and_original_incoming_conversion(nodes, budget):
    z = np.array(nodes)
    a, b, T, tx = .4, -.03, .7, 1.2
    coefficient = a+b*T
    chi, velocity = coefficient*z, b*z
    R, A = np.tanh(z[-1]), 1/np.cosh(z[-1])**2
    incoming = (velocity[-1]+coefficient)/(A*R)-chi[-1]/R**2
    row = evaluate(budget, z_nodes=z, chi=chi, velocity=velocity,
                   incoming_w_plus=incoming, T_x=tx)
    np.testing.assert_allclose(row['rates']['chi'], tx*velocity, rtol=0, atol=1e-16)
    np.testing.assert_allclose(row['rates']['velocity'], 0, rtol=0, atol=3e-14)
    expected = b*b*z[-1]**3/12+row['positive_trace_gap']*coefficient**2*z[-1]**2/4
    assert row['mass'] == pytest.approx(expected, rel=2e-15, abs=1e-17)
    expected_rate = tx*b*coefficient*row['positive_trace_gap']*z[-1]**2/2
    assert row['mass_rate'] == pytest.approx(expected_rate, abs=1e-16)
    assert abs(row['wall']['incoming_discrepancy']) < 1e-15
    assert row['centre'] == pytest.approx(dict(f=coefficient, Pi=b, Phi=0.))


def test_independent_dense_element_matrices_and_mass_gradient(budget):
    row = evaluate(budget)
    z, chi, u = row['z_nodes'], row['chi'], row['velocity']
    M, K = np.zeros((len(z), len(z))), np.zeros((len(z), len(z)))
    budget.check(16*(len(z)-1))
    for j, h in enumerate(np.diff(z)):
        M[j:j+2, j:j+2] += h*np.array([[2., 1.], [1., 2.]])/6
        K[j:j+2, j:j+2] += np.array([[1., -1.], [-1., 1.]])/h
    S = K.copy()
    S[-1, -1] -= row['wall']['A']/row['radius']
    expected = (u@M@u+chi@S@chi)/4
    assert row['mass'] == pytest.approx(expected, abs=2e-17)
    np.testing.assert_allclose(row['gradient']['chi'], (S@chi)[1:]/2, rtol=1e-14, atol=1e-16)
    np.testing.assert_allclose(row['gradient']['velocity'], (M@u)[1:]/2, rtol=1e-14, atol=1e-16)
    rhs = -(S@chi)[1:]
    rhs[-1] += row['wall']['A']*row['radius']*.4-u[-1]
    acceleration = np.linalg.solve(M[1:, 1:], rhs)
    np.testing.assert_allclose(row['rates']['velocity'][1:], .7*acceleration, rtol=1e-14, atol=1e-15)
    assert abs(row['balance_defect']) < 2e-17
    assert row['linear_solve_defect'] < 2e-16
    # Weak wall flux is exposed separately from the reconstructed derivative.
    assert abs(row['wall']['incoming_discrepancy']) > .1
    assert row['wall']['incoming_discrepancy'] == pytest.approx(
        row['wall']['derivative_discrepancy']/(row['wall']['A']*row['radius']))


def test_independent_original_radial_density_integral(budget):
    row = evaluate(budget)
    z, chi, u = row['z_nodes'], row['chi'], row['velocity']
    nodes, weights = leggauss(16)
    total = 0.
    for j, (lo, hi) in enumerate(zip(np.tanh(z[:-1]), np.tanh(z[1:]))):
        budget.check(9*len(nodes))
        r = (lo+hi)/2+(hi-lo)*nodes/2
        zr, A = np.arctanh(r), 1-r*r
        h = z[j+1]-z[j]
        slope = (chi[j+1]-chi[j])/h
        field = chi[j]+slope*(zr-z[j])
        velocity = u[j]+(u[j+1]-u[j])*(zr-z[j])/h
        f = field/r
        pi = velocity/(A*r)
        phi = (slope/A-f)/r
        density = r*r*f*f/2+r*r*A*(pi*pi+phi*phi)/4
        total += (hi-lo)*(weights@density)/2
    assert total == pytest.approx(row['mass'], rel=3e-13, abs=1e-16)


def test_homogeneous_incoming_mass_decay_and_positive_trace_control(budget):
    row = evaluate(budget, incoming_w_plus=0.)
    u = row['wall']['u']
    assert row['mass_rate'] == pytest.approx(-row['T_x']*u*u/2, abs=2e-17)
    assert row['incoming_flux_bound'] == 0
    assert row['mass'] > 0
    assert 0 < row['coercivity_factor'] < 1
    gradient_energy = np.sum(np.diff(row['chi'])**2/np.diff(row['z_nodes']))/4
    assert row['mass'] >= row['kinetic_mass']+row['coercivity_factor']*gradient_energy
    driven = evaluate(budget)
    assert driven['mass_rate'] <= driven['incoming_flux_bound']


def test_regular_centre_and_trace_gap_by_independent_symbolic_limits(budget):
    budget.check(16)
    r = s.symbols('r', positive=True)
    c, d = s.symbols('c d')
    A = 1-r*r
    chi, u = c*s.atanh(r), d*s.atanh(r)
    f, pi = chi/r, u/(A*r)
    phi = (s.diff(chi, r)-f)/r
    assert s.limit(f, r, 0) == c
    assert s.limit(pi, r, 0) == d
    assert s.limit(phi, r, 0) == 0
    assert s.limit(phi/r, r, 0) == 2*c/3
    numerator = r-A*s.atanh(r)
    assert s.limit(numerator, r, 0) == 0
    assert s.simplify(s.diff(numerator, r)-2*r*s.atanh(r)) == 0


def test_incoming_mass_estimate_is_an_exact_completed_square(budget):
    budget.check(8)
    u, drive = s.symbols('u drive', real=True)
    flux = u*(drive-u)/2
    assert s.expand(flux-drive**2/8+(u-drive/2)**2/2) == 0


def test_quiet_cap_has_exact_zero_rhs_mass_and_boundary_flux(budget):
    row = evaluate(budget, chi=np.zeros(4), velocity=np.zeros(4), incoming_w_plus=0.)
    assert row['mass'] == row['mass_rate'] == row['physical_mass_flux'] == 0
    np.testing.assert_array_equal(row['rates']['chi'], np.zeros(4))
    np.testing.assert_array_equal(row['rates']['velocity'], np.zeros(4))
    assert row['wall']['incoming_discrepancy'] == 0


@pytest.mark.parametrize("z", [1e-8, .5, .75])
def test_positive_gap_against_high_precision_at_distinct_arithmetic_regimes(z, budget):
    budget.check(32)
    with mp.workdps(65):
        value = mp.mpf(z)
        exact = float(1/value-2/mp.sinh(2*value))
    assert model._positive_gap(z) == pytest.approx(exact, rel=5e-16, abs=0.)


def test_small_radius_linear_shape_retains_mass_without_catastrophic_subtraction(budget):
    z = np.array([0., 1e-8])
    c = .4
    row = evaluate(budget, z_nodes=z, chi=c*z, velocity=np.zeros(2), incoming_w_plus=0.)
    expected = c*c*z[-1]**3/6  # leading term; next relative term is O(z^2)
    assert row['mass'] == pytest.approx(expected, rel=5e-16, abs=0.)
    assert row['mass'] > 0


def test_inputs_and_returned_arrays_do_not_alias(budget):
    z, chi, u = np.array([0., .2, .5]), np.array([0., .1, .3]), np.array([0., .2, .1])
    before = copy.deepcopy((z, chi, u))
    row = evaluate(budget, z_nodes=z, chi=chi, velocity=u)
    for original, saved in zip((z, chi, u), before):
        np.testing.assert_array_equal(original, saved)
    for key, original in [('z_nodes', z), ('chi', chi), ('velocity', u)]:
        assert not np.shares_memory(row[key], original)


@pytest.mark.parametrize("changed", [
    dict(z_nodes=[0.]), dict(z_nodes=[0., .2, .2, .5]),
    dict(z_nodes=[.01, .07, .21, .5]), dict(z_nodes=[0., .07, .21, 100.]),
    dict(chi=[0., 1.]), dict(chi=[.01, .03, .12, .13]),
    dict(velocity=[.01, -.02, .07, .08]), dict(velocity=[0., np.nan, .07, .08]),
    dict(T_x=0.), dict(incoming_w_plus=True),
])
def test_invalid_data_stops_before_the_linear_solve(changed, budget, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('invalid data reached the linear solve')
    monkeypatch.setattr(model, 'solve_banded', forbidden)
    with pytest.raises(ValueError):
        evaluate(budget, **changed)


def test_shared_budget_stops_before_the_linear_solve(budget, monkeypatch):
    previous = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    def forbidden(*args, **kwargs):
        raise AssertionError('exhausted budget reached the linear solve')
    monkeypatch.setattr(model, 'solve_banded', forbidden)
    try:
        with pytest.raises(RuntimeError, match='budget exhausted'):
            evaluate(budget)
    finally:
        budget.maximum_evaluations = previous
