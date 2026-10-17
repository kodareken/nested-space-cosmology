"""Manufactured point rows for the compiled fixed-cap coefficient kernel."""
import inspect

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground, SMALL_Z
from recursive_horizons.nsc_atlas_cap_compiled_coefficients import CompiledCapCoefficients
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import _Fields, fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_panel_jets import history_piece_indices, panel_profile_jet
from recursive_horizons.nsc_atlas_cap_potential_pair import reference_potential_pair
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

CONTROLS = ControlBudget(maximum_evaluations=2_000_000, cpu_seconds=120.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


def reference_row(fields, alpha, count, cell, xi, mass, budget, background=None, time=0.,
                  indices=None):
    row = fields.evaluate(cell, xi, budget)
    if background is not None:
        jet = background.reduced_jets(time, row['z'], budget=budget, polynomial_indices=indices)
        for key in ('f', 'p', 'h', 'chi', 'k', 'chi_z'):
            row[key] += jet[key]
        row['background'] = jet
    partials = np.zeros(2*count+(background is not None))
    if row['r'] == 0:
        return (row['r'], fields.width[cell]*row['A0'], 1., 0., 0., partials)
    f, p, h = row['f'], row['p'], row['h']
    P, U = reference_potential_pair(f, alpha, budget=budget)
    A = row['A0']-2*alpha*mass/row['r']
    square = p*p+h*h
    a = alpha*square/(2*row['r'])
    rho = row['r']**2*P/2+A*square/4
    idx = row['indices']
    partials[idx] = row['r']**2*U*row['fb']/2+A*h*row['hb']/2
    partials[count+idx] = A*p*row['pb']/2
    if background is not None:
        jet = row['background']
        partials[-1] = row['r']**2*U*jet['f_T']/2+A*(p*jet['p_T']+h*jet['h_T'])/2
    return (row['r'], fields.width[cell]*row['A0'], A, a, rho, partials)


def cache_for(basis, chi, density, breaks, owners, alpha, budget, background=None, time=0.,
              panel_indices=None):
    fields = _Fields(basis, chi, density, budget)
    return fields, CompiledCapCoefficients(
        fields, breaks, owners, cut=0.5, alpha=alpha, count=basis.dimension, background=background,
        time=time, panel_indices=panel_indices, budget=budget)


def test_default_keeps_the_python_callback():
    parameter = inspect.signature(fixed_cap_reference).parameters['compiled_coefficients']
    assert parameter.default is False
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY


def test_interior_rows_match_sparse_zero_and_positive_amplitude(budget):
    basis = FixedCapBasis([0, .3, .9], budget=budget)
    chi = np.zeros(basis.dimension)
    density = np.zeros(basis.dimension)
    chi[0] = .02
    density[4] = -.01
    breaks = np.array([1., 2.])
    owners = np.array([1])
    xi = .4
    for alpha in (0., .003):
        fields, compiled = cache_for(basis, chi, density, breaks, owners, alpha, budget)
        got = compiled.evaluate(0, 1.+xi, .2, budget)
        expected = reference_row(fields, alpha, basis.dimension, 1, xi, .2, budget)
        assert got[0] == pytest.approx(expected[0], rel=0, abs=1e-15)
        assert got[1] == pytest.approx(expected[1], rel=0, abs=1e-15)
        assert got[2] == pytest.approx(expected[2], rel=0, abs=1e-15)
        assert got[3] == pytest.approx(expected[3], rel=0, abs=1e-15)
        assert got[4] == pytest.approx(expected[4], rel=0, abs=1e-15)
        assert got[5] == pytest.approx(expected[5], rel=0, abs=1e-14)
        assert np.count_nonzero(got[5]) < got[5].size
        assert got[5][4] != 0 or got[5][basis.dimension+4] != 0


def test_background_interior_and_knot_side_keep_panel_ownership(budget):
    basis = FixedCapBasis([0, .2, 1.], budget=budget)
    chi = np.zeros(basis.dimension)
    density = np.zeros(basis.dimension)
    chi[3] = .01
    labels = np.array([-1., 0., .45, 1.])
    values = np.array([0., .1, -.2, .3])
    first = np.array([0., .2, -.4, .1])
    second = np.array([0., 0., 1., -.5])
    profile = CharacteristicProfile(labels, values, first, second)
    background = RetainedCapBackground(labels, values, first, second)
    time = 0.
    z = float(np.nextafter(.45, 0.))
    xi = (z-.2)/.8
    point = 1.+xi
    actual = history_piece_indices(profile, np.array([-z, z]), budget=budget)
    selected = actual.copy()
    selected[1] = actual[1]+1
    fields, compiled = cache_for(
        basis, chi, density, np.array([1., 2.]), np.array([1]), .002, budget,
        background=background, time=time, panel_indices=[selected])
    got = compiled.evaluate(0, point, .05, budget)
    expected = reference_row(
        fields, .002, basis.dimension, 1, xi, .05, budget, background=background, time=time,
        indices=selected)
    assert got[5] == pytest.approx(expected[5], rel=0, abs=2e-14)
    assert got[2] == pytest.approx(expected[2], rel=0, abs=2e-14)
    opposite = selected.copy()
    opposite[1] = actual[1]
    sided = panel_profile_jet(profile, np.array([-z, z]), selected, budget=budget)
    plain = panel_profile_jet(profile, np.array([-z, z]), opposite, budget=budget)
    assert sided['F_third'][1] != pytest.approx(plain['F_third'][1], rel=0, abs=0)
    assert z > SMALL_Z


def test_near_centre_background_stays_on_the_reference_path(budget):
    basis = FixedCapBasis([0, .3, .9], budget=budget)
    labels = np.array([-1., 0., 1.])
    zeros = np.zeros(3)
    background = RetainedCapBackground(labels, zeros, zeros, zeros.copy())
    fields, compiled = cache_for(
        basis, np.zeros(basis.dimension), np.zeros(basis.dimension), np.array([0., 1.]),
        np.array([0]), 0., budget, background=background, time=0.)
    assert compiled.evaluate(0, .01, 0., budget) is None
    xi = .01/.3
    assert fields.left[0]+fields.width[0]*xi <= SMALL_Z


def test_batch_matches_scalar_rows_and_rejects_bad_input(budget):
    basis = FixedCapBasis([0, .3, .9], budget=budget)
    chi = np.zeros(basis.dimension)
    chi[0] = -.004
    density = np.zeros(basis.dimension)
    density[2] = .006
    breaks = np.array([1., 1.2, 2.])
    owners = np.array([1, 1])
    fields, compiled = cache_for(basis, chi, density, breaks, owners, .001, budget)
    panels = np.array([0, 1, 1])
    points = np.array([1.05, 1.4, 1.8])
    masses = np.array([.1, -.2, .3])
    batch = compiled.evaluate_batch(panels, points, masses, budget)
    for row, panel, point, mass in zip(batch, panels, points, masses):
        scalar = compiled.evaluate(int(panel), float(point), float(mass), budget)
        assert row[0] == scalar[0]
        assert row[5:] == pytest.approx(scalar[5], rel=0, abs=0)
    with pytest.raises(ValueError, match='escaped'):
        compiled.evaluate(2, 1.1, 0., budget)
    with pytest.raises(ValueError, match='positive metric'):
        compiled.evaluate(1, 1.5, 1e6, budget)
    with pytest.raises(ValueError, match='escaped'):
        compiled.evaluate(0, np.nan, 0., budget)
    with pytest.raises(ValueError, match='nonfinite'):
        compiled.evaluate_batch(panels, points, np.array([.1, np.nan, .3]), budget)
    with pytest.raises(ValueError, match='nonfinite'):
        compiled.evaluate(1, 1.5, np.nan, budget)
    with pytest.raises(ValueError, match='integer'):
        compiled.evaluate_batch(panels+.1, points, masses, budget)
    labels = np.array([-1., 0., 1.])
    zeros = np.zeros(3)
    background = RetainedCapBackground(labels, zeros, zeros, zeros.copy())
    _, centre = cache_for(
        basis, chi, density, np.array([0., 1.]), np.array([0]), .001, budget,
        background=background, time=0.)
    with pytest.raises(ValueError, match='near-centre'):
        centre.evaluate_batch(np.array([0]), np.array([.01]), np.array([0.]), budget)
    before = budget.evaluations
    compiled.evaluate(1, 1.5, .2, budget)
    assert budget.evaluations > before
    del fields
