"""Synthetic endpoint options, grid legality and time partition; no physical solve."""
import math
from fractions import Fraction as Q

import numpy as np
import pytest

import run_nsc_atlas_fixed_endpoint as runner
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_endpoint_comparison import jordan_radius


X = np.linspace(runner.X_ANCHOR, runner.X_LAST, runner.SOURCE_ROWS)
FRONT = .42269179059454887
BACK = .4631925883576189


def test_levels_refine_space_and_time_together():
    assert [runner.level_spec(level)['cap_cells'] for level in range(3)] == [2, 4, 8]
    assert [runner.level_spec(level)['exterior_cells'] for level in range(3)] == [32, 64, 128]
    assert [runner.level_spec(level)['step_budget'] for level in range(3)] == [2000, 4000, 8000]
    assert [runner.level_spec(level)['max_dx'] for level in range(3)] == [1e-4, 1e-4/2, 1e-4/4]
    assert runner.LAST_WIDTHS == (Q('0.8'), Q('0.73'), Q('0.66'))
    assert runner.CHARGED_EVALUATIONS == 500_000_000_000
    assert runner.CPU_SECONDS == 10000.
    assert runner.ORIGINAL_ALPHA == .00012983127142551443
    with pytest.raises(ValueError):
        runner.level_spec(3)


def test_last_cell_covers_the_reference_radius_with_a_finite_guard():
    budget = Budget(maximum_evaluations=2_000_000, cpu_seconds=30.)
    for level in range(3):
        grid = runner.grid_report(level)
        knots, widths = grid['knots'], grid['widths']
        assert knots[0] == 0
        assert runner.last_cell_cut(knots, runner.REFERENCE_RADIUS_ANCHOR) == runner.INITIAL_WALL_CUT
        assert knots[-1]-knots[-2] == runner.LAST_WIDTHS[level]
        assert len(set(widths[:-1])) == 1
        assert grid['width_max'] > grid['width_min']
        assert grid['uniform_h_family'] is False
        assert grid['remesh'] is grid['adapt'] is grid['activation'] is False
        assert all(Q(1, 4) < cut < Q(3, 4) for cut in grid['reference_cuts'])
        assert grid['cut_guard'] > 0
        basis = FixedCapBasis(knots, budget=budget)
        for radius in (runner.REFERENCE_RADIUS_ANCHOR, runner.REFERENCE_RADIUS_ENDPOINT):
            basis.gram(float(math.atanh(radius)), budget=budget)
    ratios = [runner.LAST_WIDTHS[i+1]/runner.LAST_WIDTHS[i] for i in range(2)]
    assert ratios != [Q(1, 2), Q(1, 2)]


def test_partition_uses_rows_and_fronts_without_ulp_steps():
    pieces, total = runner.time_partition(X, FRONT, BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
                                          max_dx=1e-4, step_budget=2000)
    assert pieces[0][0] == runner.X_ANCHOR and pieces[-1][1] == runner.X_LAST
    assert total <= 2000
    assert min(right-left for left, right, _ in pieces) > 1e-8
    for left, right, count in pieces:
        assert (right-left)/count <= 1e-4
        nodes = runner.step_nodes(left, right, count, 1e-4)
        assert nodes[0] == left and nodes[-1] == right
    bounds = [pieces[0][0]]+[right for _, right, _ in pieces]
    for cut in runner.common_source_cuts(X):
        assert any(abs(float(cut)-bound) <= runner._ulp_gap(float(cut), bound) for bound in bounds)
    shifted, shifted_total = runner.time_partition(
        X, float(X[40])+math.ulp(float(X[40])), BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
        max_dx=1e-4, step_budget=2000)
    exact, exact_total = runner.time_partition(
        X, float(X[40]), BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST, max_dx=1e-4, step_budget=2000)
    assert shifted_total == exact_total
    assert [(a, b) for a, b, _ in shifted] == [(a, b) for a, b, _ in exact]
    midpoint = (float(X[40])+float(X[41]))/2
    split, split_total = runner.time_partition(X, midpoint, BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
                                               max_dx=1e-4, step_budget=2000)
    assert split_total >= exact_total
    assert len(split) > len(exact)
    rows = list(X)
    rows[-1] = runner.X_LAST+math.ulp(runner.X_LAST)
    endpoint, _ = runner.time_partition(rows, FRONT, BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
                                        max_dx=1e-4, step_budget=2000)
    assert endpoint[-1][1] == runner.X_LAST


def test_levels_share_the_nine_source_cuts_and_fit_their_step_budgets():
    cuts = runner.common_source_cuts(X)
    assert len(cuts) == 9
    assert np.array_equal(cuts, X[::32])
    assert float(cuts[0]) == runner.X_ANCHOR and float(cuts[-1]) == runner.X_LAST
    totals = []
    for level in range(3):
        spec = runner.level_spec(level)
        _, total = runner.time_partition(X, FRONT, BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
                                         max_dx=spec['max_dx'], step_budget=spec['step_budget'])
        totals.append(total)
    assert totals[0] < totals[1] < totals[2]
    assert totals[2] <= 8000


def test_over_budget_schedule_is_rejected_before_any_step():
    with pytest.raises(ValueError, match='exceeds admitted steps'):
        runner.time_partition(X, FRONT, BACK, x_left=runner.X_ANCHOR, x_right=runner.X_LAST,
                              max_dx=1e-4, step_budget=3)


def test_quadratures_and_transport_gap_do_not_correct_a_state():
    width = .2
    assert runner.ssprk3_increment(width, 3., 3., 3.) == pytest.approx(.6)
    assert runner.ssprk3_increment(width, 0., width/2, width) == pytest.approx(width*width/2)
    assert runner.ssprk3_increment(1., 0., .25, 1.) == pytest.approx(1/3)
    trap = runner.RunningTrapezoid()
    trap.add(0., dict(physical_temporal_mass_actual=0.))
    trap.add(1., dict(physical_temporal_mass_actual=2.))
    assert trap.signed['physical_temporal_mass_actual'] == pytest.approx(1.)
    assert trap.absolute['physical_temporal_mass_actual'] == pytest.approx(1.)
    comparison = runner.TransportComparison()
    comparison.start(1., 0.)
    comparison.advance(.5, .25, 1.4, .2)
    assert comparison.last['mass_signed_gap'] == pytest.approx(-.1)
    assert comparison.last['mass_absolute_gap'] == pytest.approx(.1)
    assert comparison.last['scalar_signed_gap'] == pytest.approx(-.05)
    assert comparison.last['scalar_absolute_gap'] == pytest.approx(.05)
    assert comparison.last['state_corrected'] is False
    rate = runner.matched_mass_rate(radius=.5, radius_rate=.2, alpha=.1, moving_density=.3,
                                    clock_rate=1.5, metric_A=.8, lapse=0., Pi=.4, Phi=.5)
    transport = 1.5*(.8**2*math.exp(0.)/2)*.5**2*.4*.5
    assert rate == pytest.approx(1.5*.5**2*.2+.1*(.3+transport))
    assert runner.physical_mass(.5, .1, .01) == pytest.approx(.5**3/2+.01*.1)


def test_jordan_radius_uses_the_invariant_scalar_and_flags_a_cosh_mismatch():
    assert runner.jordan_radius is jordan_radius
    radius, canonical_f, alpha = .5, .2, runner.ORIGINAL_ALPHA
    q = math.sqrt(alpha/6)*canonical_f
    scalar_c = radius*canonical_f*math.sinh(q)/q
    observables = runner.frame_observables(radius=radius, scalar_c=scalar_c, canonical_f=canonical_f, alpha=alpha,
                                           normalized_mass=.1, central_T=1.25, einstein_s=.4)
    assert observables['jordan_radius'] == pytest.approx(radius*math.cosh(q))
    assert observables['frame_formula_needs_check'] is False
    assert observables['wall_proper_time_frame'] == 'Einstein'
    assert observables['einstein_wall_proper_time_s'] == .4
    assert observables['central_proper_clock_T'] == 1.25
    assert 'jordan_wall_proper_time_s' not in observables
    assert observables['units'] == 'dimensionless'
    assert observables['physical_mass'] == pytest.approx(radius**3/2+alpha*.1)
    mismatched = runner.frame_observables(radius=radius, scalar_c=scalar_c+1e-3, canonical_f=canonical_f,
                                          alpha=alpha, normalized_mass=.1, central_T=1.25, einstein_s=.4)
    assert mismatched['frame_formula_needs_check'] is True


def test_label_margin_rejects_extrapolation_and_physical_gaps_exclude_weak_balance():
    assert runner.require_label_inside(1., 0., .4, .2, quiet_extension=False) == pytest.approx(.6)
    with pytest.raises(ValueError, match='extrapolate past'):
        runner.require_label_inside(1., 0., 1.1, .2, quiet_extension=True)
    with pytest.raises(ValueError, match='below'):
        runner.require_label_inside(1., 0., .4, -.1, quiet_extension=False)
    assert runner.require_label_inside(1., 0., .4, -.1, quiet_extension=True) == pytest.approx(.6)
    gaps = runner.physical_gaps(dict(physical_temporal_mass_actual=1., physical_temporal_mass_matched=-2.,
                                     conformal_wall_flux=3., weak_balance_defect=99.))
    assert gaps == dict(mass_actual=1., mass_matched=-2., flux=3.)
    assert 'weak_balance_defect' not in gaps


def test_summary_guard_and_vacuum_initial_state_reject_raw_arrays():
    with pytest.raises(ValueError, match='raw array'):
        runner.assert_summary(dict(final=np.zeros(4)))
    with pytest.raises(ValueError, match='raw numeric series'):
        runner.assert_summary(list(range(20)))
    runner.assert_summary(dict(cuts=[float(v) for v in X[::32]], notes=dict(frame='Einstein')))
    budget = Budget(maximum_evaluations=100000, cpu_seconds=30.)
    basis = FixedCapBasis(runner.cap_knots(0), budget=budget)
    state = runner.initial_state(basis, 32)
    assert set(state) == {'departure', 'wall', 'cap'}
    assert all(np.all(value == 0.) for value in state['departure'].values())
    assert all(value == 0. for value in state['wall'].values())
    assert all(np.all(value == 0.) for value in state['cap'].values())
    assert list(runner.eta_patch(32)) == list(runner.COMMON_ETA_OVER_32)
    assert list(runner.eta_patch(128)) == [4*k for k in runner.COMMON_ETA_OVER_32]


def test_existing_output_is_preserved_before_any_calculation(tmp_path, monkeypatch):
    path = tmp_path/'saved.json'
    path.write_bytes(b'preserve this receipt')
    monkeypatch.setattr(runner.sys, 'argv', ['runner', '--level', '1', '--output', str(path)])

    def forbidden():
        raise AssertionError('source discovery must not run')

    monkeypatch.setattr(runner, 'frozen_sources', forbidden)
    with pytest.raises(ValueError, match='preserve the existing result'):
        runner.main()
    assert path.read_bytes() == b'preserve this receipt'


def test_level_and_output_are_required_before_source_discovery(monkeypatch):
    def forbidden():
        raise AssertionError('source discovery must not run')

    monkeypatch.setattr(runner, 'frozen_sources', forbidden)
    monkeypatch.setattr(runner.sys, 'argv', ['runner', '--output', 'results/summaries/new.json'])
    with pytest.raises(SystemExit) as missing:
        runner.main()
    assert missing.value.code == 2
    monkeypatch.setattr(runner.sys, 'argv', ['runner', '--level', '3', '--output', 'results/summaries/new.json'])
    with pytest.raises(SystemExit) as rejected:
        runner.main()
    assert rejected.value.code == 2
