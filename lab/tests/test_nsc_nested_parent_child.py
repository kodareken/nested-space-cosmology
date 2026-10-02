"""Canonical, source, and inherited-functional checks for the nested pair.

Manufactured admissible columns make the reciprocal source-force tests
nonvacuous; the physical T0 sigma2-eigenstate preparation can have zero
conformal source force. No saved scientific trajectory is regenerated.
"""
import numpy as np
import pytest

from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_feedback_action as action
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture(scope="module")
def pair():
    return nested.build_pair(64)


def manufactured(pair):
    base = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    x = pair.grid.xi_g
    angle = 2 * np.pi * x / pair.grid.length
    base.Q = 1.1 + .05 * np.cos(angle)
    base.r = 1.2 + .04 * np.sin(2 * angle)
    base.chi = .03 * np.cos(2 * angle)
    base.p_Q = .07 + .03 * np.sin(angle)
    base.p_r = .06 * np.cos(2 * angle)
    base.p_chi = .04 + .02 * np.sin(2 * angle)
    rng = np.random.default_rng(4701)
    columns, _ = np.linalg.qr(rng.normal(size=(2 * pair.grid.nf, 6)) +
                              1j * rng.normal(size=(2 * pair.grid.nf, 6)))
    base.phi0, base.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    return nested.encode_state(pair, base), base


def test_new_preparation_preserves_child_and_separates_source_from_observer(pair):
    assert pair.parent_interval[0] < pair.child_interval[0] < pair.child_interval[1] < pair.parent_interval[1]
    assert np.array_equal(pair.source_columns[:, 2:4], pair.reference_columns[:, 2:4])
    assert not np.array_equal(pair.source_columns, pair.reference_columns)
    assert pair.source_metadata["layout"] == "separated"
    assert pair.source_metadata["source_support_closures"] == [[0., 1.], [1., 3.], [3., 4.]]
    assert pair.source_metadata["exact_spatial_support"] is False
    assert max(pair.source_metadata["fine_outside_support_power_fraction"]) > 0
    assert np.max(abs(pair.source_columns.conj().T @ pair.source_columns - np.eye(6))) < 1e-12
    phase = np.exp(1j * coupling.CALIBRATION["phase"])
    for left in (0, 2, 4):
        assert np.max(abs(pair.source_columns[:, left + 1] - phase * pair.source_columns[:, left].conj())) < 2e-12
    assert pair.weights.sum() == 3
    # Fixed observer/source arrays cannot be accidentally reset by a driver.
    with pytest.raises(ValueError):
        pair.reference_columns[0, 0] = 0


def test_geometry_frame_is_canonical_complete_and_child_nodally_local(pair):
    matrix = pair.geometry_map
    assert np.isrealobj(matrix)
    assert matrix.shape == (pair.grid.ng, pair.grid.ng)
    assert np.max(abs(matrix.T @ matrix - np.eye(pair.grid.ng))) < 2e-11
    child = matrix[:, pair.geometry_child_indices]
    outside = (pair.grid.xi_g < 1) | (pair.grid.xi_g > 3)
    assert np.max(abs(child[outside])) == 0
    assert pair.geometry_metadata["child_coarse_moment_max"] < 2e-11
    assert pair.geometry_metadata["exact_spatial_support"] is False
    assert np.intersect1d(pair.geometry_parent_indices, pair.geometry_child_indices).size == 0
    assert np.union1d(pair.geometry_parent_indices, pair.geometry_child_indices).size == pair.grid.ng
    state, base = manufactured(pair)
    recovered = nested.reconstruct_state(pair, state)
    for name in nested.STATE_NAMES:
        assert np.max(abs(getattr(recovered, name) - getattr(base, name))) < 5e-12
    rng = np.random.default_rng(44)
    for name, momentum in zip(nested.GEOMETRY_NAMES, nested.MOMENTUM_NAMES):
        direction = rng.normal(size=pair.grid.ng)
        old_form = pair.grid.dx_g * np.dot(getattr(base, momentum), matrix @ direction)
        new_form = np.dot(getattr(state, momentum), direction)
        assert abs(old_form - new_form) < 1e-12


def test_rates_are_full_action_pullback_and_one_step_matches_owner(pair):
    state, base = manufactured(pair)
    result, bundle = nested.rates(pair, state, return_bundle=True)
    original = galerkin.rates(pair.grid, base)
    for name in nested.GEOMETRY_NAMES:
        assert np.max(abs(pair.geometry_map @ getattr(result, name) - getattr(original, name))) < 2e-11
    for name in nested.MOMENTUM_NAMES:
        assert np.max(abs(pair.geometry_map @ getattr(result, name) / pair.grid.dx_g - getattr(original, name))) < 2e-9
    assert np.max(abs(result.phi0 - original.phi0)) < 2e-12
    assert abs(result.fieldwork_power - original.fieldwork_power) < 1e-10
    assert np.array_equal(bundle["lapse_dot"], bundle["lifted_Q"])
    expected = galerkin.rk4_step(pair.grid, base, 2e-5)
    actual = nested.reconstruct_state(pair, nested.rk4_step(pair, state, 2e-5))
    for name in nested.STATE_NAMES:
        assert np.max(abs(getattr(actual, name) - getattr(expected, name))) < 2e-10


def test_canonical_hamiltonian_directional_derivatives_in_both_groups(pair):
    state, _base = manufactured(pair)
    rate = nested.rates(pair, state)
    eps = 1e-6
    for group in (pair.geometry_parent_indices, pair.geometry_child_indices):
        direction = np.zeros(pair.grid.ng)
        direction[group[0]] = 1.
        for coordinate, momentum in zip(nested.GEOMETRY_NAMES, nested.MOMENTUM_NAMES):
            for name, expected in ((coordinate, -np.dot(getattr(rate, momentum), direction)),
                                   (momentum, np.dot(getattr(rate, coordinate), direction))):
                plus, minus = state.copy(), state.copy()
                setattr(plus, name, getattr(plus, name) + eps * direction)
                setattr(minus, name, getattr(minus, name) - eps * direction)
                numerical = (nested.energy(pair, plus) - nested.energy(pair, minus)) / (2 * eps)
                assert abs(numerical - expected) < 3e-6 * max(1., abs(expected))


def test_source_force_and_cross_energy_are_not_omitted(pair):
    state, _base = manufactured(pair)
    source = nested.source_geometry_forces(pair, state)
    forced = nested.rates(pair, state)
    gravity_only = nested.rates(pair, state, include_matter_force=False)
    assert np.max(abs(forced.p_Q - gravity_only.p_Q - source["p_Q_source_rate"])) < 2e-9
    for group in (pair.geometry_parent_indices, pair.geometry_child_indices):
        direction = np.zeros(pair.grid.ng)
        direction[group[0]] = 1.
        eps = 1e-6
        plus, minus = state.copy(), state.copy()
        plus.Q += eps * direction
        minus.Q -= eps * direction
        numerical = (nested.energy_accounting(pair, plus)["field"] -
                     nested.energy_accounting(pair, minus)["field"]) / (2 * eps)
        expected = np.dot(source["Q_energy_gradient"], direction)
        assert abs(expected) > 1e-5
        assert abs(numerical - expected) < 2e-7
    accounting = nested.energy_accounting(pair, state)
    assert abs(accounting["field_closure_error"]) < 2e-10
    assert abs(accounting["total"] - nested.energy(pair, state)) < 2e-10
    cross = sum(abs(value) for key, value in accounting.items() if key.endswith("_cross"))
    assert cross > 1e-3
    # Removing the cross terms is a real mutation, not another normalization.
    diagonal = sum(accounting[f"field_{name}"] for name in ("child", "parent_detail", "ambient"))
    assert abs(diagonal - accounting["field"]) > 1e-4


def test_full_representative_hamiltonian_and_nested_blocks(pair):
    state, _base = manufactured(pair)
    matrix = nested.hamiltonian(pair, state)
    assert matrix.shape == (2 * pair.grid.nf, 2 * pair.grid.nf)
    assert np.max(abs(matrix - matrix.conj().T)) < 2e-11
    columns = np.vstack((state.phi0, state.phi1))
    assert np.max(abs(matrix @ columns - nested.apply_hamiltonian(pair, state, columns))) < 2e-11
    rate = nested.rates(pair, state)
    assert np.max(abs(np.vstack((rate.phi0, rate.phi1)) + 1j * matrix @ columns)) < 2e-11
    blocks = nested.live_blocks(pair, state)
    child = pair.reference_columns[:, pair.child_indices]
    detail = pair.reference_columns[:, [0, 1, 4, 5]]
    assert np.max(abs(blocks["child_parent_link"] - child.conj().T @ matrix @ detail)) < 1e-12
    assert np.linalg.norm(blocks["child_parent_link"]) > 1e-3
    assert np.linalg.norm(blocks["ambient_image"]) > 1e-3


def test_affine_inherited_action_and_clock_map_without_running_parameters():
    jets = dict(L=1.2, L_x=.17, Q=1.1, Q_t=.2, Q_x=-.13,
                r=1.3, r_t=.14, r_x=.21, chi=.09, chi_t=-.16, chi_x=.08,
                beta=.12, beta_x=-.07)
    parameters = dict(A=.04, C_W=-.001, C_F=.01, flux=4.)
    density = action.first_order_density(**jets, **parameters)
    for scale in (.5, 2., 4.):
        local = nested.functional_pullback(jets, scale)
        actual = action.first_order_density(**local, **parameters)
        assert abs(actual - scale ** 2 * density) < 2e-12
        assert local["r"] * local["L"] == pytest.approx(scale * jets["r"] * jets["L"])
        # Wrong momenta or half-density weights would not preserve the action.
        extra = nested.functional_pullback(dict(p_Q=.3, p_r=.4, p_chi=.5,
                                                 phi=np.array([1.+2j])), scale)
        assert extra["p_Q"] == .3
        assert extra["p_r"] == pytest.approx(scale * .4)
        assert abs(extra["phi"][0]) ** 2 == pytest.approx(scale * 5.)


def test_full_metric_restrictions_and_original_clock_locations(pair):
    state, base = manufactured(pair)
    report = nested.metrics(pair, state, nested.rates(pair, state))
    fine = galerkin.prolong_state(pair.grid, base)
    assert report["one_physical_metric"] is True
    assert report["clock_locations"] == [1., 2., 3.]
    assert report["parent_proper_length"] > report["child_proper_length"] > 0
    assert report["parent_annulus_proper_length"] == pytest.approx(
        report["parent_proper_length"] - report["child_proper_length"])
    # Independent analytic product integral for rQ, rather than the same FFT helper.
    k = 2 * np.pi / pair.grid.length
    def primitive(x):
        return 1.32 * x + .06 * np.sin(k * x) / k - .044 * np.cos(2*k*x) / (2*k) - \
               .001 * np.cos(3*k*x) / (3*k) - .001 * np.cos(k*x) / k
    for interval, key in ((pair.parent_interval, "parent_proper_length"),
                          (pair.child_interval, "child_proper_length")):
        assert report[key] == pytest.approx(primitive(interval[1]) - primitive(interval[0]), abs=2e-11)
    expected_clocks = (1.2 + .04 * np.sin(2*k*np.array(pair.clock_locations))) * \
                      (1.1 + .05 * np.cos(k*np.array(pair.clock_locations)))
    assert np.max(abs(np.array(report["clock_rates"]) - expected_clocks)) < 2e-11
    assert fine.r.min() > 0


def test_api_rejects_invalid_dimensions_scale_and_chart(pair):
    with pytest.raises(ValueError, match="nf"):
        nested.build_pair(63)
    with pytest.raises(ValueError, match="nullspace|fit"):
        nested.build_pair(32, child_details=12, source_layout="original")
    with pytest.raises(ValueError, match="source_layout"):
        nested.build_pair(64, source_layout="unknown")
    with pytest.raises(ValueError, match="positive"):
        nested.functional_pullback({"Q": 1.}, 0)
    with pytest.raises(ValueError, match="unsupported"):
        nested.functional_pullback({"unknown": 1.}, 1)
    state, _base = manufactured(pair)
    bad = state.copy()
    bad.Q = np.ones(pair.grid.ng + 1)
    with pytest.raises(ValueError, match="ng-vector"):
        nested.rates(pair, bad)
    bad = state.copy()
    bad.Q = -np.sqrt(pair.grid.ng) * np.eye(pair.grid.ng)[0]
    with pytest.raises(coupling.PositiveChartExit, match="Q"):
        nested.metrics(pair, bad)


@pytest.mark.parametrize("nf", (128, 256))
def test_successor_initial_radius_is_source_derived_at_measured_arithmetic_floor(nf):
    local_pair = nested.build_pair(nf)
    state, report = nested.initial_state(local_pair)
    assert report["converged"]
    assert report["correction_at_arithmetic_floor"]
    assert report["newton_correction_inf"] <= report["radius_correction_floor"]
    assert report["reordered_correction_inf"] <= report["radius_correction_floor"]
    assert report["linear_backward_error"] < 2e-13
    assert np.isfinite(report["jacobian_condition_inf"])
    assert report["bracket_positive"]
    assert report["rho_independent_of_r_max"] == 0
    assert report["continuum_initial_state_certified"] is False
    assert report["arithmetic_indicator_certified"] is False
    assert report["imposed_radius"] is False
    assert max(step["iterations"] for step in report["steps"]) < 24
    assert np.array_equal(state.phi0, local_pair.source_phi0)
    assert np.array_equal(state.phi1, local_pair.source_phi1)
    actual = nested.metrics(local_pair, state)
    assert actual["r_min"] > 0
    assert actual["child_r_proper_mean"] > 1
    assert abs(report["source_current_mean"]) < 2e-12
    # A source change independently changes its constraint-derived radius.
    if nf == 128:
        other_pair = nested.build_pair(nf, occupations=(.85, .85, .5, .5, .15, .15))
        other_state, other_report = nested.initial_state(other_pair)
        assert other_report["converged"]
        assert np.linalg.norm(other_state.r - state.r) > .01
