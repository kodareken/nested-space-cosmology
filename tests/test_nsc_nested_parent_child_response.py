"""Independent nested trace, Schur, and same-source response controls."""
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from recursive_horizons import nsc_nested_parent_child_response as response


WEIGHTS = np.array([.75, .75, .5, .5, .25, .25])


def problem():
    """Different observer/source and noncommuting live generator; no record run."""
    dimension = 10
    reference = np.eye(dimension, dtype=complex)[:, :6]
    reference *= np.exp(1j * np.arange(6) * .13)
    hierarchy = response.nested_frame(reference)
    rng = np.random.default_rng(49)
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
    h0 = (raw + raw.conj().T) / (2 * dimension)
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
    derivative = (raw + raw.conj().T) / (2 * dimension)
    columns = expm(-.7j * h0) @ reference

    def hamiltonian(t):
        return h0 + (.15 + .2 * np.sin(1.7 * t)) * derivative

    return hierarchy, columns, hamiltonian, derivative


def covariance(columns):
    return (columns * WEIGHTS) @ columns.conj().T


def test_fixed_nested_observers_preserve_phases_and_source_is_independent():
    hierarchy, columns, _, _ = problem()
    reference = hierarchy.parent
    np.testing.assert_array_equal(hierarchy.child, reference[:, [2, 3]])
    np.testing.assert_array_equal(hierarchy.detail, reference[:, [0, 1, 4, 5]])
    np.testing.assert_allclose(hierarchy.frame.conj().T @ hierarchy.frame, np.eye(10), atol=1e-14)
    child_projector = hierarchy.child @ hierarchy.child.conj().T
    parent_projector = hierarchy.parent @ hierarchy.parent.conj().T
    np.testing.assert_allclose(parent_projector @ child_projector, child_projector, atol=1e-14)
    assert np.linalg.norm(hierarchy.ambient.conj().T @ columns) > .1
    assert not np.array_equal(hierarchy.child, columns[:, [2, 3]])
    with pytest.raises(ValueError, match="orthonormal"):
        bad = reference.copy()
        bad[:, 0] += .2 * bad[:, 2]
        response.nested_frame(bad)
    with pytest.raises(ValueError, match="distinct"):
        response.nested_frame(reference, child_indices=(2, 2))


def test_full_covariance_energy_and_reciprocal_currents_match_ambient_commutator():
    hierarchy, columns, operator, _ = problem()
    h = operator(.2)
    c = covariance(columns)
    rate = -1j * (h @ c - c @ h)
    report = response.block_accounting(h, columns, WEIGHTS, hierarchy, multiplicity=4)
    np.testing.assert_allclose(report["source_energy"], 4 * np.trace(c @ h).real, atol=2e-14)
    assert report["energy_closure_gap"] < 2e-14
    assert report["current_closure_gap"] < 2e-14
    np.testing.assert_allclose(report["pair_current_into_row"], -report["pair_current_into_row"].T, atol=1e-14)
    for i, basis in enumerate((hierarchy.child, hierarchy.detail, hierarchy.ambient)):
        np.testing.assert_allclose(report["population"][i], np.trace(basis.conj().T @ c @ basis).real, atol=1e-14)
        np.testing.assert_allclose(report["population_derivative"][i], np.trace(basis.conj().T @ rate @ basis).real, atol=1e-14)
    assert abs(report["population_derivative"].sum()) < 1e-14
    assembled = (report["parent_internal_energy"] + report["ambient_internal_energy"]
                 + report["parent_ambient_interaction_energy"])
    assert assembled == pytest.approx(report["source_energy"], abs=2e-14)
    # The diagonal-only source is a real mutation: it discards actual cross work.
    assert abs(report["source_energy"] - sum(report["diagonal_energy"])) > 1e-4
    assert report["child_ambient_covariance_norm"] > .01


def test_variational_source_contains_link_derivatives_not_only_onsite_force():
    hierarchy, columns, operator, derivative = problem()
    report = response.force_accounting({"parent_coordinate": derivative, "child_coordinate": 2 * derivative},
                                      columns, WEIGHTS, hierarchy, multiplicity=4)
    h, step = operator(.2), 1e-5
    c = covariance(columns)
    finite_difference = 4 * (np.trace(c @ (h + step * derivative)).real - np.trace(c @ (h - step * derivative)).real) / (2 * step)
    assert report["parent_coordinate"]["source_energy"] == pytest.approx(finite_difference, abs=2e-11)
    assert report["child_coordinate"]["source_energy"] == pytest.approx(2 * finite_difference, abs=4e-11)
    assert report["parent_coordinate"]["energy_closure_gap"] < 1e-14
    assert abs(report["parent_coordinate"]["source_energy"] - sum(report["parent_coordinate"]["diagonal_energy"])) > 1e-4
    with pytest.raises(ValueError, match="Hermitian"):
        response.force_accounting({"wrong_link": derivative + .1j * np.eye(10)}, columns, WEIGHTS, hierarchy)


def test_nested_schur_retains_ambient_and_agrees_with_direct_full_resolvent():
    hierarchy, columns, operator, _ = problem()
    z = .3 + .7j
    h = operator(.2)
    report = response.nested_schur(h, hierarchy, z)
    ambient_direct = np.linalg.solve(z * np.eye(10) - h, hierarchy.child)
    direct = hierarchy.child.conj().T @ ambient_direct
    np.testing.assert_allclose(report["direct_child_response"], direct, atol=2e-14)
    assert report["nested_direct_gap"] < 2e-14
    assert report["joined_direct_gap"] < 2e-14
    assert report["nested_joined_inverse_gap"] < 2e-14
    assert report["ambient_omission_response_gap"] > 1e-3
    # Mutation: the independent ambient-only omission has a different response.
    blocks = response.nested_blocks(h, hierarchy)
    assert np.linalg.norm(blocks["child_ambient"]) > .01
    with pytest.raises(ValueError, match="upper half-plane"):
        response.nested_schur(h, hierarchy, z=1.)


@pytest.fixture(scope="module")
def reduced():
    hierarchy, columns, operator, _ = problem()
    times = np.linspace(.2, .3, 21)

    def rhs(t, flat):
        return (-1j * operator(t) @ flat.reshape(columns.shape)).ravel()

    solution = solve_ivp(rhs, (times[0], times[-1]), columns.ravel(), t_eval=times,
                         method="DOP853", rtol=1e-12, atol=1e-14)
    assert solution.success
    actual = np.stack([solution.y[:, i].reshape(columns.shape) for i in range(len(times))])
    # solve_ivp's first sample is the identical supplied actual segment state.
    np.testing.assert_array_equal(actual[0], columns)
    result = response.retained_response(operator, times, columns, WEIGHTS, hierarchy,
                                        autonomous_columns=actual, cpu_limit_s=10)
    return hierarchy, columns, operator, times, actual, result


def test_same_source_full_reduced_controls_and_nested_marginal(reduced):
    hierarchy, columns, operator, times, actual, result = reduced
    child_actual = np.einsum("ij,tjk->tik", hierarchy.child.conj().T, actual)
    covariance_actual = np.einsum("tak,k,tbk->tab", child_actual, WEIGHTS, child_actual.conj())
    np.testing.assert_allclose(result["child"]["full_covariance"], covariance_actual, atol=2e-9)
    assert result["child"]["reduction_error"]["within_one_percent"]
    assert result["child"]["full_midpoint_indicator"]["within_one_percent"]
    assert result["child"]["conditional_versus_autonomous"]["within_one_percent"]
    assert result["parent"]["reduction_error"]["within_one_percent"]
    assert result["nested_full_covariance_gap"] < 1e-14
    assert result["initial_child_detail_cross_norm"] > .01
    assert result["initial_child_ambient_cross_norm"] > .01
    for name, control in result["child"]["controls"].items():
        assert control["occupation_movement"] > 1e-4
        assert control["own_reference_error"]["within_one_percent"]
        assert control["reference_midpoint_indicator"]["within_one_percent"]
        child = result["child"]
        discrepancy = ((child[name] - child["occupation"])
                       - (child[name + "_full"] - child["occupation_full"]))
        assert control["own_reference_error"]["movement"] == pytest.approx(np.max(abs(discrepancy)), abs=1e-15)
        midpoint = ((child[name + "_full"] - child["occupation_full"])
                    - (child[name + "_full_refined"] - child["occupation_full_refined"]))
        assert control["reference_midpoint_indicator"]["movement"] == pytest.approx(np.max(abs(midpoint)), abs=1e-15)
    assert result["controls_additive"] is False
    assert result["all_source_components_retained"] is True
    assert result["observer_preserved"] is True
    for allocation in result["allocation"].values():
        assert allocation["time_indexed_exterior_propagator_bytes"] == 0
        assert allocation["propagator_storage"] == "none"
    assert result["allocation"]["child"]["history_shape"] == [21, 8, 18]
    assert result["allocation"]["parent"]["history_shape"] == [21, 4, 6]


def test_response_domain_mutations_rejected_before_evolution(reduced):
    hierarchy, columns, operator, times, actual, _ = reduced
    changed = actual.copy()
    changed[0, 0, 0] += .001
    with pytest.raises(ValueError, match="actual autonomous segment"):
        response.retained_response(operator, times, columns, WEIGHTS, hierarchy, autonomous_columns=changed)
    changed = times.copy()
    changed[1] += .001
    with pytest.raises(ValueError, match="uniform"):
        response.retained_response(operator, changed, columns, WEIGHTS, hierarchy)
    with pytest.raises(ValueError, match="more midpoint substeps"):
        response.retained_response(operator, times, columns, WEIGHTS, hierarchy, reference_substeps=2)
    with pytest.raises(ValueError, match="admissible"):
        response.block_accounting(operator(.2), 2 * columns, WEIGHTS, hierarchy)


def test_allocation_plan_scales_with_source_columns_not_dense_time_propagators():
    hierarchy, _, _, _ = problem()
    plan = response.response_allocation_plan(hierarchy, 21)
    assert plan["child_history_bytes"] == 21 * 8 * 18 * 16
    assert plan["parent_history_bytes"] == 21 * 4 * 6 * 16
    assert plan["dense_propagator_history_bytes"] == 0
    assert plan["cpu_forecast_available"] is False
