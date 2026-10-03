"""Finite in-memory rank, local induced norm, and sealed-RHS checks."""
from dataclasses import replace
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_parent_step_control as step
from recursive_horizons import nsc_discovery_backend as backend


@pytest.fixture(autouse=True)
def one_thread():
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        yield


def fixture(rank, carrier="dense"):
    leading = step.leading
    grid = leading.galerkin.build_grid(32, gauge="conformal")
    # The immutable FFT factory retains its six-weight metadata contract;
    # resolve carriers before installing the actual source. No dummy columns.
    grid = backend.resolve_grid(grid, carrier)
    weights = np.linspace(.08, .12, rank)
    grid.fine = replace(grid.fine, occupations=weights.copy())
    pair = SimpleNamespace(grid=grid, geometry_map=np.eye(grid.ng), weights=weights)
    rng = np.random.default_rng(1200 + rank)
    columns = np.zeros((2 * grid.nf, rank), complex)
    if rank:
        columns = np.linalg.qr(rng.normal(size=columns.shape) + 1j * rng.normal(size=columns.shape))[0]
        columns = columns @ (np.eye(rank) + .17 * np.triu(np.ones((rank, rank)), 1))
    x = 2 * np.pi * grid.xi_g / grid.length
    nodal = leading.State(.82 + .05 * np.cos(x), 1.8 + .04 * np.cos(2 * x),
                          .03 + .02 * np.sin(x), .06 * np.cos(2 * x),
                          columns[:grid.nf], columns[grid.nf:])
    return pair, leading.encode(pair, nodal)


def real_vector(state):
    return np.r_[[getattr(state, key)[0] for key in step.leading.REAL_FIELDS],
                 state.phi0[0].real, state.phi0[0].imag,
                 state.phi1[0].real, state.phi1[0].imag]


@pytest.mark.parametrize("rank", [0, 1, 2, 6])
def test_all_component_row_norm_against_unscaled_unit_column_jacobian(rank):
    pair, state = fixture(rank)
    _, bundle = step.leading.rates(pair, state, return_bundle=True)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    rows, scales, *_ = step.scaled_local_reaction(pair.grid, fine, system)
    point = 13
    one = step.leading.State(*(getattr(fine, name)[point:point + 1].copy()
                               for name in step.leading.FIELDS))
    local = step.local_system(system)
    local.length_density, local.shift = one.Q, np.zeros(1)
    dimension = 4 + 4 * rank
    jacobian = np.zeros((dimension, dimension))
    for column in range(dimension):
        tangent = step.leading.State(*(np.zeros_like(getattr(one, name)) for name in step.leading.FIELDS))
        if column < 4:
            getattr(tangent, step.leading.REAL_FIELDS[column])[:] = 1
        else:
            group, index = divmod(column - 4, rank)
            field = tangent.phi0 if group < 2 else tangent.phi1
            field[0, index] = 1 if group % 2 == 0 else 1j
        jacobian[:, column] = real_vector(step.leading.fine_jvp(local, one, tangent))
    physical = np.r_[scales[point], np.full(4 * rank, np.sqrt(pair.grid.dx_q))]
    weighted = jacobian * physical[None, :] / physical[:, None]
    np.testing.assert_allclose(rows[point], np.sum(abs(weighted), axis=1), rtol=2e-13, atol=2e-12)
    if rank > 1:
        gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
        assert np.max(abs(gram - np.eye(rank))) > .1
        assert np.linalg.eigvalsh(np.sqrt(pair.weights[:, None] * pair.weights[None, :]) * gram).max() < 1


@pytest.mark.parametrize("rank", [0, 1, 2, 6])
@pytest.mark.parametrize("carrier", ["dense", "fft"])
def test_admission_preserves_state_weights_and_rhs_and_uses_projected_scale_rate(rank, carrier):
    pair, state = fixture(rank, carrier)
    before = state.copy()
    functions = step.leading.rates, step.leading.fine_jvp, step.leading.step_restriction
    rate = step.leading.rates(pair, state)
    record = step.scaled_step_admission(pair, state, .001)
    after = step.leading.rates(pair, state)
    for key in step.leading.FIELDS:
        np.testing.assert_array_equal(getattr(state, key), getattr(before, key))
        np.testing.assert_array_equal(getattr(rate, key), getattr(after, key))
    assert functions == (step.leading.rates, step.leading.fine_jvp, step.leading.step_restriction)
    assert record["source_rank"] == rank and record["local_real_variables"] == 4 + 4 * rank
    assert 0 < record["dt"] <= .001 and not record["physical_growth_removed"]
    assert record["reaction_majorant"] == record["frozen_scaled_reaction_norm"] + record["absolute_moving_scale_rate"]
    fine = step.leading.fine_state(pair, state)
    qdot = step.leading.galerkin.prolong_geometry(pair.grid, pair.geometry_map @ rate.Q)
    rdot = step.leading.galerkin.prolong_geometry(pair.grid, pair.geometry_map @ rate.r)
    expected = np.max(abs(np.column_stack((qdot / fine.Q, rdot / fine.r,
                                          2 * rdot / fine.r - qdot / fine.Q, rdot / fine.r))))
    assert record["absolute_moving_scale_rate"] == expected
    np.testing.assert_array_equal(pair.grid.fine.occupations, pair.weights)
    assert step.step_restriction(pair, state, .001) == (record["dt"], record)
    assert step.scaled_step_restriction(pair, state, .001) == (record["dt"], record)
    assert step.stable_timestep(pair, state, .001) == record["dt"]
    assert record["column_rank"] == rank and not record["raw_coordinate_norm_used"]
    frozen = step.scaled_step_admission(pair, state, .001, control_mode="frozen_geometry")
    assert frozen["absolute_moving_scale_rate"] == 0


def test_rank_two_local_jvp_is_directional_derivative_of_owned_source_rate():
    pair, state = fixture(2)
    fine = step.leading.fine_state(pair, state)
    point = 17
    one = step.leading.State(*(getattr(fine, key)[point:point + 1].copy() for key in step.leading.FIELDS))
    local = step.local_system(step.leading.active_system(pair.grid, fine))
    local.length_density, local.shift = one.Q, np.zeros(1)
    rng = np.random.default_rng(741)
    tangent = step.leading.State(*(rng.normal(size=getattr(one, key).shape)
        + (1j * rng.normal(size=getattr(one, key).shape) if key.startswith("phi") else 0)
        for key in step.leading.FIELDS))
    analytic = step.leading.fine_jvp(local, one, tangent)
    epsilon = 2e-6
    images = []
    for sign in (1, -1):
        varied = step.leading.combine(one, tangent, sign * epsilon)
        active = step.local_system(local)
        active.length_density = varied.Q
        images.append(step.leading.fine_rates(active, varied))
    finite = step.leading.State(*((plus - minus) / (2 * epsilon)
                                  for plus, minus in zip(*images)))
    np.testing.assert_allclose(real_vector(analytic), real_vector(finite), rtol=2e-8, atol=2e-9)


def test_rank_six_reproduces_unchanged_saved_legacy_fixture():
    legacy = step.legacy
    sealed = [Path(step.leading.__file__), Path(legacy.__file__), Path(backend.__file__)]
    hashes = [hashlib.sha256(path.read_bytes()).hexdigest() for path in sealed]
    record, arrays = step.leading.episode.load_checkpoint(legacy.SOURCE, "nf128_coherent_dt0.0005")
    pair = backend.make_fft_pair(step.leading.episode.pair_from_arrays(arrays, record))
    state = step.leading.state_from_arrays(arrays)
    old_dt, old = legacy.step_restriction(pair, state, record["step_cap"])
    new = step.scaled_step_admission(pair, state, record["step_cap"])
    for key, value in old.items():
        if isinstance(value, (int, float)):
            assert new[key] == value
    assert new["dt"] == old_dt
    assert hashes == [hashlib.sha256(path.read_bytes()).hexdigest() for path in sealed]


def test_invalid_rank_weight_spacing_and_cap_rejected():
    pair, state = fixture(2)
    _, bundle = step.leading.rates(pair, state, return_bundle=True)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    for weights in (np.ones(3), np.array([-.1, .2]), np.array([np.nan, .2])):
        with pytest.raises(ValueError, match="weight"):
            step.scaled_local_reaction(pair.grid, fine, replace(system, occupations=weights))
    with pytest.raises(ValueError, match="spacing"):
        step.scaled_local_reaction(pair.grid, fine, replace(system, dx=2 * system.dx))
    for cap in (0, -.1, np.inf, np.nan):
        with pytest.raises(ValueError, match="cap"):
            step.scaled_step_admission(pair, state, cap)
