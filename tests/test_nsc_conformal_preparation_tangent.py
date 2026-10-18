"""Manufactured controls for one paired-preparation tangent. No archive is read."""
import numpy as np
import pytest
from scipy.linalg import expm

from recursive_horizons.nsc_conformal_preparation_tangent import (
    RELATIVE_SOLVE_GATE,
    _factor_once,
    apply_preparation,
    preparation_generator,
    preparation_tangent,
    spinor_tangent,
)
from recursive_horizons.nsc_covariant_operator import SIGMA2
from recursive_horizons.nsc_spherical_coupling import (
    CALIBRATION,
    OCCUPATIONS,
    PERIOD,
    hamilton_constraint,
    shift_constraint,
    source_from_columns,
)
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    _constraint_directional,
    _projected_radius_jacobian,
    blank_state,
    build_grid,
    prolong_columns,
    prolong_geometry,
    prolong_state,
    pull_geometry,
    solve_initial_radius,
)


def _chiral_columns(nf, phase=None):
    if phase is None:
        phase = float(CALIBRATION["phase"])
    xi = np.arange(nf, dtype=float) * (PERIOD / nf)
    phi0 = np.zeros((nf, 6), dtype=complex)
    phi1 = np.zeros((nf, 6), dtype=complex)
    recorded = np.exp(1j * phase)
    for region in range(3):
        mode = 0.5 + region
        carrier = np.exp(1j * 2.0 * np.pi * mode * xi / PERIOD)
        carrier = carrier * (1.0 + 0.35 * np.cos(2.0 * np.pi * xi / PERIOD))
        carrier = carrier / np.linalg.norm(carrier)
        plus0 = carrier / np.sqrt(2.0)
        plus1 = 1j * plus0
        phi0[:, 2 * region] = plus0
        phi1[:, 2 * region] = plus1
        phi0[:, 2 * region + 1] = recorded * np.conjugate(plus0)
        phi1[:, 2 * region + 1] = recorded * np.conjugate(plus1)
    return phi0, phi1


def _grid():
    return build_grid(10, gauge="conformal")


def _blank():
    grid = _grid()
    phi0, phi1 = _chiral_columns(grid.nf)
    return grid, blank_state(grid, phi0, phi1)


def _covariance(phi0, phi1):
    packed = np.vstack((phi0, phi1))
    return (packed * OCCUPATIONS) @ packed.conj().T


def _gram(phi0, phi1):
    return phi0.conj().T @ phi0 + phi1.conj().T @ phi1


@pytest.fixture(scope="module")
def solved_tangent():
    grid, guess = _blank()
    solved, info = solve_initial_radius(grid, guess)
    if not info["converged"]:
        raise AssertionError(info["blocker"])
    radius = solved.r.copy()
    phi0 = solved.phi0.copy()
    tangent = preparation_tangent(grid, solved)
    if not np.array_equal(radius, solved.r) or not np.array_equal(phi0, solved.phi0):
        raise AssertionError("preparation tangent mutated the supplied state")
    return grid, solved, tangent


def test_component_formula_matches_sigma2_exponential():
    generator = np.array([0.3, -0.7, 1.0])
    epsilon = 0.2
    phi0 = np.array([[1.0 + 0.2j, -0.4], [0.5j, 0.2], [-0.3 + 0.1j, 0.7]])
    phi1 = np.array([[0.3j, 1.0], [-0.2 + 0.1j, -0.5j], [0.4, 0.25 - 0.3j]])
    expected0 = np.empty_like(phi0)
    expected1 = np.empty_like(phi1)
    for node, value in enumerate(generator):
        unitary = expm(1j * epsilon * value * SIGMA2)
        rotated = unitary @ np.vstack((phi0[node], phi1[node]))
        expected0[node] = rotated[0]
        expected1[node] = rotated[1]
    got0, got1 = apply_preparation(phi0, phi1, generator, epsilon)
    assert got0 == pytest.approx(expected0, abs=1e-12)
    assert got1 == pytest.approx(expected1, abs=1e-12)
    identity0, identity1 = apply_preparation(phi0, phi1, generator, 0.0)
    assert identity0 == pytest.approx(phi0, abs=1e-15)
    assert identity1 == pytest.approx(phi1, abs=1e-15)
    step = 1e-7
    plus0, plus1 = apply_preparation(phi0, phi1, generator, step)
    minus0, minus1 = apply_preparation(phi0, phi1, generator, -step)
    tangent0, tangent1 = spinor_tangent(phi0, phi1, generator)
    assert (plus0 - minus0) / (2 * step) == pytest.approx(tangent0, abs=1e-8)
    assert (plus1 - minus1) / (2 * step) == pytest.approx(tangent1, abs=1e-8)


def test_prolonged_family_preserves_conjugacy_gram_and_cancels_current():
    grid, state = _blank()
    generator = preparation_generator(grid.xi_f, grid.length)
    original_gram = _gram(state.phi0, state.phi1)
    recorded = np.exp(1j * float(CALIBRATION["phase"]))
    tangent0, tangent1 = spinor_tangent(state.phi0, state.phi1, generator)
    for epsilon in (0.0, 0.2, -0.5, 1.3):
        rotated0, rotated1 = apply_preparation(state.phi0, state.phi1, generator, epsilon)
        assert _gram(rotated0, rotated1) == pytest.approx(original_gram, abs=1e-10)
        for region in range(3):
            plus0, plus1 = rotated0[:, 2 * region], rotated1[:, 2 * region]
            minus0, minus1 = rotated0[:, 2 * region + 1], rotated1[:, 2 * region + 1]
            assert plus1 == pytest.approx(1j * plus0, abs=1e-10)
            assert minus0 == pytest.approx(recorded * np.conjugate(plus0), abs=1e-10)
            assert minus1 == pytest.approx(recorded * np.conjugate(plus1), abs=1e-10)
        trial = state.copy()
        trial.phi0, trial.phi1 = rotated0, rotated1
        fine = prolong_state(grid, trial)
        source = source_from_columns(grid.fine, fine)
        assert np.max(np.abs(source["S1"])) < 1e-8
        assert np.max(np.abs(source["Pmom"])) < 1e-8
        assert np.max(np.abs(source["force_beta"])) < 1e-8
    assert tangent0 == pytest.approx(generator[:, None] * state.phi1)
    assert tangent1 == pytest.approx(-generator[:, None] * state.phi0)
    # The highest AP mode times this cosine aliases on nf. Prolonging that
    # product is not the fine-grid product of the two prolongations.
    mode = grid.nf / 2 - 0.5
    carrier = np.exp(1j * 2.0 * np.pi * mode * grid.xi_f / grid.length)
    aliased = prolong_columns(grid, (generator * carrier)[:, None])[:, 0]
    pointwise = preparation_generator(grid.xi_q, grid.length) * prolong_columns(
        grid, carrier[:, None],
    )[:, 0]
    assert np.max(np.abs(aliased - pointwise)) > 1e-3


def test_constant_generator_leaves_covariance_and_source_unchanged():
    grid, state = _blank()
    generator = np.ones(grid.nf)
    rotated0, rotated1 = apply_preparation(state.phi0, state.phi1, generator, 0.37)
    assert _covariance(rotated0, rotated1) == pytest.approx(_covariance(state.phi0, state.phi1), abs=1e-10)
    original = source_from_columns(grid.fine, prolong_state(grid, state))
    trial = state.copy()
    trial.phi0, trial.phi1 = rotated0, rotated1
    rotated = source_from_columns(grid.fine, prolong_state(grid, trial))
    for name in ("K", "S1", "Pmom", "force_L", "force_Q", "force_beta"):
        assert rotated[name] == pytest.approx(original[name], abs=1e-10)


def test_nonconstant_source_derivative_matches_centered_difference():
    grid, state = _blank()
    generator = preparation_generator(grid.xi_f, grid.length)
    fine = prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    tangent0, tangent1 = spinor_tangent(state.phi0, state.phi1, generator)
    zero = np.zeros(grid.nq)
    _h, _d, h_dot, d_dot = _constraint_directional(
        grid.fine,
        fine,
        (
            zero, zero, zero, zero, zero, zero,
            prolong_columns(grid, tangent0),
            prolong_columns(grid, tangent1),
        ),
        source,
    )

    def forces(epsilon):
        rotated0, rotated1 = apply_preparation(state.phi0, state.phi1, generator, epsilon)
        trial = state.copy()
        trial.phi0, trial.phi1 = rotated0, rotated1
        return source_from_columns(grid.fine, prolong_state(grid, trial))

    step = 1e-6
    forward, backward = forces(step), forces(-step)
    force_l = (forward["force_L"] - backward["force_L"]) / (2 * step)
    force_beta = (forward["force_beta"] - backward["force_beta"]) / (2 * step)
    assert np.max(np.abs(force_l)) > 1e-4
    assert force_l / grid.dx_q == pytest.approx(np.real(h_dot) / fine.Q, abs=1e-6)
    assert force_beta / grid.dx_q == pytest.approx(np.real(d_dot), abs=1e-7)
    assert np.max(np.abs(force_beta)) < 1e-7


def _projected_density(grid, state):
    fine = prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    density = hamilton_constraint(grid.fine, fine) + source["force_L"] / grid.dx_q
    momentum = shift_constraint(grid.fine, fine) + source["force_beta"] / grid.dx_q
    return pull_geometry(grid, density), density, momentum


def test_radius_jacobian_matches_independent_constraint_difference(solved_tangent):
    grid, solved, _tangent = solved_tangent
    direction = np.cos(2.0 * np.pi * grid.xi_g / grid.length)
    step = 1e-6

    def projected(radius):
        trial = solved.copy()
        trial.r = radius
        values, _density, _momentum = _projected_density(grid, trial)
        return values

    difference = (
        projected(solved.r + step * direction) - projected(solved.r - step * direction)
    ) / (2 * step)
    jacobian = _projected_radius_jacobian(grid, solved.r)
    assert difference == pytest.approx(jacobian @ direction, abs=1e-7, rel=1e-6)


def test_tangent_cancels_projected_constraint_and_reports_nonzero_held_out(solved_tangent):
    grid, solved, tangent = solved_tangent
    generator = preparation_generator(grid.xi_f, grid.length)
    assert tangent.phi0_tangent == pytest.approx(generator[:, None] * solved.phi1)
    assert tangent.phi1_tangent == pytest.approx(-generator[:, None] * solved.phi0)
    assert tangent.report["geometry_variations_per_epsilon"] == {
        "Q": 0.0, "chi": 0.0, "p_Q": 0.0, "p_r": 0.0, "p_chi": 0.0,
    }
    assert tangent.report["solve"]["certified_invertibility"] is False
    assert tangent.report["solve"]["relative_residual"] < RELATIVE_SOLVE_GATE
    assert tangent.report["source_tangent"]["projected_source_max_abs"] > 1e-6
    assert np.max(np.abs(tangent.delta_r)) > 1e-8
    # The cancelling projected pieces are O(1). A 1e-6 step loses the
    # unresolved remainder under that cancellation; 1e-4 keeps it.
    step = 1e-4

    def sampled(epsilon):
        rotated0, rotated1 = apply_preparation(solved.phi0, solved.phi1, generator, epsilon)
        trial = solved.copy()
        trial.r = solved.r + epsilon * tangent.delta_r
        trial.phi0, trial.phi1 = rotated0, rotated1
        return _projected_density(grid, trial)

    projected_plus, density_plus, momentum_plus = sampled(step)
    projected_minus, density_minus, momentum_minus = sampled(-step)
    projected_difference = (projected_plus - projected_minus) / (2 * step)
    density_difference = (density_plus - density_minus) / (2 * step)
    momentum_difference = (momentum_plus - momentum_minus) / (2 * step)
    assert projected_difference == pytest.approx(0.0, abs=1e-6)
    assert momentum_difference == pytest.approx(0.0, abs=1e-6)
    held = density_difference - prolong_geometry(grid, pull_geometry(grid, density_difference))
    held_max = float(np.max(np.abs(held)))
    assert held_max > 1e-6
    reported = tangent.report["directional_check"]["held_out_max_abs"]
    assert reported == pytest.approx(held_max, rel=1e-2, abs=1e-7)
    assert reported > 1e-6
    observer = tangent.report["observer"]
    assert observer["x"] == 1.0
    assert observer["delta_N_over_N_per_epsilon"] == pytest.approx(
        observer["delta_r_per_epsilon"] / observer["original_radius"]
    )
    assert observer["light_cone_evolution"] is False


def test_invalid_premises_stop_before_factorization(monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("factorization ran before a premise was accepted")

    monkeypatch.setattr(
        "recursive_horizons.nsc_conformal_preparation_tangent._factor_once", forbidden,
    )
    grid, state = _blank()
    prescribed = build_grid(10)
    with pytest.raises(ValueError, match="conformal"):
        preparation_tangent(prescribed, state)

    short = state.copy()
    short.phi0 = state.phi0[:, :5]
    with pytest.raises(ValueError, match="column arrays"):
        preparation_tangent(grid, short)

    shifted = state.copy()
    shifted.Q = state.Q + 1e-6
    with pytest.raises(ValueError, match="b0/a0"):
        preparation_tangent(grid, shifted)

    chiral = state.copy()
    chiral.chi = np.full(grid.ng, 1e-6)
    with pytest.raises(ValueError, match="chi"):
        preparation_tangent(grid, chiral)

    moving = state.copy()
    moving.p_r = np.full(grid.ng, 1e-6)
    with pytest.raises(ValueError, match="p_r"):
        preparation_tangent(grid, moving)

    unpaired = state.copy()
    unpaired.phi1 = state.phi1.copy()
    unpaired.phi1[:, 0] *= -1
    with pytest.raises(ValueError, match="phi1 = i\\*phi0"):
        preparation_tangent(grid, unpaired)

    reweighted = build_grid(10, gauge="conformal")
    reweighted.fine.occupations = np.array([0.2, 0.9, 0.5, 0.5, 0.25, 0.25])
    with pytest.raises(ValueError, match="occupations"):
        preparation_tangent(reweighted, state)

    with pytest.raises(ValueError, match="finite"):
        apply_preparation(state.phi0, state.phi1, np.ones(grid.nf), np.nan)
    with pytest.raises(ValueError, match="every fermion node"):
        apply_preparation(state.phi0, state.phi1, np.ones(grid.nf - 1), 0.1)
    with pytest.raises(ValueError, match="positive period"):
        preparation_generator(grid.xi_f, 0.0)
    for period in (np.nan, np.inf):
        with pytest.raises(ValueError, match="positive period"):
            preparation_generator(grid.xi_f, period)
    with pytest.raises(ValueError, match="singular"):
        _factor_once(np.zeros((3, 3)), np.ones(3))


def test_factor_residual_is_not_an_invertibility_certificate():
    matrix = np.array([[2.0, 1.0], [0.5, 3.0]])
    rhs = np.array([1.0, -2.0])
    solution, info = _factor_once(matrix, rhs)
    assert matrix @ solution == pytest.approx(rhs, abs=1e-12)
    assert info["certified_invertibility"] is False
    assert info["factorizations"] == 1
    assert info["rcond"] > 1e-8


def test_budget_exhaustion_after_matrix_assembly_prevents_factorization(solved_tangent, monkeypatch):
    grid, state, _ = solved_tangent
    exhausted = False

    def assembled(*_args):
        nonlocal exhausted
        exhausted = True
        return np.eye(grid.ng)

    def checkpoint():
        if exhausted:
            raise RuntimeError("manufactured CPU budget exhausted")

    def forbidden(*_args):
        raise AssertionError("factorization ran after the budget was exhausted")

    prefix = "recursive_horizons.nsc_conformal_preparation_tangent."
    monkeypatch.setattr(prefix + "_projected_radius_jacobian", assembled)
    monkeypatch.setattr(prefix + "_factor_once", forbidden)
    with pytest.raises(RuntimeError, match="CPU budget"):
        preparation_tangent(grid, state, checkpoint=checkpoint)
