"""Focused finite static identities and creation-only control checks."""
from pathlib import Path
import tempfile
import json

import numpy as np
import pytest
from scipy.linalg import eigh

from recursive_horizons import nsc_discovery_stationary as stationary
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture
def tmp_path():
    with tempfile.TemporaryDirectory(prefix="nsc-stationary-test-", dir="/tmp") as directory:
        yield Path(directory)


@pytest.fixture(scope="module")
def preview():
    return stationary.run_stationary()


def test_radial_seed_monotonicity_and_positive_zero():
    seed = stationary.radial_seed(points=31)
    scales = [0., seed["scale"] / 2, seed["scale"], 2 * seed["scale"], 1.]
    levels = [eigh(stationary.radial_operator(scale=s)[0], eigvals_only=True)[0] for s in scales]
    assert np.all(np.diff(levels) > 0)
    assert levels[0] < 0 < levels[-1]
    assert abs(levels[2]) < 1e-12
    assert np.min(seed["Q"]) > 0 and np.min(seed["r_shape"]) > 0
    assert np.mean(seed["r_shape"]) == pytest.approx(1.)
    assert seed["radial_residual_max"] < 1e-11
    assert not seed["full_stationary_balance_claimed"]


def test_static_integral_identities_and_constant_q_obstruction():
    seed = stationary.radial_seed(points=63)
    assert seed["integral_Q_square"] == pytest.approx(seed["integral_3_logr_x_square"], abs=1e-11)
    assert seed["integral_r_x_logQ_x"] == pytest.approx(seed["minus_integral_Q_square_r"], abs=1e-11)
    assert seed["integral_r_x_logQ_x"] < 0
    # Constant Q and positive r cannot integrate 3r_xx=Q^2 r to zero.
    r = seed["r_shape"]
    D = stationary.periodic_odd_derivative(len(r), 8.).real
    assert abs(np.sum(3 * (D @ D @ r))) < 1e-10
    assert np.sum(.25**2 * r) > 0


def test_finite_car_commutant_current_and_eigenphase_invariance():
    pair = stationary.make_pair()
    x, _ = stationary.initial_unknown(pair)
    state, spec, _, _, bundle, constraints = stationary.evaluate(pair, x)
    C = spec["covariance"]
    eig = np.linalg.eigvalsh(C)
    assert eig[0] >= -1e-12 and eig[-1] <= 1 + 1e-12
    np.testing.assert_allclose(eig[-6:], np.sort(stationary.WEIGHTS), atol=1e-12)
    assert spec["commutator_gap"] < 1e-11
    assert spec["retained_eigen_residual"] < 1e-11
    assert constraints["current_max"] < 1e-11
    altered = state.copy()
    phase = np.exp(1j * np.arange(6))
    altered.phi0 *= phase
    altered.phi1 *= phase
    _, varied_bundle = galerkin.compose_fine_hamiltonian(pair.grid, altered)
    for name in ("force_L", "force_Q", "force_beta"):
        np.testing.assert_allclose(varied_bundle["source"][name], bundle["source"][name], atol=1e-12)
    assert np.min(spec["eigenvalues"]) > 0


def test_residual_is_original_discrete_rate_and_constraints():
    pair = stationary.make_pair()
    x, _ = stationary.initial_unknown(pair)
    state, _, raw, _, bundle, constraints = stationary.evaluate(pair, x)
    original = galerkin.rates(pair.grid, state)
    for name in stationary.BLOCKS[:3]:
        np.testing.assert_array_equal(raw[name], getattr(original, name))
    source = bundle["source"]
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    lapse = coupling.hamilton_constraint(system, fine) + source["force_L"] / pair.grid.dx_q
    shift = coupling.shift_constraint(system, fine) + source["force_beta"] / pair.grid.dx_q
    np.testing.assert_array_equal(raw["lapse"], galerkin.pull_geometry(pair.grid, lapse))
    np.testing.assert_array_equal(raw["shift"], galerkin.pull_geometry(pair.grid, shift))
    assert constraints["current_max"] < 1e-11
    # Mean Q and physical radius remain independent coordinates.
    moved = x.copy()
    moved[:pair.grid.ng] += .1
    moved[pair.grid.ng:2*pair.grid.ng] += .2
    decoded = stationary.nodal_state(pair.grid, moved)
    np.testing.assert_allclose(decoded.Q, state.Q * np.exp(.1))
    np.testing.assert_allclose(decoded.r, state.r * np.exp(.2))
    assert np.all(decoded.p_Q == 0)


def test_preview_keeps_unsatisfied_relations_and_source_tail(preview):
    record, arrays = preview
    assert record["status"] == "RADIAL_SEED_PREVIEW"
    assert not record["optimizer"]["performed"]
    assert record["normalized_residual_max"] > record["numerical_target"]
    assert record["next_unsatisfied_relation"] in stationary.BLOCKS
    assert not record["scope"]["nonexistence_inferred_from_optimizer"]
    assert not record["measurements"]["rho_positivity_certified"]
    np.testing.assert_allclose(arrays["fine_rho"] - arrays["positive_carrier_rho"], arrays["pointwise_source_tail"])
    hf = record["measurements"]["hellmann_feynman"]
    assert hf["absolute_gap"] < 1e-7


def test_creation_only_and_read_only_authentication(preview, tmp_path):
    record, arrays = preview
    path = tmp_path / "stationary"
    stationary.write_record(record, arrays, path)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    check = stationary.check_record(path)
    assert check["ok"] and check["bytes_written"] == 0 and not check["solver_executed"]
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    with pytest.raises(FileExistsError):
        stationary.write_record(record, arrays, path)
    with pytest.raises(ValueError):
        stationary.output_paths(stationary.LAB / "results/development/old-evidence")
    # Record summaries must agree with payload arrays even if payload hashes match.
    json_path, _ = stationary.output_paths(path)
    changed = json.loads(json_path.read_text())
    changed["measurements"]["quadrature_rho_min"] += 1
    json_path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="source summary"):
        stationary.check_record(path)


def test_existing_payload_alone_is_preserved(tmp_path):
    stem = tmp_path / "partial"
    _, payload = stationary.output_paths(stem)
    payload.write_bytes(b"prior bytes")
    with pytest.raises(FileExistsError):
        stationary.validate_new_output(stem)
    assert payload.read_bytes() == b"prior bytes"


def test_bounded_solver_controls_rejected_before_numerics():
    with pytest.raises(ValueError, match="CPU"):
        stationary.run_stationary(solve=True, cpu_limit=121)


def test_bounded_stop_returns_primary_unsatisfied_measurements(monkeypatch):
    def stop_after_one_evaluation(objective, unknown, **kwargs):
        assert kwargs["max_nfev"] == 1
        objective(unknown)
        raise stationary.BudgetReached("test bounded stop")
    monkeypatch.setattr(stationary, "least_squares", stop_after_one_evaluation)
    report, arrays = stationary.run_stationary(solve=True, max_nfev=1)
    assert report["status"] == "UNSATISFIED_STATIONARY_RELATIONS"
    assert report["optimizer"]["termination"] == "test bounded stop"
    assert report["optimizer"]["actual_residual_evaluations"] == 1
    assert all("residual_" + name in arrays for name in stationary.BLOCKS)
    assert not report["scope"]["nonexistence_inferred_from_optimizer"]


def test_radial_identity_is_euler_variation_of_the_owned_action():
    import sympy as sp
    from recursive_horizons import nsc_spherical_feedback_action as action
    x = sp.symbols("x", real=True)
    A, CW, CF, flux = sp.symbols("A CW CF flux", nonzero=True, real=True)
    r, u, chi = (sp.Function(name)(x) for name in ("r", "u", "chi"))
    F = action.feedback_F(r, chi, A, CW)
    V = action.feedback_V(r, chi, A, CW, CF, flux)
    # Original static conformal SBP energy's continuum density. Vary first;
    # neither chi's equation nor an independently prescribed source is used.
    h = action.feedback_Z(A) * sp.diff(r, x)**2 - sp.exp(2*u) * V + 2 * sp.diff(F, x) * sp.diff(u, x)
    Er = -sp.diff(h, r) + sp.diff(sp.diff(h, sp.diff(r, x)), x)
    static_relation = 3 * sp.diff(r, x, 2) / r + sp.diff(u, x, 2) - sp.exp(2*u)
    assert sp.simplify(Er + 16 * sp.pi * A * r * static_relation) == 0
    R4 = 2 * static_relation / r**2
    assert sp.simplify(Er + 8 * sp.pi * A * r**3 * R4) == 0
