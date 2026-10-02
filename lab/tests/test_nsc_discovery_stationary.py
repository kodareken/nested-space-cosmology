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


def test_exact_discrete_virial_with_noncritical_geometry():
    pair = stationary.make_pair()
    x, _ = stationary.initial_unknown(pair)
    ng = pair.grid.ng
    angle = 2 * np.pi * pair.grid.xi_g / pair.grid.length
    x[ng:2*ng] += .07 * np.cos(2 * angle)
    x[2*ng:] += .4 * np.sin(3 * angle)
    state, spec, raw, _, bundle, constraints = stationary.evaluate(pair, x)
    virial = stationary.virial_diagnostic(pair, state, raw, bundle, constraints, spec)
    assert abs(virial['gap']) > 1
    assert abs(virial['fine_identity_defect']) < 1e-10
    assert abs(virial['coarse_identity_defect']) < 1e-10
    assert abs(virial['adjoint_pairing_defect']) < 1e-10
    assert abs(virial['source_owned_field_energy_gap']) < 1e-11
    assert abs(virial['source_spectral_field_energy_gap']) < 1e-10
    assert virial['necessary_only'] and not virial['stationary_balance_claimed']
    state.p_Q[0] = .1
    with pytest.raises(ValueError, match='zero momenta'):
        stationary.virial_diagnostic(pair, state, raw, bundle, constraints, spec)


def test_finite_radius_affinity_and_unchanged_source():
    pair = stationary.make_pair()
    x, _ = stationary.initial_unknown(pair)
    state, spectral, *_ = stationary.evaluate(pair, x)
    original = state.copy()
    report, arrays = stationary.select_radius_affine(pair, state, spectral)
    assert report['finite_affinity_defect_max'] < 1e-10
    assert report['source_change_max'] == 0
    predicted = arrays['radius_pQ_remainder'] + 2.25 * arrays['radius_pQ_coefficient']
    np.testing.assert_allclose(arrays['radius_third_pQ'], predicted, atol=1e-10)
    for name in stationary.nested.STATE_NAMES:
        np.testing.assert_array_equal(getattr(state, name), getattr(original, name))
    assert not report['stationary_balance_claimed']
    if report['least_squares_amplitude_square'] <= 0:
        assert report['branch'] == 'nonpositive_radius_square_branch_mismatch'
        assert 'radius_selected_r' not in arrays
    else:
        assert report['branch'] == 'positive_conditional_radius'
        assert report['radius_amplitude']**2 == pytest.approx(report['least_squares_amplitude_square'])
        assert set(report['raw_residual_max']) == set(stationary.BLOCKS)


def test_higher_harmonic_trial_uses_actual_fine_virial():
    # One small-band unit query exercises the same trial helper. It does not
    # execute the production NF128 family or a joint geometry solve.
    pair = stationary.make_pair()
    report, arrays = stationary.balance_trial(pair, harmonic=3, amplitude=.35)
    Q, chi = arrays['fine_Q'], arrays['fine_chi']
    system = pair.grid.fine
    expected = (float(coupling.alpha_of(system.C_W)) * pair.grid.dx_q * np.sum(Q**2 * chi**2)
                - pair.grid.dx_q * np.dot(Q, arrays['fine_rho'])
                - 2*np.pi*system.C_F*system.flux**2*pair.grid.dx_q*np.sum(Q**2))
    assert report['virial_gap'] == pytest.approx(expected, abs=1e-11)
    assert report['seed']['harmonic'] == 3
    assert report['measurements']['fine_chi_seed_relation_defect_max'] >= 0
    assert abs(report['measurements']['virial']['fine_identity_defect']) < 1e-10
    assert set(report['measurements']['fine_momentum_residual_max']) == set(stationary.BLOCKS[:3])
    assert not report['measurements']['hellmann_feynman']['performed']


def test_balance_brackets_never_skip_failed_endpoint():
    measured = lambda amp, gap: {'harmonic':6, 'amplitude':amp, 'status':'MEASURED', 'virial_gap':gap}
    first, last = measured(.025, -1.), measured(.35, 1.)
    failed = {'harmonic':6, 'amplitude':.1, 'status':'FAILED'}
    assert stationary.sign_change_brackets([first, failed, last]) == []
    middle = measured(.1, .2)
    assert stationary.sign_change_brackets([first, middle, last]) == [(first, middle)]


def test_balance_preflight_runs_no_numerical_trials(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('preflight must not construct matrices or execute trials')
    monkeypatch.setattr(stationary, 'make_pair', forbidden)
    monkeypatch.setattr(stationary, 'balance_trial', forbidden)
    report, arrays = stationary.run_balance_seed()
    assert report['status'] == 'BALANCE_SEED_PREFLIGHT'
    assert report['nf'] == 128 and report['preflight']['endpoint_trials'] == 16
    assert report['preflight']['harmonics'] == [6,8,10,12]
    assert report['preflight']['endpoint_amplitudes'] == [.025,.1,.35,.65]
    assert report['actual_trials'] == 0 and arrays == {}


def test_balance_execution_cache_and_brent_admission_without_campaign(monkeypatch):
    # Synthetic scalar trial oracle tests scheduling/caching only. No model
    # state or scientific record is created from this oracle.
    calls = []
    monkeypatch.setattr(stationary, 'make_pair', lambda nf: object())
    def scalar_trial(pair, harmonic, amplitude):
        calls.append((harmonic, amplitude))
        gap = amplitude-.3 if harmonic == 6 else 1+amplitude
        return {'harmonic':harmonic, 'amplitude':amplitude, 'status':'MEASURED',
                'virial_gap':gap, 'measurements':{'virial':{'gap':gap}}}, {'test_array':np.ones(1)}
    monkeypatch.setattr(stationary, 'balance_trial', scalar_trial)
    monkeypatch.setattr(stationary, 'least_squares', lambda *args, **kwargs: pytest.fail('joint solve forbidden'))
    report, _ = stationary.run_balance_seed(execute=True)
    assert len(calls) == len(set(calls))
    assert len(report['endpoint_assessment']) == 16
    assert len(report['brackets']) == 1
    root = report['roots'][0]
    assert root['status'] == 'VIRIAL_SELECTED_SHAPE' and root['amplitude'] == pytest.approx(.3)
    assert report['status'] == 'VIRIAL_SELECTED_SHAPE_ONLY'
    assert not report['scope']['stationary_balance_claimed']
    assert report['actual_trials'] <= stationary.MAX_BALANCE_TRIALS


def test_balance_cost_forecast_preserves_unadmitted_and_failed_trials(monkeypatch):
    clock = [0.]
    monkeypatch.setattr(stationary.time, 'process_time', lambda: clock[0])
    monkeypatch.setattr(stationary, 'make_pair', lambda nf: object())
    def costly_scalar_trial(pair, harmonic, amplitude):
        clock[0] += .25
        if amplitude == .1:
            raise ValueError('synthetic failed endpoint')
        return {'harmonic':harmonic, 'amplitude':amplitude, 'status':'MEASURED',
                'virial_gap':1., 'measurements':{'virial':{'gap':1.}}}, {'test_array':np.ones(1)}
    monkeypatch.setattr(stationary, 'balance_trial', costly_scalar_trial)
    report, _ = stationary.run_balance_seed(execute=True, cpu_limit=1.)
    assert len(report['endpoint_assessment']) == 16
    assert report['actual_trials'] == 3
    assert report['trials'][1]['status'] == 'FAILED'
    assert any(row['status'] == 'NOT_ADMITTED' for row in report['endpoint_assessment'])
    assert report['preflight']['forecast_endpoint_cpu_seconds'] == 4.
    assert report['cpu_seconds'] == .75 and not report['cpu_budget_exceeded']
    assert report['status'] == 'NO_COMPLETED_VIRIAL_BRACKET'


def test_sealed_v1_authenticates_historical_bytes_and_diagnoses_without_solve(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('saved diagnosis must not optimize or select a new eigensource')
    monkeypatch.setattr(stationary, 'least_squares', forbidden)
    monkeypatch.setattr(stationary, 'spectral_source', forbidden)
    check = stationary.check_record(stationary.LEGACY_OUTPUT)
    assert check['ok'] and check['schema'] == stationary.LEGACY_SCHEMA
    assert check['source_authentication'] == 'immutable_producing_commit'
    report, arrays = stationary.run_saved_diagnostic()
    measurements = report['measurements']
    assert report['status'] == 'SAVED_STATIONARY_DIAGNOSTIC'
    assert not report['optimizer']['executed_now']
    assert measurements['virial']['gap'] == pytest.approx(-12.9592469545, abs=1e-8)
    assert measurements['fine_momentum_residual_max']['p_r'] == pytest.approx(.45172235, abs=1e-7)
    assert measurements['fine_momentum_residual_max']['p_chi'] == pytest.approx(.00296554, abs=1e-7)
    assert abs(measurements['virial']['coarse_identity_defect']) < 1e-10
    assert 'fine_p_r' in arrays and 'fine_p_chi' in arrays


def test_successor_paths_remain_creation_only(tmp_path):
    successor = stationary.LAB / 'results/development/nsc-discovery-stationary-balance-seed-v2'
    assert stationary.output_paths(successor)[0].name == 'nsc-discovery-stationary-balance-seed-v2.json'
    assert stationary.output_paths('lab/results/development/nsc-discovery-stationary-v1')[0] == stationary.LEGACY_OUTPUT.with_suffix('.json')
    with pytest.raises(FileExistsError):
        stationary.validate_new_output(stationary.LEGACY_OUTPUT)


def test_saved_diagnostic_successor_round_trip(tmp_path):
    report, arrays = stationary.run_saved_diagnostic()
    stem = tmp_path / 'sealed-diagnostic-v2'
    stationary.write_record(report, arrays, stem)
    check = stationary.check_record(stem)
    assert check['ok'] and check['schema'] == stationary.SCHEMA
    assert check['status'] == 'SAVED_STATIONARY_DIAGNOSTIC'


def test_balance_preflight_round_trip(tmp_path):
    report, arrays = stationary.run_balance_seed()
    stem = tmp_path / 'balance-preflight-v2'
    stationary.write_record(report, arrays, stem)
    assert stationary.check_record(stem)['status'] == 'BALANCE_SEED_PREFLIGHT'
