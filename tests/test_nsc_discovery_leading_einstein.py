"""New finite action, canonical maps, gauge, principal symbol and short controls."""
from dataclasses import replace
import json
import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_leading_einstein as le


def manufactured(nf=32, *, vacuum=False):
    grid = le.galerkin.build_grid(nf, gauge="conformal")
    phi0, phi1 = le.galerkin.manufactured_columns(nf)
    pair = le.nested.build_pair(nf, coarse_modes=1, child_details=2, source_layout="override",
                               columns_override=(phi0, phi1))
    x = grid.xi_g
    nodal = le.State(.3+.01*np.cos(2*np.pi*x/grid.length),
                    1.4+.03*np.sin(2*np.pi*x/grid.length),
                    .02*np.cos(2*np.pi*x/grid.length), .03*np.sin(2*np.pi*x/grid.length),
                    np.zeros_like(phi0) if vacuum else phi0,
                    np.zeros_like(phi1) if vacuum else phi1)
    return pair, le.encode(pair, nodal)


def block_exact_auxiliary(monkeypatch):
    for module, name in ((le.galerkin, "compose_fine_hamiltonian"), (le.galerkin, "rates"),
                         (le.coupling, "geometric_rates"), (le.coupling, "hamilton_constraint"),
                         (le.tidal, "analytic_accelerations")):
        monkeypatch.setattr(module, name, lambda *_args, **_kwargs: pytest.fail("exact auxiliary dynamics called"))


def test_actual_finite_energy_gradient_and_constraint_identity(monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1):
        pair, state = manufactured()
        rate = le.rates(pair, state)
        rng = np.random.default_rng(18)
        for name, conjugate, sign in (("Q", "p_Q", -1), ("r", "p_r", -1),
                                      ("p_Q", "Q", 1), ("p_r", "r", 1)):
            direction = rng.normal(size=pair.grid.ng)
            direction /= np.linalg.norm(direction)
            tangent = state.copy()
            for key in le.FIELDS:
                setattr(tangent, key, np.zeros_like(getattr(tangent, key)))
            setattr(tangent, name, direction)
            numeric = (le.energy(pair, le.combine(state, tangent, 1e-6))["total"]
                       -le.energy(pair, le.combine(state, tangent, -1e-6))["total"])/2e-6
            analytic = sign*np.dot(getattr(rate, conjugate), direction)
            assert numeric == pytest.approx(analytic, rel=2e-7, abs=1e-6)
        fine = le.fine_state(pair, state)
        system = le.active_system(pair.grid, fine)
        C, _ = le.constraint_arrays(pair.grid, fine, system)
        assert le.energy(pair, state)["total"] == pytest.approx(pair.grid.dx_q*np.sum(fine.Q*C), abs=1e-9)
        nodal = le.decode(pair, state)
        assert np.allclose(le.encode(pair, nodal).p_Q, state.p_Q, atol=1e-12)
        assert not hasattr(state, "chi") and not hasattr(state, "p_chi")


def test_new_Jv_projected_jets_and_fft_match_dense(monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1), le.backend.fft_thread_limit(1):
        pair, state = manufactured()
        rng = np.random.default_rng(29)
        tangent = le.State(*(rng.normal(size=getattr(state, name).shape)*.02 for name in le.FIELDS))
        tangent.phi0 = tangent.phi0+1j*rng.normal(size=tangent.phi0.shape)*.01
        tangent.phi1 = tangent.phi1+1j*rng.normal(size=tangent.phi1.shape)*.01
        analytic = le.jvp(pair, state, tangent)
        plus, minus = le.rates(pair, le.combine(state, tangent, 1e-6)), le.rates(pair, le.combine(state, tangent, -1e-6))
        for name in le.FIELDS:
            np.testing.assert_allclose(getattr(analytic, name),
                (getattr(plus, name)-getattr(minus, name))/2e-6, rtol=2e-5, atol=2e-6)
        fft_pair = le.backend.make_fft_pair(pair)
        dense, fft = le.rates(pair, state), le.rates(fft_pair, state)
        for name in le.FIELDS:
            np.testing.assert_allclose(getattr(dense, name), getattr(fft, name), rtol=1e-9, atol=1e-9)
        curvature, _velocity, acceleration = le.metric_jets(pair, state)
        along = le.rates(pair, state)
        before = le.rates(pair, le.combine(state, along, -1e-6))
        after = le.rates(pair, le.combine(state, along, 1e-6))
        numeric = le.prolong(pair.grid, le.decode(pair, le.State(*(
            (getattr(after, name)-getattr(before, name))/2e-6 for name in le.FIELDS))))
        np.testing.assert_allclose(acceleration.Q, numeric.Q, rtol=1e-5, atol=1e-6)
        assert np.isfinite(curvature["base"]["K"]).all()


def test_conformal_guard_and_principal_AB_identity():
    with threadpool_limits(limits=1):
        pair, state = manufactured(vacuum=True)
        nodal = le.decode(pair, state)
        fine = le.prolong(pair.grid, nodal)
        bad_grid = replace(pair.grid, gauge="prescribed")
        with pytest.raises(ValueError, match="conformal"):
            le.active_system(bad_grid, fine)
        with pytest.raises(ValueError, match="L=Q"):
            le.fine_rates(pair.grid.fine, fine)
        a, Z, _ = le.coefficients(pair.grid.fine)
        q, r, wave = .3, 1.4, 2*np.pi/pair.grid.length
        b = a*r
        A = np.array([[-Z*q*q/(2*b*b), q/(2*b)], [q/(2*b), 0.]])
        B = np.array([[0., 2*b/q], [2*b/q, 2*Z]])
        np.testing.assert_allclose(A@B, np.eye(2), atol=1e-12)
        # Execute the new SBP Jv on a cosine and subtract its zero-mode
        # reaction to recover the derivative momentum block, not just algebra.
        phi = np.zeros_like(nodal.phi0)
        constant = le.State(np.full(pair.grid.ng, q), np.full(pair.grid.ng, r),
                            np.zeros(pair.grid.ng), np.zeros(pair.grid.ng), phi, phi.copy())
        fixed = le.encode(pair, constant)
        mode = np.cos(wave*pair.grid.xi_g)
        for column, field in enumerate(("Q", "r")):
            images = []
            for pattern in (mode, np.ones_like(mode)):
                t = le.State(*(pattern if name==field else np.zeros_like(getattr(constant, name)) for name in le.FIELDS))
                image = le.decode(pair, le.jvp(pair, fixed, le.encode(pair, t)))
                images.append(np.column_stack((image.p_Q, image.p_r)))
            expected = -wave**2*mode[:, None]*B[:, column][None, :]
            np.testing.assert_allclose(images[0]-mode[:, None]*images[1], expected, atol=1e-8)
        dt, info = le.step_restriction(pair, fixed, .001)
        assert 0 < dt <= .001 and not info["field_cfl_alone"]


def test_short_magnetic_BR_control_and_chart_stop(monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1):
        pair, template = manufactured(vacuum=True)
        q0 = pair.grid.fine.calibration["b0"]/pair.grid.fine.calibration["a0"]
        a, _Z, _mag = le.coefficients(pair.grid.fine)
        phi = np.zeros_like(template.phi0)
        zero = np.zeros(pair.grid.ng)
        nodal = le.State(np.full(pair.grid.ng, q0), np.ones(pair.grid.ng), zero, zero.copy(), phi, phi.copy())
        state = le.encode(pair, nodal)
        for _ in range(5):
            state = le.rk4_step(pair, state, .001)
        actual = le.decode(pair, state)
        t = .005
        q = q0/np.cosh(q0*t)
        D = -q0*np.tanh(q0*t)
        np.testing.assert_allclose(actual.Q, q, atol=1e-10)
        np.testing.assert_allclose(actual.r, 1., atol=1e-10)
        np.testing.assert_allclose(actual.p_r, 2*a*D, atol=1e-9)
        assert abs(le.energy(pair, state)["total"]) < 1e-8
        assert le.constraints(pair, state)["h_c_max"] < 1e-7
        curvature, _, _ = le.metric_jets(pair, state)
        np.testing.assert_allclose(curvature["tides"]["R_0101"], -1., atol=1e-7)
        np.testing.assert_allclose(curvature["base"]["K"], 8., atol=1e-6)
        bad = state.copy()
        bad.Q = -bad.Q
        with pytest.raises(le.coupling.PositiveChartExit):
            le.rk4_step(pair, bad, .001)


def test_existing_runner_adapter_checkpoint_and_restoration(tmp_path, monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1), le.backend.fft_thread_limit(1):
        pair, state = manufactured()
        saved = le.episode.execute_case
        with le.episode_adapter():
            arrays = {key: getattr(state, name) for name, key in
                      (("Q", "Q"), ("r", "r"), ("p_Q", "pi_Q"), ("p_r", "pi_r"), ("phi0", "phi0"), ("phi1", "phi1"))}
            arrays.update(W=pair.geometry_map, source_phi0=pair.source_phi0, source_phi1=pair.source_phi1,
                observer_columns=pair.reference_columns, source_weights=pair.weights,
                normal_clocks=np.zeros(3), clock_rates=le.clock_rates(pair, state))
            record = {"case_id": "tiny", "nf": 32, "coordinate_time": 0., "steps": 0,
                "stations": [.002], "step_cap": .001, "control_mode": "coupled", "snapshot_kind": "handoff",
                "source_pins": le.episode._live_pins(arrays, state), "status": "PREPARED",
                "coarse_indices": pair.geometry_coarse_indices.tolist(), "child_indices": pair.geometry_child_indices.tolist(),
                "parent_indices": pair.geometry_parent_indices.tolist()}
            le.episode.commit_checkpoint(tmp_path, record, arrays)
            outcome = le.execute_case(tmp_path, "tiny", cpu_allowance=5., max_steps=2)
            assert outcome["coordinate_time"] == pytest.approx(.002)
            last, values = le.episode.load_checkpoint(tmp_path, "tiny")
            assert "chi" not in values and "pi_chi" not in values
            assert values["normal_clocks"].min() > 0
            assert le.episode.check(tmp_path) if (tmp_path/"manifest.json").exists() else last["immutable"]
        assert le.episode.execute_case is saved
    preflight = le.prepare()
    assert preflight["status"] == "PREFLIGHT" and len(preflight["cases"]) == 8
    assert not preflight["evolved"]


def test_matched_frozen_controller_jets_metric_and_field_conservation(monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1), le.backend.fft_thread_limit(1):
        pair, state = manufactured()
        original = state.copy()
        before_energy = le.energy(pair, state)
        before_gram = state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1
        rate, bundle = le.rates(pair, state, control_mode="frozen_geometry", return_bundle=True)
        coupled = le.rates(pair, state)
        assert max(np.max(abs(getattr(coupled, name))) for name in le.REAL_FIELDS) > 0
        for name in le.REAL_FIELDS:
            assert np.array_equal(getattr(rate, name), np.zeros_like(getattr(rate, name)))
        assert rate.fieldwork_power == 0
        assert np.max(abs(rate.phi0)) > 0
        initial_curvature, velocity, acceleration = le.metric_jets(pair, state, rate, bundle)
        for name in ("Q", "r"):
            assert np.max(abs(getattr(velocity, name))) == 0
            assert np.max(abs(getattr(acceleration, name))) == 0
        for _ in range(5):
            state = le.frozen_geometry_step(pair, state, .001)
        for name in le.REAL_FIELDS:
            assert np.array_equal(getattr(state, name), getattr(original, name))
        assert not np.array_equal(state.phi0, original.phi0)
        after_gram = state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1
        np.testing.assert_allclose(after_gram, before_gram, atol=1e-10)
        assert le.energy(pair, state)["field"] == pytest.approx(before_energy["field"], abs=1e-8)
        row = le.observe(pair, state, .005, "frozen_geometry")
        assert row["control_mode"] == "frozen_geometry" and row["controlled_metric"]
        assert not row["constraint_preservation_claim"]
        assert row["inactive_geometry_forcing"]["p_Q"]["max_abs"] > 0
        assert row["coordinate_fieldwork"] == row["pressure_work"] == 0
        assert not row["coordinate_readout_is_equal_proper_time"]
        np.testing.assert_array_equal(le.clock_rates(pair, state), le.clock_rates(pair, original))
        later_curvature, _, _ = le.metric_jets(pair, state, control_mode="frozen_geometry")
        for name in ("R_0101", "R_0202"):
            np.testing.assert_array_equal(later_curvature["tides"][name], initial_curvature["tides"][name])
        cases = le.plan()
        assert [row["case_id"] for row in cases[:6]] == [
            f"nf{nf}_{name}_dt{cap}" for nf, name in ((128, "uniform"), (128, "coherent"), (256, "coherent"))
            for cap in ("0.001", "0.0005")]
        assert all(row["initial_case"]=="coherent" and row["nf"]==128
                   and row["control_mode"]=="frozen_geometry" for row in cases[6:])


def test_existing_advance_uses_frozen_four_variable_step_and_exact_clocks(tmp_path, monkeypatch):
    block_exact_auxiliary(monkeypatch)
    with threadpool_limits(limits=1), le.backend.fft_thread_limit(1), le.episode_adapter():
        pair, state = manufactured()
        arrays = {key: getattr(state, name) for name, key in
                  (("Q", "Q"), ("r", "r"), ("p_Q", "pi_Q"), ("p_r", "pi_r"), ("phi0", "phi0"), ("phi1", "phi1"))}
        initial_rates = le.clock_rates(pair, state)
        arrays.update(W=pair.geometry_map, source_phi0=pair.source_phi0, source_phi1=pair.source_phi1,
            observer_columns=pair.reference_columns, source_weights=pair.weights,
            normal_clocks=np.zeros(3), clock_rates=initial_rates)
        record = {"case_id": "tiny_frozen", "nf": 32, "coordinate_time": 0., "steps": 0,
            "stations": [.002], "step_cap": .001, "control_mode": "frozen_geometry", "snapshot_kind": "handoff",
            "source_pins": le.episode._live_pins(arrays, state), "status": "PREPARED",
            "coarse_indices": pair.geometry_coarse_indices.tolist(), "child_indices": pair.geometry_child_indices.tolist(),
            "parent_indices": pair.geometry_parent_indices.tolist()}
        le.episode.commit_checkpoint(tmp_path, record, arrays)
        result = le.execute_case(tmp_path, "tiny_frozen", cpu_allowance=5., max_steps=2)
        assert result["coordinate_time"] == pytest.approx(.002)
        last, final = le.episode.load_checkpoint(tmp_path, "tiny_frozen")
        for name in ("Q", "r", "pi_Q", "pi_r"):
            np.testing.assert_array_equal(final[name], arrays[name])
        np.testing.assert_allclose(final["normal_clocks"], .002*initial_rates, atol=1e-14)
        assert last["control_mode"] == "frozen_geometry"
        assert last["channel_sample"]["controlled_metric"]
