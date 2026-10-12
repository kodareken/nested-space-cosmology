"""Same-action conformal gauge and the actual projected constraint forcing.

These checks independently differentiate the Hamiltonian and constraints.
The forcing estimate controls constraints, not a physical renewal verdict.
"""
from dataclasses import replace

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def _state(nf=16, nq=64, *, top=False):
    grid = galerkin.build_grid(nf, quadrature=nq, length=2 * np.pi, gauge="conformal")
    x, xf = grid.xi_g, grid.xi_f
    phi0 = np.stack([np.exp(1j * k * xf) for k in (-2.5, -1.5, -.5, .5, 1.5, 2.5)], axis=1)
    phi1 = np.stack([(.3 + .2j) * np.exp(1j * (k + 1) * xf)
                     for k in (-2.5, -1.5, -.5, .5, 1.5, 2.5)], axis=1)
    phi0 += .15 * np.roll(phi0, 1, axis=1)
    phi0, phi1 = phi0 / np.sqrt(nf), phi1 / np.sqrt(nf)
    state = galerkin.blank_state(grid, phi0, phi1)
    mode = grid.ng // 2 if top else 1
    state.Q = 1.2 + .07 * np.cos(mode * x)
    state.r = 1.1 + .08 * np.sin(2 * x)
    state.chi = .05 * np.cos(x)
    state.p_Q = .1 + .05 * np.sin(2 * x)
    state.p_r = .08 * np.cos(x)
    state.p_chi = .1 + .03 * np.sin(x)
    return grid, state


def _constraints(grid, state):
    fine = galerkin.prolong_state(grid, state)
    system = replace(grid.fine, length_density=fine.Q, shift=np.zeros(grid.nq))
    source = coupling.source_from_columns(system, fine)
    c = coupling.hamilton_constraint(system, fine) + source["force_L"] / grid.dx_q
    d = coupling.shift_constraint(system, fine) + source["force_beta"] / grid.dx_q
    return fine.Q * c, d


def _displace(state, rate, epsilon):
    out = state.copy()
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        setattr(out, name, getattr(state, name) + epsilon * getattr(rate, name))
    return out


def test_default_path_and_initial_constraints_remain_the_same():
    conformal, state = _state()
    prescribed = galerkin.build_grid(conformal.nf, quadrature=conformal.nq, length=conformal.length)
    assert prescribed.gauge == "prescribed"
    assert galerkin.active_fine_system(prescribed, galerkin.prolong_state(prescribed, state)) is prescribed.fine
    fine = galerkin.prolong_state(conformal, state)
    source_old = coupling.source_from_columns(prescribed.fine, fine)
    active = galerkin.active_fine_system(conformal, fine)
    source_new = coupling.source_from_columns(active, fine)
    assert np.array_equal(source_old["force_L"], source_new["force_L"])
    assert np.array_equal(source_old["force_beta"], source_new["force_beta"])
    for name in ("full_hamilton_max", "full_momentum_max", "projected_hamilton_max", "projected_momentum_max"):
        assert galerkin.constraint_diagnostics(conformal, state)[name] == galerkin.constraint_diagnostics(prescribed, state)[name]
    old_rate, old_bundle = galerkin.compose_fine_hamiltonian(prescribed, state)
    raw = coupling.geometric_rates(prescribed.fine, fine)
    assert np.array_equal(old_rate.p_Q, galerkin.pull_geometry(prescribed, raw[3] - source_old["force_Q"] / prescribed.dx_q))
    assert np.array_equal(old_bundle["lapse_dot"], np.zeros(prescribed.nq))
    assert np.array_equal(prescribed.fine.length_density, conformal.fine.length_density)
    assert not np.array_equal(active.length_density, conformal.fine.length_density)


def test_gauge_fixed_energy_derivative_needs_lapse_chain_rule_and_dynamic_work():
    grid, state = _state(top=True)
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    assert galerkin.hamiltonian_fd_errors(grid, state, eps=2e-6)["relative_worst"] < 2e-6
    fine, system = bundle["fine_state"], bundle["fine_system"]
    source = bundle["source"]
    fixed = coupling.geometric_rates(system, fine)
    required = fixed[3] - coupling.hamilton_constraint(system, fine) - (source["force_Q"] + source["force_L"]) / grid.dx_q
    assert np.max(np.abs(rate.p_Q - galerkin.pull_geometry(grid, required))) < 1e-10
    wrong = galerkin.pull_geometry(grid, fixed[3] - source["force_Q"] / grid.dx_q)
    assert np.max(np.abs(rate.p_Q - wrong)) > 1
    assert np.array_equal(bundle["lapse_dot"], bundle["lifted_Q"])
    dynamic_work = float(np.sum((source["force_Q"] + source["force_L"]) * bundle["lifted_Q"]))
    assert abs(rate.fieldwork_power - dynamic_work) < 1e-10
    legacy_work = float(np.sum(source["force_Q"] * bundle["lifted_Q"]))
    assert abs(dynamic_work - legacy_work) > 1e-3
    eps = 1e-6

    def field_energy(trial):
        lifted = galerkin.prolong_state(grid, trial)
        active = replace(grid.fine, length_density=lifted.Q, shift=np.zeros(grid.nq))
        return coupling.field_energy(active, lifted)

    numerical = (field_energy(_displace(state, rate, eps)) - field_energy(_displace(state, rate, -eps))) / (2 * eps)
    assert abs(numerical - dynamic_work) < 2e-8
    assert galerkin.work_balance(grid, state)["lifted_identity_error"] < 2e-7


def test_actual_constraint_derivative_and_sbp_forcing_budget_are_independent_of_fd_step():
    grid, state = _state(top=True)
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    report = galerkin.conformal_constraint_transport(grid, state, rate, bundle, include_vectors=True)
    for epsilon in (1e-6, 5e-7):
        plus = _constraints(grid, _displace(state, rate, epsilon))
        minus = _constraints(grid, _displace(state, rate, -epsilon))
        for index, name in enumerate(("h_dot", "d_dot")):
            numerical = (plus[index] - minus[index]) / (2 * epsilon)
            assert np.max(np.abs(numerical - report[name])) < 2e-6
    assert report["projection_forcing_norm"] > 1e-3
    assert report["forcing_norm"] > 1e-3
    assert report["sbp_identity_error"] < 1e-9
    for suffix in ("h", "d"):
        assert np.max(np.abs(report[f"forcing_{suffix}"] - report[f"quadrature_{suffix}"] - report[f"projection_{suffix}"])) < 1e-12
    norm_rate = report["constraint_energy_rate"] / (2 * report["constraint_norm"])
    assert abs(norm_rate) <= report["forcing_norm"] + 1e-10
    assert report["continuum_constraint_certified"] is False
    assert report["time_integral_certified"] is False


def test_resolved_smooth_data_approach_the_continuum_transport():
    grid, state = _state(nf=32, nq=128)
    report = galerkin.conformal_constraint_transport(grid, state)
    assert report["quadrature_forcing_norm"] < 1e-8
    assert report["projection_forcing_norm"] < 1e-7
    assert report["forcing_norm"] < 1e-7
    expected_omega = max(abs(grid.modes_f)) + grid.fine.kappa * np.max(galerkin.prolong_geometry(grid, state.Q))
    assert abs(galerkin.subspace_frequency(grid, state)[0] - expected_omega) < 1e-12


def test_api_rejects_unknown_gauge_wrong_dimensions_and_nonpositive_chart():
    with pytest.raises(ValueError, match="gauge"):
        galerkin.build_grid(16, gauge="unknown")
    grid, state = _state()
    with pytest.raises(ValueError, match="conformal gauge"):
        galerkin.conformal_constraint_transport(replace(grid, gauge="prescribed"), state)
    bad = state.copy()
    bad.Q = np.ones(grid.nf)
    with pytest.raises(ValueError, match="shape"):
        galerkin.compose_fine_hamiltonian(grid, bad)
    bad = state.copy()
    bad.Q[:] = -1
    with pytest.raises((ValueError, coupling.PositiveChartExit)):
        galerkin.compose_fine_hamiltonian(grid, bad)


def test_full_continuum_action_has_exact_constraint_transport_with_angular_mass():
    """Independent symbolic Euler derivatives of Q C_total, with free jets.

    F=ar²/2+f chi, Z=3a and V=-ar²-f(4chi+chi²)/2+v0 reproduce
    the owned action; one weighted real/imaginary column suffices because
    arbitrary column weights add linearly.
    """
    x = sp.symbols("x")
    a, f, v0, k, weight = sp.symbols("a f v0 k weight", nonzero=True)
    q, r, chi, p_q, p_r, p_chi, u, v, w, z = fields = [sp.Function(name)(x) for name in
        ("q", "r", "chi", "p_q", "p_r", "p_chi", "u", "v", "w", "z")]
    dx = lambda expr: sp.diff(expr, x)
    ff = a * r ** 2 / 2 + f * chi
    potential = -a * r ** 2 - f * (4 * chi + chi ** 2) / 2 + v0
    pi = p_r - a * r * p_chi / f
    cg = p_q * p_chi / (2 * f) + pi ** 2 / (12 * a * q) + 3 * a * dx(r) ** 2 / q - q * potential - 2 * dx(dx(ff) / q)
    kinetic = -u * dx(w) - v * dx(z) + w * dx(u) + z * dx(v)
    scalar = 2 * (u * w + v * z)
    current = u * dx(v) - v * dx(u) + w * dx(z) - z * dx(w)
    h = sp.expand(q * cg + weight * (kinetic + k * q * scalar))
    d = p_r * dx(r) + p_chi * dx(chi) - q * dx(p_q) - weight * current

    def jet_partial(expr, field, order):
        jet = sp.diff(field, x, order)
        symbol = sp.Dummy()
        return sp.diff(expr.xreplace({jet: symbol}), symbol).xreplace({symbol: jet})

    def euler(field):
        return sp.simplify(sum((-1) ** order * sp.diff(jet_partial(h, field, order), x, order) for order in range(3)))

    rates = {q: euler(p_q), r: euler(p_r), chi: euler(p_chi),
             p_q: -euler(q), p_r: -euler(r), p_chi: -euler(chi),
             u: euler(v) / (2 * weight), v: -euler(u) / (2 * weight),
             w: euler(z) / (2 * weight), z: -euler(w) / (2 * weight)}

    def variation(expr):
        return sum(jet_partial(expr, field, order) * sp.diff(rates[field], x, order)
                   for field in fields for order in range(4))

    assert sp.factor(sp.expand(variation(h) - dx(d)).doit()) == 0
    assert sp.factor(sp.expand(variation(d) - dx(h)).doit()) == 0
