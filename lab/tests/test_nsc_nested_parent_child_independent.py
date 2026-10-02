"""Independent finite-action and nesting checks; no trajectory is produced.

Generic admissible data exercise nonzero reciprocal source forces. Compact
SOURCE supports and overlapping OBSERVER links are deliberately distinct.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture(scope="module")
def pair():
    return nested.build_pair(64)


@pytest.fixture(scope="module")
def driver():
    path = Path(__file__).resolve().parents[1] / "scripts/derive_nsc_nested_parent_child.py"
    spec = importlib.util.spec_from_file_location("nested_pair_independent_driver", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _state(pair):
    x = pair.grid.xi_g
    angle = 2 * np.pi * x / pair.grid.length
    local = pair.geometry_map[:, pair.geometry_child_indices]
    rng = np.random.default_rng(878721)
    columns, _ = np.linalg.qr(rng.normal(size=(2 * pair.grid.nf, 6)) +
                              1j * rng.normal(size=(2 * pair.grid.nf, 6)))
    nodal = coupling.CauchyState(
        Q=.3 + .015 * np.cos(angle) + .01 * local[:, 0],
        r=3.2 + .08 * np.cos(2 * angle) + .01 * local[:, 1],
        chi=.02 * np.sin(angle) + .01 * local[:, 2],
        p_Q=.002 + .001 * np.sin(angle) + .001 * local[:, 2],
        p_r=.002 * np.cos(2 * angle) + .002 * local[:, 1],
        p_chi=.0004 + .0002 * np.sin(2 * angle) + .001 * local[:, 0],
        phi0=columns[:pair.grid.nf], phi1=columns[pair.grid.nf:])
    return nested.encode_state(pair, nodal)


def _columns(state):
    return np.vstack((state.phi0, state.phi1))


def _field_energy(pair, state):
    columns = _columns(state)
    image = nested.apply_hamiltonian(pair, state, columns)
    return float(4 * pair.grid.fine.kappa * np.dot(
        pair.weights, np.sum(columns.conj() * image, axis=0).real))


def _directional_energy(pair, state, rate, step=2e-7):
    plus, minus = state.copy(), state.copy()
    for name in nested.STATE_NAMES:
        setattr(plus, name, getattr(state, name) + step * getattr(rate, name))
        setattr(minus, name, getattr(state, name) - step * getattr(rate, name))
    return (nested.energy(pair, plus) - nested.energy(pair, minus)) / (2 * step)


def test_parent_intervention_is_spatially_exterior_and_preserves_child_covariance(pair):
    source, observer = pair.source_columns, pair.reference_columns
    child = observer[:, pair.child_indices]
    x = pair.grid.xi_f
    in_child = (x >= 1) & (x <= 3)
    exterior = source[:, [0, 1, 4, 5]]
    spinor_mask = np.r_[in_child, in_child]
    assert np.max(abs(exterior[spinor_mask])) < 1e-13
    assert np.max(abs(child.conj().T @ exterior)) < 1e-13
    initial = (child.conj().T @ source * pair.weights) @ (source.conj().T @ child)
    changed = pair.weights + np.array([.025, .025, 0., 0., -.025, -.025])
    control = (child.conj().T @ source * changed) @ (source.conj().T @ child)
    assert changed.sum() == pytest.approx(pair.weights.sum())
    assert np.max(abs(control - initial)) < 1e-14
    # The complete state really changes; invariant child data is not a no-op.
    difference = (source * (changed - pair.weights)) @ source.conj().T
    assert np.linalg.norm(difference) > .04
    fine0, fine1 = pair.grid.U_f @ pair.source_phi0, pair.grid.U_f @ pair.source_phi1
    power = abs(fine0) ** 2 + abs(fine1) ** 2
    measured = []
    for packet, (left, right) in enumerate(((0., 1.), (1., 3.), (3., 4.))):
        outside = (pair.grid.xi_q < left) | (pair.grid.xi_q > right)
        for index in (2 * packet, 2 * packet + 1):
            measured.append(np.sum(power[outside, index]) / np.sum(power[:, index]))
    np.testing.assert_allclose(measured, pair.source_metadata["fine_outside_support_power_fraction"],
                               rtol=2e-14, atol=1e-16)
    assert pair.source_metadata["exact_spatial_support"] is False


def test_live_observer_link_is_not_a_disjoint_source_tail_link(pair):
    state = _state(pair)
    observer = pair.reference_columns
    child, detail = observer[:, [2, 3]], observer[:, [0, 1, 4, 5]]
    actual = child.conj().T @ nested.apply_hamiltonian(pair, state, detail)
    np.testing.assert_allclose(nested.live_blocks(pair, state)["child_parent_link"],
                               actual, rtol=1e-12, atol=1e-12)
    assert np.linalg.norm(actual) > .01
    # A local differential stencil has zero overlap on these disjoint compact
    # source supports. Its result is NOT used as the production Hamiltonian.
    outer = pair.source_columns[:, [0, 1, 4, 5]]
    f0, f1 = outer[:pair.grid.nf], outer[pair.grid.nf:]
    d0 = (np.roll(f0, -1, axis=0) - np.roll(f0, 1, axis=0)) / (2 * pair.grid.dx_f)
    d1 = (np.roll(f1, -1, axis=0) - np.roll(f1, 1, axis=0)) / (2 * pair.grid.dx_f)
    image = np.vstack((-d1 + .3 * f1, d0 + .3 * f0))
    assert np.max(abs(child.conj().T @ image)) < 1e-13
    # Inserting the measured observer link into this source split is a mutation.
    assert np.linalg.norm(actual - child.conj().T @ image) > .01


def test_autonomous_full_energy_and_omitted_source_mutation(pair):
    state = _state(pair)
    full = nested.rates(pair, state)
    without_source = nested.rates(pair, state, include_matter_force=False)
    derivative = _directional_energy(pair, state, full)
    bad_derivative = _directional_energy(pair, state, without_source)
    assert abs(derivative) < 2e-6
    assert abs(full.fieldwork_power) > 1e-5
    assert bad_derivative == pytest.approx(full.fieldwork_power, rel=2e-4, abs=2e-6)
    assert abs(bad_derivative) > 10 * max(abs(derivative), 1e-7)


def test_actual_canonical_forces_include_the_cross_complete_field_trace(pair):
    state = _state(pair)
    gradient = nested.source_geometry_forces(pair, state)["Q_energy_gradient"]
    for group in (pair.geometry_parent_indices, pair.geometry_child_indices):
        index = group[np.argmax(abs(gradient[group]))]
        assert abs(gradient[index]) > 1e-4
        plus, minus = state.copy(), state.copy()
        plus.Q[index] += 1e-6
        minus.Q[index] -= 1e-6
        numerical = (_field_energy(pair, plus) - _field_energy(pair, minus)) / 2e-6
        assert numerical == pytest.approx(gradient[index], rel=2e-6, abs=3e-8)
    full = nested.rates(pair, state)
    without = nested.rates(pair, state, include_matter_force=False)
    np.testing.assert_allclose(full.p_Q - without.p_Q, -gradient, rtol=1e-11, atol=1e-9)
    # pi carries dx_g; using the old nodal momentum convention here is wrong.
    assert np.linalg.norm((full.p_Q - without.p_Q) / pair.grid.dx_g + gradient) > .01


def test_full_link_energy_is_counted_once_and_parent_already_contains_child(pair):
    state = _state(pair)
    columns = _columns(state)
    covariance = (columns * pair.weights) @ columns.conj().T
    matrix = nested.hamiltonian(pair, state)
    factor = 4 * pair.grid.fine.kappa
    direct = float(factor * np.trace(covariance @ matrix).real)
    accounting = nested.energy_accounting(pair, state)
    assert accounting["field"] == pytest.approx(direct, abs=3e-11)
    names = ("child", "parent_detail", "ambient")
    stored = sum(accounting["field_" + name] for name in names)
    stored += sum(accounting["field_" + left + "_" + right + "_cross"]
                  for index, left in enumerate(names) for right in names[index + 1:])
    assert stored == pytest.approx(direct, abs=3e-11)
    parent = pair.reference_columns @ pair.reference_columns.conj().T
    parent_only = float(factor * np.trace(parent @ covariance @ parent @ matrix).real)
    assert accounting["field_parent"] == pytest.approx(parent_only, abs=3e-11)
    diagonal = sum(accounting["field_" + name] for name in names)
    assert abs(diagonal - direct) > 1e-3
    # Parent aggregate is a report, not an additional independent energy store.
    assert abs(stored + parent_only - direct) > 1e-3


def test_driver_full_normal_shell_slope_and_disjoint_window_accounting(pair, driver):
    state = _state(pair)
    rate, bundle = nested.rates(pair, state, return_bundle=True)
    lifted = driver.lifted_rate(pair, rate, bundle)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    analytic = driver.normal_energy_slope(system, fine, lifted, bundle["source"])
    step = 2e-7
    plus, minus = fine.copy(), fine.copy()
    for name in nested.STATE_NAMES:
        setattr(plus, name, getattr(fine, name) + step * getattr(lifted, name))
        setattr(minus, name, getattr(fine, name) - step * getattr(lifted, name))
    def shell_density(point):
        return coupling.source_from_columns(system, point)["force_L"] / point.r
    numerical = (shell_density(plus) - shell_density(minus)) / (2 * step)
    np.testing.assert_allclose(analytic, numerical, rtol=3e-6, atol=3e-8)
    hierarchy = driver.response.nested_frame(pair.reference_columns)
    row = driver.observe(pair, state, 0., hierarchy)
    assert row["field_energy"] == pytest.approx(_field_energy(pair, state), abs=3e-11)
    windows = row["windows"]
    for field in ("normal_shell", "probability", "boundary_inflow", "pressure_work",
                  "lapse_gradient_work", "balance_defect"):
        disjoint = sum(windows[name][field] for name in
                       ("child", "left_annulus", "right_annulus"))
        assert windows["parent"][field] == pytest.approx(disjoint, abs=2e-9)


def test_driver_control_declarations_preserve_the_named_initial_data(pair, driver):
    baseline = np.array(driver.WEIGHTS["baseline"])
    parent = np.array(driver.WEIGHTS["parent_only"])
    child = np.array(driver.WEIGHTS["child_only"])
    assert baseline.sum() == pytest.approx(3.)
    assert parent.sum() == pytest.approx(3.)
    assert child.sum() == pytest.approx(3.2)
    assert np.array_equal(parent[2:4], baseline[2:4])
    assert np.array_equal(child[[0, 1, 4, 5]], baseline[[0, 1, 4, 5]])
    retained = pair.reference_columns[:, pair.child_indices]
    overlap = retained.conj().T @ pair.source_columns
    local = lambda weights: (overlap * weights) @ overlap.conj().T
    np.testing.assert_allclose(local(parent), local(baseline), rtol=0., atol=1e-14)
    assert np.linalg.norm(local(child) - local(baseline)) > .1


def test_two_way_assessment_uses_increment_effect_and_own_refinement_fraction(driver):
    times = np.array([0., .5, 1.])
    results = {}
    coarse, fine = driver.COARSE_NF, driver.FINE_NF
    for nf in (coarse, fine):
        for cap in (.001, .0005):
            scale = 1 + (.005 if nf == coarse else 0) + (.002 if cap == .001 else 0)
            for control in ("baseline", "parent_only", "child_only"):
                parent_delta = scale if control == "parent_only" else 0.
                child_delta = scale if control == "child_only" else 0.
                # Source-dependent initial geometry offsets must not count as
                # dynamic feedback. A response outside the child aggregate is
                # assessed on the parent ANNULUS and parent-detail modes.
                offset = 100. if control != "baseline" else 0.
                rows = [{"time": mark,
                         "child_occupation": [.5 + (-.1 + .015 * parent_delta) * mark] * 2,
                         "parent_detail_occupation": [.25 + (-.02 + .025 * child_delta) * mark] * 4,
                         "windows": {"child": {"probability": .8 + (-.1 + .03 * parent_delta) * mark},
                                     "left_annulus": {"probability": .4 + (-.01 + .02 * child_delta) * mark},
                                     "right_annulus": {"probability": .5 + (-.01 + .02 * child_delta) * mark}},
                         "metrics": {"child_r_proper_mean": 3.2 + offset + .02 * parent_delta * mark,
                                     "parent_annulus_r_proper_mean": 3.4 + offset + .05 * child_delta * mark}}
                        for mark in times]
                results[f"nf{nf}_{control}_dt{cap:g}"] = {"rows": rows}
    assessment = driver.compare_cases(results)
    for name, effect in (("parent_to_child_state", .03),
                         ("child_to_parent_state", .04),
                         ("parent_to_child_modal_state", .015),
                         ("child_to_parent_detail_modal_state", .025),
                         ("parent_to_child_geometry", .02),
                         ("child_to_parent_annulus_geometry", .05)):
        assert assessment[name]["effect_max"] == pytest.approx(effect, abs=2e-13)
        assert assessment[name]["time_fraction"] == pytest.approx(.002, abs=2e-11)
        assert assessment[name]["space_fraction"] == pytest.approx(.005, abs=2e-11)
        assert assessment[name]["resolved_at_one_percent"] is True
        assert assessment[name]["matched_proper_time_claimed"] is False
    for nf in (coarse, fine):
        for cap in (.001, .0005):
            result = results[f"nf{nf}_parent_only_dt{cap:g}"]
            baseline = results[f"nf{nf}_baseline_dt{cap:g}"]
            for row, base in zip(result["rows"], baseline["rows"]):
                row["windows"]["child"]["probability"] = base["windows"]["child"]["probability"]
    unresolved = driver.compare_cases(results)["parent_to_child_state"]
    assert unresolved["effect_max"] == 0.
    assert unresolved["time_fraction"] is None and unresolved["space_fraction"] is None
    assert unresolved["resolved_at_one_percent"] is False


@pytest.mark.parametrize("nf", [128, 256])
def test_successor_initial_preparation_reports_actual_returned_source_and_correction(nf):
    pair = nested.build_pair(nf)
    state, report = nested.initial_state(pair)
    actual = nested.reconstruct_state(pair, state)
    fine = galerkin.prolong_state(pair.grid, actual)
    system = galerkin.active_fine_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    rho = source["force_L"] / pair.grid.dx_q
    residual, _full = galerkin.projected_radius_operator(pair.grid, actual.r, rho)
    jacobian = galerkin._projected_radius_jacobian(pair.grid, actual.r)
    delta = np.linalg.solve(jacobian, -residual)
    correction = report["post_canonical_map_radius_correction"]
    assert correction["newton_correction_inf"] == pytest.approx(np.max(abs(delta)), rel=1e-10, abs=1e-17)
    assert correction["newton_correction_l2"] == pytest.approx(
        np.sqrt(pair.grid.dx_g) * np.linalg.norm(delta), rel=1e-10, abs=1e-17)
    assert report["newton_correction_inf"] <= report["radius_correction_floor"]
    assert report["reordered_correction_inf"] <= report["radius_correction_floor"]
    assert report["radius_correction_floor"] == pytest.approx(
        report["radius_ulp_floor"] + report["arithmetic_correction_indicator"])
    assert report["arithmetic_indicator_certified"] is False
    assert report["continuum_initial_state_certified"] is False
    assert report["bracket_positive"] and report["converged"]
    assert np.min(fine.r) > 0 and np.min(fine.Q) > 0
    assert np.array_equal(state.phi0, pair.source_phi0)
    assert np.array_equal(state.phi1, pair.source_phi1)
    full_constraint = coupling.hamilton_constraint(system, fine) + rho
    assert report["diagnostics"]["full_hamilton_max"] == pytest.approx(np.max(abs(full_constraint)), abs=1e-12)
    shift = coupling.shift_constraint(system, fine) + source["force_beta"] / pair.grid.dx_q
    assert report["shift_residual_mean"] == pytest.approx(np.mean(shift), abs=1e-15)
    assert report["source_current_mean"] == pytest.approx(
        np.mean(source["force_beta"] / pair.grid.dx_q), abs=1e-15)
    # Check the retained Jacobian against a separate directional residual,
    # including both broad-parent and localized-child radius directions.
    direction = pair.geometry_map[:, pair.geometry_parent_indices[0]] + \
                pair.geometry_map[:, pair.geometry_child_indices[0]]
    step = 1e-5
    plus = galerkin.projected_radius_operator(pair.grid, actual.r + step * direction, rho)[0]
    minus = galerkin.projected_radius_operator(pair.grid, actual.r - step * direction, rho)[0]
    numerical = (plus - minus) / (2 * step)
    expected = jacobian @ direction
    assert np.linalg.norm(numerical - expected) / np.linalg.norm(expected) < 2e-6


def _density(j):
    A, f, Z, CF = .04, -8 * np.pi * -.001 / 3, -24 * np.pi * .04, .01
    Fr = -8 * np.pi * A * j["r"]
    Fx = Fr * j["r_x"] + f * j["chi_x"]
    D = (j["Q_t"] - j["beta_x"] * j["Q"] - j["beta"] * j["Q_x"]) / j["L"]
    u, v = j["r_t"] - j["beta"] * j["r_x"], j["chi_t"] - j["beta"] * j["chi_x"]
    V = 8 * np.pi * A * j["r"] ** 2 - f / 2 * (4 * j["chi"] + j["chi"] ** 2) - 2 * np.pi * CF * 16
    return (2 * D * (Fr * u + f * v) - 2 * Fx * j["L_x"] / j["Q"] +
            Z * (j["Q"] / j["L"] * u ** 2 - j["L"] / j["Q"] * j["r_x"] ** 2) +
            j["L"] * j["Q"] * V)


def _dirac_action(j):
    sigma1 = np.array([[0., 1.], [1., 0.]])
    sigma2 = np.array([[0., -1j], [1j, 0.]])
    c = j["L"] / j["Q"]
    cx = (j["L_x"] * j["Q"] - j["L"] * j["Q_x"]) / j["Q"] ** 2
    return (-1j * sigma2 @ (c * j["phi_x"] + cx / 2 * j["phi"]) +
            j["L"] * sigma1 @ j["phi"] +
            1j * (j["beta"] * j["phi_x"] + j["beta_x"] / 2 * j["phi"]))


@pytest.mark.parametrize("scale", [2., 4.])
def test_affine_pullback_action_spinor_generator_symplectic_form_and_clock(scale):
    jets = dict(L=1.2, L_x=.17, Q=.7, Q_t=.2, Q_x=-.13,
                r=1.3, r_t=.14, r_x=.21, chi=.09, chi_t=-.16, chi_x=.08,
                beta=.12, beta_x=-.07, p_Q=.3, p_r=.4, p_chi=.5,
                phi=np.array([.2 + .3j, -.4 + .7j]),
                phi_x=np.array([.7 - .2j, .1 + .5j]))
    mapped = nested.functional_pullback(jets, scale)
    assert _density(mapped) == pytest.approx(scale ** 2 * _density(jets), abs=2e-13)
    np.testing.assert_allclose(_dirac_action(mapped), scale ** 1.5 * _dirac_action(jets),
                               rtol=2e-14, atol=2e-14)
    old_form = sum(jets[p] * jets[v] for p, v in
                   (("p_Q", "Q_t"), ("p_r", "r_t"), ("p_chi", "chi_t")))
    new_form = sum(mapped[p] * mapped[v] for p, v in
                   (("p_Q", "Q_t"), ("p_r", "r_t"), ("p_chi", "chi_t")))
    assert new_form == pytest.approx(scale ** 2 * old_form)
    assert np.vdot(mapped["phi"], mapped["phi"]).real / scale == pytest.approx(
        np.vdot(jets["phi"], jets["phi"]).real)
    assert mapped["r"] * mapped["L"] / scale == pytest.approx(jets["r"] * jets["L"])
    wrong = dict(mapped, phi=jets["phi"], phi_x=scale * jets["phi_x"])
    assert np.linalg.norm(_dirac_action(wrong) - scale ** 1.5 * _dirac_action(jets)) > .1
    wrong = dict(mapped, r=scale * jets["r"])
    assert abs(_density(wrong) - scale ** 2 * _density(jets)) > .1
