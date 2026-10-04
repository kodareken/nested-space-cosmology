"""Manufactured checks for the coupled RN PG state. No production trajectory."""
import math
import os

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import numpy as np
import pytest

from recursive_horizons import nsc_rn_observables as observables
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_source as source


def _family():
    return 1.04, 1.0


def test_sbp_closure_and_polynomial_jets():
    pack = source.radial_sbp(12, 0.0, 2.2)
    assert pack["closure_residual"] < 1e-12
    assert pack["periodic"] is False
    radius = np.linspace(1.0, 3.4, 13)
    sample = radius ** 4
    first, second = pg._differentiate(sample, radius[1] - radius[0])
    assert np.max(np.abs(first - 4.0 * radius ** 3)) < 1e-9
    assert np.max(np.abs(second - 12.0 * radius ** 2)) < 1e-8


def test_sourcefree_rn_is_stationary_with_the_pg_shift():
    mass, r_m = _family()
    grid, state, rates, lapse_target, shift_target = pg.sourcefree_rates(
        mass, r_m=r_m, points=33, killing_frequency=1.5,
    )
    assert grid["r_out"] == pytest.approx(32.0 * r_m)
    assert grid["r_m"] != pytest.approx(grid["r_minus"])
    assert grid["child_interval"][1] <= grid["parent_interval"][0]
    assert grid["same_areal_radius"] is False
    assert np.max(np.abs(rates["lapse"] - lapse_target)) < 1e-10
    assert np.max(np.abs(rates["shift"] - shift_target)) < 1e-10
    assert np.max(np.abs(rates["mass"] - mass)) < 1e-10
    assert np.max(np.abs(rates["phi_rate"])) < 1e-10
    assert abs(rates["inner_mass_rate"]) < 1e-10
    assert np.max(np.abs(rates["lapse_rate"])) < 1e-10
    assert np.max(np.abs(rates["shift_rate"])) < 1e-10
    assert np.max(np.abs(rates["balance_residual"])) < 1e-10
    assert rates["admission"]["excision_pure_outflow"] is True
    assert rates["admission"]["outer_incoming_only"] is True
    assert rates["electric_current"] == 0.0
    assert rates["fitted_restoring_force"] is False
    assert rates["lapse_held_at_one"] is False
    assert state["occupations_fixed"] is True
    assert state["multiplicity"] == 4
    outgoing = -rates["shift"] + rates["lapse"]
    ingoing = -rates["shift"] - rates["lapse"]
    assert outgoing[0] < 0.0 and ingoing[0] < 0.0
    assert outgoing[-1] > 0.0 and ingoing[-1] < 0.0
    corrected = observables.charged_mass(
        grid["radius"], rates["lapse"], rates["shift"], grid["magnetic_r2"],
    )
    assert np.max(np.abs(corrected - mass)) < 1e-8


def test_reference_curvature_and_nodal_jets():
    mass, r_m = _family()
    grid, _state, rates, _lapse, shift_target = pg.sourcefree_rates(
        mass, r_m=r_m, points=257,
    )
    invariants = grid["model"].invariants(grid["radius"])
    assert np.max(np.abs(invariants["R4"])) == pytest.approx(0.0)
    jets = grid["model"].pg_metric_jets(grid["radius"])
    assert np.max(np.abs(rates["jets"]["shift_r"] - jets["beta_r"])) <= 1e-2
    assert np.max(np.abs(rates["jets"]["shift_rr"] - jets["beta_rr"])) <= 1e-2
    assert np.max(np.abs(rates["shift"] - shift_target)) < 1e-10
    assert rates["jets"]["jet_provenance"] == pg.JET_PROVENANCE


def test_manufactured_mass_slope_does_not_apply_the_shell_twice():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33)
    state = pg.make_state(grid, inner_mass=mass)
    bare = 0.01
    zeros = np.zeros(grid["points"])
    densities = {
        "eps_N": np.full(grid["points"], bare),
        "eps_beta": zeros,
        "angular": zeros,
        "electric_current": zeros,
    }
    rates = pg.stage_rates(state, grid, densities=densities)
    slope = pg.NEWTON_G * bare
    expected = mass + slope * (grid["radius"] - grid["radius"][0])
    assert source.shell_multiplicity(1) == 4
    with pytest.raises(ValueError):
        source.shell_multiplicity(2)
    assert np.max(np.abs(rates["mass"] - expected)) < 1e-9
    assert np.max(np.abs(rates["lapse"] - 1.0)) < 1e-10
    derivative = grid["derivative"] @ rates["mass"]
    assert np.max(np.abs(derivative - slope)) < 1e-9


def test_shell_force_scales_once():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33, r_out=8.0)
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    wave = 2.0 * math.pi / (grid["radius"][-1] - grid["radius"][0])
    phi[0, :, 0] = np.exp(1j * wave * grid["radius"])
    occupation = np.array([0.25])
    state = pg.make_state(grid, phi, occupations=occupation, killing_frequency=wave)
    channel = source.mode_quadratic_form(
        phi[:, :, 0], np.ones(grid["points"]), np.zeros(grid["points"]),
        np.ones(grid["points"]), grid["radius"], grid["weights"], grid["derivative"],
    )
    evaluated = source.source_evaluate(
        phi, occupation, np.ones(grid["points"]), np.zeros(grid["points"]),
        np.ones(grid["points"]), grid["radius"], grid["weights"], grid["derivative"],
    )
    assert evaluated["multiplicity"] == 4
    assert evaluated["killing_energy"] == pytest.approx(4.0 * 0.25 * channel["energy"], rel=1e-10)
    assert state["multiplicity"] == 4


def test_momentum_density_moves_the_lapse_and_keeps_outer_normalization():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33, r_out=6.0)
    state = pg.make_state(grid, inner_mass=mass)
    zeros = np.zeros(grid["points"])
    densities = {
        "eps_N": zeros,
        "eps_beta": np.full(grid["points"], -1e-3),
        "angular": zeros,
        "electric_current": zeros,
    }
    rates = pg.stage_rates(state, grid, densities=densities)
    assert rates["lapse"][-1] == pytest.approx(1.0)
    assert np.max(np.abs(rates["lapse"] - 1.0)) > 1e-4
    assert rates["lapse_rate"][-1] == pytest.approx(0.0, abs=1e-12)
    assert np.all(rates["lapse"] > 0.0)


def test_chart_refuses_an_inadmissible_boundary():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33)
    with pytest.raises(pg.ChartAdmissionError, match="pure outflow"):
        pg.admit_chart(grid, np.ones(grid["points"]), np.full(grid["points"], 0.2))
    with pytest.raises(ValueError, match="Killing frequency"):
        pg.make_state(grid, killing_frequency=-1.0)


def test_occupations_stay_fixed_and_are_not_the_gram():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=17, r_out=8.0)
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    phi[0, :, 0] = 0.3
    occupations = np.array([0.2])
    state = pg.make_state(grid, phi, occupations=occupations, killing_frequency=0.4)
    gram = source.weighted_gram(state["phi"], grid["weights"])
    assert state["occupations"][0] == pytest.approx(0.2)
    assert abs(gram[0, 0] - 0.2) > 1e-3
    updated = pg.rk4_step(state, grid, 0.2 * pg.stage_rates(state, grid)["cfl_dt"])
    assert np.array_equal(updated["occupations"], state["occupations"])
    assert updated["clocks"]["t"] > 0.0
    assert updated["electric_current"] == 0.0


def test_plane_wave_energy_and_angular_coupling():
    radius = np.linspace(4.0, 8.0, 129)
    pack = source.radial_sbp(radius.size, float(radius[0]), float(radius[-1]))
    grid = {
        "radius": radius,
        "weights": pack["weights"],
        "derivative": pack["derivative"],
        "spacing": pack["step"],
        "points": radius.size,
        "charge": 0.2,
        "mass": 1.0,
    }
    wave = 2.0 * np.pi / (radius[-1] - radius[0])
    phi = np.zeros((2, radius.size, 1), dtype=complex)
    phi[0, :, 0] = np.exp(1j * wave * radius)
    occupation = np.array([0.5])
    evaluated = source.source_evaluate(
        phi, occupation, np.ones(radius.size), np.full(radius.size, 0.25),
        np.ones(radius.size), radius, pack["weights"], pack["derivative"],
    )
    omega = (1.0 - 0.25) * wave
    energy = evaluated["coordinate_energy"]
    expected = 4.0 * 0.5 * omega
    interior = slice(6, -6)
    assert np.max(np.abs(energy[interior] - expected)) <= 1e-2
    lapse = np.ones(radius.size)
    shift = np.full(radius.size, 0.25)
    sample = np.zeros_like(phi)
    sample[0, :, 0] = 1.0
    rates, _debit = pg.dirac_rates(grid, sample, lapse, shift, kappa=1)
    assert rates[1, radius.size // 2, 0] == pytest.approx(1.0 / radius[radius.size // 2], abs=1e-9)


def _principal_rates(grid, phi, lapse, shift):
    plus = phi[0, :, 0]
    minus = phi[1, :, 0]
    plus_rate = -pg._split(grid, lapse, plus) + pg._split(grid, shift, plus)
    minus_rate = pg._split(grid, lapse, minus) + pg._split(grid, shift, minus)
    return plus_rate, minus_rate


def test_ingoing_angular_sign_is_the_hermitian_generator():
    from recursive_horizons.nsc_covariant_operator import SIGMA1, SIGMA2
    basis = source.characteristic_transform()
    assert np.allclose(basis.conj().T @ SIGMA1 @ basis, SIGMA2, atol=1e-12)
    radius = np.linspace(4.0, 8.0, 65)
    pack = source.radial_sbp(radius.size, float(radius[0]), float(radius[-1]))
    grid = {
        "radius": radius, "weights": pack["weights"], "derivative": pack["derivative"],
        "spacing": pack["step"], "points": radius.size,
    }
    phi = np.zeros((2, radius.size, 1), dtype=complex)
    phi[0, :, 0] = np.exp(1j * radius)
    phi[1, :, 0] = 0.3 * np.exp(-0.5j * radius)
    lapse = np.ones(radius.size)
    shift = np.zeros(radius.size)
    rates, _debit = pg.dirac_rates(grid, phi, lapse, shift, kappa=1)
    plus_transport, minus_transport = _principal_rates(grid, phi, lapse, shift)
    alpha = 1.0 / radius
    assert np.max(np.abs(rates[0, :, 0] - (plus_transport - alpha * phi[1, :, 0]))) < 1e-12
    interior = slice(2, -2)
    assert np.max(np.abs(rates[1, interior, 0] - (minus_transport[interior] + alpha[interior] * phi[0, interior, 0]))) < 1e-12
    angular_rate_plus = -alpha * phi[1, :, 0]
    angular_rate_minus = alpha * phi[0, :, 0]
    overlap = np.sum(pack["weights"] * (
        np.conj(phi[0, :, 0]) * angular_rate_plus + np.conj(phi[1, :, 0]) * angular_rate_minus
    ))
    assert abs(overlap.real) < 1e-12


def test_prepared_stationary_mode_tracks_the_killing_rate():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=81, r_out=8.0)
    energy = 1.0
    solved = source.scattering_mode(grid["radius"], energy, grid["model"])
    phi = solved["mode"][:, :, None]
    lapse = np.asarray(grid["model"].vacuum_lapse(grid["radius"]), dtype=float)
    shift = np.asarray(grid["model"].beta_rn(grid["radius"]), dtype=float)
    rates, _debit = pg.dirac_rates(grid, phi, lapse, shift, kappa=1)
    expected = -1j * energy * solved["mode"]
    interior = slice(6, -6)
    scale = np.max(np.abs(expected[:, interior]))
    assert np.max(np.abs(rates[:, interior, 0] - expected[:, interior])) / scale < 0.05


def test_integer_flux_matches_the_shared_action_and_is_not_fitted():
    mass, r_m = _family()
    model = pg.action_model(mass, r_m)
    audited = source.action_reference(mass)
    assert model.flux == 1.0
    assert model.flux != pytest.approx(0.28209479177387814)
    assert model.G_N == pytest.approx(audited.G_N)
    assert model.G_N != pytest.approx(1.0)
    assert model.A == pytest.approx(audited.A)
    assert model.magnetic_r2 == pytest.approx(audited.magnetic_r2)
    assert math.sqrt(model.magnetic_r2) == pytest.approx(r_m)
    assert pg.ACTION_FLUX == 1.0
    with pytest.raises(pg.ChartAdmissionError, match="fit"):
        pg.action_model(mass, 2.0)


def test_sourcefree_rates_match_the_independent_reference():
    mass, r_m = _family()
    grid, _state, rates, _lapse, _shift = pg.sourcefree_rates(
        mass, r_m=r_m, points=33, killing_frequency=1.5,
    )
    sample = source_reference_static(grid)
    jets = sample["metric_jets"]
    assert np.max(np.abs(rates["lapse"] - jets["N"])) < 1e-10
    assert np.max(np.abs(rates["shift"] - jets["beta"])) < 1e-10
    assert np.max(np.abs(rates["phi_rate"])) < 1e-10
    assert np.max(np.abs(rates["lapse_rate"])) < 1e-10
    assert np.max(np.abs(rates["shift_rate"])) < 1e-10
    assert np.max(np.abs(rates["mass_rate"])) < 1e-10
    assert np.max(np.abs(sample["invariants"]["R4"])) == pytest.approx(0.0)
    assert rates["dynamic_contract"]["mass_definition_residual"] < 1e-10
    assert rates["dynamic_contract"]["static_R4_used"] is False
    assert rates["time_derivative_constraints"]["N_outer"] == pytest.approx(1.0)
    assert rates["time_derivative_constraints"]["N_t_outer"] == pytest.approx(0.0, abs=1e-12)
    assert rates["time_derivative_constraints"]["N_t_outer_fixed_by_normalization"] is True
    assert rates["jets"]["spatial_only"] is True
    assert rates["jets"]["dynamic_R4_not_claimed"] is True
    assert rates["time_jets"]["static_R4_used"] is False
    observer = pg.geometry_for_observer(rates, grid, stationary=True)
    assert observer["dynamic_R4_from_these_jets"] is False
    assert observer["time_derivatives_included"] is False
    directional = pg.directional_metric_check(
        pg.make_state(grid, inner_mass=mass, killing_frequency=1.5), grid, epsilon=1.0e-6,
    )
    assert directional["static_R4_used"] is False
    assert directional["lapse_residual"] < 1e-8
    assert directional["shift_residual"] < 1e-8
    assert directional["mass_residual"] < 1e-8


def source_reference_static(grid):
    from recursive_horizons import nsc_rn_reference as reference
    return reference.sourcefree_static_rn(
        grid["radius"], grid["mass"], A=grid["A"], C_F=grid["C_F"], flux=grid["flux"],
    )


def test_weak_coherent_packet_does_not_freeze_the_lapse_and_inner_flux_is_positive():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33, r_out=8.0)
    model = source.action_reference(mass)
    radius = grid["radius"]
    packet = _exterior_column(model)
    phi = np.zeros((2, radius.size, 1), dtype=complex)
    for component in range(2):
        phi[component, :, 0] = _sample_on_grid(
            packet["radius"], packet["phi"][component], radius,
        )
    state = pg.make_state(
        grid, phi, occupations=np.array([packet["nu"]]),
        killing_frequency=packet["omega"], inner_mass=mass,
    )
    assert state["killing_frequency_metadata"]["metadata_only"] is True
    assert state["killing_frequency_metadata"]["proves_exterior_killing_frequency"] is False
    assert state["multiplicity_in_covariance"] is False
    assert state["covariance_divided_by_rank"] is False
    assert state["multiplicity_over_rank_used"] is False
    assert state["covariance"].shape == (1, 1)
    rates = pg.stage_rates(state, grid)
    assert rates["shell_factor"] == 1.0
    assert rates["G_N"] == pytest.approx(model.G_N)
    assert rates["G_N"] == pytest.approx(packet["G_N"])
    assert rates["dynamic_contract"]["intrinsic_stress_cached_at_unit_lapse"] is True
    assert rates["dynamic_contract"]["stress_iteration"] is False
    assert rates["dynamic_contract"]["q_adm"] == 1.0
    assert rates["lapse_held_at_one"] is False
    assert np.max(np.abs(rates["lapse"] - 1.0)) > 1e-8
    assert rates["dynamic_contract"]["second_v_fn_term_max"] > 0.0
    assert "probability_stocks" in state
    assert set(state["stocks"]) == {"excision_flux", "outer_flux", "sat_debit"}
    control = _ingoing_control(grid)
    control_rates = pg.stage_rates(control, grid)
    assert control_rates["excision_flux_rate"] > 0.0
    assert control_rates["admission"]["components_clamped"] is False
    assert control_rates["admission"]["excision_pure_outflow"] is True
    assert state["probability_stocks"]["excision_outflow"] == 0.0
    assert control_rates["sat_debit_rate"] != control_rates["excision_flux_rate"]


def _exterior_column(model):
    r_m = source.magnetic_radius(model)
    request = source.packet_request(r_m)
    left = float(model.horizons["r_plus"]) + 0.25 * r_m
    right = request["center"] + 4.0 * request["width"]
    pack = source.radial_sbp(33, left, right)
    chart = source.vacuum_chart(model, pack["grid"])
    solved = source.scattering_mode(pack["grid"], request["omega"], model)
    mode = solved["mode"]
    norm = np.sqrt(np.real(np.sum(pack["weights"] * np.sum(np.abs(mode) ** 2, axis=0))))
    mode = mode / norm
    measured = source.mode_quadratic_form(
        mode, chart["lapse"], chart["shift"], chart["q_adm"],
        chart["radius"], pack["weights"], pack["derivative"],
    )
    filling, _actual = source.select_packet_filling(
        measured["energy"], model.mass, 1.0e-3, model.G_N,
    )
    return {
        "radius": pack["grid"],
        "phi": mode,
        "nu": filling,
        "omega": request["omega"],
        "G_N": float(model.G_N),
    }


def _sample_on_grid(source_radius, values, target_radius):
    sampled = np.zeros(target_radius.shape, dtype=complex)
    mask = (target_radius >= source_radius[0]) & (target_radius <= source_radius[-1])
    if not np.any(mask):
        return sampled
    sampled[mask] = (
        np.interp(target_radius[mask], source_radius, np.real(values))
        + 1j * np.interp(target_radius[mask], source_radius, np.imag(values))
    )
    return sampled


def _ingoing_control(grid):
    radius = grid["radius"]
    envelope = np.exp(-0.5 * ((radius - radius[0]) / grid["spacing"]) ** 2)
    phi = np.zeros((2, radius.size, 1), dtype=complex)
    phi[1, :, 0] = envelope * np.exp(-1j * radius)
    return pg.make_state(
        grid, phi, occupations=np.array([1.0e-4]), killing_frequency=1.0, inner_mass=grid["mass"],
    )


def test_weighted_columns_keep_the_physical_car_unscaled():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=17, r_out=6.0)
    phi = np.zeros((2, grid["points"], 2), dtype=complex)
    phi[0, :, 0] = 0.2
    phi[1, :, 1] = 0.3
    coefficients = np.array([[0.2, 0.05], [0.05, 0.1]], dtype=complex)
    state = pg.make_state(grid, phi, occupations=coefficients, killing_frequency=0.7)
    gram = source.weighted_gram(phi, grid["weights"])
    car = source.physical_car(gram, coefficients)
    assert state["occupation_matrix"].shape == (2, 2)
    assert np.max(np.abs(state["covariance"] - car)) < 1e-10
    assert state["covariance_is_physical_car"] is True
    assert state["multiplicity_in_covariance"] is False
    assert state["multiplicity"] == 4
    rates = pg.stage_rates(state, grid)
    assert rates["shell_factor"] == 1.0
    assert rates["dynamic_contract"]["mass_slope_residual"] < 1e-8


def test_vacuum_rk4_step_does_not_move_the_inner_stock():
    mass, r_m = _family()
    grid = pg.build_grid(mass, r_m=r_m, points=33)
    state = pg.make_state(grid, inner_mass=mass, killing_frequency=1.0)
    cfl = pg.stage_rates(state, grid)["cfl_dt"]
    updated = pg.rk4_step(state, grid, 0.5 * cfl)
    assert updated["inner_mass"] == pytest.approx(mass, abs=1e-12)
    assert np.max(np.abs(updated["phi"])) == pytest.approx(0.0)
    assert updated["stocks"]["sat_debit"] == pytest.approx(0.0)
    with pytest.raises(pg.ChartAdmissionError, match="CFL"):
        pg.rk4_step(state, grid, 2.0 * cfl)


def _column(grid, scale):
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    center = float(grid["radius"][0] + 0.35 * (grid["radius"][-1] - grid["radius"][0]))
    width = 0.15 * (grid["radius"][-1] - grid["radius"][0])
    envelope = np.exp(-0.5 * ((grid["radius"] - center) / width) ** 2)
    phi[0, :, 0] = scale * envelope
    phi[1, :, 0] = 0.25 * scale * envelope * np.exp(-1j * grid["radius"])
    return phi


def _compare_layouts(points, r_out, scale, steps=1):
    from scipy import sparse
    dense = pg.build_grid(1.04, r_m=1.0, points=points, r_out=r_out, layout="dense")
    csr = pg.build_grid(1.04, r_m=1.0, points=points, r_out=r_out, layout="csr")
    assert csr["derivative_layout"] == "csr" and csr["production_operator"] is True
    assert csr["derivative_family"] == "diagonal_norm_sbp_d21"
    assert sparse.isspmatrix_csr(csr["derivative"])
    assert isinstance(dense["derivative"], np.ndarray)
    assert np.array_equal(csr["radius"], dense["radius"])
    assert np.max(np.abs(csr["derivative"].toarray() - dense["derivative"])) == 0.0
    phi = _column(dense, scale)
    occupation = np.array([1.0e-4])
    left = pg.make_state(dense, phi, occupations=occupation, inner_mass=dense["mass"])
    right = pg.make_state(csr, phi, occupations=occupation, inner_mass=csr["mass"])
    left_rates = pg.stage_rates(left, dense)
    right_rates = pg.stage_rates(right, csr)
    assert np.max(np.abs(left_rates["phi_rate"] - right_rates["phi_rate"])) == 0.0
    assert np.max(np.abs(left_rates["lapse"] - right_rates["lapse"])) == 0.0
    assert np.max(np.abs(left_rates["shift"] - right_rates["shift"])) == 0.0
    assert np.max(np.abs(left_rates["mass_slope"] - right_rates["mass_slope"])) == 0.0
    dt = min(0.5 * left_rates["cfl_dt"], 0.001)
    moved_left = left
    moved_right = right
    for _ in range(steps):
        moved_left = pg.rk4_step(moved_left, dense, dt)
        moved_right = pg.rk4_step(moved_right, csr, dt)
    assert np.max(np.abs(moved_left["phi"] - moved_right["phi"])) == 0.0
    assert moved_left["inner_mass"] == pytest.approx(moved_right["inner_mass"], abs=0.0)
    assert moved_left["stocks"]["excision_flux"] == pytest.approx(
        moved_right["stocks"]["excision_flux"], abs=0.0,
    )
    return dense, csr, dt


def test_csr_matches_dense_field_source_metric_and_rk4():
    production = pg.build_grid(1.04, r_m=1.0, points=33, r_out=16.0)
    assert production["derivative_layout"] == "csr"
    _compare_layouts(33, 16.0, 1.0e-4)


def test_saved_state_agrees_and_current_car_is_not_the_cached_field():
    from pathlib import Path
    saved = Path(
        "/Users/admin/Documents/BlackHoles-Infinity/lab/results/development/"
        "nsc-rn-neutral-continuation-v1/t20.0000.npz"
    )
    if not saved.is_file():
        pytest.skip("saved T20 state is not on this machine")
    stamp = saved.stat().st_mtime_ns
    size = saved.stat().st_size
    with np.load(saved) as payload:
        radius = np.array(payload["radius"])
        phi = np.array(payload["coupled_phi"])
        occupation = np.array(payload["occupations"])
        inner = float(payload["coupled_inner_mass"][0])
    assert saved.stat().st_mtime_ns == stamp and saved.stat().st_size == size
    dense = pg.build_grid(1.04, r_m=1.0, points=radius.size, r_out=float(radius[-1]), layout="dense")
    csr = pg.build_grid(1.04, r_m=1.0, points=radius.size, r_out=float(radius[-1]), layout="csr")
    assert np.array_equal(dense["radius"], radius)
    state = pg.make_state(csr, phi, occupations=occupation, inner_mass=inner)
    assert np.array_equal(state["occupations"], occupation)
    cached = np.array(state["covariance"])
    loaded = pg.current_occupation_car(state["phi"], csr["weights"], state["occupations"])
    assert loaded["psd"] is True
    assert loaded["min_eigenvalue"] == pytest.approx(float(occupation[0]) * float(np.real(
        source.weighted_gram(phi, csr["weights"])[0, 0]
    )))
    rates = pg.stage_rates(state, csr)
    dense_state = pg.make_state(dense, phi, occupations=occupation, inner_mass=inner)
    dense_rates = pg.stage_rates(dense_state, dense)
    assert np.max(np.abs(rates["phi_rate"] - dense_rates["phi_rate"])) == 0.0
    assert np.max(np.abs(rates["lapse"] - dense_rates["lapse"])) == 0.0
    dt = min(0.5 * rates["cfl_dt"], 0.001)
    moved = pg.rk4_step(state, csr, dt)
    assert np.array_equal(moved["occupations"], occupation)
    assert np.array_equal(moved["covariance"], cached)
    assert moved["covariance_describes"] == "load_time_phi"
    assert moved["covariance_is_current_eigenvalue"] is False
    current = pg.current_occupation_car(moved["phi"], csr["weights"], moved["occupations"])
    assert current["psd"] is True
    assert current["min_eigenvalue"] != pytest.approx(float(np.linalg.eigvalsh(
        0.5 * (cached + cached.conj().T)
    )[0]), abs=0.0)
    assert saved.stat().st_mtime_ns == stamp and saved.stat().st_size == size


def test_csr_full_pair_at_3921_is_at_least_one_and_a_half_times_dense():
    import resource
    import sys
    import time
    import derive_nsc_rn_parent as driver
    from scipy import sparse
    dense = pg.build_grid(1.04, r_m=1.0, points=3921, r_out=32.0, layout="dense")
    csr = pg.build_grid(1.04, r_m=1.0, points=3921, r_out=32.0, layout="csr")
    assert sparse.isspmatrix_csr(csr["derivative"])
    assert np.max(np.abs(csr["derivative"].toarray() - dense["derivative"])) == 0.0
    phi = _column(dense, 1.0e-4)
    occupation = np.array([1.0e-4])
    probe = pg.make_state(dense, phi, occupations=occupation, inner_mass=dense["mass"])
    rates = pg.stage_rates(probe, dense)
    metric = driver.freeze_sourced_metric(rates)
    dt = min(0.5 * rates["cfl_dt"], 0.001)

    def once(grid):
        coupled = pg.make_state(grid, phi, occupations=occupation, inner_mass=grid["mass"])
        fixed = pg.make_state(grid, phi, occupations=occupation, inner_mass=grid["mass"])
        started = time.process_time()
        marched = driver.march_pair(coupled, fixed, grid, metric, dt=dt, steps=1)
        elapsed = time.process_time() - started
        return marched, elapsed

    left, dense_cpu = once(dense)
    right, csr_cpu = once(csr)
    assert np.max(np.abs(left["coupled"]["phi"] - right["coupled"]["phi"])) == 0.0
    assert np.max(np.abs(left["fixed"]["phi"] - right["fixed"]["phi"])) == 0.0
    assert left["coupled"]["inner_mass"] == pytest.approx(right["coupled"]["inner_mass"], abs=0.0)
    assert csr_cpu > 0.0 and dense_cpu / csr_cpu >= 1.5
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    memory = rss if sys.platform == "darwin" else rss * 1024
    assert memory < 8 * 1024 ** 3
    print(
        f"N3921 dense_cpu={dense_cpu:.6f} csr_cpu={csr_cpu:.6f} "
        f"speedup={dense_cpu / csr_cpu:.3f} rss={memory}"
    )
