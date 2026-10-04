"""Coupled RN response tangent. No production march and no new result record."""
import hashlib
import math
import os
from pathlib import Path
import time

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import numpy as np
import pytest

from recursive_horizons import nsc_rn_observables as observables
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

import derive_nsc_rn_parent as driver


LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
ARCHIVE_T0 = LAB / "archive/rn-first-infall-pilot-2026-10-04/nsc-rn-pilot-rank9/t0.npz"
EVIDENCE = []
CPU_LIMIT_SECONDS = 60.0
FROZEN_SHA256 = {
    "lab/src/recursive_horizons/nsc_rn_pg.py":
        "d90fd73eb951cfb03582d459723993e707617be17bda5da0383e93c8f59285ec",
    "lab/src/recursive_horizons/nsc_rn_source.py":
        "eca21554c70d4fb0acdfef60977ff986221394927a550227b29517b0074f213b",
    "lab/src/recursive_horizons/nsc_rn_observables.py":
        "cd2e3b58129acd92a513f7e03b4db283ff15011293c2ed0ae12091764e6a3307",
    "lab/src/recursive_horizons/nsc_rn_reference.py":
        "4e61faaa18ba7f82461a6e2ac7e19590786a24d7cce79cf30241079a22bbbecc",
    "lab/scripts/derive_nsc_rn_parent.py":
        "8207f06e83a2dc8118a3dc366706a07c2f0dbc3d987406afbd7cc92272211cdf",
    "lab/docs/nsc-rn-parent-episode.md":
        "afc5b7729c31967308b8e7a83f6d72e42443979c96aba1e97af95baec0b10283",
}


def _record(name, **values):
    EVIDENCE.append({"check": name, **values})
    print(f"EVIDENCE {name} {values}")


def _max_abs(value):
    return float(np.max(np.abs(np.asarray(value))))


def _rel(error, signal):
    return error / max(signal, 1.0e-30)


@pytest.fixture(scope="session", autouse=True)
def _cpu_budget():
    started = time.process_time()
    yield
    used = time.process_time() - started
    print(f"EVIDENCE cpu_seconds {used}")
    assert used < CPU_LIMIT_SECONDS


def _grid(points=17, layout="dense", r_out=16.0):
    return pg.build_grid(1.04, r_m=1.0, points=points, r_out=r_out, layout=layout)


def _state(grid, filling=1.0e-4):
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    radius = grid["radius"]
    envelope = np.exp(-0.5 * ((radius - 12.0) / 2.0) ** 2)
    phi[0, :, 0] = 0.2 * envelope * np.exp(0.4j * (radius - 12.0))
    phi[1, :, 0] = (0.15 + 0.05j) * envelope
    return pg.make_state(
        grid, phi, occupations=np.array([filling]), inner_mass=float(grid["mass"]),
        killing_frequency=1.0,
    )


def _sourced_direction(state):
    radius = np.linspace(0.0, 1.0, state["phi"].shape[1])
    dphi = np.zeros_like(state["phi"])
    dphi[0, :, 0] = 0.05 * np.exp(-((radius - 0.45) / 0.15) ** 2)
    dphi[1, :, 0] = 0.03j * np.exp(-((radius - 0.55) / 0.18) ** 2)
    filling = float(state["occupations"][0])
    return response.state_direction(dphi, np.array([0.25 * filling]))


def _halving(error, half_error, *, floor):
    if error <= floor and half_error <= floor:
        return True
    return half_error < error / 3.0


def _within_signal(error, signal, *, relative=1.0e-5, absolute=1.0e-9, signal_floor=1.0e-6):
    """Relative check on a visible signal. Roundoff is an absolute bound."""
    if signal < signal_floor:
        return error <= absolute
    return error <= relative * signal


def test_frozen_producer_files_are_untouched():
    for relative, digest in FROZEN_SHA256.items():
        found = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert found == digest


def test_baseline_window_multiplicity_and_production_step():
    declared = response.baseline()
    assert declared["new_dynamical_relation"] is False
    assert declared["held_forecast"]["status"] == "not_run"
    assert declared["step_cap_over_rm"] == 0.001
    assert declared["production_half_step_over_rm"] == 0.0005
    assert declared["rejected_half_cfl_over_rm"] == 0.008742
    assert declared["reconstruction_weight"] == 1.0
    assert declared["multiplicity"] == 4
    assert declared["shell_insertions"] == 1
    assert source.shell_multiplicity(1) == 4
    with pytest.raises(ValueError):
        source.shell_multiplicity(0)
    with pytest.raises(observables.ObservableError):
        observables.shell_factor(4, 4)
    window = response.baseline_window()
    assert window["truncated"] is True
    assert window["energies_match_omega_rule"] is True
    assert window["plus_five_percent_width_truncated"] is True
    assert window["plus_five_percent_frequency_truncated"] is True
    assert np.allclose(window["energies"], window["width_shifted_energies"])
    assert np.max(np.abs(window["frequency_shifted_energies"] / window["energies"] - 1.05)) < 1e-12
    assert window["span"] == pytest.approx(1.2784926817352849, rel=1e-8)
    assert abs(window["ingoing_speed"]) == pytest.approx(1.4079079416840141, rel=1e-8)
    grid = _grid(points=17, layout="csr", r_out=32.0)
    state = _state(grid)
    rates = pg.stage_rates(state, grid)
    choice = response.production_step(rates["cfl_dt"], grid["r_m"])
    assert choice["limited_by"] == "step_cap"
    assert choice["dt"] == pytest.approx(0.001 * grid["r_m"])
    assert choice["production_half_step"] == pytest.approx(0.0005 * grid["r_m"])
    assert choice["dt"] < choice["cfl_half_step"]
    assert abs(choice["dt"] - response.REJECTED_HALF_CFL_OVER_RM) > 1e-3
    _record(
        "production_step",
        dt=choice["dt"],
        half_step=choice["production_half_step"],
        cfl_half=choice["cfl_half_step"],
        limited_by=choice["limited_by"],
        span=window["span"],
        ingoing_speed=window["ingoing_speed"],
    )


def test_filling_slope_normalization_and_car_initials():
    newton = pg.NEWTON_G
    low, actual_low = source.select_packet_filling(0.2, 1.04, 1.0e-3, newton)
    high, actual_high = source.select_packet_filling(0.2, 1.04, 1.0e-2, newton)
    slope = (high - low) / (1.0e-2 - 1.0e-3)
    assert slope == pytest.approx(low / 1.0e-3)
    assert actual_low == pytest.approx(1.0e-3)
    assert actual_high == pytest.approx(1.0e-2)
    derived = response.filling_derivative(
        0.2, 1.04, 1.0e-3, newton, d_epsilon=1.0e-3, d_packet_energy=0.0,
    )
    assert derived["filling"] == pytest.approx(low)
    assert derived["relative"] == pytest.approx(1.0)
    assert derived["multiplicity_applied_once"] is True
    grid = _grid(points=33, r_out=20.0)
    prepared = _synthetic_packet(grid)
    energy = response.energy_initial_tangent(prepared, d_epsilon=prepared["energy_target"])
    assert _max_abs(energy["direction"]["dphi"]) == 0.0
    assert energy["direction"]["d_inner_mass"] == 0.0
    assert energy["direction"]["d_clocks"] == response._zero_clocks()
    assert energy["direction"]["d_stocks"] == response._zero_stocks()
    assert energy["d_actual_energy_target"] == pytest.approx(prepared["energy_target"])
    assert energy["direction"]["d_occupations"][0] == pytest.approx(prepared["packet_filling"])
    assert energy["real_overlap"] == 0.0
    assert _max_abs(energy["d_car"] - energy["direction"]["d_occupations"][0] / prepared["packet_filling"] * prepared["car"]["physical"]) < 1e-12
    width = response.width_initial_tangent(prepared, d_width=0.05)
    assert abs(width["real_overlap"]) < 1e-9
    assert _max_abs(width["d_sample_energies"]) == 0.0
    assert width["modes_move"] is False
    assert width["d_actual_energy_target"] == 0.0
    relative = width["direction"]["d_occupations"][0] / prepared["packet_filling"]
    assert relative == pytest.approx(-width["d_packet_energy"] / prepared["packet_mode_energy"])
    assert width["multiplicity_in_car"] is False
    readout = response.current_car_tangent(
        prepared["phi"], width["direction"]["dphi"], grid["weights"],
        prepared["packet_filling"], width["direction"]["d_occupations"][0],
    )
    assert readout["multiplicity_included"] is False
    assert abs(readout["d_gram"]) < 1e-8
    state = pg.make_state(
        grid, prepared["phi"], occupations=prepared["occupations"],
        inner_mass=grid["mass"], killing_frequency=1.0,
    )
    other = pg.make_state(
        grid, prepared["phi"], occupations=prepared["occupations"],
        inner_mass=grid["mass"], killing_frequency=4.0,
    )
    left = pg.stage_rates(state, grid)
    right = pg.stage_rates(other, grid)
    assert _max_abs(left["phi_rate"] - right["phi_rate"]) == 0.0
    assert state["killing_frequency_metadata"]["proves_exterior_killing_frequency"] is False
    _record(
        "initials",
        filling_slope=slope,
        width_real_overlap=width["real_overlap"],
        width_relative_filling=relative,
        car_multiplicity=False,
    )


def _synthetic_packet(grid, rank=9):
    chart = source.vacuum_chart(grid["model"], grid["radius"])
    request = source.packet_request(grid["r_m"])
    speed = float(np.interp(request["center"], chart["radius"], chart["speed_minus"]))
    spectral_width = abs(speed) / request["width"]
    energies, truncated, _span = source._positive_frequency_window(
        request["omega"], spectral_width, rank,
    )
    assert truncated
    weights = grid["weights"]
    modes = []
    for index, energy in enumerate(energies):
        mode = np.zeros((2, grid["points"]), dtype=complex)
        center = request["center"] + 0.15 * (index - 0.5 * (rank - 1))
        envelope = np.exp(-0.5 * ((grid["radius"] - center) / request["width"]) ** 2)
        phase = np.exp(-1j * float(energy) * (grid["radius"] - request["center"]))
        mode[1] = envelope * phase
        mode[0] = 0.2 * envelope * phase * np.exp(0.2j * index)
        normalized, _norm = source._normalize_mode(mode, weights)
        aligned, _index = source._phase_align_ingoing(normalized, grid["radius"], request["center"])
        modes.append(aligned)
    spectral = np.stack(modes, axis=-1)
    amplitudes = source._gaussian_amplitudes(energies, request["omega"], spectral_width)
    superposition = spectral @ amplitudes
    norm = math.sqrt(float(np.real(np.sum(weights * np.sum(np.abs(superposition) ** 2, axis=0)))))
    packet = superposition / norm
    measured = source.mode_quadratic_form(
        packet, chart["lapse"], chart["shift"], chart["q_adm"], chart["radius"],
        weights, grid["derivative"],
    )
    if measured["energy"] <= 0.0:
        raise AssertionError("synthetic packet energy is not positive")
    filling, actual = source.select_packet_filling(
        measured["energy"], grid["mass"], 1.0e-3, grid["G_N"],
    )
    alpha = amplitudes / norm
    coefficient = filling * np.outer(alpha, np.conj(alpha))
    gram = source.weighted_gram(spectral, weights)
    car = source.physical_car(gram, coefficient)
    return {
        "phi": packet[:, :, None],
        "occupations": np.array([filling]),
        "packet_filling": filling,
        "packet_mode_energy": measured["energy"],
        "energy_target": 1.0e-3,
        "actual_energy_target": actual,
        "G_N": float(grid["G_N"]),
        "spectral_phi": spectral,
        "spectrum": {
            "requested": energies,
            "spectral_width": spectral_width,
            "ingoing_group_velocity": speed,
        },
        "radial_moments": {
            "requested_width": request["width"],
            "requested_omega": request["omega"],
            "requested_center": request["center"],
        },
        "background": {
            "weights": weights,
            "lapse": chart["lapse"],
            "shift": chart["shift"],
            "radius": chart["radius"],
            "derivative": grid["derivative"],
            "mass": float(grid["mass"]),
            "grid": chart["radius"],
        },
        "car": {"physical": car},
    }


def test_width_tangent_matches_owner_assembly_steps():
    grid = _grid(points=33, r_out=20.0)
    prepared = _synthetic_packet(grid)
    tangent = response.width_initial_tangent(prepared, d_width=1.0)
    request = prepared["radial_moments"]
    speed = prepared["spectrum"]["ingoing_group_velocity"]
    spectral = prepared["spectral_phi"]

    def assemble(width):
        spectral_width = abs(speed) / width
        energies, truncated, _span = source._positive_frequency_window(
            request["requested_omega"], spectral_width, spectral.shape[2],
        )
        assert truncated
        amplitudes = source._gaussian_amplitudes(energies, request["requested_omega"], spectral_width)
        superposition = spectral @ amplitudes
        norm = math.sqrt(float(np.real(np.sum(
            grid["weights"] * np.sum(np.abs(superposition) ** 2, axis=0)
        ))))
        packet = superposition / norm
        measured = source.mode_quadratic_form(
            packet, prepared["background"]["lapse"], prepared["background"]["shift"],
            np.ones(grid["points"]), grid["radius"], grid["weights"], grid["derivative"],
        )
        filling, actual = source.select_packet_filling(
            measured["energy"], grid["mass"], 1.0e-3, grid["G_N"],
        )
        return packet, filling, measured["energy"], energies, actual

    def difference(step):
        plus = assemble(request["requested_width"] + step)
        minus = assemble(request["requested_width"] - step)
        return tuple((left - right) / (2.0 * step) for left, right in zip(plus, minus))

    coarse = difference(1.0e-4)
    fine = difference(5.0e-5)
    phi_error = _max_abs(coarse[0] - tangent["direction"]["dphi"][:, :, 0])
    phi_error_half = _max_abs(fine[0] - tangent["direction"]["dphi"][:, :, 0])
    filling_error = abs(coarse[1] - tangent["direction"]["d_occupations"][0])
    energy_error = abs(coarse[2] - tangent["d_packet_energy"])
    energy_error_half = abs(fine[2] - tangent["d_packet_energy"])
    signal = _max_abs(tangent["direction"]["dphi"])
    assert _max_abs(coarse[3]) < 1e-9
    assert abs(coarse[4]) < 1e-9
    assert _rel(phi_error, signal) < 1e-6
    assert _halving(phi_error, phi_error_half, floor=1e-10)
    assert _rel(energy_error, abs(tangent["d_packet_energy"])) < 1e-5
    assert _halving(energy_error, energy_error_half, floor=1e-10)
    assert filling_error / max(abs(tangent["direction"]["d_occupations"][0]), 1e-30) < 1e-5
    _record(
        "width_fd",
        phi_rel=_rel(phi_error, signal),
        phi_rel_half=_rel(phi_error_half, signal),
        energy_rel=_rel(energy_error, abs(tangent["d_packet_energy"])),
        filling_abs=filling_error,
        signal=signal,
    )


def _force_bundle(grid):
    radius = grid["radius"]
    phase = np.linspace(0.0, math.pi, radius.size)
    force_n = 1.0e-3 * np.sin(phase)
    force_beta = -4.0e-4 * np.cos(phase)
    force_r = 2.0e-4 * np.sin(2.0 * phase)
    return force_n, force_beta, force_r


def test_reconstruction_and_stress_match_two_centered_steps():
    errors = []
    for layout in ("dense", "csr"):
        grid = _grid(points=17, layout=layout, r_out=16.0)
        state = _state(grid)
        direction = _sourced_direction(state)
        stress = response.stress_tangent(
            grid, state["phi"], state["occupations"], direction["dphi"], direction["d_occupations"],
        )

        def stress_of(scale):
            moved = response.displace_state(state, direction, scale)
            densities = pg.matter_densities(moved["phi"], grid, moved["occupations"])
            return densities["eps_N"], densities["eps_beta"], densities["force_r_unit"]

        def stress_fd(step):
            plus = stress_of(step)
            minus = stress_of(-step)
            return tuple((left - right) / (2.0 * step) for left, right in zip(plus, minus))

        coarse = stress_fd(1.0e-5)
        fine = stress_fd(5.0e-6)
        names = ("d_force_n", "d_force_beta", "d_force_r_unit")
        for name, coarse_value, fine_value in zip(names, coarse, fine):
            signal = _max_abs(stress[name])
            error = _max_abs(coarse_value - stress[name])
            half = _max_abs(fine_value - stress[name])
            assert _rel(error, signal) < 1e-6
            assert _halving(error, half, floor=1.0e-9)
            errors.append((layout, name, _rel(error, signal), _rel(half, signal)))
        force_n, force_beta, force_r = _force_bundle(grid)
        base_inner = float(grid["mass"])
        tangent = response.reconstruction_tangent(
            grid, base_inner, force_n, force_beta, force_r,
            0.0, np.cos(np.linspace(0.0, math.pi, grid["points"])),
            -0.5 * np.sin(np.linspace(0.0, math.pi, grid["points"])),
            0.2 * np.ones(grid["points"]),
        )
        d_inner = 0.0
        d_force_n = np.cos(np.linspace(0.0, math.pi, grid["points"]))
        d_force_beta = -0.5 * np.sin(np.linspace(0.0, math.pi, grid["points"]))
        d_force_r = 0.2 * np.ones(grid["points"])

        def geometry_of(step):
            geometry = pg.reconstruct_geometry(
                grid,
                base_inner + step * d_inner,
                force_n + step * d_force_n,
                force_beta + step * d_force_beta,
                -(force_r + step * d_force_r) * grid["radius"],
                weight=1.0,
                force_r_unit=force_r + step * d_force_r,
            )
            return geometry["mass"], geometry["lapse"], geometry["shift"], geometry["q_equation_mass_rate"]

        def geometry_fd(step):
            plus = geometry_of(step)
            minus = geometry_of(-step)
            return tuple((left - right) / (2.0 * step) for left, right in zip(plus, minus))

        coarse = geometry_fd(1.0e-5)
        fine = geometry_fd(5.0e-6)
        analytic = (
            tangent["d_mass"], tangent["d_lapse"], tangent["d_shift"], tangent["d_q_equation_mass_rate"],
        )
        labels = ("d_mass", "d_lapse", "d_shift", "d_flux")
        for label, coarse_value, fine_value, exact in zip(labels, coarse, fine, analytic):
            signal = _max_abs(exact)
            error = _max_abs(coarse_value - exact)
            half = _max_abs(fine_value - exact)
            assert _rel(error, signal) < 1e-6
            assert _halving(error, half, floor=1.0e-9)
            errors.append((layout, label, _rel(error, signal), _rel(half, signal)))
        assert tangent["dN_out"] == pytest.approx(0.0, abs=1e-14)
        assert tangent["reconstruction_weight"] == 1.0
        assert tangent["second_v_fn_coefficient"] == 2.0
        horizon = response.trapping_horizon_tangent(
            grid["radius"], tangent["lapse"], tangent["shift"], tangent["d_lapse"], tangent["d_shift"],
        )

        def horizon_radius(step):
            geometry = pg.reconstruct_geometry(
                grid,
                base_inner,
                force_n + step * d_force_n,
                force_beta + step * d_force_beta,
                -(force_r + step * d_force_r) * grid["radius"],
                weight=1.0,
                force_r_unit=force_r + step * d_force_r,
            )
            return observables.trapping_horizon(
                grid["radius"], geometry["lapse"], geometry["shift"],
            )["trapping_horizon"]

        horizon_fd = (horizon_radius(1.0e-5) - horizon_radius(-1.0e-5)) / 2.0e-5
        horizon_half = (horizon_radius(5.0e-6) - horizon_radius(-5.0e-6)) / 1.0e-5
        horizon_error = abs(horizon_fd - horizon["d_radius"])
        horizon_error_half = abs(horizon_half - horizon["d_radius"])
        assert _within_signal(horizon_error, abs(horizon["d_radius"]))
        assert _halving(horizon_error, horizon_error_half, floor=1.0e-9)
        errors.append((layout, "d_trapping_radius", _rel(horizon_error, abs(horizon["d_radius"])), abs(horizon["d_radius"])))
        radial_only = response.reconstruction_tangent(
            grid, base_inner, force_n, force_beta, force_r,
            0.0, np.zeros(grid["points"]), np.zeros(grid["points"]), np.ones(grid["points"]),
        )
        assert _max_abs(radial_only["d_mass"]) < 1e-10
        assert _max_abs(radial_only["d_q_equation_mass_rate"]) > 1e-8
    dense = _grid(points=17, layout="dense")
    csr = _grid(points=17, layout="csr")
    state = _state(dense)
    direction = _sourced_direction(state)
    # The same nodal samples on both operator layouts.
    state_csr = pg.make_state(
        csr, state["phi"], occupations=state["occupations"], inner_mass=state["inner_mass"],
        killing_frequency=state["killing_frequency"],
    )
    left = response.stage_tangent(state, dense, direction)
    right = response.stage_tangent(state_csr, csr, direction)
    assert _max_abs(left["dphi_rate"] - right["dphi_rate"]) < 1e-10
    assert _max_abs(left["d_mass"] - right["d_mass"]) < 1e-10
    assert _max_abs(left["d_lapse"] - right["d_lapse"]) < 1e-10
    _record("reconstruction_fd", rows=errors, dN_out=0.0, dense_csr_phi_rate_gap=_max_abs(left["dphi_rate"] - right["dphi_rate"]))


def _output_fields(updated):
    remaining = float(updated["probability_stocks"]["remaining"])
    excision = float(updated["probability_stocks"]["excision_outflow"])
    outer = float(updated["probability_stocks"]["outer_outflow"])
    sat = float(updated["stocks"]["sat_debit"])
    return {
        "phi": np.asarray(updated["phi"]),
        "inner_mass": float(updated["inner_mass"]),
        "t": float(updated["clocks"]["t"]),
        "killing": float(updated["clocks"]["killing"]),
        "normal": float(updated["clocks"]["normal"]),
        "excision_flux": float(updated["stocks"]["excision_flux"]),
        "outer_flux": float(updated["stocks"]["outer_flux"]),
        "sat_debit": sat,
        "remaining": remaining,
        "excision_outflow": excision,
        "outer_outflow": outer,
        "closure": remaining + excision + outer + 2.0 * sat,
    }


def _tangent_fields(bundle):
    probability = bundle["d_probability_stocks"]
    stocks = bundle["d_stocks"]
    clocks = bundle["d_clocks"]
    return {
        "phi": bundle["dphi"],
        "inner_mass": float(bundle["d_inner_mass"]),
        "t": float(clocks["t"]),
        "killing": float(clocks["killing"]),
        "normal": float(clocks["normal"]),
        "excision_flux": float(stocks["excision_flux"]),
        "outer_flux": float(stocks["outer_flux"]),
        "sat_debit": float(stocks["sat_debit"]),
        "remaining": float(probability["remaining"]),
        "excision_outflow": float(probability["excision_outflow"]),
        "outer_outflow": float(probability["outer_outflow"]),
        "closure": float(bundle.get(
            "d_closure",
            probability["remaining"] + probability["excision_outflow"]
            + probability["outer_outflow"] + 2.0 * stocks["sat_debit"],
        )),
    }


def _compare_step(plus, minus, analytic, step):
    found = {}
    worst_large = 0.0
    worst_large_name = None
    for name, exact in analytic.items():
        difference = (plus[name] - minus[name]) / (2.0 * step)
        error = _max_abs(difference - exact)
        signal = max(_max_abs(exact), _max_abs(difference))
        assert _within_signal(error, signal)
        if signal >= 1.0e-6:
            relative = _rel(error, signal)
            found[name] = relative
            if relative >= worst_large:
                worst_large = relative
                worst_large_name = name
        else:
            found[name] = error
    return found, worst_large, worst_large_name


def test_rk4_and_frozen_tangents_match_the_paired_step():
    report = []
    for layout in ("dense", "csr"):
        grid = _grid(points=17, layout=layout, r_out=16.0)
        state = _state(grid)
        direction = _sourced_direction(state)
        for dt in (0.001, 0.0005):
            analytic = _tangent_fields(response.rk4_tangent(state, grid, dt, direction))
            frozen = response.frozen_arm_tangent(state, grid, dt, direction)
            assert frozen["metric_variation_included"] is True
            assert frozen["mass_flux_advanced"] is False
            assert _max_abs(frozen["d_lapse"]) > 0.0
            assert frozen["dN_out"] == pytest.approx(0.0, abs=1e-12)
            assert frozen["d_stocks"]["excision_flux"] == 0.0
            assert frozen["d_inner_mass"] == 0.0

            def coupled(scale):
                return _output_fields(pg.rk4_step(
                    response.displace_state(state, direction, scale), grid, dt,
                ))

            def fixed(scale):
                moved = response.displace_state(state, direction, scale)
                rates = pg.stage_rates(moved, grid)
                updated = driver.fixed_rk4(moved, grid, driver.freeze_sourced_metric(rates), dt)
                return _output_fields(updated)

            coupled_errors, coupled_worst, coupled_name = _compare_step(
                coupled(1.0e-5), coupled(-1.0e-5), analytic, 1.0e-5,
            )
            coupled_half, coupled_half_worst, _coupled_half_name = _compare_step(
                coupled(5.0e-6), coupled(-5.0e-6), analytic, 5.0e-6,
            )
            frozen_errors, frozen_worst, frozen_name = _compare_step(
                fixed(1.0e-5), fixed(-1.0e-5), _tangent_fields(frozen), 1.0e-5,
            )
            frozen_half, frozen_half_worst, _frozen_half_name = _compare_step(
                fixed(5.0e-6), fixed(-5.0e-6), _tangent_fields(frozen), 5.0e-6,
            )
            assert _halving(coupled_errors["phi"], coupled_half["phi"], floor=1.0e-8)
            assert _halving(frozen_errors["phi"], frozen_half["phi"], floor=1.0e-8)
            report.append({
                "layout": layout,
                "dt": dt,
                "coupled": coupled_worst,
                "coupled_name": coupled_name,
                "coupled_half": coupled_half_worst,
                "frozen": frozen_worst,
                "frozen_name": frozen_name,
                "frozen_half": frozen_half_worst,
                "coupled_fields": coupled_errors,
                "frozen_fields": frozen_errors,
                "coupled_half_fields": coupled_half,
                "frozen_half_fields": frozen_half,
            })
    _record("rk4_fd", rows=[{
        key: row[key] for key in (
            "layout", "dt", "coupled", "coupled_name", "coupled_half", "frozen", "frozen_name", "frozen_half",
        )
    } for row in report])


def test_matched_proper_time_is_the_outer_static_protocol():
    protocol = response.OUTER_STATIC_PROTOCOL
    assert response.matched_proper_derivative(
        0.2, 0.5, 0.1, 2.0, protocol=protocol,
    ) == pytest.approx(0.175)
    assert response.matched_proper_derivative(
        0.3, 4.0, 0.3, 4.0, protocol=protocol,
    ) == pytest.approx(0.0)
    grid = _grid()
    state = _state(grid)
    direction = _sourced_direction(state)
    bundle = response.stage_tangent(state, grid, direction)
    proper = response.declared_proper_clock(bundle)
    assert proper["proper_clock"] is True
    assert proper["protocol"] == protocol
    assert proper["node"] == "outer"
    assert proper["legacy_normal_is_proper"] is False
    assert bundle["clock_semantics"]["normal"]["proper_clock"] is False
    assert bundle["clock_semantics"]["normal"]["protocol"] == response.LEGACY_EXCISION_LAPSE
    assert bundle["clock_semantics"]["normal"]["formula"] == "N[0]"
    assert bundle["clock_semantics"]["killing"]["formula"] == "sqrt(N[-1]**2 - shift[-1]**2)"
    assert bundle["clocks_move_the_metric"] is False
    assert bundle["d_normal_rate"] == pytest.approx(float(bundle["d_lapse"][0]))
    assert bundle["rates"]["normal_rate"] == pytest.approx(float(bundle["lapse"][0]))
    killing = math.sqrt(bundle["lapse"][-1] ** 2 - bundle["shift"][-1] ** 2)
    assert bundle["rates"]["killing_rate"] == pytest.approx(killing)
    matched = response.matched_proper_derivative(
        bundle["d_coordinate_energy"],
        0.0,
        proper["d_rate"],
        proper["rate"],
        protocol=protocol,
    )
    assert matched == pytest.approx(bundle["d_coordinate_energy"])
    delayed = response.matched_proper_derivative(
        proper["d_rate"], proper["rate"], proper["d_rate"], proper["rate"],
        protocol=protocol,
    )
    assert delayed == pytest.approx(0.0)
    for refused, message in (
        ("normal", "not a proper clock"),
        (response.LEGACY_EXCISION_LAPSE, "not a proper clock"),
        (response.TRANSPORTED_NORMAL_PROTOCOL, "not supplied"),
    ):
        with pytest.raises(ValueError, match=message):
            response.declared_proper_clock(bundle, protocol=refused)
        with pytest.raises(ValueError, match=message):
            response.matched_proper_derivative(
                bundle["d_coordinate_energy"], 0.0,
                bundle["d_normal_rate"], bundle["rates"]["normal_rate"],
                protocol=refused,
            )
    stage = response.observer_metric_stage(bundle, grid)
    assert stage["transported_clock_computed"] is False
    assert stage["legacy_normal_not_a_sample"] is True
    assert stage["field_rk4_stage_fractions"] == (0.0, 0.5, 0.5, 1.0)
    assert np.shares_memory(stage["lapse"], bundle["lapse"]) or _max_abs(
        stage["lapse"] - bundle["lapse"]
    ) == 0.0
    assert stage["radius"][0] == pytest.approx(float(grid["radius"][0]))
    assert stage["excision_index"] == 0 and stage["outer_index"] == -1
    _record(
        "proper_clock_contract",
        killing_self_match=delayed,
        legacy_normal_proper=False,
        transported_normal_available=bundle["clock_semantics"]["transported_normal_available"],
    )


def test_saved_t0_stage_and_one_rk4_are_cheap():
    if not ARCHIVE_T0.is_file():
        pytest.skip("archived pilot t0 is not on this host")
    payload = np.load(ARCHIVE_T0)
    grid = pg.build_grid(1.04, r_m=1.0, points=int(payload["radius"].shape[0]), r_out=32.0, layout="csr")
    assert _max_abs(grid["radius"] - payload["radius"]) == 0.0
    state = pg.make_state(
        grid, payload["coupled_phi"], occupations=payload["occupations"],
        inner_mass=float(payload["coupled_inner_mass"][0]), killing_frequency=1.0,
    )
    assert float(state["clocks"]["t"]) == 0.0
    direction = _sourced_direction(state)
    direction["d_occupations"] = np.array([0.2 * float(state["occupations"][0])])
    rates = pg.stage_rates(state, grid)
    choice = response.production_step(rates["cfl_dt"], grid["r_m"])
    assert choice["dt"] == pytest.approx(0.001)
    assert choice["limited_by"] == "step_cap"
    analytic = response.rk4_tangent(state, grid, choice["dt"], direction)
    frozen = response.frozen_arm_tangent(state, grid, choice["dt"], direction)
    assert _max_abs(frozen["d_lapse"]) > 0.0

    def coupled(scale):
        return _output_fields(pg.rk4_step(response.displace_state(state, direction, scale), grid, choice["dt"]))

    coarse, worst, worst_name = _compare_step(
        coupled(1.0e-6), coupled(-1.0e-6), _tangent_fields(analytic), 1.0e-6,
    )
    _half, half_worst, _half_name = _compare_step(
        coupled(5.0e-7), coupled(-5.0e-7), _tangent_fields(analytic), 5.0e-7,
    )
    assert _halving(coarse["phi"], _half["phi"], floor=1.0e-8)
    _record(
        "saved_t0",
        points=int(grid["points"]),
        dt=choice["dt"],
        cfl_half=choice["cfl_half_step"],
        coupled_rel=worst,
        coupled_name=worst_name,
        coupled_rel_half=half_worst,
        frozen_metric_variation=_max_abs(frozen["d_lapse"]),
        fields=coarse,
    )


def test_frequency_preparation_uses_two_phase_aligned_steps():
    grid = _grid(points=17, r_out=20.0)
    request = source.packet_request(grid["r_m"])
    started = time.process_time()
    prepared = response.frequency_preparation_derivative(
        1.0e-3,
        mass_over_rm=1.04,
        energy_target=1.0e-3,
        rank=9,
        model=grid["model"],
        center=request["center"],
        width=request["width"],
        omega=request["omega"],
        radii=grid["radius"],
        weights=grid["weights"],
        derivative=grid["derivative"],
    )
    elapsed = time.process_time() - started
    assert prepared["phase_aligned"] is True
    assert prepared["modes_move"] is True
    assert prepared["variational_ode"] is False
    assert prepared["scattering_campaign"] is False
    assert abs(prepared["relative_filling_plus_relative_energy"]) < 1e-6
    ell = np.linspace(-1.0, 1.0, prepared["d_sample_energies"].shape[0])
    assert _max_abs(prepared["d_sample_energies"] - (1.0 + 0.9 * ell)) < 1.0e-6
    assert prepared["phi_step_discrepancy"] / prepared["phi_signal"] < 5.0e-2
    assert elapsed < 30.0
    _record(
        "frequency_fd",
        cpu_seconds=elapsed,
        phi_step_gap=prepared["phi_step_discrepancy"] / prepared["phi_signal"],
        filling_energy_identity=prepared["relative_filling_plus_relative_energy"],
        phase_aligned=prepared["phase_aligned"],
    )
