"""Manufactured leading parent-observer checks. No production trajectory is loaded."""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_extent as extent
from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_parent_observer as parent
from recursive_horizons import nsc_discovery_regions as regions
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def _block_chi(monkeypatch):
    for module, name in (
        (leading.galerkin, "compose_fine_hamiltonian"),
        (leading.galerkin, "rates"),
        (leading.coupling, "geometric_rates"),
        (leading.tidal, "analytic_accelerations"),
        (regional, "geometric_rates"),
        (regional, "compose_fine_hamiltonian"),
    ):
        monkeypatch.setattr(module, name, lambda *_args, **_kwargs: pytest.fail("chi carrier called"))


def _envelope(grid, focus):
    distance = parent.signed_distance(grid.xi_q, grid.length)
    bump = np.exp(-0.5 * ((np.abs(distance) - focus) / 0.22) ** 2)
    phase = np.exp(1j * np.pi * np.asarray(grid.xi_q, dtype=float) / float(grid.length))
    fine = bump * phase
    coefficient = grid.U_f.conj().T @ fine
    return coefficient, grid.U_f @ coefficient


def _holder(focus):
    grid = galerkin.build_grid(32, gauge="conformal")
    rank = 4
    generator = np.random.default_rng(19)
    raw = generator.normal(size=(grid.nf, rank)) + 1j * generator.normal(size=(grid.nf, rank))
    orthogonal, _ = np.linalg.qr(raw)
    bump, reconstructed = _envelope(grid, focus)
    phi0 = orthogonal / np.sqrt(2.0)
    phi0[:, 0] = bump / max(np.linalg.norm(bump), 1.0e-30) / np.sqrt(2.0)
    phi1 = 1j * phi0
    weights = np.array([0.82, 0.11, 0.07, 0.04])
    grid.fine = replace(grid.fine, occupations=weights.copy())
    length = float(grid.length)
    centre = 0.5 * length
    coordinate = np.asarray(grid.xi_g, dtype=float)
    target = centre + focus
    radius = 1.35 + 0.08 * np.cos(2.0 * np.pi * coordinate / length)
    lapse = 0.9 + 0.04 * np.sin(2.0 * np.pi * coordinate / length)
    geometry = np.eye(grid.ng)
    geometry.setflags(write=False)
    pair = SimpleNamespace(
        grid=grid,
        geometry_map=geometry,
        weights=weights.copy(),
        reference_columns=np.vstack((phi0, phi1)),
        source_columns=np.vstack((phi0, phi1)),
        source_metadata={"layout": "variable_rank", "rank": rank, "focus": focus},
        geometry_metadata={"W": "identity", "all_geometry_degrees_retained": True},
        child_interval=(centre - 0.5, centre + 0.5),
        parent_interval=(centre - 3.0, centre + 3.0),
        clock_locations=(centre - 1.0, centre, centre + 1.0),
    )
    nodal = leading.State(
        lapse, radius,
        0.02 * np.cos(2.0 * np.pi * coordinate / length),
        0.015 * np.sin(2.0 * np.pi * coordinate / length),
        phi0, phi1,
    )
    return pair, leading.encode(pair, nodal), reconstructed


def test_signed_distance_regions_are_disjoint_and_cover_the_carrier():
    grid = galerkin.build_grid(32, gauge="conformal")
    partition = parent.region_partition(grid.length)
    masks = parent.region_masks(grid)
    assert partition["centre"] == pytest.approx(4.0)
    assert partition["protected_collar_accounted_separately"] is False
    assert partition["ambient"]["wrapped"] is True
    ones = np.ones(grid.nq)
    total = sum(parent._piece_integral(grid, ones, partition[name]) for name in parent.REGION_ORDER)
    assert total == pytest.approx(float(grid.length), abs=1.0e-12)
    covered = np.zeros(grid.nq, dtype=int)
    for name in parent.REGION_ORDER:
        assert np.count_nonzero(masks[name]) > 0
        covered += masks[name].astype(int)
    assert np.all(covered == 1)
    assert abs(float(masks["signed_distance"][int(np.argmin(np.abs(grid.xi_q - 4.0)))])) <= 0.5


def test_observe_uses_leading_rates_actual_weights_and_wrapped_ambient(monkeypatch):
    _block_chi(monkeypatch)
    with threadpool_limits(limits=1):
        pair, state, _reconstructed = _holder(0.0)
        row = parent.observe(pair, state, 0.2)
    assert row["chi_rate_used"] is False
    assert row["passed_bundle"] is False
    assert row["rank"] == 4
    assert row["weights"] == pytest.approx([0.82, 0.11, 0.07, 0.04])
    assert row["positive_proper_probability"] is True
    assert row["per_proper_length_minimum"] >= 0.0
    assert row["piston_work_included"] is False
    assert row["independent_energy_parcel"] is False
    assert row["column_ancestry_used_as_energy"] is False
    assert row["held_population_measured"] is False
    assert row["flags"]["dominance_forced"] is False
    assert row["holder_intervals"]["matches_signed_distance"] is True
    assert row["holder_intervals"]["declared_intervals_used_as_measurement_cuts"] is False
    assert row["partition"]["accounted_parent_is_annulus"] is True
    assert row["partition"]["overlapping_parent_sum_used"] is False
    probability = row["channels"]["probability"]
    assert probability["total"] == pytest.approx(row["probability_total_nodal"], abs=1.0e-8)
    assert sum(probability[name] for name in parent.REGION_ORDER) == pytest.approx(probability["total"])
    for name in parent.REGION_ORDER:
        ledger = row["ledgers"][name]
        assert ledger["piston_work_included"] is False
        assert ledger["surface_work"] == 0.0
        assert ledger["physical_wall"] is False
        assert ledger["boundary_velocity"] == [0.0, 0.0]
        assert abs(ledger["closure_residual"]) < 1.0e-8
        assert row["channels"]["normal_energy_rate"][name] == pytest.approx(ledger["rate"])
    assert row["ledgers"]["ambient"]["wrapped"] is True
    assert row["null_crossings"]["coordinate_speed"] == [-1.0, 1.0]
    assert row["null_crossings"]["proper_clocks_integrated"] is False
    assert row["proper_clocks"]["definition"] == "dτ/dT = rQ"
    assert row["proper_clocks"]["rates"] == pytest.approx(leading.clock_rates(pair, state).tolist())
    assert row["metric"]["chi_substituted_for_curvature"] is False
    assert np.isfinite(row["metric"]["R"]["child"]["max"])
    assert np.isfinite(row["metric"]["W"]["parent"]["max"])
    assert np.isfinite(row["metric"]["tides"]["R_0101"]["inner"]["max_abs"])
    assert np.isfinite(row["metric"]["base"]["normal_tidal_radial"]["ambient"]["max_abs"])
    assert row["contours"]["continuum_optimized_peak"] is False
    assert row["contours"]["child_peak_replaced_global"] is False
    assert row["state_token_unchanged"] is True
    regions.jsonable(row)


def test_passed_leading_jets_are_used_and_chi_bundles_are_refused(monkeypatch):
    _block_chi(monkeypatch)
    with threadpool_limits(limits=1):
        pair, state, _reconstructed = _holder(0.0)
        rate, bundle = leading.rates(pair, state, return_bundle=True)
        curvature, velocity, acceleration = leading.metric_jets(pair, state, rate, bundle)
    tides = dict(curvature["tides"])
    tides["owned_W"] = np.ones_like(curvature["tides"]["owned_W"])
    sent = dict(curvature)
    sent["tides"] = tides
    passed = dict(bundle)
    passed["leading_rate"] = rate
    passed["metric_jets"] = (sent, velocity, acceleration)
    monkeypatch.setattr(leading, "rates", lambda *_args, **_kwargs: pytest.fail("rates recomputed"))
    monkeypatch.setattr(leading, "metric_jets", lambda *_args, **_kwargs: pytest.fail("jets recomputed"))
    row = parent.observe(pair, state, 0.2, bundle=passed)
    assert row["passed_bundle"] is True
    assert row["chi_rate_used"] is False
    for name in parent.REGION_ORDER:
        assert row["metric"]["W"][name]["min"] == pytest.approx(1.0)
        assert row["metric"]["W"][name]["max"] == pytest.approx(1.0)
    missing = dict(bundle)
    with pytest.raises(ValueError, match="leading_rate"):
        parent.observe(pair, state, 0.2, bundle=missing)
    chi_bundle = dict(passed)
    chi_bundle["fine_state"] = SimpleNamespace(chi=np.zeros(1))
    with pytest.raises(ValueError, match="chi"):
        parent.observe(pair, state, 0.2, bundle=chi_bundle)
    chi_state = state.copy()
    chi_state.chi = np.zeros(1)
    with pytest.raises(ValueError, match="chi"):
        parent.observe(pair, chi_state, 0.2)


def test_parent_peak_does_not_replace_the_child_peak_and_frozen_pressure_vanishes(monkeypatch):
    _block_chi(monkeypatch)
    with threadpool_limits(limits=1):
        parent_pair, parent_state, parent_bump = _holder(2.0)
        child_pair, child_state, child_bump = _holder(0.0)
        dominated = parent.observe(parent_pair, parent_state, 0.4)
        centred = parent.observe(child_pair, child_state, 0.4)
        frozen = parent.observe(child_pair, child_state, 0.4, control_mode="frozen_geometry")
    assert float(np.max(np.abs(parent_bump))) > 0.0
    assert float(np.max(np.abs(child_bump))) > 0.0
    assert dominated["flags"]["dominance_forced"] is False
    assert dominated["flags"]["parent_peak_dominates"] is True
    assert dominated["flags"]["global_peak_region"] == "parent"
    assert dominated["child_peak"]["region"] == "child"
    assert dominated["child_peak"]["replaces_global_contours"] is False
    assert dominated["child_peak"]["is_global_peak"] is False
    assert dominated["contours"]["global_peak_region"] == "parent"
    assert centred["flags"]["parent_peak_dominates"] is False
    assert centred["flags"]["global_peak_region"] == "child"
    assert centred["child_peak"]["is_global_peak"] is True
    assert centred["flags"]["dominance_forced"] is False
    assert frozen["rates"]["geometry_jets_zero"] is True
    assert frozen["rates"]["Q_dot_max"] == 0.0
    assert frozen["rates"]["r_dot_max"] == 0.0
    assert frozen["channels"]["pressure"]["child"] == pytest.approx(0.0, abs=1.0e-12)
    assert frozen["ledgers"]["child"]["boundary_velocity"] == [0.0, 0.0]
    assert frozen["piston_work_included"] is False


def _reduced_geometry_map(count, seed=14):
    generator = np.random.default_rng(seed)
    coordinate = np.arange(count, dtype=float)
    columns = [np.ones(count)]
    for mode in range(1, 4):
        angle = 2.0 * np.pi * mode * coordinate / count
        columns.extend((np.cos(angle), np.sin(angle)))
    while len(columns) < count:
        columns.append(generator.normal(size=count))
    matrix, _ = np.linalg.qr(np.column_stack(columns))
    matrix.setflags(write=False)
    return matrix


def _high_mode_holder():
    """Generic AP phases, a mode-14 lapse, and a non-identity geometry frame."""
    grid = galerkin.build_grid(32, gauge="conformal")
    rank = 3
    generator = np.random.default_rng(14)
    phase = generator.normal(size=(grid.nf, rank)) + 1j * generator.normal(size=(grid.nf, rank))
    other = generator.normal(size=(grid.nf, rank)) + 1j * generator.normal(size=(grid.nf, rank))
    phi0, _ = np.linalg.qr(phase)
    phi1, _ = np.linalg.qr(other)
    phi0 = phi0 / np.sqrt(2.0)
    phi1 = phi1 / np.sqrt(2.0)
    weights = np.array([0.23, 0.51, 0.74])
    grid.fine = replace(grid.fine, occupations=weights.copy())
    length = float(grid.length)
    centre = 0.5 * length
    coordinate = np.asarray(grid.xi_g, dtype=float)
    wave = 2.0 * np.pi / length
    lapse = 0.9 + 0.12 * np.cos(14.0 * wave * coordinate)
    radius = 1.28 + 0.04 * np.sin(wave * coordinate)
    pair = SimpleNamespace(
        grid=grid,
        geometry_map=_reduced_geometry_map(grid.ng),
        weights=weights.copy(),
        reference_columns=np.vstack((phi0, phi1)),
        source_columns=np.vstack((phi0, phi1)),
        source_metadata={"layout": "generic_phase", "rank": rank},
        geometry_metadata={"W": "reduced_low_mode_completion", "all_geometry_degrees_retained": True},
        child_interval=(centre - 0.5, centre + 0.5),
        parent_interval=(centre - 3.0, centre + 3.0),
        clock_locations=(centre - 1.0, centre, centre + 1.0),
    )
    nodal = leading.State(
        lapse, radius,
        0.03 * np.cos(wave * coordinate),
        0.02 * np.sin(2.0 * wave * coordinate),
        phi0, phi1,
    )
    return pair, leading.encode(pair, nodal)


def _per_proper(pair, fine):
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * pair.weights, axis=1)
    return mass / (float(pair.grid.dx_q) * fine.r * fine.Q)


def _child_shell(pair, fine):
    system = leading.active_system(pair.grid, fine)
    nodal = regional.matter_ledger(system, fine)["normal_energy_nodal"]
    partition = parent.region_partition(pair.grid.length)
    return parent._piece_integral(pair.grid, nodal / float(pair.grid.dx_q), partition["child"])


def test_projected_flow_matches_the_decoded_state_and_not_the_raw_rate(monkeypatch):
    _block_chi(monkeypatch)
    with threadpool_limits(limits=1):
        pair, state = _high_mode_holder()
        rate, bundle = leading.rates(pair, state, return_bundle=True)
        projected, diagnostic = parent.projected_state_rate(pair, rate, bundle)
        raw = leading.State(*bundle["unprojected_rates"])
        fine = bundle["fine_state"]
        system = bundle["fine_system"]
        row = parent.observe(pair, state, 0.15, bundle={**bundle, "leading_rate": rate})
    assert diagnostic["unprojected_rate_is_state_flow"] is False
    assert diagnostic["field_projection_gap_max"] > 1.0e-3
    assert row["state_flow"] == parent.STATE_FLOW
    assert row["rate_projection"]["unprojected_rate_is_state_flow"] is False
    assert row["rate_projection"]["field_projection_gap_max"] == pytest.approx(diagnostic["field_projection_gap_max"])
    step = 1.0e-6
    ahead = leading.fine_state(pair, leading.combine(state, rate, step))
    behind = leading.fine_state(pair, leading.combine(state, rate, -step))
    probability_change = (_per_proper(pair, ahead) - _per_proper(pair, behind)) / (2.0 * step)
    projected_probability = regions._probability_rate(fine, system, projected)
    raw_probability = regions._probability_rate(fine, system, raw)
    scale = max(float(np.max(np.abs(probability_change))), 1.0e-12)
    assert float(np.max(np.abs(probability_change - projected_probability))) / scale < 2.0e-4
    assert float(np.max(np.abs(raw_probability - projected_probability))) > 1.0e-3
    peak = int(row["contours"]["peak_index"])
    assert row["contours"]["peak_rate"] == pytest.approx(float(probability_change[peak]), abs=1.0e-5)
    assert abs(float(raw_probability[peak]) - float(projected_probability[peak])) > 1.0e-4
    shell_change = (_child_shell(pair, ahead) - _child_shell(pair, behind)) / (2.0 * step)
    projected_slope = regions.analytic_normal_energy_slope(system, fine, projected, bundle["source"])
    raw_slope = regions.analytic_normal_energy_slope(system, fine, raw, bundle["source"])
    window = parent._piece_integral(pair.grid, projected_slope / float(pair.grid.dx_q), parent.region_partition(pair.grid.length)["child"])
    raw_window = parent._piece_integral(pair.grid, raw_slope / float(pair.grid.dx_q), parent.region_partition(pair.grid.length)["child"])
    assert row["channels"]["normal_energy_rate"]["child"] == pytest.approx(shell_change, rel=2.0e-4, abs=1.0e-7)
    assert row["channels"]["normal_energy_rate"]["child"] == pytest.approx(window, abs=1.0e-8)
    assert abs(raw_window - window) > 1.0e-4
    assert row["channels"]["normal_energy_rate"]["child"] * shell_change >= -1.0e-10
    raw_pressure = regional.proper_balance_terms(system, fine, raw, regional.matter_ledger(system, fine))
    assert row["channels"]["pressure"]["child"] == pytest.approx(
        parent._fixed_ledger(pair.grid, {
            "energy": np.zeros(pair.grid.nq),
            "flux": np.zeros(pair.grid.nq),
            "pressure": regional.proper_balance_terms(system, fine, projected)["proper_pressure_work"],
            "lapse": np.zeros(pair.grid.nq),
            "slope": np.zeros(pair.grid.nq),
        }, parent.region_partition(pair.grid.length)["child"])["pressure"]
    )
    assert row["ledgers"]["child"]["piston_work_included"] is False
    del raw_pressure


def test_source_free_zero_weights_keep_the_label_and_report_no_matter(monkeypatch):
    _block_chi(monkeypatch)
    with pytest.raises(ValueError, match="control_mode"):
        episode.normalize_control_mode("source_free")
    with threadpool_limits(limits=1):
        grid = galerkin.build_grid(32, gauge="conformal")
        rank = 3
        weights = np.zeros(rank)
        grid.fine = replace(grid.fine, occupations=weights.copy())
        phi0 = np.zeros((grid.nf, rank), dtype=complex)
        phi1 = np.zeros((grid.nf, rank), dtype=complex)
        length = float(grid.length)
        centre = 0.5 * length
        coordinate = np.asarray(grid.xi_g, dtype=float)
        wave = 2.0 * np.pi / length
        geometry = np.eye(grid.ng)
        geometry.setflags(write=False)
        pair = SimpleNamespace(
            grid=grid,
            geometry_map=geometry,
            weights=weights.copy(),
            reference_columns=np.vstack((phi0, phi1)),
            source_columns=np.vstack((phi0, phi1)),
            source_metadata={"control_mode": "source_free", "source_reprepared": True},
            geometry_metadata={"W": "identity"},
            child_interval=(centre - 0.5, centre + 0.5),
            parent_interval=(centre - 3.0, centre + 3.0),
            clock_locations=(centre,),
        )
        nodal = leading.State(
            0.9 + 0.05 * np.cos(wave * coordinate),
            1.2 + 0.03 * np.sin(wave * coordinate),
            0.04 * np.cos(wave * coordinate),
            0.05 * np.ones(grid.ng),
            phi0, phi1,
        )
        state = leading.encode(pair, nodal)
        row = parent.observe(pair, state, 0.0, control_mode="source_free")
        invalid = state.copy()
        invalid.phi0 = np.array(invalid.phi0, copy=True)
        invalid.phi0[0, 0] = np.nan
    assert row["control_mode"] == "source_free"
    assert row["source_free"] is True
    assert row["source_free_label_preserved"] is True
    assert row["geometry_rate_mode"] == "coupled"
    assert row["rates"]["geometry_jets_zero"] is False
    assert row["rates"]["Q_dot_max"] > 0.0 or row["rates"]["r_dot_max"] > 0.0
    assert row["matter_present"] is False
    assert row["matter_curves"] == "flat"
    assert row["positive_proper_probability"] is False
    assert row["proper_probability_nonnegative"] is True
    assert row["contours"]["flags"]["flat"] is True
    assert row["contours"]["levels"] == []
    assert row["channels"]["probability"]["total"] == pytest.approx(0.0, abs=1.0e-12)
    assert row["weights"] == pytest.approx([0.0, 0.0, 0.0])
    assert row["piston_work_included"] is False
    with pytest.raises(ValueError, match="nonfinite"):
        parent.observe(pair, invalid, 0.0, control_mode="source_free")
    pair.weights = np.array([-0.2, 0.2, 0.2])
    with pytest.raises(ValueError, match="\\[0, 1\\]"):
        parent.observe(pair, state, 0.0, control_mode="source_free")
