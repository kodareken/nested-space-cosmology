"""Stage-5 first-wave directional primitive. Centred differences are the oracle only."""
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_discovery_response as discovery
from recursive_horizons import nsc_influence as influence
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


_LAB = Path(__file__).resolve().parents[1]
_V5 = _LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
_RATE_FIELDS = (
    "Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1",
    "force_L", "force_Q", "force_beta", "fieldwork_power",
)


def _assert_match(analytic, numeric, name):
    gap = float(np.max(np.abs(np.asarray(analytic) - np.asarray(numeric))))
    scale = max(
        1e-8,
        float(np.max(np.abs(np.asarray(analytic)))),
        float(np.max(np.abs(np.asarray(numeric)))),
    )
    assert gap <= 2e-5 * scale + 2e-8, (name, gap, scale, gap / scale)


def _manufactured_state(nf=16, quadrature=32, seed=4):
    grid = galerkin.build_grid(nf, quadrature=quadrature, gauge="conformal")
    phi0, phi1 = galerkin.manufactured_columns(nf, seed=seed)
    state = galerkin.blank_state(grid, phi0, phi1)
    angle = 2 * np.pi * grid.xi_g / grid.length
    state.Q = state.Q * (1 + 0.03 * np.sin(angle))
    state.r = 1 + 0.04 * np.cos(angle)
    state.chi = 0.02 * np.sin(2 * angle)
    state.p_Q = 0.01 * np.cos(angle)
    state.p_r = 0.015 * np.sin(angle)
    state.p_chi = 0.012 * np.cos(2 * angle)
    # A pure σ2 column has vanishing S1, so the conformal sum F_L+F_Q
    # does not see an occupation tangent. A relative phase makes S1 nonzero.
    state.phi1 = np.exp(-0.5j) * state.phi1
    return grid, state


def _direction(grid, name, seed=7):
    tangent = discovery.zero_tangent(grid)
    generator = np.random.default_rng(seed)
    angle = 2 * np.pi * grid.xi_g / grid.length
    if name == "Q":
        tangent.Q = np.sin(angle)
    elif name == "r":
        tangent.r = np.cos(angle)
    elif name == "chi":
        tangent.chi = np.sin(2 * angle)
    elif name == "p_Q":
        tangent.p_Q = np.cos(angle)
    elif name == "p_r":
        tangent.p_r = np.sin(angle)
    elif name == "p_chi":
        tangent.p_chi = np.cos(2 * angle)
    elif name == "phi":
        raw0 = generator.normal(size=tangent.phi0.shape) + 1j * generator.normal(size=tangent.phi0.shape)
        raw1 = generator.normal(size=tangent.phi1.shape) + 1j * generator.normal(size=tangent.phi1.shape)
        tangent.phi0 = raw0 / np.linalg.norm(raw0)
        tangent.phi1 = raw1 / np.linalg.norm(raw1)
    elif name == "occupations":
        tangent.occupations[1] = 1.0
    else:
        raise AssertionError(name)
    return tangent


def _displace(state, tangent, step):
    return coupling.CauchyState(
        *(getattr(state, name) + step * getattr(tangent, name) for name in nested.STATE_NAMES)
    )


def _centred_galerkin(grid, state, tangent, h=1e-6):
    base = np.array(grid.fine.occupations, copy=True)

    def evaluate(sign):
        grid.fine.occupations = base + sign * h * tangent.occupations
        try:
            return galerkin.rates(grid, _displace(state, tangent, sign * h))
        finally:
            grid.fine.occupations = base

    plus, minus = evaluate(1), evaluate(-1)
    return {
        name: (getattr(plus, name) - getattr(minus, name)) / (2 * h)
        for name in _RATE_FIELDS
    }


def test_contract_names_the_missing_stress_and_does_not_call_newton():
    contract = discovery.specification()
    assert contract["schema"] == "NSC-DISCOVERY-RESPONSE-v1"
    assert contract["method"] == "analytic_jacobian_vector"
    assert contract["finite_difference_used_as_analytic"] is False
    assert contract["newton_production"] is False
    assert contract["primary_readout"] == "child_regional_content"
    assert contract["secondary_readout"] == "child_proper_mean_r"
    assert contract["kernel_2x2_is_spatial_outside"] is False
    assert "ctp_future_effective_stress" in contract["missing_primitives"]
    assert "retarded_delta_state" in contract["missing_primitives"]
    assert "SpectralInducedSource" in contract["induced_owner"]
    source = Path(discovery.__file__).read_text()
    assert "solve_initial_radius" not in source
    assert "(plus - minus) / (2" not in source


def test_manufactured_galerkin_directions_match_centred_differences():
    grid, state = _manufactured_state()
    assert grid.fine.multiplicity == 4 * grid.fine.kappa
    zero = discovery.rate_jacobian_vector(grid, state, discovery.zero_tangent(grid))
    for name in _RATE_FIELDS:
        assert np.max(np.abs(getattr(zero, name))) < 1e-10
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi", "occupations"):
        tangent = _direction(grid, name)
        analytic = discovery.rate_jacobian_vector(grid, state, tangent)
        numeric = _centred_galerkin(grid, state, tangent)
        for field in _RATE_FIELDS:
            _assert_match(getattr(analytic, field), numeric[field], f"{name}.{field}")
    geometry = discovery.rate_jacobian_vector(grid, state, _direction(grid, "Q"))
    assert np.max(np.abs(geometry.phi0)) > 1e-8
    occupation = discovery.rate_jacobian_vector(grid, state, _direction(grid, "occupations"))
    assert np.max(np.abs(occupation.phi0)) < 1e-10
    assert np.max(np.abs(occupation.r)) < 1e-10
    assert np.max(np.abs(occupation.p_Q)) > 1e-8
    assert np.max(np.abs(occupation.force_L)) > 0
    assert np.max(np.abs(occupation.force_Q)) > 0
    bare = discovery.rate_jacobian_vector(
        grid, state, _direction(grid, "phi"), include_matter_force=False,
    )
    full = discovery.rate_jacobian_vector(grid, state, _direction(grid, "phi"))
    pulled = galerkin.pull_geometry(grid, (full.force_Q + full.force_L) / grid.dx_q)
    _assert_match(full.p_Q, bare.p_Q - pulled, "conformal force_L+force_Q")


def test_nested_canonical_map_matches_centred_differences():
    pair = nested.build_pair(64, quadrature=64)
    phi0, phi1 = galerkin.manufactured_columns(pair.grid.nf, seed=5)
    nodal = galerkin.blank_state(pair.grid, phi0, phi1)
    nodal.phi1 = np.exp(-0.5j) * nodal.phi1
    angle = 2 * np.pi * pair.grid.xi_g / pair.grid.length
    nodal.Q = nodal.Q * (1 + 0.02 * np.sin(angle))
    nodal.r = 1.1 + 0.03 * np.cos(angle)
    nodal.chi = 0.01 * np.sin(2 * angle)
    nodal.p_r = 0.02 * np.sin(angle)
    nodal.p_chi = 0.01 * np.cos(angle)
    state = nested.encode_state(pair, nodal)
    for name in ("Q", "r", "p_r", "phi", "occupations"):
        nodal_tangent = _direction(pair.grid, name, seed=9)
        increment = coupling.CauchyState(
            nodal_tangent.Q, nodal_tangent.r, nodal_tangent.chi,
            nodal_tangent.p_Q, nodal_tangent.p_r, nodal_tangent.p_chi,
            nodal_tangent.phi0, nodal_tangent.phi1,
        )
        canonical = nested.encode_state(pair, increment)
        tangent = discovery.StateTangent(
            canonical.Q, canonical.r, canonical.chi,
            canonical.p_Q, canonical.p_r, canonical.p_chi,
            canonical.phi0, canonical.phi1, nodal_tangent.occupations.copy(),
        )
        analytic = discovery.nested_rate_jacobian_vector(pair, state, tangent)
        nodal_analytic = discovery.rate_jacobian_vector(pair.grid, nodal, nodal_tangent)
        _assert_match(analytic.r, pair.geometry_map.T @ nodal_analytic.r, f"{name}.canonical r")
        _assert_match(
            analytic.p_r,
            pair.grid.dx_g * pair.geometry_map.T @ nodal_analytic.p_r,
            f"{name}.canonical pi",
        )
        base = np.array(pair.grid.fine.occupations, copy=True)

        def evaluate(sign, step_tangent=tangent, stored=base):
            pair.grid.fine.occupations = stored + sign * 1e-6 * step_tangent.occupations
            trial = nested.NestedState(*(
                getattr(state, field) + sign * 1e-6 * getattr(step_tangent, field)
                for field in nested.STATE_NAMES
            ))
            try:
                return nested.rates(pair, trial)
            finally:
                pair.grid.fine.occupations = stored

        plus, minus = evaluate(1), evaluate(-1)
        for field in _RATE_FIELDS:
            numeric = (getattr(plus, field) - getattr(minus, field)) / 2e-6
            _assert_match(getattr(analytic, field), numeric, f"nested {name}.{field}")


def test_linearized_radius_tangent_cancels_the_source_residual():
    grid = galerkin.build_grid(16, quadrature=32, gauge="conformal")
    phi0, phi1 = galerkin.manufactured_columns(16, seed=6)
    state = galerkin.blank_state(grid, phi0, phi1)
    tangent = _direction(grid, "phi", seed=8)
    prepared = discovery.prepared_radius_tangent(grid, state, tangent)
    assert prepared.available is True
    assert prepared.newton_used is False
    assert prepared.missing_primitive is None
    assert prepared.linear_residual_max < 1e-9
    assert prepared.jacobian_owner.endswith("_projected_radius_jacobian")

    def residual(radius, column0, column1, weights):
        grid.fine.occupations = weights
        trial = state.copy()
        trial.r = radius
        trial.phi0 = column0
        trial.phi1 = column1
        fine = galerkin.prolong_state(grid, trial)
        source = coupling.source_from_columns(galerkin.active_fine_system(grid, fine), fine)
        projected, _full = galerkin.projected_radius_operator(
            grid, radius, source["force_L"] / grid.dx_q,
        )
        return projected

    base = np.array(grid.fine.occupations, copy=True)
    step = 1e-6
    source_only = (
        residual(state.r, phi0 + step * tangent.phi0, phi1 + step * tangent.phi1, base)
        - residual(state.r, phi0 - step * tangent.phi0, phi1 - step * tangent.phi1, base)
    ) / (2 * step)
    joint = (
        residual(state.r + step * prepared.delta_r, phi0 + step * tangent.phi0, phi1 + step * tangent.phi1, base)
        - residual(state.r - step * prepared.delta_r, phi0 - step * tangent.phi0, phi1 - step * tangent.phi1, base)
    ) / (2 * step)
    grid.fine.occupations = base
    source_scale = max(1e-8, float(np.max(np.abs(source_only))))
    assert float(np.max(np.abs(joint))) <= 2e-5 * source_scale + 2e-8
    occupation = _direction(grid, "occupations")
    occupation_radius = discovery.prepared_radius_tangent(grid, state, occupation)
    assert occupation_radius.available is True
    assert np.max(np.abs(occupation_radius.delta_r)) > 0
    with pytest.raises(ValueError, match="zero chi"):
        moved = state.copy()
        moved.chi = np.full(grid.ng, 1e-3)
        discovery.prepared_radius_tangent(grid, moved, tangent)


def test_readouts_clock_and_matched_tau_correction():
    grid, state = _manufactured_state()
    phi = _direction(grid, "phi", seed=2)
    radius = _direction(grid, "r")
    content, delta_content = discovery.child_regional_content(grid, state, phi)
    _frozen, frozen_content = discovery.child_regional_content(grid, state, radius)
    assert frozen_content == pytest.approx(0.0, abs=1e-12)
    fine = galerkin.prolong_state(grid, state)
    weights = grid.fine.occupations
    density = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights, axis=1) / grid.dx_q
    assert content == pytest.approx(nested.interval_integral(grid, density, nested.CHILD_INTERVAL))
    step = 1e-6
    plus = discovery.child_regional_content(grid, _displace(state, phi, step))
    minus = discovery.child_regional_content(grid, _displace(state, phi, -step))
    _assert_match(delta_content, (plus - minus) / (2 * step), "child content")
    mean, delta_mean = discovery.child_proper_mean_r(grid, state, radius)
    plus = discovery.child_proper_mean_r(grid, _displace(state, radius, step))
    minus = discovery.child_proper_mean_r(grid, _displace(state, radius, -step))
    _assert_match(delta_mean, (plus - minus) / (2 * step), "proper mean r")
    tau_dot, delta_tau_dot, delta_tau = discovery.child_clock_rate(
        grid, state, radius, coordinate_duration=0.2,
    )
    assert delta_tau == pytest.approx(delta_tau_dot * 0.2)

    def clock(sign):
        return discovery.child_clock_rate(grid, _displace(state, radius, sign * step))

    _assert_match(delta_tau_dot, (clock(1) - clock(-1)) / (2 * step), "clock rate")
    report = discovery.galerkin_readout(grid, state, phi, coordinate_duration=0.2)
    assert report["primary_readout"] == "child_regional_content"
    assert report["secondary_readout"] == "child_proper_mean_r"
    assert report["evolved_clock_integral"] is False
    assert report["delta_child_regional_content_tau"] == pytest.approx(
        discovery.matched_tau_correction(
            report["delta_child_regional_content_t"],
            report["child_regional_content_dot"],
            report["delta_tau"],
            report["tau_dot"],
        )
    )
    assert tau_dot > 0
    assert mean > 0

    def observable(eps, tau_target=1.0):
        rate = 2.0 + 0.5 * eps
        return 3.0 * eps + (1.0 + 0.2 * eps) * (tau_target / rate)

    numeric = (observable(step) - observable(-step)) / (2 * step)
    analytic = discovery.matched_tau_correction(3.1, 1.0, 0.25, 2.0)
    assert analytic == pytest.approx(2.975)
    _assert_match(analytic, numeric, "matched tau")


def test_full_gaussian_connected_term_is_nonzero_and_stress_stays_open():
    _grid, state = _manufactured_state()
    weights = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25])
    covariance = discovery.gaussian_covariance(state.phi0, state.phi1, weights)
    columns = np.vstack((state.phi0, state.phi1))
    first = columns[:, 0]
    vertex = np.outer(first, first.conj())
    value = discovery.connected_bilinear(covariance, vertex, vertex)
    assert value == pytest.approx(influence.connected(covariance, vertex, vertex))
    assert value.real == pytest.approx(weights[0] * (1.0 - weights[0]))
    assert abs(value) > 1e-8
    mixed0 = (columns[:, 2] + columns[:, 0]) / np.sqrt(2)
    mixed1 = (columns[:, 3] + columns[:, 1]) / np.sqrt(2)
    modal = discovery.child_modal_kernel(covariance, np.column_stack((mixed0, mixed1)))
    assert modal["kernel"].shape == (2, 2)
    assert modal["spatial_outside"] is False
    assert modal["spatial_outside_by_definition"] is False
    assert modal["initial_cross_is_spatial_exterior"] is False
    assert modal["initial_cross_frobenius"] > 1e-6
    stress = discovery.future_effective_stress_specification(
        covariance, np.column_stack((mixed0, mixed1)),
    )
    assert stress["implemented"] is False
    assert stress["effective_stress"] is None
    assert stress["parts"]["complementary_force"]["implemented"] is False
    assert stress["parts"]["initial_cross"]["implemented"] is True
    assert stress["parts"]["retarded_delta_state"]["implemented"] is False
    assert stress["parts"]["retarded_delta_state"]["value"] is None
    assert stress["induced_force_included"] is False
    assert stress["kernel_2x2_is_spatial_outside"] is False
    assert stress["finite_difference_used_as_analytic"] is False


def _saved_initial_state():
    with np.load(_V5, allow_pickle=False) as data:
        fermions = int(data["nf"])
        quadrature = int(data["nq"])
        length = float(data["length"])
        grid = galerkin.build_grid(fermions, quadrature=quadrature, length=length, gauge="conformal")
        state = coupling.CauchyState(
            Q=np.array(data["geometry_Q"], dtype=float, copy=True),
            r=np.array(data["geometry_r"], dtype=float, copy=True),
            chi=np.array(data["geometry_chi"], dtype=float, copy=True),
            p_Q=np.array(data["geometry_p_Q"], dtype=float, copy=True),
            p_r=np.array(data["geometry_p_r"], dtype=float, copy=True),
            p_chi=np.array(data["geometry_p_chi"], dtype=float, copy=True),
            phi0=np.array(data["columns_phi0"], copy=True),
            phi1=np.array(data["columns_phi1"], copy=True),
        )
        modes = np.array(data["geometry_modes"], dtype=float, copy=True)
        occupations = np.array(data["occupations"], dtype=float, copy=True)
    if not np.array_equal(modes, grid.modes_g):
        raise AssertionError("saved geometry modes differ from the owned grid")
    if not np.allclose(occupations, grid.fine.occupations):
        raise AssertionError("saved occupations differ from the live system")
    return grid, state


def test_saved_initial_state_selected_directions_match_centred_differences():
    grid, state = _saved_initial_state()
    assert np.min(state.r) > 0 and np.min(state.Q) > 0
    for name in ("Q", "phi", "occupations"):
        tangent = _direction(grid, name, seed=11)
        analytic = discovery.rate_jacobian_vector(grid, state, tangent)
        numeric = _centred_galerkin(grid, state, tangent)
        for field in _RATE_FIELDS:
            _assert_match(getattr(analytic, field), numeric[field], f"saved {name}.{field}")
    phi = _direction(grid, "phi", seed=11)
    prepared = discovery.prepared_radius_tangent(grid, state, phi)
    assert prepared.available is True and prepared.newton_used is False
    assert prepared.linear_residual_max < 1e-8

    def saved_residual(radius, column0, column1):
        trial = state.copy()
        trial.r = radius
        trial.phi0 = column0
        trial.phi1 = column1
        fine = galerkin.prolong_state(grid, trial)
        source = coupling.source_from_columns(galerkin.active_fine_system(grid, fine), fine)
        projected, _full = galerkin.projected_radius_operator(
            grid, radius, source["force_L"] / grid.dx_q,
        )
        return projected

    step = 1e-6
    source_only = (
        saved_residual(state.r, state.phi0 + step * phi.phi0, state.phi1 + step * phi.phi1)
        - saved_residual(state.r, state.phi0 - step * phi.phi0, state.phi1 - step * phi.phi1)
    ) / (2 * step)
    joint = (
        saved_residual(
            state.r + step * prepared.delta_r,
            state.phi0 + step * phi.phi0, state.phi1 + step * phi.phi1,
        )
        - saved_residual(
            state.r - step * prepared.delta_r,
            state.phi0 - step * phi.phi0, state.phi1 - step * phi.phi1,
        )
    ) / (2 * step)
    source_scale = max(1e-8, float(np.max(np.abs(source_only))))
    assert float(np.max(np.abs(joint))) <= 2e-5 * source_scale + 2e-8
    content = discovery.child_regional_content(grid, state)
    mean = discovery.child_proper_mean_r(grid, state)
    assert 0 < content < float(np.sum(grid.fine.occupations))
    assert mean > 0
    columns = np.vstack((state.phi0, state.phi1))
    covariance = discovery.gaussian_covariance(state.phi0, state.phi1, grid.fine.occupations)
    child = columns[:, list(nested.CHILD_INDICES)]
    modal = discovery.child_modal_kernel(covariance, child)
    assert modal["kernel"].shape == (2, 2)
    assert modal["spatial_outside_by_definition"] is False
    stress = discovery.future_effective_stress_specification(covariance, child)
    assert stress["effective_stress"] is None
    assert stress["missing_primitive"] == "ctp_future_effective_stress"
