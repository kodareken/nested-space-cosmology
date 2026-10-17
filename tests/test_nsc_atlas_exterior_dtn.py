"""Retarded exterior Dirichlet map: kernel, corners, normal, and l=0 recovery."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_exterior_dtn import (
    Q0, contraction_bound, kruskal_factor, mode_kernel, outward_normal_derivative,
    proper_null_speeds, radius_null_derivatives, solve_exterior_dirichlet,
)
from recursive_horizons.nsc_atlas_field_preparation import (
    M_OVER_ELL, RADIAL_FLOOR, WALL_UPPER, coefficient_from_radius_derivatives,
    exterior_coefficient, preparation_record, radius_from_null, radius_product, source_profile,
)
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _geometry, _solve


def _quadratic_wall(Va, Vb):
    span = Vb-Va

    def trace(V):
        coordinate = np.asarray(V, dtype=float)
        slope = (coordinate-Va)/span
        return slope**2, 2*slope/span

    return trace


def test_mode_kernel_comes_from_radius_curvature_and_decreases():
    radii = np.linspace(RADIAL_FLOOR, WALL_UPPER, 9)
    curvature = 16*M_OVER_ELL**4*np.exp(-radii/(2*M_OVER_ELL))/radii**4
    centrifugal = 8*M_OVER_ELL**3*np.exp(-radii/(2*M_OVER_ELL))/radii**3
    assert np.all(np.diff(curvature) < 0)
    assert np.all(np.diff(centrifugal) < 0)
    base = coefficient_from_radius_derivatives(radii)
    for degree in range(4):
        expected = base+kruskal_factor(radii)*degree*(degree+1)/(4*radii**2)
        assert np.max(np.abs(mode_kernel(radii, degree)-expected)) < 1e-12
    floor = float(exterior_coefficient(RADIAL_FLOOR))
    for degree in range(4):
        sharp = 1+.5*degree*(degree+1)
        loose = 1+degree*(degree+1)*WALL_UPPER/(2*M_OVER_ELL)
        assert mode_kernel(RADIAL_FLOOR, degree) == pytest.approx(floor*sharp)
        assert contraction_bound(degree)["q_l"] == pytest.approx(Q0*sharp)
        assert contraction_bound(degree)["q_l"] < 1
        if degree:
            assert sharp < loose
    assert Q0*(1+.5*4*5) > 1
    with pytest.raises(ValueError, match="0 through 3"):
        solve_exterior_dirichlet(_quadratic_wall(-2.5, -2.3), angular_degree=4)


def test_constant_phi_has_zero_outward_derivative_and_speeds_share_one_slope():
    geometry = _geometry(129)
    V = .5*(geometry["Va"]+geometry["Vb"])
    wall = float(geometry["wall"](V))
    prime = float(geometry["wall"](V, 1))
    radius = float(radius_from_null(wall, V))
    radius_U, radius_V = radius_null_derivatives(wall, V, radius)
    step = 1e-6
    numerical_U = float(radius_from_null(wall+step, V)-radius_from_null(wall-step, V))/(2*step)
    numerical_V = float(radius_from_null(wall, V+step)-radius_from_null(wall, V-step))/(2*step)
    assert radius_U == pytest.approx(numerical_U, rel=1e-8)
    assert radius_V == pytest.approx(numerical_V, rel=1e-8)
    udot, vdot = proper_null_speeds(prime, radius)
    assert udot*vdot*float(kruskal_factor(radius)) == pytest.approx(1.)
    assert udot/vdot == pytest.approx(prime)
    scaled_udot, scaled_vdot = proper_null_speeds(4*prime, radius)
    assert scaled_udot == pytest.approx(2*udot)
    assert scaled_vdot == pytest.approx(vdot/2)
    phi = 2.5
    normal = outward_normal_derivative(phi*radius_U, phi*radius_V, phi*radius, wall, V, prime)
    assert float(normal["phi"]) == pytest.approx(phi)
    assert float(normal["outward_dphi"]) == pytest.approx(0., abs=1e-12)
    omitted = udot*(phi*radius_U/radius)-vdot*(phi*radius_V/radius)
    assert omitted != pytest.approx(0., abs=1e-8)


def test_corners_reject_incompatible_traces_and_zero_kernel_is_exact(monkeypatch):
    geometry = _geometry(129)
    Va, Vb = geometry["Va"], geometry["Vb"]
    budget = PreparedFieldBudget()
    with pytest.raises(ValueError, match="U0"):
        solve_exterior_dirichlet(
            _quadratic_wall(Va, Vb),
            source=lambda V: (np.ones(np.shape(V)), np.zeros(np.shape(V))),
            budget=budget)
    with pytest.raises(ValueError, match="Ua"):
        solve_exterior_dirichlet(
            lambda V: (np.ones(np.shape(V)), np.zeros(np.shape(V))), budget=budget)
    with pytest.raises(ValueError, match="C1"):
        solve_exterior_dirichlet(
            lambda V: (np.asarray(V, dtype=float)-Va, np.ones(np.shape(V))), budget=budget)
    import recursive_horizons.nsc_atlas_exterior_dtn as dtn
    cached = _geometry(129)

    def geometry_without_independent_vdot(knots):
        if knots != 129:
            raise AssertionError("unexpected geometry knots")
        curves = dict(cached)
        curves["Vdot"] = lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("independent Vdot"))
        return curves

    monkeypatch.setattr(dtn, "_geometry", geometry_without_independent_vdot)
    monkeypatch.setattr(
        dtn, "mode_kernel",
        lambda radius, angular_degree=0: np.zeros(np.shape(radius), dtype=float))
    wall = _quadratic_wall(Va, Vb)
    solution = solve_exterior_dirichlet(wall, budget=budget)
    sample = solution.boundary(np.linspace(Va, Vb, 9))
    prime = np.array([float(geometry["wall"](value, 1)) for value in sample["V"]])
    _, derivative = wall(sample["V"])
    assert sample["wall_integral"] == pytest.approx(0., abs=1e-12)
    assert sample["chi_V"] == pytest.approx(0., abs=1e-12)
    assert sample["chi_U"] == pytest.approx(derivative/prime, rel=0, abs=1e-9)


def test_linear_response_of_one_prescribed_wall():
    geometry = _geometry(129)
    Va, Vb = geometry["Va"], geometry["Vb"]
    wall = _quadratic_wall(Va, Vb)

    def doubled(V):
        value, derivative = wall(V)
        return 2*value, 2*derivative

    budget = PreparedFieldBudget()
    first = solve_exterior_dirichlet(wall, angular_degree=1, budget=budget)
    second = solve_exterior_dirichlet(doubled, angular_degree=1, budget=budget)
    left = first.boundary(np.linspace(Va, Vb, 7))
    right = second.boundary(left["V"])
    assert left["phi"] == pytest.approx(.5*right["phi"], rel=1e-6, abs=1e-8)
    assert left["outward_dphi"] == pytest.approx(.5*right["outward_dphi"], rel=1e-6, abs=1e-8)
    assert budget.kernel_values < budget.kernel_values_max


def test_prescribed_prepared_wall_recovers_vanishing_chi_U():
    geometry = _geometry(129)
    prepared = preparation_record()
    for key in ("Va", "Vb", "Ua", "Ub"):
        assert geometry[key] == pytest.approx(prepared["geometry"][key], abs=1e-12)
    budget = PreparedFieldBudget()
    extract, _, _, _ = _solve(129, geometry, prepared, budget, include_dense=True)
    Va, Vb = geometry["Va"], geometry["Vb"]

    def wall_trace(V):
        coordinate = np.asarray(V, dtype=float)
        row = extract((coordinate-Va)/(Vb-Va))
        return (np.reshape(row["wall_chi_over_epsilon"], coordinate.shape),
                np.reshape(row["wall_dchi_dV_over_epsilon"], coordinate.shape))

    solution = solve_exterior_dirichlet(
        wall_trace, source=lambda V: source_profile(V, Va, Vb), bottom=None, budget=budget)
    sample = solution.boundary(np.linspace(Va, Vb, 17))
    reference = extract((sample["V"]-Va)/(Vb-Va))
    assert sample["chi_U"] == pytest.approx(0., abs=1e-6)
    assert sample["chi_V"] == pytest.approx(reference["wall_dchi_dV_over_epsilon"], rel=1e-5, abs=1e-6)
    assert budget.kernel_values <= 30_000_000


def test_constant_physical_field_and_prescribed_node_entry_values():
    # phi=1 is an exact l=0 solution in the Ricci-flat exterior. This checks
    # the full curved solver, independent of the old prepared-pulse solution.
    from scipy.optimize import brentq
    geometry = _geometry(129)
    Va, Vb, Ua = (geometry[key] for key in ('Va', 'Vb', 'Ua'))
    U0 = float(radius_product(RADIAL_FLOOR))/Va

    def wall(V):
        U = geometry['wall'](V)
        radius = radius_from_null(U, V)
        rU, rV = radius_null_derivatives(U, V, radius)
        return radius, rU*geometry['wall'](V, 1)+rV

    def source(V):
        radius = radius_from_null(U0, V)
        return radius, radius_null_derivatives(U0, V, radius)[1]

    def bottom(U):
        radius = radius_from_null(U, Va)
        return radius, radius_null_derivatives(U, Va, radius)[0]

    budget = PreparedFieldBudget(cpu_seconds=30, kernel_values_max=10_000_000)
    solution = solve_exterior_dirichlet(wall, source=source, bottom=bottom, budget=budget)
    sample = solution.boundary(np.linspace(Va, Vb, 17))
    assert sample['phi'] == pytest.approx(1., abs=2e-15)
    # One fixed 129-node numerical control; not a continuum error enclosure.
    assert sample['outward_dphi'] == pytest.approx(0., abs=1e-5)
    entering = np.flatnonzero(solution.grid > Ua)
    for index in entering[::max(1, len(entering)//3)]:
        coordinate = solution.grid[index]
        entry = brentq(lambda v: geometry['wall'](v)-coordinate, Va, Vb, xtol=1e-14)
        state = solution._state(np.array([(entry-Va)/(Vb-Va)]))[index, 0]
        prescribed = wall(entry)[0]-source(entry)[0]
        assert state == pytest.approx(prescribed, rel=0, abs=1e-13)
