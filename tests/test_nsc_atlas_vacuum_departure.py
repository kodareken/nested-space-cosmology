"""Exact and manufactured algebra; no physical source or evolution is run."""
import copy
import time

import numpy as np
import pytest
from scipy.integrate import cumulative_simpson, cumulative_trapezoid

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_exterior_characteristic_slice import exterior_characteristic_slice
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator
from recursive_horizons.nsc_atlas_vacuum_departure import departure_slice, vacuum_slice


TUBE = {"r": (.2, .9), "L": (-3., 3.), "c": (-4., 4.)}


@pytest.fixture(scope="module")
def shared_budget():
    budget = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield budget
    assert budget.evaluations < 20000
    assert time.process_time()-budget.started < 30.
    print("vacuum departure test charges", budget.evaluations,
          "process CPU", time.process_time()-budget.started)


def _vacuum_data(budget, nodes=17, scalar=True):
    U = np.linspace(0., .4, nodes)
    x = .45
    base = vacuum_slice(U, x, budget=budget)
    c = .2+U+.3*U*U if scalar else np.zeros(nodes)
    states = {"r": base["r"].copy(), "L": base["L"].copy(), "c": c}
    jets = {"r_Ubar": np.zeros(nodes), "log_C_Ubar": np.zeros(nodes),
            "c_Ubar": 1+.6*U if scalar else np.zeros(nodes)}
    source = {"r": x, "log_C": 0., "c": float(c[0]),
              "r_x": 1., "log_C_x": 0., "c_x": -.3 if scalar else 0.}
    return U, states, jets, source, x


def _solve(data, budget, *, alpha=0., order=2):
    U, states, jets, source, x = data
    return departure_slice(U, states, jets, source, x=x, alpha=alpha,
                           tube=TUBE, budget=budget, quadrature_order=order)


@pytest.mark.parametrize("order", [2, 4])
@pytest.mark.parametrize("alpha,scalar", [(0., True), (0., False), (.7, False)])
def test_vacuum_geometry_preserved_with_active_test_field_or_quiet_finite_field(shared_budget, order, alpha, scalar):
    data = _vacuum_data(shared_budget, scalar=scalar)
    row = _solve(data, shared_budget, alpha=alpha, order=order)
    for name, base_name in (("r_x", "r_x"), ("log_C_x", "log_C_x"),
                            ("r_Ubar", "r_Ubar"), ("log_C_Ubar", "log_C_Ubar")):
        np.testing.assert_array_equal(row[name], row["vacuum"][base_name])
    np.testing.assert_array_equal(row["delta_r_x"], 0.)
    np.testing.assert_array_equal(row["delta_log_C_x"], 0.)
    if scalar:
        assert np.max(abs(row["f_Ubar"])) > .5
        assert np.max(abs(row["c_Ubar_x"])) > .01
    else:
        np.testing.assert_array_equal(row["c_x"], 0.)
    # A nonzero old quadrature residual is recorded, not used for evolution.
    assert np.max(abs(row["removed_radial_quadrature_error"])) > 1e-11
    assert row["scope"]["vacuum_quadrature_drives_evolution"] is False


def test_exact_implicit_radius_jets_mass_and_mixed_equations(shared_budget):
    U, x = np.array([0., .1, .3, .7]), np.array([.4, .43, .46, .49])
    v = vacuum_slice(U, x, budget=shared_budget)
    r, p, R, C = (v[key] for key in ("r", "r_Ubar", "r_x", "C"))
    a = .8
    F = lambda z: (1-z/a)*np.exp(z/a)
    Fp = lambda z: -z*np.exp(z/a)/(a*a)
    Fpp = lambda z: -(1+z/a)*np.exp(z/a)/(a*a)
    np.testing.assert_allclose(F(r), (1-U/(4*a))*F(x), atol=3e-16, rtol=0)
    np.testing.assert_allclose(Fp(r)*p, -F(x)/(4*a), rtol=3e-15)
    np.testing.assert_allclose(Fp(r)*R, (1-U/(4*a))*Fp(x), rtol=3e-15)
    np.testing.assert_allclose(Fpp(r)*p*p+Fp(r)*v["r_Ubar_Ubar"], 0., atol=3e-16)
    np.testing.assert_allclose(Fpp(r)*R*R+Fp(r)*v["r_xx"], (1-U/(4*a))*Fpp(x), rtol=3e-15)
    np.testing.assert_allclose(Fpp(r)*p*R+Fp(r)*v["r_Ubar_x"], -Fp(x)/(4*a), rtol=3e-15)
    np.testing.assert_allclose(r*(1+4*p*R/C)/2, .4, atol=3e-16, rtol=0)
    rhs = bulk_operator({"r": r, "C": C, "c": np.zeros_like(r)},
                        {"r_U": p, "r_V": R, "c_U": 0., "c_V": 0.}, alpha=.3, region="exterior")
    np.testing.assert_allclose(v["r_Ubar_x"], rhs["r_UV"], atol=3e-16, rtol=2e-15)
    np.testing.assert_allclose(v["log_C_Ubar_x"], rhs["log_C_UV"], atol=3e-16, rtol=2e-15)


def test_lapse_pure_jets_from_independent_directional_changes(shared_budget):
    U, x, h = .27, .46, 2e-5
    v = vacuum_slice(U, x, budget=shared_budget)
    for coordinate, suffix in (("U", "Ubar"), ("x", "x")):
        plus = vacuum_slice(U+h if coordinate == "U" else U,
                            x+h if coordinate == "x" else x, budget=shared_budget)
        minus = vacuum_slice(U-h if coordinate == "U" else U,
                             x-h if coordinate == "x" else x, budget=shared_budget)
        half_plus = vacuum_slice(U+h/2 if coordinate == "U" else U,
                                 x+h/2 if coordinate == "x" else x, budget=shared_budget)
        half_minus = vacuum_slice(U-h/2 if coordinate == "U" else U,
                                  x-h/2 if coordinate == "x" else x, budget=shared_budget)
        first = "log_C_"+suffix
        second = "log_C_Ubar_Ubar" if coordinate == "U" else "log_C_xx"
        assert float((plus["L"]-minus["L"])/(2*h)) == pytest.approx(float(v[first]), abs=3e-9)
        derivative = (4*(half_plus[first]-half_minus[first])/h-(plus[first]-minus[first])/(2*h))/3
        assert float(derivative) == pytest.approx(float(v[second]), abs=1e-8)


def test_exact_source_anchor_and_beta_zero_are_regular(shared_budget):
    x = .45
    r_cross = (4/25)**(1/3)
    F = lambda r: (1-r/.8)*np.exp(r/.8)
    U_cross = 3.2*(1-F(r_cross)/F(x))
    v = vacuum_slice([0., U_cross], x, budget=shared_budget)
    assert v["r"][0] == x and v["C"][0] == 1 and v["L"][0] == 0
    assert v["r_x"][0] == 1 and v["log_C_x"][0] == 0
    assert v["r"][1] == pytest.approx(r_cross, abs=2e-16)
    assert abs(1/(5*v["r"][1]**2)-5*v["r"][1]/4) < 6e-16
    assert all(np.all(np.isfinite(value)) for value in v.values())


def _manufactured(nodes, budget, alpha=.3):
    # Constant physical r,C and linear canonical f give explicit polynomial
    # physical prefixes; the nonzero geometric departure is not a solution.
    U, x, r, C = np.linspace(0., .4, nodes), .45, .65, 1.3
    v = vacuum_slice(U, x, budget=budget)
    f0, b, R0, S0, L0 = .3, 1.2, .8, -.4, .15
    f = f0+b*U
    scale = np.sqrt(alpha/6)
    u = np.sinh(scale*f)/scale
    root_D = np.cosh(scale*f)
    states = {"r": np.full(nodes, r), "L": np.full(nodes, np.log(C)), "c": r*u}
    jets = {"r_Ubar": -v["r_Ubar"], "log_C_Ubar": -v["log_C_Ubar"],
            "c_Ubar": r*root_D*b}
    source = {"r": r, "log_C": float(np.log(C)), "c": float(states["c"][0]),
              "r_x": R0, "log_C_x": L0, "c_x": float(u[0]*R0+r*root_D[0]*S0)}
    R = R0-C*U/(4*r)
    S = S0-b*(R0*U-C*U*U/(8*r))/r
    exact = {"r_x": R, "f_x": S, "c_x": u*R+r*root_D*S,
             "log_C_x": L0+C*U/(2*r*r)-alpha*b*(S0*U-b*R0*U*U/(2*r)+b*C*U**3/(24*r*r))}
    return (U, states, jets, source, x), exact


@pytest.mark.parametrize("order,ratio", [(2, 3.7), (4, 12.)])
def test_finite_departures_converge_to_independent_polynomial_solution(shared_budget, order, ratio):
    errors, old_differences = [], []
    for nodes in (9, 17, 33):
        data, exact = _manufactured(nodes, shared_budget)
        row = _solve(data, shared_budget, alpha=.3, order=order)
        physical_jets = {key: row[key] for key in ("r_Ubar", "log_C_Ubar", "c_Ubar")}
        old = exterior_characteristic_slice(data[0], data[1], physical_jets,
            {key: data[3][key] for key in ("r_x", "log_C_x", "c_x")},
            alpha=.3, quadrature_order=order, tube=TUBE, budget=shared_budget)
        errors.append(max(float(np.max(abs(row[key]-value))) for key, value in exact.items()))
        old_differences.append(max(float(np.max(abs(row[key]-old[key]))) for key in exact))
        np.testing.assert_allclose(row["r_x"]-old["r_x"],
                                   row["removed_radial_quadrature_error"], atol=6e-16, rtol=1e-8)
    print("departure manufactured", order, "errors", errors, "old differences", old_differences)
    assert np.all(np.array(errors[:-1])/errors[1:] > ratio)
    assert np.all(np.array(old_differences[:-1])/old_differences[1:] > ratio)


def test_subtracted_lapse_algebra_and_full_kinetic_equation(shared_budget):
    data, _ = _manufactured(17, shared_budget)
    row = _solve(data, shared_budget, alpha=.3, order=4)
    U, states = data[:2]
    rhs = bulk_operator({"r": states["r"], "C": np.exp(states["L"]), "c": states["c"]},
        {"r_U": row["r_Ubar"], "r_V": row["r_x"],
         "c_U": row["c_Ubar"], "c_V": row["c_x"]}, alpha=.3, region="exterior")
    for key, other in (("r_Ubar_x", "r_UV"), ("log_C_Ubar_x", "log_C_UV"), ("c_Ubar_x", "c_UV")):
        np.testing.assert_allclose(row[key], rhs[other], rtol=3e-15, atol=4e-16)
    np.testing.assert_allclose(row["lapse_difference_integrand"],
        rhs["log_C_UV"]-row["vacuum"]["log_C_Ubar_x"], rtol=3e-15, atol=8e-16)
    ordinary_lapse = data[3]["log_C_x"]+cumulative_simpson(rhs["log_C_UV"], x=U, initial=0.)
    np.testing.assert_allclose(row["log_C_x"]-ordinary_lapse,
        row["removed_lapse_quadrature_error"], atol=9e-16, rtol=1e-8)
    assert np.max(abs(.3*row["f_Ubar"]*row["f_x"])) > .1


def test_original_source_and_nonzero_source_offsets_are_preserved(shared_budget):
    data, _ = _manufactured(5, shared_budget)
    data[3]["h"] = 99.  # Extra snapshot fields are accepted, never substituted.
    original = copy.deepcopy(data)
    row = _solve(data, shared_budget, alpha=.3)
    for name in ("r_x", "log_C_x", "c_x"):
        assert row[name][0] == data[3][name]
    assert row["delta_r"][0] == data[3]["r"]-data[4]
    assert row["delta_log_C"][0] == data[3]["log_C"]
    for group, before in zip(data[1:4], original[1:4]):
        for name in group:
            np.testing.assert_array_equal(group[name], before[name])
    row["r_x"][0] = 99
    assert data[3]["r_x"] == .8


def test_zero_source_trace_does_not_erase_transverse_field_memory(shared_budget):
    data = _vacuum_data(shared_budget, nodes=9)
    U, states, jets, source, x = data
    states["c"] = U+.3*U*U
    source["c"], source["c_x"] = 0., 0.
    row = _solve(data, shared_budget)
    assert row["c_x"][0] == 0.
    assert row["c_Ubar"][0] == 1.
    assert row["f_Ubar"][0] == 1/x
    assert np.max(abs(row["c_x"][1:])) > 1e-4


def test_declared_counts_and_zero_budget_reset(shared_budget):
    data = _vacuum_data(shared_budget, nodes=5)
    before, started = shared_budget.evaluations, shared_budget.started
    row = _solve(data, shared_budget)
    assert shared_budget.evaluations-before == row["settings"]["charged_evaluations"] == 35
    assert shared_budget.started == started
    low = Budget(maximum_evaluations=8, cpu_seconds=30., evaluations=2)
    clock = low.started
    with pytest.raises(RuntimeError, match="budget exhausted"):
        _solve(data, low)
    assert low.evaluations == 12 and low.started == clock


@pytest.mark.parametrize("kind", ["missing_source", "source_shape", "source_mismatch", "source_boolean",
    "not_uniform", "not_anchored", "too_few", "too_many", "x_shape", "x_domain", "U_domain",
    "missing_jet", "state_shape", "state_boolean", "nonfinite", "tube", "alpha", "order", "budget"])
def test_invalid_departure_inputs_stop_before_charging(shared_budget, kind):
    data = list(_vacuum_data(shared_budget, nodes=5))
    U, states, jets, source, x = data
    kwargs = dict(x=x, alpha=0., tube=copy.deepcopy(TUBE), budget=shared_budget, quadrature_order=2)
    if kind == "missing_source":
        del source["r_x"]
    elif kind == "source_shape":
        source["c_x"] = [1.]
    elif kind == "source_mismatch":
        source["r"] += 1e-12
    elif kind == "source_boolean":
        source["r_x"] = True
    elif kind == "not_uniform":
        U[1] += .001
    elif kind == "not_anchored":
        U += .01
    elif kind == "too_few":
        U = U[:2]
    elif kind == "too_many":
        U = np.linspace(0., .4, 514)
    elif kind == "x_shape":
        kwargs["x"] = [x]
    elif kind == "x_domain":
        kwargs["x"] = .8
    elif kind == "U_domain":
        U *= 10
    elif kind == "missing_jet":
        del jets["r_Ubar"]
    elif kind == "state_shape":
        states["r"] = .45
    elif kind == "state_boolean":
        states["L"] = [True, 0., 0., 0., 0.]
    elif kind == "nonfinite":
        jets["c_Ubar"][1] = np.nan
    elif kind == "tube":
        kwargs["tube"]["r"] = (.7, .9)
    elif kind == "alpha":
        kwargs["alpha"] = -.1
    elif kind == "order":
        kwargs["quadrature_order"] = True
    else:
        kwargs["budget"] = None
    before = shared_budget.evaluations
    with pytest.raises(ValueError):
        departure_slice(U, states, jets, source, **kwargs)
    assert shared_budget.evaluations == before


@pytest.mark.parametrize("U,x", [(True, .45), (0., True), (-.1, .45), (3.2, .45),
                                  (0., 0.), (0., .8), ([[0.]], .45), ([0., .1], [.4, .5, .6])])
def test_invalid_vacuum_coordinates(shared_budget, U, x):
    before = shared_budget.evaluations
    with pytest.raises(ValueError):
        vacuum_slice(U, x, budget=shared_budget)
    assert shared_budget.evaluations == before
