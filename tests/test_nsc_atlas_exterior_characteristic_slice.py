"""Manufactured and exact-vacuum controls; no prepared source is constructed."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_exterior_characteristic_slice import exterior_characteristic_slice
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator


TUBE = {"r": (.2, 1.), "L": (-3., 3.), "c": (-4., 4.)}


@pytest.fixture(scope="module")
def shared_budget():
    budget = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield budget
    assert budget.evaluations <= 20000
    assert time.process_time()-budget.started <= 30.


def _minkowski(nodes=17, u=.7):
    U = np.linspace(0., .15, nodes)
    r = .6-U/2
    states = {"r": r, "L": np.zeros(nodes), "c": u*r}
    jets = {"r_Ubar": np.full(nodes, -.5), "log_C_Ubar": np.zeros(nodes),
            "c_Ubar": np.full(nodes, -u/2)}
    source = {"r_x": .5, "log_C_x": 0., "c_x": u/2}
    return U, states, jets, source


def _solve(data, budget, alpha=.2, order=4, tube=TUBE):
    return exterior_characteristic_slice(*data, alpha=alpha, quadrature_order=order,
                                         tube=tube, budget=budget)


@pytest.mark.parametrize("alpha", [0., .2, 2.])
@pytest.mark.parametrize("order", [2, 4])
def test_constant_canonical_field_in_flat_geometry(shared_budget, alpha, order):
    data = _minkowski()
    row = _solve(data, shared_budget, alpha, order)
    np.testing.assert_allclose(row["r_x"], .5, atol=2e-16, rtol=0)
    np.testing.assert_allclose(row["c_x"], .35, atol=2e-16, rtol=0)
    np.testing.assert_allclose(row["log_C_x"], 0., atol=2e-16, rtol=0)
    np.testing.assert_allclose(row["f_Ubar"], 0., atol=2e-16, rtol=0)
    assert row["scope"] == {"source_integrations": 0, "time_steps": 0,
                            "pure_constraints_evaluated": False,
                            "quadrature_error_is_enclosed": False}


def _canonical_manufactured(nodes, alpha):
    # An arbitrary supplied slice, not a claim that its pure constraints hold.
    U = np.linspace(0., .2, nodes)
    r, C, f0, b, R0, S0, L0 = .65, 1.3, .3, 1.2, .8, -.4, .15
    f = f0+b*U
    if alpha:
        q = np.sqrt(alpha/6)*f
        u, root_D = np.sinh(q)/np.sqrt(alpha/6), np.cosh(q)
    else:
        u, root_D = f, np.ones_like(f)
    states = {"r": np.full(nodes, r), "L": np.full(nodes, np.log(C)), "c": r*u}
    jets = {"r_Ubar": np.zeros(nodes), "log_C_Ubar": np.zeros(nodes), "c_Ubar": r*root_D*b}
    source = {"r_x": R0, "log_C_x": L0, "c_x": u[0]*R0+r*root_D[0]*S0}
    R = R0-C*U/(4*r)
    S = S0-b*(R0*U-C*U*U/(8*r))/r
    exact = {"r_x": R, "f_x": S, "c_x": u*R+r*root_D*S,
             "log_C_x": L0+C*U/(2*r*r)-alpha*b*(S0*U-b*R0*U*U/(2*r)+b*C*U**3/(24*r*r))}
    return (U, states, jets, source), exact


@pytest.mark.parametrize("alpha", [0., .3])
def test_nonzero_canonical_profile_has_exact_polynomial_prefixes(shared_budget, alpha):
    data, exact = _canonical_manufactured(17, alpha)
    row = _solve(data, shared_budget, alpha)
    for key, value in exact.items():
        np.testing.assert_allclose(row[key], value, rtol=1e-14, atol=4e-16)
    if alpha:
        # Detect omission or reversal of the scalar kinetic lapse term.
        scalar_free_lapse = data[3]["log_C_x"]+1.3*data[0]/(2*.65**2)
        assert abs(row["log_C_x"][-1]-scalar_free_lapse[-1]) > .03


def test_pointwise_rhs_agrees_with_independent_mixed_operator(shared_budget):
    data, _ = _canonical_manufactured(19, .7)
    row = _solve(data, shared_budget, .7, order=2)
    _, states, jets, _ = data
    rhs = bulk_operator({"r": states["r"], "C": np.exp(states["L"]), "c": states["c"]},
                        {"r_U": jets["r_Ubar"], "r_V": row["r_x"],
                         "c_U": jets["c_Ubar"], "c_V": row["c_x"]},
                        alpha=.7, region="exterior")
    for name, other in (("r_Ubar_x", "r_UV"), ("log_C_Ubar_x", "log_C_UV"), ("c_Ubar_x", "c_UV")):
        np.testing.assert_allclose(row[name], rhs[other], rtol=2e-15, atol=2e-15)


def test_affine_spatial_system_and_canonical_reduction_agree_for_arbitrary_jets():
    # Independent pointwise chain rule, before any prefix approximation.
    r, C, u, p, h, R, q, alpha = .51, 1.2, .8, .31, -.27, .77, .19, .4
    D = 1+alpha*u*u/6
    z = (h-u*p)/r
    M = np.array([[-p/r, 0, 0],
                  [2*p/r**2+alpha*u*z/(D*r), 0, -alpha*z/(D*r)],
                  [-u*p/r-alpha*u*u*z/(6*D), 0, alpha*u*z/(6*D)]])
    b = np.array([-C/(4*r), C/(2*r*r), -u*C/(4*r)])
    affine = M@np.array([R, -9., q])+b
    rhs = bulk_operator({"r": r, "C": C, "c": r*u},
                        {"r_U": p, "r_V": R, "c_U": h, "c_V": q}, alpha=alpha, region="exterior")
    np.testing.assert_allclose(affine, [rhs["r_UV"], rhs["log_C_UV"], rhs["c_UV"]], rtol=2e-15)
    fU, fx = z/np.sqrt(D), (q-u*R)/(r*np.sqrt(D))
    ux = (q-u*R)/r
    uUx = (rhs["c_UV"]-ux*p-u*rhs["r_UV"]-z*R)/r
    fUx = uUx/np.sqrt(D)-alpha*u*z*ux/(6*D**1.5)
    assert abs(r*fUx+p*fx+R*fU) < 2e-15


def _schwarzschild(nodes, u=.8):
    U = np.linspace(0., .4, nodes)
    x, a = .45, .8
    product = (1-U/(4*a))*(1-x/a)*np.exp(x/a)
    r = a*(1+lambertw(-product/np.e, 0).real)
    L = np.log(x/r)+(x-r)/a
    C = np.exp(L)
    p, R = C*(a-x)/(4*x), (1-U/(4*a))*C
    LU, Lx = -(1/r+1/a)*p, 1/x+1/a-(1/r+1/a)*R
    states = {"r": r, "L": L, "c": u*r}
    jets = {"r_Ubar": p, "log_C_Ubar": LU, "c_Ubar": u*p}
    source = {"r_x": 1., "log_C_x": 0., "c_x": u}
    return (U, states, jets, source), {"r_x": R, "log_C_x": Lx, "c_x": u*R}


@pytest.mark.parametrize("order,minimum_ratio", [(2, 3.8), (4, 13.)])
def test_trapped_schwarzschild_prefix_quadrature_converges(shared_budget, order, minimum_ratio):
    errors = []
    for nodes in (9, 17, 33):
        data, exact = _schwarzschild(nodes)
        row = _solve(data, shared_budget, alpha=.6, order=order)
        errors.append(max(float(np.max(abs(row[key]-value))) for key, value in exact.items()))
        assert np.all(data[2]["r_Ubar"] > 0)
        assert np.all(data[1]["r"] < .8)
    ratios = np.array(errors[:-1])/errors[1:]
    print("Schwarzschild prefix", order, "errors", errors, "ratios", ratios.tolist())
    assert np.all(ratios > minimum_ratio)


def test_affine_outgoing_coordinate_rescaling_and_source_anchor(shared_budget):
    data, _ = _canonical_manufactured(17, .3)
    row = _solve(data, shared_budget, alpha=.3)
    scale = 2.3
    U, states, jets, source = copy.deepcopy(data)
    U *= scale
    states["L"] -= np.log(scale)
    jets = {key: value/scale for key, value in jets.items()}
    other = _solve((U, states, jets, source), shared_budget, alpha=.3)
    for key in ("r_x", "log_C_x", "c_x", "f_x"):
        np.testing.assert_allclose(other[key], row[key], rtol=2e-15, atol=3e-16)
    assert all(row[key][0] == source[key] for key in source)


def test_declared_charge_count_and_input_immutability(shared_budget):
    data = _minkowski(5)
    before = copy.deepcopy(data)
    charged = shared_budget.evaluations
    row = _solve(data, shared_budget)
    assert shared_budget.evaluations-charged == 20 == row["settings"]["charged_evaluations"]
    np.testing.assert_array_equal(data[0], before[0])
    for mapping, original in zip(data[1:], before[1:]):
        for key in mapping:
            np.testing.assert_array_equal(mapping[key], original[key])
    row["r_x"][0] = 99
    assert data[3]["r_x"] == .5


@pytest.mark.parametrize("bad_alpha", [True, -1., np.inf, 1j, [0.]])
def test_invalid_amplitude(shared_budget, bad_alpha):
    with pytest.raises(ValueError):
        _solve(_minkowski(5), shared_budget, bad_alpha)


@pytest.mark.parametrize("kind", ["too_few", "too_many", "nonuniform", "wrong_source", "reversed",
                                  "state_shape", "jet_missing", "rate_shape", "boolean", "nonfinite",
                                  "tube", "tube_lapse", "zero_radius", "order", "budget"])
def test_malformed_inputs_stop_without_charging(shared_budget, kind):
    data = list(_minkowski(5))
    kwargs = dict(alpha=.2, quadrature_order=4, tube=copy.deepcopy(TUBE), budget=shared_budget)
    if kind == "too_few":
        data = _minkowski(2)
    elif kind == "too_many":
        data = _minkowski(514)
    elif kind == "nonuniform":
        data[0][1] += .001
    elif kind == "wrong_source":
        data[0] += .1
    elif kind == "reversed":
        data[0] = data[0][::-1]
    elif kind == "state_shape":
        data[1]["r"] = .6
    elif kind == "jet_missing":
        del data[2]["log_C_Ubar"]
    elif kind == "rate_shape":
        data[3]["r_x"] = [.5]
    elif kind == "boolean":
        data[1]["L"] = [True, 0., 0., 0., 0.]
    elif kind == "nonfinite":
        data[2]["c_Ubar"][2] = np.nan
    elif kind == "tube":
        kwargs["tube"]["r"] = (.7, .8)
    elif kind == "tube_lapse":
        kwargs["tube"]["L"] = (-1000., 1000.)
    elif kind == "zero_radius":
        kwargs["tube"]["r"] = (0., 1.)
    elif kind == "order":
        kwargs["quadrature_order"] = True
    elif kind == "budget":
        kwargs["budget"] = None
    before = shared_budget.evaluations
    with pytest.raises(ValueError):
        exterior_characteristic_slice(*data, **kwargs)
    assert shared_budget.evaluations == before


def test_budget_exhaustion_keeps_caller_clock_and_counts():
    budget = Budget(maximum_evaluations=7, cpu_seconds=30., evaluations=2)
    started = budget.started
    with pytest.raises(RuntimeError, match="budget exhausted"):
        _solve(_minkowski(5), budget)
    assert budget.started == started
    assert budget.evaluations == 12  # input pass succeeds; first prefix charge stops


def test_finite_input_that_overflows_derived_jets_stops(shared_budget):
    data = list(_minkowski(5))
    data[2]["c_Ubar"][:] = 1e308
    with pytest.raises(ValueError, match="nonfinite"):
        _solve(data, shared_budget)
