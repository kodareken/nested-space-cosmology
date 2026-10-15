"""Independent cap identities and manufactured/quiet metric controls."""
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_atlas_cap_response as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget, background
from recursive_horizons.nsc_atlas_cap_wall_response import (
    field_mass_ingredients, reference_coefficients, wall_response_coefficients,
)
from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _geometry


class ManufacturedScalar:
    def __init__(self, constant=0., time_slope=0.):
        self.constant, self.time_slope = constant, time_slope
        self.budget = PreparedFieldBudget()
        self.edges = np.array([-2., -1., 0., 1.])
        self.prepared_field = object()

    def jet(self, T, r):
        T, r = np.broadcast_arrays(T, r)
        p = self.constant+self.time_slope*T
        return {"p": p, "p_T": np.full_like(T, self.time_slope),
                "p_r": np.zeros_like(T), "p_center": p}


def test_radial_mass_flux_conservation_is_the_scalar_wave_equation():
    T, r = sp.symbols("T r", real=True)
    p = sp.Function("p")(T, r)
    A = 1-r*r
    pT, pr = sp.diff(p, T), sp.diff(p, r)
    density = r*r*(pT*pT/A+A*pr*pr)/4+r*r*p*p/2
    flux = r*r*A*pT*pr/2
    wave = sp.diff(p, T, 2)/A-sp.diff(r*r*A*pr, r)/(r*r)+2*p
    assert sp.simplify(sp.diff(density, T)-sp.diff(flux, r)-r*r*pT*wave/2) == 0
    wrong = density-r*r*p*p/2
    assert sp.simplify(sp.diff(wrong, T)-sp.diff(flux, r)-r*r*pT*wave/2) == -r*r*p*pT


def test_optical_combination_cancels_explicit_kinetic_terms():
    r = sp.symbols("r", positive=True)
    p, pT, pr = sp.symbols("p pT pr", real=True)
    mass, lapse = sp.Function("mass")(r), sp.Function("lapse")(r)
    A = 1-r*r
    density = r*r*(pT*pT/A+A*pr*pr)/4+r*r*p*p/2
    lapse_rate = r*(pr*pr+pT*pT/A**2)/2
    D = lapse-2*mass/(r*A)
    derivative = sp.diff(D, r).subs({sp.diff(mass, r): density,
                                     sp.diff(lapse, r): lapse_rate})
    expected = 2*mass*(1-3*r*r)/(r*r*A*A)-r*p*p/A
    assert sp.simplify(derivative-expected) == 0
    wrong = sp.diff(D, r).subs({sp.diff(mass, r): density-r*r*p*p/2,
                               sp.diff(lapse, r): lapse_rate})
    assert sp.simplify(wrong-expected-r*p*p/A) == 0


def test_constant_manufactured_scalar_retains_potential_and_central_gauge():
    field = ManufacturedScalar(constant=.7)
    radii = np.array([0., .1, .47])
    metric = model.radial_metric(field, .2, radii, radial_order=8)
    assert metric["nu_E"] == pytest.approx(.7**2*radii**3/6, abs=2e-16)
    assert metric["nu_J"] == pytest.approx(.7**2*radii**3/12, abs=2e-16)
    assert metric["delta_E"] == pytest.approx(np.full(3, -.7**2/12))
    assert np.count_nonzero(metric["mass_kinetic"]) == 0
    assert metric["mass_potential"][-1] > 0
    assert np.count_nonzero(metric["central_clock_normalization_defect"]) == 0
    # Omitting the transformed cap potential would reverse this Jordan mass.
    wrong = metric["nu_J"]-metric["mass_potential"]
    assert wrong[-1] < 0 < metric["nu_J"][-1]


def test_time_dependent_manufactured_scalar_matches_closed_radial_integrals():
    c, T, r = .4, .3, .5
    metric = model.radial_metric(ManufacturedScalar(time_slope=c), T, r, radial_order=16)
    expected_mass = c*c*(np.arctanh(r)-r)/4+c*c*T*T*r**3/6
    expected_lapse = -c*c*T*T/12+c*c*(1/(1-r*r)-1)/4
    assert float(metric["nu_E"]) == pytest.approx(expected_mass, abs=2e-15)
    assert float(metric["delta_E"]) == pytest.approx(expected_lapse, abs=2e-15)


def test_wall_normalization_detects_independent_cap_mass_mismatch():
    x, rho, mu, nu = .51, .13, -.2, .05
    g = reference_coefficients(x)
    scalar = field_mass_ingredients(x, .2, .3)
    response = wall_response_coefficients(x, rho, mu, nu)
    wall = {**{key: float(value) for key, value in g.items()},
            **{key: float(value) for key, value in scalar.items()},
            **{key: float(value) for key, value in response.items()},
            "rho": rho, "rho_s": float(response["displacement_rate"])}
    nuE = nu-float(scalar["Q_minus"])
    valid = model.wall_clock_coefficients(wall, nuE, .7)
    changed = model.wall_clock_coefficients(wall, nuE+.01, .7)
    assert valid["wall_normalization_defect"] == pytest.approx(0., abs=2e-15)
    assert changed["wall_normalization_defect"] == pytest.approx(.02/(x*(1-x*x)), abs=2e-15)


def test_ray_quadrature_matches_a_closed_manufactured_metric():
    field = ManufacturedScalar(constant=.7)
    x = .5
    ray = model.ray_integral(field, .4, x, radial_order=8, ray_order=16)
    expected = .7**2*(x/(6*(1-x*x))-np.arctanh(x)/12)
    assert ray["integral"] == pytest.approx(expected, abs=2e-15)
    assert ray["positive_integral"]+ray["negative_integral"] == ray["integral"]
    assert ray["negative_integral"] == 0
    assert 0 < ray["integrand_sample_min"] <= ray["integrand_sample_max"]
    assert model.null_integrand(0., 0., -.7**2/12) == pytest.approx(.7**2/12)
    with pytest.raises(ValueError, match="zero enclosed"):
        model.null_integrand(0., 1., 0.)


@pytest.fixture(scope="module")
def geometry():
    return source_geometry()


class QuietCap(ManufacturedScalar):
    def __init__(self, geometry):
        super().__init__()
        self.geometry = _geometry(129)
        qa, qb = self.geometry["arrival"]([geometry["Va"], geometry["Vb"]])
        Ta, Tb = qa-np.arctanh(.45), qb-np.arctanh(.55)
        self.domain = {key: geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
        self.domain.update(Ta=float(Ta), Tb=float(Tb), qa=float(qa), qb=float(qb),
                           qmin=float(Ta-np.arctanh(.45)), xa=.45, xb=.55)
        self.edges = np.array([qa, qa+(qb-qa)/4, qa+3*(qb-qa)/4, qb])

    def jet(self, T, r):
        T, r = np.broadcast_arrays(T, r)
        if np.any(T < self.domain["Ta"]-2e-13) or np.any(T+np.arctanh(r) > self.domain["qb"]+2e-13):
            raise ValueError("quiet test query exits the cap domain")
        return super().jet(T, r)


class QuietWall:
    def __init__(self, cap, geometry, label_shift=0.):
        self.geometry, self.cap, self.label_shift = geometry, cap, label_shift
        self.budget = BulkBudget(started=cap.budget.started, external=cap.budget)
        self.wedge = SimpleNamespace(prepared_field=cap.prepared_field)

    def _wall(self, t):
        V = self.geometry["Va"]+t*(self.geometry["Vb"]-self.geometry["Va"])
        U = self.cap.geometry["wall"](V)
        x = min(.55, max(.45, float(background(U, V)["r"])))
        g = {key: float(value) for key, value in reference_coefficients(x).items()}
        Vs = -V/(1.6*(g["v"]-g["beta_plus"]))
        return {**g, "p": 0., "A": 0., "B_minus": 0., "rho": 0., "rho_s": 0.,
                "delta_beta_minus": 0., "integral_L_minus": 0., "nu_minus": 0.,
                "Q_minus": 0., "V": V, "Vdot": Vs,
                "ds_dt": (self.geometry["Vb"]-self.geometry["Va"])/Vs}

    def sample(self, t):
        return {"V1": self.label_shift*t}


def test_quiet_cap_recovers_zero_metric_wall_clock_and_both_arrivals(geometry):
    cap = QuietCap(geometry)
    solution = model.solve_cap_response(cap, QuietWall(cap, geometry),
                                        radial_order=4, ray_order=4, time_panels=4)
    summary = solution.summary()
    assert summary["wall_metric_evaluations"] == 25
    assert summary["ray_evaluations"] == 9
    for row in summary["cap_samples"]:
        for key in ("nu_E_wall", "nu_J_wall", "delta_E_wall", "T1", "T1_s",
                    "Delta_C_fixed_s", "Delta_C_fixed_V", "ray_positive_integral",
                    "ray_negative_integral", "ray_integrand_sample_min", "ray_integrand_sample_max"):
            assert row[key] == 0
    for key, value in summary["maximum_absolute_defects"].items():
        if key != "reference_arrival_rate_defect":
            assert value == 0
    assert summary["maximum_absolute_defects"]["reference_arrival_rate_defect"] < 2e-8


def test_fixed_exterior_label_conversion_has_the_inverse_embedding_sign(geometry):
    cap = QuietCap(geometry)
    solution = model.solve_cap_response(cap, QuietWall(cap, geometry, label_shift=.3),
                                        radial_order=4, ray_order=4, time_panels=4)
    row = solution.sample(.5)
    assert row["Delta_C_fixed_s"] == 0
    assert row["Delta_C_fixed_V"] == pytest.approx(-row["q0_s"]*.15/row["V0_s"])
    assert row["Delta_C_fixed_V"] < 0
    assert row["fixed_V_conversion_defect"] == 0


def test_changed_prepared_state_and_unshared_budget_are_rejected(geometry):
    cap = QuietCap(geometry)
    wall = QuietWall(cap, geometry)
    wall.wedge.prepared_field = object()
    with pytest.raises(ValueError, match="same prepared field solution"):
        model.solve_cap_response(cap, wall)
    wall.wedge.prepared_field = cap.prepared_field
    with pytest.raises(ValueError, match="shared"):
        model.solve_cap_response(cap, wall, budget=BulkBudget())
    with pytest.raises(ValueError, match="cumulative"):
        model.solve_cap_response(cap, wall,
            budget=BulkBudget(started=cap.budget.started, external=cap.budget))
    with pytest.raises(ValueError, match="radial order"):
        model.solve_cap_response(cap, wall, radial_order=16.)


def test_domain_exit_and_nonfinite_queries_fail_before_quadrature(geometry):
    cap = QuietCap(geometry)
    with pytest.raises(ValueError, match="domain"):
        model.radial_metric(cap, cap.domain["qb"]+.1, .2)
    with pytest.raises(ValueError, match="finite"):
        model.radial_metric(cap, np.nan, .2)
    with pytest.raises(ValueError, match="admitted wall"):
        model.ray_integral(cap, cap.domain["qa"], .8)
