"""Independent algebra, manufactured jets and quiet-source wall controls."""
import time

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad

from recursive_horizons import nsc_atlas_wall_evolution as model
from recursive_horizons.nsc_atlas_cap_wall_response import (
    reference_coefficients, field_mass_ingredients, wall_response_coefficients,
)
from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _geometry


@pytest.fixture(scope="module")
def geometry():
    return source_geometry()


class ManufacturedWedge:
    def __init__(self, geometry, gradient=0.):
        self.domain = dict(geometry)
        self.geometry = _geometry(129)
        self.grid = np.linspace(geometry["U0"], geometry["Ub"], 17)
        self.budget = PreparedFieldBudget()
        self.gradient = gradient

    def bulk(self, U, t, order=32):
        V = self.domain["Va"]+t*(self.domain["Vb"]-self.domain["Va"])
        bg = model.background(U, V)
        p = self.gradient*(U-self.domain["U0"])
        pU, pV = np.full_like(U, self.gradient), np.zeros_like(U)
        return {"p": p, "p_U": pU, "p_V": pV,
                "chi": bg["r"]*p, "chi_U": bg["r_U"]*p+bg["r"]*pU,
                "chi_V": bg["r_V"]*p,
                "chi_U_lower_limit_bound": np.zeros_like(U)}


def test_clock_ratio_follows_from_normalization_including_shift():
    alpha, V, V1, a, d, dd, zeta = sp.symbols("alpha V V1 a d dd zeta", nonzero=True)
    # Z0 depends on the moved V, even at fixed proper wall time.
    exact = -(V+alpha*V1)/(2*a*(d+alpha*dd)*(1+alpha*zeta))
    V0_s = -V/(2*a*d)
    V1_s = sp.diff(exact, alpha).subs(alpha, 0)
    eta_s = sp.factor(V1_s/V-V1*V0_s/V**2)
    assert sp.factor(eta_s-(zeta+dd/d)/(2*a*d)) == 0
    incorrectly_frozen_V_s = -V*(1-alpha*zeta)/(2*a*(d+alpha*dd))
    wrong = sp.diff(incorrectly_frozen_V_s, alpha).subs(alpha, 0)/V-V1*V0_s/V**2
    assert sp.factor(wrong-eta_s) != 0


def test_wall_frame_lapse_requires_nonzero_physical_p_U():
    r, ru, p, zeta = sp.symbols("r ru p zeta", nonzero=True)
    pU = -p*ru/r  # chi_U=(r*p)_U=0, not p_U=0.
    generic = zeta+p*p/12-r*p*pU/(6*ru)
    assert sp.factor(generic-zeta-p*p/4) == 0
    omitted = zeta+p*p/12
    assert sp.factor(omitted-generic) != 0


def test_wall_frame_mass_matches_independent_first_jet_trace():
    x, V, chi, chiV = .51, -.7, .6, -.9
    from recursive_horizons.nsc_atlas_field_preparation import radius_product
    U = float(radius_product(x))/V
    bg, g = model.background(U, V), reference_coefficients(x)
    v, bp = float(g["v"]), float(g["beta_plus"])
    Vdot = -V/(2*model.A_RADIUS*(v-bp))
    p, pU = chi/x, -chi*bg["r_U"]/x**2
    pV = (chiV-p*bg["r_V"])/x
    generic = model.M*p*p/12+x*x*p*(bg["r_U"]*pV+bg["r_V"]*pU)/(3*bg["C"])
    trace = field_mass_ingredients(x, chi, chiV*Vdot)["Q_plus"]
    assert float(generic) == pytest.approx(float(trace), abs=2e-15)
    assert pU != 0


def test_radial_constraints_match_independent_integrals_and_generic_frame(geometry):
    field = ManufacturedWedge(geometry, gradient=.3)
    t = .6
    row = model.radial_response(field, geometry, t, upper_U=geometry["Ua"],
                                radial_order=4, source_scale=0.)
    V = row["V"]
    def mass_integrand(U):
        g = model.background(U, V)
        return -float(g["r"]**2*g["r_V"]*.3**2/g["C"])
    def lapse_integrand(U):
        g = model.background(U, V)
        return float(g["r"]*.3**2/(2*g["r_U"]))
    assert row["mu_E"] == pytest.approx(quad(mass_integrand, geometry["U0"], geometry["Ua"])[0], abs=2e-14)
    assert row["zeta_E"] == pytest.approx(quad(lapse_integrand, geometry["U0"], geometry["Ua"])[0], abs=2e-14)
    assert row["mu_E"] < 0 < row["zeta_E"]
    assert row["chi_U"] != 0
    assert abs(row["zeta_J"]-row["zeta_E"]-row["p"]**2/4) > 1e-8


def test_gauss_partial_integrals_are_exact_for_cubic_rates():
    nodes, basis = model._gauss_antiderivatives()
    for degree in range(4):
        for upper in (-1., -.3, .6, 1.):
            computed = sum(poly(upper)*node**degree for node, poly in zip(nodes, basis))
            expected = (upper**(degree+1)-(-1.)**(degree+1))/(degree+1)
            assert computed == pytest.approx(expected, abs=2e-15)


def test_clock_is_regular_at_signed_beta_zero_and_respects_null_boost():
    x = .16**(1/3)
    g = reference_coefficients(x)
    assert abs(g["beta_plus"]) < 1e-14
    rho, mu, nu, zeta = .2, -.4, -.3, .12
    rate = model.clock_rate(x, rho, mu, nu, zeta)
    assert np.isfinite(rate)
    response = wall_response_coefficients(x, rho, mu, nu)
    d = float(g["v"]-g["beta_plus"])
    E = zeta+float(response["displacement_rate"]-response["delta_beta_plus"])/d
    for V in (-.8, -.8/3):
        eta = .07
        Vs = -V/(2*model.A_RADIUS*d)
        V1s = Vs*eta+V*rate
        assert V1s == pytest.approx(Vs*(eta-E))


def test_quiet_source_evolution_is_identically_zero_with_real_reference(geometry):
    solution = model.solve_wall_evolution(ManufacturedWedge(geometry), geometry,
                                         source_scale=0., time_panels=4)
    result = solution.summary()
    assert result["radial_evaluations"] == 25
    assert result["clock_time_order"] == 4
    for row in result["wall_samples"]:
        for key in ("mu_plus", "nu_minus", "rho", "rho_s", "eta", "eta_s", "V1", "V1_s"):
            assert row[key] == 0
    for key, value in result["maximum_absolute_defects"].items():
        if key != "reference_clock_normalization_defect":
            assert value == 0
    # The 129-knot reference spline has its own finite interpolation error.
    assert result["maximum_absolute_defects"]["reference_clock_normalization_defect"] < 2e-8


def test_actual_source_cannot_be_combined_with_quiet_field(geometry):
    with pytest.raises(ValueError, match="incompatible"):
        model.solve_wall_evolution(ManufacturedWedge(geometry), geometry)


def test_domain_orders_and_shared_budget_fail_closed(geometry):
    field = ManufacturedWedge(geometry)
    with pytest.raises(ValueError, match="retained"):
        model.radial_response(field, geometry, .5, upper_U=0., source_scale=0.)
    with pytest.raises(ValueError, match="radial order"):
        model.solve_wall_evolution(field, geometry, radial_order=3, source_scale=0.)
    with pytest.raises(ValueError, match="time panels"):
        model.solve_wall_evolution(field, geometry, time_panels=8., source_scale=0.)
    with pytest.raises(ValueError, match="shared"):
        model.solve_wall_evolution(field, geometry, budget=model.BulkBudget(), source_scale=0.)
    field.budget.started = time.process_time()-121
    with pytest.raises(RuntimeError, match="budget"):
        model.solve_wall_evolution(field, geometry, source_scale=0.)


@pytest.mark.parametrize("invalid", [True, np.nan, -1, 1.1])
def test_time_queries_fail_closed(geometry, invalid):
    with pytest.raises(ValueError):
        model.radial_response(ManufacturedWedge(geometry), geometry, invalid, source_scale=0.)
