"""Independent local checks for formal cap mass and material-wall response."""
import time

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_atlas_cap_wall_response as model


def test_exact_stress_frame_and_junction_certificate():
    started = time.process_time()
    assert model.exact_algebra_checks(started) == model.expected_certificate()
    assert time.process_time()-started < model.CPU_SECONDS


def test_independent_dimensionful_cap_hilbert_transport():
    R, H, field, tang, norm, speed, beta, acceleration = sp.symbols(
        "R H field tang norm speed beta acceleration")
    Htt, Htn = sp.symbols("Htt Htn")
    Hnn = Htt+2*(speed*tang-beta*norm)/R+2*H*H*field
    box_phi_squared = 2*(-tang*tang+norm*norm)+4*H*H*field*field
    # Start from the canonical Hilbert tensor rather than the module formulas.
    Ttn = tang*norm-(2*tang*norm+2*field*Htn)/6
    Tnn = (tang*tang+norm*norm)/2+(
        -3*H*H*field*field+box_phi_squared-2*norm*norm-2*field*Hnn)/6
    jordan = R*R*(beta*Ttn-speed*Tnn)/2
    endpoint = R*R*field*(speed*tang-beta*norm)/6-H*H*R**3*field*field/12
    rates = {R: speed, field: tang, speed: acceleration*beta+H*H*R,
             beta: acceleration*speed, tang: Htt+acceleration*norm,
             norm: Htn+acceleration*tang}
    endpoint_rate = sum(sp.diff(endpoint, variable)*rate for variable, rate in rates.items())
    minimal = R*R*(beta*tang*norm-speed*(tang*tang+norm*norm)/2)/2
    potential = H*H*R*R*speed*field*field/2
    assert sp.factor(minimal+potential+endpoint_rate-jordan) == 0
    assert sp.factor(minimal+endpoint_rate-jordan) != 0


def test_first_derivative_ingredients_match_independent_trace_forms():
    x, w, j = np.array([.45, .5, .55]), np.array([.1, -.2, .3]), np.array([1., -2., 3.])
    fields = model.field_mass_ingredients(x, w, j)
    g = model.reference_coefficients(x)
    assert np.allclose(fields["B_plus"]-fields["B_minus"], 2*fields["p"], rtol=0, atol=1e-14)
    for side, f in (("minus", 1-x*x), ("plus", 1-.8/x)):
        beta, v, p = g["beta_"+side], g["v"], w/x
        # Expand both squared null derivatives independently of the adapter.
        expected = (beta-v)*j*j/2-f*j*p/2+v*f*p*p/4
        if side == "minus":
            expected += x*x*v*p*p/2
        assert np.allclose(fields["L_"+side], expected, rtol=2e-14, atol=1e-14)
        assert np.allclose(fields["Q_"+side], fields["Q_"+side+"_trace"], rtol=1e-14, atol=1e-14)


def test_quiet_data_and_field_amplitude_scaling():
    for value in model.field_mass_ingredients(.5, 0., 0.).values():
        assert value == 0
    quiet = model.wall_response_coefficients(.5, 0., 0., 0.)
    assert all(value == 0 for value in quiet.values())
    base = model.field_mass_ingredients(.5, .7, -.2)
    scaled = model.field_mass_ingredients(.5, 3*.7, 3*(-.2))
    for key in ("L_minus", "L_plus", "Q_minus", "Q_plus", "cap_potential_density"):
        assert scaled[key] == pytest.approx(9*base[key], rel=2e-14)


def test_actual_cap_radius_and_linear_constraints_are_retained():
    x, rho, mu, nu = .51, .3, -.2, .4
    state = model.reference_coefficients(x)
    result = model.wall_response_coefficients(x, rho, mu, nu)
    assert result["cap_total_mass"] == pytest.approx(nu+1.5*x*x*rho)
    assert model.cap_total_mass_coefficient(x, rho, nu) == result["cap_total_mass"]
    assert result["displacement_rate"] == pytest.approx(result["cap_side_displacement_rate"])
    assert result["delta_beta_minus"]-result["delta_beta_plus"] == pytest.approx(2*rho)
    v, bp, bm = state["v"], state["beta_plus"], state["beta_minus"]
    assert bp*result["delta_beta_plus"]-v*result["displacement_rate"]-.4*rho/x**2+mu/x == pytest.approx(0, abs=2e-14)
    assert bm*result["delta_beta_minus"]-v*result["displacement_rate"]+x*rho+nu/x == pytest.approx(0, abs=2e-14)
    mean_beta = (bp+bm)/2
    mean_delta_beta = (result["delta_beta_plus"]+result["delta_beta_minus"])/2
    assert result["delta_normal_acceleration"]+2*mean_delta_beta/x-2*mean_beta*rho/x**2 == pytest.approx(0, abs=2e-14)


def test_beta_plus_zero_is_regular_without_changing_wall_branch():
    x = (.16)**(1/3)
    state = model.reference_coefficients(x)
    assert abs(state["beta_plus"]) < 1e-14
    response = model.wall_response_coefficients(x, .2, -.3, .4)
    ingredients = model.field_mass_ingredients(x, .3, -.7)
    assert all(np.isfinite(value) for value in (*response.values(), *ingredients.values()))
    assert response["displacement_forcing"] == pytest.approx(-.3/(x*state["v"]), abs=1e-14)


@pytest.mark.parametrize("invalid", [True, False, np.nan, np.inf, -np.inf, 1j, "0.5", [], [True]])
def test_malformed_inputs_are_rejected(invalid):
    with pytest.raises(ValueError):
        model.field_mass_ingredients(.5, invalid, 0.)
    with pytest.raises(ValueError):
        model.wall_response_coefficients(.5, 0., invalid, 0.)
    with pytest.raises(ValueError):
        model.reference_coefficients(invalid)


@pytest.mark.parametrize("outside", [.449, .551, 0., -1.])
def test_radius_domain_exit_is_rejected(outside):
    with pytest.raises(ValueError, match="wall strip"):
        model.reference_coefficients(outside)


def test_incompatible_broadcast_shapes_fail():
    with pytest.raises(ValueError, match="broadcast"):
        model.field_mass_ingredients([.45, .55], [1., 2., 3.], 0.)


def test_changed_cap_vacuum_formula_fails_independent_check(monkeypatch):
    monkeypatch.setitem(model.FORMULAS, "cap_endpoint", model.FORMULAS["cap_endpoint"]+"+h**2*x**3*p**2")
    with pytest.raises(ValueError, match="exact identity failed"):
        model.exact_algebra_checks()


def test_wrong_symbolic_runtime_and_exhausted_budget_fail(monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(sp, "__version__", "unexpected")
        with pytest.raises(RuntimeError, match="pinned SymPy"):
            model.exact_algebra_checks()
    monkeypatch.setattr(model.time, "process_time", lambda: 61.)
    with pytest.raises(RuntimeError, match="CPU budget"):
        model.exact_algebra_checks(started=0.)
