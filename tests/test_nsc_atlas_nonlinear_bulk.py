"""Independent canonical conversion, vacuum and constraint controls."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons import nsc_atlas_nonlinear_bulk as model


@pytest.fixture(autouse=True)
def bounded_test_cpu():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def inputs():
    return dict(r=.5, C=1.3, c=.2), dict(r_U=-.4, r_V=.6, c_U=.3, c_V=-.2)


@pytest.mark.parametrize("region,h_squared", [("cap", 1), ("exterior", 0)])
def test_canonical_action_and_independent_mixed_chain_rule(region, h_squared):
    alpha = np.array([.03, .2, .7])
    r = np.array([[.4], [.8]])
    C, rU, rV = 1.2, -.3, .7
    psi = np.array([.12, -.2, .31])
    psiU, psiV = np.array([[.23], [-.17]]), np.array([-.1, .21, .07])
    w = np.sqrt(6)*np.sinh(psi/np.sqrt(6))
    w_psi = np.cosh(psi/np.sqrt(6))
    state = dict(r=r, C=C, c=r*w/np.sqrt(alpha))
    jets = dict(r_U=rU, r_V=rV, c_U=(rU*w+r*w_psi*psiU)/np.sqrt(alpha),
                c_V=(rV*w+r*w_psi*psiV)/np.sqrt(alpha))
    result = model.bulk_operator(state, jets, alpha=alpha, region=region)
    # Differentiate the canonical action's potential, before reducing to c.
    P = 3*h_squared*np.cosh(psi/np.sqrt(6))**4
    P_psi = 2*np.sqrt(6)*h_squared*np.cosh(psi/np.sqrt(6))**3*np.sinh(psi/np.sqrt(6))
    rUV = -rU*rV/r-C/(4*r)+C*r*P/4
    LUV = 2*rU*rV/r**2+C/(2*r**2)-psiU*psiV
    psiUV = -(rU*psiV+rV*psiU)/r-C*P_psi/4
    cUV = (rUV*w+rU*w_psi*psiV+rV*w_psi*psiU
           +r*(w_psi*psiUV+w*psiU*psiV/6))/np.sqrt(alpha)
    for key, expected in (("r_UV", rUV), ("log_C_UV", LUV), ("c_UV", cUV)):
        np.testing.assert_allclose(result[key], expected, rtol=3e-14, atol=3e-14)
    np.testing.assert_allclose(np.sqrt(alpha)*result["psi_U_normalized"], psiU+np.zeros((2, 3)), atol=1e-15)
    np.testing.assert_allclose(np.sqrt(alpha)*result["psi_V_normalized"], psiV+np.zeros((2, 3)), atol=1e-15)
    np.testing.assert_allclose(result["potential"], P+np.zeros((2, 3)), rtol=2e-15)


@pytest.mark.parametrize("geometry", ["de_sitter", "schwarzschild", "minkowski"])
def test_exact_static_vacuum_geometries_and_nonzero_normalized_scalar(geometry):
    if geometry == "de_sitter":
        r, region = np.array([.2, .45, .55]), "cap"
        f, fp, fpp = 1-r*r, -2*r, np.full_like(r, -2.)
        mass = r**3/2
    elif geometry == "schwarzschild":
        r, region = np.array([.9, 1.2, 2.]), "exterior"
        f, fp, fpp = 1-.8/r, .8/r**2, -1.6/r**3
        mass = np.full_like(r, .4)
    else:
        r, region = np.array([.3, .8, 1.4]), "exterior"
        f, fp, fpp = np.ones_like(r), np.zeros_like(r), np.zeros_like(r)
        mass = np.zeros_like(r)
    rU, rV = -f/2, f/2
    # Schwarzschild/Minkowski constant scalar has reduced field c=r.
    # The cap's constant reduced field is a solution of its free null equation.
    c = np.ones_like(r) if geometry == "de_sitter" else r
    cU, cV = (np.zeros_like(r), np.zeros_like(r)) if geometry == "de_sitter" else (rU, rV)
    result = model.bulk_operator(dict(r=r, C=f, c=c),
        dict(r_U=rU, r_V=rV, c_U=cU, c_V=cV), alpha=0., region=region)
    np.testing.assert_allclose(result["r_UV"], -f*fp/4, atol=2e-15)
    np.testing.assert_allclose(result["log_C_UV"], -f*fpp/4, atol=2e-15)
    np.testing.assert_allclose(result["m"], mass, atol=2e-15)
    expected_K = f*.4/(2*r**3) if geometry == "schwarzschild" else np.zeros_like(r)
    np.testing.assert_allclose(result["K"], expected_K, atol=2e-15)
    np.testing.assert_allclose(result["c_UV"], -expected_K*c, atol=2e-15)
    massU = 1.5*r*r*rU if geometry == "de_sitter" else np.zeros_like(r)
    massV = 1.5*r*r*rV if geometry == "de_sitter" else np.zeros_like(r)
    np.testing.assert_allclose(result["mass_U_on_constraints"], massU, atol=2e-15)
    np.testing.assert_allclose(result["mass_V_on_constraints"], massV, atol=2e-15)


def test_trapped_schwarzschild_in_the_original_affine_source_radius_domain():
    r, C = np.array([.4, .45, .5, .55]), 1.7
    # Affine source radius x=r has r_V=1 and C=C0>0. Both null radius
    # derivatives are positive inside r<2M; no exterior static chart is used.
    rU, rV = C*(.8/r-1)/4, np.ones_like(r)
    assert np.all(rU > 0) and np.all(rV > 0)
    row = model.bulk_operator(dict(r=r, C=C, c=r),
        dict(r_U=rU, r_V=rV, c_U=rU, c_V=rV), alpha=0., region="exterior")
    for key, expected in (("r_UV", -.2*C/r**2), ("log_C_UV", .4*C/r**3),
                           ("K", .2*C/r**3), ("m", np.full_like(r, .4)),
                           ("c_UV", -.2*C/r**2)):
        np.testing.assert_allclose(row[key], expected, rtol=2e-15, atol=2e-15)
    np.testing.assert_array_equal(row["mass_U_on_constraints"], np.zeros_like(r))
    np.testing.assert_array_equal(row["mass_V_on_constraints"], np.zeros_like(r))


@pytest.mark.parametrize("region", ["cap", "exterior"])
def test_mass_derivative_retains_nonzero_null_constraint_residuals(region):
    state, jets = inputs()
    alpha = .27
    row = model.bulk_operator(state, jets, alpha=alpha, region=region)
    r, C, c = (state[key] for key in ("r", "C", "c"))
    rU, rV, cU, cV = (jets[key] for key in ("r_U", "r_V", "c_U", "c_V"))
    psiU = np.sqrt(alpha)*(cU/r-c*rU/r**2)/np.sqrt(1+alpha*c*c/(6*r*r))
    psiV = np.sqrt(alpha)*(cV/r-c*rV/r**2)/np.sqrt(1+alpha*c*c/(6*r*r))
    LU, LV, EU, EV = .31, -.21, .12, -.08
    rUU, rVV = LU*rU-r*psiU*psiU/2+EU, LV*rV-r*psiV*psiV/2+EV
    assert row["constraint_U_rhs"] == pytest.approx(rUU-LU*rU-EU)
    assert row["constraint_V_rhs"] == pytest.approx(rVV-LV*rV-EV)
    # Directly differentiate m=r/2+2*r*r_U*r_V/C with arbitrary nonzero E.
    mU = rU/2+2*(rU*rU*rV+r*rUU*rV+r*rU*row["r_UV"]-r*rU*rV*LU)/C
    mV = rV/2+2*(rV*rU*rV+r*row["r_UV"]*rV+r*rU*rVV-r*rU*rV*LV)/C
    assert mU == pytest.approx(row["mass_U_on_constraints"]+row["mass_U_residual_coefficient"]*EU)
    assert mV == pytest.approx(row["mass_V_on_constraints"]+row["mass_V_residual_coefficient"]*EV)
    assert abs(mU-row["mass_U_on_constraints"]) > .01
    assert abs(mV-row["mass_V_on_constraints"]) > .01


def test_null_coordinate_rescaling_and_transport_weights():
    state, jets = inputs()
    original = model.bulk_operator(state, jets, alpha=.13, region="cap")
    a, b = 2.3, .7
    state["C"] /= a*b
    jets = {key: value/(a if key.endswith("_U") else b) for key, value in jets.items()}
    changed = model.bulk_operator(state, jets, alpha=.13, region="cap")
    factors = {"r_UV": 1/(a*b), "log_C_UV": 1/(a*b), "c_UV": 1/(a*b), "K": 1/(a*b),
               "psi_U_normalized": 1/a, "psi_V_normalized": 1/b,
               "constraint_U_rhs": 1/a**2, "constraint_V_rhs": 1/b**2,
               "constraint_U_transport_coefficient": 1/b, "constraint_V_transport_coefficient": 1/a,
               "mass_U_on_constraints": 1/a, "mass_V_on_constraints": 1/b,
               "mass_U_residual_coefficient": a, "mass_V_residual_coefficient": b}
    for key, value in original.items():
        assert changed[key] == pytest.approx(factors.get(key, 1)*value, rel=2e-14, abs=2e-14)
    # These coefficients preserve rE along the respective null direction.
    assert original["constraint_U_transport_coefficient"]*state["r"]+inputs()[1]["r_V"] == pytest.approx(0.)
    assert original["constraint_V_transport_coefficient"]*state["r"]+inputs()[1]["r_U"] == pytest.approx(0.)


def test_small_alpha_increment_is_not_lost_to_vacuum_subtraction():
    # Exact binary cap data make the vacuum kernel cancel exactly.
    state = dict(r=.5, C=.75, c=.5)
    jets = dict(r_U=-.375, r_V=.375, c_U=-.375, c_V=.375)
    zero = model.bulk_operator(state, jets, alpha=0., region="cap")
    tiny = model.bulk_operator(state, jets, alpha=1e-18, region="cap")
    assert zero["K"] == 0 and zero["potential_increment"] == 0
    assert tiny["D"] == 1 and tiny["potential"] == 3
    assert tiny["potential_increment"] == pytest.approx(1e-18, rel=1e-15, abs=0.)
    assert tiny["K"] == pytest.approx(-.75e-18/12, rel=1e-15, abs=0.)
    assert tiny["c_UV"] != 0


def test_alpha_zero_keeps_normalized_field_and_zero_null_radius_derivative():
    state, jets = inputs()
    jets["r_U"] = 0.
    zero = model.bulk_operator(state, jets, alpha=0., region="exterior")
    small = model.bulk_operator(state, jets, alpha=1e-12, region="exterior")
    assert zero["psi_U_normalized"] != 0 and zero["c_UV"] != 0
    assert zero["constraint_U_rhs"] == zero["constraint_V_rhs"] == 0.
    for key in zero:
        assert small[key] == pytest.approx(zero[key], rel=1e-11, abs=1e-11)
    quiet_state, quiet_jets = copy.deepcopy(state), copy.deepcopy(jets)
    quiet_state["c"] = quiet_jets["c_U"] = quiet_jets["c_V"] = 0.
    quiet = model.bulk_operator(quiet_state, quiet_jets, alpha=.7, region="exterior")
    for key in ("c_UV", "psi_U_normalized", "psi_V_normalized", "constraint_U_rhs", "constraint_V_rhs"):
        assert quiet[key] == 0.


@pytest.mark.parametrize("invalid", [True, False, np.nan, np.inf, -np.inf, 1j, "0.5", [], [True], [True, 1.], None])
def test_malformed_numeric_inputs_are_rejected(invalid):
    state, jets = inputs()
    jets["c_U"] = invalid
    with pytest.raises(ValueError, match="finite real"):
        model.bulk_operator(state, jets, alpha=.2, region="cap")


@pytest.mark.parametrize("section,key,value", [("state", "r", 0.), ("state", "r", -.1),
    ("state", "C", 0.), ("state", "C", -1.), ("parameter", "alpha", -.1)])
def test_invalid_physical_domains(section, key, value):
    state, jets = inputs()
    alpha = .2
    if section == "state":
        state[key] = value
    else:
        alpha = value
    with pytest.raises(ValueError, match="require"):
        model.bulk_operator(state, jets, alpha=alpha, region="cap")


@pytest.mark.parametrize("region", [True, 0, None, "de_sitter", "CAP", ["cap"]])
def test_only_the_fixed_model_regions_are_accepted(region):
    with pytest.raises(ValueError, match="region"):
        model.bulk_operator(*inputs(), alpha=.2, region=region)


def test_mapping_shapes_and_no_input_mutation():
    state, jets = inputs()
    before = copy.deepcopy((state, jets))
    model.bulk_operator(state, jets, alpha=.2, region="cap")
    assert (state, jets) == before
    with pytest.raises(ValueError, match="state"):
        model.bulk_operator({**state, "m": .4}, jets, alpha=.2, region="cap")
    with pytest.raises(ValueError, match="first_jets"):
        model.bulk_operator(state, {"r_U": 1.}, alpha=.2, region="cap")
    with pytest.raises(ValueError, match="state"):
        model.bulk_operator(None, jets, alpha=.2, region="cap")
    state["r"], jets["c_U"] = [.4, .5], [.2, .3, .4]
    with pytest.raises(ValueError, match="broadcast"):
        model.bulk_operator(state, jets, alpha=.2, region="cap")


@pytest.mark.parametrize("section,key,value", [("state", "c", 1e308), ("state", "r", 1e-200),
    ("state", "C", 1e-320), ("jets", "c_U", 1e308)])
def test_nonrepresentable_derived_quantities_are_rejected(section, key, value):
    state, jets = inputs()
    (state if section == "state" else jets)[key] = value
    with pytest.raises(ValueError, match="overflow|finite"):
        model.bulk_operator(state, jets, alpha=1., region="cap")
