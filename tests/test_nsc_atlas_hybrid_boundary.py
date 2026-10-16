"""Independent polar boost, canonical action and exterior gauge controls."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary


@pytest.fixture(autouse=True)
def bounded_test_cpu():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def inputs():
    return (dict(r=.5, c=.3, eta_plus=.2),
            dict(r_V=.4, c_V=-.2, log_C_V=.17),
            dict(A=.75, delta=-.03, w_minus=.4),
            dict(C_wall=1.3, log_C_s=.09))


def test_signed_vacuum_reference_through_beta_plus_zero():
    state, outgoing, cap, gauge = inputs()
    x = np.array([.45, .5, .16**(1/3), .55])
    # Independent signed junction solution for m_plus=.4, h_minus=1, b=2.
    bp, bm = .2/x**2-1.25*x, .2/x**2+.75*x
    v = np.sqrt(bm**2-(1-x*x))
    Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["C_wall"])
    state.update(r=x, c=0.)
    outgoing.update(r_V=(v+bp)/(2*Vp), c_V=0.)
    cap.update(A=1-x*x, delta=0., w_minus=0.)
    row = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=0.)
    for key, expected in (("r_s", v), ("beta_plus", bp), ("beta_minus", bm),
                          ("k_plus", -.4/x**3-1.25), ("k_minus", -.4/x**3+.75),
                          ("T_s", bm/(1-x*x))):
        np.testing.assert_allclose(row[key], expected, rtol=2e-14, atol=2e-14)
    assert abs(row["beta_plus"][2]) < 2e-15
    for key in ("c_s", "n_c", "w_plus", "Pi", "Phi", "f_s", "psi_s_normalized"):
        np.testing.assert_array_equal(row[key], np.zeros(4))
    assert all(np.all(np.isfinite(value)) for value in row.values())


def test_finite_amplitude_recovers_independent_canonical_polar_boost():
    state, outgoing, cap, gauge = inputs()
    alpha = np.array([.01, .2, .8])
    r, f = np.array([[.42], [.57]]), np.array([.3, -.4, .1])
    A, boost = np.array([[.7], [.8]]), np.array([[-.3], [.6]])
    Pi, Phi = np.array([.2, -.1, .3]), np.array([-.4, .2, .1])
    q = np.sqrt(alpha/6)*f
    u = np.sqrt(6/alpha)*np.sinh(q)
    cosh = np.cosh(q)
    v, bm = np.sqrt(A)*np.sinh(boost), np.sqrt(A)*np.cosh(boost)
    tension = 2*cosh**3
    bp = bm-tension*r
    # Canonical field gradient in the independently chosen polar orthonormal
    # frame, boosted to the wall tangent and outward normal.
    f_s, Nminus = bm*Pi+v*Phi, v*Pi+bm*Phi
    cs, nc = v*u+r*cosh*f_s, bm*u+r*cosh*Nminus
    Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["C_wall"])
    state.update(r=r, c=r*u)
    outgoing.update(r_V=(v+bp)/(2*Vp), c_V=(cs+nc)/(2*Vp))
    cap.update(A=A, w_minus=Pi-Phi)
    row = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=alpha)
    for key, expected in (("r_s", v+np.zeros_like(u)), ("beta_minus", bm+np.zeros_like(u)),
                          ("beta_plus", bp), ("Pi", Pi+np.zeros_like(r)),
                          ("Phi", Phi+np.zeros_like(r)), ("f_s", f_s),
                          ("psi_s_normalized", f_s), ("N_minus_normalized", Nminus),
                          ("c_s", cs), ("n_c", nc)):
        np.testing.assert_allclose(row[key], expected, rtol=3e-14, atol=3e-14)
    # The common value rate equals f_T*T_s+f_r*r_s. This is the value
    # transported into radial points newly admitted by an expanding wall.
    f_T = A*np.exp(cap["delta"])*Pi
    np.testing.assert_allclose(f_T*row["T_s"]+Phi*row["r_s"], f_s, atol=2e-14)
    np.testing.assert_allclose(A*np.exp(2*cap["delta"])*row["T_s"]**2
                               -row["r_s"]**2/A, 1., atol=2e-14)
    assert np.all(row["T_s"] > 0) and np.all(row["H"] > 0)
    assert np.any(row["r_s"] < 0) and np.any(row["r_s"] > 0)
    # Exact transformed action, differentiated independently in Psi.
    tension_psi = np.sqrt(6)*cosh**2*np.sinh(q)
    potential_jump = -3*cosh**4
    Np = np.sqrt(alpha)*(nc-u*bp)/(r*cosh)
    Nm = np.sqrt(alpha)*Nminus
    np.testing.assert_allclose(Np-Nm, 2*tension_psi+np.zeros_like(r), atol=3e-14)
    np.testing.assert_allclose(row["beta_plus"]-row["beta_minus"], -tension*r, atol=2e-14)
    np.testing.assert_allclose(row["k_plus"]-row["k_minus"], -tension+np.zeros_like(r), atol=2e-14)
    force = tension*(row["k_plus"]+row["k_minus"]+2*(bp+bm)/r)
    force += tension_psi*(Np+Nm)-potential_jump
    np.testing.assert_allclose(force, 0., atol=3e-14)
    assert row["r_s"].shape == (2, 3)


def test_exterior_first_jets_and_both_coordinate_acceleration_equations():
    state, outgoing, cap, gauge = inputs()
    row = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.3)
    U, V = row["Udot_plus"], row["Vdot_plus"]
    assert gauge["C_wall"]*U*V == pytest.approx(1.)
    assert U*row["r_U_plus"]+V*outgoing["r_V"] == pytest.approx(row["r_s"])
    assert -U*row["r_U_plus"]+V*outgoing["r_V"] == pytest.approx(row["beta_plus"])
    assert U*row["c_U_plus"]+V*outgoing["c_V"] == pytest.approx(row["c_s"])
    assert -U*row["c_U_plus"]+V*outgoing["c_V"] == pytest.approx(row["n_c"])
    assert U*row["log_C_U_plus"]+V*outgoing["log_C_V"] == pytest.approx(gauge["log_C_s"])
    assert -row["eta_plus_s"]-gauge["log_C_s"]/2 == pytest.approx(
        -row["k_plus"]-U*row["log_C_U_plus"])
    assert row["eta_plus_s"]-gauge["log_C_s"]/2 == pytest.approx(
        row["k_plus"]-V*outgoing["log_C_V"])


def test_constant_exterior_null_rescaling_preserves_physical_map():
    state, outgoing, cap, gauge = inputs()
    old = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.3)
    a, b = 2.3, .7
    state["eta_plus"] += .5*np.log(b/a)
    gauge["C_wall"] /= a*b
    outgoing = {key: value/b for key, value in outgoing.items()}
    new = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.3)
    scales = {"Udot_plus": a, "Vdot_plus": b, "r_U_plus": 1/a,
              "c_U_plus": 1/a, "log_C_U_plus": 1/a}
    for key in old:
        assert new[key] == pytest.approx(scales.get(key, 1)*old[key], rel=2e-14, abs=2e-14)


def test_alpha_zero_with_nonzero_scalar_and_stationary_wall():
    state, outgoing, cap, gauge = inputs()
    Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["C_wall"])
    outgoing["r_V"] = (np.sqrt(cap["A"])-2*state["r"])/(2*Vp)
    zero = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=0.)
    small = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=1e-12)
    assert zero["r_s"] == pytest.approx(0., abs=2e-16)
    assert zero["c_s"] != 0 and zero["f_s"] != 0
    assert zero["f_s"] == pytest.approx(zero["psi_s_normalized"])
    assert zero["F"] == 1.
    for key in zero:
        assert small[key] == pytest.approx(zero[key], rel=1e-11, abs=1e-11)


def test_constraint_residual_selection_has_no_stationary_or_beta_zero_degeneracy():
    # From retained Codazzi/Gauss residuals, cap constraints=0 and exterior
    # outgoing constraint=0: J_plus=-v*d and J_plus=beta_minus*d, where
    # d=actual cap acceleration minus the predicted one. Their determinant
    # is -(v+beta_minus)=-H, not v or beta_plus.
    state, outgoing, cap, gauge = inputs()
    for v in (-.3, 0., .4):
        bm = np.sqrt(cap["A"]+v*v)
        bp = bm-2*state["r"]
        Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["C_wall"])
        outgoing["r_V"] = (v+bp)/(2*Vp)
        row = hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=0.)
        matrix = np.array([[1., row["r_s"]], [1., -row["beta_minus"]]])
        assert np.linalg.det(matrix) == pytest.approx(-row["H"])
        np.testing.assert_array_equal(np.linalg.solve(matrix, [0., 0.]), [0., 0.])


@pytest.mark.parametrize("invalid", [True, False, np.nan, np.inf, -np.inf, 1j,
                                     "0.5", [], [True], [True, 1.], None])
def test_malformed_numeric_inputs(invalid):
    state, outgoing, cap, gauge = inputs()
    cap["w_minus"] = invalid
    with pytest.raises(ValueError):
        hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.2)


@pytest.mark.parametrize("location,key,value", [
    (0, "r", 0.), (0, "r", -.1), (2, "A", 0.), (2, "A", -.1),
    (3, "C_wall", 0.), (3, "C_wall", -1.), (1, "r_V", -100.)])
def test_nonpositive_physical_margins_are_rejected(location, key, value):
    arguments = inputs()
    arguments[location][key] = value
    with pytest.raises(ValueError, match="require"):
        hybrid_wall_boundary(*arguments, alpha=.2)


@pytest.mark.parametrize("location,key,value", [
    (0, "eta_plus", 1000.), (0, "eta_plus", -1000.), (0, "c", 1e308),
    (2, "delta", 1000.), (2, "delta", -1000.), (1, "c_V", 1e308)])
def test_derived_overflow_and_underflow_are_rejected(location, key, value):
    arguments = inputs()
    arguments[location][key] = value
    with pytest.raises(ValueError, match="overflow|finite|positive"):
        hybrid_wall_boundary(*arguments, alpha=1.)


def test_timelike_margin_underflow_is_rejected():
    state, outgoing, cap, gauge = inputs()
    cap["A"], outgoing["r_V"] = 1e-300, 1e300
    with pytest.raises(ValueError, match="timelike margin"):
        hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=0.)


def test_exact_mapping_keys_shapes_alpha_and_no_input_mutation():
    for index, label in enumerate(("state", "outgoing", "cap", "gauge")):
        arguments = list(inputs())
        arguments[index] = {**arguments[index], "extra": 1.}
        with pytest.raises(ValueError, match=label):
            hybrid_wall_boundary(*arguments, alpha=.2)
        arguments[index] = None
        with pytest.raises(ValueError, match=label):
            hybrid_wall_boundary(*arguments, alpha=.2)
    state, outgoing, cap, gauge = inputs()
    before = copy.deepcopy((state, outgoing, cap, gauge))
    hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.2)
    assert (state, outgoing, cap, gauge) == before
    with pytest.raises(ValueError, match="require"):
        hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=-.1)
    state["r"], cap["A"] = [.5, .6], [.7, .8, .9]
    with pytest.raises(ValueError, match="broadcast"):
        hybrid_wall_boundary(state, outgoing, cap, gauge, alpha=.2)
