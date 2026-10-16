"""Independent wall-action, null-frame and observer conversion controls."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons import nsc_atlas_nonlinear_boundary as model


@pytest.fixture(autouse=True)
def bounded_test_cpu():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def inputs():
    state = dict(r=.5, c=.3, eta_plus=.2, eta_minus=-.3)
    outgoing = dict(plus=dict(r_V=.4, c_V=-.2, log_C_V=.17),
                    minus=dict(r_Qminus=-.3, c_Qminus=.7, log_C_Qminus=-.11))
    gauge = dict(plus=dict(C_wall=1.3, log_C_s=.09),
                 minus=dict(C_wall=.8, log_C_s=-.07))
    return state, outgoing, gauge


def test_exact_reference_wall_and_beta_plus_zero():
    state, outgoing, gauge = inputs()
    x = np.array([.45, .5, .16**(1/3), .55])
    # Independent signed vacuum angular solution and its proper acceleration.
    bp, bm = .2/x**2-1.25*x, .2/x**2+.75*x
    v = np.sqrt(bm*bm-(1-x*x))
    state.update(r=x, c=0.)
    Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["plus"]["C_wall"])
    Qm = np.exp(-state["eta_minus"])/np.sqrt(gauge["minus"]["C_wall"])
    outgoing["plus"].update(r_V=(v+bp)/(2*Vp), c_V=0.)
    outgoing["minus"].update(r_Qminus=(v-bm)/(2*Qm), c_Qminus=0.)
    result = model.wall_boundary_operator(state, outgoing, gauge, alpha=0.)
    for name, expected in (("r_s", v), ("beta_plus", bp), ("beta_minus", bm),
                           ("k_plus", -.4/x**3-1.25), ("k_minus", -.4/x**3+.75)):
        np.testing.assert_allclose(result[name], expected, rtol=2e-14, atol=2e-14)
    assert abs(result["beta_plus"][2]) < 2e-15
    assert all(np.all(np.isfinite(value)) for value in result.values())
    for name in ("c_s", "n_c", "psi_s_normalized", "N_plus_normalized", "N_minus_normalized"):
        np.testing.assert_array_equal(result[name], np.zeros(4))


def test_finite_amplitude_matches_canonical_action_jumps_and_force():
    state, outgoing, gauge = inputs()
    state.update(r=np.array([[.42], [.57]]), c=np.array([.3, -.4, .1]))
    alpha = np.array([.01, .2, .8])
    row = model.wall_boundary_operator(state, outgoing, gauge, alpha=alpha)
    r, c = state["r"], state["c"]
    psi = np.sqrt(6)*np.arcsinh(np.sqrt(alpha)*c/(np.sqrt(6)*r))
    cosh, sinh = np.cosh(psi/np.sqrt(6)), np.sinh(psi/np.sqrt(6))
    # Differentiate the transformed wall action S=2 cosh(Psi/sqrt6)^3.
    tension = 2*cosh**3
    tension_psi = np.sqrt(6)*cosh*cosh*sinh
    potential_jump = -3*cosh**4
    canonical_normals = []
    for side in ("plus", "minus"):
        normal = np.sqrt(alpha)*(row["n_c"]/r-c*row["beta_"+side]/r**2)/cosh
        canonical_normals.append(normal)
        np.testing.assert_allclose(np.sqrt(alpha)*row["N_"+side+"_normalized"], normal, atol=1e-14)
    Np, Nm = canonical_normals
    np.testing.assert_allclose(Np-Nm, 2*tension_psi, rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(row["beta_plus"]-row["beta_minus"], -tension*r, atol=2e-14)
    np.testing.assert_allclose(row["k_plus"]-row["k_minus"], -tension, atol=2e-14)
    force = tension*(row["k_plus"]+row["k_minus"]
                     +2*(row["beta_plus"]+row["beta_minus"])/r)
    force += tension_psi*(Np+Nm)-potential_jump
    np.testing.assert_allclose(force, 0., atol=2e-14)
    tangent = np.sqrt(alpha)*(row["c_s"]/r-c*row["r_s"]/r**2)/cosh
    np.testing.assert_allclose(np.sqrt(alpha)*row["psi_s_normalized"], tangent, atol=1e-14)
    assert row["r_s"].shape == (2, 3)


def test_incoming_jets_share_radius_scalar_lapse_and_tangent_normalization():
    state, outgoing, gauge = inputs()
    row = model.wall_boundary_operator(state, outgoing, gauge, alpha=.3)
    for side, first_rate, second_rate, r_first, r_second, c_first, c_second, L_first, L_second in (
        ("plus", row["Udot_plus"], row["Vdot_plus"], row["r_U_plus"], outgoing["plus"]["r_V"],
         row["c_U_plus"], outgoing["plus"]["c_V"], row["log_C_U_plus"], outgoing["plus"]["log_C_V"]),
        ("minus", row["Qminus_dot_minus"], row["Qplus_dot_minus"], outgoing["minus"]["r_Qminus"],
         row["r_Qplus_minus"], outgoing["minus"]["c_Qminus"], row["c_Qplus_minus"],
         outgoing["minus"]["log_C_Qminus"], row["log_C_Qplus_minus"]),
    ):
        assert gauge[side]["C_wall"]*first_rate*second_rate == pytest.approx(1.)
        assert first_rate*r_first+second_rate*r_second == pytest.approx(row["r_s"])
        assert -first_rate*r_first+second_rate*r_second == pytest.approx(row["beta_"+side])
        assert first_rate*c_first+second_rate*c_second == pytest.approx(row["c_s"])
        assert -first_rate*c_first+second_rate*c_second == pytest.approx(row["n_c"])
        assert first_rate*L_first+second_rate*L_second == pytest.approx(gauge[side]["log_C_s"])
        eta_s = row["eta_"+side+"_s"]
        # Both independently signed coordinate acceleration equations.
        assert -eta_s-gauge[side]["log_C_s"]/2 == pytest.approx(-row["k_"+side]-first_rate*L_first)
        assert eta_s-gauge[side]["log_C_s"]/2 == pytest.approx(row["k_"+side]-second_rate*L_second)


def test_constant_null_coordinate_rescaling_preserves_physical_map():
    state, outgoing, gauge = inputs()
    original = model.wall_boundary_operator(state, outgoing, gauge, alpha=.3)
    transforms = {"plus": (2.3, .7), "minus": (.6, 1.8)}
    for side, (a, b) in transforms.items():
        state["eta_"+side] += .5*np.log(b/a)
        gauge[side]["C_wall"] /= a*b
        divisor = b if side == "plus" else a
        outgoing[side] = {key: value/divisor for key, value in outgoing[side].items()}
    changed = model.wall_boundary_operator(state, outgoing, gauge, alpha=.3)
    coordinate_scales = {"Udot_plus": 2.3, "Vdot_plus": .7,
                         "Qminus_dot_minus": .6, "Qplus_dot_minus": 1.8,
                         "r_U_plus": 1/2.3, "c_U_plus": 1/2.3, "log_C_U_plus": 1/2.3,
                         "r_Qplus_minus": 1/1.8, "c_Qplus_minus": 1/1.8,
                         "log_C_Qplus_minus": 1/1.8}
    for key, value in original.items():
        assert changed[key] == pytest.approx(coordinate_scales.get(key, 1)*value, rel=2e-14, abs=2e-14)


def test_normalized_scalar_survives_alpha_zero_and_zero_radial_speed():
    state, outgoing, gauge = inputs()
    Vp = np.exp(state["eta_plus"])/np.sqrt(gauge["plus"]["C_wall"])
    Qm = np.exp(-state["eta_minus"])/np.sqrt(gauge["minus"]["C_wall"])
    outgoing["plus"]["r_V"] = -.3/(2*Vp)
    outgoing["minus"]["r_Qminus"] = -.7/(2*Qm)
    zero = model.wall_boundary_operator(state, outgoing, gauge, alpha=0.)
    small = model.wall_boundary_operator(state, outgoing, gauge, alpha=1e-12)
    assert zero["r_s"] == pytest.approx(0., abs=2e-16)
    assert zero["c_s"] != 0 and zero["N_plus_normalized"] != 0
    for key in zero:
        assert small[key] == pytest.approx(zero[key], abs=1e-11, rel=1e-11)


def test_central_observer_conversion_from_independent_regular_axis_data():
    alpha, p = .3, np.array([-.7, 0., .4])
    F = 1-alpha*p*p/6
    C, qm, qp = 1.4, .7, 1.8
    # Construct a regular axis in a separate parameter lambda. Its tangent
    # is timelike, dr/dlambda=0 and dc/dlambda=0, with invariant radial norm1.
    r_plus = .5*np.sqrt(C*qm/qp)
    c_plus = p*r_plus/np.sqrt(F)
    expected_rate = np.sqrt(C*qm*qp/F)
    row = model.central_readout(alpha=alpha, C=C, q_minus_rate=qm,
                                q_plus_rate=qp, c_Qplus=c_plus)
    np.testing.assert_allclose(row["T_rate"], expected_rate, rtol=2e-15)
    np.testing.assert_allclose(row["p"], p, rtol=2e-15, atol=1e-15)
    np.testing.assert_allclose(row["F"], F, rtol=2e-15)
    # Constant changes of both null labels and the axis parameter.
    a, b, rate = 2.3, .6, 1.7
    changed = model.central_readout(alpha=alpha, C=C/(a*b),
        q_minus_rate=a*qm/rate, q_plus_rate=b*qp/rate, c_Qplus=c_plus/b)
    np.testing.assert_allclose(changed["T_rate"], row["T_rate"]/rate, rtol=2e-15)
    for key in ("p", "F"):
        np.testing.assert_allclose(changed[key], row[key], rtol=2e-15, atol=1e-15)


def test_anchored_central_gauge_quiet_and_alpha_zero():
    alpha, cp = .2, np.array([-.4, 0., .7])
    F = 1-2*alpha*cp*cp/3
    row = model.central_readout(alpha=alpha, C=F, q_minus_rate=1., q_plus_rate=1., c_Qplus=cp)
    np.testing.assert_allclose(row["T_rate"], 1., atol=1e-15)
    np.testing.assert_allclose(row["p"], 2*cp, atol=1e-15)
    np.testing.assert_allclose(row["F"], F, atol=1e-15)
    quiet = model.central_readout(alpha=3., C=2., q_minus_rate=.5, q_plus_rate=4., c_Qplus=0.)
    assert quiet == dict(T_rate=2., p=0., F=1.)
    zero = model.central_readout(alpha=0., C=2., q_minus_rate=.5, q_plus_rate=4., c_Qplus=.3)
    assert zero == dict(T_rate=2., p=1.2, F=1.)


@pytest.mark.parametrize("invalid", [True, False, np.nan, np.inf, -np.inf, 1j, "0.5", [], [True], [True, 1.], None])
def test_non_real_or_nonfinite_inputs_are_rejected(invalid):
    state, outgoing, gauge = inputs()
    outgoing["minus"]["c_Qminus"] = invalid
    with pytest.raises(ValueError):
        model.wall_boundary_operator(state, outgoing, gauge, alpha=.2)
    with pytest.raises(ValueError):
        model.central_readout(alpha=.2, C=1., q_minus_rate=1., q_plus_rate=1., c_Qplus=invalid)


@pytest.mark.parametrize("section,key,value", [("state", "r", 0.), ("state", "r", -.1),
    ("plus", "C_wall", 0.), ("minus", "C_wall", -1.)])
def test_invalid_wall_domains(section, key, value):
    state, outgoing, gauge = inputs()
    (state if section == "state" else gauge[section])[key] = value
    with pytest.raises(ValueError, match="require"):
        model.wall_boundary_operator(state, outgoing, gauge, alpha=.2)


@pytest.mark.parametrize("key,value", [("alpha", -.1), ("C", 0.), ("C", -1.),
    ("q_minus_rate", 0.), ("q_plus_rate", -1.)])
def test_invalid_central_domains(key, value):
    kwargs = dict(alpha=.2, C=1., q_minus_rate=1., q_plus_rate=1., c_Qplus=.1)
    kwargs[key] = value
    with pytest.raises(ValueError, match="require"):
        model.central_readout(**kwargs)


def test_malformed_mappings_and_incompatible_shapes():
    state, outgoing, gauge = inputs()
    with pytest.raises(ValueError, match="state"):
        model.wall_boundary_operator({**state, "extra": 1}, outgoing, gauge, alpha=.2)
    with pytest.raises(ValueError, match="outgoing.plus"):
        model.wall_boundary_operator(state, {**outgoing, "plus": None}, gauge, alpha=.2)
    with pytest.raises(ValueError, match="gauge.minus"):
        model.wall_boundary_operator(state, outgoing, {**gauge, "minus": {}}, alpha=.2)
    state["r"], state["c"] = [.5, .6], [.1, .2, .3]
    with pytest.raises(ValueError, match="broadcast"):
        model.wall_boundary_operator(state, outgoing, gauge, alpha=.2)
    with pytest.raises(ValueError, match="broadcast"):
        model.central_readout(alpha=[0., .1], C=[1., 2., 3.], q_minus_rate=1., q_plus_rate=1., c_Qplus=.1)
    with pytest.raises(ValueError, match="require"):
        model.wall_boundary_operator(*inputs(), alpha=-.1)


@pytest.mark.parametrize("location,key,value", [("state", "eta_plus", 1000.),
    ("state", "c", 1e308), ("outgoing", "c_V", 1e308)])
def test_wall_overflow_is_rejected(location, key, value):
    state, outgoing, gauge = inputs()
    target = state if location == "state" else outgoing["plus"]
    target[key] = value
    with pytest.raises(ValueError, match="overflow|finite|positive"):
        model.wall_boundary_operator(state, outgoing, gauge, alpha=1.)


def test_central_overflow_and_input_mutation_are_controlled():
    with pytest.raises(ValueError, match="overflow|finite"):
        model.central_readout(alpha=1., C=1., q_minus_rate=1., q_plus_rate=1., c_Qplus=1e308)
    state, outgoing, gauge = inputs()
    before = copy.deepcopy((state, outgoing, gauge))
    model.wall_boundary_operator(state, outgoing, gauge, alpha=.2)
    assert (state, outgoing, gauge) == before
