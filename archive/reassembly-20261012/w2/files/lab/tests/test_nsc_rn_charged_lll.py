"""Unit tests for the kappa = 0 charged lowest-Landau column."""

import numpy as np

from recursive_horizons import nsc_rn_charged_lll as lll
from recursive_horizons import nsc_rn_charged_window as window

# Manufactured dyonic sample from nsc-rn-charged-source.md, fixed digits.
A = 0.125
C_F = 0.5
Q_MAG = 1.0
Q_E = 4.0
MASS_RATIO = 2.0
SIGN = 1
OMEGA = 0.04063643378570953 / 2.0  # Phi_H / 2
PHI_H = 0.04063643378570953
EXPECTED_WINDOW_VALUE = -0.020318216892854766


def _reference():
    return lll.dyonic_reference(A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E,
                                mass_over_total_radius=MASS_RATIO)


def test_dyonic_reference_matches_note_digits():
    ref = _reference()
    assert abs(ref["Phi_H"] - PHI_H) < 1e-14
    assert abs(ref["r_m_tot"] - 1.0494385087475768) < 1e-12
    window_entry = window.window_gate(
        A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E,
        mass_over_total_radius=MASS_RATIO, sign=SIGN,
    )
    assert window_entry["window_open"] is True
    # The note's own manufactured sample: the window is open but every B_<
    # frequency is narrower than the declared budget width, so the gate is
    # closed. This test pins the *recorded* reading, not production physics.
    assert window_entry["admitted"] is False
    assert window_entry["reason"] == "every frequency in B_< is narrower than the declared budget width"
    assert min(1e-4, window_entry["cap_room"]) > 0.0
    assert ref["horizons"]["r_plus"] == window_entry["r_plus"]
    assert ref["horizons"]["r_minus"] == window_entry["r_minus"]


def test_electric_potential_closed_form():
    radius = np.array([1.0, 2.0, 4.0])
    phi = lll.electric_potential(radius, Q_E=Q_E, C_F=C_F)
    expected = Q_E / (16.0 * np.pi * C_F * radius)
    assert np.allclose(phi, expected, rtol=0.0, atol=1e-15)


def test_charged_grid_ends_and_cap():
    ref = _reference()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    assert abs(grid["outer"] - 32.0 * ref["r_m_tot"]) < 1e-12
    assert grid["inner_edge"] > ref["horizons"]["r_minus"]
    assert grid["radius"][0] == grid["inner_edge"]
    assert np.all(np.diff(grid["radius"]) > 0.0)


def test_horizon_jet_is_regular_and_in_window():
    ref = _reference()
    jet = lll.horizon_jet(model=ref["model"], omega=OMEGA, sign=SIGN,
                          Q_E=Q_E, C_F=C_F)
    assert jet["plus"] == 0.0
    assert abs(jet["effective_frequency_at_horizon"] - EXPECTED_WINDOW_VALUE) < 1e-14
    assert jet["window_member"] is True
    assert np.isfinite(jet["minus_derivative"])


def test_charged_lll_mode_carried_by_minus_only():
    ref = _reference()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=400)
    result = lll.charged_lll_mode(
        model=ref["model"], radius=grid["radius"], omega=OMEGA, sign=SIGN,
        Q_E=Q_E, C_F=C_F,
    )
    assert result["kappa"] == 0
    assert result["plus_identically_zero"] is True
    assert result["mode"].shape == (2, 400)
    assert np.isfinite(result["mode"]).all()
    assert np.max(np.abs(result["mode"][1])) > 1e-4


def test_charged_mode_form_drops_angular_row():
    ref = _reference()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=200)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_lll_mode(
        model=ref["model"], radius=grid["radius"], omega=OMEGA, sign=SIGN,
        Q_E=Q_E, C_F=C_F,
    )
    form = lll.charged_mode_form(
        result["mode"], chart=chart, weights=grid["weights"],
        derivative=grid["derivative"], radius=grid["radius"],
        Q_E=Q_E, C_F=C_F, sign=SIGN, omega=OMEGA,
    )
    assert form["angular_row_dropped_by_kappa"] is True
    assert form["window_entered"] is True
    assert np.isfinite(form["energy"])
    assert form["weighted_norm"] > 0.0


def test_normalized_column_carries_unit_charge_norm():
    ref = _reference()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=200)
    result = lll.charged_lll_mode(
        model=ref["model"], radius=grid["radius"], omega=OMEGA, sign=SIGN,
        Q_E=Q_E, C_F=C_F,
    )
    normed = lll.normalized_column(result["mode"], grid["weights"])
    weight_row = np.asarray(grid["weights"], dtype=float)
    norm = np.sqrt(np.real(np.sum(weight_row * np.sum(np.abs(normed["column"]) ** 2, axis=0))))
    assert abs(norm - 1.0) < 1e-12
    assert normed["phase_index"] == int(np.argmax(np.abs(result["mode"][1])))


def test_infall_winding_is_finite_and_ordered():
    ref = _reference()
    far = lll.infall_winding(
        model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F,
        r_outer=32.0 * ref["r_m_tot"], points=1200,
    )
    near = lll.infall_winding(
        model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F,
        r_outer=4.0 * ref["r_m_tot"], points=600,
    )
    assert far["T_total"] > near["T_total"] > 0.0
    assert 0.0 < far["Phi_mean"] < far["horizon_potential"]
    assert abs(far["omega_star"] - SIGN * far["Phi_mean"]) < 1e-15


def test_infall_phase_vanishes_at_omega_star():
    ref = _reference()
    entry = lll.infall_winding(
        model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F,
        r_outer=32.0 * ref["r_m_tot"], points=2400,
    )
    entry = lll.infall_winding(
        model=ref["model"], omega=entry["omega_star"], sign=SIGN, Q_E=Q_E, C_F=C_F,
        r_outer=32.0 * ref["r_m_tot"], points=2400,
    )
    assert abs(entry["total_phase_at_omega"]) < 1e-8


def test_infall_winding_quadrature_converges():
    ref = _reference()
    args = dict(model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F)
    coarse = lll.infall_winding(r_outer=32.0 * ref["r_m_tot"], points=600, **args)
    fine = lll.infall_winding(r_outer=32.0 * ref["r_m_tot"], points=4800, **args)
    # The near-horizon weight has a square-root endpoint, so the trapezoid
    # converges at first order there; pin the measured rate, not h^2.
    assert abs(coarse["Phi_mean"] - fine["Phi_mean"]) < 5e-4 * fine["Phi_mean"]


def test_infall_winding_omega_zero_is_negative():
    ref = _reference()
    entry = lll.infall_winding(
        model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F,
        r_outer=32.0 * ref["r_m_tot"],
    )
    assert entry["net_winding"] == "negative"
    assert entry["total_phase_at_omega"] < 0.0
    # At omega = 0 the effective frequency is negative all the way to the
    # cap, so the whole fall is inside the window and the fraction is 1.
    assert 0.99 < entry["negative_window_fraction"] <= 1.0 + 1e-12


def test_admitted_slice_winding_constants():
    """Pin the first admitted slice's phase-averaged potential."""
    ref = lll.dyonic_reference(
        A=0.04501936182826115, C_F=0.01125484045706527, q_mag=1.0, Q_E=4.0,
        mass_over_total_radius=1.01,
    )
    entry = lll.infall_winding(
        model=ref["model"], omega=0.0, sign=1, Q_E=4.0,
        C_F=0.01125484045706527, r_outer=32.0 * ref["r_m_tot"], points=2400,
    )
    assert abs(entry["T_total"] - 58.186762) < 1e-3
    assert abs(entry["omega_star"] - 0.237125) < 5e-4
    assert entry["omega_star"] < 0.8843561919019439  # admitted band centre
    assert entry["omega_star"] < 0.28216036056980787  # neutral frequency


ADMITTED_A = 0.04501936182826115
ADMITTED_C_F = 0.01125484045706527
ADMITTED_Q_E = 4.0
ADMITTED_ETA = 1.01
ADMITTED_OMEGA = 0.8843561919019439
TWIST_KAPPA = 2.0 ** 0.5


def _admitted():
    return lll.dyonic_reference(
        A=ADMITTED_A, C_F=ADMITTED_C_F, q_mag=1.0, Q_E=ADMITTED_Q_E,
        mass_over_total_radius=ADMITTED_ETA,
    )


def test_twisted_column_populates_both_rows():
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=400)
    result = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    assert result["kappa"] == TWIST_KAPPA
    assert result["plus_populated"] is True
    assert np.max(np.abs(result["mode"][0])) > 1e-6
    assert np.isfinite(result["mode"]).all()


def test_twisted_reflection_is_finite_and_subunity_here():
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    ratio = lll.reflection_ratio(result["mode"], chart=chart)
    assert ratio["ingoing_flux"] > 0.0
    assert ratio["outgoing_flux"] > 0.0
    assert 0.0 < ratio["reflection_ratio"] < 1.0


def test_kappa_zero_limit_gives_vanishing_reflection():
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=1.0e-9,
    )
    ratio = lll.reflection_ratio(result["mode"], chart=chart)
    assert ratio["reflection_ratio"] < 1e-12


def test_charge_flux_is_conservative_between_windows():
    """No amplification claim: reflection stays below unity at band centre."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=800)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    ratio = lll.reflection_ratio(result["mode"], chart=chart)
    # Sub-unity at the admitted band centre: the hole absorbs, nothing emits.
    assert ratio["reflection_ratio"] < 1.0


def test_opposite_sign_window_never_opens():
    """s = -1: E_eff = omega + Phi > 0 everywhere, so no window shell exists."""
    ref = _admitted()
    import numpy as _np
    radius = _np.linspace(ref["horizons"]["r_plus"], 32.0 * ref["r_m_tot"], 400)
    phi = lll.electric_potential(radius, Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F)
    effective = ADMITTED_OMEGA + phi
    assert _np.all(effective > 0.0)


def test_opposite_sign_pivot_is_smaller_at_band_centre():
    """Regular-root anatomy: |E_eff| = omega + Phi_H at the horizon, larger
    then the s = +1 value |omega - Phi_H|, so the pivot ratio is smaller."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=400)
    plus = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    minus = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=-1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    assert abs(minus["horizon_pivot"]["plus"]) < abs(plus["horizon_pivot"]["plus"])


def test_opposite_sign_reflection_finite_and_controlled():
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    result = lll.charged_twisted_column(
        model=ref["model"], radius=grid["radius"], omega=ADMITTED_OMEGA, sign=-1,
        Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
    )
    ratio = lll.reflection_ratio(result["mode"], chart=chart)
    assert 0.0 < ratio["reflection_ratio"] < 1.0


def test_window_floor_anomaly_reabsorbs():
    """Round-6 structure: at the window floor the s=+1 reflection collapses
    far below its mid-window plateau; no channel ever crosses unity."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    plateaus = {}
    for omega in (0.250000, 0.036600):
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=omega, sign=1,
            Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
        )
        plateaus[omega] = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
    assert plateaus[0.250000] > 0.5
    assert plateaus[0.036600] < 0.06
    assert plateaus[0.036600] < plateaus[0.250000]


def test_kappa_sign_symmetry_of_reflection():
    """Reflection is a flux ratio after column exchange, so +kappa and
    -kappa must reproduce the same ratio to integration accuracy."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    ratios = {}
    for kappa in (TWIST_KAPPA, -TWIST_KAPPA):
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=0.10, sign=1,
            Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=kappa,
        )
        ratios[kappa] = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
    assert abs(ratios[TWIST_KAPPA] - ratios[-TWIST_KAPPA]) < 2e-4


def test_reflection_resonance_is_bounded_and_inside_window():
    """The scanned maximum must stay sub-unity and sit inside [0.06, 0.14]."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    best = None
    for omega in (0.06, 0.08, 0.10, 0.12, 0.14):
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=omega, sign=1,
            Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=TWIST_KAPPA,
        )
        value = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
        if best is None or value > best[1]:
            best = (omega, value)
    assert 0.06 <= best[0] <= 0.14
    assert best[1] < 1.0


def test_omega_star_decreases_away_from_extremality():
    """phi_mean is horizon-weighted, so a farther horizon lowers omega_star."""
    ref_near = _admitted()
    ref_far = lll.dyonic_reference(
        A=ADMITTED_A, C_F=ADMITTED_C_F, q_mag=1.0, Q_E=ADMITTED_Q_E,
        mass_over_total_radius=1.5,
    )
    near = lll.infall_winding(
        model=ref_near["model"], omega=0.0, sign=1, Q_E=ADMITTED_Q_E,
        C_F=ADMITTED_C_F, r_outer=32.0 * ref_near["r_m_tot"], points=2400,
    )
    far = lll.infall_winding(
        model=ref_far["model"], omega=0.0, sign=1, Q_E=ADMITTED_Q_E,
        C_F=ADMITTED_C_F, r_outer=32.0 * ref_far["r_m_tot"], points=2400,
    )
    assert far["omega_star"] < near["omega_star"]


def test_higher_angular_ladder_step_stays_subunity():
    """kappa = sqrt(6) (n = 2) is the next multiplet; its reflection stays
    finite and below unity across the scan window."""
    ref = _admitted()
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=600)
    chart = lll.charged_chart(ref["model"], grid["radius"])
    best = None
    for omega in (0.06, 0.08, 0.10, 0.12, 0.14):
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=omega, sign=1,
            Q_E=ADMITTED_Q_E, C_F=ADMITTED_C_F, kappa=6.0 ** 0.5,
        )
        value = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
        if best is None or value > best[1]:
            best = (omega, value)
    assert best[1] < 1.0
    assert np.isfinite(best[1])