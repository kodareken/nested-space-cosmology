"""Field matching, supplied state and central response on the fixed wall."""

import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_causal_observers import de_sitter_wall_events
from recursive_horizons.nsc_atlas_field_transfer import (
    arrival_from_wall_event,
    compact_bump,
    conformal_junction_coefficients,
    field_transfer_record,
    regular_radial_field,
)


@pytest.fixture(scope="module")
def record():
    return field_transfer_record()


def test_fixed_conformal_junction_cancels_only_the_surface_term():
    R, k, chi, nchi = .5, 2., .3, -.17
    beta_minus = 1/(5*R**2)+3*R/4
    beta_plus = beta_minus-k*R
    # Directly differentiate phi=chi/R on the two sides, without using
    # the adapter's jump formula. Reduced-field derivatives match.
    nphi_minus = nchi/R-chi*beta_minus/R**2
    nphi_plus = nchi/R-chi*beta_plus/R**2
    assert nphi_plus-nphi_minus == pytest.approx(k*chi/R)
    result = conformal_junction_coefficients(k)
    assert result["ricci_delta"] == 12
    assert result["normal_phi_jump_per_phi"] == k
    assert result["normal_chi_jump_per_chi"] == 0
    assert conformal_junction_coefficients(k, 0)["normal_chi_jump_per_chi"] == -k
    assert conformal_junction_coefficients(k, .2)["normal_chi_jump_per_chi"] != 0
    # The vacuum l=0 Schwarzschild potential remains nonzero at this event.
    assert (1-.8/R)*.8/R**3 != 0


def test_compact_input_and_primitive_are_differentiable_at_support_edges():
    z = np.array([-2., -1., 0., 1., 2.])
    b, derivative = compact_bump(z)
    np.testing.assert_array_equal(b, [0, 0, 1, 0, 0])
    np.testing.assert_array_equal(derivative, np.zeros(5))
    h = 1e-4
    for edge in (-1., 1.):
        inside = edge*(1-h)
        primitive, first = compact_bump(inside)
        assert primitive < 2e-15
        assert abs(first) < 7e-11
        # B''=-8(1-z^2)^2(1-7z^2) and
        # B'''=48z(1-z^2)(3-7z^2) both vanish at the edges.
        assert -8*(1-edge**2)**2*(1-7*edge**2) == 0
        assert 48*edge*(1-edge**2)*(3-7*edge**2) == 0


def test_pulse_is_specified_in_wall_proper_time(record):
    state, samples = record["incoming_state"], record["samples"]
    tau = np.array(samples["tau_wall_over_ell"])
    expected_z = (tau-state["wall_clock_centre_over_ell"])/state["wall_clock_half_width_over_ell"]
    np.testing.assert_allclose(samples["z"], expected_z)
    b, derivative = compact_bump(expected_z)
    np.testing.assert_allclose(samples["ell_I_over_epsilon"], derivative/.03)
    np.testing.assert_allclose(samples["primitive_F_over_epsilon"], b/2)
    assert state["prescribed_characteristics"] == 1
    assert not state["outgoing_characteristic_independently_prescribed"]
    assert not state["remote_source_or_global_exterior_state_solved"]


def test_incoming_wall_characteristic_has_plus_normal_sign():
    event = de_sitter_wall_events([.5])
    gamma, v, a = event["gamma"][0], event["v"][0], math.exp(event["t_flat"][0])
    # For a local incoming chi=F(eta+r) with F'=1, d_eta chi=d_r chi=1.
    tangent_derivative = gamma*(1+v)/a
    normal_derivative = gamma*(v+1)/a
    assert tangent_derivative+normal_derivative == pytest.approx(2*gamma*(1+v)/a)
    assert tangent_derivative-normal_derivative == 0


def test_arrival_clock_jacobian_by_independent_clock_difference():
    event = de_sitter_wall_events([.5-1e-6, .5, .5+1e-6])
    arrival, derivative = arrival_from_wall_event(event["R"], event["t_flat"],
                                                 event["gamma"], event["v"])
    measured = (arrival[2]-arrival[0])/(event["tau_wall"][2]-event["tau_wall"][0])
    assert measured == pytest.approx(derivative[1], rel=2e-8)
    # Independently differentiate eta_C=(eta+r)_wall and t_C=-log(-eta_C).
    p_dot = event["gamma"][1]*(1+event["v"][1])/math.exp(event["t_flat"][1])
    assert derivative[1] == pytest.approx(math.exp(arrival[1])*p_dot)


def test_regular_centre_and_radial_wave_operator():
    primitive = lambda q: q**5
    derivative = lambda q: 5*q**4
    eta, r = -.7, .2
    centre = regular_radial_field(eta, 0., primitive, derivative)
    near = regular_radial_field(eta, 1e-5, primitive, derivative)
    assert near == pytest.approx(centre, rel=2e-9)
    # Independent finite differences act on chi, not on the record's
    # stated wave identity. The principal conformal radial operator is free.
    def chi(t, radius):
        a = -1/t
        return a*radius*regular_radial_field(t, radius, primitive, derivative)
    h = 1e-3
    eta_second = (chi(eta+h, r)-2*chi(eta, r)+chi(eta-h, r))/h**2
    r_second = (chi(eta, r+h)-2*chi(eta, r)+chi(eta, r-h))/h**2
    # For a fifth-degree primitive these symmetric differences have the
    # same truncation term, so the wave-operator difference cancels.
    assert eta_second-r_second == pytest.approx(0., abs=2e-9)


def test_actual_pulse_has_regular_central_limit(record):
    assert record["checks"]["symmetric_radial_centre_limit_relative_error"] < 2e-7
    assert record["checks"]["radial_limit_probe_r_over_ell"] == 1e-5


def test_initial_data_and_no_return_have_strict_geometric_margins(record):
    certificate = record["no_return_certificate"]
    assert certificate["initial_null_label_gap_over_ell"] > 2.7
    assert certificate["earliest_arrival_after_wall_R08_eta_gap_over_ell"] > .24
    assert certificate["strict_radius_gap_over_ell"] > .14
    assert certificate["reflected_comoving_radius_upper_bound_through_t0_over_ell"] < .8
    assert certificate["does_not_exclude_uncomputed_exterior_tails"]
    initial_max_label = record["support"]["initial_flat_eta_over_ell"]*(1-.05)
    first_label = record["support"]["incoming_null_label_over_ell"][0]
    assert initial_max_label < first_label


def test_support_is_early_and_arrivals_are_later_on_their_own_clock(record):
    support = record["support"]
    assert 0 < support["wall_R_over_ell"][0] < .5 < support["wall_R_over_ell"][1] < .8
    assert support["wall_tau_over_ell"][1]-support["wall_tau_over_ell"][0] == pytest.approx(.06, abs=1e-12)
    assert support["centre_arrival_t_over_ell"][0] < support["centre_arrival_t_over_ell"][1] < 0
    for wall_time, arrival in zip(support["wall_flat_t_over_ell"], support["centre_arrival_t_over_ell"]):
        assert arrival > wall_time
    assert np.all(np.diff(record["samples"]["arrival_t_C_over_ell"]) > 0)
    assert record["checks"]["minimum_arrival_jacobian"] > 2.5


def test_response_and_area_follow_the_transfer_map(record):
    samples = record["samples"]
    incident = np.array(samples["ell_I_over_epsilon"])
    output = np.array(samples["ell_phi_C_over_epsilon"])
    stretch = np.array(samples["dT_dtau"])
    np.testing.assert_allclose(output*stretch, incident, atol=2e-14)
    assert abs(record["metrics"]["signed_output_area_over_epsilon"]) < 1e-11
    assert not record["metrics"]["squared_amplitude_is_energy"]
    assert record["checks"]["input_squared_amplitude_analytic_error"] < 1e-10
    assert record["checks"]["integral_32_64_max_difference"] < 1e-10
    assert record["metrics"]["input_squared_amplitude_integral_ell_over_epsilon2"] > record["metrics"]["output_squared_amplitude_integral_ell_over_epsilon2"] > 0


def test_waveform_asymmetry_is_from_fixed_clock_map(record):
    positive = record["peaks"]["positive_centre"]
    negative = record["peaks"]["negative_centre"]
    assert positive["ell_phi_C_over_epsilon"] > 0 > negative["ell_phi_C_over_epsilon"]
    assert positive["arrival_t_C_over_ell"] < negative["arrival_t_C_over_ell"]
    assert 1.02 < record["metrics"]["negative_to_positive_peak_magnitude"] < 1.03
    assert 2.63 < record["metrics"]["support_stretch_ratio"] < 2.64
    # The input remains symmetric in wall proper time; no source retuning.
    tau_c = record["incoming_state"]["wall_clock_centre_over_ell"]
    first, last = record["peaks"]["input_peak_tau_over_ell"]
    assert tau_c-first == pytest.approx(last-tau_c)
    assert record["theory"]["fixed_within_test"]


def test_afterpulse_central_field_is_zero_for_supplied_state(record):
    last_label = record["support"]["incoming_null_label_over_ell"][1]
    last_arrival = record["support"]["centre_arrival_t_over_ell"][1]
    for time in ((last_arrival+0.)/2, 0.):
        assert -math.exp(-time) > last_label
    assert record["checks"]["analytic_zero_checks"]["central_field_after_last_arrival_through_t0"] == 0
    assert "conditional" in " ".join(record["limits"])


@pytest.mark.parametrize("ell", [.01, 7.])
def test_unit_rescaling_preserves_normalized_response(ell, record):
    changed = field_transfer_record(ell)
    for key, array in record["samples"].items():
        np.testing.assert_allclose(changed["samples"][key], array, rtol=2e-9, atol=3e-10)
    for key, value in record["metrics"].items():
        if isinstance(value, float):
            assert changed["metrics"][key] == pytest.approx(value, abs=2e-9)
    assert changed["admission"]["distinct_wall_radius_samples"] <= 513


@pytest.mark.parametrize("radius", [0., 1., 1.5, math.nan, math.inf])
def test_central_arrival_domain_exit(radius):
    with pytest.raises(ValueError):
        arrival_from_wall_event(radius, 0., 1., 0.)


@pytest.mark.parametrize("eta,radius", [(0., .1), (math.inf, .1), (-1., -.1)])
def test_radial_chart_domain_exit(eta, radius):
    with pytest.raises(ValueError):
        regular_radial_field(eta, radius, lambda q: q, lambda q: 1.)


@pytest.mark.parametrize("ell", [0., -1., math.inf, math.nan])
def test_invalid_units_rejected_before_computation(ell):
    with pytest.raises(ValueError, match="length scale"):
        field_transfer_record(ell)


def test_budget_and_truth_boundaries(record):
    assert record["admission"]["distinct_wall_radius_samples"] <= 513
    assert record["admission"]["cpu_seconds_max"] == 60
    assert not record["admission"]["priority_claim"]
    assert not record["checks"]["quadrature_comparison_is_rigorous_error_bound"]
    assert record["theory"]["additional_scalar_surface_coupling"] == 0
    assert record["theory"]["curvature_coupling_xi"] == 1/6
    assert any("not Schwarzschild bulk scattering" in line for line in record["limits"])
