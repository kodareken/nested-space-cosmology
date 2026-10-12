"""Mutation controls on freshly evolved fields and their actual state tangent."""
import numpy as np

from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons import nsc_ks_source_envelope as E
from test_nsc_ks_difference_envelope import inputs


def evolved(amplitude, tangents):
    return evolve_ks_difference_envelope(
        *inputs(amplitude), E.computational_z_grid(32),
        np.linspace(*E.physical_incoming_interval(), 9),
        np.pi / 2, np.sqrt(5), 1.03, axial_support=E.usual_axial_support(),
        rtol=2e-12, atol=2e-14, max_step=0.001, tangents=tangents)


def contract(state, *, drop_response=False, diagonal=False, weight_factor=1.0, label_momentum=False):
    # column_weights already contains sqrt(dE/(2*pi)); apply it once.
    weights = np.asarray(state.column_weights) * np.sqrt(weight_factor)
    f, fz = state.columns * weights, state.axial_columns * weights
    df, dfz = state.column_tangents * weights, state.axial_tangents * weights
    if drop_response:
        df, dfz = np.zeros_like(df), np.zeros_like(dfz)
    if label_momentum:
        fz, dfz = -1j * state.source_energies * f, -1j * state.source_energies * df
    covariance = np.asarray(state.source_covariance)
    if diagonal:
        covariance = np.diag(np.diag(covariance))
    return source_column_matter(
        f, fz, covariance, df, dfz, mass=state.mass, angular=state.angular,
        axial_scale=E.reference_chart(1.0)[2], radius=np.sqrt(2), multiplicity=7.0)


def test_omitted_state_response_fails_derivative_of_fresh_evolution():
    h = 1e-4
    base = evolved(0.002, "all")
    plus, minus = evolved(0.002 + h, "zero"), evolved(0.002 - h, "zero")
    assert base.fixed_preparation_digest == plus.fixed_preparation_digest == minus.fixed_preparation_digest
    fd = (contract(plus)["action_gradient"] - contract(minus)["action_gradient"]) / (2*h)
    correct = contract(base)["action_gradient_tangent"][0]
    removed = contract(base, drop_response=True)["action_gradient_tangent"][0]
    np.testing.assert_allclose(correct, fd, rtol=0, atol=3e-8)
    # Both constraints must notice the missing retarded derivative.
    assert np.all(np.max(np.abs(removed-fd), axis=0) > 3e-8)


def test_evolved_constraint_values_detect_coherence_weights_and_wrong_momentum():
    state = evolved(0.002, "all")
    honest = contract(state)
    for options in ({"diagonal": True}, {"weight_factor": 2.0}, {"label_momentum": True}):
        changed = contract(state, **options)
        assert np.max(np.abs(changed["action_gradient"]-honest["action_gradient"])) > 1e-8
        assert np.max(np.abs(changed["action_gradient_tangent"]-honest["action_gradient_tangent"])) > 1e-8


def test_omitting_negative_energy_sector_changes_the_two_constraint_sum():
    from test_nsc_ks_signed_state import signed_sources, A_up, evolve, ANGULAR

    positive_source, negative_source = signed_sources()
    positive = evolve(positive_source, A_up(3), ANGULAR)
    negative = negative_angular_partner(positive, negative_source)
    assert np.all(positive.source_energies > 0)
    assert np.all(negative.source_energies < 0)
    pos, neg = contract(positive), contract(negative)
    summed = pos["action_gradient"] + neg["action_gradient"]
    dropped = pos["action_gradient"]
    assert np.all(np.max(np.abs(summed-dropped), axis=0) > 1e-8)
    complete_tangent = pos["action_gradient_tangent"] + neg["action_gradient_tangent"]
    assert np.max(np.abs(complete_tangent-pos["action_gradient_tangent"])) > 1e-8
