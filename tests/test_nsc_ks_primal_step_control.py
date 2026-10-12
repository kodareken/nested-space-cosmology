"""Primal-block DOP853 control is independent of concatenated tangent slots."""
import numpy as np
from scipy.integrate import DOP853

from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection)
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_primal_step_control import (
    PrimalBlockDOP853, evolve_primal_controlled, primal_initial_step)
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from tests.test_nsc_ks_difference_envelope import evolve, inputs


def _integrate(factory, n_extra):
    def fun(t, y):
        out = np.empty_like(y)
        out[0] = -y[0]
        if y.shape[0] > 1:
            out[1:] = -30.0 * y[1:] + 4.0
        return out

    y0 = np.zeros(1 + n_extra, complex)
    y0[0] = 1.0 + 0.2j
    solver = factory(fun, 0.0, y0, 1.0, rtol=1e-5, atol=1e-7, max_step=0.2)
    steps = 0
    while solver.status == 'running':
        assert solver.step() is None
        steps += 1
    return solver.y[0], steps


def test_extra_components_move_joint_steps_and_not_the_primal_controller():
    joint_plain, joint_steps = _integrate(DOP853, 0)
    joint_wide, joint_wide_steps = _integrate(DOP853, 80)
    plain, plain_steps = _integrate(lambda *args, **kwargs: PrimalBlockDOP853(*args, n_primal=1, **kwargs), 0)
    wide, wide_steps = _integrate(lambda *args, **kwargs: PrimalBlockDOP853(*args, n_primal=1, **kwargs), 80)
    assert joint_wide_steps != joint_steps
    assert plain_steps == joint_steps
    np.testing.assert_array_equal(plain, joint_plain)
    assert wide_steps == plain_steps
    np.testing.assert_allclose(wide, plain, rtol=0, atol=1e-15)


def test_full_block_initial_step_matches_the_unsplit_rule():
    def fun(t, y):
        return -y

    y0 = np.array([0.4 + 0.1j, -0.2])
    f0 = fun(0.0, y0)
    got = primal_initial_step(fun, 1.03, y0, 1.0, 0.002, f0, -1.0, 7, 2e-11, 2e-13, 2)
    direct = DOP853(fun, 1.03, y0, 1.0, rtol=2e-11, atol=2e-13, max_step=0.002)
    np.testing.assert_allclose(got, direct.h_abs, rtol=0, atol=0)


def test_controlled_difference_envelope_matches_across_tangent_count():
    before = evolve_ks_difference_envelope.__globals__['DOP853']
    one = evolve(owner=evolve_primal_controlled)
    source, initial, _family = inputs()
    w = LocalAxialFunction((0., 1.), E.axial_center())
    U = LocalAxialFunction((0.2, -0.1), E.axial_center())
    base = CompatibleRadiusDirection(w, U, .007, .03)
    extra = CompatibleRadiusDirection(
        LocalAxialFunction((0., 0., 1.), E.axial_center()),
        LocalAxialFunction((0.,), E.axial_center()), .007, .03)
    directions = (base,) + (extra,) * 15
    wide_family = CompatibleIncomingMetric((0.002,) + (0.,) * 15, directions)
    wide = evolve_primal_controlled(
        source, initial, wide_family, E.computational_z_grid(32),
        np.linspace(*E.physical_incoming_interval(), 9), np.pi / 2, np.sqrt(5), 1.03,
        axial_support=E.usual_axial_support(), rtol=2e-12, atol=2e-14, max_step=.001)
    assert wide.diagnostics['accepted_steps'] == one.diagnostics['accepted_steps']
    np.testing.assert_array_equal(wide.envelope_difference, one.envelope_difference)
    np.testing.assert_array_equal(wide.reference_amplitudes, one.reference_amplitudes)
    assert evolve_ks_difference_envelope.__globals__['DOP853'] is before
    assert one.diagnostics['step_control'] == 'primal'
    assert one.diagnostics['n_primal'] == wide.diagnostics['n_primal']
    assert wide.diagnostics['tangent_directions'] > one.diagnostics['tangent_directions']


def test_tangent_accuracy_control_cannot_relax_primal_steps():
    def fun(t, y):
        out = np.empty_like(y)
        out[0] = -y[0]
        if y.shape[0] > 1:
            out[1:] = -30.0 * y[1:] + 4.0
        return out

    y0 = np.zeros(81, complex)
    y0[0] = 1.0 + 0.2j
    y0[1:] = 0.05
    common = dict(fun=fun, t0=0.0, y0=y0, t_bound=1.0, n_primal=1, rtol=1e-5, atol=1e-7, max_step=0.2)

    def run(solver):
        steps = 0
        while solver.status == 'running':
            assert solver.step() is None
            steps += 1
        return solver.y[0], steps

    primal, primal_steps = run(PrimalBlockDOP853(**common))
    loose, loose_steps = run(PrimalBlockDOP853(tangent_rtol=1e-1, tangent_atol=1e-1, **common))
    tight, tight_steps = run(PrimalBlockDOP853(tangent_rtol=1e-12, tangent_atol=1e-14, **common))
    assert loose_steps == primal_steps
    np.testing.assert_allclose(loose, primal, rtol=0, atol=1e-13)
    assert tight_steps > primal_steps
