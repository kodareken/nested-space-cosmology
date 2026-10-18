"""Constraint, regular-frame, geometric-derivative and intrinsic-operator controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_slice import (
    cap_slice, constant_rate_bridge, outgoing_slice, source_anchor_slice,
    spatial_dirac_coefficients,
)


def test_anchor_boost_and_shell_tensor_constraints():
    w = source_anchor_slice()
    p = w['parameters']
    bm = .4/(2*.45**2)+.75*.45
    bp = .4/(2*.45**2)-1.25*.45
    speed = math.sqrt(bm*bm-(1-.45**2))
    assert w['beta_minus'] == pytest.approx(bm)
    assert w['beta_plus'] == pytest.approx(bp)
    assert w['Rdot'] == pytest.approx(speed)
    boost = np.array([[w['gamma'], -w['gamma_v']], [-w['gamma_v'], w['gamma']]])
    assert boost @ np.diag([-1., 1.]) @ boost.T == pytest.approx(np.diag([-1., 1.]))
    assert boost @ np.array([speed, bm]) == pytest.approx([0., w['Rprime_minus']])
    assert boost @ np.array([speed, bp]) == pytest.approx([w['V_plus'], w['Rprime_plus']])
    # Project S_ab=-sigma*(g_ab-n_a*n_b), including the surface delta factor.
    n_cov = np.array([-w['gamma_v'], w['gamma']])
    stress = -p.sigma*(np.diag([-1., 1.])-np.outer(n_cov, n_cov))/w['gamma']
    jump_r = w['Rprime_plus']-w['Rprime_minus']
    jump_k = -w['V_plus']/w['R']
    assert -4*jump_r/w['R'] == pytest.approx(16*math.pi*stress[0, 0])
    assert -2*jump_k == pytest.approx(-8*math.pi*stress[0, 1])


def test_exact_cap_and_regular_minimum_differ_from_horizon():
    w = source_anchor_slice()
    cap = cap_slice(np.array([0., w['wall_distance']]), w)
    assert cap['R'] == pytest.approx([0., w['R']])
    assert cap['Rprime'][0] == 1
    assert cap['mass'][0] == 0
    assert cap['scalar_curvature'] == pytest.approx(6.)
    qh = math.acosh(math.sqrt(2*w['parameters'].mass/w['minimum_radius']))
    b = constant_rate_bridge(np.array([w['q_wall'], 0., qh]), w)
    assert b['Rprime'][0] < 0
    assert b['Rprime'][1] == 0
    assert b['Rsecond'][1] > 0
    assert b['R'][1] < w['R'] < b['R'][2]
    assert b['R'][2] == pytest.approx(.8)
    assert b['Rprime'][2] == pytest.approx(w['V_plus'])
    # Direct contraction in outgoing Eddington--Finkelstein coordinates,
    # including the point where Schwarzschild static time is singular.
    for i in range(3):
        f = 1-.8/b['R'][i]
        metric = np.array([[-f, -1.], [-1., 0.]])
        tangent = np.array([b['outgoing_EF_u_l'][i], b['Rprime'][i]])
        normal = np.array([b['outgoing_EF_normal_u'][i], b['V'][i]])
        assert tangent @ metric @ tangent == pytest.approx(1.)
        assert normal @ metric @ normal == pytest.approx(-1.)
        assert tangent @ metric @ normal == pytest.approx(0., abs=3e-15)


def test_constraints_by_independent_numerical_derivatives():
    w = source_anchor_slice()
    # Differentiate the parametric embedding rather than the algebraic RHS.
    q = np.array([-.16, .12, .35, .7])
    step = 1e-4
    def derivative(field):
        return (constant_rate_bridge(q-2*step, w)[field]
            -8*constant_rate_bridge(q-step, w)[field]
            +8*constant_rate_bridge(q+step, w)[field]
            -constant_rate_bridge(q+2*step, w)[field])/(12*step)
    b = constant_rate_bridge(q, w)
    lq = derivative('l')
    assert derivative('R')/lq == pytest.approx(b['Rprime'], abs=5e-12)
    rpp = derivative('Rprime')/lq
    scalar = -4*rpp/b['R']+2*(1-b['Rprime']**2)/b['R']**2
    hamiltonian = scalar+4*b['Kll']*b['Ktheta']+2*b['Ktheta']**2
    momentum = 2*b['Rprime']/b['R']*(b['Kll']-b['Ktheta'])-2*derivative('Ktheta')/lq
    assert hamiltonian == pytest.approx(0., abs=5e-10)
    assert momentum == pytest.approx(0., abs=5e-10)
    R = np.array([1.035, 1.08, 1.125, 1.165])
    data = outgoing_slice(R, w)
    def radial_derivative(field):
        return (outgoing_slice(R-2*step, w)[field]-8*outgoing_slice(R-step, w)[field]
            +8*outgoing_slice(R+step, w)[field]-outgoing_slice(R+2*step, w)[field])/(12*step)
    rpp = radial_derivative('Rprime')*data['Rprime']
    scalar = -4*rpp/R+2*(1-data['Rprime']**2)/R**2
    assert scalar+4*data['Kll']*data['Ktheta']+2*data['Ktheta']**2 == pytest.approx(0., abs=4e-8)
    momentum = 2*data['Rprime']/R*(data['Kll']-data['Ktheta'])-2*radial_derivative('Ktheta')*data['Rprime']
    assert momentum == pytest.approx(0., abs=4e-8)


def test_exact_tail_and_taper_overlap():
    w = source_anchor_slice()
    qend = math.acosh(math.sqrt(w['taper_start']/w['minimum_radius']))
    b = constant_rate_bridge(np.array([qend]), w)
    out = outgoing_slice(np.array([1., 1.2, 4.]), w)
    for key in ('R', 'Rprime', 'Rsecond', 'V', 'V_R', 'Kll', 'Ktheta'):
        assert b[key][0] == pytest.approx(out[key][0], abs=1e-14)
    for key in ('V', 'V_R', 'Kll', 'Ktheta'):
        assert out[key][1:] == pytest.approx(0.)
    assert out['Rprime'][1:] == pytest.approx(np.sqrt(1-.8/np.array([1.2, 4.])))
    assert out['mass'] == pytest.approx(.4)


def test_intrinsic_spinor_rescaling_square_and_wall_step():
    l = s.symbols('l', positive=True)
    r = s.Function('r')(l)
    chi = s.Matrix([l+l**3, 1+2*l*l])
    assert all(s.simplify(v)==0 for v in
        r*((chi/r).diff(l)+s.diff(r, l)/r*(chi/r))-chi.diff(l))
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    sigma2 = s.Matrix([[0, -s.I], [s.I, 0]])
    W = 1/r
    def D(v): return -s.I*sigma2*v.diff(l)+sigma1*W*v
    square = -chi.diff(l, 2)+s.diag(W*W-s.diff(W, l), W*W+s.diff(W, l))*chi
    assert all(s.simplify(v)==0 for v in D(D(chi))-square)
    w = source_anchor_slice()
    for k in (-2, 1, 3):
        inside = spatial_dirac_coefficients(w['R'], w['Rprime_minus'], angular_k=k)
        outside = spatial_dirac_coefficients(w['R'], w['Rprime_plus'], angular_k=k)
        assert outside['W'] == inside['W']
        assert outside['W_l']-inside['W_l'] == pytest.approx(k*2*w['gamma']/w['R'])
        flat = spatial_dirac_coefficients(.3, 1., angular_k=k)
        assert flat['V_plus'] == pytest.approx(k*(k-1)/.3**2)


def test_length_scale_covariance_of_joined_geometry():
    base, scaled = source_anchor_slice(), source_anchor_slice(7.)
    for key in ('R', 'minimum_radius', 'minimum_distance', 'wall_distance', 'taper_start', 'taper_stop'):
        assert scaled[key] == pytest.approx(7*base[key])
    for key in ('Rdot', 'gamma', 'gamma_v', 'Rprime_plus', 'Rprime_minus', 'V_plus'):
        assert scaled[key] == pytest.approx(base[key])
    q = np.array([-.1, 0., .4])
    b, c = constant_rate_bridge(q, base), constant_rate_bridge(q, scaled)
    assert c['R'] == pytest.approx(7*b['R'])
    assert c['Ktheta'] == pytest.approx(b['Ktheta']/7)
    assert c['scalar_curvature'] == pytest.approx(b['scalar_curvature']/49)


@pytest.mark.parametrize('scale', [0., -1., math.inf, math.nan])
def test_reject_invalid_length(scale):
    with pytest.raises(ValueError): source_anchor_slice(scale)


def test_reject_outside_declared_domains():
    w = source_anchor_slice()
    for q in (w['q_wall']-.01, 2., math.inf):
        with pytest.raises(ValueError): constant_rate_bridge(q, w)
    with pytest.raises(ValueError): cap_slice(-.1, w)
    with pytest.raises(ValueError): outgoing_slice(.9, w)
    with pytest.raises(ValueError): spatial_dirac_coefficients(0., 1.)
    with pytest.raises(ValueError): spatial_dirac_coefficients(.5, 1., angular_k=True)
