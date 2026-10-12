"""Charged light action versus the independently checked homogeneous tensor."""
import numpy as np
import pytest

from recursive_horizons.nsc_angular_stress import profile_jets
from recursive_horizons.nsc_horizon_source import conformal_stress
from recursive_horizons.nsc_magnetic_light_reference import magnetic_light_spectrum, homogeneous_magnetic_restoration
from recursive_horizons.nsc_light_restoration_action import (
    LightRestorationAction, CHANNELS, conformal_metric_jets, weyl_euler_average,
)
from recursive_horizons.nsc_spherical_local_history import SphericalHistoryGrid


def product(left, right):
    x, dx, ddx = left; y, dy, ddy = right
    return (x*y, dx*y[..., None]+x[..., None]*dy,
            ddx*y[..., None, None]+x[..., None, None]*ddy
            +dx[..., :, None]*dy[..., None, :]+dy[..., :, None]*dx[..., None, :])


def bump_axis(x, left, right):
    u = (2*x-left-right)/(right-left); scale = 2/(right-left); d = 1-u*u
    value = np.exp(1-1/d)
    log1 = -2*u/d**2*scale
    log2 = (-2/d**2-8*u*u/d**3)*scale**2
    return value, value*log1, value*(log2+log1*log1)


def homogeneous_variation(grid, amplitudes):
    """Old geometry in t=-theta, same future KS foliation; compact probes only."""
    t, z = grid.time, grid.radius_coordinate
    theta = -t
    W, W1, W2, _, _ = np.array([profile_jets(q) for q in theta]).T
    r = 1/np.sin(theta); s1 = -np.cos(theta)/np.sin(theta); s2 = r*r
    N, a = r/np.sqrt(W), r*np.sqrt(W)
    lnN1, lna1 = s1-W1/(2*W), s1+W1/(2*W)
    lnN2, lna2 = s2-.5*(W2/W-(W1/W)**2), s2+.5*(W2/W-(W1/W)**2)
    shape = (len(t), len(z)); dtype = np.result_type(amplitudes, float)
    def scalar(value, first, second):
        v = np.broadcast_to(value[:, None], shape).astype(dtype).copy()
        d = np.zeros((*shape, 2), dtype=dtype); dd = np.zeros((*shape, 2, 2), dtype=dtype)
        d[..., 0] = first[:, None]; dd[..., 0, 0] = second[:, None]
        return v, d, dd
    fields = [scalar(N, -N*lnN1, N*(lnN2+lnN1**2)),
              scalar(np.zeros_like(N), np.zeros_like(N), np.zeros_like(N)),
              scalar(a, -a*lna1, a*(lna2+lna1**2)),
              scalar(r, -r*s1, r*(s2+s1*s1))]
    # Endpoints are the declared test integration interval, not a physical
    # duration or selected geometry. Both the probe and all jets vanish there.
    time_range, space_range = (-2.4, -1.3), (-.7, .7)
    bt, dt, ddt = bump_axis(t, *time_range)
    bz, dz, ddz = bump_axis(z, *space_range)
    bz, dz, ddz = bz*(1+.2*z), dz*(1+.2*z)+.2*bz, ddz*(1+.2*z)+.4*dz
    envelope = bt[:, None]*bz
    de = np.stack((dt[:, None]*bz, bt[:, None]*dz), axis=-1)
    dde = np.empty((*shape, 2, 2))
    dde[..., 0, 0] = ddt[:, None]*bz
    dde[..., 0, 1] = dde[..., 1, 0] = dt[:, None]*dz
    dde[..., 1, 1] = bt[:, None]*ddz
    for i, amplitude in enumerate(amplitudes):
        fields[i] = tuple(v+amplitude*e for v, e in zip(fields[i], (envelope, de, dde)))
    NJ, BJ, AJ, RJ = fields
    nn, aa, bb = product(NJ, NJ), product(AJ, AJ), product(BJ, BJ)
    aabb, aab = product(aa, bb), product(aa, BJ)
    metric_jets = [[tuple(x-y for x, y in zip(nn, aabb)), tuple(-v for v in aab)],
                   [tuple(-v for v in aab), tuple(-v for v in aa)]]
    g = np.empty((*shape, 2, 2), dtype=dtype)
    dg = np.empty((*shape, 2, 2, 2), dtype=dtype)
    ddg = np.empty((*shape, 2, 2, 2, 2), dtype=dtype)
    for i in range(2):
        for j in range(2):
            g[..., i, j], dg[..., :, i, j], ddg[..., :, :, i, j] = metric_jets[i][j]
    return g, RJ[0], (dg, ddg, RJ[1], RJ[2]), envelope, (N, a, r)


def full_geometric_tensor(theta, spectrum):
    nonzero = homogeneous_magnetic_restoration(theta, spectrum)
    W, W1, W2, _, _ = profile_jets(theta)
    r = 1/np.sin(theta); s1 = -np.cos(theta)/np.sin(theta)
    A = -r*r*W; Ap = -W1-2*W*s1
    App = -(W2+2*W1*s1+2*W*r*r)/(r*r)
    uu, uv, vv = conformal_stress(A, Ap, App, 0., 0., central_charge=spectrum.charge)
    rho_lll = (uu-2*uv+vv)/(-A)/(4*np.pi*r*r)
    parallel_lll = (uu+2*uv+vv)/(-A)/(4*np.pi*r*r)
    return nonzero['rho']+rho_lll, nonzero['p_parallel']+parallel_lll, nonzero['p_sphere']


@pytest.mark.parametrize('charge', [0, 4])
def test_all_compact_metric_variations_match_homogeneous_restoration(charge):
    # Independent time/spatial refinement resolves the second-derivative
    # compact-bump quadrature; 64x48 has a measured O(1e-7) pairing error.
    grid = SphericalHistoryGrid.gaussian(128, 128, (-2.4, -1.3), (-.7, .7))
    spectrum = magnetic_light_spectrum(charge); owner = LightRestorationAction(spectrum)
    g, rgrid, jets, envelope, (N, a, r) = homogeneous_variation(grid, np.zeros(4))
    tensor = np.array([full_geometric_tensor(q, spectrum) for q in -grid.time])
    rho, parallel, sphere = tensor.T
    density = np.stack((-4*np.pi*a*r*r*rho, np.zeros_like(a),
                        4*np.pi*N*r*r*parallel, 8*np.pi*N*a*r*sphere), axis=-1)
    expected = grid.integral(density[:, None, :]*envelope[..., None])
    gradient = np.empty(4); step = 1e-28
    for B in range(4):
        amplitudes = np.zeros(4, complex); amplitudes[B] = 1j*step
        G, R, J, _, _ = homogeneous_variation(grid, amplitudes)
        actions, data = owner.actions(grid, G, R, J, canonical_coordinates='KS_time_reparametrization')
        assert tuple(actions) == CHANNELS
        gradient[B] = sum(actions.values()).imag/step
        assert data['scope']['LLL_geometry_included_once'] is True
        assert data['scope']['LLL_state_included'] is False
        assert data['scope']['compact_induced_action_included'] is False
    assert np.max(abs(gradient-expected)) < 3e-11


def test_weyl_path_jets_and_formal_quadrature_keep_spatial_dependence():
    grid = SphericalHistoryGrid.gaussian(12, 12, (-2.4, -1.3), (-.7, .7))
    g, r, jets, _, _ = homogeneous_variation(grid, np.array([.001, .002, .003, .001]))
    G, R, J = conformal_metric_jets(g, r, jets, 1.)
    assert np.max(abs(G-g)) == 0 and np.max(abs(R-r)) == 0
    assert all(np.max(abs(a-b)) == 0 for a, b in zip(J, jets))
    G, R, J = conformal_metric_jets(g, r, jets, 0.)
    assert np.max(abs(G*r[..., None, None]**2-g)) < 3e-14
    assert np.max(abs(R-1)) == 0 and np.max(abs(J[2])) == 0 and np.max(abs(J[3])) == 0
    three = weyl_euler_average(grid, g, r, jets)
    five = weyl_euler_average(grid, g, r, jets, quadrature_order=5)
    assert np.max(abs(three-five)/(1+abs(five))) < 3e-12
    owner = LightRestorationAction(magnetic_light_spectrum(4))
    with pytest.raises(ValueError, match='foliation'):
        owner.actions(grid, g, r, jets, canonical_coordinates='PG')
    with pytest.raises(ValueError, match='analytic'):
        owner.actions(grid, g, r, None, canonical_coordinates='KS')
