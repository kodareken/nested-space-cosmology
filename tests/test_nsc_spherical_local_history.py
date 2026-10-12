"""Focused checks of the extension of the existing local action domain."""
import numpy as np
import pytest

from recursive_horizons.nsc_spherical_local_history import (
    SphericalHistoryGrid, spherical_invariants, OwnedCompactKSHistory,
)


def limit_control(kind):
    """Apply the already owned static/homogeneous curvature formulas."""
    grid = SphericalHistoryGrid.gaussian(12, 12, (0., 1.), (-.5, .5))
    s = grid.time[:, None]+np.zeros((1, 12)) if kind == 'homogeneous' else np.zeros((12, 1))+grid.radius_coordinate
    N, q, r = 1+.1*s, 1+.2*s+.05*s*s, 1+.125*s+.025*s*s
    Np, qp, rp = .1+0*s, .2+.1*s, .125+.05*s
    Npp, qpp, rpp = 0*s, .1+0*s, .05+0*s
    g = np.zeros((*s.shape, 2, 2)); g[..., 0, 0] = N*N; g[..., 1, 1] = -q*q
    dg = np.zeros((*s.shape, 2, 2, 2)); ddg = np.zeros((*s.shape, 2, 2, 2, 2))
    dr = np.zeros((*s.shape, 2)); ddr = np.zeros((*s.shape, 2, 2))
    axis = 0 if kind == 'homogeneous' else 1
    dg[..., axis, 0, 0], dg[..., axis, 1, 1] = 2*N*Np, -2*q*qp
    ddg[..., axis, axis, 0, 0], ddg[..., axis, axis, 1, 1] = 2*(Np*Np+N*Npp), -2*(qp*qp+q*qpp)
    dr[..., axis], ddr[..., axis, axis] = rp, rpp
    got = spherical_invariants(grid, g, r, (dg, ddg, dr, ddr))
    if kind == 'homogeneous':
        A = (qpp-qp*Np/N)/(N*N*q)
        B = qp*rp/(N*N*q*r); C = (rpp-rp*Np/N)/(N*N*r); D = (1+(rp/N)**2)/r**2
        expected = {'R': -2*(A+2*B+2*C+D), 'C2': 4*(A-B-C+D)**2/3, 'E4': 8*(A*D+2*B*C)}
    else:
        A = (Npp-Np*qp/q)/(N*q*q)
        B = Np*rp/(N*q*q*r); C = (rpp-rp*qp/q)/(q*q*r); D = (1-(rp/q)**2)/r**2
        expected = {'R': 2*(A+2*B+2*C-D), 'C2': 4*(A-B-C-D)**2/3, 'E4': 8*(-A*D+2*B*C)}
    return {k: float(np.max(abs(got[k]-v))) for k, v in expected.items()}


@pytest.mark.parametrize('kind', ['static', 'homogeneous'])
def test_same_owned_curvature_limits(kind):
    assert max(limit_control(kind).values()) < 3e-11


def test_compact_variation_rejects_truncated_support():
    grid = SphericalHistoryGrid.gaussian(12, 12, (0., .025), (-1., 1.))
    with pytest.raises(ValueError, match='complete compact support'):
        OwnedCompactKSHistory(grid)


def test_analytic_metric_value_matches_value_only_path():
    grid = SphericalHistoryGrid.gaussian(12, 12)
    history = OwnedCompactKSHistory(grid)
    amplitudes = np.array([.002, .001, .003, .001])
    g, r, _ = history.fields_and_jets(amplitudes)
    G, R = history.fields(amplitudes)
    assert np.max(abs(g-G)) < 3e-14
    assert np.max(abs(r-R)) < 3e-14
