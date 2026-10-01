"""Actual projected metric jets and an independent Ricci contraction."""
import json
import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation
import derive_nsc_spherical_conformal_curvature as curvature
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_episode_assessment as metric
from recursive_horizons.nsc_spherical_coupling import _combine


def test_Qddot_is_the_actual_projected_Qrate_directional_derivative():
    grid = galerkin.build_grid(16, quadrature=64, length=2 * np.pi, gauge="conformal")
    phi0, phi1 = galerkin.manufactured_columns(16)
    state = galerkin.blank_state(grid, phi0, phi1)
    x = grid.xi_g
    state.Q = 1.2 + .1 * np.cos(7 * x)
    state.r = 1.1 + .08 * np.sin(2 * x)
    state.p_chi = .2 * np.sin(7 * x)
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    qt, qtt, before_projection = curvature.actual_q_second_rate(grid, state, rate, bundle)
    assert np.max(abs(qt - galerkin.prolong_geometry(grid, rate.Q))) < 1e-14
    for epsilon in (1e-6, 2e-7):
        plus = galerkin.rates(grid, _combine(state, rate, epsilon)).Q
        minus = galerkin.rates(grid, _combine(state, rate, -epsilon)).Q
        directional = galerkin.prolong_geometry(grid, (plus - minus) / (2 * epsilon))
        assert np.max(abs(qtt - directional)) < 2e-7
    assert np.max(abs(qtt - before_projection)) > 1e-3


def test_metric_scalar_matches_an_independent_Christoffel_contraction():
    count, length = 64, 8.
    x = np.arange(count) * length / count
    q = 1.1 + .1 * np.cos(2 * np.pi * x / length)
    qt = .07 * np.sin(4 * np.pi * x / length)
    qtt = -.03 * np.cos(2 * np.pi * x / length)
    zero = np.zeros(count)
    rh = metric.direct_rh_grid(q, q, zero, qt, qtt, length, L_dot=qt, beta_dot=zero)["R_h"]
    qx = metric.spectral_dx(q, length)
    qxx = metric.spectral_dx(qx, length)
    qtx = metric.spectral_dx(qt, length)
    for node in (0, 7, 19, 31, 47):
        signs = np.array([1., -1.])
        g = np.diag(signs * q[node] ** 2)
        first = np.zeros((2, 2, 2))
        second = np.zeros((2, 2, 2, 2))
        for i, sign in enumerate(signs):
            first[i, i, 0] = 2 * sign * q[node] * qt[node]
            first[i, i, 1] = 2 * sign * q[node] * qx[node]
            second[0, i, i, 0] = 2 * sign * (qt[node] ** 2 + q[node] * qtt[node])
            second[0, i, i, 1] = second[1, i, i, 0] = 2 * sign * (qt[node] * qx[node] + q[node] * qtx[node])
            second[1, i, i, 1] = 2 * sign * (qx[node] ** 2 + q[node] * qxx[node])
        assert abs(metric.ricci_scalar_from_christoffel(g, first, second) - rh[node]) < 1e-13


def test_saved_metric_curvature_is_bound_and_auxiliary_is_only_compared_afterward():
    record = json.loads(curvature.OUT.read_text())
    assert record["source_unchanged"] is True
    assert record["geometry_evolved"] is False
    assert record["auxiliary_substituted"] is False
    assert record["source_bindings"]["payload_v2"] == episode.sha256(continuation.NPZ)
    with np.load(curvature.NPZ, allow_pickle=False) as saved:
        for spec in episode.RUN_PLAN:
            name = spec["name"]
            actual = saved[name + "_final_R_h"]
            r = saved[name + "_final_r"]
            c2 = saved[name + "_final_weyl_C2"]
            expected = (actual - 2) ** 2 / (3 * r ** 4)
            assert np.max(abs(c2 - expected)) < 1e-14
            assert abs(np.max(c2) - record["results"][name]["final"]["weyl_C2_max"]) < 1e-14
            assert record["results"][name]["final"]["positive_chart"] is True
            assert record["results"][name]["maximum_auxiliary_gap"] > 0
    primary = record["results"]["nf512_dt_0.0005"]["final"]["weyl_C2_max"]
    time = record["results"]["nf512_dt_0.0010"]["final"]["weyl_C2_max"]
    space = record["results"]["nf256_dt_0.0005"]["final"]["weyl_C2_max"]
    assert abs(primary - time) / primary < .01
    assert abs(primary - space) / primary < .01
