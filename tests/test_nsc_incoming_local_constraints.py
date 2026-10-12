"""Local Euler forces against inherited tensors and independent weak variation."""
import json
from math import comb
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_angular_stress import profile_jets
from recursive_horizons.nsc_horizon_source import conformal_stress
from recursive_horizons.nsc_incoming_cauchy_jets import (
    incoming_cauchy_jets, IncomingNormalJetChange,
)
from recursive_horizons.nsc_incoming_local_constraints import (
    incoming_local_constraints, local_euler_gradients, local_action_densities,
    scalar_taylor_coefficients, taylor_two_jets,
)
from recursive_horizons.nsc_light_restoration_action import LightRestorationAction
from recursive_horizons.nsc_magnetic_light_reference import (
    magnetic_light_spectrum, homogeneous_magnetic_restoration,
)
from recursive_horizons.nsc_spherical_local_history import LockedSphericalLocalAction
from recursive_horizons.nsc_spatial_reference_symbol import (
    SymbolJet, INDEX, homogeneous_seed_metric_jets,
)


def make_owners():
    root = Path(__file__).resolve().parents[1]
    ledger = json.loads((root/'results/development/nsc-subgap-history-response.json').read_text())['locked_inputs']
    return LightRestorationAction(magnetic_light_spectrum(4)), LockedSphericalLocalAction(ledger)


@pytest.fixture(scope="module")
def owners():
    return make_owners()


def geometric_tensor(theta, spectrum):
    """The already owned homogeneous restoration plus physical LLL geometry."""
    source = homogeneous_magnetic_restoration(theta, spectrum)
    W, W1, W2, _, _ = profile_jets(theta)
    r = 1/np.sin(theta); s1 = -np.cos(theta)/np.sin(theta)
    A = -r*r*W; Ap = -W1-2*W*s1
    App = -(W2+2*W1*s1+2*W*r*r)/(r*r)
    uu, uv, vv = conformal_stress(A, Ap, App, 0., 0., central_charge=spectrum.charge)
    factor = 1/(-A)/(4*np.pi*r*r)
    return np.array([source['rho']+(uu-2*uv+vv)*factor,
                     source['p_parallel']+(uu+2*uv+vv)*factor,
                     source['T01'], source['p_sphere']])


def tensor_gradients(fields, tensor):
    N, _, a, r = [float(f.value[0, 0, 0].real) for f in fields]
    rho, parallel, current, sphere = tensor
    return np.array([-4*np.pi*a*r*r*rho, -4*np.pi*a*a*r*r*current,
                     4*np.pi*N*r*r*parallel, 8*np.pi*N*a*r*sphere])


def inherited_light_check(owners):
    light, compact = owners
    domain = incoming_cauchy_jets()
    before = [f.data.copy() for f in domain.fields]
    values = local_euler_gradients(domain.fields, light, compact,
                                  varied=('N', 'beta', 'a', 'r'), samples=16)
    got = sum(v[0] for k, v in values.items() if k.startswith('light/'))
    expected = tensor_gradients(domain.fields, geometric_tensor(3*np.pi/4, light.spectrum))
    np.testing.assert_allclose(got, expected, rtol=0, atol=3e-11)
    assert all(np.array_equal(a, f.data) for a, f in zip(before, domain.fields))
    result = incoming_local_constraints(domain, light, compact)
    assert result['numerically_resolved']
    assert result['maximum_numerical_indicator'] < 3e-11
    assert result['Euler_bulk_identity_residual'] < 3e-11
    np.testing.assert_allclose(result['light_action_gradient'], expected[:2], rtol=0, atol=3e-11)
    np.testing.assert_array_equal(result['local_force'], -result['local_action_gradient'])
    assert not result['scope']['physical_state_reference_or_band_included']
    assert not result['scope']['constraint_roots_or_four_metric_equations_solved']
    return {"maximum_error": float(np.max(abs(got-expected))),
            "expected_four_gradients": expected.tolist(), "observed_four_gradients": got.real.tolist(),
            "local_result": result}


def test_incoming_local_forces_match_inherited_light_tensor_before_gauge_fixing(owners):
    inherited_light_check(owners)


def compact_neck_check(owners):
    light, compact = owners
    fields = homogeneous_seed_metric_jets()
    values = local_euler_gradients(fields, light, compact,
                                  varied=('N', 'beta', 'a', 'r'), samples=16)
    # nsc_charged_ctp_neck.charged_ctp_neck_source's unchanged compact tensor.
    c = compact.ledger['C_Weyl']
    tensor = c*np.array([-48., 112/3+16*np.pi, 0., -128/3-8*np.pi])
    expected = tensor_gradients(fields, tensor)
    np.testing.assert_allclose(values['compact/weyl_bulk'][0], expected, rtol=0, atol=3e-11)
    assert np.max(abs(values['compact/euler_bulk_diagnostic'])) < 3e-11
    return {"maximum_error": float(np.max(abs(values['compact/weyl_bulk'][0]-expected))),
            "expected_four_gradients": expected.tolist(),
            "observed_four_gradients": values['compact/weyl_bulk'][0].real.tolist()}


def test_locked_weyl_at_original_neck_matches_committed_compact_tensor(owners):
    compact_neck_check(owners)


def shifted_fields(coefficients, T, Z):
    """Translate a polynomial test germ exactly; no physical evolution."""
    c = coefficients[0]
    fields = []
    for f in range(4):
        out = SymbolJet.constant(np.zeros(len(T)))
        for dt in range(5):
            for dz in range(5-dt):
                coefficient = sum(c[f, t, z]*comb(t, dt)*comb(z, dz)
                                  *T**(t-dt)*Z**(z-dz)
                                  for t in range(dt, 5) for z in range(dz, 5-t))
                out.data[INDEX[dt, dz, 0]] = coefficient[:, None, None]*np.eye(2)
        fields.append(out)
    return tuple(fields)


def probe(x, length):
    # Value and first derivative vanish at both boundaries, as required by
    # the second-derivative action. This is a test function, not a history.
    u = x/length
    return ((1-u*u)**3, -6*u*(1-u*u)**2/length,
            (-6+36*u*u-30*u**4)/length**2)


def weak_variation_check(owners):
    light, compact = owners
    domain = incoming_cauchy_jets((
        IncomingNormalJetChange('a', 1, 1, .03),
        IncomingNormalJetChange('r', 1, 1, -.02),
        IncomingNormalJetChange('r', 2, 1, .01),
        IncomingNormalJetChange('a', 2, 0, .04),
    ))
    coefficients = scalar_taylor_coefficients(domain.fields)
    x, weights = np.polynomial.legendre.leggauss(10)
    length = .08
    t, z = np.meshgrid(length*x, length*x, indexing='ij')
    wt = length*length*weights[:, None]*weights
    ft, dt, ddt = probe(t, length); fz, dz, ddz = probe(z, length)
    test_jets = np.stack((ft*fz, dt*fz, ft*dz, ddt*fz, dt*dz, ft*ddz), axis=-1)
    raw = taylor_two_jets(coefficients, t, z)[0]
    direct = {}
    for B in range(2):
        varied = raw.copy()
        varied[..., B, :] += 1e-25j*test_jets
        densities = local_action_densities(varied, light, compact)
        for k, value in densities.items():
            direct.setdefault(k, np.zeros(2))[B] = np.sum(wt*value.imag/1e-25)
    # Integrate the independently extracted pointwise Euler coefficients.
    pointwise = {}
    flat_t, flat_z = t.ravel(), z.ravel()
    for start in range(0, t.size, 5):
        end = start+5
        fields = shifted_fields(coefficients, flat_t[start:end], flat_z[start:end])
        values = local_euler_gradients(fields, light, compact, samples=12)
        for k, value in values.items():
            pointwise.setdefault(k, []).append(value)
    errors = {}
    for k, parts in pointwise.items():
        values = np.concatenate(parts).reshape(*t.shape, 2)
        paired = np.sum(wt[..., None]*test_jets[..., :1]*values, axis=(0, 1))
        errors[k] = float(np.max(abs(paired-direct[k])))
        np.testing.assert_allclose(paired, direct[k], rtol=0, atol=3e-11, err_msg=k)
    # This development control resolves a nonzero local momentum force;
    # homogeneous zero-shift checks alone cannot certify the adapter.
    result = incoming_local_constraints(domain, light, compact)
    assert result['numerically_resolved']
    assert abs(result['local_action_gradient'][1]) > 1e-5
    return {"maximum_error": max(errors.values()), "channel_errors": errors,
            "quadrature": [10, 10], "test_interval": [-length, length],
            "changed_normal_entries": domain.changed_normal_entries(), "local_result": result}


def test_mixed_normal_jet_forces_equal_independent_compact_weak_variation(owners):
    weak_variation_check(owners)


def radius_check(owners):
    light, compact = owners
    fields = incoming_cauchy_jets().fields
    a = local_euler_gradients(fields, light, compact, samples=16)
    b = local_euler_gradients(fields, light, compact, samples=16,
                             coordinate_radius=.04, jet_radius=.02)
    errors = {k: float(np.max(abs(a[k]-b[k]))) for k in a}
    for k in a:
        np.testing.assert_allclose(a[k], b[k], rtol=0, atol=3e-11)
    return {"maximum_error": max(errors.values()), "channel_errors": errors,
            "numerical_radii": [[.03125, .015625], [.04, .02]]}


def test_contour_radius_change_and_input_contract(owners):
    radius_check(owners)
    light, compact = owners
    fields = incoming_cauchy_jets().fields
    with pytest.raises(ValueError, match='same intrinsic'):
        incoming_local_constraints(fields, light, compact)
    with pytest.raises(ValueError, match='distinct'):
        local_euler_gradients(fields, light, compact, varied=('N', 'N'))
    with pytest.raises(ValueError, match='eight'):
        local_euler_gradients(fields, light, compact, samples=4)
