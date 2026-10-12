"""Regular mu-jets of the auxiliary PST projector; not a source or state run."""
import resource

import numpy as np
import pytest

pytest.importorskip('flint')
from flint import acb, arb, ctx

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_ks_spacetime_geometry_jets import spacetime_geometry_jets
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_scaled_reference_projector import (
    MatrixJet, ScalarJet, momentum_D, positive_inv, projector_derivative,
    scaled_reference_projector, scaled_star, symbol_metric_from_geometry,
)
from recursive_horizons.nsc_spatial_reference_symbol import reference_projector, star_order


S1 = (0, 1, 1, 0)
S3 = (1, 0, 0, -1)


def _as2(vec):
    return np.array([[complex(vec[0]), complex(vec[1])],
                     [complex(vec[2]), complex(vec[3])]], complex)


def _close(vec, mat, atol=2e-9):
    err = np.max(np.abs(_as2(vec) - np.asarray(mat, complex)))
    assert err < atol, err


def _cpu():
    return resource.getrusage(resource.RUSAGE_SELF).ru_utime


def _family(alpha=.001):
    w = LocalAxialFunction((.2, .3, -.1, .07), 1.2)
    U = LocalAxialFunction((.1, -.03, .05), 1.2)
    return CompatibleIncomingMetric((alpha,), (CompatibleRadiusDirection(w, U, .007, .03),))


def test_directed_inverse_does_not_round_upper():
    x = arb(2, .1)
    interior_small = float(x.mid()) - .5 * float(x.rad())
    interior_large = float(x.mid()) + .5 * float(x.rad())
    assert x.contains(interior_small) and x.contains(interior_large)
    inv = positive_inv(x)
    assert inv.contains(1 / interior_small) and inv.contains(1 / interior_large)
    wrong = arb(1 / float(x.upper()))
    assert not wrong.contains(1 / interior_small)
    with pytest.raises(ValueError, match='positive'):
        positive_inv(arb(0, .2))


def test_mixed_jet_product_keeps_cross_derivatives():
    T = ScalarJet.variable(arb(0), 0, 3, 2)
    z = ScalarJet.variable(arb(0), 1, 3, 2)
    mu = ScalarJet.variable(arb(0), 2, 3, 2)
    product = T * z * mu
    assert product.coefficient(1, 1, 1) == 1
    assert product.derivative(t=1, z=1, mu=1).value == 1
    with pytest.raises(ValueError, match='invent|exceeds'):
        product.trim(4, 2)
    with pytest.raises(ValueError, match='exceeds'):
        T.derivative(t=4)


def test_momentum_D_matches_exact_chain_rule_on_constant_symbol():
    B = MatrixJet.constant([acb(c) for c in S1], 2, 4)
    mu0 = arb('1/5')
    d1 = momentum_D(B, 2, 1, 1, mu0)
    # D[2,1]B = -e (2 B) = -2 S1 for e=+1 when d_mu B=0.
    _close(d1.value, -2 * np.array([[0, 1], [1, 0]], float))
    d2 = momentum_D(B, 2, 2, 1, mu0)
    # D[2,2]B = p(p+1) B = 6 S1.
    _close(d2.value, 6 * np.array([[0, 1], [1, 0]], float))
    d1m = momentum_D(B, 2, 1, -1, mu0)
    _close(d1m.value, 2 * np.array([[0, 1], [1, 0]], float))


def test_constant_geometry_higher_B_vanish_and_mu0_is_finite():
    kwargs = dict(mass=arb(1), angular=arb(2), axial=arb('11/10'), radius=arb('7/5'),
                  physical=4, momentum=4, bits=80)
    for mu, sign in ((arb('1/4'), 1), (arb(0), 1), (arb('1/3'), -1), (arb(0), -1)):
        jets = scaled_reference_projector(mu=mu, sign=sign, **kwargs)
        for entry in jets.B[0].value:
            assert entry.is_finite()
        for j in range(1, 5):
            for entry in jets.B[j].value:
                assert entry.is_finite()
                assert entry.contains(0)
        if mu.is_zero():
            _close(jets.P0.value, _as2(jets.Pinf.value), atol=1e-18)
            for j in range(1, 5):
                _close(jets.reconstruct_P_value(j), np.zeros((2, 2)), atol=1e-18)


def test_rejects_callable_metric_and_signed_mu_interval():
    with pytest.raises(TypeError, match='callables|ball'):
        scaled_reference_projector(mu=arb('1/5'), sign=1, mass=1, angular=1,
                                   axial=1, radius=lambda z, n: 1)
    with pytest.raises(ValueError, match='nonnegative'):
        scaled_reference_projector(mu=arb(0, .2), sign=1, mass=1, angular=1,
                                   axial=1, radius=2)
    with pytest.raises(ValueError, match=r'\+1 or -1'):
        scaled_reference_projector(mu=arb('1/5'), sign=True, mass=1, angular=1,
                                   axial=1, radius=2)
    with pytest.raises(TypeError, match='KSGeometryJets'):
        scaled_reference_projector(mu=arb('1/5'), sign=1, mass=1, angular=1,
                                   geometry={'a': 1})


def _point_geometry(order=4):
    return spacetime_geometry_jets(_family(), 1.012, 1.245, order=order, bits=120)


def _old_defects(row, n):
    P, H = row['jets'], row['hamiltonian']
    G = 0
    for i in range(n):
        for j in range(n):
            ell = n - i - j
            if ell >= 0:
                G = G + star_order(P[i], P[j], ell)
    F = 1j * P[n - 1].derivative(t=1)
    for j in range(n):
        F = F - (star_order(H, P[j], n - j) - star_order(P[j], H, n - j))
    return G, F


def test_finite_k_matches_old_owner_and_pst_powers_signs_factorials():
    geometry = _point_geometry(4)
    fields = symbol_metric_from_geometry(geometry)
    mass, angular = float(np.pi / 2), float(np.sqrt(5))
    momenta = [-2.0, -0.5, 0.4, 1.7]
    old = reference_projector(*fields, momenta, [mass] * 4, [angular] * 4)
    with ctx.workprec(120):
        for index, k in enumerate(momenta):
            sign = 1 if k > 0 else -1
            mu = arb(1 / abs(k))
            jets = scaled_reference_projector(
                mu=mu, sign=sign, mass=arb(mass), angular=arb(angular),
                geometry=geometry, physical=4, momentum=4, bits=120)
            for j in range(5):
                _close(jets.reconstruct_P_value(j), old['orders'][j, index], atol=5e-9)
            # k-derivative scaling against the unscaled owner, including P0.
            for h in (1, 2):
                _close(projector_derivative(jets, 0, k=h),
                       old['jets'][0].derivative(k=h).value[index], atol=5e-8)
            _close(projector_derivative(jets, 1, z=1),
                   old['jets'][1].derivative(z=1).value[index], atol=5e-8)
            _close(projector_derivative(jets, 1, t=1),
                   old['jets'][1].derivative(t=1).value[index], atol=5e-8)
            for n in range(1, 5):
                G, F = _old_defects(old, n)
                scale_g = mu ** (n + 2)
                scale_f = mu ** n
                _close([scale_g * c for c in jets.Gtilde[n - 1].value], G.value[index], atol=2e-8)
                _close([scale_f * c for c in jets.Ftilde[n - 1].value], F.value[index], atol=2e-8)
                for i in range(n):
                    for j in range(n):
                        ell = n - i - j
                        if ell < 0:
                            continue
                        unscaled = star_order(old['jets'][i], old['jets'][j], ell).value[index]
                        scaled = scaled_star(jets.B[i], jets.B[j], i, j, ell, sign, mu).value
                        power = mu ** (i + j + ell + 2)
                        _close([power * c for c in scaled], unscaled, atol=2e-8)


def test_scalar_trace_is_the_berry_term_not_forced_zero():
    geometry = _point_geometry(4)
    mass, angular, k = float(np.pi / 2), float(np.sqrt(5)), .35
    jets = scaled_reference_projector(
        mu=arb(1 / k), sign=1, mass=arb(mass), angular=arb(angular),
        geometry=geometry, physical=4, momentum=4, bits=120)
    fields = symbol_metric_from_geometry(geometry)
    old = reference_projector(*fields, [k], [mass], [angular])
    observed = complex(jets.trace_P[1])
    assert abs(observed - np.trace(old['orders'][1, 0])) < 5e-9
    a = float(geometry.coefficients['a'][(0, 0)])
    r = float(geometry.coefficients['r'][(0, 0)])
    rz = float(geometry.derivative('r', 0, 1))
    gap = np.sqrt(mass * mass + (angular / r) ** 2 + (k / a) ** 2)
    expected = -.5 * mass * angular * rz / (a * r * r * gap ** 3)
    assert abs(observed - expected) < 5e-9
    assert abs(expected) > 1e-8
    assert abs(complex(jets.trace_B[1])) > 1e-8


def test_actual_geometry_pilot_retains_nine_by_five(record_property):
    geometry = spacetime_geometry_jets(_family(), 1.012, 1.245, order=9, bits=160)
    mass, angular = arb.pi() / 2, arb(5).sqrt()
    start = _cpu()
    jets = scaled_reference_projector(
        mu=arb('1/4'), sign=1, mass=mass, angular=angular,
        geometry=geometry, physical=9, momentum=5, bits=160)
    cpu = _cpu() - start
    record_property('pilot_cpu_seconds', cpu)
    assert cpu < 90
    assert jets.retained == ((9, 5), (8, 4), (7, 3), (6, 2), (5, 1))
    for j, B in enumerate(jets.B):
        assert (B.physical, B.momentum) == jets.retained[j]
        for t, z, m in B.indices:
            for entry in B.coefficient(t, z, m):
                assert entry.is_finite()
    mixed = jets.B[0].coefficient(4, 5, 0)
    assert all(entry.is_finite() for entry in mixed)
    assert any(not entry.contains(0) for entry in mixed)
    p4 = projector_derivative(jets, 4, t=1, z=4)
    assert all(entry.is_finite() for entry in p4)
    assert any(not entry.contains(0) for entry in p4)
    assert any(not entry.contains(0) for entry in jets.B[4].coefficient(1, 4, 0))
    assert abs(complex(jets.trace_P[1])) > 1e-10
    assert jets.auxiliary and 'P_g + D_g' in jets.state_law
    assert jets.history_identity == geometry.history_identity


def test_actual_saved_control_box_reaches_infinite_momentum_without_a_pole():
    from recursive_horizons.nsc_ks_source_envelope import axial_center
    from recursive_horizons.nsc_ks_profile_identity import profile_identity
    center=axial_center()
    w=LocalAxialFunction((0.,1.),center);zero=LocalAxialFunction((0.,),center)
    family=CompatibleIncomingMetric((.001,),(CompatibleRadiusDirection(w,zero,.007,.03),))
    with ctx.workprec(160):
        mu=arb((1,-7),(1,-7)) # exact [0,1/64], not a rounded negative endpoint.
        geometry=spacetime_geometry_jets(family,arb(1.012,1e-7),arb(center+.045,1e-7),order=9,bits=160)
        point_geometry=spacetime_geometry_jets(family,1.012,center+.045,order=9,bits=160)
        for sign in (-1,1):
            box=scaled_reference_projector(mu=mu,sign=sign,mass=float(np.pi/2),angular=float(np.sqrt(5)),geometry=geometry)
            assert box.history_identity==profile_identity(family)
            for mu_point in (arb(0),arb(1)/128,arb(1)/64):
                point=scaled_reference_projector(mu=mu_point,sign=sign,mass=float(np.pi/2),angular=float(np.sqrt(5)),geometry=point_geometry)
                for j,t,z in ((0,0,0),(0,0,4),(4,1,4)):
                    a=box.B[j].derivative(t=t,z=z).value
                    b=point.B[j].derivative(t=t,z=z).value
                    assert all(x.is_finite() and x.contains(y) for x,y in zip(a,b))
