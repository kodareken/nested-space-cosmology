import json
from pathlib import Path

import numpy as np
import sympy as sp

from recursive_horizons.nsc_incoming_finite_preparation import finite_endpoint_maps, symbolic_proof
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance


ROOT = Path(__file__).resolve().parents[1]


def test_finite_maps_agree_with_independently_owned_physical_tangent_coefficients():
    parent = json.loads((ROOT/'results/development/nsc-incoming-jet-rigidity.json').read_text())
    context = {name: sp.Symbol(name) for name in ('a1', 'r1', 'm', 'ell', 'A_k', 'R_k')}
    context.update(Delta_a_k=context['A_k'], Delta_r_k=context['R_k'], I=sp.I)
    for finite, tangent in zip(finite_endpoint_maps(), parent['result']['endpoint_maps']):
        for actual, old in zip(finite['finite_endpoint_coefficient'], (tangent['limit01'], tangent['limit10'])):
            assert sp.simplify(sp.sympify(actual, locals=context)-sp.sympify(old, locals=context)) == 0


def test_finite_second_normal_difference_keeps_nonzero_lower_jets():
    t = sp.Symbol('t')
    # Same nonzero first jets and unequal finite second jets; no small parameter.
    a0, r0 = 2+3*t+5*t*t/2, 7+11*t+13*t*t/2
    a1, r1 = a0+17*t*t/2, r0-19*t*t/2
    p1 = lambda a,r: a*(23+sp.I*29/r)/2
    difference = sp.diff(p1(a1,r1)-p1(a0,r0),t,2).subs(t,0)
    expected = ((23+sp.I*29/7)*17-sp.I*29*2*(-19)/49)/2
    assert sp.simplify(difference-expected) == 0


def test_actual_signed_source_preserves_complement_on_open_and_closed_fibers():
    for energy in (0.2, 9.0):
        positive = source_covariance(energy, 0.238, 3.97, 1.5)
        negative = source_covariance(-energy, 0.238, 3.97, 1.5)
        np.testing.assert_allclose(negative, np.eye(3)-positive.conj(), atol=2e-15, rtol=0)
    proof = symbolic_proof()
    assert proof['exact_scalar_residual_count'] == 33
    assert all(set(row['exact_residuals']) == {'0'} for row in proof['finite_endpoint_maps'])
