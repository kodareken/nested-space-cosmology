"""Independent projector geometry and incoming-vertex enclosure checks."""
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from recursive_horizons.nsc_incoming_projector_energy_bound import (
    rank_one_energy_identity, local_cross_product_coefficients,
    projected_energy_bound, incoming_riccati_interval_coefficients, RefinedDefectGeometry,
)
from recursive_horizons.nsc_incoming_source_tail import _riccati_at_one
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi
from recursive_horizons.nsc_transmitting_dirac_domain import S1, S2, S3

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']


def test_rank_one_bound_against_independent_matrix_pairing():
    assert rank_one_energy_identity()['transverse_squared_residual'] == '0'
    rng = np.random.default_rng(917)
    for _ in range(80):
        v = rng.normal(size=2)+1j*rng.normal(size=2); v /= np.linalg.norm(v)
        w = rng.normal(size=2)+1j*rng.normal(size=2); w /= np.linalg.norm(w)
        P, Q = np.outer(v, v.conj()), np.outer(w, w.conj())
        n = np.array([np.trace(P@S).real for S in (S1, S2, S3)])
        h = rng.normal(size=3); H = sum(x*S for x, S in zip(h, (S1, S2, S3)))
        d = np.linalg.norm(Q-P, 2)
        actual = abs(np.trace(H@(Q-P)))
        bound = 2*np.linalg.norm(np.cross(h, n))*d+2*abs(h@n)*d*d
        assert actual <= bound+3e-14


def test_local_directed_coefficients_enclose_existing_arbitrary_precision_owner():
    for group in (1, 13, 14, 32):
        c = CHANNELS[group]
        with _precision(40):
            intervals = incoming_riccati_interval_coefficients(c)
            values = _riccati_at_one(mp.mpf(c['compact_mass']), mp.mpf(c['angular_eigenvalue']))
            for interval, value in zip(intervals, values):
                assert _lo(interval.real) <= value.real <= _hi(interval.real)
                assert _lo(interval.imag) <= value.imag <= _hi(interval.imag)


def test_cross_polynomial_bounds_actual_owned_Hamiltonian_vertex():
    for group in (1, 13, 14, 32):
        c = CHANNELS[group]; bound = local_cross_product_coefficients(c)
        with mp.workdps(70):
            a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
            m, ell = mp.mpf(c['compact_mass']), mp.mpf(c['angular_eigenvalue'])
            positive = _riccati_at_one(m, ell)
            for sign in (1, -1):
                cs = _riccati_at_one(m, sign*ell)
                for j, value in enumerate(cs, 1):
                    assert abs(value-(-1)**j*mp.conj(positive[j-1])) < mp.mpf('1e-48')*max(1, abs(value)) if sign == -1 else True
                for E in (16, 40, 160):
                    x = mp.mpf(1)/(2*E); S = sum(value*x**j for j, value in enumerate(cs, 1))
                    u = a*a*abs(S)**2
                    n = mp.matrix([2*a*S.imag, -2*a*S.real, 1-u])/(1+u)
                    h = mp.matrix([-m, sign*ell/r, -E/a])
                    cross = [h[1]*n[2]-h[2]*n[1], h[2]*n[0]-h[0]*n[2], h[0]*n[1]-h[1]*n[0]]
                    upper = sum(mp.mpf(v)*x**j for j, v in enumerate(bound['cross_numerator_coefficient_upper'], 1))
                    assert mp.sqrt(sum(v*v for v in cross)) <= upper
        assert bound['constant_coefficient_exact'] == 0.


def test_projected_bound_improves_generic_group14_energy_bound():
    old = json.loads((ROOT/'results/development/nsc-incoming-vacuum-tail-bound.json').read_text())['groups'][13]
    previous = json.loads((ROOT/'results/development/nsc-incoming-middle-bound.json').read_text())['groups'][13]
    result = projected_energy_bound(CHANNELS[14], old, local_cross_product_coefficients(CHANNELS[14]), 16., 160.)
    assert result['density_error_upper'] < previous['vacuum_stress_error_upper'][0]/100
    assert result['quadratic_longitudinal_density_upper'] > 0
    assert result['vacuum_current_error'] == 0.
    assert not result['thermal_quadrature_low_band_errors_included']


def test_new_partition_preserves_directed_full_collar_coverage():
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    owner = RefinedDefectGeometry(config, intervals=32)
    with _precision(40):
        assert _hi(owner.cells[0][0]) >= _hi(owner.horizon)
        assert _lo(owner.cells[-1][0]) <= 1
        for earlier, later in zip(owner.cells, owner.cells[1:]):
            assert _hi(later[0]) >= _lo(earlier[0])
