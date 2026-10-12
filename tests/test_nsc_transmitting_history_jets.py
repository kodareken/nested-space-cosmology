"""New nonlinear metric chain rule and state-dependent Gaussian contraction."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_transmitting_history_jets import (
    FourthOrderModePropagator, raw_KS_amplitude_directions, generator_direction,
    evolve_field_jets, relative_branch_ctp_control,
)
from recursive_horizons.nsc_transmitting_dirac_domain import MODE_TO_CURRENT


def test_actual_metric_generator_jets_and_flux_identity():
    x = np.linspace(-2., 1.2, 81)
    owner = FourthOrderModePropagator(x, 1.3, 2.)
    amplitudes = np.array([.002, .001, .003, .001])
    provider = SuppliedKSHarmonicMetric(tuple(amplitudes), support_time=(0., .05))
    metric = provider.values(.025, x)
    directions = raw_KS_amplitude_directions(provider, .025, x, metric)
    assert owner.sbp_residual < 3e-11
    for field in range(4):
        tangent, residual = generator_direction(owner, metric, directions[field])
        assert residual < 3e-11
        h = 1e-5
        direction = h*np.eye(4)[field]
        Lp, _ = owner.generator(SuppliedKSHarmonicMetric(tuple(amplitudes+direction), support_time=(0., .05)).values(.025, x))
        Lm, _ = owner.generator(SuppliedKSHarmonicMetric(tuple(amplitudes-direction), support_time=(0., .05)).values(.025, x))
        diff = (Lp-Lm)/(2*h)-tangent
        assert (np.max(abs(diff.data)) if diff.nnz else 0.) < 3e-8
    wrong = directions[0].copy(); wrong[0, -1] = 1.
    with pytest.raises(ValueError, match='reference exterior'):
        generator_direction(owner, metric, wrong)


def test_nonlinear_field_tangent_against_independent_history_perturbation():
    mode = json.loads((ROOT/'results/development/nsc-pg-ctp-mode-jets.json').read_text())
    history = json.loads((ROOT/'results/development/nsc-transmitting-history-modes.json').read_text())
    name = '14_1'
    with np.load(ROOT/mode['payload']['path'], allow_pickle=False) as a:
        d = {k[len(name)+1:]: a[k] for k in a.files if k.startswith(name+'/')}
    with np.load(ROOT/history['payload']['path'], allow_pickle=False) as a:
        x = a[name+'/x']; char = a[name+'/reference_field']
    fields = char.reshape(2, len(x), len(d['energies']), 3).transpose(1, 2, 0, 3)
    fields = np.einsum('ab,nebc->neac', MODE_TO_CURRENT, fields)[::8]
    x = x[::8]
    owner = FourthOrderModePropagator(x, float(d['mass'].item()), float(d['angular'].item()))
    amplitudes = np.array([.002, .001, .003, .001]); times = np.linspace(0., .05, 9)
    provider = SuppliedKSHarmonicMetric(tuple(amplitudes), support_time=(0., .05))
    analytic = evolve_field_jets(owner, d['energies'], fields, times, provider)
    h = 1e-5; direction = h*np.eye(4)[2]
    parts = [evolve_field_jets(owner, d['energies'], fields, times,
              SuppliedKSHarmonicMetric(tuple(amplitudes+sign*direction), support_time=(0., .05)),
              with_jets=False)['evolved_field'] for sign in (1., -1.)]
    assert np.max(abs((parts[0]-parts[1])/(2*h)-analytic['tangent_fields'][2])) < 3e-8
    assert analytic['physical_endpoint_family'] is None
    assert analytic['full_spectral_covariance'] is None


def test_relative_branch_contraction_requires_the_physical_source_support():
    C = np.array([np.diag([.2, .8, 0.])])
    P = np.array([np.diag([1., 1., 0.])])
    A = np.array([1j*np.diag([1., -2., 0.])])
    result = relative_branch_ctp_control(np.array([.3]), np.array([.2]), C, P, A)
    assert result[0]['trace_residual'] < 3e-11
    A[0, 2, 2] = 1j
    with pytest.raises(ValueError, match='closed source fibers'):
        relative_branch_ctp_control(np.array([.3]), np.array([.2]), C, P, A)
