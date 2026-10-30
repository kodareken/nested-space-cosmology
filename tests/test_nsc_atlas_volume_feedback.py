"""Causal, energetic and readout controls for the inherited volume response."""
import json
import hashlib
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_horizon_response import solve_source
from recursive_horizons.nsc_atlas_inherited_charge import harmonic_sector
from recursive_horizons.nsc_atlas_volume_feedback import (
    calculation,gram_matrix,solve_feedback,source_coefficients)


@pytest.fixture(scope='module')
def record():return calculation()


def test_two_methods_and_independently_integrated_energy(record):
    assert record['numerical_control']['coefficient_maximum_absolute_gap']<2e-10
    for r in record['integrated_energy_balances']:
        assert r['final_energy_over_epsilon_fourth']>0
        assert r['absolute_residual']<1e-7
        assert r['quadrature_error_estimate']<1e-6


def test_nonzero_future_memory_survives_analytic_tail(record):
    m=record['causal_memory']['rescaled_quadrature_and_tail']
    a=np.array(m['A']);b=np.array(m['B']);d=np.array(m['each_quadrature_future_error_upper_bound'])
    assert np.any(np.abs(a)>d) and np.any(np.abs(b)>d)
    assert m['phase_averaged_late_energy_estimate']>m['late_energy_future_error_upper_bound']>0


def test_gram_positive_and_shared_eigenvalue_block():
    G=gram_matrix()
    assert np.linalg.eigvalsh(G)[0]>0
    assert G[1,2]==pytest.approx(-1/24)
    assert G[1,1]==pytest.approx(1/16)
    assert G[2,2]==pytest.approx(1/12)


def wrapped_source(source,amplitude):
    def state(u):
        z=source.state(u);z[2:6]*=amplitude
        return z
    return SimpleNamespace(parent=source.parent,endpoint_u=.3,state=state)


def test_zero_source_and_quadratic_source_scaling():
    source=solve_source(endpoint_u=2.)
    zero=solve_feedback(source=wrapped_source(source,0.))
    assert np.array_equal(zero.solution.y,np.zeros_like(zero.solution.y))
    one=solve_feedback(source=wrapped_source(source,1.))
    two=solve_feedback(source=wrapped_source(source,2.))
    assert np.max(np.abs(two.coefficients(.3)[0]-4*one.coefficients(.3)[0]))<1e-10


def test_initial_acceleration_uses_actual_turning_point_not_phase_average():
    source=solve_source(endpoint_u=2.);r=solve_feedback(source=wrapped_source(source,1.))
    du=1e-5;c=r.coefficients(du)[0]
    expected=r.d*source.parent.mass_squared*source_coefficients(source,0.)/(2*r.kinetic*source.parent.H**2)
    assert np.max(np.abs(2*c/du**2-expected))<2e-5


def test_volume_changes_physical_mass_ratio_and_coupling_in_exact_lockstep(record):
    h=harmonic_sector();s0=h['s_star']
    def values(s):
        mu2=400*math.exp(-12*s);gap2=(15/11)*math.exp(-14*s)
        g=math.exp(-7*s)
        return g,math.sqrt(gap2/mu2)
    g0,m0=values(s0)
    for row in (record['samples'][2],record['samples'][4]):
        shift=.05**2*row['central_volume_shift_over_epsilon_squared']
        g,m=values(s0+shift)
        assert abs(math.log(g/g0)-7*math.log(m/m0))<2e-14


def test_retarded_energy_is_not_static_minimum(record):
    r=record['samples'][4]
    assert r['induced_volume_energy_over_epsilon_fourth']>0
    assert r['static_relaxation_energy_contribution_over_epsilon_fourth']<0
    assert abs(r['static_relaxation_energy_contribution_over_epsilon_fourth'])<abs(r['zero_gradient_local_relaxation_energy_over_epsilon_fourth'])
    # A positive source does not force a monotone or positive instantaneous response.
    assert r['spatial_mean_volume_shift_over_epsilon_squared']<0


def test_published_record_source_bindings():
    root=Path(__file__).resolve().parents[1]
    record=json.loads((root/'results/summaries/nsc-atlas-volume-feedback-v1.json').read_text(encoding='utf-8'))
    for name,expected in record['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
