"""Full-PDE, known-limit and physical-unit checks for finite parent data."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

import recursive_horizons.nsc_atlas_finite_parent_data as m
from recursive_horizons.nsc_atlas_parent_circulation import Parent,VOLUME,initial_constraint


@pytest.fixture(scope='module')
def record():return m.calculation()


def test_full_unprojected_residual_contracts_across_three_resolutions(record):
    rows=record['levels'];r=[x['full_conformal_PDE_RMS_residual'] for x in rows]
    assert r[0]/r[1]>10 and r[1]/r[2]>10
    assert r[-1]<2e-4
    assert all(x['continuation'][-1]['projected_residual_norm']<1e-10 for x in rows)
    # A tiny weak residual cannot replace the residual outside the truncation.
    assert r[0]>1e7*rows[0]['continuation'][-1]['projected_residual_norm']


def test_positive_metric_and_finite_geometry_comparison(record):
    for r in record['levels']:
        assert 0<r['analytic_positivity_lower_bound']<r['sampled_Psi_min']
        assert r['sampled_Psi_max']<1
        v=r['proper_volume_over_reference_round_volume']
        assert r['volume_equivalent_length_over_reference_length']**3==pytest.approx(v)
        assert r['mean_excitation_density_over_vacuum']*v==pytest.approx(r['excitation_energy_over_reference_round_vacuum_energy'])
    assert record['convergence']['fine_middle_volume_spread']<2e-8


def test_exact_round_vacuum_limit():
    with threadpool_limits(limits=2):data=m.solve_initial_data(4,epsilon=0.)
    np.testing.assert_allclose(data.basis.evaluate(data.coefficients),1,atol=2e-14,rtol=0)


def test_homogeneous_source_recovers_independent_exact_radius(monkeypatch):
    f0=.03
    monkeypatch.setattr(m,'source_profile',lambda points,parent:(np.full(len(points),f0),np.zeros(len(points))))
    with threadpool_limits(limits=2):data=m.solve_initial_data(4,epsilon=.7)
    p=data.parent
    exact=(6/(6+p.kappa*.7**2*p.mass_squared/p.H**2*f0*f0))**.25
    np.testing.assert_allclose(data.basis.evaluate(data.coefficients),exact,atol=3e-12,rtol=0)


def test_small_amplitude_recovers_existing_quadratic_geometry():
    errors=[]
    with threadpool_limits(limits=2):
        for eps in (.04,.02):
            data=m.solve_initial_data(6,epsilon=eps)
            x,y,_,w=data.basis.points.T;c=initial_constraint(data.parent)['coefficients']
            u=c[0]+c[1]*x*x+c[2]*(y*y+w*w)+c[3]*y*y*w*w+c[4]*x*y*w
            errors.append(np.max(abs(data.basis.evaluate(data.coefficients)-(1+eps*eps*u))))
    assert 8<errors[0]/errors[1]<25


def test_existence_certificate_has_strict_nonconcentration_margin(record):
    c=record['finite_amplitude_existence']
    assert c['mass_squared_over_H_squared']>33
    assert c['constant_trial_over_bubbling_threshold']<1
    assert float(c['mean_A_at_epsilon1'])>0
    assert c['maximum_B_at_epsilon1']<30.5


def test_initial_expansion_is_not_future_collapse(record):
    r=record['levels'][-1]
    assert r['initial_normal_expansion_derivative_over_H_squared_at_anchor']==pytest.approx(3*(1+r['scalar_maximum_anchor_density_over_vacuum']))
    assert r['initial_normal_expansion_derivative_over_H_squared_at_anchor']>0


def test_source_bindings_are_current():
    root=Path(__file__).resolve().parents[1]
    r=json.loads((root/'results/summaries/nsc-atlas-finite-parent-data-v1.json').read_text(encoding='utf-8'))
    for name,expected in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name


def test_initial_generalized_harmonic_gauge_from_connection_definition():
    g,Pi,Phi=m.gh_initial_geometry(.81,np.array([.13,-.21,.07]))
    inverse=np.linalg.inv(g);derivatives=np.r_[-Pi[None,:,:],Phi]
    def contracted_connection(D):
        C=np.zeros((4,4,4))
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    C[a,b,c]=sum(inverse[a,d]*(D[b,d,c]+D[c,d,b]-D[d,b,c])/2 for d in range(4))
        return np.einsum('bc,abc->a',inverse,C)
    assert np.max(np.abs(contracted_connection(derivatives)))<1e-14
    assert np.array_equal(Pi[1:,1:],np.zeros((3,3)))
    # Zero extrinsic curvature does not mean zero mixed time-metric data.
    derivatives[0]=0.
    assert np.linalg.norm(contracted_connection(derivatives))>.1
