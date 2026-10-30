"""Constraints, matter work and anchored-observer controls for the parent."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_cubic_operators import numeric_operator_matrices
from recursive_horizons.nsc_atlas_horizon_response import solve_source
from recursive_horizons.nsc_atlas_parent_backreaction import (
    C1,C2,anchor_point,calculation,energy_and_work,initial_metric,matter_source,
    solve_backreaction,tensor_gram)
from recursive_horizons.nsc_atlas_parent_circulation import Parent,VOLUME,initial_constraint
from recursive_horizons.nsc_atlas_s3_polynomials import X,gradient_dot,mean
from recursive_horizons.nsc_atlas_volume_feedback import gram_matrix,solve_feedback


@pytest.fixture(scope='module')
def record():return calculation()


def test_full_hamiltonian_scalar_and_vector_momentum_constraints(record):
    c=record['numerical_controls']
    assert c['Hamiltonian_max_absolute_residual']<1e-8
    assert c['scalar_momentum_max_absolute_residual']<1e-9
    assert c['vector_momentum_max_absolute_residual']<1e-10
    assert c['metric_two_integrator_max_gap']<1e-9
    assert c['cubic_two_integrator_max_gap']<1e-7


def test_complete_excitation_work_including_measure_and_volume(record):
    for r in record['energy_balances']:
        scale=abs(r['initial_excitation_energy_coefficient'])
        assert r['conservation_absolute_residual']/scale<1e-8
        assert abs(r['integrated_work_from_evolution']-r['independent_integrated_work'])/scale<1e-9
        assert r['parts']['volume_mode_energy']>0
        assert abs(sum(r['parts'].values())-r['final_excitation_energy_coefficient'])<1e-9


def test_physical_tensor_gram_and_trace_of_gradient_stress():
    G=gram_matrix();GT=tensor_gram(G);ops=numeric_operator_matrices()
    assert np.linalg.eigvalsh(GT)[0]>0
    z=Parent.retained().initial();z[2:6]=[.13,.27,-.41,.07]
    _,_,D,Q,_=matter_source(z)
    assert np.max(np.abs(ops['trace']@Q-D))<1e-14
    # A gradient dyad has squared norm (|Dphi|^2)^2 pointwise.
    assert abs(float(Q@GT@Q)-float(D@G@D))<1e-15


def test_initial_energy_agrees_with_independent_conformal_measure_expansion():
    p=Parent.retained();volume=solve_feedback(source=solve_source(endpoint_u=2.))
    ops=numeric_operator_matrices();G=gram_matrix()
    from recursive_horizons.nsc_atlas_cubic_operators import cubic_basis
    from recursive_horizons.nsc_atlas_s3_polynomials import inner_product
    P=cubic_basis();GP=np.array([[float(inner_product(a,b)) for b in P] for a in P])
    y=np.zeros(55);y[:11]=initial_metric(p)
    E=energy_and_work(volume,ops,G,GP,tensor_gram(G),0.,y)[0]
    c=initial_constraint(p)['coefficients'];z=p.initial();x1,x2,_,x4=X
    u=c[0]+c[1]*x1*x1+c[2]*(x2*x2+x4*x4)+c[3]*x2*x2*x4*x4+c[4]*x1*x2*x4
    phi=C1*z[2]*x1+C2*z[4]*x2*x4
    independent=VOLUME/p.H**3*float(mean(u*(3*p.mass_squared*phi**2+p.H**2*gradient_dot(phi,phi))))
    assert abs(E-independent)<1e-10


def test_anchored_scalar_phase_event_and_independent_integrator(record):
    r=record['observer']
    assert r['phase_speed_over_H']>0
    assert r['two_integrator_absolute_H_phase_shift_gap']<1e-8
    total=r['metric_feedback']['H_phase_event_proper_shift_over_epsilon_squared']+r['volume_feedback']['H_phase_event_proper_shift_over_epsilon_squared']
    assert total==pytest.approx(r['total_H_phase_event_proper_shift_over_epsilon_squared'])


def test_anchor_is_initial_scalar_maximum():
    p=Parent.retained();x=anchor_point(p);q=p.initial()
    a=C1*q[2];b=C2*q[4]
    ambient=np.array([a,b*x[3],0.,b*x[1]])
    intrinsic=ambient-x*np.dot(x,ambient)
    assert abs(np.dot(x,x)-1)<1e-14
    assert np.max(np.abs(intrinsic))<1e-14
    assert a*x[0]+b*x[1]*x[3]==pytest.approx(b/2+a*a/(2*b))


def test_same_two_mode_preparation_generates_additional_degrees(record):
    row=record['harmonic_generation'][0]['cubic_field_L2_squared_by_degree']
    assert all(row[str(l)]>1e-8 for l in (3,4,5,6))


def test_cross_protocol_relation_recovers_existing_guide_signal(record):
    root=Path(__file__).resolve().parents[1]
    prior=json.loads((root/'results/summaries/nsc-atlas-parent-circulation-v1.json').read_text(encoding='utf-8'))
    old=prior['light_comparisons'][0]['H_delay_over_epsilon_squared']
    r=record['cross_observable']
    assert r['H_guided_differential_light_delay_over_epsilon_squared']==pytest.approx(old,rel=2e-8)
    assert r['scalar_clock_shift_over_guided_light_delay']*r['H_guided_differential_light_delay_over_epsilon_squared']==pytest.approx(r['H_scalar_phase_event_shift_over_epsilon_squared'])


class EmptyParent(Parent):
    def initial(self,eigenvalues=(3.,8.)):
        z=super().initial(eigenvalues);z[2:6]=0.
        return z


def test_vacuum_recovers_zero_metric_volume_and_cubic_response():
    p=Parent.retained();empty=EmptyParent(p.vacuum,p.mass_squared,p.offset)
    volume=solve_feedback(source=solve_source(parent=empty,endpoint_u=2.))
    r=solve_backreaction(volume=volume)
    assert np.array_equal(r.solution.y,np.zeros_like(r.solution.y))


def test_compact_record_binds_current_sources():
    root=Path(__file__).resolve().parents[1]
    r=json.loads((root/'results/summaries/nsc-atlas-parent-backreaction-v1.json').read_text(encoding='utf-8'))
    for name,expected in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
