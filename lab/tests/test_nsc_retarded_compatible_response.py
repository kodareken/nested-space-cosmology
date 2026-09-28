"""New support/geometry controls and replay of the three saved response solves."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('retarded_response',ROOT/'scripts/derive_nsc_retarded_compatible_response.py')
PILOT=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(PILOT)


def test_compact_axial_bump_has_the_declared_coherent_derivatives():
    u=sp.Symbol('u');symbolic=sp.exp(1-1/(1-u*u))
    functions=[sp.lambdify(u,sp.diff(symbolic,u,n),'numpy') for n in range(4)]
    bump=PILOT.axial_bump(1.2,.06)
    for point in (-.8,-.4,0.,.3,.75):
        for n,function in enumerate(functions):
            expected=float(function(point))/.06**n
            assert abs(bump(1.2+.06*point,n)-expected) < 2e-11*(1+abs(expected))
    for point in (.9,1.,1.4,1.5):
        assert all(bump(point,n)==0 for n in range(4))


def test_retarded_support_and_cached_zero_geometry_match_owned_provider():
    guard=PILOT.support_guard();lower,upper=guard['radial_support'];start,stop=guard['PG_time_support']
    assert -2<lower<1<upper<1.2 and 0<start<stop<.3
    assert guard['maximum_checked_characteristic_speed']<0
    direction=PILOT.CompatibleRadiusDirection(PILOT.axial_bump(guard['axial_center']),lambda z,n:0.,.007,.03)
    x=np.array([-2.,lower,.99,1.,1.01,upper,1.2])
    provider=PILOT.CachedZeroRadiusProvider(x,direction)
    for time in (0.,.11,.15,.21,.3):
        metric=provider.values(time,x)
        assert np.array_equal(metric,provider.original.values(time,x))
        assert np.max(abs(provider.log_directions(time,x,metric)-provider.original.log_directions(time,x,metric)))<2e-16
    assert np.max(abs(provider.log_directions(0.,x,provider.metric)))==0
    assert np.max(abs(provider.log_directions(.3,x,provider.metric)))==0


def test_receipt_replays_arrays_without_field_or_radial_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('scientific producer called during receipt replay')
    monkeypatch.setattr(PILOT,'one_case',forbidden)
    monkeypatch.setattr(PILOT,'evolve_prepared_field_jets',forbidden)
    record=PILOT.check()
    assert len(record['control']['cases_completed'])<=3
    assert not record['scope']['four_energy_contribution_decides_C0']
    assert record['scope']['continuum_error_bound'] is None
    assert record['control']['amplitude']==0
    assert not record['scope']['metric_timestep']
    if not record['runtime']['cpu_budget_exceeded']:
        assert record['control']['cases_completed']==['401_64','801_64','801_128']
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:
        data={k:a[k] for k in a.files};meta=json.loads(data['metadata_json'].tobytes())
        summaries,comparisons,passed,_,_=PILOT.analyze(data,meta)
        assert record['status'].startswith('PASS' if passed else 'OPEN')
        assert comparisons==record['comparison']
        source=data['source/source'];P=data['source/projector'];weights=data['source/weights']
        assert np.array_equal(P,np.broadcast_to(np.diag([1,1,0]),(4,3,3)))
        assert np.max(abs(source[:,:,2]))==0 and np.all(weights>0)
        C=PILOT.block_diag(*source);root=np.repeat(np.sqrt(weights/(2*np.pi)),3)
        for case in record['control']['cases_completed']:
            F=(data[case+'/trace']*root).reshape(-1,12)
            dF=(data[case+'/trace_tangent']*root).reshape(-1,12)
            expected=dF@C@F.conj().T+F@C@dF.conj().T
            assert np.max(abs(expected-data[case+'/covariance_tangent']))<1e-15
            assert np.max(abs(data[case+'/trace_tangent'][0]))==0
            assert summaries[case]['diagnostics']['initial_tangent_norm']==0
            assert summaries[case]['diagnostics']['incident_tangent_norm']==0
            assert summaries[case]['diagnostics']['source_covariance_tangent_norm']==0
            assert not summaries[case]['diagnostics']['stationary_phase_reapplied']


def test_fixed_surface_constructor_correction_preserves_executed_maps_and_arrays():
    identity=PILOT.frame_constructor_correction_identity()
    assert identity['restriction_residual']==0
    assert identity['trace_map_triples_bitwise_equal']
    assert identity['restriction_map_bitwise_equal']
    record=PILOT.check()
    assert record['postprocessing_correction']['field_solves_repeated']==0
    assert record['execution_source_hashes']!=record['verification_source_hashes']
    assert not identity['coordinate_normal_or_metric_equality_claimed']
