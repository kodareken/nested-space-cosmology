"""The probe inherits the actual dense background and canonical mode scale."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_boson_hamiltonian import potential_hessian
from recursive_horizons.nsc_atlas_boson_parent import ParentBosonProbe
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import ParentGH,evolve_parent


@pytest.fixture(scope='module')
def vacuum_history():
    history=DenseGHHistory()
    result=evolve_parent(degree=2,epsilon=0.,threads=1,proper_endpoint=.16,
        budget_seconds=90,dense_sink=history.record)
    point=result['snapshots'][-1]
    history.metadata=dict(degree=2,epsilon=0.,gauge='lapse_controlled',backend='hopf',quadrature_order=None,
        physical_coordinate_endpoint=point['coordinate_t'],proper_endpoint=point['observer_proper_duration'])
    return history


def test_adapter_never_recalculates_initial_data_and_uses_the_exact_grid(vacuum_history,monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('initial data must not be recalculated')
    monkeypatch.setattr(ParentGH,'initial',forbidden)
    probe=ParentBosonProbe(vacuum_history,2)
    np.testing.assert_array_equal(probe.basis.points,probe.system.basis.points)
    np.testing.assert_array_equal(probe.basis.weights,probe.system.basis.weights)
    np.testing.assert_allclose(probe.Y.T@(probe.basis.weights[:,None]*probe.Y),np.eye(probe.n),atol=3e-15,rtol=0)
    assert np.all(probe.fields[probe.s_indices]==0) and np.all(probe.fields[probe.chi_indices]==1)
    assert probe.phase_dimension==2*probe.n


def test_round_vacuum_matches_known_mode_hamiltonian_without_cross_field_blocks(vacuum_history):
    probe=ParentBosonProbe(vacuum_history,2)
    masses=potential_hessian([0.,0.],sstar=probe.system.sstar).diagonal()
    for t in (0.,.07,probe.coordinate_endpoint):
        G=probe.hamiltonian(t);a=math.cosh(probe.system.r*t)
        degrees=np.array([label[0] for label in probe.labels])
        expectedK=a*probe.system.r**2*degrees*(degrees+2)+a**3*masses[probe.fields]
        np.testing.assert_allclose(G[:probe.n,:probe.n],np.diag(expectedK),rtol=3e-11,atol=3e-10)
        np.testing.assert_allclose(G[probe.n:,probe.n:],a**-3*np.eye(probe.n),rtol=3e-12,atol=3e-12)
        np.testing.assert_allclose(G[:probe.n,probe.n:],0.,atol=2e-13,rtol=0)


def test_unit_kicks_and_observer_rows_keep_the_physical_harmonic_normalization(vacuum_history):
    probe=ParentBosonProbe(vacuum_history,2);rows=probe.observer_rows(0.)
    point=vacuum_history.evaluate(0.)[probe.system.dof:probe.system.dof+4]
    expected=[math.sqrt(2)/math.pi*point[0],math.sqrt(12)/math.pi*point[1]*point[3]]
    kicks=[probe.source_kick(number) for number in (1,2)]
    for number,kick,profile in zip((1,2),kicks,expected):
        np.testing.assert_array_equal(kick[:probe.n],0.)
        np.testing.assert_array_equal(kick[probe.n+probe.s_indices],0.)
        assert np.linalg.norm(kick)==pytest.approx(1.,abs=2e-16)
        assert rows['chi']@kick[probe.n:]==pytest.approx(probe.system.r**1.5*profile,abs=2e-16)
        assert rows['S']@kick[probe.n:]==0.
        assert probe.source_projection_residuals[number]<2e-15
    assert kicks[0]@kicks[1]==pytest.approx(0.,abs=2e-16)
    np.testing.assert_array_equal(rows['log_g_KK'],-7/math.sqrt(84)*rows['S'])
    with pytest.raises(ValueError):probe.source_kick(3)


def test_proper_checkpoints_and_physical_endpoint_do_not_extrapolate(vacuum_history):
    probe=ParentBosonProbe(vacuum_history,2)
    proper=np.array([point['observer_proper_duration'] for point in probe.proper_checkpoints])
    np.testing.assert_allclose(proper,np.linspace(probe.proper_start,probe.proper_endpoint,5),rtol=0,atol=2e-12)
    assert np.all(np.diff(probe.checkpoint_times)>0)
    for t in (-1e-6,probe.coordinate_endpoint+1e-9,float('nan')):
        with pytest.raises(ValueError,match='physical background interval'):probe.hamiltonian(t)
    if vacuum_history.bounds[-1][1]>probe.coordinate_endpoint:
        in_archive=(vacuum_history.bounds[-1][1]+probe.coordinate_endpoint)/2
        vacuum_history.evaluate(in_archive)
        with pytest.raises(ValueError):probe.observer_rows(in_archive)


def test_decoupled_control_removes_only_cross_field_stiffness(vacuum_history):
    original=ParentBosonProbe(vacuum_history,2);state=vacuum_history.evaluate(0.).copy()
    coefficients=original.system.unpack(state)
    index=np.flatnonzero(np.all(original.system.basis.characters==[-1,-1],axis=1))[0]
    coefficients[21,index]=.01
    state[:original.system.dof]=original.system.pack(coefficients)
    history=DenseGHHistory();history.bounds=[(0.,.1)];history.old=[state]
    F=np.zeros((7,len(state)));F[0,original.system.dof+8]=.1
    history.coefficients=[F];history.metadata=dict(vacuum_history.metadata,
        physical_coordinate_endpoint=.1,proper_endpoint=.1)
    coupled=ParentBosonProbe(history,2);control=ParentBosonProbe(history,2,decoupled_control=True)
    A=coupled.hamiltonian(0.);B=control.hamiltonian(0.)
    cross=np.zeros_like(A,dtype=bool)
    cross[:coupled.n,:coupled.n]=coupled.fields[:,None]!=coupled.fields[None,:]
    assert np.max(abs(A[cross]))>.1
    np.testing.assert_array_equal(B[cross],0.)
    np.testing.assert_array_equal(B[~cross],A[~cross])


def test_missing_model_metadata_and_unsupported_probe_cutoffs_are_explicit(vacuum_history):
    for degree in (True,1,2.5):
        with pytest.raises(ValueError):ParentBosonProbe(vacuum_history,degree)
    empty=DenseGHHistory()
    with pytest.raises(ValueError):ParentBosonProbe(empty,2)
    copied=DenseGHHistory();copied.bounds=vacuum_history.bounds;copied.old=vacuum_history.old
    copied.coefficients=vacuum_history.coefficients;copied.metadata=dict(vacuum_history.metadata,backend='dense')
    with pytest.raises(ValueError,match='Hopf'):ParentBosonProbe(copied,2)
    copied.metadata=dict(vacuum_history.metadata,physical_coordinate_endpoint=vacuum_history.bounds[-1][1]+1.)
    with pytest.raises(ValueError,match='endpoint'):ParentBosonProbe(copied,2)
