import numpy as np
import pytest
from flint import arb,ctx
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_source_occupation_enclosure import occupation_vacuum_distance


def test_exact_horizon_block_identity():
    import sympy as sp
    x=sp.symbols('x',positive=True);f=1/(1+x*x);s=x/(1+x*x)
    D=sp.Matrix([[f,-sp.I*s],[sp.I*s,-f]])
    assert sp.simplify(D*D-f*sp.eye(2))==sp.zeros(2)


def test_complete_source_and_normalized_sewings_with_negative_partners():
    for E,m,k,w in ((.5,.1,1.,4.),(.2,1.,.3,4.),(2.071557525252559,np.pi/2,.23832579963401956,3.973074368754331)):
        C=source_covariance(E,k,w,m);P=np.diag([0.,1.,0.])
        # Form the analytic difference directly: subtracting stored 1-f from
        # 1 loses f at tiny occupations and is a separate rounding error.
        f=C[0,0].real;s=abs(C[0,1]);n=C[2,2].real
        D=np.array([[f,-1j*s,0],[1j*s,-f,0],[0,0,n]],complex)
        bound=float(occupation_vacuum_distance(E,m,k,w)['operator_distance_upper'])
        assert np.linalg.norm(D,2)<=bound*(1+2e-14)
        for probability in (0.,.3,1.):
            R=np.sqrt(probability)*np.exp(.37j)
            A=np.array([[0,1,0],[R,0,np.sqrt(1-probability)]],complex)
            assert np.max(abs(A@A.conj().T-np.eye(2)))<1e-15
            assert np.linalg.norm(A@D@A.conj().T,2)<=bound*(1+2e-14)
        negative=source_covariance(-E,k,w,m)-(np.eye(3)-P)
        assert np.max(abs(negative+D.conj()))<2e-16


def test_coherence_and_open_incoming_terms_cannot_be_dropped():
    E,m,k,w=.5,.1,1.,4.
    C=source_covariance(E,k,w,m);P=np.diag([0.,1.,0.]);D=C-P
    f=C[0,0].real;n=C[2,2].real
    assert np.linalg.norm(D[:2,:2],2)>f
    assert n>np.sqrt(f)
    assert np.linalg.norm(D,2)>np.sqrt(f)
    result=occupation_vacuum_distance(E,m,k,w)
    assert result['incoming_gap_used'] is False
    assert float(result['operator_distance_upper'])>=n*(1-2e-15)


def test_original_gap_and_interval_crossing_are_distinct():
    with ctx.workprec(192):
        r=occupation_vacuum_distance(2.071557525252559,arb.pi()/2,.23832579963401956,3.973074368754331)
        assert r['incoming_gap_used'] is True and r['incoming_occupation_upper']==0
        assert r['operator_distance_upper']<arb('1.383e-12')
        crossing=occupation_vacuum_distance(arb('.4','.01'),.1,1,4)
        assert crossing['incoming_gap_used'] is False and crossing['incoming_occupation_upper']>0
        for missing in (None,False):
            with pytest.raises(ValueError):occupation_vacuum_distance(missing,1,1,1)
        with pytest.raises(ValueError):occupation_vacuum_distance(1,-1,1,1)


def test_analytic_bound_is_not_a_rounding_bound_for_stored_covariance():
    E,m,k,w=2.071557525252559,np.pi/2,.23832579963401956,3.973074368754331
    C=source_covariance(E,k,w,m)
    bound=float(occupation_vacuum_distance(E,m,k,w)['operator_distance_upper'])
    assert C[1,1]==1 and C[0,0]>0  # binary 1-f rounded to one
    assert np.linalg.norm(C-np.diag([0.,1.,0.]),2)>bound
