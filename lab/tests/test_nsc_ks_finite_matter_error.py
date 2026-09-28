"""Operator-norm and coherent insertion perturbation checks."""
import pytest
pytest.importorskip('flint')
import numpy as np
from flint import arb,ctx
from recursive_horizons.nsc_ks_finite_matter_error import source_norm_upper,finite_matter_error
from recursive_horizons.nsc_ks_finite_matter_error import endpoint_norms
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter


def test_source_norm_retains_off_diagonal_coherence():
    C=np.array([[.6,.2j],[-.2j,.5]])
    with ctx.workprec(120):
        bound=source_norm_upper(C)
        assert float(bound)>=np.linalg.norm(C,2)
        assert bound>source_norm_upper(np.diag(np.diag(C)))


def test_finite_source_bound_covers_field_perturbations():
    rng=np.random.default_rng(7)
    C=np.array([[.6,.2j],[-.2j,.5]])
    F=rng.normal(size=(3,2,2))+1j*rng.normal(size=(3,2,2))
    Z=rng.normal(size=F.shape)+1j*rng.normal(size=F.shape)
    dF=1e-5*(rng.normal(size=F.shape)+1j*rng.normal(size=F.shape))
    dZ=1e-5*(rng.normal(size=F.shape)+1j*rng.normal(size=F.shape))
    empty=np.zeros((0,*F.shape),complex)
    parameters=dict(mass=1.5,angular=-2.,axial_scale=.8,radius=1.4,multiplicity=3.)
    direct=source_column_matter(F+dF,Z+dZ,C,empty,empty,**parameters)['action_gradient']-source_column_matter(F,Z,C,empty,empty,**parameters)['action_gradient']
    norm=lambda v:float(np.nextafter(np.linalg.norm(v,axis=(1,2)).max(),np.inf))
    with ctx.workprec(120):
        error=finite_matter_error(norm(F),norm(Z),norm(dF),norm(dZ),source_norm_upper(C),
            mass=1.5,absolute_angular=2.,axial_lower=.8,radius_lower=1.4,multiplicity=3.)
        assert np.all(np.max(abs(direct),axis=0)<[float(restored_upper(error['N'])),float(restored_upper(error['beta']))])
        assert not error['source_accuracy_included']


def test_endpoint_norm_retains_the_actual_envelope_momentum():
    harmonic=np.array([1,1j,-1,-1j]*2,complex)
    start=np.r_[np.zeros(2),np.stack((harmonic,2*harmonic)).ravel()]
    segment=TrajectorySegment(1.03,1.02,start,2*start,np.zeros((6,len(start))))
    field=BallFourierSegment(segment,[1.],8,1.,bits=120)
    with ctx.workprec(120):
        norm,derivative=endpoint_norms(field,[.7])
        assert norm>=arb(20).sqrt().upper()
        assert derivative>=((4*arb.pi()-arb(.7))*arb(20).sqrt()).upper()


@pytest.mark.parametrize('size',[1,3,12])
def test_source_norm_accepts_complete_coherent_source_sizes(size):
    with ctx.workprec(120):
        assert source_norm_upper(np.ones((size,size),complex))==arb(size)
        assert source_norm_upper(np.eye(size,dtype=complex))==arb(1)
