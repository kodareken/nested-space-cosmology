"""Finite propagation is a hypothesis; numerical leakage is never subtracted."""
import pytest
pytest.importorskip('flint')
import numpy as np
from flint import arb,ctx
from recursive_horizons.nsc_ks_reference_error import spectator_reference_error,homogeneous_reference_norms
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper


def test_reference_bound_retains_numerical_leakage():
    with ctx.workprec(120):
        leak=np.array([[1e-8,0],[0,2e-8]],complex)
        bound=spectator_reference_error(1.03,0.,.4,.14,.26,leak,[1.,1.],[.5,2.],arb('1e-10'))
        assert restored_upper(bound['row_reference_error_upper'])>=arb(float(2e-8))+arb('1e-10')
        assert not bound['stress_drift_subtracted']
        assert not bound['physical_source_error_included']


def test_causal_contact_is_rejected_and_reference_derivative_is_carrier_only():
    with pytest.raises(ValueError,match='backward cone'):
        spectator_reference_error(1.03,0.,.4,.01,.26,np.zeros((2,1)),[1.],[1.],0.)
    with ctx.workprec(120):
        F,Fz=homogeneous_reference_norms(np.array([[3.],[4.]]),[1.],[2.])
        assert F==5 and Fz==10


def test_cached_reference_shift_is_included_with_directed_arithmetic():
    with ctx.workprec(120):
        joint=np.array([[1.],[0.]],complex)
        cached=joint+1e-7
        bound=spectator_reference_error(1.03,0.,.4,.14,.26,np.zeros((2,1)),[1.],[1.],0.,
            reference_pair=(joint,cached))
        assert restored_upper(bound['row_reference_error_upper'])>=arb(float(1e-7))
        assert bound['reference_pair_difference_retained']
