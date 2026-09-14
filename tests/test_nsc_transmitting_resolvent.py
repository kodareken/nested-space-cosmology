"""Checks of the domain-correct response, not legacy source calculations."""
import numpy as np
import pytest
from recursive_horizons.nsc_transmitting_resolvent import (
    TransmittingCrossResolvent, probes, quadrature, inverse_response,
)


def test_spatial_packet_embedding_is_isometric_and_has_disjoint_support():
    x,w=quadrature(8)
    gram=sum(weight*probes(rho).conj().T@probes(rho) for rho,weight in zip(x,w))
    assert np.max(abs(gram-np.eye(4))) < 2e-14


def test_real_energy_and_nontrapped_deformation_rejected():
    owner=TransmittingCrossResolvent(0.,0.)
    with pytest.raises(ValueError,match='nonreal'):owner.solve(.2)
    with pytest.raises(ValueError,match='trapped'):owner.solve(.2+.25j,parameters=[1.,0.,0.,0.])


def test_inverse_response_derivative_is_not_a_unitary_history_jet():
    from types import SimpleNamespace
    from recursive_horizons.nsc_transmitting_resolvent import PacketResolvent
    r=np.array([[2.,.2j],[-.2j,1.]])+1j*np.eye(2)
    dr=np.array([np.array([[.1,.2],[.2,-.3]])])
    k,dk=inverse_response(SimpleNamespace(response=r,directional_jets=dr))
    assert np.max(abs(dk[0]+k@dr[0]@k)) < 1e-14
    with pytest.raises(ValueError,match='not dU'):
        PacketResolvent(None,1j,None,{},None,None,0.,0).as_endpoint_branch_jets()
