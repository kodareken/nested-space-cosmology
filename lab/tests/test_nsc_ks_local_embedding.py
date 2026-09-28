"""Strict cone geometry and distinction from numerical-grid causality."""
from fractions import Fraction as Q
import pytest
from recursive_horizons.nsc_ks_local_embedding import continuum_local_embedding


def arguments():
    return dict(rho_up=Q(103,100),rho_sigma=1,axial_lower=Q(4,5),
        period_left=Q(-205,1024),period_length=Q(205,512),
        physical_interval=(Q(-3,100),Q(3,100)),axial_support=(Q(-3,50),Q(3,50)),radius_positive=True)


def test_owned_local_interval_has_strict_continuum_buffer():
    r=continuum_local_embedding(**arguments())
    assert Q(r['characteristic_distance_upper'])==Q(3,64)
    assert all(Q(v)>0 for v in r['two_distance_buffer_margin'])
    assert r['continuum_local_embedding']
    assert not r['finite_Fourier_causality_certified']
    assert r['physical_local_gate']=='OPEN'
    assert not r['global_matching_required']


def test_insufficient_buffer_and_unproved_radius_are_rejected():
    args=arguments();args.update(period_left=Q(-1,10),period_length=Q(1,5))
    with pytest.raises(ValueError,match='buffer'):continuum_local_embedding(**args)
    args=arguments();args['radius_positive']=False
    with pytest.raises(ValueError,match='positivity'):continuum_local_embedding(**args)
    args=arguments();args['axial_lower']=Q(9,10)
    with pytest.raises(ValueError,match='conservative'):continuum_local_embedding(**args)
