"""Independent derivatives and canonical masses in four log-radius coordinates."""
import numpy as np
import pytest
from scipy.linalg import eigh

from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    combined_record,charged_potential,shape_mass_squared,
    log_radius_potential,log_radius_gradient,log_radius_hessian,log_radius_kinetic_metric,
)


def direct(u,q):
    s=float(np.mean(u));beta=np.asarray(u)-s
    model=combined_record().scalar_sector
    return charged_potential(s,model.minimum(q),q,beta)


@pytest.mark.parametrize('u,q',[
    ([-.03,-.02,-.04,-.05],1.5),([.04,.08,-.01,.03],1.),([.14,.03,.2,.09],2.),
])
def test_coordinate_potential_gradient_and_hessian_against_independent_differences(u,q):
    u=np.array(u);value=lambda v:direct(v,q)
    assert log_radius_potential(u,q)==pytest.approx(value(u),rel=3e-13,abs=2e-14)
    hg=1e-5;gradient=np.array([(value(u+hg*e)-value(u-hg*e))/(2*hg) for e in np.eye(4)])
    np.testing.assert_allclose(log_radius_gradient(u,q),gradient,rtol=2e-7,atol=2e-8)
    def hessian(h):
        H=np.empty((4,4));zero=value(u)
        for i in range(4):
            ei=h*np.eye(4)[i];H[i,i]=(value(u+ei)-2*zero+value(u-ei))/h**2
            for j in range(i):
                ej=h*np.eye(4)[j]
                H[i,j]=H[j,i]=(value(u+ei+ej)-value(u+ei-ej)-value(u-ei+ej)+value(u-ei-ej))/(4*h*h)
        return H
    extrapolated=(4*hessian(2e-4)-hessian(4e-4))/3
    np.testing.assert_allclose(log_radius_hessian(u,q),extrapolated,rtol=3e-6,atol=5e-7)
    np.testing.assert_array_equal(log_radius_hessian(u,q),log_radius_hessian(u,q).T)


def test_full_positive_kinetic_metric_matches_independent_Einstein_reduction():
    G=log_radius_kinetic_metric()
    np.testing.assert_allclose(np.linalg.eigvalsh(G),[3.,3.,3.,21.],rtol=0,atol=8e-15)
    for velocity in (np.array([.2,-.1,.03,.04]),np.ones(4)):
        direct_norm=3*np.sum(velocity**2)+.5*(3*np.sum(velocity))**2
        assert velocity@G@velocity==pytest.approx(direct_norm,rel=2e-15)
    common=np.ones(4);shape=np.array([.3,-.1,-.15,-.05])
    assert common@G@common==84
    assert common@G@shape==pytest.approx(0.,abs=2e-15)
    assert shape@G@shape==pytest.approx(3*shape@shape,rel=2e-15)


@pytest.mark.parametrize('root',VolumeParent().vacuum_roots(),ids=lambda r:r['kind'])
def test_reference_false_vacuum_and_barrier_canonical_masses(root):
    s=root['s'];u=np.full(4,s)
    np.testing.assert_allclose(log_radius_gradient(u),0.,rtol=0,atol=5e-13)
    assert log_radius_potential(u)==pytest.approx(root['vacuum'],abs=3e-14)
    masses=eigh(log_radius_hessian(u),log_radius_kinetic_metric(),eigvals_only=True)
    expected=sorted([root['canonical_mass_squared']]+[shape_mass_squared(s)]*3)
    np.testing.assert_allclose(masses,expected,rtol=4e-13,atol=3e-13)


def test_batch_evaluation_keeps_all_coordinate_axes():
    u=np.array([[-.03,-.02,-.04,-.05],[.04,.08,-.01,.03],[.14,.03,.2,.09]])
    q=np.array([1.5,1.,2.])
    assert log_radius_gradient(u,q).shape==(3,4)
    assert log_radius_hessian(u,q).shape==(3,4,4)
    for function in (log_radius_potential,log_radius_gradient,log_radius_hessian):
        np.testing.assert_allclose(function(u,q),[function(row,float(flux)) for row,flux in zip(u,q)],rtol=3e-14,atol=3e-14)
    with pytest.raises(ValueError,match='four finite'):log_radius_potential([0.,0.,0.])
