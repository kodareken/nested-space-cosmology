"""Directional SO(3) RHS controls without any time evolution."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve


@pytest.fixture(scope='module')
def systems():
    return (SphericalValve(16),SphericalValve(16,reflection_symmetry=True),
            SphericalValve(16,reflection_symmetry=False))


def test_default_initialization_and_rhs_equal_explicit_old_symmetry(systems):
    default,symmetric,_=systems
    assert default.reflection_symmetry is True
    np.testing.assert_array_equal(default.initial,symmetric.initial)
    np.testing.assert_array_equal(default.rhs(0.,default.initial.ravel()),
                                  symmetric.rhs(0.,symmetric.initial.ravel()))


def test_symmetric_preparation_has_same_full_and_projected_equations(systems):
    _,symmetric,directional=systems
    y=symmetric.initial.copy();u=symmetric.u
    # Nonzero regular symmetric geometric and field momenta, with no ODE step.
    y[1]=.002*(1+u*u);y[2]=.01*(1-.1*u*u);y[3]=.003*(1+u*u)
    y[4]=.0002*(1+u*u);y[6]=.0003*(1-.2*u*u);y[7]=.004*u
    a=symmetric.rhs(.1,y.ravel());b=directional.rhs(.1,y.ravel())
    np.testing.assert_allclose(a,b,rtol=3e-13,atol=3e-12)


def test_odd_curvature_and_even_normal_momentum_are_not_discarded(systems):
    _,symmetric,directional=systems;u=directional.u
    # This is a local equation check, not a claim of constraint-solved data.
    y=np.zeros_like(directional.initial);y[0]=1/directional.H
    k0=.03;k1=.02;y[2]=k0+k1*u;y[5]=.02*u;y[7]=.004*(1+.3*u*u)
    y[7]=.5*(y[7]+y[7,::-1])  # Exact even input despite cosine-node roundoff.
    raw=directional.rhs(0.,y.ravel()).reshape(y.shape)
    projected=symmetric.rhs(0.,y.ravel()).reshape(y.shape)
    np.testing.assert_array_equal(raw[0],-y[2]*y[0])
    np.testing.assert_array_equal(raw[5],y[7])
    np.testing.assert_allclose(raw[2]-raw[2,::-1],12*k0*k1*u,rtol=2e-11,atol=2e-12)
    assert np.linalg.norm(raw[0]-raw[0,::-1])>.01
    assert np.linalg.norm(raw[5]+raw[5,::-1])>.01
    np.testing.assert_array_equal(projected[5],np.zeros_like(u))
    np.testing.assert_array_equal(projected[2],projected[2,::-1])
    # The default still applies exactly the original row-parity projection.
    parity=np.array([1,1,1,1,1,-1,1,-1])[:,None]
    np.testing.assert_array_equal(projected,.5*(raw+parity*raw[:,::-1]))


@pytest.mark.parametrize('method',('diagnostic','trapped_band'))
def test_reflection_specific_diagnostics_reject_directional_mode(systems,method):
    directional=systems[2]
    with pytest.raises(ValueError,match='reflection_symmetry=True'):
        if method=='diagnostic':directional.diagnostic(0.,directional.initial)
        else:directional.trapped_band(directional.initial)


@pytest.mark.parametrize('value',(None,0,1,'False',[],object()))
def test_reflection_option_requires_a_boolean(value):
    with pytest.raises(ValueError,match='boolean reflection_symmetry'):
        SphericalValve(16,reflection_symmetry=value)
