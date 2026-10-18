"""Exact nonconstant-speed integrating-factor and partial-cell controls."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dirac_wall_integral import WallPrimitive


@pytest.mark.parametrize('order',(8,12))
def test_full_triangular_solution_with_corner_masses_and_nonconstant_speed(order):
    start,end=.1,.8
    a,d,rho0=.2,.3,.4
    def masses(t):
        mu=np.exp(t)*(t*t+a)
        return np.array([mu,mu-t**3-d])
    def rate(t):
        mu=np.exp(t)*(t*t+2*t+a)
        return np.array([mu,mu-3*t*t])
    def geometry(t):
        R=np.exp(t)
        return dict(x=R,v=R,radial_acceleration=R,beta_plus=2*R*R,k_wall=2.)
    def primitive(t):
        return -np.exp(-2*t)*(.5*t**3+1.25*t*t+1.25*t+.625+(a+d)/2)
    initial=[*masses(start),rho0]
    model=WallPrimitive(geometry,rate,start=start,end=end,breakpoints=[.2,.4,.6],
        initial=initial,order=order)
    queries=np.array([start,.17,.2,.36,.4,.65,end])
    result=model(queries)
    assert result[:2]==pytest.approx(masses(queries),abs=3e-15)
    exact=np.exp(queries)*(rho0/np.exp(start)+primitive(queries)-primitive(start))
    assert result[2]==pytest.approx(exact,abs=3e-15)
    # A shared smooth-control worldtube needs no new radius equation.
    shared=WallPrimitive(geometry,rate,start=start,end=end,breakpoints=[.2,.4,.6],
        initial=[*masses(start),0.],with_displacement=False,order=order)
    assert shared(queries)[:2]==pytest.approx(masses(queries),abs=3e-15)
    assert np.all(shared(queries)[2]==0.)
    with pytest.raises(ValueError,match='interval'):
        model(end+.01)
