"""Direct 4D Riemann and observer-frame checks of the tidal readout."""
import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_dirac_tidal import tidal_contrasts,mean_tidal_response


def metric_curvature(metric,coords,point):
    """Independent coordinate Christoffel/Riemann contraction at one point."""
    g=metric.subs(point)
    inv=g.inv()
    dg=[metric.diff(x).subs(point) for x in coords]
    ddg=[[metric.diff(x,y).subs(point) for y in coords] for x in coords]
    di=[-inv@d@inv for d in dg]
    gamma=[[[sum(inv[a,d]*(dg[b][d,c]+dg[c][d,b]-dg[d][b,c])/2
        for d in range(4)) for c in range(4)] for b in range(4)] for a in range(4)]
    dgamma=[[[[sum((di[m][a,d]*(dg[b][d,c]+dg[c][d,b]-dg[d][b,c])
        +inv[a,d]*(ddg[m][b][d,c]+ddg[m][c][d,b]-ddg[m][d][b,c]))/2
        for d in range(4)) for c in range(4)] for b in range(4)] for a in range(4)] for m in range(4)]
    upper=[[[[dgamma[c][a][b][d]-dgamma[d][a][b][c]
        +sum(gamma[a][c][e]*gamma[e][b][d]-gamma[a][d][e]*gamma[e][b][c]
             for e in range(4)) for d in range(4)] for c in range(4)]
             for b in range(4)] for a in range(4)]
    ricci=sp.Matrix(4,4,lambda b,d:sum(upper[a][b][a][d] for a in range(4)))
    scalar=sum(inv[a,b]*ricci[a,b] for a in range(4) for b in range(4))
    einstein=ricci-g*scalar/2
    scales=[sp.sqrt(-g[0,0]),sp.sqrt(g[1,1]),sp.sqrt(g[2,2]),sp.sqrt(g[3,3])]
    lower=np.empty((4,4,4,4),dtype=complex)
    for a,b,c,d in np.ndindex(4,4,4,4):
        lower[a,b,c,d]=complex(sum(g[a,e]*upper[e][b][c][d] for e in range(4))
            /(scales[a]*scales[b]*scales[c]*scales[d]))
    G=np.array([[complex(einstein[a,b]/(scales[a]*scales[b]))
        for b in range(4)] for a in range(4)])
    return g,lower,G


def test_contrasts_from_full_four_dimensional_curvature_and_radial_boost():
    t,x,theta,phi=sp.symbols('t x theta phi',real=True)
    coords=(t,x,theta,phi)
    point={t:0,x:0,theta:sp.pi/2,phi:0}
    L=sp.log(sp.Rational(3,2))+2*t/7-x/5+t*t/17+t*x/11-x*x/23
    r=sp.Rational(4,5)+t/3+x/4+t*t/13+t*x/19+x*x/29
    metric=sp.diag(-sp.exp(L),sp.exp(L),r*r,r*r*sp.sin(theta)**2)
    g,R,G=metric_curvature(metric,coords,point)
    lam=3.
    T=(G+lam*np.diag([-1,1,1,1])).real
    radius=.8
    mass_total=radius*(1+(float(r.diff(t).subs(point))**2
        -float(r.diff(x).subs(point))**2)/float(sp.exp(L).subs(point)))/2
    mass=mass_total-lam*radius**3/6
    boost=np.eye(4);boost[0,0]=boost[1,1]=1.25;boost[0,1]=boost[1,0]=.75
    # Both an unboosted and a different physical radial observer.
    for frame in (np.eye(4),boost):
        curvature=np.einsum('ae,bf,cg,dh,efgh->abcd',frame,frame,frame,frame,R).real
        matter=frame@T@frame.T
        expected=tidal_contrasts(radius=radius,mass_excess=mass,
            energy_density=matter[0,0],radial_pressure=matter[1,1],angular_pressure=matter[2,2])
        assert expected['temporal']==pytest.approx(
            curvature[0,1,0,1]-curvature[0,2,0,2],abs=2e-13)
        assert expected['spatial']==pytest.approx(
            curvature[1,2,1,2]-curvature[2,3,2,3],abs=2e-13)
    assert abs(T[0,1])>.01  # This control includes radial energy flux.


def test_de_sitter_Newtonian_metric_gives_tidal_potential_combinations():
    eta,chi,theta,phi=sp.symbols('eta chi theta phi',real=True)
    coords=(eta,chi,theta,phi)
    point={eta:-2,chi:sp.Rational(2,5),theta:sp.pi/2,phi:0}
    a=-1/eta
    Phi=eta*chi**4/7
    Psi=(eta**2+1)*chi**3/11
    e=sp.I*sp.Rational(1,10**7)
    metric=a*a*sp.diag(-(1+2*e*Phi),1-2*e*Psi,
        chi*chi*(1-2*e*Psi),chi*chi*sp.sin(theta)**2*(1-2*e*Psi))
    _,R,_=metric_curvature(metric,coords,point)
    A=R[0,1,0,1]-R[0,2,0,2]
    B=R[1,2,1,2]-R[2,3,2,3]
    for actual,potential in ((A,Phi),(B,Psi)):
        expected=float(((potential.diff(chi,2)-potential.diff(chi)/chi)/(a*a)).subs(point))
        assert actual.imag/1e-7==pytest.approx(expected,abs=3e-11)


def test_mean_radius_term_and_ratio_against_exact_charged_family():
    r0,M,lam=.53,.0001,3.
    s=.02
    e=1e-7
    def family(q):
        r=r0+q*s
        # Equal-charge exact mass and Maxwell stress, coupling=1.
        return tidal_contrasts(radius=r,mass_excess=M-q/(2*r),
            energy_density=q/r**4,radial_pressure=-q/r**4,angular_pressure=q/r**4)
    plus,minus=family(e),family(-e)
    exact=mean_tidal_response(radius=r0,background_mass=M,mass_change=-1/(2*r0),
        radius_change=s,stress=dict(energy_density=1/r0**4,
            radial_pressure=-1/r0**4,angular_pressure=1/r0**4))
    assert exact['temporal_change']==pytest.approx((plus['temporal']-minus['temporal'])/(2*e),rel=2e-10)
    assert exact['spatial_change']==pytest.approx((plus['spatial']-minus['spatial'])/(2*e),rel=2e-10)
    # The ratio is differentiated at q=0; complex step avoids a near-pole error.
    analytic=family(1e-12j)
    ratio=analytic['spatial']/analytic['temporal']
    assert exact['ratio_change']==pytest.approx(ratio.imag/1e-12,rel=2e-12)
