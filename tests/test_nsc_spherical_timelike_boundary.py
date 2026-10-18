"""Independent covariant and canonical controls; no archived state is read."""
import time

import pytest
import sympy as sp

from recursive_horizons.nsc_spherical_feedback_action import (
    alpha_of, feedback_F, feedback_V, feedback_Z, first_order_density,
)
from recursive_horizons.nsc_spherical_timelike_boundary import (
    timelike_completion, timelike_gravitational_response,
)


ACTIVE_CPU = [0.0]


@pytest.fixture(autouse=True)
def bounded_controls():
    started = time.process_time()
    yield
    ACTIVE_CPU[0] += time.process_time() - started
    assert ACTIVE_CPU[0] < 60, "timelike manufactured controls exhausted their CPU allowance"


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0


@pytest.fixture(scope="module")
def covariant():
    t, x = sp.symbols("t x", real=True)
    L, Q, r, chi = [sp.Function(name)(t, x) for name in ("L", "Q", "r", "chi")]
    A, CW, CF, flux = sp.symbols("A C_W C_F flux", nonzero=True, real=True)
    g = sp.diag(L**2, -Q**2)
    inv = g.inv()
    coords = (t, x)
    gamma = [[[sum(inv[a, d]*(sp.diff(g[d, c], coords[b]) + sp.diff(g[d, b], coords[c])
                              - sp.diff(g[b, c], coords[d])) for d in range(2))/2
               for c in range(2)] for b in range(2)] for a in range(2)]
    ricci = sp.Matrix(2, 2, lambda a, b: sum(
        sp.diff(gamma[c][a][b], coords[c]) - sp.diff(gamma[c][a][c], coords[b])
        + sum(gamma[c][c][d]*gamma[d][a][b] - gamma[c][b][d]*gamma[d][a][c]
              for d in range(2)) for c in range(2)))
    scalar = sp.simplify(sum(inv[a, b]*ricci[a, b] for a in range(2) for b in range(2)))
    F, Z = feedback_F(r, chi, A, CW), feedback_Z(A)
    bulk = L*Q*(F*scalar + Z*(sp.diff(r,t)**2/L**2 - sp.diff(r,x)**2/Q**2)
                + feedback_V(r,chi,A,CW,CF,flux))
    first = first_order_density(L,sp.diff(L,x),Q,sp.diff(Q,t),sp.diff(Q,x),
                               r,sp.diff(r,t),sp.diff(r,x),chi,sp.diff(chi,t),
                               sp.diff(chi,x),0,0,A,CW,CF,flux)
    return t,x,L,Q,r,chi,A,CW,F,scalar,sp.expand(bulk),first


def test_independent_ricci_and_both_boundary_fluxes(covariant):
    t,x,L,Q,r,chi,A,CW,F,scalar,bulk,first = covariant
    curvature_density = -2*sp.diff(sp.diff(Q,t)/L,t) + 2*sp.diff(sp.diff(L,x)/Q,x)
    zero(L*Q*scalar-curvature_density)
    zero(sp.diff(bulk,chi)-2*alpha_of(CW)*L*Q*(scalar-2-chi))
    zero(bulk-first-sp.diff(-2*F*sp.diff(Q,t)/L,t)-sp.diff(2*F*sp.diff(L,x)/Q,x))
    # Independently compute the outward-normal divergence on a timelike face.
    s=sp.symbols("s", real=True)
    K=sp.diff(L*Q*(s/Q),x)/(L*Q)
    zero(-2*L*F*K + 2*s*F*sp.diff(L,x)/Q)


@pytest.mark.parametrize("orientation", [-1,1])
def test_complete_covariant_side_variation_cancels_normal_derivatives(covariant,orientation):
    t,x,L,Q,r,chi,A,CW,F,scalar,bulk,first=covariant
    # Differentiate the second-order covariant action, not the Legendre flux.
    completion=-2*F*sp.diff(L,x)/Q
    coefficients={}
    for a in (L,Q,r,chi):
        ax,axx=sp.diff(a,x),sp.diff(a,x,2)
        higher=sp.diff(bulk,axx)
        normal=sp.diff(bulk,ax)-sp.diff(higher,x)
        zero(higher+sp.diff(completion,ax))
        coefficients[a]=sp.simplify(normal+sp.diff(completion,a))
        zero(coefficients[a]-sp.diff(first,ax))
    lp,qp,rp=sp.symbols("l q R",positive=True)
    cp,lx,rx,cx=sp.symbols("c lx rx cx",real=True)
    substitutions={sp.diff(L,x):lx,sp.diff(Q,x):0,sp.diff(r,x):rx,sp.diff(chi,x):cx,
                   L:lp,Q:qp,r:rp,chi:cp}
    response=timelike_gravitational_response(lp,qp,rp,cp,lx,rx,cx,A,CW,orientation=orientation)
    for a,key in ((L,"L"),(Q,"Q"),(r,"r"),(chi,"chi")):
        zero(orientation*coefficients[a].subs(substitutions)-response["coordinate"][key])
    zero(timelike_completion(lp,qp,rp,cp,lx,A,CW,orientation=orientation)
         -orientation*completion.subs(substitutions))


def test_physical_clock_one_form_weyl_trace_and_smooth_cut_additivity():
    L,Q,r=sp.symbols("L Q r",positive=True)
    A,CW=sp.symbols("A C_W",nonzero=True,real=True)
    chi,lx,rx,cx,dN,dr,dc=sp.symbols("chi lx rx cx dN dr dc",real=True)
    plus=timelike_gravitational_response(L,Q,r,chi,lx,rx,cx,A,CW)
    c,p=plus["coordinate"],plus["physical"]
    zero(c['L']*(dN/r-L*dr/r)+c['r']*dr+c['chi']*dc-p['N']*dN-p['r']*dr-p['chi']*dc)
    zero(-plus['unreferenced_energy']*dN + r*L*plus['unreferenced_area_pressure']*8*sp.pi*r*dr
         +p['chi']*dc-(p['N']*dN+p['r']*dr+p['chi']*dc))
    # The auxiliary part is exactly traceless under physical conformal rescaling.
    energy_aux=sp.expand(plus['unreferenced_energy']).coeff(CW)*CW
    pressure_aux=sp.expand(plus['unreferenced_area_pressure']).coeff(CW)*CW
    zero(energy_aux-8*sp.pi*r**2*pressure_aux)
    assert sp.simplify(sp.diff(energy_aux,cx)) != 0
    minus=timelike_gravitational_response(L,Q,r,chi,lx,rx,cx,A,CW,orientation=-1)
    for key in p:zero(plus['physical'][key]+minus['physical'][key])
    # Free normal jets are not an on-shell exterior response: the one-form
    # need not be the differential of a local potential on this jet space.
    zero(sp.diff(c['r'],L)-sp.diff(c['L'],r)-32*sp.pi*A*rx/Q)


def test_einstein_sector_schwarzschild_and_declared_flat_reference():
    R=sp.symbols("R",positive=True)
    mass,A,CW=sp.symbols("M A C_W",positive=True)
    chi,cx=sp.symbols("chi cx",real=True)
    # Extract the A coefficient after deriving in the nonzero-Weyl chart.
    L,Q,r,lx,rx=sp.symbols("L Q r lx rx",positive=True)
    response=timelike_gravitational_response(L,Q,r,chi,lx,rx,cx,A,CW)
    f=1-2*mass/R
    # The static-face geometric interpretation assumes R>2M.
    subst={L:sp.sqrt(f)/R,Q:1/(R*sp.sqrt(f)),r:R,lx:sp.diff(sp.sqrt(f)/R,R),rx:1}
    e=sp.expand(response['unreferenced_energy']).coeff(A)*A
    p=sp.expand(response['unreferenced_area_pressure']).coeff(A)*A
    reference_action=-16*sp.pi*A*sp.Symbol('N',positive=True)*R
    referenced_e=e.subs(subst)+16*sp.pi*A*R
    referenced_p=p.subs(subst)-2*A/R
    zero(referenced_e-16*sp.pi*A*R*(1-sp.sqrt(f)))
    zero(referenced_p-2*A/R*((1-mass/R)/sp.sqrt(f)-1))
    zero(referenced_p+sp.diff(referenced_e,R)/(8*sp.pi*R))
    zero(referenced_e.subs(mass,0));zero(referenced_p.subs(mass,0))
    zero(-sp.diff(reference_action,sp.Symbol('N',positive=True))-16*sp.pi*A*R)


def test_clock_and_length_unit_changes_and_invalid_domains():
    L,Q,r,a,scale=sp.symbols("L Q r a scale",positive=True)
    A,CW=sp.symbols("A C_W",nonzero=True,real=True)
    chi,lx,rx,cx=sp.symbols("chi lx rx cx",real=True)
    response=timelike_gravitational_response(L,Q,r,chi,lx,rx,cx,A,CW)
    for key,power in (('unreferenced_energy',1),('unreferenced_area_pressure',3),
                      ('auxiliary_response_per_proper_time',1)):
        value=response[key]
        zero(value.subs({L:L/a,lx:lx/a},simultaneous=True)-value)
        zero(value.subs({r:scale*r,rx:scale*rx,A:A/scale**2},simultaneous=True)-value/scale**power)
    with pytest.raises(ValueError,match='Q'):
        timelike_gravitational_response(1,0,2,0,0,0,0,1,.1)
    with pytest.raises(ValueError,match='C_W=0'):
        timelike_gravitational_response(1,1,2,0,0,0,0,1,0)
    with pytest.raises(ValueError,match='orientation'):
        timelike_gravitational_response(1,1,2,0,0,0,0,1,.1,orientation=0)
    with pytest.raises(ValueError,match='chi_x'):
        timelike_gravitational_response(1,1,2,0,0,0,float('nan'),1,.1)
